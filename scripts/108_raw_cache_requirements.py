"""
Raw-cache requirements for reproducing v2 from scratch.

The raw caches are gitignored. This script records, in a committed manifest, every raw
file the v2 rebuild reads, where it came from and how to verify it, and can restore and
verify a teammate's local caches.

Caches
    data/raw/mmcif/<PDB>.cif.gz                     1,777 release files (frozen v1 cache)
    data/interim/raw_mmcif_v2_additions/<PDB>.cif.gz  1KB2, 1KB4, 1KB6 (approved v2 additions;
                                                    kept separate so the release cache stays 1,777)
    data/raw/api/uniprot/<ACC>.json                 48 UniProt entries

Manifest: data/rebuild_v2/inputs/raw_cache_requirements_v2.csv
    mmCIF: source URL, retrieval time, SHA-256 of the gzip file and of the decompressed
           mmCIF (from mmcif_download_manifest_v2.csv);
    UniProt: source URL, retrieval time, SHA-256 of the cached JSON (uniprot_api_manifest.csv),
           plus entry version, sequence version, sequence length and SHA-256 of the sequence.

Verification levels (a re-download can differ in bytes even when the content is the same):
    EXACT            file SHA-256 equals the recorded one
    CONTENT_MATCH    mmCIF: decompressed SHA-256 equal (gzip re-encoded);
                     UniProt: sequence SHA-256 equal (entry metadata updated)
    CONTENT_CHANGED  the content itself differs: v2 may not reproduce exactly; report it
    MISSING

Usage:
    python scripts/108_raw_cache_requirements.py --write-manifest   (from verified local caches)
    python scripts/108_raw_cache_requirements.py --check
    python scripts/108_raw_cache_requirements.py --download         (network: RCSB, UniProt REST)
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import time
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MMCIF_MANIFEST = PROJECT_ROOT / "data" / "rebuild_v2" / "inputs" / "mmcif_download_manifest_v2.csv"
MMCIF_ADDITIONS = PROJECT_ROOT / "data" / "rebuild_v2" / "inputs" / "mmcif_download_additions_v2.csv"
UNIPROT_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "uniprot_api_manifest.csv"
OUT = PROJECT_ROOT / "data" / "rebuild_v2" / "inputs" / "raw_cache_requirements_v2.csv"
RELEASE_CACHE = "data/raw/mmcif"
ADDITIONS_CACHE = "data/interim/raw_mmcif_v2_additions"


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def uniprot_content(raw: bytes) -> dict:
    e = json.loads(raw)
    seq = e["sequence"]["value"]
    return {"entry_version": e.get("entryAudit", {}).get("entryVersion"),
            "sequence_version": e.get("entryAudit", {}).get("sequenceVersion"),
            "sequence_length": len(seq), "sha256_content": sha(seq.encode())}


def write_manifest() -> pd.DataFrame:
    m = pd.read_csv(MMCIF_MANIFEST, dtype={"pdb_id": str}, keep_default_na=False)
    additions = set(pd.read_csv(MMCIF_ADDITIONS, dtype={"pdb_id": str})["pdb_id"])
    rows = []
    for r in m.itertuples():
        cache = ADDITIONS_CACHE if r.pdb_id in additions else RELEASE_CACHE
        rows.append({"cache": "mmcif", "id": r.pdb_id, "local_path": f"{cache}/{r.pdb_id}.cif.gz",
                     "source_url": r.download_url, "retrieved_at_utc": r.retrieved_at_utc,
                     "sha256_file": r.sha256_gzip, "sha256_content": r.sha256_decompressed,
                     "content_definition": "sha256 of the decompressed mmCIF",
                     "entry_version": "", "sequence_version": "", "sequence_length": "",
                     "provenance": ("approved v2 addition (VDR DBD entries admitted by the final span), downloaded "
                                    "2026-10-21" if r.pdb_id in additions else "v1 release cache (script 18)")})
    u = pd.read_csv(UNIPROT_MANIFEST, keep_default_na=False)
    for r in u.itertuples():
        path = PROJECT_ROOT / r.raw_file
        raw = path.read_bytes()
        if sha(raw) != r.sha256:
            raise SystemExit(f"local {r.raw_file} does not match uniprot_api_manifest.csv; refusing to record it")
        rows.append({"cache": "uniprot", "id": r.requested_accession, "local_path": r.raw_file,
                     "source_url": r.request_url, "retrieved_at_utc": r.retrieved_at_utc,
                     "sha256_file": r.sha256, "content_definition": "sha256 of the canonical sequence string",
                     "provenance": "UniProt REST (script 02)", **uniprot_content(raw)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False, lineterminator="\n")
    return df


def verify(row) -> str:
    p = PROJECT_ROOT / row.local_path
    if not p.exists():
        return "MISSING"
    raw = p.read_bytes()
    if sha(raw) == row.sha256_file:
        return "EXACT"
    try:
        content = sha(gzip.decompress(raw)) if row.cache == "mmcif" else uniprot_content(raw)["sha256_content"]
    except Exception:  # unreadable file counts as changed content
        return "CONTENT_CHANGED"
    return "CONTENT_MATCH" if content == row.sha256_content else "CONTENT_CHANGED"


def check() -> pd.Series:
    req = pd.read_csv(OUT, dtype={"id": str}, keep_default_na=False)
    req["status"] = [verify(r) for r in req.itertuples()]
    print(req.groupby(["cache", "status"]).size().to_string())
    bad = req[req["status"].isin(["MISSING", "CONTENT_CHANGED"])]
    if len(bad):
        print(bad[["cache", "id", "local_path", "status"]].to_string(index=False))
    return req["status"]


def download() -> None:
    import requests

    from scripts.utils.mmcif_client import download_one_mmcif, save_validated_download
    req = pd.read_csv(OUT, dtype={"id": str}, keep_default_na=False)
    for r in req.itertuples():
        p = PROJECT_ROOT / r.local_path
        if p.exists():
            continue
        if r.cache == "mmcif":  # project client: retries, gzip/mmCIF validation, atomic write
            res = download_one_mmcif(r.id)
            if res.download_status != "SUCCESS":
                print(f"FAILED {r.id}: {res.error_message}")
                continue
            save_validated_download(res, p.parent)
        else:
            resp = requests.get(r.source_url, timeout=60)
            if resp.status_code != 200:
                print(f"FAILED {r.id}: HTTP {resp.status_code}")
                continue
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(resp.content)
        time.sleep(0.1)
    check()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write-manifest", action="store_true")
    g.add_argument("--check", action="store_true")
    g.add_argument("--download", action="store_true")
    a = ap.parse_args()
    if a.write_manifest:
        df = write_manifest()
        print(f"wrote {OUT.relative_to(PROJECT_ROOT)}: {df.groupby('cache').size().to_dict()}")
        check()
    elif a.check:
        check()
    else:
        download()


if __name__ == "__main__":
    main()
