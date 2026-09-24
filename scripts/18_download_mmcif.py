"""
Stage 3B, step 1: define the coordinate-audit selection pool and download
one mmCIF file per selected unique PDB entry.

Reads:
    data/manifests/rcsb_precoordinate_qc.csv

Writes:
    data/manifests/stage3b_coordinate_pool.csv   (2072 rows)
    data/raw/mmcif/{PDB_ID}.cif.gz               (gitignored)
    data/manifests/mmcif_download_manifest.csv

Selection rule (Section 1): method_is_xray == True AND mapping_overlaps_lbd
== True, using the Stage 3A STANDARDIZED LBD interval. No resolution,
R-free, taxonomy, fusion, or prefilter restriction is applied here.

Idempotent: a cached, re-validated .cif.gz is reused (CACHED_VALID) unless
--refresh. A failed re-validation triggers a fresh download, never a silent
trust of a possibly-corrupt cache.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.mmcif_client import (  # noqa: E402
    RCSB_FILES_BASE_URL,
    download_one_mmcif,
    load_cached_mmcif,
    save_validated_download,
)

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"
PRECOORDINATE_QC_CSV = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"

POOL_CSV = MANIFESTS_DIR / "stage3b_coordinate_pool.csv"
DOWNLOAD_MANIFEST_CSV = MANIFESTS_DIR / "mmcif_download_manifest.csv"

MAX_CONCURRENT_DOWNLOADS = 4
POOL_COLUMNS = [
    "uniprot_id", "nr_code", "common_name", "group",
    "pdb_id", "polymer_entity_id", "entity_id",
    "method_is_xray", "lbd_mapping_category", "mapping_overlaps_lbd", "lbd_mapping_coverage",
    "resolution_project_range_status", "r_free_status",
    "passes_project_metadata_prefilter",
    "coordinate_audit_selected", "coordinate_audit_reason",
]
DOWNLOAD_COLUMNS = [
    "pdb_id", "download_url", "http_status", "retrieved_at_utc",
    "compressed_size_bytes", "decompressed_size_bytes",
    "sha256_gzip", "sha256_decompressed",
    "gzip_valid", "mmcif_parseable", "parsed_entry_id", "entry_id_matches_request",
    "raw_file", "download_status", "error_message",
]


def build_pool(qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in qc.iterrows():
        is_xray = bool(row["method_is_xray"])
        overlaps_lbd = bool(row["mapping_overlaps_lbd"])
        selected = is_xray and overlaps_lbd
        if selected:
            reason = "SELECTED_XRAY_LBD_OVERLAP"
        elif not is_xray and not overlaps_lbd:
            reason = "NON_XRAY_AND_NO_LBD_OVERLAP"
        elif not is_xray:
            reason = "NON_XRAY"
        else:
            reason = "NO_LBD_OVERLAP"

        rows.append(
            {
                "uniprot_id": row["uniprot_id"], "nr_code": row["nr_code"],
                "common_name": row["common_name"], "group": row["group"],
                "pdb_id": row["pdb_id"], "polymer_entity_id": row["polymer_entity_id"], "entity_id": row["entity_id"],
                "method_is_xray": is_xray, "lbd_mapping_category": row["lbd_mapping_category"],
                "mapping_overlaps_lbd": overlaps_lbd, "lbd_mapping_coverage": row["lbd_mapping_coverage"],
                "resolution_project_range_status": row["resolution_project_range_status"],
                "r_free_status": row["r_free_status"],
                "passes_project_metadata_prefilter": bool(row["passes_project_metadata_prefilter"]),
                "coordinate_audit_selected": selected,
                "coordinate_audit_reason": reason,
            }
        )
    return pd.DataFrame(rows)[POOL_COLUMNS]


def resolve_one(pdb_id: str, refresh: bool) -> dict:
    cached_path = RAW_MMCIF_DIR / f"{pdb_id.upper()}.cif.gz"

    if not refresh and cached_path.exists():
        gzip_valid, mmcif_parseable, parsed_entry_id, err, decompressed = load_cached_mmcif(cached_path)
        if gzip_valid and mmcif_parseable and parsed_entry_id.upper() == pdb_id.upper():
            gzip_bytes = cached_path.read_bytes()
            return {
                "pdb_id": pdb_id, "download_url": f"{RCSB_FILES_BASE_URL}/{pdb_id}.cif.gz",
                "http_status": 200,
                "retrieved_at_utc": datetime.fromtimestamp(cached_path.stat().st_mtime, tz=timezone.utc).isoformat(),
                "compressed_size_bytes": len(gzip_bytes), "decompressed_size_bytes": len(decompressed),
                "sha256_gzip": hashlib.sha256(gzip_bytes).hexdigest(),
                "sha256_decompressed": hashlib.sha256(decompressed).hexdigest(),
                "gzip_valid": True, "mmcif_parseable": True,
                "parsed_entry_id": parsed_entry_id, "entry_id_matches_request": True,
                "raw_file": str(cached_path.relative_to(PROJECT_ROOT)),
                "download_status": "CACHED_VALID", "error_message": "",
            }
            # else: fall through to a fresh download — never silently trust
            # a cache file that failed re-validation.

    result = download_one_mmcif(pdb_id)
    if result.download_status == "SUCCESS":
        raw_file = save_validated_download(result, RAW_MMCIF_DIR)
        return {
            "pdb_id": pdb_id, "download_url": result.download_url, "http_status": result.http_status,
            "retrieved_at_utc": result.retrieved_at_utc,
            "compressed_size_bytes": result.compressed_size_bytes,
            "decompressed_size_bytes": result.decompressed_size_bytes,
            "sha256_gzip": result.sha256_gzip, "sha256_decompressed": result.sha256_decompressed,
            "gzip_valid": result.gzip_valid, "mmcif_parseable": result.mmcif_parseable,
            "parsed_entry_id": result.parsed_entry_id, "entry_id_matches_request": result.entry_id_matches_request,
            "raw_file": str(Path(raw_file).relative_to(PROJECT_ROOT)),
            "download_status": "SUCCESS", "error_message": "",
        }

    return {
        "pdb_id": pdb_id, "download_url": result.download_url, "http_status": result.http_status or "",
        "retrieved_at_utc": result.retrieved_at_utc,
        "compressed_size_bytes": result.compressed_size_bytes, "decompressed_size_bytes": result.decompressed_size_bytes,
        "sha256_gzip": result.sha256_gzip, "sha256_decompressed": result.sha256_decompressed,
        "gzip_valid": result.gzip_valid, "mmcif_parseable": result.mmcif_parseable,
        "parsed_entry_id": result.parsed_entry_id, "entry_id_matches_request": result.entry_id_matches_request,
        "raw_file": "", "download_status": "FAILED", "error_message": result.error_message,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--only", type=str, default=None, help="Comma-separated PDB IDs to restrict to (pilot mode).")
    args = parser.parse_args()

    qc = pd.read_csv(PRECOORDINATE_QC_CSV, dtype={"pdb_id": str})
    pool_df = build_pool(qc)
    pool_df.to_csv(POOL_CSV, index=False)
    print(f"Wrote {POOL_CSV.relative_to(PROJECT_ROOT)} ({len(pool_df)} rows)")
    print(pool_df["coordinate_audit_reason"].value_counts().to_string())

    selected = pool_df[pool_df["coordinate_audit_selected"]]
    unique_pdb_ids = sorted(selected["pdb_id"].unique())
    if args.only:
        only_set = set(args.only.split(","))
        unique_pdb_ids = [p for p in unique_pdb_ids if p in only_set]

    print(f"\nSelected candidate relationships: {len(selected)}")
    print(f"Selected unique polymer entities: {selected['polymer_entity_id'].nunique()}")
    print(f"Selected unique PDB entries: {selected['pdb_id'].nunique()}")
    print(f"PDB entries to process this run: {len(unique_pdb_ids)}")

    RAW_MMCIF_DIR.mkdir(parents=True, exist_ok=True)

    download_rows = []
    start_time = time.monotonic()
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENT_DOWNLOADS) as executor:
        futures = {executor.submit(resolve_one, pdb_id, args.refresh): pdb_id for pdb_id in unique_pdb_ids}
        completed = 0
        for future in as_completed(futures):
            row = future.result()
            download_rows.append(row)
            completed += 1
            if completed % 50 == 0 or completed == len(unique_pdb_ids):
                elapsed = time.monotonic() - start_time
                print(f"  [{completed}/{len(unique_pdb_ids)}] elapsed={elapsed:.0f}s "
                      f"last={row['pdb_id']} -> {row['download_status']}")

    new_df = pd.DataFrame(download_rows)
    if DOWNLOAD_MANIFEST_CSV.exists() and args.only:
        existing = pd.read_csv(DOWNLOAD_MANIFEST_CSV)
        existing = existing[~existing["pdb_id"].isin(new_df["pdb_id"])]
        new_df = pd.concat([existing, new_df], ignore_index=True)
    new_df = new_df.sort_values("pdb_id")[DOWNLOAD_COLUMNS]
    new_df.to_csv(DOWNLOAD_MANIFEST_CSV, index=False)
    print(f"\nWrote {DOWNLOAD_MANIFEST_CSV.relative_to(PROJECT_ROOT)} ({len(new_df)} rows)")

    status_counts = new_df["download_status"].value_counts().to_dict()
    print(f"download_status counts: {status_counts}")
    total_compressed = int(new_df.loc[new_df["download_status"] != "FAILED", "compressed_size_bytes"].sum())
    total_decompressed = int(new_df.loc[new_df["download_status"] != "FAILED", "decompressed_size_bytes"].sum())
    print(f"total compressed bytes: {total_compressed}  total decompressed bytes: {total_decompressed}")
    failed = new_df[new_df["download_status"] == "FAILED"]
    if len(failed):
        print(f"FAILED PDB IDs: {failed['pdb_id'].tolist()}")


if __name__ == "__main__":
    main()
