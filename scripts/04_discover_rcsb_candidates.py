"""
Stage 2, step 1: discover current EXPERIMENTAL RCSB polymer entities for
every validated receptor, via the official RCSB Search API.

Reads:
    config/nr_metadata.py   (validated 48-record working master)

Writes:
    data/raw/api/rcsb/search/{uniprot_id}.json   (raw responses, gitignored)
    data/manifests/rcsb_search_manifest.csv       (one row per receptor query)
    data/manifests/rcsb_discovery_by_receptor.csv (one row per receptor x entity)

Idempotent by default: a cached raw response is reused unless --refresh is
passed. Use --only P11473,P10275 to restrict to specific accessions (e.g.
for the Stage 2 pilot).

This script performs NO Data API calls and downloads NO coordinate files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402
from scripts.utils.rcsb_search_client import (  # noqa: E402
    RCSB_SEARCH_URL,
    RCSBSearchResult,
    search_polymer_entities_by_uniprot,
)
from scripts.utils.rcsb_extract import parse_polymer_entity_identifier  # noqa: E402

RAW_SEARCH_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "search"
MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"

SEARCH_MANIFEST_PATH = MANIFESTS_DIR / "rcsb_search_manifest.csv"
DISCOVERY_BY_RECEPTOR_PATH = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"

SEARCH_MANIFEST_COLUMNS = [
    "record_index", "uniprot_id", "nr_code", "common_name", "group",
    "request_url", "request_method", "query_return_type", "experimental_only",
    "http_status", "retrieved_at_utc",
    "result_count_reported", "result_count_parsed",
    "raw_response_file", "raw_response_sha256",
    "success", "error_message",
]

DISCOVERY_COLUMNS = [
    "uniprot_id", "nr_code", "common_name", "group",
    "query_uniprot_id", "polymer_entity_id", "pdb_id", "entity_id",
    "query_retrieved_at_utc",
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def utc_mtime_iso(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def load_cached_response(accession: str) -> dict | None:
    raw_path = RAW_SEARCH_DIR / f"{accession}.json"
    if not raw_path.exists():
        return None
    try:
        return json.loads(raw_path.read_text())
    except json.JSONDecodeError:
        return None


def resolve_one(accession: str, refresh: bool) -> dict:
    """Fetch (or reuse cached) search results for one accession. The saved
    wrapper preserves BOTH the exact request payload and the raw response,
    per Stage 2 spec option A."""
    raw_path = RAW_SEARCH_DIR / f"{accession}.json"

    if not refresh:
        cached = load_cached_response(accession)
        if cached is not None and "response" in cached and "request" in cached:
            raw_bytes = raw_path.read_bytes()
            response = cached["response"]
            identifiers = [r["identifier"] for r in response.get("result_set", [])]
            return {
                "request_url": cached["request"]["url"],
                "http_status": cached.get("http_status", 200),
                "retrieved_at_utc": utc_mtime_iso(raw_path),
                "raw_file": str(raw_path.relative_to(PROJECT_ROOT)),
                "sha256": sha256_bytes(raw_bytes),
                "success": True,
                "error_message": "",
                "total_count": response.get("total_count", len(identifiers)),
                "identifiers": identifiers,
                "_from_cache": True,
            }

    result: RCSBSearchResult = search_polymer_entities_by_uniprot(accession)

    if result.success:
        RAW_SEARCH_DIR.mkdir(parents=True, exist_ok=True)
        wrapper = {
            "request": {"url": result.request_url, "method": "POST", "payload": result.request_payload},
            "http_status": result.http_status,
            "retrieved_at_utc": result.retrieved_at_utc,
            "response": result.raw_json,
        }
        raw_bytes = json.dumps(wrapper, indent=2, ensure_ascii=False).encode("utf-8")
        raw_path.write_bytes(raw_bytes)
        return {
            "request_url": result.request_url,
            "http_status": result.http_status,
            "retrieved_at_utc": utc_mtime_iso(raw_path),
            "raw_file": str(raw_path.relative_to(PROJECT_ROOT)),
            "sha256": sha256_bytes(raw_bytes),
            "success": True,
            "error_message": "",
            "total_count": result.total_count,
            "identifiers": result.identifiers,
            "_from_cache": False,
        }

    return {
        "request_url": result.request_url,
        "http_status": result.http_status if result.http_status is not None else "",
        "retrieved_at_utc": result.retrieved_at_utc,
        "raw_file": "",
        "sha256": "",
        "success": False,
        "error_message": result.error_message,
        "total_count": None,
        "identifiers": [],
        "_from_cache": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Force re-querying every accession from the live API.")
    parser.add_argument(
        "--only", type=str, default=None,
        help="Comma-separated list of UniProt accessions to restrict processing to (e.g. for a pilot run).",
    )
    args = parser.parse_args()

    only_set = set(args.only.split(",")) if args.only else None
    records = [r for r in NR_METADATA if only_set is None or r["uniprot_id"] in only_set]

    RAW_SEARCH_DIR.mkdir(parents=True, exist_ok=True)
    MANIFESTS_DIR.mkdir(parents=True, exist_ok=True)

    search_manifest_rows = []
    discovery_rows = []

    print(f"Discovering current experimental RCSB polymer entities for {len(records)} receptor(s) "
          f"via {RCSB_SEARCH_URL} ...")

    # Map original record_index (position in the full 48-record NR_METADATA)
    # so a --only pilot run still reports a meaningful index.
    full_index_by_uniprot = {r["uniprot_id"]: i for i, r in enumerate(NR_METADATA)}

    for record in records:
        accession = record["uniprot_id"]
        record_index = full_index_by_uniprot[accession]
        result = resolve_one(accession, refresh=args.refresh)
        source = "cache" if result.get("_from_cache") else ("live" if result["success"] else "FAILED")
        n_results = result["total_count"] if result["success"] else "N/A"
        print(f"  [{record_index+1:2d}/48] {accession:10s} common_name={record['common_name']:10s} "
              f"-> {source}  n_results={n_results}")

        search_manifest_rows.append(
            {
                "record_index": record_index,
                "uniprot_id": accession,
                "nr_code": record["nr_code"],
                "common_name": record["common_name"],
                "group": record["group"],
                "request_url": result["request_url"],
                "request_method": "POST",
                "query_return_type": "polymer_entity",
                "experimental_only": True,
                "http_status": result["http_status"],
                "retrieved_at_utc": result["retrieved_at_utc"],
                "result_count_reported": result["total_count"] if result["success"] else "",
                "result_count_parsed": len(result["identifiers"]) if result["success"] else "",
                "raw_response_file": result["raw_file"],
                "raw_response_sha256": result["sha256"],
                "success": result["success"],
                "error_message": result["error_message"],
            }
        )

        if result["success"]:
            for identifier in result["identifiers"]:
                pdb_id, entity_id = parse_polymer_entity_identifier(identifier)
                discovery_rows.append(
                    {
                        "uniprot_id": accession,
                        "nr_code": record["nr_code"],
                        "common_name": record["common_name"],
                        "group": record["group"],
                        "query_uniprot_id": accession,
                        "polymer_entity_id": identifier,
                        "pdb_id": pdb_id,
                        "entity_id": entity_id,
                        "query_retrieved_at_utc": result["retrieved_at_utc"],
                    }
                )

    # When --only restricts the run, merge into any existing full manifests
    # rather than truncating them, so a pilot run doesn't destroy prior
    # full-run results. On a full run (no --only), this naturally rewrites
    # the whole table since `records` covers all 48.
    search_df = pd.DataFrame(search_manifest_rows)[SEARCH_MANIFEST_COLUMNS]
    if only_set is not None and SEARCH_MANIFEST_PATH.exists():
        existing = pd.read_csv(SEARCH_MANIFEST_PATH)
        existing = existing[~existing["uniprot_id"].isin(only_set)]
        search_df = pd.concat([existing, search_df], ignore_index=True).sort_values("record_index")
    search_df.to_csv(SEARCH_MANIFEST_PATH, index=False)
    print(f"\nWrote {SEARCH_MANIFEST_PATH.relative_to(PROJECT_ROOT)} ({len(search_df)} rows)")

    discovery_df = pd.DataFrame(discovery_rows, columns=DISCOVERY_COLUMNS)
    if only_set is not None and DISCOVERY_BY_RECEPTOR_PATH.exists():
        existing = pd.read_csv(DISCOVERY_BY_RECEPTOR_PATH)
        existing = existing[~existing["uniprot_id"].isin(only_set)]
        discovery_df = pd.concat([existing, discovery_df], ignore_index=True)
    discovery_df.to_csv(DISCOVERY_BY_RECEPTOR_PATH, index=False)
    print(f"Wrote {DISCOVERY_BY_RECEPTOR_PATH.relative_to(PROJECT_ROOT)} ({len(discovery_df)} rows)")

    n_success = int(search_df["success"].sum())
    n_failed = len(search_df) - n_success
    n_zero = int(((search_df["success"]) & (search_df["result_count_parsed"] == 0)).sum())
    print(f"\n--- Discovery summary ({'PILOT: ' + args.only if only_set else 'FULL 48'}) ---")
    print(f"queried: {len(search_df)}  successful: {n_success}  failed: {n_failed}  zero-result: {n_zero}")
    print(f"total receptor->entity relationships: {len(discovery_df)}")
    print(f"unique polymer entities: {discovery_df['polymer_entity_id'].nunique() if len(discovery_df) else 0}")
    print(f"unique PDB entries: {discovery_df['pdb_id'].nunique() if len(discovery_df) else 0}")


if __name__ == "__main__":
    main()
