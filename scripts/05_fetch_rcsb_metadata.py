"""
Stage 2, step 2: enrich the discovered polymer entities and their parent PDB
entries with metadata from the official RCSB Data API (GraphQL), in
conservative batches.

Reads:
    data/manifests/rcsb_discovery_by_receptor.csv

Writes:
    data/raw/api/rcsb/data/polymer_entities_batch_XXXX.json   (gitignored)
    data/raw/api/rcsb/data/entries_batch_XXXX.json            (gitignored)
    data/manifests/rcsb_data_api_manifest.csv

Idempotent by default: if a valid cached batch file set already covers all
required IDs, it is reused. Pass --refresh to force re-fetching everything.

This script downloads metadata only — no coordinate files.
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

from scripts.utils.rcsb_data_client import (  # noqa: E402
    DEFAULT_BATCH_SIZE,
    RCSB_DATA_API_URL,
    fetch_batches,
)

RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "data"
MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
DISCOVERY_BY_RECEPTOR_PATH = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
DATA_API_MANIFEST_PATH = MANIFESTS_DIR / "rcsb_data_api_manifest.csv"

DATA_API_MANIFEST_COLUMNS = [
    "batch_type", "batch_number", "number_requested",
    "request_url", "http_status", "retrieved_at_utc",
    "raw_response_file", "raw_response_sha256", "success", "error_message",
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_manifest_of_cached_ids(batch_type: str) -> set[str]:
    """IDs already covered by successful cached batches of this type,
    recovered by reading back the saved raw batch files."""
    covered: set[str] = set()
    prefix = "polymer_entities" if batch_type == "polymer_entity" else "entries"
    for path in sorted(RAW_DATA_DIR.glob(f"{prefix}_batch_*.json")):
        try:
            wrapper = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        covered.update(wrapper.get("ids_requested", []))
    return covered


def run_batches_for(batch_type: str, ids: list[str], refresh: bool) -> list[dict]:
    prefix = "polymer_entities" if batch_type == "polymer_entity" else "entries"

    if not refresh:
        cached_ids = load_manifest_of_cached_ids(batch_type)
        missing = [i for i in ids if i not in cached_ids]
        if not missing:
            print(f"  {batch_type}: all {len(ids)} IDs already covered by cached batches — reusing.")
            manifest_rows = []
            for path in sorted(RAW_DATA_DIR.glob(f"{prefix}_batch_*.json")):
                wrapper = json.loads(path.read_text())
                raw_bytes = path.read_bytes()
                manifest_rows.append(
                    {
                        "batch_type": batch_type,
                        "batch_number": wrapper["batch_number"],
                        "number_requested": len(wrapper["ids_requested"]),
                        "request_url": RCSB_DATA_API_URL,
                        "http_status": wrapper["http_status"],
                        "retrieved_at_utc": datetime.fromtimestamp(
                            path.stat().st_mtime, tz=timezone.utc
                        ).isoformat(),
                        "raw_response_file": str(path.relative_to(PROJECT_ROOT)),
                        "raw_response_sha256": sha256_bytes(raw_bytes),
                        "success": True,
                        "error_message": "",
                    }
                )
            return manifest_rows
        print(f"  {batch_type}: {len(missing)} of {len(ids)} IDs not yet cached — fetching those.")
        ids_to_fetch = missing
    else:
        ids_to_fetch = ids
        # Clear old batch files for this type so we never mix stale and
        # freshly refreshed batches under the same numbering scheme.
        for path in RAW_DATA_DIR.glob(f"{prefix}_batch_*.json"):
            path.unlink()

    RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)
    results = fetch_batches(batch_type, ids_to_fetch, batch_size=DEFAULT_BATCH_SIZE)

    manifest_rows = []
    for r in results:
        raw_path = RAW_DATA_DIR / f"{prefix}_batch_{r.batch_number:04d}.json"
        wrapper = {
            "batch_type": batch_type,
            "batch_number": r.batch_number,
            "ids_requested": r.ids_requested,
            "retrieved_at_utc": r.retrieved_at_utc,
            "http_status": r.http_status,
            "success": r.success,
            "error_message": r.error_message,
            "response": r.raw_json,
        }
        raw_bytes = json.dumps(wrapper, indent=2, ensure_ascii=False).encode("utf-8")
        raw_path.write_bytes(raw_bytes)
        manifest_rows.append(
            {
                "batch_type": batch_type,
                "batch_number": r.batch_number,
                "number_requested": len(r.ids_requested),
                "request_url": RCSB_DATA_API_URL,
                "http_status": r.http_status if r.http_status is not None else "",
                "retrieved_at_utc": r.retrieved_at_utc,
                "raw_response_file": str(raw_path.relative_to(PROJECT_ROOT)),
                "raw_response_sha256": sha256_bytes(raw_bytes),
                "success": r.success,
                "error_message": r.error_message,
            }
        )
        status_label = "OK" if r.success else "FAILED"
        print(f"    batch {r.batch_number:3d} ({len(r.ids_requested):3d} ids) -> {status_label}"
              + (f"  [{r.error_message}]" if not r.success else ""))

    if not refresh and load_manifest_of_cached_ids(batch_type):
        # Merge with any previously-cached batches that already existed
        # before this run (covers the "some missing" partial-refresh path).
        pass  # batches are additive on disk; manifest is rebuilt fully below.

    # Rebuild the full manifest rows from ALL batch files on disk for this
    # type, so a partial fetch run still yields a complete, non-duplicated
    # manifest.
    all_rows = []
    for path in sorted(RAW_DATA_DIR.glob(f"{prefix}_batch_*.json")):
        wrapper = json.loads(path.read_text())
        raw_bytes = path.read_bytes()
        all_rows.append(
            {
                "batch_type": batch_type,
                "batch_number": wrapper["batch_number"],
                "number_requested": len(wrapper["ids_requested"]),
                "request_url": RCSB_DATA_API_URL,
                "http_status": wrapper["http_status"] if wrapper["http_status"] is not None else "",
                "retrieved_at_utc": wrapper["retrieved_at_utc"],
                "raw_response_file": str(path.relative_to(PROJECT_ROOT)),
                "raw_response_sha256": sha256_bytes(raw_bytes),
                "success": wrapper["success"],
                "error_message": wrapper["error_message"],
            }
        )
    return all_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Force re-fetching every batch from the live API.")
    args = parser.parse_args()

    if not DISCOVERY_BY_RECEPTOR_PATH.exists():
        raise RuntimeError(
            f"STOPPING: {DISCOVERY_BY_RECEPTOR_PATH} not found — run "
            f"scripts/04_discover_rcsb_candidates.py first."
        )
    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR_PATH)

    unique_entity_ids = sorted(discovery_df["polymer_entity_id"].unique().tolist())
    unique_pdb_ids = sorted(discovery_df["pdb_id"].unique().tolist())

    print(f"Unique polymer entities to enrich: {len(unique_entity_ids)}")
    print(f"Unique PDB entries to enrich: {len(unique_pdb_ids)}")

    RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("\nFetching polymer entity metadata...")
    entity_manifest_rows = run_batches_for("polymer_entity", unique_entity_ids, refresh=args.refresh)

    print("\nFetching entry metadata...")
    entry_manifest_rows = run_batches_for("entry", unique_pdb_ids, refresh=args.refresh)

    manifest_df = pd.DataFrame(entity_manifest_rows + entry_manifest_rows)[DATA_API_MANIFEST_COLUMNS]
    manifest_df = manifest_df.sort_values(["batch_type", "batch_number"]).reset_index(drop=True)
    manifest_df.to_csv(DATA_API_MANIFEST_PATH, index=False)
    print(f"\nWrote {DATA_API_MANIFEST_PATH.relative_to(PROJECT_ROOT)} ({len(manifest_df)} batch rows)")

    n_failed = int((~manifest_df["success"]).sum())
    if n_failed:
        print(f"\nWARNING: {n_failed} batch(es) failed and could not be recovered even at minimum batch size.")
        print("Stage 2 Data API enrichment is INCOMPLETE — do not treat downstream outputs as a full PASS.")
    else:
        print("\nAll Data API batches succeeded.")


if __name__ == "__main__":
    main()
