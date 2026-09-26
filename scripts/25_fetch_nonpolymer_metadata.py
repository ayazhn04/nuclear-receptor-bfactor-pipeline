"""
Final Data Release, Step 5-6 (fetch phase): pull official RCSB Data API
metadata for every non-polymer (ligand/ion/solvent/additive) entity
belonging to a PDB entry with Stage 3B coordinate data (all 1840 coordinate candidates) — chemical identity (formula,
name, type), atom counts, and any RCSB-computed nonpolymer-entity
annotations. Introspected against the live GraphQL schema before writing
(CoreNonpolymerEntity / CoreChemComp / RcsbNonpolymerEntityContainerIdentifiers)
rather than guessing field names.

Water entities are intentionally excluded here (RCSB's own schema tracks
water separately via water_entity_ids, not non_polymer_entity_ids); water
instances are enumerated directly from the cached mmCIF atom records in the
next step instead.

Reads:
    data/manifests/final_structure_selection_audit.csv

Writes (raw, gitignored, cached/idempotent):
    data/raw/api/rcsb/data/entry_container_identifiers_batch_XXXX.json
    data/raw/api/rcsb/data/nonpolymer_entities_batch_XXXX.json

Writes (tracked):
    data/manifests/rcsb_nonpolymer_data_api_manifest.csv
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

from scripts.utils.rcsb_data_client import DEFAULT_BATCH_SIZE, RCSB_DATA_API_URL, fetch_batches  # noqa: E402

RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "data"
MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
SELECTION_AUDIT_CSV = MANIFESTS_DIR / "final_structure_selection_audit.csv"
NONPOLYMER_DATA_API_MANIFEST_CSV = MANIFESTS_DIR / "rcsb_nonpolymer_data_api_manifest.csv"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def cached_ids(prefix: str) -> set[str]:
    covered: set[str] = set()
    for path in sorted(RAW_DATA_DIR.glob(f"{prefix}_batch_*.json")):
        try:
            wrapper = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        covered.update(wrapper.get("ids_requested", []))
    return covered


def run_batches_for(batch_type: str, prefix: str, ids: list[str], refresh: bool) -> list[dict]:
    if not refresh:
        already = cached_ids(prefix)
        missing = [i for i in ids if i not in already]
        if not missing:
            print(f"  {batch_type}: all {len(ids)} IDs already cached — reusing.")
            return _manifest_rows_from_disk(batch_type, prefix)
        print(f"  {batch_type}: {len(missing)} of {len(ids)} IDs not yet cached — fetching those.")
        ids_to_fetch = missing
    else:
        ids_to_fetch = ids
        for path in RAW_DATA_DIR.glob(f"{prefix}_batch_*.json"):
            path.unlink()

    RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)
    # Continue numbering after any existing cached batches for this prefix.
    existing_batches = sorted(RAW_DATA_DIR.glob(f"{prefix}_batch_*.json"))
    start_number = len(existing_batches)

    results = fetch_batches(batch_type, ids_to_fetch, batch_size=DEFAULT_BATCH_SIZE)
    for r in results:
        batch_number = start_number + r.batch_number
        raw_path = RAW_DATA_DIR / f"{prefix}_batch_{batch_number:04d}.json"
        wrapper = {
            "batch_type": batch_type, "batch_number": batch_number, "ids_requested": r.ids_requested,
            "retrieved_at_utc": r.retrieved_at_utc, "http_status": r.http_status,
            "success": r.success, "error_message": r.error_message, "response": r.raw_json,
        }
        raw_path.write_bytes(json.dumps(wrapper, indent=2, ensure_ascii=False).encode("utf-8"))
        status_label = "OK" if r.success else "FAILED"
        print(f"    batch {batch_number:4d} ({len(r.ids_requested):4d} ids) -> {status_label}"
              + (f"  [{r.error_message}]" if not r.success else ""))
        if not r.success:
            raise RuntimeError(f"STOPPING: {batch_type} batch {batch_number} failed: {r.error_message}")

    return _manifest_rows_from_disk(batch_type, prefix)


def _manifest_rows_from_disk(batch_type: str, prefix: str) -> list[dict]:
    rows = []
    for path in sorted(RAW_DATA_DIR.glob(f"{prefix}_batch_*.json")):
        wrapper = json.loads(path.read_text())
        raw_bytes = path.read_bytes()
        rows.append(
            {
                "batch_type": batch_type, "batch_number": wrapper["batch_number"],
                "number_requested": len(wrapper["ids_requested"]), "request_url": RCSB_DATA_API_URL,
                "http_status": wrapper["http_status"] if wrapper["http_status"] is not None else "",
                "retrieved_at_utc": wrapper["retrieved_at_utc"],
                "raw_response_file": str(path.relative_to(PROJECT_ROOT)),
                "raw_response_sha256": sha256_bytes(raw_bytes),
                "success": wrapper["success"], "error_message": wrapper["error_message"],
            }
        )
    return rows


def load_nonpolymer_entity_ids_from_cache(pdb_ids: list[str]) -> dict[str, list[str]]:
    """Read back cached entry_container_identifiers batches and return
    {pdb_id: [non_polymer_entity_id_full, ...]} for the requested PDB IDs."""
    by_pdb: dict[str, list[str]] = {}
    for path in sorted(RAW_DATA_DIR.glob("entry_container_identifiers_batch_*.json")):
        wrapper = json.loads(path.read_text())
        response = wrapper.get("response") or {}
        for entry in (response.get("data", {}) or {}).get("entries", []) or []:
            pdb_id = entry["rcsb_id"]
            ids = entry.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids") or []
            by_pdb[pdb_id] = [f"{pdb_id}_{eid}" for eid in ids]
    return {pdb_id: by_pdb.get(pdb_id, []) for pdb_id in pdb_ids}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()

    audit = pd.read_csv(SELECTION_AUDIT_CSV, dtype={"pdb_id": str})
    # Target population: ALL coordinate-assessable candidates (QC is metadata, not a filter).
    primary_pdb_ids = sorted(audit.loc[audit["stage3b_coordinate_selected"], "pdb_id"].unique().tolist())
    print(f"Coordinate-candidate PDB entries: {len(primary_pdb_ids)}")

    print("\nFetching entry container identifiers (to enumerate non-polymer entity IDs)...")
    id_manifest_rows = run_batches_for("entry_container_identifiers", "entry_container_identifiers", primary_pdb_ids, args.refresh)

    nonpolymer_entity_ids_by_pdb = load_nonpolymer_entity_ids_from_cache(primary_pdb_ids)
    all_nonpolymer_entity_ids = sorted({eid for ids in nonpolymer_entity_ids_by_pdb.values() for eid in ids})
    print(f"Non-polymer entity IDs to fetch: {len(all_nonpolymer_entity_ids)}")
    print(f"PDB entries with zero non-polymer entities: {sum(1 for v in nonpolymer_entity_ids_by_pdb.values() if not v)}")

    print("\nFetching non-polymer entity metadata (chem_comp identity, formula, annotations)...")
    if all_nonpolymer_entity_ids:
        entity_manifest_rows = run_batches_for("nonpolymer_entity", "nonpolymer_entities", all_nonpolymer_entity_ids, args.refresh)
    else:
        entity_manifest_rows = []

    manifest_df = pd.DataFrame(id_manifest_rows + entity_manifest_rows)
    manifest_df.to_csv(NONPOLYMER_DATA_API_MANIFEST_CSV, index=False)
    print(f"\nWrote {NONPOLYMER_DATA_API_MANIFEST_CSV.relative_to(PROJECT_ROOT)} ({len(manifest_df)} rows)")


if __name__ == "__main__":
    main()
