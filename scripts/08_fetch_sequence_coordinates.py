"""
Stage 2.2, step 1: query the official RCSB Sequence Coordinates API
(UNIPROT -> PDB_ENTITY) for each of the 48 validated project accessions,
to independently obtain official reference/entity sequence coverage.

Reads:
    config/nr_metadata.py
    data/manifests/rcsb_discovery_by_receptor.csv   (Stage 2 candidate set, for matching)

Writes:
    data/raw/api/rcsb/sequence_coordinates/{uniprot_id}.json   (gitignored)
    data/manifests/rcsb_sequence_coordinates_manifest.csv
    data/manifests/rcsb_sequence_coordinate_alignments.csv

One request per UniProt accession (not per polymer entity) — the API
returns all PDB_ENTITY targets for a single UNIPROT query. Computed
structure model targets (e.g. AlphaFold "AF_..." IDs) are excluded, since
Stage 2's candidate set is experimental-only; they are counted but not
retained as rows.

Idempotent by default: a cached raw response is reused unless --refresh.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402
from scripts.utils.rcsb_sequence_coordinates_client import (  # noqa: E402
    SEQUENCE_COORDINATES_URL,
    SequenceCoordinatesResult,
    fetch_alignments_for_uniprot,
)

RAW_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "sequence_coordinates"
MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
DISCOVERY_BY_RECEPTOR_PATH = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"

MANIFEST_PATH = MANIFESTS_DIR / "rcsb_sequence_coordinates_manifest.csv"
ALIGNMENTS_PATH = MANIFESTS_DIR / "rcsb_sequence_coordinate_alignments.csv"

# Standard 4-character PDB ID: digit followed by 3 alphanumerics.
PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def utc_mtime_iso(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def parse_target_id(target_id: str) -> tuple[str, str] | None:
    """Returns (pdb_id, entity_id) for an experimental PDB_ENTITY target,
    or None if `target_id` is a computed structure model (e.g. AF_/MA_)."""
    if "_" not in target_id:
        return None
    pdb_id, entity_id = target_id.rsplit("_", 1)
    if not PDB_ID_RE.match(pdb_id):
        return None
    return pdb_id, entity_id


def resolve_one(accession: str, refresh: bool) -> dict:
    raw_path = RAW_DIR / f"{accession}.json"

    if not refresh and raw_path.exists():
        try:
            wrapper = json.loads(raw_path.read_text())
            raw_bytes = raw_path.read_bytes()
            return {
                "request_url": wrapper["request"]["url"],
                "http_status": wrapper.get("http_status", 200),
                "retrieved_at_utc": utc_mtime_iso(raw_path),
                "raw_file": str(raw_path.relative_to(PROJECT_ROOT)),
                "sha256": sha256_bytes(raw_bytes),
                "success": True,
                "error_message": "",
                "raw_json": wrapper["response"],
                "_from_cache": True,
            }
        except (json.JSONDecodeError, KeyError):
            pass  # fall through to a live fetch if the cache is corrupt

    result: SequenceCoordinatesResult = fetch_alignments_for_uniprot(accession)

    if result.success:
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        wrapper = {
            "request": {
                "url": result.request_url,
                "method": "POST",
                "payload": {"from": "UNIPROT", "to": "PDB_ENTITY", "queryId": accession},
            },
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
            "raw_json": result.raw_json,
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
        "raw_json": None,
        "_from_cache": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Force re-querying every accession from the live API.")
    args = parser.parse_args()

    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR_PATH, dtype={"pdb_id": str, "entity_id": str})
    stage2_entities_by_uniprot = discovery_df.groupby("uniprot_id")["polymer_entity_id"].apply(set).to_dict()

    RAW_DIR.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    alignment_rows = []

    print(f"Querying Sequence Coordinates API (UNIPROT -> PDB_ENTITY) for {len(NR_METADATA)} receptors "
          f"via {SEQUENCE_COORDINATES_URL} ...")

    for idx, record in enumerate(NR_METADATA):
        accession = record["uniprot_id"]
        result = resolve_one(accession, refresh=args.refresh)
        source = "cache" if result.get("_from_cache") else ("live" if result["success"] else "FAILED")

        stage2_targets_expected = stage2_entities_by_uniprot.get(accession, set())

        experimental_targets: list[tuple[str, dict]] = []
        computed_model_count = 0
        query_length = None
        if result["success"] and result["raw_json"]:
            alignments_obj = result["raw_json"]["data"]["alignments"]
            query_length_seq = alignments_obj.get("query_sequence")
            query_length = len(query_length_seq) if query_length_seq else None
            for t in alignments_obj.get("target_alignments") or []:
                parsed = parse_target_id(t["target_id"])
                if parsed is None:
                    computed_model_count += 1
                    continue
                experimental_targets.append((t["target_id"], t))

        targets_returned_ids = {tid for tid, _ in experimental_targets}
        matched = targets_returned_ids & stage2_targets_expected
        # Targets the Sequence Coordinates API returns that are not in the
        # Stage 2 Search API candidate set. This does NOT imply an error or
        # an indexing-timing discrepancy — the two RCSB services apply
        # different scope filters (Stage 2's Search query restricts to
        # structure_determination_methodology == "experimental"; Sequence
        # Coordinates applies no such filter). See Stage 2.3
        # (rcsb_sequence_coordinates_extra_targets.csv) for the
        # per-target, metadata-based classification of why each one falls
        # outside Stage 2's scope.
        outside_scope = targets_returned_ids - stage2_targets_expected
        missing = stage2_targets_expected - targets_returned_ids

        print(f"  [{idx+1:2d}/48] {accession:10s} common_name={record['common_name']:10s} -> {source}  "
              f"experimental_targets={len(experimental_targets)}  matched={len(matched)}  "
              f"missing={len(missing)}  outside_stage2_scope={len(outside_scope)}  computed_models_excluded={computed_model_count}")

        manifest_rows.append(
            {
                "uniprot_id": accession,
                "nr_code": record["nr_code"],
                "common_name": record["common_name"],
                "http_status": result["http_status"],
                "retrieved_at_utc": result["retrieved_at_utc"],
                "targets_returned": len(experimental_targets) if result["success"] else "",
                "stage2_targets_expected": len(stage2_targets_expected),
                "stage2_targets_matched": len(matched) if result["success"] else "",
                "targets_outside_stage2_scope": json.dumps(sorted(outside_scope)) if result["success"] else "",
                "missing_stage2_targets": json.dumps(sorted(missing)) if result["success"] else "",
                "raw_response_file": result["raw_file"],
                "raw_response_sha256": result["sha256"],
                "success": result["success"],
                "error_message": result["error_message"],
            }
        )

        if result["success"]:
            for target_id, t in experimental_targets:
                pdb_id, entity_id = parse_target_id(target_id)
                coverage = t.get("coverage") or {}
                alignment_rows.append(
                    {
                        "uniprot_id": accession,
                        "nr_code": record["nr_code"],
                        "common_name": record["common_name"],
                        "polymer_entity_id": target_id,
                        "pdb_id": pdb_id,
                        "entity_id": entity_id,
                        "query_length": coverage.get("query_length"),
                        "target_length": coverage.get("target_length"),
                        "rcsb_reference_sequence_coverage": coverage.get("query_coverage"),
                        "rcsb_entity_sequence_coverage": coverage.get("target_coverage"),
                        "aligned_regions_json": json.dumps(t.get("aligned_regions") or [], separators=(",", ":")),
                        "is_stage2_candidate": target_id in stage2_targets_expected,
                    }
                )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(MANIFEST_PATH, index=False)
    print(f"\nWrote {MANIFEST_PATH.relative_to(PROJECT_ROOT)} ({len(manifest_df)} rows)")

    # Flag duplicate alignment objects for the same (uniprot_id, polymer_entity_id)
    # pair rather than silently collapsing them.
    alignments_df = pd.DataFrame(alignment_rows)
    if len(alignments_df):
        dup_mask = alignments_df.duplicated(subset=["uniprot_id", "polymer_entity_id"], keep=False)
        alignments_df["duplicate_alignment_flagged"] = dup_mask
    alignments_df.to_csv(ALIGNMENTS_PATH, index=False)
    print(f"Wrote {ALIGNMENTS_PATH.relative_to(PROJECT_ROOT)} ({len(alignments_df)} rows)")

    n_success = int(manifest_df["success"].sum())
    n_failed = len(manifest_df) - n_success
    print(f"\n--- Sequence Coordinates summary ---")
    print(f"queried: {len(manifest_df)}  successful: {n_success}  failed: {n_failed}")
    if n_success:
        print(f"total experimental target alignments: {len(alignments_df)}")
        print(f"Stage 2 candidates matched: {int(manifest_df['stage2_targets_matched'].fillna(0).sum())}")


if __name__ == "__main__":
    main()
