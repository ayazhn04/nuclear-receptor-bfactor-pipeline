"""
Stage 1, step 1: ingest and locally validate the professor-supplied NR
metadata snapshot, and export it to tracked, auditable manifest files.

Reads:
    config/nr_metadata_professor.py   (immutable professor snapshot)

Writes:
    data/manifests/nr_metadata_professor.csv
    data/manifests/nr_metadata_professor_provenance.json
    data/manifests/nr_metadata_local_validation.csv
    data/manifests/nr_metadata_local_validation_summary.json

This script performs NO network access. It never modifies
config/nr_metadata_professor.py or config/nr_metadata.py.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata_professor import NR_METADATA  # noqa: E402
from scripts.utils.local_validation import build_local_validation_summary, validate_record  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
PROFESSOR_SNAPSHOT_PATH = PROJECT_ROOT / "config" / "nr_metadata_professor.py"


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def export_professor_csv() -> Path:
    rows = []
    for idx, record in enumerate(NR_METADATA):
        rows.append(
            {
                "record_index": idx,
                "uniprot_id": record.get("uniprot_id"),
                "nr_code": record.get("nr_code"),
                "common_name": record.get("common_name"),
                # Deterministic JSON representation: stable key order not
                # needed here since it's a flat list, but we fix separators
                # so the same input always serializes identically.
                "search_terms": json.dumps(
                    record.get("search_terms"), ensure_ascii=False, separators=(",", ":")
                ),
                "group": record.get("group"),
            }
        )
    df = pd.DataFrame(rows)
    out_path = MANIFESTS_DIR / "nr_metadata_professor.csv"
    df.to_csv(out_path, index=False)
    return out_path


def write_provenance(csv_path: Path) -> Path:
    provenance = {
        "project_name": "Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families",
        "description": (
            "Tabular export of the professor-supplied human nuclear receptor "
            "metadata snapshot (config/nr_metadata_professor.py), produced "
            "for downstream Stage 1 validation. This is a transcription of "
            "the supplied data, not an independently verified master list."
        ),
        "source": "Professor-supplied NR_METADATA",
        "number_of_records": len(NR_METADATA),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "sha256_nr_metadata_professor_py": sha256_of_file(PROFESSOR_SNAPSHOT_PATH),
        "sha256_nr_metadata_professor_csv": sha256_of_file(csv_path),
        "python_version": sys.version,
        "generating_script": "scripts/01_ingest_professor_metadata.py",
    }
    out_path = MANIFESTS_DIR / "nr_metadata_professor_provenance.json"
    out_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return out_path


def export_local_validation() -> tuple[Path, Path]:
    rows = []
    for idx, record in enumerate(NR_METADATA):
        checks = validate_record(record)
        row = {"record_index": idx, "uniprot_id": record.get("uniprot_id"), "nr_code": record.get("nr_code")}
        row.update(checks)
        rows.append(row)
    df = pd.DataFrame(rows)
    csv_path = MANIFESTS_DIR / "nr_metadata_local_validation.csv"
    df.to_csv(csv_path, index=False)

    summary = build_local_validation_summary(NR_METADATA)
    summary_path = MANIFESTS_DIR / "nr_metadata_local_validation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    return csv_path, summary_path


def main() -> None:
    print(f"Loaded {len(NR_METADATA)} professor-supplied records from "
          f"{PROFESSOR_SNAPSHOT_PATH.relative_to(PROJECT_ROOT)}")

    summary = build_local_validation_summary(NR_METADATA)
    print("\n--- Local validation summary ---")
    print(f"total_records: {summary['total_records']} (expected 48)")
    print(f"record_count_matches_expected: {summary['record_count_matches_expected']}")
    print(f"unique_uniprot_id_count: {summary['unique_uniprot_id_count']}")
    print(f"unique_nr_code_count: {summary['unique_nr_code_count']}")
    print(f"group_counts: {summary['group_counts']}")
    print(f"duplicate_uniprot_ids: {summary['duplicate_uniprot_ids']}")
    print(f"duplicate_nr_codes: {summary['duplicate_nr_codes']}")
    print(f"duplicate_common_names: {summary['duplicate_common_names']}")
    print(f"exact_duplicate_record_indexes: {summary['exact_duplicate_record_indexes']}")
    print(f"records_missing_required_fields: {summary['records_missing_required_fields']}")
    print(f"all_records_locally_valid: {summary['all_records_locally_valid']}")

    if not summary["record_count_matches_expected"]:
        raise RuntimeError(
            f"STOPPING: expected exactly 48 professor records, found "
            f"{summary['total_records']}. The source list was NOT modified. "
            f"Report this to the user rather than proceeding."
        )
    if summary["duplicate_uniprot_ids"] or summary["duplicate_nr_codes"]:
        raise RuntimeError(
            f"STOPPING: duplicate identifiers detected — "
            f"duplicate_uniprot_ids={summary['duplicate_uniprot_ids']}, "
            f"duplicate_nr_codes={summary['duplicate_nr_codes']}. "
            f"The source list was NOT modified. Report this rather than proceeding."
        )

    csv_path = export_professor_csv()
    print(f"\nWrote {csv_path.relative_to(PROJECT_ROOT)}")

    provenance_path = write_provenance(csv_path)
    print(f"Wrote {provenance_path.relative_to(PROJECT_ROOT)}")

    local_csv_path, local_summary_path = export_local_validation()
    print(f"Wrote {local_csv_path.relative_to(PROJECT_ROOT)}")
    print(f"Wrote {local_summary_path.relative_to(PROJECT_ROOT)}")

    print("\nLOCAL INGESTION: PASS" if summary["all_records_locally_valid"] else "\nLOCAL INGESTION: COMPLETED WITH FLAGS (see summary)")


if __name__ == "__main__":
    main()
