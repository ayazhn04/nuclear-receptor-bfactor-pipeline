"""
Stage 1.5: freeze the Stage 1 UniProt-validated professor metadata into the
project's working receptor master.

Reads:
    config/nr_metadata_professor.py           (immutable professor snapshot)
    data/manifests/nr_metadata_uniprot_audit.csv   (Stage 1 audit, must all be PASS)

Writes:
    config/nr_metadata.py                          (generated, not hand-typed)
    data/manifests/nr_metadata_validated.csv
    data/manifests/nr_metadata_validated_provenance.json

This script performs NO network access. It fails loudly (raises, does not
print a false PASS) if:
  - the professor snapshot and the Stage 1 audit disagree in row count,
  - the professor snapshot and audit disagree in per-record identity
    (uniprot_id / nr_code / common_name / group),
  - any audited receptor is not audit_status == PASS,
  - any of the specific Stage-1 assertions listed in the Stage 1.5 spec fail.

It never edits config/nr_metadata_professor.py or the raw UniProt responses.
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

from config.nr_metadata_professor import NR_METADATA as PROFESSOR_NR_METADATA  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
CONFIG_DIR = PROJECT_ROOT / "config"

AUDIT_CSV_PATH = MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv"
PROFESSOR_SNAPSHOT_PATH = CONFIG_DIR / "nr_metadata_professor.py"
PROFESSOR_CSV_PATH = MANIFESTS_DIR / "nr_metadata_professor.csv"
WORKING_METADATA_PATH = CONFIG_DIR / "nr_metadata.py"
VALIDATED_CSV_PATH = MANIFESTS_DIR / "nr_metadata_validated.csv"
VALIDATED_PROVENANCE_PATH = MANIFESTS_DIR / "nr_metadata_validated_provenance.json"

EXPECTED_RECORD_COUNT = 48


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_and_verify_audit() -> pd.DataFrame:
    if not AUDIT_CSV_PATH.exists():
        raise RuntimeError(f"STOPPING: Stage 1 audit not found at {AUDIT_CSV_PATH}")

    audit_df = pd.read_csv(AUDIT_CSV_PATH)

    if len(PROFESSOR_NR_METADATA) != EXPECTED_RECORD_COUNT:
        raise RuntimeError(
            f"STOPPING: professor snapshot has {len(PROFESSOR_NR_METADATA)} records, "
            f"expected {EXPECTED_RECORD_COUNT}."
        )
    if len(audit_df) != len(PROFESSOR_NR_METADATA):
        raise RuntimeError(
            f"STOPPING: audit row count ({len(audit_df)}) does not match professor "
            f"snapshot record count ({len(PROFESSOR_NR_METADATA)})."
        )

    # Row-by-row identity agreement between the immutable professor snapshot
    # and the audit table (defends against the audit having been generated
    # from a different/stale snapshot).
    audit_by_index = audit_df.set_index("record_index")
    for idx, record in enumerate(PROFESSOR_NR_METADATA):
        if idx not in audit_by_index.index:
            raise RuntimeError(f"STOPPING: audit is missing record_index {idx}")
        row = audit_by_index.loc[idx]
        mismatches = []
        if row["professor_uniprot_id"] != record["uniprot_id"]:
            mismatches.append("uniprot_id")
        if row["professor_nr_code"] != record["nr_code"]:
            mismatches.append("nr_code")
        if row["professor_common_name"] != record["common_name"]:
            mismatches.append("common_name")
        if row["professor_group"] != record["group"]:
            mismatches.append("group")
        if mismatches:
            raise RuntimeError(
                f"STOPPING: professor snapshot and audit disagree at record_index {idx} "
                f"on fields {mismatches}. Snapshot and audit are out of sync — re-run "
                f"Stage 1 before freezing."
            )

    assertions = {
        "row_count_48": len(audit_df) == EXPECTED_RECORD_COUNT,
        "unique_professor_uniprot_id_48": audit_df["professor_uniprot_id"].nunique() == EXPECTED_RECORD_COUNT,
        "unique_professor_nr_code_48": audit_df["professor_nr_code"].nunique() == EXPECTED_RECORD_COUNT,
        "all_audit_status_pass": (audit_df["audit_status"] == "PASS").all(),
        "all_taxon_9606": (audit_df["organism_taxon_id"] == 9606).all(),
        "all_reviewed_true": (audit_df["entry_is_reviewed"] == True).all(),  # noqa: E712
        "resolved_accession_nonnull": audit_df["resolved_primary_accession"].notna().all()
        and (audit_df["resolved_primary_accession"].astype(str) != "").all(),
        "resolved_equals_professor": (
            audit_df["resolved_primary_accession"] == audit_df["professor_uniprot_id"]
        ).all(),
        "no_fail": (audit_df["audit_status"] == "FAIL").sum() == 0,
        "no_review_needed": (audit_df["audit_status"] == "REVIEW_NEEDED").sum() == 0,
    }

    failed = {k: v for k, v in assertions.items() if not v}
    print("--- Stage 1 audit verification ---")
    for k, v in assertions.items():
        print(f"  {k}: {v}")
    if failed:
        raise RuntimeError(f"STOPPING: Stage 1 audit verification failed: {list(failed.keys())}")

    print("ALL STAGE 1 ASSERTIONS PASS — safe to freeze.\n")
    return audit_df


def format_python_literal(value) -> str:
    """repr() with a stable, human-readable style for the generated config file."""
    return repr(value)


def write_working_metadata(audit_df: pd.DataFrame) -> None:
    audit_by_index = audit_df.set_index("record_index")

    lines: list[str] = []
    lines.append('"""')
    lines.append("Validated working human nuclear receptor metadata master.")
    lines.append("")
    lines.append(
        "GENERATED FILE — do not hand-edit. Produced by "
        "scripts/03_freeze_validated_metadata.py from the professor-supplied "
        "snapshot (config/nr_metadata_professor.py) after Stage 1 UniProt "
        "identity validation."
    )
    lines.append("")
    lines.append("Provenance:")
    lines.append("  - Source: professor-supplied 48-receptor NR_METADATA list.")
    lines.append(
        "  - Every supplied UniProt accession was independently validated "
        "against the current official UniProt REST API "
        "(https://rest.uniprot.org/) during Stage 1."
    )
    lines.append(
        "  - All 48 accessions resolved successfully as reviewed (Swiss-Prot) "
        "Homo sapiens (taxon 9606) entries."
    )
    lines.append("  - No primary-accession changes were detected for any record.")
    lines.append(
        "  - Professor-supplied `group` labels (project-level groupings, not "
        "an authoritative taxonomy) are retained unchanged."
    )
    lines.append(
        "  - Detailed audit and provenance: "
        "data/manifests/nr_metadata_uniprot_audit.csv, "
        "data/manifests/nr_metadata_uniprot_audit_summary.json, "
        "data/manifests/nr_metadata_validated.csv, "
        "data/manifests/nr_metadata_validated_provenance.json."
    )
    lines.append(
        "  - Do not hand-edit this list. Any future correction must go "
        "through a documented, re-run validation step, not a silent edit here."
    )
    lines.append('"""')
    lines.append("")
    lines.append("NR_METADATA = [")
    for idx, record in enumerate(PROFESSOR_NR_METADATA):
        row = audit_by_index.loc[idx]
        validated_uniprot_id = row["resolved_primary_accession"]
        lines.append("    {")
        lines.append(f'        "uniprot_id": {format_python_literal(validated_uniprot_id)},')
        lines.append(f'        "nr_code": {format_python_literal(record["nr_code"])},')
        lines.append(f'        "common_name": {format_python_literal(record["common_name"])},')
        lines.append(f'        "search_terms": {format_python_literal(record["search_terms"])},')
        lines.append(f'        "group": {format_python_literal(record["group"])},')
        lines.append("    },")
    lines.append("]")
    lines.append("")

    WORKING_METADATA_PATH.write_text("\n".join(lines))


def write_validated_csv(audit_df: pd.DataFrame) -> None:
    audit_by_index = audit_df.set_index("record_index")
    generated_at_utc = datetime.now(timezone.utc).isoformat()

    rows = []
    for idx, record in enumerate(PROFESSOR_NR_METADATA):
        row = audit_by_index.loc[idx]
        rows.append(
            {
                "record_index": idx,
                "uniprot_id": row["resolved_primary_accession"],
                "nr_code": record["nr_code"],
                "common_name": record["common_name"],
                "search_terms": json.dumps(record["search_terms"], ensure_ascii=False, separators=(",", ":")),
                "group": record["group"],
                "uniprot_entry_name": row["uniprot_entry_name"],
                "primary_gene_name": row["primary_gene_name"],
                "recommended_protein_name": row["recommended_protein_name"],
                "organism_taxon_id": row["organism_taxon_id"],
                "reviewed_status": row["reviewed_status"],
                "uniprot_entry_version": row["entry_version"],
                "uniprot_sequence_version": row["sequence_version"],
                "uniprot_last_annotation_update_date": row["last_annotation_update_date"],
                "validation_status": row["audit_status"],
                "validated_at_utc": generated_at_utc,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(VALIDATED_CSV_PATH, index=False)


def write_validated_provenance(audit_df: pd.DataFrame) -> None:
    status_counts = audit_df["audit_status"].value_counts().to_dict()

    provenance = {
        "project_name": "Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families",
        "description": (
            "Provenance record for the Stage 1.5 freeze: the professor-supplied "
            "nuclear receptor metadata, after independent Stage 1 UniProt identity "
            "validation, converted into the project's working receptor master "
            "(config/nr_metadata.py) and a parallel validated CSV."
        ),
        "professor_source_record_count": len(PROFESSOR_NR_METADATA),
        "validated_record_count": len(audit_df),
        "pass_count": int(status_counts.get("PASS", 0)),
        "review_needed_count": int(status_counts.get("REVIEW_NEEDED", 0)),
        "fail_count": int(status_counts.get("FAIL", 0)),
        "validation_source": "UniProt REST API (https://rest.uniprot.org/)",
        "stage1_audit_file": str(AUDIT_CSV_PATH.relative_to(PROJECT_ROOT)),
        "professor_snapshot_file": str(PROFESSOR_SNAPSHOT_PATH.relative_to(PROJECT_ROOT)),
        "validated_master_csv_file": str(VALIDATED_CSV_PATH.relative_to(PROJECT_ROOT)),
        "working_metadata_file": str(WORKING_METADATA_PATH.relative_to(PROJECT_ROOT)),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "sha256": {
            "config/nr_metadata_professor.py": sha256_of_file(PROFESSOR_SNAPSHOT_PATH),
            "config/nr_metadata.py": sha256_of_file(WORKING_METADATA_PATH),
            "data/manifests/nr_metadata_professor.csv": sha256_of_file(PROFESSOR_CSV_PATH),
            "data/manifests/nr_metadata_uniprot_audit.csv": sha256_of_file(AUDIT_CSV_PATH),
            "data/manifests/nr_metadata_validated.csv": sha256_of_file(VALIDATED_CSV_PATH),
        },
        "python_version": sys.version,
        "generating_script": "scripts/03_freeze_validated_metadata.py",
    }
    VALIDATED_PROVENANCE_PATH.write_text(json.dumps(provenance, indent=2) + "\n")


def main() -> None:
    audit_df = load_and_verify_audit()

    write_working_metadata(audit_df)
    print(f"Wrote {WORKING_METADATA_PATH.relative_to(PROJECT_ROOT)}")

    write_validated_csv(audit_df)
    print(f"Wrote {VALIDATED_CSV_PATH.relative_to(PROJECT_ROOT)}")

    write_validated_provenance(audit_df)
    print(f"Wrote {VALIDATED_PROVENANCE_PATH.relative_to(PROJECT_ROOT)}")

    # Final self-check: re-import the freshly generated config/nr_metadata.py
    # and confirm it matches expectations before declaring success.
    import importlib
    import config.nr_metadata as working_module

    importlib.reload(working_module)
    working_records = working_module.NR_METADATA

    assert len(working_records) == EXPECTED_RECORD_COUNT, (
        f"Generated config/nr_metadata.py has {len(working_records)} records, "
        f"expected {EXPECTED_RECORD_COUNT}"
    )
    assert len({r["uniprot_id"] for r in working_records}) == EXPECTED_RECORD_COUNT
    assert len({r["nr_code"] for r in working_records}) == EXPECTED_RECORD_COUNT

    print(f"\nGenerated NR_METADATA record count: {len(working_records)}")
    print("FREEZE: PASS")


if __name__ == "__main__":
    main()
