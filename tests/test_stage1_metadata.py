"""
Stage 1 tests: professor metadata ingestion and UniProt audit.

These tests operate ONLY on already-generated manifest/audit files and the
immutable professor snapshot. They perform NO network access — if the audit
outputs don't exist yet, the tests that need them are skipped with a clear
message rather than triggering fresh API calls.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"

PROFESSOR_CSV = MANIFESTS_DIR / "nr_metadata_professor.csv"
LOCAL_VALIDATION_CSV = MANIFESTS_DIR / "nr_metadata_local_validation.csv"
AUDIT_CSV = MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv"
AUDIT_SUMMARY_JSON = MANIFESTS_DIR / "nr_metadata_uniprot_audit_summary.json"
MANIFEST_CSV = MANIFESTS_DIR / "uniprot_api_manifest.csv"

VALID_AUDIT_STATUSES = {"PASS", "REVIEW_NEEDED", "FAIL"}


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet — run the Stage 1 scripts first")


def test_professor_snapshot_has_exactly_48_records():
    from config.nr_metadata_professor import NR_METADATA

    assert len(NR_METADATA) == 48


def test_professor_snapshot_records_have_required_fields():
    from config.nr_metadata_professor import NR_METADATA

    required = {"uniprot_id", "nr_code", "common_name", "search_terms", "group"}
    for record in NR_METADATA:
        assert required.issubset(record.keys())


def test_nr_metadata_working_file_reflects_current_stage():
    """Historical note: during Stage 1 (audit only, no freeze yet) this test
    asserted `NR_METADATA == []`, since Stage 1 was explicitly forbidden from
    populating the working file. Stage 1.5 (scripts/03_freeze_validated_metadata.py)
    intentionally promotes the audited professor data into config/nr_metadata.py,
    so this test now asserts the post-freeze size instead. Detailed content
    checks for the frozen master live in tests/test_stage1_5_frozen_metadata.py.
    """
    from config.nr_metadata import NR_METADATA

    assert len(NR_METADATA) == 48, (
        "config/nr_metadata.py should hold the 48 Stage-1-validated records "
        "once Stage 1.5 has run; found a different count"
    )


def test_professor_csv_has_exactly_48_rows():
    _require(PROFESSOR_CSV)
    df = pd.read_csv(PROFESSOR_CSV)
    assert len(df) == 48


def test_professor_csv_uniprot_ids_are_unique():
    _require(PROFESSOR_CSV)
    df = pd.read_csv(PROFESSOR_CSV)
    assert df["uniprot_id"].is_unique


def test_professor_csv_nr_codes_are_unique():
    _require(PROFESSOR_CSV)
    df = pd.read_csv(PROFESSOR_CSV)
    assert df["nr_code"].is_unique


def test_local_validation_all_records_pass():
    _require(LOCAL_VALIDATION_CSV)
    df = pd.read_csv(LOCAL_VALIDATION_CSV)
    assert len(df) == 48
    assert df["record_locally_valid"].all()


def test_audit_has_exactly_one_row_per_professor_record():
    _require(AUDIT_CSV)
    df = pd.read_csv(AUDIT_CSV)
    assert len(df) == 48
    assert df["record_index"].nunique() == 48


def test_every_audit_status_is_valid():
    _require(AUDIT_CSV)
    df = pd.read_csv(AUDIT_CSV)
    assert set(df["audit_status"].unique()).issubset(VALID_AUDIT_STATUSES)


def test_non_pass_records_have_a_concrete_reason():
    _require(AUDIT_CSV)
    df = pd.read_csv(AUDIT_CSV)
    non_pass = df[df["audit_status"] != "PASS"]
    for _, row in non_pass.iterrows():
        reason = row["audit_reason"]
        assert isinstance(reason, str) and len(reason.strip()) > 10


def test_every_accession_has_response_or_documented_error():
    _require(MANIFEST_CSV)
    df = pd.read_csv(MANIFEST_CSV)
    assert len(df) == 48
    for _, row in df.iterrows():
        if row["success"]:
            assert isinstance(row["sha256"], str) and len(row["sha256"]) == 64
        else:
            assert isinstance(row["error_message"], str) and len(row["error_message"].strip()) > 0


def test_professor_values_unchanged_in_audit():
    """The audit table must carry the professor's original values verbatim,
    never replaced by whatever UniProt returned."""
    _require(AUDIT_CSV)
    from config.nr_metadata_professor import NR_METADATA

    audit_df = pd.read_csv(AUDIT_CSV).set_index("record_index")
    for idx, record in enumerate(NR_METADATA):
        row = audit_df.loc[idx]
        assert row["professor_uniprot_id"] == record["uniprot_id"]
        assert row["professor_nr_code"] == record["nr_code"]
        assert row["professor_common_name"] == record["common_name"]
        assert row["professor_group"] == record["group"]
        assert json.loads(row["professor_search_terms"]) == record["search_terms"]


def test_audit_summary_counts_are_internally_consistent():
    _require(AUDIT_SUMMARY_JSON)
    _require(AUDIT_CSV)
    summary = json.loads(AUDIT_SUMMARY_JSON.read_text())
    audit_df = pd.read_csv(AUDIT_CSV)

    assert summary["professor_record_count"] == 48
    assert summary["pass_count"] == (audit_df["audit_status"] == "PASS").sum()
    assert summary["review_needed_count"] == (audit_df["audit_status"] == "REVIEW_NEEDED").sum()
    assert summary["fail_count"] == (audit_df["audit_status"] == "FAIL").sum()
    assert (
        summary["pass_count"] + summary["review_needed_count"] + summary["fail_count"]
        == summary["professor_record_count"]
    )
