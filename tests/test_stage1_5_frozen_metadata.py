"""
Stage 1.5 tests: the frozen, UniProt-validated working receptor master.

Operate only on already-generated files (config/nr_metadata.py,
data/manifests/nr_metadata_validated.csv, the Stage 1 audit, and the
immutable professor snapshot). No network access, no raw UniProt responses
required.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
VALIDATED_CSV = MANIFESTS_DIR / "nr_metadata_validated.csv"
VALIDATED_PROVENANCE_JSON = MANIFESTS_DIR / "nr_metadata_validated_provenance.json"
AUDIT_CSV = MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv"
PROFESSOR_SNAPSHOT_PATH = PROJECT_ROOT / "config" / "nr_metadata_professor.py"


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet — run scripts/03_freeze_validated_metadata.py first")


def test_working_metadata_has_exactly_48_records():
    from config.nr_metadata import NR_METADATA

    assert len(NR_METADATA) == 48


def test_working_metadata_uniprot_ids_unique():
    from config.nr_metadata import NR_METADATA

    ids = [r["uniprot_id"] for r in NR_METADATA]
    assert len(set(ids)) == 48


def test_working_metadata_nr_codes_unique():
    from config.nr_metadata import NR_METADATA

    codes = [r["nr_code"] for r in NR_METADATA]
    assert len(set(codes)) == 48


def test_working_records_match_professor_for_non_accession_fields():
    from config.nr_metadata import NR_METADATA as WORKING
    from config.nr_metadata_professor import NR_METADATA as PROFESSOR

    assert len(WORKING) == len(PROFESSOR)
    for working, professor in zip(WORKING, PROFESSOR):
        assert working["nr_code"] == professor["nr_code"]
        assert working["common_name"] == professor["common_name"]
        assert working["search_terms"] == professor["search_terms"]
        assert working["group"] == professor["group"]


def test_working_accession_equals_stage1_resolved_primary_accession():
    _require(AUDIT_CSV)
    from config.nr_metadata import NR_METADATA as WORKING

    audit_df = pd.read_csv(AUDIT_CSV).set_index("record_index")
    for idx, record in enumerate(WORKING):
        resolved = audit_df.loc[idx, "resolved_primary_accession"]
        assert record["uniprot_id"] == resolved


def test_validated_csv_has_exactly_48_rows():
    _require(VALIDATED_CSV)
    df = pd.read_csv(VALIDATED_CSV)
    assert len(df) == 48


def test_validated_csv_every_status_is_pass():
    _require(VALIDATED_CSV)
    df = pd.read_csv(VALIDATED_CSV)
    assert (df["validation_status"] == "PASS").all()


def test_validated_csv_search_terms_are_valid_json_matching_professor():
    _require(VALIDATED_CSV)
    from config.nr_metadata_professor import NR_METADATA as PROFESSOR

    df = pd.read_csv(VALIDATED_CSV).sort_values("record_index").reset_index(drop=True)
    for idx, professor in enumerate(PROFESSOR):
        parsed = json.loads(df.loc[idx, "search_terms"])
        assert parsed == professor["search_terms"]


def test_professor_immutable_file_unchanged_by_freeze():
    """The freeze must never rewrite the immutable professor snapshot."""
    _require(VALIDATED_PROVENANCE_JSON)
    import hashlib

    provenance = json.loads(VALIDATED_PROVENANCE_JSON.read_text())
    recorded_sha256 = provenance["sha256"]["config/nr_metadata_professor.py"]
    actual_sha256 = hashlib.sha256(PROFESSOR_SNAPSHOT_PATH.read_bytes()).hexdigest()
    assert actual_sha256 == recorded_sha256, (
        "config/nr_metadata_professor.py has changed since the Stage 1.5 freeze "
        "recorded its checksum — the immutable snapshot must never be edited"
    )


def test_stage1_audit_unchanged_by_freeze():
    """The freeze reads the Stage 1 audit but must never modify it."""
    _require(VALIDATED_PROVENANCE_JSON)
    import hashlib

    provenance = json.loads(VALIDATED_PROVENANCE_JSON.read_text())
    recorded_sha256 = provenance["sha256"]["data/manifests/nr_metadata_uniprot_audit.csv"]
    actual_sha256 = hashlib.sha256(AUDIT_CSV.read_bytes()).hexdigest()
    assert actual_sha256 == recorded_sha256, (
        "data/manifests/nr_metadata_uniprot_audit.csv has changed since the "
        "Stage 1.5 freeze recorded its checksum"
    )


def test_validated_provenance_counts_are_consistent():
    _require(VALIDATED_PROVENANCE_JSON)
    provenance = json.loads(VALIDATED_PROVENANCE_JSON.read_text())
    assert provenance["professor_source_record_count"] == 48
    assert provenance["validated_record_count"] == 48
    assert provenance["pass_count"] == 48
    assert provenance["review_needed_count"] == 0
    assert provenance["fail_count"] == 0
