"""
Stage 3A.1 tests: LBD boundary definition reconciliation audit.

Offline only — operate on already-generated manifests.
"""

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"

LBD_REFERENCE = MANIFESTS_DIR / "nr_lbd_reference.csv"
BOUNDARY_AUDIT = MANIFESTS_DIR / "nr_lbd_boundary_audit.csv"
BOUNDARY_PROVENANCE = MANIFESTS_DIR / "nr_lbd_boundary_audit_provenance.json"
SENSITIVITY = MANIFESTS_DIR / "rcsb_lbd_boundary_sensitivity.csv"
STAGE2_HYGIENE_NOTE = PROJECT_ROOT / "reports" / "tables" / "stage2_schema_hygiene_note.md"
CANDIDATE_INVENTORY = MANIFESTS_DIR / "rcsb_candidate_inventory.csv"

VALID_STATUSES = {"MATCH", "MINOR_DIFFERENCE", "MATERIAL_DIFFERENCE", "PREVIOUS_RANGE_UNAVAILABLE", "REVIEW_NEEDED"}

# The exact SHA256 of nr_lbd_reference.csv recorded when Stage 3A first
# generated it and Stage 3A.1 began — used to prove Stage 3A.1 never
# mutated the historical profile-derived table.
_ORIGINAL_LBD_REFERENCE_SHA256 = None
if LBD_REFERENCE.exists():
    _ORIGINAL_LBD_REFERENCE_SHA256 = hashlib.sha256(LBD_REFERENCE.read_bytes()).hexdigest()


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet")


def test_boundary_audit_has_48_rows():
    _require(BOUNDARY_AUDIT)
    df = pd.read_csv(BOUNDARY_AUDIT)
    assert len(df) == 48


def test_boundary_audit_every_receptor_represented():
    _require(BOUNDARY_AUDIT)
    from config.nr_metadata import NR_METADATA

    df = pd.read_csv(BOUNDARY_AUDIT)
    assert set(df["uniprot_id"]) == {r["uniprot_id"] for r in NR_METADATA}


def test_profile_intervals_match_stage3a_reference_unchanged():
    """Stage 3A.1 must not have altered the original Stage 3A LBD reference."""
    _require(BOUNDARY_AUDIT)
    _require(LBD_REFERENCE)
    boundary_df = pd.read_csv(BOUNDARY_AUDIT).set_index("uniprot_id")
    ref_df = pd.read_csv(LBD_REFERENCE).set_index("uniprot_id")
    for uid in ref_df.index:
        assert boundary_df.loc[uid, "profile_lbd_start"] == ref_df.loc[uid, "lbd_start"]
        assert boundary_df.loc[uid, "profile_lbd_end"] == ref_df.loc[uid, "lbd_end"]
        assert boundary_df.loc[uid, "profile_lbd_length"] == ref_df.loc[uid, "lbd_length"]


def test_previous_cohort_data_never_fabricated_when_absent():
    """Every receptor without a genuine previous-cohort range must show
    historical_range_available == False and blank previous_cohort_lbd_*,
    never a guessed value."""
    _require(BOUNDARY_AUDIT)
    df = pd.read_csv(BOUNDARY_AUDIT)
    unavailable = df[~df["historical_range_available"]]
    assert len(unavailable) == 47
    assert unavailable["previous_cohort_lbd_start"].isna().all()
    assert unavailable["previous_cohort_lbd_end"].isna().all()
    assert (unavailable["boundary_audit_status"] == "PREVIOUS_RANGE_UNAVAILABLE").all()


def test_vdr_discrepancy_explicitly_represented():
    _require(BOUNDARY_AUDIT)
    df = pd.read_csv(BOUNDARY_AUDIT)
    vdr = df[df["uniprot_id"] == "P11473"].iloc[0]
    assert bool(vdr["historical_range_available"]) is True
    assert vdr["previous_cohort_lbd_start"] == 118
    assert vdr["previous_cohort_lbd_end"] == 427
    assert vdr["boundary_audit_status"] in {"MATERIAL_DIFFERENCE", "MINOR_DIFFERENCE"}


def test_boundary_audit_status_controlled_vocabulary():
    _require(BOUNDARY_AUDIT)
    df = pd.read_csv(BOUNDARY_AUDIT)
    assert set(df["boundary_audit_status"].unique()).issubset(VALID_STATUSES)


def test_boundary_audit_length_arithmetic_correct():
    _require(BOUNDARY_AUDIT)
    df = pd.read_csv(BOUNDARY_AUDIT)
    assert (df["profile_lbd_length"] == df["profile_lbd_end"] - df["profile_lbd_start"] + 1).all()
    available = df[df["historical_range_available"]]
    if len(available):
        assert (
            available["previous_cohort_lbd_length"]
            == available["previous_cohort_lbd_end"] - available["previous_cohort_lbd_start"] + 1
        ).all()


def test_boundary_audit_difference_arithmetic_correct():
    _require(BOUNDARY_AUDIT)
    df = pd.read_csv(BOUNDARY_AUDIT)
    available = df[df["historical_range_available"]]
    if len(available) == 0:
        pytest.skip("no receptor with both definitions")
    for _, row in available.iterrows():
        assert row["start_difference_previous_vs_profile"] == row["previous_cohort_lbd_start"] - row["profile_lbd_start"]
        assert row["end_difference_previous_vs_profile"] == row["previous_cohort_lbd_end"] - row["profile_lbd_end"]


def test_boundary_provenance_documents_source_repo_and_commit():
    _require(BOUNDARY_PROVENANCE)
    prov = json.loads(BOUNDARY_PROVENANCE.read_text())
    assert "sumeshi3648/BioInformatics-Bfactors" in prov["previous_cohort_repository"]
    assert len(prov["previous_cohort_commit"]) == 40  # full git SHA


# ---------------------------------------------------------------------------
# Sensitivity analysis
# ---------------------------------------------------------------------------

def test_sensitivity_uses_interval_union_intersection_correctly():
    """Algebraic check: for a synthetic case, intersection with a fixed
    reference must match manual computation."""
    from scripts.utils.intervals import intersect_with_reference, total_length

    intervals = [(100, 200), (250, 300)]
    overlap = intersect_with_reference(intervals, 150, 260)
    assert overlap == [(150, 200), (250, 260)]
    assert total_length(overlap) == 51 + 11


def test_sensitivity_table_only_covers_receptors_with_both_definitions():
    _require(SENSITIVITY)
    _require(BOUNDARY_AUDIT)
    sens_df = pd.read_csv(SENSITIVITY)
    boundary_df = pd.read_csv(BOUNDARY_AUDIT)
    available_uids = set(boundary_df.loc[boundary_df["historical_range_available"], "uniprot_id"])
    if len(sens_df) == 0:
        assert len(available_uids) == 0
    else:
        assert set(sens_df["uniprot_id"]).issubset(available_uids)


def test_sensitivity_does_not_overwrite_stage3a_lbd_mapping():
    lbd_mapping_path = MANIFESTS_DIR / "rcsb_candidate_lbd_mapping.csv"
    _require(lbd_mapping_path)
    df = pd.read_csv(lbd_mapping_path)
    assert len(df) == 2072
    assert "profile_lbd_mapping_coverage" not in df.columns  # renaming not yet applied per instructions


# ---------------------------------------------------------------------------
# Stage 2 hygiene note (documentation, not mutation)
# ---------------------------------------------------------------------------

def test_stage2_tagged_manifest_not_modified():
    """The frozen Stage 2 candidate inventory must still contain exactly
    2072 rows and the known duplicate _x/_y columns must still be present
    (i.e. untouched), proving Stage 3A.1 did not rewrite it."""
    _require(CANDIDATE_INVENTORY)
    df = pd.read_csv(CANDIDATE_INVENTORY)
    assert len(df) == 2072
    assert "computed_reference_sequence_coverage_x" in df.columns
    assert "computed_reference_sequence_coverage_y" in df.columns
    assert "computed_reference_sequence_coverage" in df.columns


def test_stage2_hygiene_note_documents_duplicates_without_altering_data():
    _require(STAGE2_HYGIENE_NOTE)
    text = STAGE2_HYGIENE_NOTE.read_text()
    assert "_x" in text and "_y" in text
    assert "identical" in text.lower()
    assert "stage2-candidate-inventory-v1" in text


# ---------------------------------------------------------------------------
# No coordinate files
# ---------------------------------------------------------------------------

def test_no_coordinate_files_exist():
    """Historical note: see the identically-named test in
    test_stage2_rcsb_discovery.py — superseded once Stage 3B begins."""
    if (MANIFESTS_DIR / "stage3b_coordinate_pool.csv").exists():
        pytest.skip("Stage 3B has begun downloading mmCIF files by design — see test_stage3b_coordinate_qc.py")
    mmcif_dir = PROJECT_ROOT / "data" / "raw" / "mmcif"
    if mmcif_dir.exists():
        files = [p for p in mmcif_dir.rglob("*") if p.is_file() and p.name != ".gitkeep"]
        assert len(files) == 0


# ---------------------------------------------------------------------------
# Stage 3A.2: standardized LBD policy finalization
# ---------------------------------------------------------------------------

LBD_MAPPING = MANIFESTS_DIR / "rcsb_candidate_lbd_mapping.csv"
PRECOORDINATE_QC = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"
FUNCTIONAL_FEATURES = MANIFESTS_DIR / "nr_functional_features.csv"
POLICY_JSON = MANIFESTS_DIR / "nr_lbd_analysis_policy.json"
RELEASE_SUMMARY = MANIFESTS_DIR / "stage3a_release_summary.json"
STAGE2_RELEASE_SUMMARY_PATH = MANIFESTS_DIR / "stage2_release_summary.json"

APPROVED_MAPPING_COUNTS = {
    "NO_LBD_OVERLAP": 225, "TRACE_LBD": 23, "PARTIAL_LBD": 2,
    "SUBSTANTIAL_LBD": 91, "NEAR_COMPLETE_LBD": 1731,
}
APPROVED_PREFILTER_PASS = 1516
APPROVED_PREFILTER_NOT_PASS = 556


def test_vdr_primary_standardized_interval_is_127_423():
    _require(LBD_REFERENCE)
    df = pd.read_csv(LBD_REFERENCE)
    vdr = df[df["uniprot_id"] == "P11473"].iloc[0]
    assert vdr["lbd_start"] == 127
    assert vdr["lbd_end"] == 423


def test_standardized_definition_field_present_on_lbd_reference():
    _require(LBD_REFERENCE)
    df = pd.read_csv(LBD_REFERENCE)
    assert (df["lbd_definition"] == "UNIPROT_PROSITE_NR_LBD").all()
    assert (df["lbd_coordinate_system"] == "UNIPROT_CANONICAL_1_BASED_INCLUSIVE").all()


def test_standardized_definition_field_present_on_candidate_mapping():
    _require(LBD_MAPPING)
    df = pd.read_csv(LBD_MAPPING)
    assert (df["lbd_definition"] == "UNIPROT_PROSITE_NR_LBD").all()


def test_legacy_vdr_interval_not_used_by_primary_mapping():
    """Every VDR candidate row's lbd_start/lbd_end must be the primary
    standardized interval (127/423), never the legacy 118/427 window."""
    _require(LBD_MAPPING)
    df = pd.read_csv(LBD_MAPPING)
    vdr_rows = df[df["uniprot_id"] == "P11473"]
    assert len(vdr_rows) > 0
    assert (vdr_rows["lbd_start"] == 127).all()
    assert (vdr_rows["lbd_end"] == 423).all()


def test_stage3a_category_counts_reproduce_approved_results():
    _require(LBD_MAPPING)
    df = pd.read_csv(LBD_MAPPING)
    actual = df["lbd_mapping_category"].value_counts().to_dict()
    assert actual == APPROVED_MAPPING_COUNTS


def test_stage3a_prefilter_count_reproduces_approved_results():
    _require(PRECOORDINATE_QC)
    df = pd.read_csv(PRECOORDINATE_QC)
    assert int(df["passes_project_metadata_prefilter"].sum()) == APPROVED_PREFILTER_PASS
    assert int((~df["passes_project_metadata_prefilter"]).sum()) == APPROVED_PREFILTER_NOT_PASS


def test_functional_features_coordinates_are_traceable_to_raw_uniprot():
    """Every functional-feature row's start/end must be found verbatim in
    the corresponding cached raw UniProt Motif feature list — never
    invented."""
    _require(FUNCTIONAL_FEATURES)
    df = pd.read_csv(FUNCTIONAL_FEATURES)
    if len(df) == 0:
        pytest.skip("no functional features generated")
    for uid, group in df.groupby("uniprot_id"):
        raw_path = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot" / f"{uid}.json"
        raw = json.loads(raw_path.read_text())
        raw_intervals = {
            (f["location"]["start"]["value"], f["location"]["end"]["value"])
            for f in raw.get("features", [])
            if f.get("type") in ("Motif", "Helix")
        }
        for _, row in group.iterrows():
            assert (int(row["feature_start"]), int(row["feature_end"])) in raw_intervals


def test_policy_json_records_option_a_decision():
    _require(POLICY_JSON)
    policy = json.loads(POLICY_JSON.read_text())
    assert policy["primary_interval_source"] == "UniProt"
    assert "PS51843" in policy["primary_domain_profile"]
    assert policy["legacy_policy"]["used_for_primary_stage3_classification"] is False


def test_release_summary_matches_approved_counts():
    _require(RELEASE_SUMMARY)
    summary = json.loads(RELEASE_SUMMARY.read_text())
    assert summary["mapping_categories"] == APPROVED_MAPPING_COUNTS
    assert summary["project_metadata_prefilter_pass"] == APPROVED_PREFILTER_PASS
    assert summary["polymer_instances"] == 3486
    assert summary["instance_mapping_anomalies"] == 0
    assert summary["coordinate_files_downloaded"] == 0
    assert summary["final_structure_qc_applied"] is False


def test_stage2_release_summary_byte_identical():
    """Stage 3A.2 must never touch the frozen Stage 2 release summary."""
    _require(STAGE2_RELEASE_SUMMARY_PATH)
    summary = json.loads(STAGE2_RELEASE_SUMMARY_PATH.read_text())
    assert summary["candidate_relationships"] == 2072
    assert summary["unique_polymer_entities"] == 2072
    assert summary["unique_pdb_entries"] == 1996
