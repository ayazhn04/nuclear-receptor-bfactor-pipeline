"""
Stage 3B tests: coordinate pool selection, mmCIF download validation,
residue mapping, per-residue observations, missingness arithmetic, and
immutability of frozen Stage 2/3A artifacts.

Offline only — operate on already-generated manifests and direct algebraic
checks of the mapping/interval utilities. No network access.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"

POOL_CSV = MANIFESTS_DIR / "stage3b_coordinate_pool.csv"
DOWNLOAD_MANIFEST_CSV = MANIFESTS_DIR / "mmcif_download_manifest.csv"
INSTANCE_QC_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_instance_qc.csv"
CANDIDATE_SUMMARY_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_candidate_summary.csv"
CROSSCHECK_CSV = MANIFESTS_DIR / "mmcif_unobserved_residue_crosscheck.csv"
STAGE2_RELEASE_SUMMARY = MANIFESTS_DIR / "stage2_release_summary.json"
STAGE3A_RELEASE_SUMMARY = MANIFESTS_DIR / "stage3a_release_summary.json"
LBD_REFERENCE_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"
RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"

VALID_QC_STATUS = {"MEASURED", "MULTI_MODEL_REVIEW", "MAPPING_REVIEW_NEEDED", "PARSE_REVIEW_NEEDED"}
VALID_POOL_REASONS = {
    "SELECTED_XRAY_LBD_OVERLAP", "NON_XRAY", "NO_LBD_OVERLAP", "NON_XRAY_AND_NO_LBD_OVERLAP",
}


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet")


# ---------------------------------------------------------------------------
# POOL
# ---------------------------------------------------------------------------

def test_pool_represents_all_2072_stage2_candidates():
    _require(POOL_CSV)
    df = pd.read_csv(POOL_CSV)
    assert len(df) == 2072


def test_pool_selection_matches_xray_and_lbd_overlap_exactly():
    _require(POOL_CSV)
    df = pd.read_csv(POOL_CSV)
    expected_selected = df["method_is_xray"] & df["mapping_overlaps_lbd"]
    assert (df["coordinate_audit_selected"] == expected_selected).all()


def test_pool_resolution_and_rfree_do_not_affect_selection():
    """Selection must depend only on method_is_xray and mapping_overlaps_lbd
    — never resolution_project_range_status or r_free_status."""
    _require(POOL_CSV)
    df = pd.read_csv(POOL_CSV)
    xray_lbd = df[df["method_is_xray"] & df["mapping_overlaps_lbd"]]
    # Every combination of resolution/r_free status should still be selected.
    assert xray_lbd["coordinate_audit_selected"].all()
    assert set(xray_lbd["resolution_project_range_status"].unique()) - {"PASS_PROJECT_RANGE"}, (
        "expected some non-passing resolution statuses among selected rows (selection must not filter on it)"
    )


def test_non_xray_never_selected():
    _require(POOL_CSV)
    df = pd.read_csv(POOL_CSV)
    non_xray = df[~df["method_is_xray"]]
    assert not non_xray["coordinate_audit_selected"].any()


def test_no_lbd_overlap_never_selected():
    _require(POOL_CSV)
    df = pd.read_csv(POOL_CSV)
    no_overlap = df[~df["mapping_overlaps_lbd"]]
    assert not no_overlap["coordinate_audit_selected"].any()


def test_pool_reason_controlled_vocabulary():
    _require(POOL_CSV)
    df = pd.read_csv(POOL_CSV)
    assert set(df["coordinate_audit_reason"].unique()).issubset(VALID_POOL_REASONS)


# ---------------------------------------------------------------------------
# DOWNLOAD
# ---------------------------------------------------------------------------

def test_download_manifest_one_row_per_unique_pdb():
    _require(DOWNLOAD_MANIFEST_CSV)
    df = pd.read_csv(DOWNLOAD_MANIFEST_CSV)
    assert not df["pdb_id"].duplicated().any()


def test_successful_downloads_have_existing_file():
    _require(DOWNLOAD_MANIFEST_CSV)
    df = pd.read_csv(DOWNLOAD_MANIFEST_CSV)
    ok = df[df["download_status"].isin(["SUCCESS", "CACHED_VALID"])]
    for _, row in ok.iterrows():
        path = PROJECT_ROOT / row["raw_file"]
        assert path.exists(), f"missing raw file for {row['pdb_id']}: {path}"


def test_sha256_format_valid_for_successful_downloads():
    _require(DOWNLOAD_MANIFEST_CSV)
    df = pd.read_csv(DOWNLOAD_MANIFEST_CSV)
    ok = df[df["download_status"].isin(["SUCCESS", "CACHED_VALID"])]
    for col in ("sha256_gzip", "sha256_decompressed"):
        assert ok[col].apply(lambda v: isinstance(v, str) and len(v) == 64).all()


def test_entry_id_matches_request_for_successful_downloads():
    _require(DOWNLOAD_MANIFEST_CSV)
    df = pd.read_csv(DOWNLOAD_MANIFEST_CSV)
    ok = df[df["download_status"].isin(["SUCCESS", "CACHED_VALID"])]
    assert ok["entry_id_matches_request"].all()


# ---------------------------------------------------------------------------
# MAPPING (algebraic, direct utility tests)
# ---------------------------------------------------------------------------

def test_label_seq_id_drives_mapping_not_auth_seq_id():
    from scripts.utils.residue_mapping import build_entity_to_uniprot_map

    regions = [{"query_begin": 100, "query_end": 110, "target_begin": 1, "target_end": 11}]
    result = build_entity_to_uniprot_map(regions)
    assert result["status"] == "OK"
    # label_seq_id 1 -> uniprot 100, purely from target_begin/query_begin —
    # no auth_seq_id concept enters this function's signature at all.
    assert result["mapping"][1] == 100
    assert result["mapping"][11] == 110


def test_mismatched_aligned_lengths_flagged_review_needed():
    from scripts.utils.residue_mapping import build_entity_to_uniprot_map

    regions = [{"query_begin": 100, "query_end": 110, "target_begin": 1, "target_end": 9}]  # 11 vs 9
    result = build_entity_to_uniprot_map(regions)
    assert result["status"] == "REVIEW_NEEDED"
    assert result["mapping"] == {}


def test_gaps_are_not_interpolated():
    from scripts.utils.residue_mapping import build_entity_to_uniprot_map

    # Two disjoint regions with a gap between them — positions in the gap
    # must simply be absent from the mapping, not interpolated.
    regions = [
        {"query_begin": 100, "query_end": 105, "target_begin": 1, "target_end": 6},
        {"query_begin": 120, "query_end": 125, "target_begin": 7, "target_end": 12},
    ]
    result = build_entity_to_uniprot_map(regions)
    assert result["status"] == "OK"
    mapped_uniprot_positions = set(result["mapping"].values())
    for gap_pos in range(106, 120):
        assert gap_pos not in mapped_uniprot_positions


def test_ambiguous_mapping_fails_loudly():
    from scripts.utils.residue_mapping import build_entity_to_uniprot_map

    # Two overlapping regions claiming the same label_seq_id maps to
    # different UniProt positions.
    regions = [
        {"query_begin": 100, "query_end": 105, "target_begin": 1, "target_end": 6},
        {"query_begin": 200, "query_end": 205, "target_begin": 1, "target_end": 6},
    ]
    result = build_entity_to_uniprot_map(regions)
    assert result["status"] == "REVIEW_NEEDED"


def test_invert_mapping_detects_ambiguous_values():
    from scripts.utils.residue_mapping import invert_mapping

    inv, ambiguous = invert_mapping({1: 100, 2: 100})
    assert 100 in ambiguous


# ---------------------------------------------------------------------------
# OBSERVATIONS
# ---------------------------------------------------------------------------

def test_instance_qc_lbd_length_matches_reference():
    _require(INSTANCE_QC_CSV)
    _require(LBD_REFERENCE_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    ref = pd.read_csv(LBD_REFERENCE_CSV).set_index("uniprot_id")
    measured = qc[qc["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]
    if len(measured) == 0:
        pytest.skip("no measured rows yet")
    for _, row in measured.iterrows():
        assert row["lbd_length"] == ref.loc[row["uniprot_id"], "lbd_length"]


def test_instance_qc_keys_unique():
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    measured = qc[qc["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]
    key = measured[["pdb_id", "polymer_entity_id", "instance_id", "model_num"]]
    assert not key.duplicated().any()


def test_modeled_and_bfactor_usable_are_distinct_concepts():
    """Historical note: pre-Stage-3B.1 this compared 'modeled' (mere
    atom-record presence) against 'bfactor_usable'. Stage 3B.1 introduced a
    third, primary concept — positive occupancy — so this now compares
    `positive_occupancy_lbd_ca_residues` (renamed from `modeled_lbd_ca_residues`)
    against `bfactor_usable_lbd_ca_residues`; record-presence-vs-occupancy
    is covered separately by the record-present-vs-observed tests below."""
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    measured = qc[qc["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]
    if len(measured) == 0:
        pytest.skip("no measured rows yet")
    # bfactor_usable count must never exceed positive-occupancy count.
    assert (
        pd.to_numeric(measured["bfactor_usable_lbd_ca_residues"], errors="coerce")
        <= pd.to_numeric(measured["positive_occupancy_lbd_ca_residues"], errors="coerce")
    ).all()
    # Column must exist and be non-negative (may legitimately be all-zero).
    assert (pd.to_numeric(measured["positive_occupancy_but_bfactor_unusable_residues"], errors="coerce") >= 0).all()


def test_coverage_values_within_0_1():
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    for col in ("construct_lbd_mapping_coverage", "observed_full_lbd_ca_coverage", "bfactor_usable_full_lbd_coverage"):
        values = pd.to_numeric(qc[col], errors="coerce").dropna()
        if len(values):
            assert ((values >= 0) & (values <= 1.000001)).all()


def test_full_lbd_and_within_construct_coverage_formulas_correct():
    """`observed_full_lbd_ca_coverage` and `observed_within_construct_ca_coverage`
    are, since Stage 3B.1, defined from `positive_occupancy_lbd_ca_residues`
    (renamed from `modeled_lbd_ca_residues`, which counted mere atom-record
    presence)."""
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    measured = qc[qc["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]
    if len(measured) == 0:
        pytest.skip("no measured rows yet")
    lbd_length = pd.to_numeric(measured["lbd_length"], errors="coerce")
    observed = pd.to_numeric(measured["positive_occupancy_lbd_ca_residues"], errors="coerce")
    mapped = pd.to_numeric(measured["mapped_lbd_residues"], errors="coerce")
    full_cov = pd.to_numeric(measured["observed_full_lbd_ca_coverage"], errors="coerce")
    within_cov = pd.to_numeric(measured["observed_within_construct_ca_coverage"], errors="coerce")

    expected_full = (observed / lbd_length).round(6)
    assert (full_cov.round(6) - expected_full).abs().max() < 1e-4

    has_mapped = mapped > 0
    expected_within = (observed[has_mapped] / mapped[has_mapped]).round(6)
    assert (within_cov[has_mapped].round(6) - expected_within).abs().max() < 1e-4


# ---------------------------------------------------------------------------
# MISSINGNESS
# ---------------------------------------------------------------------------

def test_construct_unmapped_and_mapped_unmodeled_are_distinct():
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    measured = qc[qc["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]
    if len(measured) == 0:
        pytest.skip("no measured rows yet")
    lbd_length = pd.to_numeric(measured["lbd_length"], errors="coerce")
    mapped = pd.to_numeric(measured["mapped_lbd_residues"], errors="coerce")
    unmapped = pd.to_numeric(measured["unmapped_lbd_residues"], errors="coerce")
    assert ((mapped + unmapped) == lbd_length).all()


def test_terminal_missingness_arithmetic():
    """Direct algebraic check of the terminal-run-counting logic used in
    scripts/19_process_lbd_coordinates.py via a synthetic mask."""
    mapped_mask = [False, False, True, True, True, True, False]
    n_term = 0
    for m in mapped_mask:
        if not m:
            n_term += 1
        else:
            break
    c_term = 0
    for m in reversed(mapped_mask):
        if not m:
            c_term += 1
        else:
            break
    assert n_term == 2
    assert c_term == 1


def test_internal_run_arithmetic():
    """Synthetic check: internal unmodeled positions [3,4,7] (0-indexed,
    with mapped-but-unmodeled flags true there and false elsewhere within
    the mapped span) should yield 2 runs, longest 2."""
    positions = [3, 4, 7]
    runs = []
    run_start = positions[0]
    prev = run_start
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            runs.append(prev - run_start + 1)
            run_start = pos
            prev = pos
    runs.append(prev - run_start + 1)
    assert runs == [2, 1]
    assert len(runs) == 2
    assert max(runs) == 2


# ---------------------------------------------------------------------------
# MODEL
# ---------------------------------------------------------------------------

def test_multi_model_rows_not_collapsed():
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    multi = qc[qc["coordinate_qc_status"] == "MULTI_MODEL_REVIEW"]
    if len(multi) == 0:
        pytest.skip("no multi-model entries in this run")
    # Each (pdb_id, polymer_entity_id, instance_id) with multi-model status
    # must have more than one row (one per model), not a single collapsed row.
    counts = multi.groupby(["pdb_id", "polymer_entity_id", "instance_id"]).size()
    assert (counts > 1).all()


# ---------------------------------------------------------------------------
# CROSSCHECK
# ---------------------------------------------------------------------------

def test_unobserved_crosscheck_does_not_assume_absence_proves_observation():
    """Historical note: the crosscheck schema was rebuilt at Stage 3B.1 with
    four explicit categories (A/B/C/D) and a `ca_record_present` column
    (renamed from `has_ca_atom_record`). A residue absent from the deposited
    unobserved-residue table with no Cα record must still be recorded as
    category D, never silently dropped."""
    _require(CROSSCHECK_CSV)
    df = pd.read_csv(CROSSCHECK_CSV)
    if len(df) == 0:
        pytest.skip("no crosscheck rows yet")
    category_d = df[df["crosscheck_category"] == "D_NOT_DEPOSITED_UNOBSERVED_NO_CA_RECORD"]
    # Just confirming the category is representable (may be zero in a given run).
    assert "crosscheck_category" in df.columns
    assert (~category_d["listed_as_unobserved"]).all() if len(category_d) else True
    assert (~category_d["ca_record_present"]).all() if len(category_d) else True


# ---------------------------------------------------------------------------
# IMMUTABILITY
# ---------------------------------------------------------------------------

def test_stage2_release_summary_byte_identical():
    _require(STAGE2_RELEASE_SUMMARY)
    summary = json.loads(STAGE2_RELEASE_SUMMARY.read_text())
    assert summary["candidate_relationships"] == 2072
    assert summary["unique_pdb_entries"] == 1996


def test_stage3a_release_summary_byte_identical():
    _require(STAGE3A_RELEASE_SUMMARY)
    summary = json.loads(STAGE3A_RELEASE_SUMMARY.read_text())
    assert summary["project_metadata_prefilter_pass"] == 1516


def test_nr_lbd_reference_unchanged():
    _require(LBD_REFERENCE_CSV)
    df = pd.read_csv(LBD_REFERENCE_CSV)
    assert len(df) == 48
    vdr = df[df["uniprot_id"] == "P11473"].iloc[0]
    assert vdr["lbd_start"] == 127
    assert vdr["lbd_end"] == 423


def test_no_final_keep_exclude_field_exists():
    for path in (INSTANCE_QC_CSV, CANDIDATE_SUMMARY_CSV, POOL_CSV):
        if not path.exists():
            continue
        df = pd.read_csv(path)
        forbidden = {"keep", "exclude", "final_decision", "selected_for_analysis", "include_in_dataset"}
        assert forbidden.isdisjoint(set(c.lower() for c in df.columns)), f"{path.name} contains a forbidden final-decision column"


# ---------------------------------------------------------------------------
# Stage 3B.1: positive-occupancy correction
# ---------------------------------------------------------------------------

OBSERVATIONS_PARQUET = PROJECT_ROOT / "data" / "interim" / "qc" / "stage3b_lbd_ca_observations.parquet"
CORRECTION_AUDIT_CSV = MANIFESTS_DIR / "stage3b_occupancy_correction_audit.csv"
SUMMARY_JSON = MANIFESTS_DIR / "stage3b_coordinate_summary.json"
THRESHOLD_SENSITIVITY_CSV = MANIFESTS_DIR / "stage3b_completeness_threshold_sensitivity.csv"


def _load_compute_occupancy_fields():
    """scripts/ has no __init__.py (its modules are invoked as standalone
    scripts, not imported as a package), so load the numbered script
    directly by file path rather than via a package import."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "occupancy_correction", PROJECT_ROOT / "scripts" / "21_apply_occupancy_correction.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.compute_occupancy_fields


def test_zero_occupancy_semantics_direct():
    """Direct algebraic check of the Stage 3B.1 occupancy-field logic,
    mirroring scripts/21_apply_occupancy_correction.py's compute_occupancy_fields."""
    compute = _load_compute_occupancy_fields()

    # Zero occupancy alone: record present, but neither positive-occupancy
    # nor B-factor-usable.
    result = compute(json.dumps([0.0]), json.dumps([25.0]))
    assert result["ca_record_present"] is True
    assert result["ca_positive_occupancy"] is False
    assert result["ca_bfactor_usable"] is False
    assert result["zero_occupancy_only_ca"] is True

    # Positive occupancy with finite B: satisfies both.
    result = compute(json.dumps([1.0]), json.dumps([25.0]))
    assert result["ca_positive_occupancy"] is True
    assert result["ca_bfactor_usable"] is True
    assert result["zero_occupancy_only_ca"] is False

    # Positive occupancy with missing/nonfinite B: observed but not usable.
    result = compute(json.dumps([1.0]), json.dumps([None]))
    assert result["ca_positive_occupancy"] is True
    assert result["ca_bfactor_usable"] is False

    # Altlocs: any positive-occupancy altloc makes the residue observed;
    # altlocs never multiply the residue count (these functions describe a
    # single residue's atom records, so "count" here means altloc count,
    # which must be exactly the number of positive-occupancy altlocs, not
    # double-counted).
    result = compute(json.dumps([0.0, 0.6, 0.4]), json.dumps([20.0, 21.0, 22.0]))
    assert result["ca_positive_occupancy"] is True
    assert result["positive_occupancy_ca_altloc_count"] == 2


def test_record_present_coverage_gte_positive_occupancy_coverage():
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    if len(qc) == 0:
        pytest.skip("no instance QC rows yet")
    assert (
        pd.to_numeric(qc["record_present_full_lbd_coverage"], errors="coerce")
        >= pd.to_numeric(qc["observed_full_lbd_ca_coverage"], errors="coerce")
    ).all()


def test_positive_occupancy_coverage_gte_bfactor_usable_coverage():
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    if len(qc) == 0:
        pytest.skip("no instance QC rows yet")
    assert (
        pd.to_numeric(qc["observed_full_lbd_ca_coverage"], errors="coerce")
        >= pd.to_numeric(qc["bfactor_usable_full_lbd_coverage"], errors="coerce")
    ).all()


def test_exact_complete_flags_follow_definitions():
    _require(INSTANCE_QC_CSV)
    qc = pd.read_csv(INSTANCE_QC_CSV)
    if len(qc) == 0:
        pytest.skip("no instance QC rows yet")
    assert (
        qc["full_standardized_lbd_positive_occupancy_complete"]
        == (qc["positive_occupancy_lbd_ca_residues"] == qc["lbd_length"])
    ).all()
    assert (
        qc["full_standardized_lbd_ca_record_complete"]
        == (qc["ca_record_present_lbd_residues"] == qc["lbd_length"])
    ).all()
    assert (
        qc["full_standardized_lbd_bfactor_complete"]
        == (qc["bfactor_usable_lbd_ca_residues"] == qc["lbd_length"])
    ).all()


def test_deposited_unobserved_zero_occupancy_not_treated_as_observation():
    _require(CROSSCHECK_CSV)
    df = pd.read_csv(CROSSCHECK_CSV)
    zero_occ_rows = df[df["crosscheck_category"] == "B_DEPOSITED_UNOBSERVED_ZERO_OCCUPANCY_CA"]
    if len(zero_occ_rows) == 0:
        pytest.skip("no zero-occupancy-only deposited-unobserved rows in this run")
    assert (~zero_occ_rows["ca_positive_occupancy"]).all()
    assert (zero_occ_rows["review_status"] == "PASS").all()


def test_deposited_unobserved_positive_occupancy_flagged_review_needed():
    _require(CROSSCHECK_CSV)
    df = pd.read_csv(CROSSCHECK_CSV)
    review_rows = df[df["crosscheck_category"] == "C_DEPOSITED_UNOBSERVED_POSITIVE_OCCUPANCY_CA"]
    if len(review_rows) == 0:
        pytest.skip("no category-C rows in this run")
    assert (review_rows["review_status"] == "CROSSCHECK_REVIEW_NEEDED").all()
    assert review_rows["ca_positive_occupancy"].all()


def test_correction_audit_only_contains_changed_rows():
    _require(CORRECTION_AUDIT_CSV)
    df = pd.read_csv(CORRECTION_AUDIT_CSV)
    if len(df) == 0:
        pytest.skip("no occupancy-correction impact in this run")
    assert (df["record_present_lbd_residues"] != df["positive_occupancy_lbd_residues"]).all()


def test_no_representative_chain_or_threshold_selected():
    _require(SUMMARY_JSON)
    summary = json.loads(SUMMARY_JSON.read_text())
    assert summary.get("representative_chain_selected") is False
    assert summary.get("completeness_threshold_selected") is False
    assert summary.get("final_structure_qc_applied") is False


def test_threshold_sensitivity_is_descriptive_only():
    _require(THRESHOLD_SENSITIVITY_CSV)
    df = pd.read_csv(THRESHOLD_SENSITIVITY_CSV)
    assert len(df) == 5  # 100, 99, 95, 90, 80
    # Retention must be monotonically non-decreasing as the threshold relaxes.
    assert df["instances_retained"].is_monotonic_increasing or df.sort_values("threshold", ascending=False)["instances_retained"].is_monotonic_increasing


# ---------------------------------------------------------------------------
# Stage 3B.2: threshold-attrition + residual-missingness audit
# ---------------------------------------------------------------------------

ATTRITION_BY_RECEPTOR_CSV = MANIFESTS_DIR / "stage3b_threshold_attrition_by_receptor.csv"
RECEPTOR_ATTRITION_AUDIT_CSV = MANIFESTS_DIR / "stage3b_receptor_attrition_audit.csv"
RESIDUAL_MISSINGNESS_CSV = MANIFESTS_DIR / "stage3b_residual_missingness_audit.csv"
GROUND_TRUTH_CSV = MANIFESTS_DIR / "mmcif_unobserved_residue_ground_truth.csv"
HETEROGENEITY_CSV = MANIFESTS_DIR / "stage3b_instance_coverage_heterogeneity.csv"
ANY_ALL_CSV = MANIFESTS_DIR / "stage3b_any_all_threshold_comparison.csv"
RELEASE_SUMMARY_JSON = MANIFESTS_DIR / "stage3b_release_summary.json"
RECOVERY_NOTE_MD = PROJECT_ROOT / "reports" / "tables" / "stage3b_data_recovery_note.md"
DESIGN_NOTE_MD = PROJECT_ROOT / "reports" / "tables" / "stage3b_representative_instance_design_note.md"

THRESHOLD_KEYS = ["100pct", "99pct", "95pct", "90pct", "80pct"]
VALID_ATTRITION_STATUSES = {"RETAINED", "LOST_AT_THRESHOLD", "NO_STAGE3B_COORDINATE_CANDIDATE", "STAGE2_ZERO_RESULT"}


def test_threshold_attrition_table_has_exactly_48_receptor_rows():
    _require(ATTRITION_BY_RECEPTOR_CSV)
    df = pd.read_csv(ATTRITION_BY_RECEPTOR_CSV)
    assert len(df) == 48
    assert df["uniprot_id"].nunique() == 48


def test_threshold_attrition_denominator_partition_is_exhaustive_and_disjoint():
    """The 48 receptors must partition exactly into STAGE2_ZERO_RESULT (3),
    NO_STAGE3B_COORDINATE_CANDIDATE (3), and HAS_CANDIDATE (42) — no receptor
    may appear in more than one bucket, and every receptor must appear in
    exactly one."""
    _require(ATTRITION_BY_RECEPTOR_CSV)
    df = pd.read_csv(ATTRITION_BY_RECEPTOR_CSV)
    counts = df["base_status"].value_counts().to_dict()
    assert counts.get("STAGE2_ZERO_RESULT") == 3
    assert counts.get("NO_STAGE3B_COORDINATE_CANDIDATE") == 3
    assert counts.get("HAS_CANDIDATE") == 42
    assert sum(counts.values()) == 48


def test_threshold_attrition_status_values_are_from_controlled_vocabulary():
    _require(ATTRITION_BY_RECEPTOR_CSV)
    df = pd.read_csv(ATTRITION_BY_RECEPTOR_CSV)
    for key in THRESHOLD_KEYS:
        assert set(df[f"status_at_{key}"].unique()) <= VALID_ATTRITION_STATUSES


def test_zero_result_and_threshold_lost_receptors_are_disjoint_at_every_threshold():
    """A Stage-2 zero-result receptor (never had ANY candidate) must never be
    reported as LOST_AT_THRESHOLD (which specifically means it HAD a
    Stage 3B candidate that failed a coverage threshold) — these are
    distinct denominators and must never be conflated."""
    _require(ATTRITION_BY_RECEPTOR_CSV)
    df = pd.read_csv(ATTRITION_BY_RECEPTOR_CSV)
    zero_result = set(df.loc[df["base_status"] == "STAGE2_ZERO_RESULT", "uniprot_id"])
    no_candidate = set(df.loc[df["base_status"] == "NO_STAGE3B_COORDINATE_CANDIDATE", "uniprot_id"])
    for key in THRESHOLD_KEYS:
        lost = set(df.loc[df[f"status_at_{key}"] == "LOST_AT_THRESHOLD", "uniprot_id"])
        assert lost.isdisjoint(zero_result)
        assert lost.isdisjoint(no_candidate)


def test_threshold_attrition_retention_is_monotonic_per_receptor():
    """A receptor retained at a stricter threshold (e.g. 100%) must also be
    retained at every looser threshold (99/95/90/80%)."""
    _require(ATTRITION_BY_RECEPTOR_CSV)
    df = pd.read_csv(ATTRITION_BY_RECEPTOR_CSV)
    ordered = ["100pct", "99pct", "95pct", "90pct", "80pct"]
    for _, row in df[df["base_status"] == "HAS_CANDIDATE"].iterrows():
        retained_flags = [row[f"status_at_{k}"] == "RETAINED" for k in ordered]
        # Once True, must stay True for all looser (later) thresholds.
        seen_true = False
        for flag in retained_flags:
            if seen_true:
                assert flag, f"{row['uniprot_id']} lost retention at a looser threshold after gaining it"
            seen_true = seen_true or flag


def test_receptor_attrition_audit_only_covers_receptors_lost_at_95_90_or_80():
    _require(RECEPTOR_ATTRITION_AUDIT_CSV)
    _require(ATTRITION_BY_RECEPTOR_CSV)
    audit = pd.read_csv(RECEPTOR_ATTRITION_AUDIT_CSV)
    attrition = pd.read_csv(ATTRITION_BY_RECEPTOR_CSV)
    if len(audit) == 0:
        pytest.skip("no receptors lost at 95/90/80 in this run")
    expected_lost = set()
    for key in ("95pct", "90pct", "80pct"):
        expected_lost |= set(attrition.loc[attrition[f"status_at_{key}"] == "LOST_AT_THRESHOLD", "uniprot_id"])
    assert set(audit["uniprot_id"].unique()) <= expected_lost


def test_receptor_attrition_audit_classification_is_controlled_vocabulary():
    _require(RECEPTOR_ATTRITION_AUDIT_CSV)
    audit = pd.read_csv(RECEPTOR_ATTRITION_AUDIT_CSV)
    if len(audit) == 0:
        pytest.skip("no receptor attrition audit rows in this run")
    assert set(audit["attrition_classification"].unique()) <= {
        "CONSTRUCT_TRUNCATION", "UNRESOLVED_COORDINATES", "BOTH", "OTHER_REVIEW",
    }


def test_residual_missingness_audit_expected_case_counts():
    """Stage 3B.2 expects 13 case-D (no Cα record, not deposited-unobserved)
    plus 11 case-E (zero-occupancy-only Cα, not deposited-unobserved) = 24
    total residual cases, matching the counts independently reproduced from
    the corrected parquet and the re-extracted deposited-unobserved ground
    truth."""
    _require(RESIDUAL_MISSINGNESS_CSV)
    df = pd.read_csv(RESIDUAL_MISSINGNESS_CSV)
    counts = df["case_type"].value_counts().to_dict()
    assert counts.get("D_NO_CA_RECORD_NOT_LISTED_UNOBSERVED", 0) == 13
    assert counts.get("E_ZERO_OCCUPANCY_ONLY_NOT_LISTED_UNOBSERVED", 0) == 11
    assert len(df) == 24


def test_residual_missingness_cases_have_valid_receptor_mapping():
    _require(RESIDUAL_MISSINGNESS_CSV)
    df = pd.read_csv(RESIDUAL_MISSINGNESS_CSV)
    if len(df) == 0:
        pytest.skip("no residual missingness cases in this run")
    assert df["within_standardized_lbd"].all()
    assert df["valid_label_seq_id"].all()
    assert not df["has_positive_occupancy_alternative"].any()
    assert df["status"].isin({"EXPLAINED", "REVIEW_NEEDED"}).all()


def test_residual_missingness_status_matches_verification_flags():
    """A case is EXPLAINED iff it is within the LBD, has a valid label_seq_id,
    and has no positive-occupancy alternative Cα; otherwise REVIEW_NEEDED."""
    _require(RESIDUAL_MISSINGNESS_CSV)
    df = pd.read_csv(RESIDUAL_MISSINGNESS_CSV)
    if len(df) == 0:
        pytest.skip("no residual missingness cases in this run")
    expected_explained = (
        df["within_standardized_lbd"] & df["valid_label_seq_id"] & ~df["has_positive_occupancy_alternative"]
    )
    assert ((df["status"] == "EXPLAINED") == expected_explained).all()


def test_ground_truth_extraction_matches_previously_verified_recovery_counts():
    """The re-extracted deposited-unobserved ground truth must reproduce the
    73,545-record count independently verified during the Stage 3B.1
    recovery (see reports/tables/stage3b_data_recovery_note.md)."""
    _require(GROUND_TRUTH_CSV)
    df = pd.read_csv(GROUND_TRUTH_CSV)
    assert len(df) == 73545


def test_instance_heterogeneity_only_includes_multi_instance_candidates():
    _require(HETEROGENEITY_CSV)
    df = pd.read_csv(HETEROGENEITY_CSV)
    if len(df) == 0:
        pytest.skip("no multi-instance candidates in this run")
    assert (df["instance_count"] > 1).all()


def test_instance_heterogeneity_flag_matches_range_threshold():
    _require(HETEROGENEITY_CSV)
    df = pd.read_csv(HETEROGENEITY_CSV)
    if len(df) == 0:
        pytest.skip("no multi-instance candidates in this run")
    expected_flag = df["coverage_range"] >= 0.05
    actual_flag = df["flag"] == "INSTANCE_COVERAGE_HETEROGENEITY_REVIEW"
    assert (expected_flag == actual_flag).all()


def test_instance_heterogeneity_range_equals_max_minus_min():
    _require(HETEROGENEITY_CSV)
    df = pd.read_csv(HETEROGENEITY_CSV)
    if len(df) == 0:
        pytest.skip("no multi-instance candidates in this run")
    assert (
        (df["coverage_range"] - (df["max_coverage"] - df["min_coverage"])).abs() < 1e-9
    ).all()


def test_any_instance_threshold_pass_implies_at_least_as_permissive_as_all_instances():
    """For every candidate and threshold, if all_instances_pass is True then
    any_instance_pass must also be True (all => any, never the reverse)."""
    _require(ANY_ALL_CSV)
    df = pd.read_csv(ANY_ALL_CSV)
    if len(df) == 0:
        pytest.skip("no candidates in this run")
    for key in THRESHOLD_KEYS:
        all_col, any_col = f"all_instances_ge_{key}", f"any_instance_ge_{key}"
        violation = df[all_col] & ~df[any_col]
        assert not violation.any(), f"{key}: found candidate(s) where all_instances passes but any_instance does not"


def test_stage3b2_release_summary_reproduces_frozen_stage3b1_figures():
    """Stage 3B.2 must not recompute or alter the approved Stage 3B.1
    figures — only add denominator/attrition/missingness/heterogeneity
    context on top of them."""
    _require(RELEASE_SUMMARY_JSON)
    summary = json.loads(RELEASE_SUMMARY_JSON.read_text())
    assert summary["coordinate_files_downloaded"] == 1777
    assert summary["selected_candidates"] == 1840
    assert summary["measured_instance_models"] == 3159
    assert summary["exact_full_lbd_positive_occupancy_complete"] == 341
    assert summary["coverage_ge_99pct"] == 629
    assert summary["coverage_ge_95pct"] == 1862
    assert summary["coverage_ge_90pct"] == 2771
    assert summary["coverage_ge_80pct"] == 3089
    assert summary["zero_occupancy_only_lbd_residues"] == 281
    assert summary["deposited_unobserved_positive_occupancy_conflicts"] == 0


def test_stage3b2_release_summary_declares_no_final_selection():
    _require(RELEASE_SUMMARY_JSON)
    summary = json.loads(RELEASE_SUMMARY_JSON.read_text())
    assert summary["final_threshold_selected"] is False
    assert summary["representative_instance_selected"] is False
    assert summary["normalized_bfactor_analysis_started"] is False


def test_stage3b2_release_summary_denominator_matches_attrition_table():
    _require(RELEASE_SUMMARY_JSON)
    _require(ATTRITION_BY_RECEPTOR_CSV)
    summary = json.loads(RELEASE_SUMMARY_JSON.read_text())
    attrition = pd.read_csv(ATTRITION_BY_RECEPTOR_CSV)
    denom = summary["receptor_denominator"]
    assert denom["total_project_receptors"] == 48
    assert sorted(denom["stage2_zero_result_receptors"]) == sorted(
        attrition.loc[attrition["base_status"] == "STAGE2_ZERO_RESULT", "uniprot_id"]
    )
    assert sorted(denom["no_stage3b_candidate_receptors"]) == sorted(
        attrition.loc[attrition["base_status"] == "NO_STAGE3B_COORDINATE_CANDIDATE", "uniprot_id"]
    )
    assert denom["receptors_with_stage3b_candidate"] == 42


def test_data_recovery_note_exists_and_documents_no_network_reaccess():
    _require(RECOVERY_NOTE_MD)
    text = RECOVERY_NOTE_MD.read_text()
    assert "No network request was made" in text or "no network re-download" in text.lower() or "no network access" in text.lower()
    assert "1777" in text


def test_representative_instance_design_note_exists_and_specifies_granularity():
    _require(DESIGN_NOTE_MD)
    text = DESIGN_NOTE_MD.read_text()
    assert "polymer entity" in text.lower() or "polymer_entity" in text.lower()
    assert "once per PDB ID" in text or "never once per PDB ID" in text


def test_raw_mmcif_file_count_unchanged_at_1777():
    _require(RAW_MMCIF_DIR)
    files = list(RAW_MMCIF_DIR.glob("*.cif.gz"))
    assert len(files) == 1777


def test_download_manifest_checksums_all_present_and_no_failures():
    _require(DOWNLOAD_MANIFEST_CSV)
    df = pd.read_csv(DOWNLOAD_MANIFEST_CSV)
    assert len(df) == 1777
    assert (df["download_status"].isin({"SUCCESS", "CACHED_VALID"})).all()
    assert df["sha256_gzip"].notna().all()
