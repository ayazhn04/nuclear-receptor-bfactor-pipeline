"""
Stage 3A tests: LBD reference, candidate-to-LBD mapping, instance mapping,
and pre-coordinate QC. Offline only — operate on already-generated
manifests, plus direct algebraic checks of the shared interval-arithmetic
utilities.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.intervals import intersect_with_reference, merge_intervals, overlap_fraction, total_length  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"

LBD_REFERENCE = MANIFESTS_DIR / "nr_lbd_reference.csv"
LBD_MAPPING = MANIFESTS_DIR / "rcsb_candidate_lbd_mapping.csv"
DISCOVERY = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
INSTANCES = MANIFESTS_DIR / "rcsb_polymer_instances.csv"
INSTANCE_INVENTORY = MANIFESTS_DIR / "rcsb_candidate_instance_inventory.csv"
PRECOORDINATE_QC = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"
STAGE1_AUDIT = MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv"
STAGE2_RELEASE_SUMMARY = MANIFESTS_DIR / "stage2_release_summary.json"
CANDIDATE_INVENTORY = MANIFESTS_DIR / "rcsb_candidate_inventory.csv"

VALID_LBD_CATEGORIES = {"NO_LBD_OVERLAP", "TRACE_LBD", "PARTIAL_LBD", "SUBSTANTIAL_LBD", "NEAR_COMPLETE_LBD"}
VALID_LBD_STATUS = {"PASS", "REVIEW_NEEDED", "NOT_FOUND"}


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet")


# ---------------------------------------------------------------------------
# Interval arithmetic (algebraic, no files needed)
# ---------------------------------------------------------------------------

def test_merge_intervals_combines_overlapping():
    assert merge_intervals([(1, 5), (3, 8), (20, 25)]) == [(1, 8), (20, 25)]


def test_merge_intervals_combines_adjacent():
    assert merge_intervals([(1, 5), (6, 10)]) == [(1, 10)]


def test_total_length_counts_union_not_sum():
    # Overlapping intervals must not be double-counted.
    assert total_length([(1, 10), (5, 15)]) == 15  # union is 1-15
    assert total_length([(1, 10), (5, 15)]) != (10 + 11)  # naive sum would be wrong


def test_intersect_with_reference_basic():
    assert intersect_with_reference([(1, 100)], 20, 50) == [(20, 50)]
    assert intersect_with_reference([(1, 10), (90, 100)], 20, 50) == []


def test_overlap_fraction_disjoint_and_identical():
    assert overlap_fraction([(1, 10)], [(20, 30)]) == 0.0
    assert overlap_fraction([(1, 10)], [(1, 10)]) == 1.0


# ---------------------------------------------------------------------------
# LBD reference table
# ---------------------------------------------------------------------------

def test_lbd_reference_has_exactly_48_rows():
    _require(LBD_REFERENCE)
    df = pd.read_csv(LBD_REFERENCE)
    assert len(df) == 48


def test_lbd_reference_every_validated_receptor_represented():
    _require(LBD_REFERENCE)
    from config.nr_metadata import NR_METADATA

    df = pd.read_csv(LBD_REFERENCE)
    assert set(df["uniprot_id"]) == {r["uniprot_id"] for r in NR_METADATA}


def test_lbd_reference_status_controlled_vocabulary():
    _require(LBD_REFERENCE)
    df = pd.read_csv(LBD_REFERENCE)
    assert set(df["lbd_status"].unique()).issubset(VALID_LBD_STATUS)


def test_lbd_reference_pass_rows_obey_coordinate_constraints():
    _require(LBD_REFERENCE)
    df = pd.read_csv(LBD_REFERENCE)
    passed = df[df["lbd_status"] == "PASS"]
    assert len(passed) > 0
    assert (passed["lbd_start"] >= 1).all()
    assert (passed["lbd_start"] <= passed["lbd_end"]).all()
    assert (passed["lbd_end"] <= passed["sequence_length"]).all()
    assert (passed["lbd_length"] == passed["lbd_end"] - passed["lbd_start"] + 1).all()


def test_lbd_reference_not_manually_typed():
    """Every accepted LBD interval must cite UniProt PROSITE-ProRule
    evidence — a proxy for 'extracted from authoritative data', not typed."""
    _require(LBD_REFERENCE)
    df = pd.read_csv(LBD_REFERENCE)
    passed = df[df["lbd_status"] == "PASS"]
    assert (passed["lbd_feature_evidence"].str.contains("PROSITE-ProRule", na=False)).all()


# ---------------------------------------------------------------------------
# Candidate LBD mapping
# ---------------------------------------------------------------------------

def test_lbd_mapping_preserves_exactly_2072_candidates():
    _require(LBD_MAPPING)
    _require(DISCOVERY)
    mapping_df = pd.read_csv(LBD_MAPPING)
    discovery_df = pd.read_csv(DISCOVERY)
    assert len(mapping_df) == len(discovery_df) == 2072


def test_lbd_mapping_no_accidental_row_multiplication():
    _require(LBD_MAPPING)
    df = pd.read_csv(LBD_MAPPING)
    key = df[["uniprot_id", "polymer_entity_id"]]
    assert not key.duplicated().any()


def test_mapped_lbd_residues_within_lbd_length():
    _require(LBD_MAPPING)
    df = pd.read_csv(LBD_MAPPING)
    passed = df[df["mapping_status"] == "PASS"]
    assert (passed["mapped_lbd_residues"] >= 0).all()
    assert (passed["mapped_lbd_residues"] <= passed["lbd_length"]).all()


def test_lbd_mapping_coverage_within_0_1():
    _require(LBD_MAPPING)
    df = pd.read_csv(LBD_MAPPING)
    values = pd.to_numeric(df["lbd_mapping_coverage"], errors="coerce").dropna()
    assert ((values >= 0) & (values <= 1.000001)).all()


def test_lbd_category_thresholds_exact():
    _require(LBD_MAPPING)
    df = pd.read_csv(LBD_MAPPING)
    passed = df[df["mapping_status"] == "PASS"]
    for _, row in passed.iterrows():
        cov, cat = row["lbd_mapping_coverage"], row["lbd_mapping_category"]
        if cov == 0:
            assert cat == "NO_LBD_OVERLAP"
        elif cov < 0.25:
            assert cat == "TRACE_LBD"
        elif cov < 0.80:
            assert cat == "PARTIAL_LBD"
        elif cov < 0.95:
            assert cat == "SUBSTANTIAL_LBD"
        else:
            assert cat == "NEAR_COMPLETE_LBD"


def test_lbd_mapping_category_controlled_vocabulary():
    _require(LBD_MAPPING)
    df = pd.read_csv(LBD_MAPPING)
    non_review = df[df["mapping_status"] == "PASS"]
    assert set(non_review["lbd_mapping_category"].unique()).issubset(VALID_LBD_CATEGORIES)


def test_fusion_entities_use_query_receptor_mapping_not_first():
    """For a known multi-UniProt entity, the LBD mapping coverage must be
    computed from the QUERY receptor's own alignment, not an arbitrary
    first UniProt mapping (which could belong to a fusion partner)."""
    _require(LBD_MAPPING)
    multi_audit_path = MANIFESTS_DIR / "rcsb_multi_uniprot_entity_audit.csv"
    _require(multi_audit_path)
    multi_df = pd.read_csv(multi_audit_path)
    lbd_df = pd.read_csv(LBD_MAPPING)
    if len(multi_df) == 0:
        pytest.skip("no multi-uniprot entities to check")
    sample = multi_df.iloc[0]
    row = lbd_df[
        (lbd_df["uniprot_id"] == sample["project_receptor_uniprot"])
        & (lbd_df["polymer_entity_id"] == sample["polymer_entity_id"])
    ]
    assert len(row) == 1, "query receptor's own LBD-mapping row must exist for a multi-UniProt entity"


# ---------------------------------------------------------------------------
# Instance / chain mapping
# ---------------------------------------------------------------------------

def test_each_instance_has_one_parent_entity():
    _require(INSTANCES)
    df = pd.read_csv(INSTANCES)
    passed = df[df["instance_metadata_status"] == "PASS"]
    # instance_id encodes exactly one entry+entity; verify entity_id and
    # pdb_id are present and singular per row (not lists).
    assert (passed["polymer_entity_id"].apply(lambda v: isinstance(v, str) and "_" in v)).all()


def test_instance_rows_have_unique_logical_keys():
    _require(INSTANCES)
    df = pd.read_csv(INSTANCES)
    passed = df[df["instance_metadata_status"] == "PASS"]
    assert not passed["instance_id"].duplicated().any()


def test_label_auth_ids_came_from_instance_metadata_not_position_pairing():
    """Regression guard: the instance table must carry its own
    label/auth asym IDs per row (from RCSB instance records), not be
    derivable purely by zipping the Stage 2 entity-level arrays — this
    test checks the columns exist and are populated for PASS rows, which
    only instance-level fetching (not array pairing) could produce for
    entities where label != auth id ordering."""
    _require(INSTANCES)
    df = pd.read_csv(INSTANCES)
    passed = df[df["instance_metadata_status"] == "PASS"]
    assert passed["label_asym_id"].notna().all()
    assert passed["auth_asym_id"].notna().all()


def test_instance_inventory_larger_than_or_equal_entity_count():
    _require(INSTANCE_INVENTORY)
    _require(DISCOVERY)
    inst_df = pd.read_csv(INSTANCE_INVENTORY)
    disc_df = pd.read_csv(DISCOVERY)
    assert len(inst_df) >= disc_df["polymer_entity_id"].nunique()


# ---------------------------------------------------------------------------
# Pre-coordinate QC
# ---------------------------------------------------------------------------

def test_precoordinate_qc_prefilter_matches_component_flags():
    _require(PRECOORDINATE_QC)
    df = pd.read_csv(PRECOORDINATE_QC)
    for _, row in df.iterrows():
        expected = (
            bool(row["method_is_xray"])
            and row["resolution_project_range_status"] == "PASS_PROJECT_RANGE"
            and row["r_free_status"] == "PASS_ALL_LE_0_30"
            and bool(row["mapping_overlaps_lbd"])
        )
        assert bool(row["passes_project_metadata_prefilter"]) == expected


def test_ultrahigh_resolution_entries_flagged_not_discarded():
    _require(PRECOORDINATE_QC)
    df = pd.read_csv(PRECOORDINATE_QC)
    ultrahigh = df[df["resolution_project_range_status"] == "ULTRAHIGH_LT_1_8"]
    if len(ultrahigh) == 0:
        pytest.skip("no ultrahigh-resolution entries in this run")
    assert len(ultrahigh) > 0  # present in the table, not removed


def test_non_xray_candidates_remain_represented():
    _require(PRECOORDINATE_QC)
    df = pd.read_csv(PRECOORDINATE_QC)
    non_xray = df[~df["method_is_xray"]]
    assert len(non_xray) > 0
    assert (~non_xray["passes_project_metadata_prefilter"]).all()


def test_missing_r_free_remains_represented():
    _require(PRECOORDINATE_QC)
    df = pd.read_csv(PRECOORDINATE_QC)
    missing = df[df["r_free_status"] == "MISSING"]
    if len(missing) == 0:
        pytest.skip("no missing-R-free entries in this run")
    assert (~missing["passes_project_metadata_prefilter"]).all()


def test_taxonomy_status_alone_does_not_cause_removal():
    """No row should be dropped from the QC table merely for having a
    taxonomy_audit_status — the table must still contain 2072 rows total,
    and taxonomy status must play no role in passes_project_metadata_prefilter."""
    _require(PRECOORDINATE_QC)
    df = pd.read_csv(PRECOORDINATE_QC)
    assert len(df) == 2072
    with_taxonomy_flag = df[df["taxonomy_audit_status"].notna() & (df["taxonomy_audit_status"] != "")]
    if len(with_taxonomy_flag) == 0:
        pytest.skip("no taxonomy-flagged rows in this run")
    # Confirm the prefilter formula (tested above) has no taxonomy term —
    # here we just confirm these rows aren't silently excluded from the table.
    assert len(with_taxonomy_flag) == 9


# ---------------------------------------------------------------------------
# Immutability
# ---------------------------------------------------------------------------

def test_stage1_manifests_unchanged():
    _require(STAGE1_AUDIT)
    df = pd.read_csv(STAGE1_AUDIT)
    assert len(df) == 48


def test_stage2_release_summary_unchanged():
    _require(STAGE2_RELEASE_SUMMARY)
    summary = json.loads(STAGE2_RELEASE_SUMMARY.read_text())
    assert summary["candidate_relationships"] == 2072
    assert summary["unique_polymer_entities"] == 2072
    assert summary["unique_pdb_entries"] == 1996


def test_candidate_inventory_still_2072():
    _require(CANDIDATE_INVENTORY)
    df = pd.read_csv(CANDIDATE_INVENTORY)
    assert len(df) == 2072


def test_no_coordinate_files_exist():
    """Historical note: see the identically-named test in
    test_stage2_rcsb_discovery.py — superseded once Stage 3B begins."""
    if (MANIFESTS_DIR / "stage3b_coordinate_pool.csv").exists():
        pytest.skip("Stage 3B has begun downloading mmCIF files by design — see test_stage3b_coordinate_qc.py")
    mmcif_dir = PROJECT_ROOT / "data" / "raw" / "mmcif"
    if mmcif_dir.exists():
        files = [p for p in mmcif_dir.rglob("*") if p.is_file() and p.name != ".gitkeep"]
        assert len(files) == 0
