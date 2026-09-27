"""
Structural-segment alignment tests: secondary-structure consensus, canonical
H1-H12/beta-sheet/loop labelling, segment coverage/flexibility features, and
protection of the approved normalization foundation, historical strict-QC
files and Role D outputs.

Offline; uses only tracked files (no raw mmCIF/API caches needed except the
DSSP-consensus test, which needs the already-cached data/raw/mmcif files this
session copied in from the sibling clone -- skipped if absent).
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"
REPORTS = PROJECT_ROOT / "reports" / "tables"

MASTER = PROCESSED / "structures.csv"
NORMALIZED = ANALYSIS / "normalized_lbd_bfactors_all.parquet"
CONSENSUS = ANALYSIS / "receptor_dssp_consensus.parquet"
SEGMENT_MAP = ANALYSIS / "nr_lbd_structural_segment_map.csv"
VALIDATION_CSV = REPORTS / "structural_segment_validation.csv"
COVERAGE = ANALYSIS / "pdb_segment_coverage.parquet"
PDB_FEATURES = ANALYSIS / "pdb_segment_flexibility_features.parquet"
RECEPTOR_PROFILE = ANALYSIS / "receptor_segment_flexibility_profiles.parquet"
RECEPTOR_MATRIX_MEDIAN = ANALYSIS / "receptor_segment_feature_matrix_median.parquet"
RECEPTOR_MATRIX_MEAN = ANALYSIS / "receptor_segment_feature_matrix_mean.parquet"
COVERAGE_MATRIX = ANALYSIS / "receptor_segment_coverage_matrix.parquet"
ALIGNMENT_MAP = ANALYSIS / "nr_lbd_common_alignment_map.parquet"
LBD_REFERENCE = MANIFESTS / "nr_lbd_reference.csv"

# Baseline: approved analysis-foundation commit; this stage must not alter these.
BASELINE_COMMIT = "e714bbd"
NORMALIZATION_FILES = [
    "data/processed/structures.csv",
    "data/analysis/normalized_lbd_bfactors_all.parquet",
    "data/analysis/instance_bfactor_summary.csv",
    "data/analysis/pdb_position_bfactor_profiles.parquet",
    "data/analysis/receptor_position_flexibility_profiles.parquet",
    "data/analysis/receptor_flexibility_summary.csv",
]
PROTECTED_PATHS = [
    "data/processed/structures_primary_qc_subset.csv",
    "data/processed/final_lbd_residue_map_primary_qc_subset.parquet",
    "data/processed/lbd_residue_bfactors_all.parquet",
    "scripts/31_analyze_functional_site_bfactors.py",
    "reports/ROLE_D_PROGRESS_REPORT.md",
    "reports/tables/functional_site_definitions.csv",
    "reports/tables/functional_site_entry_scores.csv",
    "reports/tables/functional_site_normalization_qc.csv",
    "reports/tables/functional_site_receptor_scores.csv",
    "tests/test_functional_site_analysis.py",
]
VDR = "P11473"
ANCHORS = {"P11473": "VDR", "P37231": "PPARgamma", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet")


@pytest.fixture(scope="module")
def master():
    return pd.read_csv(MASTER, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def consensus():
    _require(CONSENSUS)
    return pd.read_parquet(CONSENSUS)


@pytest.fixture(scope="module")
def seg_map():
    _require(SEGMENT_MAP)
    return pd.read_csv(SEGMENT_MAP, keep_default_na=False)


@pytest.fixture(scope="module")
def pdb_features():
    _require(PDB_FEATURES)
    return pd.read_parquet(PDB_FEATURES)


@pytest.fixture(scope="module")
def receptor_profile():
    _require(RECEPTOR_PROFILE)
    return pd.read_parquet(RECEPTOR_PROFILE)


# ---------------------------------------------------------------- master / normalization protection

def test_no_master_candidate_removed(master):
    assert len(master) == 2072
    assert int(master["has_stage3b_coordinate_data"].sum()) == 1840


def test_normalization_outputs_byte_identical_to_approved_baseline():
    for rel in NORMALIZATION_FILES:
        actual = _sha(PROJECT_ROOT / rel)
        baseline = hashlib.sha256(subprocess.run(
            ["git", "show", f"{BASELINE_COMMIT}:{rel}"], cwd=PROJECT_ROOT, capture_output=True, check=True
        ).stdout).hexdigest()
        assert actual == baseline, f"{rel} differs from the approved normalization-foundation baseline"


def test_historical_and_role_d_files_unchanged():
    r = subprocess.run(["git", "diff", "--quiet", BASELINE_COMMIT, "--", *PROTECTED_PATHS], cwd=PROJECT_ROOT)
    assert r.returncode == 0


# ---------------------------------------------------------------- Part 3: DSSP consensus

def test_dssp_consensus_uses_valid_uniprot_positions(consensus, master):
    lbd_ref = pd.read_csv(LBD_REFERENCE).set_index("uniprot_id")
    for uid, g in consensus.groupby("uniprot_id"):
        lo, hi = lbd_ref.loc[uid, ["lbd_start", "lbd_end"]]
        assert g["canonical_uniprot_position"].between(lo, hi).all()
    assert not consensus.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()
    assert np.isclose(consensus["helix_fraction"] + consensus["sheet_fraction"] + consensus["coil_fraction"], 1.0).all()
    assert (consensus["n_pdb_observed"] >= 1).all()


def test_dssp_consensus_covers_receptors_with_coordinate_data(consensus, master):
    coord_receptors = set(master.loc[master["has_stage3b_coordinate_data"], "uniprot_id"])
    assert set(consensus["uniprot_id"]) <= coord_receptors
    assert len(consensus["uniprot_id"].unique()) >= 40  # most receptors with structures have secondary-structure evidence


# ---------------------------------------------------------------- Part 4: segment map

def test_each_receptor_position_has_at_most_one_final_segment_assignment(seg_map):
    assert not seg_map.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()


def test_ambiguous_positions_remain_explicit_not_guessed(seg_map):
    ambiguous = seg_map[seg_map["assignment_status"] == "AMBIGUOUS"]
    assert len(ambiguous) > 0
    assert (ambiguous["structural_segment"] == "").all()
    unassigned = seg_map[seg_map["assignment_status"] == "UNASSIGNED"]
    assert (unassigned["structural_segment"] == "").all()


def test_segment_order_biologically_consistent(seg_map):
    named = seg_map[seg_map["structural_segment"].str.fullmatch(r"H\d{1,2}")]
    for uid, g in named.groupby("uniprot_id"):
        g = g.drop_duplicates("structural_segment").assign(
            n=lambda d: d["structural_segment"].str[1:].astype(int)
        ).sort_values("n")
        starts = g.groupby("structural_segment")["canonical_uniprot_position"].min().reindex(g["structural_segment"])
        # each named helix must not start before the previous-numbered named helix
        assert starts.is_monotonic_increasing


def test_no_position_assigned_a_segment_name_without_secondary_structure_evidence(seg_map, consensus):
    named = seg_map[seg_map["structural_class"].isin(["HELIX", "SHEET"])]
    have_evidence = set(zip(consensus["uniprot_id"], consensus["canonical_uniprot_position"]))
    assert set(zip(named["uniprot_id"], named["canonical_uniprot_position"])) <= have_evidence


def test_h1_and_h12_are_the_most_reliably_assigned_segments(seg_map):
    high = seg_map[seg_map["assignment_status"] == "ASSIGNED_HIGH_CONFIDENCE"]
    counts = high.groupby("structural_segment")["uniprot_id"].nunique()
    assert counts.get("H1", 0) >= 30
    assert counts.get("H12", 0) >= 30
    for seg in ["H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11"]:
        assert counts.get(seg, 0) <= counts.get("H1", 0)  # never more reliably resolved than the anchors


# ---------------------------------------------------------------- Part 5: validation

def test_anchor_receptors_present_in_validation_table():
    _require(VALIDATION_CSV)
    v = pd.read_csv(VALIDATION_CSV)
    assert set(ANCHORS) & set(v["uniprot_id"]) == set(ANCHORS) - {"P22736"} | ({"P22736"} if "P22736" in set(v["uniprot_id"]) else set())
    for uid in ANCHORS:
        assert (v["uniprot_id"] == uid).any(), f"{ANCHORS[uid]} missing from structural_segment_validation.csv"


def test_vdr_h1_and_h12_validated_and_h12_matches_known_terminal_helix():
    _require(VALIDATION_CSV)
    v = pd.read_csv(VALIDATION_CSV)
    vdr = v[v["uniprot_id"] == VDR]
    h1 = vdr[vdr["structural_segment"] == "H1"].iloc[0]
    h12 = vdr[vdr["structural_segment"] == "H12"].iloc[0]
    assert h1["start_uniprot"] <= 132  # near lbd_start=127
    # repository-documented VDR terminal helix (nr_functional_features.csv): 416-422
    assert h12["start_uniprot"] <= 420 and h12["end_uniprot"] >= 418


# ---------------------------------------------------------------- Parts 6-9: coverage/features/matrices

def test_no_missing_feature_silently_converted_to_zero(receptor_profile):
    _require(RECEPTOR_MATRIX_MEDIAN)
    matrix = pd.read_parquet(RECEPTOR_MATRIX_MEDIAN).set_index("uniprot_id")
    assert matrix.isna().any().any()
    # a real, non-trivial number of receptor-median z-scores are exactly/near zero (genuine mean flexibility) --
    # confirms 0 is a meaningful value the pipeline must not confuse with missing.
    near_zero = (receptor_profile["median_of_entry_median_z"].abs() < 0.05).sum()
    assert near_zero > 0


def test_chain_to_pdb_aggregation_is_median(pdb_features):
    pos = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")
    seg_map_df = pd.read_csv(SEGMENT_MAP, keep_default_na=False)
    named = seg_map_df[seg_map_df["structural_segment"] != ""][["uniprot_id", "canonical_uniprot_position", "structural_segment"]]
    merged = pos.merge(named, on=["uniprot_id", "canonical_uniprot_position"])
    sample = pdb_features.dropna(subset=["median_z"]).sample(min(20, len(pdb_features)), random_state=3)
    for _, r in sample.iterrows():
        vals = merged[(merged["uniprot_id"] == r["uniprot_id"]) & (merged["pdb_id"] == r["pdb_id"]) & (merged["structural_segment"] == r["structural_segment"])]["b_factor_z"]
        assert np.isclose(r["median_z"], vals.median())


def test_pdb_to_receptor_aggregation_is_median_for_primary_matrix(receptor_profile, pdb_features):
    sample = receptor_profile.dropna(subset=["median_of_entry_median_z"]).sample(min(20, len(receptor_profile)), random_state=5)
    for _, r in sample.iterrows():
        vals = pdb_features[(pdb_features["uniprot_id"] == r["uniprot_id"]) & (pdb_features["structural_segment"] == r["structural_segment"])]["median_z"]
        assert np.isclose(r["median_of_entry_median_z"], vals.median())


def test_feature_matrix_columns_correspond_to_declared_structural_segments(receptor_profile):
    _require(RECEPTOR_MATRIX_MEDIAN)
    matrix = pd.read_parquet(RECEPTOR_MATRIX_MEDIAN)
    cols = [c for c in matrix.columns if c != "uniprot_id"]
    assert set(cols) == set(receptor_profile["structural_segment"].unique())


def test_coverage_never_mixed_into_flexibility_columns(pdb_features):
    for c in ["mean_z", "median_z", "sd_z", "q25_z", "q75_z", "iqr_z"]:
        assert pdb_features[c].dropna().abs().max() < 100  # sanity: these are z-scores, not fractions in [0,1] mislabeled
    assert pdb_features["segment_coverage"].dropna().between(0, 1).all()


def test_segment_coverage_not_deleted_for_low_coverage_observations():
    _require(COVERAGE)
    cov = pd.read_parquet(COVERAGE)
    assert (cov["segment_coverage"] < 0.5).sum() > 0  # low-coverage rows retained, not filtered out


# ---------------------------------------------------------------- VDR

def test_vdr_present_in_normalization_and_segment_outputs(consensus, seg_map, receptor_profile):
    assert (consensus["uniprot_id"] == VDR).any()
    assert (seg_map["uniprot_id"] == VDR).any()
    assert (receptor_profile["uniprot_id"] == VDR).any()
    v = receptor_profile[receptor_profile["uniprot_id"] == VDR]
    assert (v["n_pdb_entries_total"] == 48).all()


def test_vdr_not_excluded_for_failing_90pct_coverage(master):
    v = master[master["uniprot_id"] == VDR]
    coord = v[v["has_stage3b_coordinate_data"]]
    assert len(coord) == 48
    assert (coord["coverage_ge_090"].astype("boolean") != True).all()  # noqa: E712 -- none pass 90%, all retained anyway
    assert not coord["previous_primary_qc_member"].any()


# ---------------------------------------------------------------- determinism / idempotency

def test_all_generated_files_deterministic_and_idempotent(tmp_path):
    files = [CONSENSUS, SEGMENT_MAP, VALIDATION_CSV, COVERAGE, PDB_FEATURES, RECEPTOR_PROFILE,
             RECEPTOR_MATRIX_MEDIAN, RECEPTOR_MATRIX_MEAN, COVERAGE_MATRIX, ALIGNMENT_MAP]
    for f in files:
        _require(f)
    before = {f: _sha(f) for f in files}
    for script in ["34_build_dssp_consensus.py", "35_build_structural_segment_map.py",
                   "36_build_segment_features.py", "37_build_segment_missingness_and_alignment.py"]:
        r = subprocess.run([sys.executable, f"scripts/{script}"], cwd=PROJECT_ROOT, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
    after = {f: _sha(f) for f in files}
    assert before == after
