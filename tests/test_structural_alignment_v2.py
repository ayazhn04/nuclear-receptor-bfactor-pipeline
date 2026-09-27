"""
Structural-segment alignment v2 tests: real 3D structural-alignment-based
label transfer (not sequential DSSP-helix counting), anchor-receptor
agreement with curated literature/repository references, and protection of
the approved normalization foundation, historical strict-QC files, Role D
outputs, and the v1 diagnostic-baseline commit.

Offline; uses only tracked files plus the locally-cached mmCIF files this
session copied in (same as test_structural_segment_alignment.py).
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
SEGMENT_MAP = ANALYSIS / "nr_lbd_structural_segment_map.csv"
TEMPLATE_CSV = ANALYSIS / "canonical_lbd_template_segments.csv"
SCORES_CSV = ANALYSIS / "lbd_structural_alignment_scores.csv"
TRANSFER_CSV = ANALYSIS / "lbd_structural_label_transfer.csv"
REF_CSV = ANALYSIS / "lbd_structural_reference_structures.csv"
VALIDATION_V2_CSV = REPORTS / "structural_segment_validation_v2.csv"
RECEPTOR_MATRIX_MEDIAN = ANALYSIS / "receptor_segment_feature_matrix_median.parquet"
PDB_FEATURES = ANALYSIS / "pdb_segment_flexibility_features.parquet"

BASELINE_COMMIT = "e714bbd"  # approved normalization foundation
V1_BASELINE_COMMIT = "9f15062"  # v1 diagnostic baseline, must remain intact in history
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
def seg_map():
    _require(SEGMENT_MAP)
    return pd.read_csv(SEGMENT_MAP, keep_default_na=False)


@pytest.fixture(scope="module")
def scores():
    _require(SCORES_CSV)
    return pd.read_csv(SCORES_CSV)


@pytest.fixture(scope="module")
def transfer():
    _require(TRANSFER_CSV)
    return pd.read_csv(TRANSFER_CSV, keep_default_na=False)


@pytest.fixture(scope="module")
def master():
    return pd.read_csv(MASTER, dtype={"pdb_id": str})


# ---------------------------------------------------------------- v1 baseline preserved

def test_v1_diagnostic_baseline_commit_still_reachable():
    r = subprocess.run(["git", "cat-file", "-e", f"{V1_BASELINE_COMMIT}^{{commit}}"], cwd=PROJECT_ROOT)
    assert r.returncode == 0, "v1 diagnostic-baseline commit 9f15062 must not be deleted or rewritten"
    r2 = subprocess.run(["git", "show", f"{V1_BASELINE_COMMIT}:data/analysis/nr_lbd_structural_segment_map.csv"],
                         cwd=PROJECT_ROOT, capture_output=True)
    assert r2.returncode == 0 and len(r2.stdout) > 0


# ---------------------------------------------------------------- labels not from ordinal DSSP counting

def test_canonical_labels_not_primarily_from_sequential_dssp_counting(seg_map):
    named = seg_map[seg_map["structural_segment"] != ""]
    methods = named["assignment_method"].value_counts()
    # v1's method name was "topology_ordinal_bracketed"; v2's primary mechanisms are
    # explicit literature anchors on the template and CE structural-alignment transfer.
    ordinal_only = methods.get("topology_ordinal_bracketed", 0)
    real_evidence = methods.get("own_canonical_template", 0) + methods.get("structural_alignment_transfer", 0) + methods.get("loop_bracketed_by_named_segments", 0)
    assert real_evidence > 0
    assert ordinal_only == 0, "v2 must not use the v1 sequential-DSSP-ordinal method at all"


def test_structural_alignment_scores_are_real_residue_correspondences(scores, transfer):
    assert {"n_aligned", "z_score", "rmsd", "coverage", "sequence_identity"} <= set(scores.columns)
    usable = scores[scores["usable"]]
    assert len(usable) > 0
    assert (usable["n_aligned"] > 0).all()
    assert (usable["rmsd"] >= 0).all()
    assert usable["sequence_identity"].between(0, 1).all()
    # the transfer table must cite an actual alignment score for every non-anchor row
    non_anchor = transfer[transfer["assignment_method"] == "structural_alignment_transfer"]
    assert len(non_anchor) > 0
    assert pd.to_numeric(non_anchor["structural_alignment_score"], errors="coerce").notna().all()


def test_ce_alignment_produces_monotonic_residue_correspondence():
    sys.path.insert(0, str(PROJECT_ROOT))
    from scripts.utils.structural_alignment import align_structures

    # two simple, clearly homologous fabricated helical tracks (a straight line of points,
    # offset) -- just exercises the real CE call end to end and checks the monotonicity
    # invariant the label-transfer step depends on.
    a = [[float(i), 0.0, 0.0] for i in range(20)]
    b = [[float(i), 0.5, 0.0] for i in range(20)]
    result = align_structures(a, b, window_size=4, max_gap=5)
    assert result is not None
    assert result.n_aligned > 0
    assert all(x < y for x, y in zip(result.idx_a, result.idx_a[1:]))
    assert all(x < y for x, y in zip(result.idx_b, result.idx_b[1:]))


# ---------------------------------------------------------------- VDR H12 / anchor agreement

def test_vdr_h12_mapping_reproduced(seg_map):
    v = seg_map[(seg_map["uniprot_id"] == VDR) & (seg_map["structural_segment"] == "H12")]
    assert len(v) > 0
    # repository-documented VDR terminal helix (nr_functional_features.csv h12/AF-2): 416-422
    assert v["canonical_uniprot_position"].min() <= 420
    assert v["canonical_uniprot_position"].max() >= 418


def test_anchor_receptors_agree_with_curated_charge_clamp_landmarks(seg_map):
    _require(VALIDATION_V2_CSV)
    fsite = pd.read_csv(REPORTS / "functional_site_definitions.csv")
    charge_clamp = fsite[fsite["site_id"] == "af2_coactivator_groove"]
    checked = 0
    for uid in ANCHORS:
        cc = charge_clamp[charge_clamp["uniprot_id"] == uid]
        for _, r in cc.iterrows():
            row = seg_map[(seg_map["uniprot_id"] == uid) & (seg_map["canonical_uniprot_position"] == r["uniprot_position"])]
            if len(row) == 0:
                continue
            seg = row["structural_segment"].iloc[0]
            if seg == "":
                continue
            assert seg in ("H3", "H12"), f"{uid} charge-clamp residue {r['uniprot_position']} mapped to {seg}, expected H3 or H12"
            checked += 1
    assert checked >= 5


def test_five_anchors_present_in_validation_v2():
    _require(VALIDATION_V2_CSV)
    v = pd.read_csv(VALIDATION_V2_CSV)
    for uid in ANCHORS:
        assert (v["uniprot_id"] == uid).any(), f"{ANCHORS[uid]} missing from structural_segment_validation_v2.csv"


# ---------------------------------------------------------------- structural sanity

def test_segment_boundaries_map_to_valid_uniprot_positions(seg_map):
    lbd_ref = pd.read_csv(MANIFESTS / "nr_lbd_reference.csv").set_index("uniprot_id")
    named = seg_map[seg_map["structural_segment"] != ""]
    for uid, g in named.groupby("uniprot_id"):
        lo, hi = lbd_ref.loc[uid, ["lbd_start", "lbd_end"]]
        assert g["canonical_uniprot_position"].between(lo, hi).all()


def test_one_uniprot_residue_has_at_most_one_final_segment(seg_map):
    assert not seg_map.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()


def test_canonical_segment_ordering_valid(seg_map):
    named = seg_map[seg_map["structural_segment"].str.fullmatch(r"H\d{1,2}")]
    for uid, g in named.groupby("uniprot_id"):
        g = g.drop_duplicates("structural_segment").assign(n=lambda d: d["structural_segment"].str[1:].astype(int)).sort_values("n")
        starts = g.groupby("structural_segment")["canonical_uniprot_position"].min().reindex(g["structural_segment"])
        assert starts.is_monotonic_increasing, f"{uid}: named helices out of N->C order"


def test_low_confidence_assignments_remain_explicit(seg_map):
    ambiguous = seg_map[seg_map["assignment_status"].isin(["AMBIGUOUS", "PENDING_EXTERNAL_VALIDATION"])]
    unassigned = seg_map[seg_map["assignment_status"] == "UNASSIGNED"]
    assert len(ambiguous) > 0 or len(unassigned) > 0
    # A low-confidence AMBIGUOUS row may still carry a segment name transferred by
    # structural alignment (Part 6's LOW/AMBIGUOUS tier: "conflicting mappings or
    # insufficient coverage" -- we know which template segment it points to, just not
    # confidently) -- but it must always be evidenced by that specific method, never guessed.
    amb = seg_map[seg_map["assignment_status"] == "AMBIGUOUS"]
    named_amb = amb[amb["structural_segment"] != ""]
    assert set(named_amb["assignment_method"].unique()) <= {"structural_alignment_transfer"}


# ---------------------------------------------------------------- protection

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


def test_master_structures_2072_rows_no_candidate_excluded(master):
    assert len(master) == 2072
    assert int(master["has_stage3b_coordinate_data"].sum()) == 1840


# ---------------------------------------------------------------- aggregation hierarchy unchanged

def test_chain_to_pdb_aggregation_still_median():
    _require(PDB_FEATURES)
    pdb_feat = pd.read_parquet(PDB_FEATURES)
    pos = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")
    seg_map_df = pd.read_csv(SEGMENT_MAP, keep_default_na=False)
    named = seg_map_df[seg_map_df["structural_segment"] != ""][["uniprot_id", "canonical_uniprot_position", "structural_segment"]]
    merged = pos.merge(named, on=["uniprot_id", "canonical_uniprot_position"])
    sample = pdb_feat.dropna(subset=["median_z"]).sample(min(15, len(pdb_feat)), random_state=11)
    for _, r in sample.iterrows():
        vals = merged[(merged["uniprot_id"] == r["uniprot_id"]) & (merged["pdb_id"] == r["pdb_id"]) & (merged["structural_segment"] == r["structural_segment"])]["b_factor_z"]
        assert np.isclose(r["median_z"], vals.median())


def test_receptor_matrix_uses_median_primary():
    _require(RECEPTOR_MATRIX_MEDIAN)
    matrix = pd.read_parquet(RECEPTOR_MATRIX_MEDIAN).set_index("uniprot_id")
    assert matrix.isna().any().any()  # missing left as NaN, never 0
    profile = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    sample = profile.dropna(subset=["median_of_entry_median_z"]).sample(min(10, len(profile)), random_state=13)
    for _, r in sample.iterrows():
        assert np.isclose(matrix.loc[r["uniprot_id"], r["structural_segment"]], r["median_of_entry_median_z"])


# ---------------------------------------------------------------- determinism

def test_v2_reference_structures_deterministic_and_idempotent():
    """Only script 39 (reference-structure selection, Part 2 -- unchanged
    input/logic between v2 and v3) is re-run here. Scripts 40-43 (v2's
    single-landmark templates/transfer/segment-map/feature builders) are
    NOT re-run: they write to the same output filenames v3's scripts 44-47
    now own (nr_lbd_structural_segment_map.csv, lbd_structural_alignment_scores.csv,
    lbd_structural_label_transfer.csv, and the downstream feature/matrix
    tables), and re-running the superseded v2 builders would silently
    overwrite the live v3 pipeline's output with v2 content. v3's own
    idempotency (scripts 44-47) is covered by
    test_structural_alignment_v3.py::test_v3_pipeline_deterministic_and_idempotent."""
    _require(REF_CSV)
    before = _sha(REF_CSV)
    r = subprocess.run([sys.executable, "scripts/39_select_structural_references.py"], cwd=PROJECT_ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    after = _sha(REF_CSV)
    assert before == after
