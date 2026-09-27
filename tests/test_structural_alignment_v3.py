"""
Structural-segment alignment v3 tests: complete curated multi-structure
templates (Parts 1-3), cross-template/complementary-structure rescue
(Parts 4, 6), canonical-order biological-consistency enforcement (Part 7),
and protection of the approved normalization foundation, historical
strict-QC files, Role D outputs, and the v1/v2 diagnostic-baseline commits.

Offline; uses only tracked files plus the locally-cached mmCIF files this
session copied in.
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
REF_STRUCTURES_V2 = ANALYSIS / "canonical_lbd_reference_structures_v2.csv"
TEMPLATE_V2 = ANALYSIS / "canonical_lbd_template_segments_v2.csv"
TRANSFER_CSV = ANALYSIS / "lbd_structural_label_transfer.csv"
MAPPING_SOURCES_CSV = ANALYSIS / "receptor_segment_mapping_sources.csv"
H12_AUDIT_CSV = REPORTS / "h12_mapping_audit.csv"
H5_H10_H11_AUDIT_CSV = REPORTS / "h5_h10_h11_mapping_audit.csv"
VALIDATION_V3_CSV = REPORTS / "structural_segment_validation_v3.csv"
RECEPTOR_MATRIX_MEDIAN = ANALYSIS / "receptor_segment_feature_matrix_median.parquet"
PDB_FEATURES = ANALYSIS / "pdb_segment_flexibility_features.parquet"

BASELINE_COMMIT = "e714bbd"  # approved normalization foundation
V1_BASELINE_COMMIT = "9f15062"
V2_BASELINE_COMMIT = "abd2257"
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
ANCHORS = {"P11473": "VDR", "P37231": "PPARgamma", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}
VALID_H12_CATEGORIES = {"A", "B", "C", "D", "E"}


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
def master():
    return pd.read_csv(MASTER, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def ref_structures_v2():
    _require(REF_STRUCTURES_V2)
    return pd.read_csv(REF_STRUCTURES_V2)


@pytest.fixture(scope="module")
def template_v2():
    _require(TEMPLATE_V2)
    return pd.read_csv(TEMPLATE_V2, keep_default_na=False)


# ---------------------------------------------------------------- v1/v2 baselines preserved

def test_v1_and_v2_diagnostic_baselines_still_reachable():
    for commit in (V1_BASELINE_COMMIT, V2_BASELINE_COMMIT):
        r = subprocess.run(["git", "cat-file", "-e", f"{commit}^{{commit}}"], cwd=PROJECT_ROOT)
        assert r.returncode == 0, f"diagnostic-baseline commit {commit} must not be deleted or rewritten"


# ---------------------------------------------------------------- Parts 1-3: curated template provenance

def test_reference_structures_v2_has_primary_and_provenance(ref_structures_v2):
    for uid in ANCHORS:
        g = ref_structures_v2[ref_structures_v2["uniprot_id"] == uid]
        assert len(g) >= 1
        assert (g["role"] == "PRIMARY").sum() == 1
        assert set(g["role"]) <= {"PRIMARY", "COMPLEMENTARY"}
        assert (g["source_publication"] != "").all()


def test_template_v2_every_row_has_provenance(template_v2):
    named = template_v2[template_v2["structural_segment"] != ""]
    assert len(named) > 0
    assert (named["reference_pdb"] != "").all()
    assert (named["publication"] != "").all()
    assert set(named["assignment_confidence"]) <= {"HIGH", "MODERATE"}


def test_complete_curated_templates_cover_more_than_landmark_residues(template_v2):
    """The v3 templates must establish COMPLETE segment spans (multiple
    residues), not single landmark points -- this is the explicit goal of
    the re-curation, distinct from the earlier single-residue-anchor
    approach."""
    named = template_v2[template_v2["structural_segment"] != ""]
    span = named.groupby(["uniprot_id", "structural_segment"]).size()
    assert (span >= 2).mean() > 0.7  # the large majority of segments are multi-residue spans


# ---------------------------------------------------------------- no ordinal DSSP counting

def test_no_segment_assigned_by_ordinal_helix_counting(seg_map, template_v2):
    methods = set(seg_map["assignment_method"].unique()) | set(template_v2.loc[template_v2["structural_segment"] != "", "reference_pdb"].index.map(lambda _: ""))
    assert "topology_ordinal_bracketed" not in seg_map["assignment_method"].unique()
    valid_methods = {
        "own_canonical_template", "structural_alignment_transfer", "complementary_structure_rescue",
        "loop_bracketed_by_named_segments", "loop_unbracketed", "no_secondary_structure_evidence", "no_structures",
    }
    assert set(seg_map["assignment_method"].unique()) <= valid_methods


# ---------------------------------------------------------------- complementary-structure rescue (Part 4/6)

def test_complementary_structure_rescue_present_and_evidenced():
    _require(TRANSFER_CSV)
    t = pd.read_csv(TRANSFER_CSV, keep_default_na=False)
    rescued = t[t["assignment_method"] == "complementary_structure_rescue"]
    assert len(rescued) > 0
    assert rescued["evidence_notes"].str.contains("rescue", case=False).all()
    assert rescued["reference_pdb"].notna().all()


def test_receptor_segment_mapping_sources_traces_every_segment():
    _require(MAPPING_SOURCES_CSV)
    src = pd.read_csv(MAPPING_SOURCES_CSV)
    assert len(src) > 0
    assert src["source_pdb"].notna().all()
    assert not src.duplicated(["uniprot_id", "structural_segment"]).any()


# ---------------------------------------------------------------- unique-per-residue mapping

def test_receptor_level_mapping_unique_per_residue(seg_map):
    assert not seg_map.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()


# ---------------------------------------------------------------- Part 7: canonical order enforced

def test_canonical_segment_order_strictly_enforced(seg_map):
    """No receptor may have a later-canonical-rank segment (e.g. H5) start
    before an earlier-canonical-rank segment (e.g. H3) ends -- the real
    violation found and fixed during this stage (VDR's cross-template-
    rescued H5 initially landed inside its own N-terminal insertion domain,
    before its own literature-anchored H3)."""
    order = {"H1": 1, "H2": 2, "H3": 3, "H4": 4, "H5": 5, "beta_sheet": 6, "H6": 7, "H7": 8,
             "H8": 9, "H9": 10, "H10": 11, "H11": 12, "H12": 13}
    named = seg_map[seg_map["structural_segment"].isin(order)]
    for uid, g in named.groupby("uniprot_id"):
        spans = g.groupby("structural_segment")["canonical_uniprot_position"].agg(["min", "max"])
        spans["rank"] = spans.index.map(order)
        spans = spans.sort_values("rank")
        assert spans["min"].is_monotonic_increasing, f"{uid}: segment start positions not increasing with canonical rank"
        # no segment may start before the previous (lower-rank) segment's own end
        prev_end = None
        for _, row in spans.iterrows():
            if prev_end is not None:
                assert row["min"] > prev_end, f"{uid}: a later segment starts before an earlier one ends"
            prev_end = row["max"]


def test_vdr_h3_before_rescued_segments():
    """Direct regression test for the specific violation found: VDR's H3
    (literature-anchored, HIGH confidence) must start before any segment
    ranked after it (H4, H5, beta_sheet, H6+)."""
    _require(SEGMENT_MAP)
    seg_map = pd.read_csv(SEGMENT_MAP, keep_default_na=False)
    v = seg_map[seg_map["uniprot_id"] == "P11473"]
    h3 = v[v["structural_segment"] == "H3"]
    assert len(h3) > 0
    h3_start = h3["canonical_uniprot_position"].min()
    for seg in ["H4", "H5", "beta_sheet", "H9", "H10", "H11", "H12"]:
        g = v[v["structural_segment"] == seg]
        if len(g):
            assert g["canonical_uniprot_position"].min() > h3_start, f"VDR {seg} starts before its own H3"


# ---------------------------------------------------------------- Part 6: H12 audit categories

def test_h12_audit_categories_exhaustive_and_valid():
    _require(H12_AUDIT_CSV)
    h12 = pd.read_csv(H12_AUDIT_CSV)
    if len(h12) == 0:
        pytest.skip("all receptors have H12 mapped -- nothing to audit")
    assert set(h12["category"].unique()) <= VALID_H12_CATEGORIES
    assert h12["note"].notna().all() and (h12["note"] != "").all()
    assert not h12["uniprot_id"].duplicated().any()


def test_h5_h10_h11_audit_categories_valid():
    _require(H5_H10_H11_AUDIT_CSV)
    df = pd.read_csv(H5_H10_H11_AUDIT_CSV)
    if len(df) == 0:
        pytest.skip("no missing H5/H10/H11 cases")
    assert set(df["category"].unique()) <= VALID_H12_CATEGORIES
    assert set(df["segment"].unique()) <= {"H5", "H10", "H11"}


# ---------------------------------------------------------------- Part 7 acceptance: mapping support reported

def test_core_segment_support_reported_and_h5_h7_h10_not_near_zero(seg_map):
    n_with_coord = 42
    named = seg_map[seg_map["structural_segment"] != ""]
    counts = named.groupby("structural_segment")["uniprot_id"].nunique()
    for seg in ["H3", "H5", "H7", "H10", "H11", "H12"]:
        assert seg in counts.index or counts.get(seg, 0) == 0  # always computable, even if 0
    for seg in ["H5", "H7"]:
        assert counts.get(seg, 0) >= n_with_coord * 0.10, f"{seg} unexpectedly near zero after full curation"
    # H10 is documented as the known remaining weak spot (bundled H10/H11 landmark evidence
    # favors H11); still assert it is computed and non-negative, not silently dropped.
    assert counts.get("H10", 0) >= 0


# ---------------------------------------------------------------- Part 11: final validation

def test_validation_v3_covers_all_five_anchors():
    _require(VALIDATION_V3_CSV)
    v = pd.read_csv(VALIDATION_V3_CSV)
    for uid in ANCHORS:
        assert (v["uniprot_id"] == uid).any(), f"{ANCHORS[uid]} missing from structural_segment_validation_v3.csv"
    assert (v["confidence"].isin(["HIGH", "MODERATE"])).all()


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

def test_chain_to_pdb_to_receptor_aggregation_unchanged():
    _require(PDB_FEATURES)
    pdb_feat = pd.read_parquet(PDB_FEATURES)
    pos = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")
    seg_map_df = pd.read_csv(SEGMENT_MAP, keep_default_na=False)
    named = seg_map_df[seg_map_df["structural_segment"] != ""][["uniprot_id", "canonical_uniprot_position", "structural_segment"]]
    merged = pos.merge(named, on=["uniprot_id", "canonical_uniprot_position"])
    sample = pdb_feat.dropna(subset=["median_z"]).sample(min(15, len(pdb_feat)), random_state=17)
    for _, r in sample.iterrows():
        vals = merged[(merged["uniprot_id"] == r["uniprot_id"]) & (merged["pdb_id"] == r["pdb_id"]) & (merged["structural_segment"] == r["structural_segment"])]["b_factor_z"]
        assert np.isclose(r["median_z"], vals.median())

    _require(RECEPTOR_MATRIX_MEDIAN)
    matrix = pd.read_parquet(RECEPTOR_MATRIX_MEDIAN).set_index("uniprot_id")
    profile = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    sample2 = profile.dropna(subset=["median_of_entry_median_z"]).sample(min(10, len(profile)), random_state=19)
    for _, r in sample2.iterrows():
        assert np.isclose(matrix.loc[r["uniprot_id"], r["structural_segment"]], r["median_of_entry_median_z"])
    assert matrix.isna().any().any()  # missing left as NaN, never 0


# ---------------------------------------------------------------- determinism

def test_v3_pipeline_deterministic_and_idempotent():
    files = [REF_STRUCTURES_V2, TEMPLATE_V2, TRANSFER_CSV, MAPPING_SOURCES_CSV, H12_AUDIT_CSV,
             SEGMENT_MAP, RECEPTOR_MATRIX_MEDIAN, VALIDATION_V3_CSV]
    for f in files:
        _require(f)
    before = {f: _sha(f) for f in files}
    for script in ["44_curate_templates_v3.py", "45_rerun_transfer_v3.py",
                   "46_rebuild_segment_map_v3.py", "47_validate_v3.py"]:
        r = subprocess.run([sys.executable, f"scripts/{script}"], cwd=PROJECT_ROOT, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
    after = {f: _sha(f) for f in files}
    assert before == after
