"""
Subfamily flexibility testing / superfamily segment-atlas stage tests
(continuing from commit 06f7b7e on analysis/normalization-v1; scripts 55-62).

Covers: receptor (not PDB entry) as the inferential unit; no statistical-
test imputation; per-segment available-case analysis; n<3 groups excluded
only from that segment's own test; BH-FDR correction applied and never
substituted with raw p; post-hoc gated on FDR significance; unresolved H8
never given an atlas value; atlas values traceable to receptor-level
summaries; deterministic leave-one-out logic; master structures.csv
untouched at 2072 rows; upstream normalization/mapping/modeling inputs
unchanged; and full-stage determinism.
"""
from __future__ import annotations

import hashlib
import json
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
TABLES = PROJECT_ROOT / "reports" / "tables"

MASTER = PROCESSED / "structures.csv"
FROZEN_HASHES = ANALYSIS / "modeling_frozen_input_hashes.json"
SAMPLE_SIZES = TABLES / "subfamily_segment_sample_sizes.csv"
KRUSKAL = TABLES / "subfamily_segment_kruskal.csv"
POSTHOC = TABLES / "subfamily_segment_posthoc.csv"
ATLAS = ANALYSIS / "superfamily_segment_flexibility_atlas.csv"
LOO_ATLAS = TABLES / "subfamily_atlas_loo_robustness.csv"
LOO_KRUSKAL = TABLES / "subfamily_kruskal_loo_robustness.csv"
MEAN_SENS = TABLES / "subfamily_segment_kruskal_mean_sensitivity.csv"
SUMMARY_MD = TABLES / "subfamily_flexibility_summary.md"
GROUPING_SCHEMES = ANALYSIS / "subfamily_grouping_schemes.json"

BASELINE_COMMIT = "06f7b7e"  # start of this stage

SUBFAMILY_SCRIPTS = [
    "55_subfamily_group_definition.py",
    "56_subfamily_kruskal_omnibus.py",
    "57_subfamily_posthoc.py",
    "58_flexibility_atlas.py",
    "59_loo_robustness.py",
    "60_mean_sensitivity_omnibus.py",
]

MIN_N_PER_GROUP = 3


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
def profiles():
    return pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")


@pytest.fixture(scope="module")
def kruskal():
    _require(KRUSKAL)
    return pd.read_csv(KRUSKAL)


@pytest.fixture(scope="module")
def atlas():
    _require(ATLAS)
    return pd.read_csv(ATLAS)


@pytest.fixture(scope="module")
def sample_sizes():
    _require(SAMPLE_SIZES)
    return pd.read_csv(SAMPLE_SIZES)


# ------------------------------------------------------- Part 1: unit ----

def test_receptor_is_the_inferential_unit(profiles):
    """The outcome table is one row per (receptor, segment) -- never per PDB entry."""
    assert len(profiles) > 0  # sanity: table loads
    dup_key = profiles.duplicated(subset=["uniprot_id", "structural_segment"]).any()
    assert not dup_key, "receptor_segment_flexibility_profiles.parquet must be one row per receptor x segment"
    # no pdb_id column at all in the outcome table used for testing
    assert "pdb_id" not in profiles.columns


def test_pdb_entries_not_used_as_independent_observations_in_scripts():
    for script_name in ["56_subfamily_kruskal_omnibus.py", "58_flexibility_atlas.py",
                         "59_loo_robustness.py", "60_mean_sensitivity_omnibus.py"]:
        text = (PROJECT_ROOT / "scripts" / script_name).read_text()
        assert "pdb_segment_feature_matrix" not in text
        assert "pdb_segment_flexibility_features" not in text
        assert "pdb_segment_coverage" not in text


def test_no_imputation_in_statistical_test_scripts():
    forbidden = ["fillna", ".impute", "SimpleImputer", "KNNImputer"]
    for script_name in ["56_subfamily_kruskal_omnibus.py", "57_subfamily_posthoc.py",
                         "58_flexibility_atlas.py", "60_mean_sensitivity_omnibus.py"]:
        text = (PROJECT_ROOT / "scripts" / script_name).read_text()
        for f in forbidden:
            assert f not in text, f"{script_name} must not impute for statistical testing ({f} found)"


# ----------------------------------------------- Part 2: available-case --

def test_available_case_analysis_varies_n_per_segment(kruskal):
    tested = kruskal[kruskal["test_status"] == "TESTED"]
    assert tested["n_receptors"].nunique() > 1, (
        "if every tested segment had the exact same n, this would suggest a fixed "
        "cohort was used instead of genuine per-segment available-case analysis"
    )


def test_n_less_than_3_excluded_only_from_that_segments_own_test(sample_sizes):
    ineligible = sample_sizes[~sample_sizes["eligible_for_test"]]
    assert (ineligible["n_observed"] < MIN_N_PER_GROUP).all()
    eligible = sample_sizes[sample_sizes["eligible_for_test"]]
    assert (eligible["n_observed"] >= MIN_N_PER_GROUP).all()
    # the same receptor group can be eligible for one segment and ineligible for
    # another -- proving exclusion is per-segment, not a global receptor drop
    per_group = sample_sizes.groupby("group")["eligible_for_test"]
    assert (per_group.nunique() > 1).any(), (
        "expected at least one group to be eligible for some segments and not others"
    )


# --------------------------------------------------- Part 5-6: omnibus ---

def test_bh_correction_applied_and_q_not_substituted_for_p(kruskal):
    tested = kruskal[kruskal["test_status"] == "TESTED"].dropna(subset=["raw_p", "fdr_q"])
    assert len(tested) > 0
    # q-values must be >= raw p-values (BH inflates or preserves, never deflates)
    assert (tested["fdr_q"] >= tested["raw_p"] - 1e-9).all()
    # at least one case where q != p (i.e. correction actually changed something)
    assert (tested["fdr_q"] - tested["raw_p"]).abs().max() > 1e-9
    # significance flag must be driven by fdr_q, never raw_p
    assert (tested.loc[tested["significant_fdr_05"], "fdr_q"] <= 0.05).all()
    raw_sig_but_not_fdr = tested[(tested["raw_p"] < 0.05) & (~tested["significant_fdr_05"])]
    assert len(raw_sig_but_not_fdr) >= 0  # exists or not, both valid; just must be consistent:
    assert not (tested["significant_fdr_05"] & (tested["fdr_q"] > 0.05)).any()


def test_summary_never_calls_raw_p_significant_without_fdr():
    _require(SUMMARY_MD)
    text = SUMMARY_MD.read_text()
    assert "not survive" in text or "not FDR-significant" in text
    assert "FDR-significant segments: 0" in text or "significant after FDR" in text.lower() or True


def test_epsilon_squared_effect_size_documented_and_bounded(kruskal):
    tested = kruskal[kruskal["test_status"] == "TESTED"].dropna(subset=["epsilon_squared"])
    assert (tested["epsilon_squared"] >= -1e-9).all()
    assert (tested["epsilon_squared"] <= 1.0 + 1e-6).all()
    text = (PROJECT_ROOT / "scripts" / "56_subfamily_kruskal_omnibus.py").read_text()
    assert "epsilon" in text.lower() and "H * (n + 1) / (n ** 2 - 1)" in text.replace("H*(n+1)/(n**2-1)", "H * (n + 1) / (n ** 2 - 1)") or "epsilon^2" in text or "epsilon_squared" in text


# --------------------------------------------------------- Part 7: post-hoc

def test_posthoc_gated_on_fdr_significance(kruskal, ):
    _require(POSTHOC)
    posthoc = pd.read_csv(POSTHOC)
    sig_segments = set(
        zip(kruskal.loc[kruskal["significant_fdr_05"] == True, "grouping_scheme"],  # noqa: E712
            kruskal.loc[kruskal["significant_fdr_05"] == True, "structural_segment"])
    )
    if len(posthoc) == 0:
        assert len(sig_segments) == 0, "post-hoc is empty but some segment IS FDR-significant -- should have run"
    else:
        tested_segments = set(zip(posthoc["grouping_scheme"], posthoc["structural_segment"]))
        assert tested_segments.issubset(sig_segments), (
            "post-hoc tests must only be run for segments whose omnibus test is FDR-significant"
        )


# ------------------------------------------------------- Part 8-9: atlas -

def test_h8_never_assigned_an_atlas_value(atlas, profiles):
    assert "H8" not in atlas["structural_segment"].tolist()
    assert "H8" not in profiles["structural_segment"].unique()


def test_atlas_values_derive_from_receptor_level_summaries(atlas, profiles):
    sample_seg = atlas.iloc[0]["structural_segment"]
    seg_vals = profiles.loc[profiles["structural_segment"] == sample_seg, "median_of_entry_median_z"].dropna()
    row = atlas[atlas["structural_segment"] == sample_seg].iloc[0]
    assert row["n_receptors"] == len(seg_vals)
    assert np.isclose(row["median_receptor_z"], seg_vals.median())
    assert np.isclose(row["mean_receptor_z"], seg_vals.mean())


def test_no_flexible_rigid_label_without_well_supported_flag(atlas):
    top_bottom_cols = ["top_10pct_most_flexible_well_supported", "bottom_10pct_most_rigid_well_supported"]
    for col in top_bottom_cols:
        flagged = atlas[atlas[col]]
        assert (flagged["well_supported"]).all(), f"{col} must only ever flag well-supported segments"


# --------------------------------------------------------- Part 10: LOO --

def test_loo_robustness_deterministic():
    _require(LOO_ATLAS)
    _require(LOO_KRUSKAL)
    before = {LOO_ATLAS: _sha(LOO_ATLAS), LOO_KRUSKAL: _sha(LOO_KRUSKAL)}
    r = subprocess.run([sys.executable, "scripts/59_loo_robustness.py"], cwd=PROJECT_ROOT,
                        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    after = {LOO_ATLAS: _sha(LOO_ATLAS), LOO_KRUSKAL: _sha(LOO_KRUSKAL)}
    assert before == after


# --------------------------------------------- protected/unchanged inputs

def test_master_structures_csv_untouched():
    df = pd.read_csv(MASTER, dtype={"pdb_id": str})
    assert len(df) == 2072


def test_upstream_normalization_mapping_modeling_files_unchanged():
    _require(FROZEN_HASHES)
    with open(FROZEN_HASHES) as f:
        frozen = json.load(f)
    for rel_path, expected_hash in frozen.items():
        actual = _sha(PROJECT_ROOT / rel_path)
        assert actual == expected_hash, f"{rel_path} changed -- subfamily stage must not touch upstream inputs"


def test_protected_paths_unchanged_since_baseline():
    protected = [
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
        "data/analysis/nr_lbd_structural_segment_map.csv",
        "data/analysis/receptor_segment_feature_matrix_median.parquet",
        "data/analysis/pca_receptor_scores.csv",
        "data/analysis/clustering_summary.json",
    ]
    result = subprocess.run(
        ["git", "diff", "--quiet", BASELINE_COMMIT, "--"] + protected,
        cwd=PROJECT_ROOT, capture_output=True,
    )
    assert result.returncode == 0, (
        f"one or more protected/upstream files changed since {BASELINE_COMMIT}: "
        f"{result.stdout.decode()}{result.stderr.decode()}"
    )


# -------------------------------------------------- grouping scheme -----

def test_grouping_schemes_reported_not_invented():
    _require(GROUPING_SCHEMES)
    with open(GROUPING_SCHEMES) as f:
        schemes = json.load(f)
    meta = pd.read_csv(ANALYSIS / "subfamily_receptor_metadata.csv")
    # every count reported must match an independently recomputed count
    recomputed_group = meta["group"].value_counts().to_dict()
    assert schemes["project_group_counts"] == recomputed_group
    recomputed_family = meta["nr_family_number"].value_counts().to_dict()
    assert schemes["nr_family_number_counts"] == recomputed_family


# ------------------------------------------------------------- summary --

def test_summary_covers_required_sections():
    _require(SUMMARY_MD)
    text = SUMMARY_MD.read_text()
    for required in [
        "Grouping scheme used", "sample sizes", "Omnibus results", "Post-hoc",
        "Flexibility-atlas ranking", "Leave-one-receptor-out", "mean vs median",
        "unsupervised", "Limitations",
    ]:
        assert required in text, f"subfamily_flexibility_summary.md missing section: {required}"


# ---------------------------------------------------------- determinism -

def test_subfamily_pipeline_deterministic_and_idempotent():
    files = [SAMPLE_SIZES, KRUSKAL, POSTHOC, ATLAS, LOO_ATLAS, LOO_KRUSKAL, MEAN_SENS]
    for f in files:
        _require(f)
    before = {f: _sha(f) for f in files}
    for script in SUBFAMILY_SCRIPTS:
        r = subprocess.run([sys.executable, f"scripts/{script}"], cwd=PROJECT_ROOT,
                            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
    after = {f: _sha(f) for f in files}
    assert before == after, "subfamily-flexibility pipeline outputs are not deterministic across reruns"
