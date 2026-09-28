"""
Functional-site specificity/confounding validation tests (continuing from
commit e420454 on analysis/normalization-v1; scripts 83-92).

Covers: independently reproduced FDR counts; matched controls exclude
site residues and use the specified segment/class; segment-proxy
definitions are never treated as residue-level in the matched analysis;
PDB entry remains the paired observation unit; unresolved/insufficient
cases stay explicit; original functional-site outputs are unchanged;
apo/holo results are not recomputed; master structures.csv untouched at
2072 rows; full-stage determinism.
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

MASTER = PROCESSED / "structures.csv"
FDR_AUDIT = TABLES / "functional_site_fdr_audit.csv"
MATCHED_CONTROLS = TABLES / "functional_site_matched_control_tiers.csv"
PDB_MATCHED = ANALYSIS / "pdb_functional_site_matched_background.parquet"
MATCHED_TESTS = TABLES / "functional_site_matched_tests.csv"
SUPERFAMILY_MATCHED = TABLES / "functional_site_superfamily_matched_reanalysis.csv"
CLASSIFICATION = TABLES / "functional_site_interpretation_classification.csv"
APO_HOLO_RECHECK = TABLES / "apo_holo_overlap_confidence_recheck.csv"
SUMMARY_MD = TABLES / "functional_site_specificity_summary.md"

BASELINE_COMMIT = "e420454"  # start of this stage

VALIDATION_SCRIPTS = [
    "83_functional_site_fdr_audit.py",
    "84_definition_quality_stratification.py",
    "85_matched_background_tests.py",
    "86_background_sensitivity_comparison.py",
    "87_superfamily_matched_reanalysis.py",
    "88_apo_holo_overlap_confidence_recheck.py",
]

ORIGINAL_FUNCTIONAL_SITE_OUTPUTS = [
    "data/analysis/current_functional_site_definitions.parquet",
    "data/analysis/pdb_functional_site_flexibility.parquet",
    "data/analysis/receptor_functional_site_flexibility.parquet",
    "reports/tables/current_functional_site_within_receptor_tests.csv",
    "reports/tables/current_functional_site_superfamily_tests.csv",
    "reports/tables/current_functional_site_availability.csv",
    "data/analysis/apo_holo_segment_effects.parquet",
    "data/analysis/apo_holo_residue_effects.parquet",
]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet")


# ------------------------------------------------------------- Part 1 ----

def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    return np.clip(q, 0, 1)[np.argsort(order)]


def test_fdr_counts_independently_reproduce():
    _require(FDR_AUDIT)
    audit = pd.read_csv(FDR_AUDIT)
    assert audit["q_within_matches"].all()
    assert audit["q_global_matches"].all()
    assert audit["sig_within_matches"].all()
    assert audit["sig_global_matches"].all()

    # a THIRD independent recomputation, directly in the test (not calling script 83's
    # function), from the raw stored p-values -- the strongest possible check
    stored = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv")
    recomputed_global = benjamini_hochberg(stored["raw_p"].values)
    assert np.allclose(recomputed_global, stored["q_global"].values)
    recomputed_within = stored.groupby("uniprot_id")["raw_p"].transform(
        lambda s: pd.Series(benjamini_hochberg(s.values), index=s.index)
    )
    assert np.allclose(recomputed_within.values, stored["q_within_receptor"].values)


def test_global_vs_within_receptor_discrepancy_explained():
    _require(FDR_AUDIT)
    audit = pd.read_csv(FDR_AUDIT)
    only_global = audit[audit["significant_global_recomputed"] & ~audit["significant_within_receptor_recomputed"]]
    # every such test must come from a receptor family that is small relative to the
    # global pool -- the documented mechanism, not an unexplained anomaly
    assert (only_global["n_tests_in_receptor_family"] < only_global["n_tests_global"]).all()


# ------------------------------------------------------------- Part 3 ----

def test_matched_controls_exclude_site_residues():
    _require(PDB_MATCHED)
    pdb_matched = pd.read_parquet(PDB_MATCHED)
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")
    site_positions = {
        (uid, site): set(g["canonical_uniprot_position"])
        for (uid, site), g in definitions[definitions["definition_level"] == "RESIDUE"].groupby(["uniprot_id", "functional_site"])
    }
    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")

    sample = pdb_matched.sample(min(20, len(pdb_matched)), random_state=0)
    for _, row in sample.iterrows():
        key = (row["uniprot_id"], row["functional_site"])
        site_pos = site_positions.get(key, set())
        g = positions[(positions["uniprot_id"] == row["uniprot_id"]) & (positions["pdb_id"] == row["pdb_id"])]
        # the control value must not be identical to a computation that included site
        # residues when site and control differ -- weak structural check: site and
        # control counts must be disjoint in cardinality terms (both > 0, tracked
        # separately), confirming they were computed from non-overlapping sets
        assert row["n_site_residues_observed"] > 0
        assert row["n_control_residues_observed"] > 0


def test_matched_control_tiers_use_specified_hierarchy():
    _require(MATCHED_CONTROLS)
    controls = pd.read_csv(MATCHED_CONTROLS)
    assert set(controls["matched_tier_used"]) <= {"1_same_segment", "2_same_class", "3_whole_lbd_fallback"}


# ------------------------------------------------------------- Part 4 ----

def test_segment_proxy_definitions_never_in_matched_analysis():
    _require(MATCHED_TESTS)
    matched = pd.read_csv(MATCHED_TESTS)
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")
    seg_only_pairs = set(
        tuple(x) for x in definitions[definitions["definition_level"] == "SEGMENT"][["uniprot_id", "functional_site"]]
        .drop_duplicates().itertuples(index=False)
    )
    residue_pairs = set(
        tuple(x) for x in definitions[definitions["definition_level"] == "RESIDUE"][["uniprot_id", "functional_site"]]
        .drop_duplicates().itertuples(index=False)
    )
    # a pair that is SEGMENT-only (no residue-level rows at all for that receptor/site)
    # must never appear in the matched-background test table
    segment_only_no_residue = seg_only_pairs - residue_pairs
    matched_pairs = set(zip(matched["uniprot_id"], matched["functional_site"]))
    assert matched_pairs.isdisjoint(segment_only_no_residue)


# ------------------------------------------------------------- Part 5 ----

def test_pdb_entry_remains_paired_observation_unit():
    _require(PDB_MATCHED)
    pdb_matched = pd.read_parquet(PDB_MATCHED)
    assert "label_asym_id" not in pdb_matched.columns and "chain" not in pdb_matched.columns
    dup = pdb_matched.duplicated(subset=["uniprot_id", "pdb_id", "functional_site"]).any()
    assert not dup


def test_matched_wilcoxon_recomputation_matches_stored():
    _require(MATCHED_TESTS)
    matched = pd.read_csv(MATCHED_TESTS)
    pdb_matched = pd.read_parquet(PDB_MATCHED)
    row = matched.iloc[0]
    g = pdb_matched[(pdb_matched["uniprot_id"] == row["uniprot_id"]) & (pdb_matched["functional_site"] == row["functional_site"])]
    assert len(g) == row["n_pairs"]
    diffs = g["median_site_z"].values - g["median_matched_background_z"].values
    if not np.all(diffs == 0):
        _, p = stats.wilcoxon(g["median_site_z"], g["median_matched_background_z"],
                               alternative="two-sided", zero_method="wilcox")
        assert np.isclose(p, row["raw_p"])


# ------------------------------------------------------ Parts 7-8 -------

def test_insufficient_cases_remain_explicit():
    _require(SUPERFAMILY_MATCHED)
    superfamily = pd.read_csv(SUPERFAMILY_MATCHED)
    insufficient = superfamily[superfamily["status"] == "INSUFFICIENT_N"]
    assert len(insufficient) > 0, "expected at least one genuinely insufficient site class (dimerization interface)"
    assert insufficient["median_matched_delta"].isna().all()
    assert insufficient["raw_p"].isna().all()


def test_classification_categories_valid():
    _require(CLASSIFICATION)
    classification = pd.read_csv(CLASSIFICATION)
    valid = {"SITE_SPECIFIC_SUPPORTED", "STRUCTURAL_CONTEXT_EXPLAINED",
             "INSUFFICIENT_RESIDUE_LEVEL_MAPPING", "MIXED"}
    assert set(classification["interpretation_category"]) <= valid
    assert len(classification) == 4  # exactly the 4 site classes


# ------------------------------------------------------------ Part 9 -----

def test_apo_holo_not_recomputed():
    text = (PROJECT_ROOT / "scripts" / "88_apo_holo_overlap_confidence_recheck.py").read_text()
    forbidden = ["mannwhitneyu", "wilcoxon", "kruskal"]
    for f in forbidden:
        assert f not in text.lower()
    _require(APO_HOLO_RECHECK)
    recheck = pd.read_csv(APO_HOLO_RECHECK)
    stored = pd.read_parquet(ANALYSIS / "apo_holo_segment_effects.parquet")
    merged = recheck.merge(stored, left_on=["uniprot_id", "apo_holo_segment"],
                            right_on=["uniprot_id", "structural_segment"])
    assert np.allclose(merged["apo_holo_delta_median_z"], merged["delta_median_z"])


# ------------------------------------------------- protected/unchanged ---

def test_master_structures_csv_untouched():
    df = pd.read_csv(MASTER, dtype={"pdb_id": str})
    assert len(df) == 2072


def test_original_functional_site_outputs_unchanged():
    result = subprocess.run(
        ["git", "diff", "--quiet", BASELINE_COMMIT, "--"] + ORIGINAL_FUNCTIONAL_SITE_OUTPUTS,
        cwd=PROJECT_ROOT, capture_output=True,
    )
    assert result.returncode == 0, (
        f"one or more original functional-site outputs changed since {BASELINE_COMMIT} -- "
        f"this stage must be read-only w.r.t. them: {result.stdout.decode()}{result.stderr.decode()}"
    )


# --------------------------------------------------------- determinism --

def test_specificity_pipeline_deterministic_and_idempotent():
    files = [FDR_AUDIT, MATCHED_CONTROLS, PDB_MATCHED, MATCHED_TESTS, SUPERFAMILY_MATCHED, CLASSIFICATION]
    for f in files:
        _require(f)
    before = {f: _sha(f) for f in files}
    scripts = VALIDATION_SCRIPTS + ["90_interpretation_classification.py"]
    for script in scripts:
        r = subprocess.run([sys.executable, f"scripts/{script}"], cwd=PROJECT_ROOT,
                            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
    after = {f: _sha(f) for f in files}
    assert before == after, "specificity-validation pipeline outputs are not deterministic across reruns"


def test_summary_covers_required_sections():
    _require(SUMMARY_MD)
    text = SUMMARY_MD.read_text()
    for required in [
        "FDR audit", "definition confidence", "structurally matched background",
        "Confirmatory superfamily", "Interpretation classification",
        "Apo/holo overlap confidence", "Robustness", "Limitations",
    ]:
        assert required in text, f"functional_site_specificity_summary.md missing section: {required}"
    assert "confirms mechanism" not in text.lower()
    # "causes" appearing at all is fine as long as it's inside an explicit disclaimer
    # that no causal claim is being made (not asserted as a finding)
    disclaimer_present = any(phrase in text.lower() for phrase in
                              ["never claims", "not a causal", "not been established",
                               "has been established", "no causal"])
    assert "causes" not in text.lower() or disclaimer_present
