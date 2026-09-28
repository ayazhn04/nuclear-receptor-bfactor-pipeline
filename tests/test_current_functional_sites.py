"""
Current full-dataset functional-site flexibility tests (continuing from
commit fb27a86 on analysis/normalization-v1; scripts 72-82).

Covers: historical Role D files remain byte-identical; master
structures.csv unchanged at 2072 rows; unresolved sites are never zero-
filled; ligand-pocket contacts exclude irrelevant HETATM species; receptor
is the superfamily inferential unit; PDB entry is the within-receptor
observational unit; paired tests use matched PDB site/background
observations; BH correction is applied; site definitions carry provenance;
AF2/H12 uses the current v3 structural mapping; apo/holo results are only
linked descriptively (never recomputed); upstream files unchanged; and
full-stage determinism.
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
DEFINITIONS = ANALYSIS / "current_functional_site_definitions.parquet"
AVAILABILITY = TABLES / "current_functional_site_availability.csv"
PDB_SITE = ANALYSIS / "pdb_functional_site_flexibility.parquet"
RECEPTOR_SITE = ANALYSIS / "receptor_functional_site_flexibility.parquet"
WITHIN_TESTS = TABLES / "current_functional_site_within_receptor_tests.csv"
SUPERFAMILY_TESTS = TABLES / "current_functional_site_superfamily_tests.csv"
APO_HOLO_OVERLAP = TABLES / "apo_holo_functional_site_overlap.csv"
POCKET_CONSENSUS = ANALYSIS / "ligand_pocket_consensus_residues.parquet"
SUMMARY_MD = TABLES / "current_functional_site_summary.md"

BASELINE_COMMIT = "fb27a86"  # start of this stage

ROLE_D_PATHS = [
    "scripts/31_analyze_functional_site_bfactors.py",
    "reports/ROLE_D_PROGRESS_REPORT.md",
    "reports/tables/functional_site_definitions.csv",
    "reports/tables/functional_site_entry_scores.csv",
    "reports/tables/functional_site_normalization_qc.csv",
    "reports/tables/functional_site_receptor_scores.csv",
    "tests/test_functional_site_analysis.py",
]

UPSTREAM_PATHS = [
    "data/analysis/nr_lbd_structural_segment_map.csv",
    "data/analysis/receptor_segment_feature_matrix_median.parquet",
    "data/analysis/pca_receptor_scores.csv",
    "data/analysis/clustering_summary.json",
    "reports/tables/subfamily_segment_kruskal.csv",
    "data/analysis/superfamily_segment_flexibility_atlas.csv",
    "data/analysis/apo_holo_segment_effects.parquet",
    "data/analysis/apo_holo_residue_effects.parquet",
    "reports/tables/apo_holo_segment_tests.csv",
]

FUNCTIONAL_SITE_SCRIPTS = [
    "73_current_functional_site_definitions.py",
    "74_pdb_site_flexibility.py",
    "75_receptor_site_summary.py",
    "76_within_receptor_site_tests.py",
    "77_superfamily_site_tests.py",
    "78_apo_holo_site_overlap.py",
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


@pytest.fixture(scope="module")
def definitions():
    _require(DEFINITIONS)
    return pd.read_parquet(DEFINITIONS)


@pytest.fixture(scope="module")
def pdb_site():
    _require(PDB_SITE)
    return pd.read_parquet(PDB_SITE)


@pytest.fixture(scope="module")
def receptor_site():
    _require(RECEPTOR_SITE)
    return pd.read_parquet(RECEPTOR_SITE)


@pytest.fixture(scope="module")
def within_tests():
    _require(WITHIN_TESTS)
    return pd.read_csv(WITHIN_TESTS)


# --------------------------------------------- Role D / protected files --

def test_role_d_files_byte_identical_to_baseline():
    result = subprocess.run(
        ["git", "diff", "--quiet", BASELINE_COMMIT, "--"] + ROLE_D_PATHS,
        cwd=PROJECT_ROOT, capture_output=True,
    )
    assert result.returncode == 0, (
        f"historical Role D file(s) changed since {BASELINE_COMMIT}: "
        f"{result.stdout.decode()}{result.stderr.decode()}"
    )


def test_new_definitions_table_is_separate_from_role_d():
    _require(DEFINITIONS)
    role_d = pd.read_csv(TABLES / "functional_site_definitions.csv")
    new = pd.read_parquet(DEFINITIONS)
    assert set(new.columns) != set(role_d.columns) or len(new) != len(role_d)
    assert len(new) > len(role_d), "the new full-dataset layer should cover far more than the 44 Role-D rows"


def test_master_structures_csv_untouched():
    df = pd.read_csv(MASTER, dtype={"pdb_id": str})
    assert len(df) == 2072


def test_upstream_files_unchanged():
    result = subprocess.run(
        ["git", "diff", "--quiet", BASELINE_COMMIT, "--"] + UPSTREAM_PATHS,
        cwd=PROJECT_ROOT, capture_output=True,
    )
    assert result.returncode == 0, (
        f"upstream clustering/subfamily/apo-holo file(s) changed since {BASELINE_COMMIT}: "
        f"{result.stdout.decode()}{result.stderr.decode()}"
    )
    _require(FROZEN_HASHES)
    with open(FROZEN_HASHES) as f:
        frozen = json.load(f)
    for rel_path, expected_hash in frozen.items():
        assert _sha(PROJECT_ROOT / rel_path) == expected_hash


# ----------------------------------------------------------- Part 1, 6 ---

def test_unresolved_sites_never_zero_filled(receptor_site):
    _require(AVAILABILITY)
    availability = pd.read_csv(AVAILABILITY)
    unresolved = availability[availability["status"] == "UNRESOLVED"]
    assert len(unresolved) > 0, "expected at least one genuinely unresolved receptor x site"
    for _, r in unresolved.iterrows():
        sub = receptor_site[(receptor_site["uniprot_id"] == r["uniprot_id"]) &
                             (receptor_site["functional_site"] == r["functional_site"])]
        assert len(sub) == 0, (
            f"{r['uniprot_id']}/{r['functional_site']} is UNRESOLVED but has a receptor-level "
            "flexibility row -- unresolved sites must never receive a value"
        )


def test_availability_statuses_are_valid_enum():
    _require(AVAILABILITY)
    availability = pd.read_csv(AVAILABILITY)
    valid = {"AVAILABLE", "PARTIAL", "AMBIGUOUS", "NOT_APPLICABLE", "UNRESOLVED"}
    assert set(availability["status"]).issubset(valid)


def test_site_definitions_have_provenance(definitions):
    required_cols = {"uniprot_id", "nr_code", "functional_site", "canonical_uniprot_position",
                      "structural_segment", "definition_level", "assignment_status",
                      "assignment_confidence", "evidence_type", "source_or_method", "notes"}
    assert required_cols.issubset(set(definitions.columns))
    assert (definitions["evidence_type"] != "").all()
    assert (definitions["source_or_method"].astype(str).str.len() > 0).all()
    assert set(definitions["definition_level"]) <= {"RESIDUE", "SEGMENT"}


def test_af2_h12_uses_current_structural_mapping(definitions):
    af2 = definitions[definitions["functional_site"] == "AF2_H12"]
    assert len(af2) > 0
    assert (af2["structural_segment"] == "H12").all()
    assert (af2["evidence_type"] == "current_structural_segment_map").all()
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    h12_map = seg_map[seg_map["structural_segment"] == "H12"]
    merged = af2.merge(h12_map, on=["uniprot_id", "canonical_uniprot_position"], suffixes=("_def", "_map"))
    assert len(merged) == len(af2), "every AF2_H12 position must trace back to the current v3 H12 mapping"


# ------------------------------------------------------ Part 3 pocket ----

def test_pocket_excludes_irrelevant_hetatm_species():
    inventory = pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "final_nonpolymer_inventory.csv",
                             dtype={"pdb_id": str}, keep_default_na=False)
    functional_ids = set(inventory.loc[inventory["nonpolymer_relevance_class"] == "FUNCTIONAL_LBD_LIGAND", "component_id"])
    # water and common crystallization additives/ions must never be classified functional
    assert "HOH" not in functional_ids
    assert "DOD" not in functional_ids
    excluded_classes = set(inventory["nonpolymer_relevance_class"].unique())
    assert "WATER" in excluded_classes
    water_rows = inventory[inventory["component_id"].isin(["HOH", "DOD"])]
    assert (water_rows["nonpolymer_relevance_class"] == "WATER").all()


def test_pocket_consensus_reports_required_columns():
    _require(POCKET_CONSENSUS)
    consensus = pd.read_parquet(POCKET_CONSENSUS)
    required = {"uniprot_id", "canonical_uniprot_position", "n_holo_pdbs_observed",
                "n_holo_pdbs_contacting", "contact_fraction"}
    assert required.issubset(set(consensus.columns))
    assert (consensus["n_holo_pdbs_contacting"] <= consensus["n_holo_pdbs_observed"]).all()
    assert ((consensus["contact_fraction"] >= 0) & (consensus["contact_fraction"] <= 1)).all()


# ----------------------------------------------------- Part 7 unit -------

def test_pdb_entry_is_within_receptor_observation_unit(pdb_site):
    assert "chain" not in pdb_site.columns and "label_asym_id" not in pdb_site.columns
    dup = pdb_site.duplicated(subset=["uniprot_id", "pdb_id", "functional_site"]).any()
    assert not dup


def test_no_site_score_without_observed_site_residues(pdb_site):
    assert (pdb_site["n_site_residues_observed"] > 0).all()


# ----------------------------------------------- Part 9 receptor unit ----

def test_receptor_is_the_superfamily_inferential_unit():
    _require(SUPERFAMILY_TESTS)
    superfamily = pd.read_csv(SUPERFAMILY_TESTS)
    receptor_site_df = pd.read_parquet(RECEPTOR_SITE)
    for _, r in superfamily.iterrows():
        n = (receptor_site_df["functional_site"] == r["functional_site"]).sum()
        assert n == r["n_receptors"]
    # superfamily test must not read the PDB-level or position-level table directly
    text = (PROJECT_ROOT / "scripts" / "77_superfamily_site_tests.py").read_text()
    assert "pdb_functional_site_flexibility" not in text
    assert "pdb_position_bfactor_profiles" not in text


# ------------------------------------------------- Part 10 paired test ---

def test_paired_test_uses_matched_pdb_observations(within_tests, pdb_site):
    row = within_tests.iloc[0]
    uid, site = row["uniprot_id"], row["functional_site"]
    g = pdb_site[(pdb_site["uniprot_id"] == uid) & (pdb_site["functional_site"] == site)]
    assert len(g) == row["n_pairs"]
    from scipy import stats
    diffs = g["median_site_z"].values - g["comparison_non_site_median_z"].values
    if not np.all(diffs == 0):
        stat, p = stats.wilcoxon(g["median_site_z"], g["comparison_non_site_median_z"],
                                  alternative="two-sided", zero_method="wilcox")
        assert np.isclose(p, row["raw_p"])


def test_bh_correction_applied_within_and_globally(within_tests):
    tested = within_tests.dropna(subset=["raw_p", "q_within_receptor", "q_global"])
    assert len(tested) > 0
    assert (tested["q_within_receptor"] >= tested["raw_p"] - 1e-9).all()
    assert (tested["q_global"] >= tested["raw_p"] - 1e-9).all()
    assert not tested["q_within_receptor"].equals(tested["q_global"])
    assert not (tested["significant_within_receptor"] & (tested["q_within_receptor"] > 0.05)).any()


# -------------------------------------------- Part 13 apo/holo link ------

def test_apo_holo_results_only_linked_not_recomputed():
    text = (PROJECT_ROOT / "scripts" / "78_apo_holo_site_overlap.py").read_text()
    forbidden = ["mannwhitneyu", "wilcoxon", "kruskal"]
    for f in forbidden:
        assert f not in text.lower(), (
            f"78_apo_holo_site_overlap.py must only read stored apo/holo results, never "
            f"recompute a statistical test ({f} found)"
        )
    assert "apo_holo_segment_effects.parquet" in text


def test_apo_holo_overlap_table_matches_stored_effects():
    _require(APO_HOLO_OVERLAP)
    overlap = pd.read_csv(APO_HOLO_OVERLAP, keep_default_na=False)
    stored = pd.read_parquet(ANALYSIS / "apo_holo_segment_effects.parquet")
    merged = overlap.drop_duplicates(["uniprot_id", "apo_holo_segment"]).merge(
        stored, left_on=["uniprot_id", "apo_holo_segment"], right_on=["uniprot_id", "structural_segment"]
    )
    assert np.allclose(merged["apo_holo_delta_median_z"], merged["delta_median_z"])
    assert np.allclose(merged["apo_holo_raw_p"], merged["raw_p"])


# ---------------------------------------------------------- determinism -

def test_functional_site_pipeline_deterministic_and_idempotent():
    files = [DEFINITIONS, AVAILABILITY, PDB_SITE, RECEPTOR_SITE, WITHIN_TESTS, SUPERFAMILY_TESTS]
    for f in files:
        _require(f)
    before = {f: _sha(f) for f in files}
    for script in FUNCTIONAL_SITE_SCRIPTS:
        r = subprocess.run([sys.executable, f"scripts/{script}"], cwd=PROJECT_ROOT,
                            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
    after = {f: _sha(f) for f in files}
    assert before == after, "functional-site pipeline outputs are not deterministic across reruns"


def test_summary_covers_required_sections():
    _require(SUMMARY_MD)
    text = SUMMARY_MD.read_text()
    for required in [
        "Site availability", "Site definitions and provenance",
        "Within-receptor site-vs-background", "Superfamily site-class results",
        "Effect sizes", "flexibility atlas", "apo/holo effects", "Robustness",
        "Limitations",
    ]:
        assert required in text, f"current_functional_site_summary.md missing section: {required}"
