"""
Analysis-foundation tests: within-instance normalization, PDB-level pseudo-
replication control, receptor x position profiles, and protection of the
master inventory, strict-QC historical files and Role D outputs.

Offline; uses only tracked files under data/ and reports/ (no raw mmCIF/API
caches needed).
"""

from __future__ import annotations

import hashlib
import io
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
REPORTS = PROJECT_ROOT / "reports" / "tables"

MASTER = PROCESSED / "structures.csv"
RESIDUES = PROCESSED / "lbd_residue_bfactors_all.parquet"
NORMALIZED = ANALYSIS / "normalized_lbd_bfactors_all.parquet"
INSTANCE_SUMMARY = ANALYSIS / "instance_bfactor_summary.csv"
PDB_PROFILE = ANALYSIS / "pdb_position_bfactor_profiles.parquet"
RECEPTOR_PROFILE = ANALYSIS / "receptor_position_flexibility_profiles.parquet"
RECEPTOR_SUMMARY = ANALYSIS / "receptor_flexibility_summary.csv"

# Commit at which the historical / Role D files were last authored by others.
BASELINE_COMMIT = "e11ddc7"
PROTECTED_PATHS = [
    "data/processed/structures_primary_qc_subset.csv",
    "data/processed/final_lbd_residue_map_primary_qc_subset.parquet",
    "data/processed/lbd_residue_bfactors_all.parquet",
    "data/processed/excluded_structures.csv",
    "scripts/31_analyze_functional_site_bfactors.py",
    "reports/ROLE_D_PROGRESS_REPORT.md",
    "reports/tables/functional_site_definitions.csv",
    "reports/tables/functional_site_entry_scores.csv",
    "reports/tables/functional_site_normalization_qc.csv",
    "reports/tables/functional_site_receptor_scores.csv",
    "tests/test_functional_site_analysis.py",
]
SUBSET_CSV_SHA256 = "5222cfd462cb558c421f2a6780952612e8b7d4bedd6a7930de287afeb90ed58b"
SUBSET_PARQUET_SHA256 = "ce118a4d049d5b4123daf460bffaca965493cb1154071a25bfe15afccadfacbc"
RESIDUES_SHA256 = "28cc0f7feff59d27cf6b724becad23391068c6ed74181db3a6b9923750f81a27"
VDR = "P11473"


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@pytest.fixture(scope="module")
def master():
    return pd.read_csv(MASTER, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def raw():
    return pd.read_parquet(RESIDUES)


@pytest.fixture(scope="module")
def norm():
    return pd.read_parquet(NORMALIZED)


@pytest.fixture(scope="module")
def inst():
    return pd.read_csv(INSTANCE_SUMMARY, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def pdb_prof():
    return pd.read_parquet(PDB_PROFILE)


@pytest.fixture(scope="module")
def rec_prof():
    return pd.read_parquet(RECEPTOR_PROFILE)


@pytest.fixture(scope="module")
def summ():
    return pd.read_csv(RECEPTOR_SUMMARY)


# ---------------------------------------------------------------- normalization

def test_normalized_row_count_equals_raw(norm, raw):
    assert len(norm) == len(raw) == 771542


def test_existing_residue_columns_preserved_and_keys_unchanged(norm, raw):
    for c in raw.columns:
        assert norm[c].equals(raw[c]), c
    for c in ["instance_raw_b_mean", "instance_raw_b_median", "instance_raw_b_sd", "instance_raw_b_q25",
              "instance_raw_b_q75", "instance_raw_b_iqr", "b_factor_z", "normalization_status"]:
        assert c in norm.columns


def test_missing_residue_rows_explicit_and_stay_missing(norm, raw):
    missing = raw["selected_ca_b_iso_or_equiv"].isna()
    assert missing.sum() == 50697
    assert norm.loc[missing, "b_factor_z"].isna().all()
    assert norm.loc[missing, "selected_ca_b_iso_or_equiv"].isna().all()
    assert (norm.loc[missing, "row_normalization_status"] == "B_MISSING").all()
    assert not raw.duplicated(["instance_id", "lbd_relative_position"]).any()


def test_finite_z_only_from_finite_usable_raw_b(norm):
    fin = norm["b_factor_z"].notna()
    assert np.isfinite(norm.loc[fin, "b_factor_z"]).all()
    assert norm.loc[fin, "selected_ca_b_iso_or_equiv"].notna().all()
    assert norm.loc[fin, "ca_bfactor_usable"].all()
    assert not np.isinf(norm["b_factor_z"].astype(float)).any()


def test_normalized_instances_mean_zero_population_sd_one(norm):
    g = norm.dropna(subset=["b_factor_z"]).groupby("instance_id")["b_factor_z"]
    assert g.mean().abs().max() < 1e-9
    assert (g.std(ddof=0) - 1).abs().max() < 1e-9


def test_z_matches_formula_population_sd(norm):
    for iid in norm["instance_id"].drop_duplicates().sample(25, random_state=7):
        d = norm[norm["instance_id"] == iid]
        b = d["selected_ca_b_iso_or_equiv"].dropna()
        if len(b) < 2 or b.std(ddof=0) == 0:
            continue
        expect = (d["selected_ca_b_iso_or_equiv"] - b.mean()) / b.std(ddof=0)
        assert np.allclose(d["b_factor_z"].dropna(), expect.dropna(), atol=1e-9)


def test_positive_occupancy_semantics_unchanged(norm, raw):
    for c in ["ca_record_present", "ca_positive_occupancy", "ca_bfactor_usable", "zero_occupancy_only_ca"]:
        assert (norm[c].values == raw[c].values).all()
    assert not (norm["ca_bfactor_usable"] & ~norm["ca_positive_occupancy"]).any()


def test_non_normalizable_instance_retained_with_explicit_status(norm, inst):
    bad = inst[inst["normalization_status"] != "NORMALIZED"]
    assert len(bad) >= 1
    assert bad["normalization_status"].str.startswith("NOT_NORMALIZABLE").all()
    for iid in bad["instance_id"]:
        rows = norm[norm["instance_id"] == iid]
        assert len(rows) > 0 and rows["b_factor_z"].isna().all()


# ---------------------------------------------------------------- instance summary / master

def test_instance_summary_one_row_per_instance_all_retained(inst, raw):
    assert len(inst) == raw["instance_id"].nunique() == 3159
    assert inst["instance_id"].is_unique
    assert (inst["n_ca_usable"] + inst["n_ca_missing"] == inst["n_lbd_rows"]).all()


def test_no_qc_threshold_removes_candidates_from_instances(inst):
    assert (inst["positive_occupancy_lbd_coverage"] < 0.5).sum() > 0     # low coverage retained
    assert (~inst["previous_primary_qc_member"]).sum() > 0                  # non-strict retained
    assert inst["resolution"].isna().sum() >= 0
    assert set(inst["ligand_state"]) >= {"APO", "HOLO"}


def test_master_still_2072_rows_with_all_candidates(master):
    assert len(master) == 2072
    assert int(master["has_stage3b_coordinate_data"].sum()) == 1840
    assert not master.duplicated(["uniprot_id", "pdb_id", "polymer_entity_id"]).any()


def test_master_existing_columns_unchanged_and_descriptive_columns_added(master):
    old = pd.read_csv(io.StringIO(subprocess.run(
        ["git", "show", f"{BASELINE_COMMIT}:data/processed/structures.csv"], cwd=PROJECT_ROOT,
        capture_output=True, text=True, check=True).stdout), dtype=str, keep_default_na=False)
    new = pd.read_csv(MASTER, dtype=str, keep_default_na=False)
    assert old.equals(new[old.columns])
    added = [c for c in new.columns if c not in old.columns]
    assert added == [c for c in added if c.startswith("representative_")] and len(added) == 8
    assert not any("z_mean" in c or "structure_mean_z" in c for c in master.columns)


def test_representative_mapping_unique_and_complete(master, inst):
    coord = master[master["has_stage3b_coordinate_data"]]
    assert coord["representative_instance_id"].is_unique and len(coord) == 1840
    reps = inst[inst["is_representative_instance"]]
    assert set(reps["instance_id"]) == set(coord["representative_instance_id"])
    assert master.loc[~master["has_stage3b_coordinate_data"], "representative_n_ca_usable"].isna().all()
    assert (master.loc[~master["has_stage3b_coordinate_data"], "coordinate_data_status"] == "NO_STAGE3B_COORDINATE_DATA").all()


# ---------------------------------------------------------------- PDB-level / receptor profiles

def test_pdb_profile_one_observation_per_receptor_pdb_position(pdb_prof):
    assert not pdb_prof.duplicated(["uniprot_id", "pdb_id", "canonical_uniprot_position"]).any()


def test_duplicate_candidates_in_one_pdb_reduced_by_median_not_double_counted(pdb_prof, norm, inst):
    multi = pdb_prof[pdb_prof["n_contributing_representative_candidates"] > 1]
    assert len(multi) > 0
    row = multi.iloc[0]
    reps = set(inst.loc[inst["is_representative_instance"], "instance_id"])
    d = norm[(norm["uniprot_id"] == row["uniprot_id"]) & (norm["pdb_id"] == row["pdb_id"])
             & (norm["canonical_uniprot_position"] == row["canonical_uniprot_position"])
             & norm["instance_id"].isin(reps)].dropna(subset=["b_factor_z"])
    assert d["instance_id"].nunique() == row["n_contributing_representative_candidates"] > 1
    assert row["b_factor_z"] == pytest.approx(d["b_factor_z"].median())


def test_receptor_counts_count_each_pdb_once(pdb_prof, rec_prof):
    c = pdb_prof.groupby(["uniprot_id", "canonical_uniprot_position"])["pdb_id"].nunique()
    m = rec_prof.set_index(["uniprot_id", "canonical_uniprot_position"])["n_pdb_entries_observed_at_position"]
    common = c.index.intersection(m.index)
    assert (c.loc[common].values == m.loc[common].values).all()
    assert (rec_prof["n_pdb_entries_observed_at_position"] <= rec_prof["n_pdb_entries_total_for_receptor"]).all()
    assert not rec_prof.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()


def test_apo_holo_summaries_exclude_ambiguous_and_not_assessable(pdb_prof, rec_prof):
    assert set(pdb_prof["ligand_state"]) <= {"APO", "HOLO", "AMBIGUOUS"}
    conflicts = pdb_prof[pdb_prof["ligand_state_conflict"]]
    assert (conflicts["ligand_state"] == "AMBIGUOUS").all()
    apo = pdb_prof[pdb_prof["ligand_state"] == "APO"].groupby(["uniprot_id", "canonical_uniprot_position"]).size()
    holo = pdb_prof[pdb_prof["ligand_state"] == "HOLO"].groupby(["uniprot_id", "canonical_uniprot_position"]).size()
    r = rec_prof.set_index(["uniprot_id", "canonical_uniprot_position"])
    assert (r.loc[apo.index, "apo_n"].values == apo.values).all()
    assert (r.loc[holo.index, "holo_n"].values == holo.values).all()
    assert ((r["apo_n"] + r["holo_n"]) <= r["n_pdb_entries_observed_at_position"]).all()
    only_one = (r["apo_n"] == 0) | (r["holo_n"] == 0)
    assert r.loc[only_one, "apo_holo_mean_delta"].isna().all()
    both = (r["apo_n"] > 0) & (r["holo_n"] > 0)
    assert np.allclose(r.loc[both, "apo_holo_mean_delta"], r.loc[both, "holo_mean_z"] - r.loc[both, "apo_mean_z"])


def test_all_48_receptors_in_summary_zero_structure_rows_have_no_invented_values(summ):
    assert len(summ) == 48 and summ["uniprot_id"].is_unique
    empty = summ[summ["n_normalizable_pdbs"] == 0]
    assert len(empty) == 6
    for c in ["median_position_median_z", "iqr_of_position_median_z", "fraction_positions_median_z_gt_1_0",
              "median_lbd_observation_fraction"]:
        assert empty[c].isna().all()
    assert set(empty["profile_status"]) == {"NO_STRUCTURES_IN_MASTER", "NO_COORDINATE_DATA"}


# ---------------------------------------------------------------- VDR

def test_vdr_present_everywhere_and_not_excluded_for_coverage(master, norm, inst, pdb_prof, rec_prof, summ):
    assert (norm["uniprot_id"] == VDR).sum() == 14256
    v = inst[inst["uniprot_id"] == VDR]
    assert len(v) == 48 and (v["normalization_status"] == "NORMALIZED").all()
    assert (v["positive_occupancy_lbd_coverage"] < 0.90).all() and not v["previous_primary_qc_member"].any()
    assert (pdb_prof["uniprot_id"] == VDR).any() and (rec_prof["uniprot_id"] == VDR).any()
    s = summ[summ["uniprot_id"] == VDR].iloc[0]
    assert s["n_normalizable_pdbs"] == 48 and s["profile_status"] == "PROFILE_AVAILABLE"
    assert (master["uniprot_id"] == VDR).sum() == 52


# ---------------------------------------------------------------- protection

def test_historical_strict_qc_files_byte_identical():
    assert _sha(PROCESSED / "structures_primary_qc_subset.csv") == SUBSET_CSV_SHA256
    assert _sha(PROCESSED / "final_lbd_residue_map_primary_qc_subset.parquet") == SUBSET_PARQUET_SHA256
    assert _sha(RESIDUES) == RESIDUES_SHA256


def test_protected_and_role_d_files_unchanged_since_baseline():
    r = subprocess.run(["git", "diff", "--quiet", BASELINE_COMMIT, "--", *PROTECTED_PATHS], cwd=PROJECT_ROOT)
    assert r.returncode == 0, "a protected historical / Role D file differs from the baseline commit"


def test_resolution_diagnostic_is_descriptive_and_leaves_z_untouched(norm):
    d = pd.read_csv(REPORTS / "raw_bfactor_resolution_diagnostic.csv")
    assert {"pearson_r", "spearman_rho", "n_pairs"} <= set(d.columns)
    assert set(d["y"]) == {"representative_raw_b_mean", "representative_raw_b_median"}
    assert (d.loc[d["scope"] == "all_representative_coordinate_instances", "n_pairs"] > 1000).all()
    # per-instance z remains exactly mean 0 irrespective of the diagnostic
    g = norm.dropna(subset=["b_factor_z"]).groupby("instance_id")["b_factor_z"].mean()
    assert g.abs().max() < 1e-9


def test_reports_exist_and_alignment_report_states_decision_pending():
    for name in ["analysis_foundation_missingness_report.md", "common_lbd_alignment_readiness.md"]:
        assert (REPORTS / name).exists()
    text = (REPORTS / "common_lbd_alignment_readiness.md").read_text().lower()
    assert "no common multi-receptor alignment coordinate exists yet" in text
    assert "none was installed" in text
