"""Synthetic tests; no network, no writes to upstream scientific outputs."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sasa = load("113_role_d_sasa.py")
analysis = load("114_role_d_sasa_iqr_analysis.py")


def test_altloc_coherent_and_occupancy_not_b_selected():
    atoms = pd.DataFrame({"seq": [1]*7, "atom": ["N", "CA", "CA", "CB", "CB", "C", "O"],
                          "alt": [".", "A", "B", "A", "B", ".", "."],
                          "occupancy": [1, .7, .3, .4, .6, 1, 0], "element": ["N", "C", "C", "C", "C", "C", "O"]})
    selected = sasa.select_atoms(atoms, {1: "B"})
    assert set(selected.atom) == {"N", "CA", "CB", "C"}
    assert set(selected.loc[selected.atom.isin(["CA", "CB"]), "alt"]) == {"B"}
    assert not selected.duplicated(["seq", "atom"]).any()


def test_missing_sidechain_not_assigned_falsely_buried_rsa():
    atoms = pd.DataFrame({"seq": [1]*4, "atom": ["N", "CA", "C", "O"], "res": ["ALA"]*4,
                          "alt": ["."]*4, "element": ["N", "C", "C", "O"], "occupancy": [1]*4,
                          "x": [0, 1.4, 2.8, 3.8], "y": [0]*4, "z": [0]*4})
    mapping = pd.DataFrame([dict(uniprot_id="TEST", pdb_id="TEST", instance_id="TEST.A", label_asym_id="A", model_num="1",
                                canonical_uniprot_position=10, entity_label_seq_id=1,
                                selected_ca_altloc="", ca_bfactor_usable=True, canonical_residue_name="ALA")])
    result = sasa.calculate_chain(atoms, mapping, 100)
    assert result.sasa_status.iloc[0] == "INCOMPLETE_HEAVY_ATOMS"
    atoms["res"] = "GLY"
    result = sasa.calculate_chain(atoms, mapping, 100)
    assert result.sasa_status.iloc[0] == "OK"
    assert np.isclose(result.rsa.iloc[0], result.sasa_a2.iloc[0] / sasa.MAX_ASA["GLY"])
    assert not result.residue_matches_canonical.iloc[0]  # observed identity used, not canonical ALA


def test_isolated_atom_area_matches_analytic_sphere():
    import biotite.structure as bs
    array = bs.AtomArray(1)
    array.coord[:] = 0
    array.element[:] = "C"
    area = bs.sasa(array, probe_radius=1.4, point_number=1000, vdw_radii=np.array([1.7]), ignore_ions=False)
    assert np.isclose(area[0], 4 * np.pi * (1.7 + 1.4)**2, rtol=1e-5)


def test_context_matching_does_not_use_cross_receptor_or_h12_controls():
    rows = []
    for uid in ["A", "B"]:
        for i in range(5):
            rows.append(dict(uniprot_id=uid, common_name=uid, eligible=True, is_lbp=i == 0,
                             is_any_functional=i in [0, 4], context="H3" if i != 4 else "H12",
                             support_band="20_PLUS", rsa_band="BURIED_LE_005", paired_iqr_z=float(i),
                             canonical_residue_name="ALA", review_required_flag=False))
    strata, receptors = analysis.matched_effects(pd.DataFrame(rows), "CONTEXT_SUPPORT_RSA")
    assert len(strata) == len(receptors) == 2
    assert (strata.n_control == 3).all()
    assert (receptors.effect == -2).all()
    assert not strata.match_stratum.str.contains("H12").any()
    rows[1]["rsa_band"] = "EXPOSED_GT_025"
    strata, receptors = analysis.matched_effects(pd.DataFrame(rows), "CONTEXT_SUPPORT_RSA")
    assert list(receptors.uniprot_id) == ["B"]  # only two remaining controls in A


def test_fdr_and_receptor_bootstrap_are_reproducible():
    assert np.allclose(analysis.bh([.01, .04, .03]), [.03, .04, .04])
    assert analysis.summarize_effects([1, 2, 3, 4, 5]) == analysis.summarize_effects([1, 2, 3, 4, 5])
    summary = analysis.summarize_effects([1, 2, 3, 4, 5])
    assert summary['n_receptors'] == 5
    assert summary['ci95_low'] <= summary['median_effect'] <= summary['ci95_high']


def test_original_v2_iqr_recomputes_from_entry_profiles():
    entries = pd.read_parquet(analysis.ANALYSIS / "pdb_position_bfactor_profiles.parquet")
    keys = ["uniprot_id", "canonical_uniprot_position"]
    q = entries.groupby(keys).b_factor_z.quantile([.25, .75]).unstack()
    q["recomputed_iqr"] = q[.75] - q[.25]
    profile = pd.read_csv(analysis.ANALYSIS / "receptor_position_professor_profile_v2.csv")
    check = profile.merge(q[["recomputed_iqr"]], on=keys, how="left", validate="one_to_one")
    assert np.allclose(check.iqr_z, check.recomputed_iqr, equal_nan=True)


def test_generated_sasa_outputs_keep_missingness_and_unique_pdb_positions():
    import pytest
    if not (analysis.OUT / "pdb_position_sasa.parquet").exists():
        pytest.skip("Full SASA run has not finished")
    instances = pd.read_parquet(analysis.OUT / "representative_instance_sasa.parquet")
    usable = instances.sasa_status.eq("OK")
    assert instances.loc[~usable, "rsa"].isna().all()
    assert (instances.loc[usable, "rsa"] >= 0).all()
    entries = pd.read_parquet(analysis.OUT / "pdb_position_sasa.parquet")
    assert not entries.duplicated(["uniprot_id", "pdb_id", "canonical_uniprot_position"]).any()
    profile, paired = analysis.build_positions()
    assert (profile.loc[profile.eligible, "n_paired_pdb"] >= 5).all()
    assert np.allclose(profile.paired_iqr_z, profile.paired_q75_z - profile.paired_q25_z, equal_nan=True)
    assert not paired.duplicated(["uniprot_id", "pdb_id", "canonical_uniprot_position"]).any()
