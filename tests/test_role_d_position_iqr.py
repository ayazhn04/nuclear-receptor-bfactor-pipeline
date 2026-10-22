"""Focused validation of the v2-only Role D position-IQR outputs."""
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_position_iqr_is_recomputable_and_canonical_v2_only():
    subprocess.run([sys.executable, "scripts/109_role_d_position_iqr.py"], cwd=ROOT, check=True)
    out = pd.read_csv(ROOT / "reports/role_d/contextual_iqr/role_d_functional_position_iqr.csv")
    observed = out.iqr_z.notna()
    assert len(out) > 0
    assert np.allclose(out.loc[observed, "iqr_z"], out.loc[observed, "q75_z"] - out.loc[observed, "q25_z"])
    assert out.functional_site.isin(["AF2_H12", "LIGAND_BINDING_POCKET"]).all()
    assert out.definition_level.eq("RESIDUE").all()
    assert not out.context_stratum_display.str.contains("loop|beta_sheet", case=False).any()


def test_primary_excess_iqr_has_declared_support_threshold():
    out = pd.read_csv(ROOT / "reports/role_d/contextual_iqr/role_d_functional_position_iqr.csv")
    eligible = out[out.iqr_eligibility.eq("ELIGIBLE_PRIMARY")]
    assert (eligible.n_pdb_entries_observed_at_position >= 5).all()
    assert np.allclose(eligible.excess_iqr_z, eligible.iqr_z - eligible.receptor_median_position_iqr_z)


def test_summary_primary_metric_excludes_sparse_positions():
    out = pd.read_csv(ROOT / "reports/role_d/contextual_iqr/role_d_functional_position_iqr.csv")
    summary = pd.read_csv(ROOT / "reports/role_d/contextual_iqr/role_d_functional_site_iqr_summary.csv")
    keys = ["uniprot_id", "nr_code", "common_name", "functional_site"]
    expected = (out[out.eligible_primary]
                .groupby(keys).iqr_z.median()
                .rename("expected").reset_index())
    checked = summary.merge(expected, on=keys, how="left")
    assert np.allclose(checked.median_position_iqr_z, checked.expected, equal_nan=True)


def test_figure_is_reproducible_headlessly():
    subprocess.run([sys.executable, "scripts/110_role_d_position_iqr_figures.py"], cwd=ROOT, check=True)
    figure = ROOT / "reports/role_d/contextual_iqr/role_d_functional_position_iqr_scatter.png"
    assert figure.is_file() and figure.stat().st_size > 0


def test_lbp_matched_analysis_has_receptor_level_inference():
    subprocess.run([sys.executable, "scripts/111_role_d_lbp_matched_iqr.py"], cwd=ROOT, check=True)
    results = pd.read_csv(ROOT / "reports/role_d/contextual_iqr/role_d_lbp_matched_iqr_results.csv")
    summary = results[results.record_type.eq("SUMMARY")].iloc[0]
    receptor = results[results.record_type.eq("RECEPTOR")]
    assert summary.analysis == "LBP_VS_NONFUNCTIONAL_EXACT_CONTEXT_SUPPORT_MATCHED"
    assert summary.n_receptors == len(receptor)
    assert summary.n_permutations == 10_000
    assert 0 <= summary.restricted_permutation_p_two_sided <= 1
