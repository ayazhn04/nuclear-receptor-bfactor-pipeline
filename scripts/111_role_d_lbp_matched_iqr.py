"""Role D v2: matched-context LBP IQR contrast.

This tests a position-level association in between-PDB crystallographic-disorder
heterogeneity.  The inferential unit is a receptor, never an individual residue
or PDB entry. AF-2/H12 is deliberately excluded: H12 itself defines that site,
so exact same-H12 non-site controls do not exist.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "data/rebuild_v2/analysis/receptor_position_professor_profile_v2.csv"
FUNCTIONAL = ROOT / "reports/role_d/contextual_iqr/role_d_functional_position_iqr.csv"
TABLES = ROOT / "reports/role_d/contextual_iqr"
MIN_PDB = 5
N_PERMUTATIONS = 10_000
RNG_SEED = 20261022


def context_stratum(data: pd.DataFrame) -> pd.Series:
    element = data["canonical_elements"].fillna("").astype(str)
    inter = "INTER:" + data["inter_helix_flank_before"].fillna("NONE") + "-" + data["inter_helix_flank_after"].fillna("NONE")
    return element.where(element.ne(""), inter)


def support_band(n: pd.Series) -> pd.Series:
    return pd.cut(n, [MIN_PDB - 1, 9, 19, np.inf], labels=["5_9", "10_19", "20_PLUS"]).astype(str)


def receptor_effect(data: pd.DataFrame) -> float:
    effects = []
    for _, stratum in data.groupby("match_stratum", sort=False):
        site = stratum.loc[stratum.is_lbp, "iqr_z"]
        control = stratum.loc[~stratum.is_lbp, "iqr_z"]
        if len(site) and len(control) >= 3:
            effects.append(site.median() - control.median())
    return float(np.median(effects)) if effects else np.nan


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    profile = pd.read_csv(PROFILE)
    functional = pd.read_csv(FUNCTIONAL)
    functional_keys = set(zip(functional.uniprot_id, functional.canonical_uniprot_position))
    lbp_keys = set(zip(functional.loc[functional.functional_site.eq("LIGAND_BINDING_POCKET"), "uniprot_id"],
                       functional.loc[functional.functional_site.eq("LIGAND_BINDING_POCKET"), "canonical_uniprot_position"]))
    data = profile.loc[profile.iqr_z.notna() & profile.n_pdb_entries_observed_at_position.ge(MIN_PDB)].copy()
    key = list(zip(data.uniprot_id, data.canonical_uniprot_position))
    data["is_lbp"] = [item in lbp_keys for item in key]
    data["is_any_functional"] = [item in functional_keys for item in key]
    data = data.loc[data.is_lbp | ~data.is_any_functional].copy()
    data["context_stratum"] = context_stratum(data)
    data["support_band"] = support_band(data.n_pdb_entries_observed_at_position)
    data["match_stratum"] = data.context_stratum + "|" + data.support_band

    rows = []
    used = []
    for receptor, receptor_data in data.groupby("uniprot_id", sort=True):
        for stratum, group in receptor_data.groupby("match_stratum", sort=True):
            n_site, n_control = int(group.is_lbp.sum()), int((~group.is_lbp).sum())
            if n_site and n_control >= 3:
                effect = group.loc[group.is_lbp, "iqr_z"].median() - group.loc[~group.is_lbp, "iqr_z"].median()
                rows.append({"uniprot_id": receptor, "nr_code": group.nr_code.iloc[0], "common_name": group.common_name.iloc[0],
                             "match_stratum": stratum, "context_stratum": group.context_stratum.iloc[0],
                             "support_band": group.support_band.iloc[0], "n_lbp_positions": n_site,
                             "n_control_positions": n_control, "lbp_median_iqr_z": group.loc[group.is_lbp, "iqr_z"].median(),
                             "control_median_iqr_z": group.loc[~group.is_lbp, "iqr_z"].median(), "matched_difference_iqr_z": effect})
                used.append(group)
    strata = pd.DataFrame(rows)
    if strata.empty:
        return strata, pd.DataFrame()
    receptor = strata.groupby(["uniprot_id", "nr_code", "common_name"], as_index=False).agg(
        n_matched_strata=("match_stratum", "size"),
        matched_difference_iqr_z=("matched_difference_iqr_z", "median"),
    )

    # Restricted null: shuffle LBP labels within every receptor × context × support stratum.
    rng = np.random.default_rng(RNG_SEED)
    # Restrict the null to strata eligible under the observed design. This preserves
    # the pre-specified matched strata rather than allowing a shuffle to create a
    # new comparison in a stratum that had no observed LBP residue.
    matched_data = pd.concat(used, ignore_index=True)
    permutation_design = []
    for _, receptor_data in matched_data.groupby("uniprot_id", sort=True):
        strata_design = []
        for _, group in receptor_data.groupby("match_stratum", sort=False):
            strata_design.append((group.iqr_z.to_numpy(), int(group.is_lbp.sum())))
        permutation_design.append(strata_design)
    simulated_receptor_effects = []
    for strata_design in permutation_design:
        stratum_effects = []
        for values, n_site in strata_design:
            # One random ordering per permutation; its first n_site values are
            # the restricted shuffled site label set.
            order = np.argsort(rng.random((N_PERMUTATIONS, len(values))), axis=1)
            shuffled = values[order]
            stratum_effects.append(np.median(shuffled[:, :n_site], axis=1) - np.median(shuffled[:, n_site:], axis=1))
        simulated_receptor_effects.append(np.median(np.column_stack(stratum_effects), axis=1))
    null = np.median(np.column_stack(simulated_receptor_effects), axis=1)
    observed = receptor.matched_difference_iqr_z.median()
    p_value = (1 + np.count_nonzero(np.abs(null) >= abs(observed))) / (N_PERMUTATIONS + 1)
    receptor_effects = receptor.matched_difference_iqr_z.to_numpy()
    bootstrap = np.median(
        receptor_effects[rng.integers(0, len(receptor_effects), size=(N_PERMUTATIONS, len(receptor_effects)))], axis=1
    )
    summary = pd.DataFrame([{
        "analysis": "LBP_VS_NONFUNCTIONAL_EXACT_CONTEXT_SUPPORT_MATCHED",
        "n_receptors": len(receptor), "n_matched_strata": len(strata),
        "receptor_median_matched_difference_iqr_z": observed,
        "bootstrap_ci95_low": np.quantile(bootstrap, 0.025),
        "bootstrap_ci95_high": np.quantile(bootstrap, 0.975),
        "receptors_with_positive_difference": int((receptor.matched_difference_iqr_z > 0).sum()),
        "restricted_permutation_p_two_sided": p_value, "n_permutations": N_PERMUTATIONS,
        "random_seed": RNG_SEED,
        "conservative_status": "INCONCLUSIVE_BOOTSTRAP_CI_CROSSES_ZERO",
        "interpretation": "POSITION_LEVEL_ASSOCIATION_BETWEEN_PDB_CRYSTALLOGRAPHIC_DISORDER_HETEROGENEITY_NOT_SOLUTION_DYNAMICS",
    }])
    return strata, pd.concat([receptor.assign(record_type="RECEPTOR"), summary.assign(record_type="SUMMARY")], sort=False)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    strata, results = build()
    strata.to_csv(TABLES / "role_d_lbp_matched_iqr_strata.csv", index=False, lineterminator="\n")
    results.to_csv(TABLES / "role_d_lbp_matched_iqr_results.csv", index=False, lineterminator="\n")
    print(f"Wrote {len(strata)} matched strata and {len(results)} matched-analysis records.")


if __name__ == "__main__":
    main()
