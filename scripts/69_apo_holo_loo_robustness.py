"""
Part 13 of the APO vs HOLO flexibility stage: leave-one-PDB-out (LOO)
robustness for every significant or nominally interesting (raw p<0.05)
receptor x segment result -- to check whether an effect is driven by a
single deposition. Deterministic (no randomness; every fold just drops
one PDB entry from the larger group, or from both if both groups have
enough members).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def main():
    effects = pd.read_parquet(ANALYSIS / "apo_holo_segment_effects.parquet")
    interesting = effects[effects["raw_p"] < 0.05].copy()
    print(f"{len(interesting)} nominal-or-better receptor x segment results to LOO-check")

    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")
    features = features[features["ligand_state"].isin(["APO", "HOLO"])]

    rows = []
    for _, r in interesting.iterrows():
        uid, seg = r["uniprot_id"], r["structural_segment"]
        g = features[(features["uniprot_id"] == uid) & (features["structural_segment"] == seg)]
        apo = g.loc[g["ligand_state"] == "APO", ["pdb_id", "median_z"]].dropna()
        holo = g.loc[g["ligand_state"] == "HOLO", ["pdb_id", "median_z"]].dropna()

        full_delta = float(holo["median_z"].median() - apo["median_z"].median())
        full_sign = np.sign(full_delta)

        loo_deltas = []
        loo_ps = []
        # drop one APO entry at a time
        for pdb_id in apo["pdb_id"]:
            a = apo.loc[apo["pdb_id"] != pdb_id, "median_z"]
            h = holo["median_z"]
            if len(a) < 2 or len(h) < 2:
                continue
            loo_deltas.append(float(h.median() - a.median()))
            _, p = stats.mannwhitneyu(h, a, alternative="two-sided")
            loo_ps.append(p)
        # drop one HOLO entry at a time
        for pdb_id in holo["pdb_id"]:
            a = apo["median_z"]
            h = holo.loc[holo["pdb_id"] != pdb_id, "median_z"]
            if len(a) < 2 or len(h) < 2:
                continue
            loo_deltas.append(float(h.median() - a.median()))
            _, p = stats.mannwhitneyu(h, a, alternative="two-sided")
            loo_ps.append(p)

        loo_deltas = np.array(loo_deltas)
        loo_ps = np.array(loo_ps)
        frac_same_sign = float((np.sign(loo_deltas) == full_sign).mean()) if len(loo_deltas) else np.nan
        frac_p_lt_05 = float((loo_ps < 0.05).mean()) if len(loo_ps) else np.nan

        rows.append({
            "uniprot_id": uid, "nr_code": r["nr_code"], "structural_segment": seg,
            "full_delta_median_z": full_delta, "full_raw_p": r["raw_p"],
            "significant_within_receptor": bool(r["significant_within_receptor"]),
            "n_loo_iterations": len(loo_deltas),
            "fraction_loo_same_direction": frac_same_sign,
            "fraction_loo_p_lt_05": frac_p_lt_05,
            "loo_delta_min": float(loo_deltas.min()) if len(loo_deltas) else np.nan,
            "loo_delta_max": float(loo_deltas.max()) if len(loo_deltas) else np.nan,
        })

    loo = pd.DataFrame(rows).sort_values("full_raw_p").reset_index(drop=True)
    out = TABLES / "apo_holo_loo_robustness.csv"
    loo.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(loo)} rows)")
    print(loo[["uniprot_id", "structural_segment", "significant_within_receptor", "full_delta_median_z",
               "fraction_loo_same_direction", "fraction_loo_p_lt_05"]].to_string(index=False))

    fragile = loo[(loo["significant_within_receptor"]) & (loo["fraction_loo_p_lt_05"] < 1.0)]
    if len(fragile):
        print(f"\n{len(fragile)} within-receptor-FDR-significant result(s) do NOT retain raw p<0.05 "
              "in every LOO fold -- flagged as at least partly driven by specific depositions:")
        print(fragile[["uniprot_id", "structural_segment", "fraction_loo_p_lt_05"]].to_string(index=False))


if __name__ == "__main__":
    main()
