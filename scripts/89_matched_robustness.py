"""
Part 10 of the functional-site specificity/validation stage: robustness
for the matched-background results, using the same fixed-seed/
deterministic conventions as every prior robustness pass in this project
(no randomness beyond the already-fixed bootstrap seed=0 used upstream).

(a) Leave-one-PDB-out within receptor, for every within-receptor-
    significant matched-background result.
(b) Leave-one-receptor-out for the confirmatory superfamily matched
    reanalysis (Part 7).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def loo_pdb():
    tests = pd.read_csv(TABLES / "functional_site_matched_tests.csv")
    pdb_matched = pd.read_parquet(ANALYSIS / "pdb_functional_site_matched_background.parquet")
    sig = tests[tests["significant_within_receptor"]]

    rows = []
    for _, r in sig.iterrows():
        uid, site = r["uniprot_id"], r["functional_site"]
        g = pdb_matched[(pdb_matched["uniprot_id"] == uid) & (pdb_matched["functional_site"] == site)]
        site_vals, bg_vals = g["median_site_z"].values, g["median_matched_background_z"].values
        full_sign = np.sign(r["median_paired_delta"])

        loo_p, loo_same_sign = [], []
        for i in range(len(g)):
            sv, bv = np.delete(site_vals, i), np.delete(bg_vals, i)
            diffs = sv - bv
            if np.all(diffs == 0) or len(sv) < 3:
                continue
            _, p = stats.wilcoxon(sv, bv, alternative="two-sided", zero_method="wilcox")
            loo_p.append(p)
            loo_same_sign.append(np.sign(np.median(diffs)) == full_sign)

        rows.append({
            "uniprot_id": uid, "nr_code": r["nr_code"], "functional_site": site,
            "n_pairs": r["n_pairs"], "full_raw_p": r["raw_p"],
            "fraction_loo_p_lt_05": float(np.mean([p < 0.05 for p in loo_p])) if loo_p else np.nan,
            "fraction_loo_same_direction": float(np.mean(loo_same_sign)) if loo_same_sign else np.nan,
        })
    return pd.DataFrame(rows)


def loo_receptor():
    superfamily = pd.read_csv(TABLES / "functional_site_superfamily_matched_reanalysis.csv")
    matched = pd.read_csv(TABLES / "functional_site_matched_tests.csv")
    receptor_site = pd.read_parquet(ANALYSIS / "receptor_functional_site_flexibility.parquet")
    matched = matched.merge(
        receptor_site[["uniprot_id", "functional_site", "site_definition_confidence"]],
        on=["uniprot_id", "functional_site"], how="left",
    )
    confident = matched[matched["site_definition_confidence"].isin(["HIGH", "MODERATE"])]

    rows = []
    for _, r in superfamily[superfamily["status"] == "TESTED"].iterrows():
        site = r["functional_site"]
        g = confident[confident["functional_site"] == site]
        deltas = g["median_paired_delta"].values
        full_sign = np.sign(r["median_matched_delta"])

        loo_p, loo_same_sign = [], []
        for i in range(len(deltas)):
            d = np.delete(deltas, i)
            if len(d) < 2 or np.all(d == d[0]):
                continue
            _, p = stats.wilcoxon(d, alternative="two-sided")
            loo_p.append(p)
            loo_same_sign.append(np.sign(np.median(d)) == full_sign)

        flips = False
        if not pd.isna(r["raw_p"]):
            was_sig_pattern = [p >= 0.05 for p in loo_p] if r["raw_p"] < 0.05 else [p < 0.05 for p in loo_p]
            flips = any(was_sig_pattern)

        rows.append({
            "functional_site": site, "n_receptors": len(deltas), "full_raw_p": r["raw_p"],
            "fraction_loo_p_lt_05": float(np.mean([p < 0.05 for p in loo_p])) if loo_p else np.nan,
            "fraction_loo_same_direction": float(np.mean(loo_same_sign)) if loo_same_sign else np.nan,
            "significance_flips_on_some_loo_fold": flips,
        })
    return pd.DataFrame(rows)


def main():
    loo_pdb_df = loo_pdb()
    out1 = TABLES / "functional_site_matched_loo_pdb_robustness.csv"
    loo_pdb_df.to_csv(out1, index=False)
    print(f"Wrote {out1.relative_to(PROJECT_ROOT)} ({len(loo_pdb_df)} within-receptor-significant "
          "matched-background results checked)")
    fragile = loo_pdb_df[loo_pdb_df["fraction_loo_p_lt_05"] < 1.0]
    print(f"{len(fragile)}/{len(loo_pdb_df)} do not retain p<0.05 in every leave-one-PDB-out fold")

    loo_receptor_df = loo_receptor()
    out2 = TABLES / "functional_site_matched_loo_receptor_robustness.csv"
    loo_receptor_df.to_csv(out2, index=False)
    print(f"\nWrote {out2.relative_to(PROJECT_ROOT)}")
    print(loo_receptor_df.to_string(index=False))


if __name__ == "__main__":
    main()
