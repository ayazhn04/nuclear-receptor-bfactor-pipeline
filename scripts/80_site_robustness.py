"""
Part 15 of the current functional-site flexibility stage: robustness.

(a) Leave-one-PDB-out, WITHIN receptor, for every within-receptor-FDR-
    significant site-vs-background result (Part 10) -- does the paired
    Wilcoxon result depend on one deposition?
(b) Leave-one-receptor-out for the superfamily site-class results
    (Part 11) -- does a site class's overall significance depend on one
    receptor?

Flags results dominated by one receptor, one PDB series, low site
coverage, or a low-confidence site definition.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def rank_biserial_effect_size(diffs: np.ndarray) -> float:
    nz = diffs[diffs != 0]
    if len(nz) == 0:
        return 0.0
    ranks = stats.rankdata(np.abs(nz))
    w_pos = ranks[nz > 0].sum()
    w_neg = ranks[nz < 0].sum()
    return float((w_pos - w_neg) / (w_pos + w_neg))


def loo_pdb_within_receptor():
    tests = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv")
    pdb_site = pd.read_parquet(ANALYSIS / "pdb_functional_site_flexibility.parquet")
    sig = tests[tests["significant_within_receptor"]]

    rows = []
    for _, r in sig.iterrows():
        uid, site = r["uniprot_id"], r["functional_site"]
        g = pdb_site[(pdb_site["uniprot_id"] == uid) & (pdb_site["functional_site"] == site)]
        site_vals, bg_vals = g["median_site_z"].values, g["comparison_non_site_median_z"].values
        full_sign = np.sign(r["median_paired_delta"])

        loo_p, loo_same_sign = [], []
        for i in range(len(g)):
            sv = np.delete(site_vals, i)
            bv = np.delete(bg_vals, i)
            diffs = sv - bv
            if np.all(diffs == 0) or len(sv) < 3:
                continue
            _, p = stats.wilcoxon(sv, bv, alternative="two-sided", zero_method="wilcox")
            loo_p.append(p)
            loo_same_sign.append(np.sign(np.median(diffs)) == full_sign)

        rows.append({
            "uniprot_id": uid, "nr_code": r["nr_code"], "functional_site": site,
            "n_pairs": r["n_pairs"], "full_raw_p": r["raw_p"],
            "n_loo_iterations": len(loo_p),
            "fraction_loo_p_lt_05": float(np.mean([p < 0.05 for p in loo_p])) if loo_p else np.nan,
            "fraction_loo_same_direction": float(np.mean(loo_same_sign)) if loo_same_sign else np.nan,
        })
    return pd.DataFrame(rows)


def loo_receptor_superfamily():
    receptor_site = pd.read_parquet(ANALYSIS / "receptor_functional_site_flexibility.parquet")
    superfamily = pd.read_csv(TABLES / "current_functional_site_superfamily_tests.csv")

    rows = []
    for _, r in superfamily.iterrows():
        site = r["functional_site"]
        g = receptor_site[receptor_site["functional_site"] == site]
        deltas = g["median_delta_site_vs_background"].values
        uids = g["uniprot_id"].values
        full_sign = np.sign(r["median_delta_vs_background"])

        loo_p, loo_same_sign = [], []
        for i in range(len(deltas)):
            d = np.delete(deltas, i)
            if len(d) < 2 or np.all(d == d[0]):
                continue
            _, p = stats.wilcoxon(d, alternative="two-sided")
            loo_p.append(p)
            loo_same_sign.append(np.sign(np.median(d)) == full_sign)

        # single-receptor dominance: does removing any one receptor flip significance?
        flips_significance = False
        if not pd.isna(r["raw_p"]):
            was_sig_pattern = [p >= 0.05 for p in loo_p] if r["raw_p"] < 0.05 else [p < 0.05 for p in loo_p]
            flips_significance = any(was_sig_pattern)

        rows.append({
            "functional_site": site, "n_receptors": len(deltas), "full_raw_p": r["raw_p"],
            "fraction_loo_p_lt_05": float(np.mean([p < 0.05 for p in loo_p])) if loo_p else np.nan,
            "fraction_loo_same_direction": float(np.mean(loo_same_sign)) if loo_same_sign else np.nan,
            "significance_flips_on_some_loo_fold": flips_significance,
        })
    return pd.DataFrame(rows)


def flag_low_confidence_or_coverage():
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")
    receptor_site = pd.read_parquet(ANALYSIS / "receptor_functional_site_flexibility.parquet")
    flags = receptor_site.copy()
    flags["low_confidence_definition"] = flags["site_definition_confidence"].isin(["LOW", "MODERATE"])
    flags["low_coverage"] = flags["median_site_coverage"] < 0.5
    flags["segment_level_only"] = flags["definition_level"] == "SEGMENT"
    return flags[["uniprot_id", "nr_code", "functional_site", "site_definition_confidence",
                  "definition_level", "median_site_coverage", "low_confidence_definition",
                  "low_coverage", "segment_level_only"]]


def main():
    loo_pdb = loo_pdb_within_receptor()
    out1 = TABLES / "functional_site_loo_pdb_robustness.csv"
    loo_pdb.to_csv(out1, index=False)
    print(f"Wrote {out1.relative_to(PROJECT_ROOT)} ({len(loo_pdb)} within-receptor-significant results checked)")
    fragile_pdb = loo_pdb[loo_pdb["fraction_loo_p_lt_05"] < 1.0]
    print(f"{len(fragile_pdb)}/{len(loo_pdb)} do not retain p<0.05 in every leave-one-PDB-out fold")

    loo_receptor = loo_receptor_superfamily()
    out2 = TABLES / "functional_site_loo_receptor_robustness.csv"
    loo_receptor.to_csv(out2, index=False)
    print(f"\nWrote {out2.relative_to(PROJECT_ROOT)}")
    print(loo_receptor.to_string(index=False))
    dominated = loo_receptor[loo_receptor["significance_flips_on_some_loo_fold"]]
    if len(dominated):
        print(f"\n{len(dominated)} site class(es) have significance that flips when a single receptor is "
              f"removed -- flagged as potentially receptor-dominated: {dominated['functional_site'].tolist()}")

    flags = flag_low_confidence_or_coverage()
    out3 = TABLES / "functional_site_confidence_coverage_flags.csv"
    flags.to_csv(out3, index=False)
    print(f"\nWrote {out3.relative_to(PROJECT_ROOT)} ({len(flags)} rows)")
    print(f"{flags['low_confidence_definition'].sum()} rows flagged low/moderate-confidence definitions; "
          f"{flags['low_coverage'].sum()} rows flagged low (<50%) median site coverage; "
          f"{flags['segment_level_only'].sum()} rows are segment-level-only definitions.")


if __name__ == "__main__":
    main()
