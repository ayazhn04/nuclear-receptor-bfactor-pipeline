"""
Parts 11-12 of the current functional-site flexibility stage.

Part 11: are particular functional-site classes consistently more flexible
or rigid across receptors? One receptor-level value per receptor/site
(receptor_functional_site_flexibility.parquet's median_delta_site_vs_
background) -- never PDB entries as independent superfamily-level
replicates. Tested against zero with a one-sample Wilcoxon signed-rank
test per site class, BH-FDR corrected across the 4 site classes. Receptors
where that site is UNRESOLVED are simply absent from that site's
receptor-level table already (Part 9) and are excluded here too.

Part 12: pairwise site-class comparison using receptors with BOTH sites
measured (a paired/repeated framework -- receptors are never treated as
contributing independent observations to both classes as if unpaired).
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 0


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def rank_biserial_effect_size(diffs: np.ndarray) -> float:
    nz = diffs[diffs != 0]
    if len(nz) == 0:
        return 0.0
    ranks = stats.rankdata(np.abs(nz))
    w_pos = ranks[nz > 0].sum()
    w_neg = ranks[nz < 0].sum()
    return float((w_pos - w_neg) / (w_pos + w_neg))


def bootstrap_ci(values: np.ndarray, seed=BOOTSTRAP_SEED, n_boot=BOOTSTRAP_N):
    if len(values) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.array([np.median(rng.choice(values, size=len(values), replace=True)) for _ in range(n_boot)])
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def main():
    receptor_site = pd.read_parquet(ANALYSIS / "receptor_functional_site_flexibility.parquet")

    # ---- Part 11: one-sample test per site class ----
    rows = []
    for site, g in receptor_site.groupby("functional_site"):
        deltas = g["median_delta_site_vs_background"].values
        n = len(deltas)
        ci_low, ci_high = bootstrap_ci(deltas)
        if n >= 2 and not np.all(deltas == deltas[0]):
            stat, p = stats.wilcoxon(deltas, alternative="two-sided")
        else:
            stat, p = np.nan, np.nan
        rows.append({
            "functional_site": site, "n_receptors": n,
            "median_site_z": float(g["median_site_z"].median()),
            "median_delta_vs_background": float(np.median(deltas)),
            "bootstrap_ci_low": ci_low, "bootstrap_ci_high": ci_high,
            "wilcoxon_statistic": stat, "raw_p": p,
            "fraction_receptors_negative_delta": float((deltas < 0).mean()),
            "fraction_receptors_positive_delta": float((deltas > 0).mean()),
        })
    site_class = pd.DataFrame(rows)
    site_class["fdr_q"] = benjamini_hochberg(site_class["raw_p"].fillna(1.0).values)
    site_class.loc[site_class["raw_p"].isna(), "fdr_q"] = np.nan
    site_class["significant_fdr_05"] = site_class["fdr_q"] <= 0.05

    out1 = TABLES / "current_functional_site_superfamily_tests.csv"
    site_class.to_csv(out1, index=False)
    print(f"Wrote {out1.relative_to(PROJECT_ROOT)} ({len(site_class)} site classes)")
    print(site_class[["functional_site", "n_receptors", "median_delta_vs_background",
                       "raw_p", "fdr_q", "significant_fdr_05"]].to_string(index=False))

    # ---- Part 12: paired site-class comparison ----
    pivot = receptor_site.pivot(index="uniprot_id", columns="functional_site", values="median_delta_site_vs_background")
    sites = sorted(pivot.columns)
    pair_rows = []
    for a, b in itertools.combinations(sites, 2):
        paired = pivot[[a, b]].dropna()
        if len(paired) < 5:
            continue
        diffs = (paired[a] - paired[b]).values
        if np.all(diffs == 0):
            continue
        stat, p = stats.wilcoxon(paired[a], paired[b], alternative="two-sided")
        effect = rank_biserial_effect_size(diffs)
        pair_rows.append({
            "site_a": a, "site_b": b, "n_receptors_paired": len(paired),
            "median_delta_a": float(paired[a].median()), "median_delta_b": float(paired[b].median()),
            "median_paired_difference_a_minus_b": float(np.median(diffs)),
            "effect_size_rank_biserial": effect, "raw_p": float(p),
        })
    site_pairs = pd.DataFrame(pair_rows)
    if len(site_pairs):
        site_pairs["fdr_q"] = benjamini_hochberg(site_pairs["raw_p"].values)
        site_pairs["significant_fdr_05"] = site_pairs["fdr_q"] <= 0.05
    out2 = TABLES / "current_functional_site_class_comparison.csv"
    site_pairs.to_csv(out2, index=False)
    print(f"\nWrote {out2.relative_to(PROJECT_ROOT)} ({len(site_pairs)} pairwise site-class comparisons)")
    if len(site_pairs):
        print(site_pairs.to_string(index=False))


if __name__ == "__main__":
    main()
