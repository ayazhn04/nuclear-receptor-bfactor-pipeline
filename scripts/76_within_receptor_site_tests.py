"""
Part 10 of the current functional-site flexibility stage: within-receptor
paired site-vs-background test. Each PDB entry contributes both a site
value and a background value (median_site_z, comparison_non_site_median_z)
-- a natural pairing, tested with the paired non-parametric Wilcoxon
signed-rank test rather than treating site and background as independent
samples. Individual residues are never treated as independent observations
here; the PDB entry is the paired unit.

Minimum >=5 paired PDB entries per receptor x site (fixed before testing).
BH-FDR applied within each receptor across its own tested site classes
(PRIMARY criterion); a secondary global BH-FDR is also reported.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

MIN_PAIRS = 5
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
    """Matched-pairs rank-biserial correlation for Wilcoxon signed-rank:
    r = (W+ - W-) / (W+ + W-), using signed ranks of the nonzero differences.
    Ranges -1 to +1; sign matches whether site > background typically."""
    nz = diffs[diffs != 0]
    if len(nz) == 0:
        return 0.0
    ranks = stats.rankdata(np.abs(nz))
    w_pos = ranks[nz > 0].sum()
    w_neg = ranks[nz < 0].sum()
    return float((w_pos - w_neg) / (w_pos + w_neg))


def bootstrap_delta_ci(diffs: np.ndarray, seed=BOOTSTRAP_SEED, n_boot=BOOTSTRAP_N):
    rng = np.random.default_rng(seed)
    boot = np.array([np.median(rng.choice(diffs, size=len(diffs), replace=True)) for _ in range(n_boot)])
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def main():
    pdb_site = pd.read_parquet(ANALYSIS / "pdb_functional_site_flexibility.parquet")

    rows = []
    for (uid, site), g in pdb_site.groupby(["uniprot_id", "functional_site"]):
        if len(g) < MIN_PAIRS:
            continue
        nr_code = g["nr_code"].iloc[0]
        site_vals = g["median_site_z"].values
        bg_vals = g["comparison_non_site_median_z"].values
        diffs = site_vals - bg_vals

        if np.all(diffs == 0):
            continue  # Wilcoxon undefined for all-zero differences
        stat, p = stats.wilcoxon(site_vals, bg_vals, alternative="two-sided", zero_method="wilcox")
        effect = rank_biserial_effect_size(diffs)
        ci_low, ci_high = bootstrap_delta_ci(diffs)

        rows.append({
            "uniprot_id": uid, "nr_code": nr_code, "functional_site": site,
            "n_pairs": len(g), "median_site_z": float(np.median(site_vals)),
            "median_background_z": float(np.median(bg_vals)), "median_paired_delta": float(np.median(diffs)),
            "bootstrap_ci_low": ci_low, "bootstrap_ci_high": ci_high,
            "wilcoxon_statistic": float(stat), "effect_size_rank_biserial": effect, "raw_p": float(p),
        })

    tests = pd.DataFrame(rows)
    tests["q_within_receptor"] = np.nan
    for uid, idx in tests.groupby("uniprot_id").groups.items():
        tests.loc[idx, "q_within_receptor"] = benjamini_hochberg(tests.loc[idx, "raw_p"].values)
    tests["significant_within_receptor"] = tests["q_within_receptor"] <= 0.05

    tests["q_global"] = benjamini_hochberg(tests["raw_p"].values)
    tests["significant_global"] = tests["q_global"] <= 0.05

    tests = tests.sort_values(["uniprot_id", "raw_p"]).reset_index(drop=True)

    out = TABLES / "current_functional_site_within_receptor_tests.csv"
    tests.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(tests)} receptor x site tests, "
          f"min {MIN_PAIRS} paired PDB entries)")

    n_nominal = (tests["raw_p"] < 0.05).sum()
    n_sig_within = tests["significant_within_receptor"].sum()
    n_sig_global = tests["significant_global"].sum()
    print(f"{n_nominal} nominal (p<0.05); {n_sig_within} significant within-receptor FDR; "
          f"{n_sig_global} significant global FDR")

    if n_sig_within:
        sig = tests[tests["significant_within_receptor"]]
        print("\nWithin-receptor-FDR-significant site-vs-background results:")
        print(sig[["nr_code", "functional_site", "n_pairs", "median_paired_delta",
                    "effect_size_rank_biserial", "raw_p", "q_within_receptor", "q_global"]].to_string(index=False))


if __name__ == "__main__":
    main()
