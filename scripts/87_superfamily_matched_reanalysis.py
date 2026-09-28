"""
Part 7 of the functional-site specificity/validation stage: repeat the
superfamily site-class test (one receptor-level value per receptor) using
the matched-background deltas instead of the original whole-LBD
background, restricted to HIGH/MODERATE-confidence RESIDUE-level
definitions only (matched-background tests are already RESIDUE-level-only
by construction; here we additionally drop any LOW-confidence residue
definition). Low-confidence segment proxies are never mixed into this
confirmatory analysis -- they remain descriptive/sensitivity evidence
only (Parts 2, 4).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 0
CONFIDENT_TIERS = {"HIGH", "MODERATE"}


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


def bootstrap_ci(values: np.ndarray, seed=BOOTSTRAP_SEED, n_boot=BOOTSTRAP_N):
    if len(values) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.array([np.median(rng.choice(values, size=len(values), replace=True)) for _ in range(n_boot)])
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def main():
    matched = pd.read_csv(TABLES / "functional_site_matched_tests.csv")
    receptor_site = pd.read_parquet(ANALYSIS / "receptor_functional_site_flexibility.parquet")

    matched = matched.merge(
        receptor_site[["uniprot_id", "functional_site", "site_definition_confidence"]],
        on=["uniprot_id", "functional_site"], how="left",
    )
    confident = matched[matched["site_definition_confidence"].isin(CONFIDENT_TIERS)]
    print(f"{len(confident)}/{len(matched)} matched-background tests have HIGH/MODERATE-confidence "
          "RESIDUE-level definitions -- only these enter the confirmatory superfamily reanalysis")

    rows = []
    for site, g in confident.groupby("functional_site"):
        deltas = g["median_paired_delta"].values
        n = len(deltas)
        if n < 3:
            print(f"  {site}: only {n} confident receptor(s) -- INSUFFICIENT for a superfamily reanalysis")
            rows.append({
                "functional_site": site, "n_receptors": n, "median_matched_delta": np.nan,
                "bootstrap_ci_low": np.nan, "bootstrap_ci_high": np.nan, "raw_p": np.nan, "status": "INSUFFICIENT_N",
            })
            continue
        ci_low, ci_high = bootstrap_ci(deltas)
        if np.all(deltas == deltas[0]):
            stat, p = np.nan, np.nan
        else:
            stat, p = stats.wilcoxon(deltas, alternative="two-sided")
        rows.append({
            "functional_site": site, "n_receptors": n, "median_matched_delta": float(np.median(deltas)),
            "bootstrap_ci_low": ci_low, "bootstrap_ci_high": ci_high, "raw_p": p, "status": "TESTED",
        })

    result = pd.DataFrame(rows)
    tested = result[result["status"] == "TESTED"]
    if len(tested):
        q = benjamini_hochberg(tested["raw_p"].values)
        result.loc[result["status"] == "TESTED", "BH_q"] = q
        result.loc[result["status"] == "TESTED", "significant_fdr_05"] = q <= 0.05
    else:
        result["BH_q"] = np.nan
        result["significant_fdr_05"] = False

    out = TABLES / "functional_site_superfamily_matched_reanalysis.csv"
    result.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)}")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
