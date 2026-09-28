"""
Part 9 of the current functional-site flexibility stage: aggregate PDB
entries to the receptor level (the approved instance-Z -> PDB-entry-median
-> receptor-summary hierarchy's final step). PDB-level rows are retained
separately (pdb_functional_site_flexibility.parquet) for within-receptor
uncertainty/robustness; they are never treated as independent superfamily-
level replicates here.
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"

BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 0


def bootstrap_ci(values: np.ndarray, seed=BOOTSTRAP_SEED, n_boot=BOOTSTRAP_N):
    if len(values) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.array([np.median(rng.choice(values, size=len(values), replace=True)) for _ in range(n_boot)])
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def main():
    pdb_site = pd.read_parquet(ANALYSIS / "pdb_functional_site_flexibility.parquet")

    rows = []
    for (uid, site), g in pdb_site.groupby(["uniprot_id", "functional_site"]):
        nr_code = g["nr_code"].iloc[0]
        site_vals = g["median_site_z"].values
        delta_vals = g["delta_site_vs_background"].values
        ci_site_low, ci_site_high = bootstrap_ci(site_vals)
        ci_delta_low, ci_delta_high = bootstrap_ci(delta_vals)

        rows.append({
            "uniprot_id": uid, "nr_code": nr_code, "functional_site": site,
            "n_pdb_entries": len(g),
            "median_site_z": float(np.median(site_vals)), "mean_site_z": float(np.mean(site_vals)),
            "iqr_site_z": float(np.percentile(site_vals, 75) - np.percentile(site_vals, 25)),
            "median_background_z": float(g["comparison_non_site_median_z"].median()),
            "median_delta_site_vs_background": float(np.median(delta_vals)),
            "bootstrap_95ci_site_median_low": ci_site_low, "bootstrap_95ci_site_median_high": ci_site_high,
            "bootstrap_95ci_delta_low": ci_delta_low, "bootstrap_95ci_delta_high": ci_delta_high,
            "site_definition_confidence": g["site_definition_confidence"].mode().iloc[0],
            "definition_level": g["definition_level"].mode().iloc[0],
            "median_site_coverage": float(g["site_coverage"].median()),
        })

    receptor_site = pd.DataFrame(rows).sort_values(["functional_site", "uniprot_id"]).reset_index(drop=True)
    out = ANALYSIS / "receptor_functional_site_flexibility.parquet"
    receptor_site.to_parquet(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(receptor_site)} rows)")
    print("\nReceptors per site class:")
    print(receptor_site["functional_site"].value_counts().to_string())
    print("\nMedian delta (site - background) by site class:")
    print(receptor_site.groupby("functional_site")["median_delta_site_vs_background"].median().to_string())


if __name__ == "__main__":
    main()
