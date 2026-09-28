"""
Part 13 of the subfamily flexibility stage: repeat the PRIMARY omnibus
Kruskal-Wallis analysis using the receptor MEAN segment summary
(`mean_of_entry_mean_z` -- the same aggregation underlying
receptor_segment_feature_matrix_mean.parquet: entry-level MEAN, then MEAN
across entries) instead of the median pipeline, as a sensitivity check.
Reuses the exact same statistical machinery as script 56 -- this does not
create a separate biological narrative, only reports agreement/
disagreement with the median-based primary analysis.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

spec = importlib.util.spec_from_file_location("kruskal_mod", PROJECT_ROOT / "scripts" / "56_subfamily_kruskal_omnibus.py")
kruskal_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kruskal_mod)


def main():
    profiles = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet").copy()
    meta = pd.read_csv(ANALYSIS / "subfamily_receptor_metadata.csv")

    # swap the outcome column: run_omnibus_for_scheme always reads
    # "median_of_entry_median_z" -- rename columns so the MEAN column takes
    # that slot for this sensitivity run only (a local copy, upstream file untouched).
    profiles_mean = profiles.rename(columns={
        "median_of_entry_median_z": "_median_of_entry_median_z_UNUSED",
        "mean_of_entry_mean_z": "median_of_entry_median_z",
    })

    all_results = []
    for scheme_name, scheme_col in kruskal_mod.SCHEME_COLUMNS.items():
        print(f"--- MEAN-sensitivity omnibus Kruskal-Wallis: {scheme_name} ---")
        result = kruskal_mod.run_omnibus_for_scheme(profiles_mean, meta, scheme_col, scheme_name)
        result = kruskal_mod.apply_fdr(result)
        result["outcome_summary"] = "mean_of_entry_mean_z"
        all_results.append(result)
        n_tested = (result["test_status"] == "TESTED").sum()
        n_sig_fdr = result["significant_fdr_05"].sum()
        print(f"  {n_tested} segments tested; {n_sig_fdr} significant after BH-FDR")

    combined = pd.concat(all_results, ignore_index=True)
    combined = combined.sort_values(["grouping_scheme", "raw_p"], na_position="last").reset_index(drop=True)
    out = TABLES / "subfamily_segment_kruskal_mean_sensitivity.csv"
    combined.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(combined)} rows)")

    # agreement/disagreement with the primary (median-based) analysis
    primary = pd.read_csv(TABLES / "subfamily_segment_kruskal.csv")
    merged = primary.merge(
        combined[["grouping_scheme", "structural_segment", "raw_p", "significant_fdr_05"]],
        on=["grouping_scheme", "structural_segment"], suffixes=("_median", "_mean"),
    )
    agreement = merged["significant_fdr_05_median"] == merged["significant_fdr_05_mean"]
    print(f"\nAgreement between median-based and mean-based FDR-significance calls: "
          f"{agreement.sum()}/{len(agreement)} segment/scheme pairs ({agreement.mean():.1%})")
    both_none_significant = (~merged["significant_fdr_05_median"]).all() and (~merged["significant_fdr_05_mean"]).all()
    if both_none_significant:
        print("Both median-based and mean-based analyses find ZERO FDR-significant segments -- "
              "fully consistent conclusions, not a separate narrative.")

    corr = merged[["raw_p_median", "raw_p_mean"]].dropna().corr(method="spearman").iloc[0, 1]
    print(f"Spearman correlation of raw p-values (median-based vs mean-based) across all "
          f"segment/scheme tests: {corr:.3f}")

    agreement_out = TABLES / "subfamily_median_vs_mean_agreement.csv"
    merged.to_csv(agreement_out, index=False)
    print(f"Wrote {agreement_out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
