"""
Part 6 of the functional-site specificity/validation stage: compare each
original whole-LBD-background effect (from the prior stage,
current_functional_site_within_receptor_tests.csv) against its
structurally matched-background counterpart (functional_site_matched_
tests.csv, this stage) for every receptor x site that has both.
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TABLES = PROJECT_ROOT / "reports" / "tables"


def main():
    original = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv")
    matched = pd.read_csv(TABLES / "functional_site_matched_tests.csv")

    merged = original.merge(
        matched, on=["uniprot_id", "nr_code", "functional_site"], suffixes=("_original", "_matched")
    )
    merged["effect_direction_preserved"] = np.sign(merged["median_paired_delta_original"]) == np.sign(merged["median_paired_delta_matched"])
    merged["significance_preserved"] = merged["significant_within_receptor_original"] & merged["significant_within_receptor_matched"]
    merged["significance_lost"] = merged["significant_within_receptor_original"] & ~merged["significant_within_receptor_matched"]
    merged["significance_gained"] = ~merged["significant_within_receptor_original"] & merged["significant_within_receptor_matched"]

    out_cols = [
        "uniprot_id", "nr_code", "functional_site", "matched_tier_used",
        "median_paired_delta_original", "median_paired_delta_matched",
        "q_within_receptor_original", "q_within_receptor_matched",
        "significant_within_receptor_original", "significant_within_receptor_matched",
        "effect_direction_preserved", "significance_preserved", "significance_lost", "significance_gained",
    ]
    sensitivity = merged.rename(columns={
        "median_paired_delta_original": "original_delta", "median_paired_delta_matched": "matched_delta",
        "q_within_receptor_original": "original_q", "q_within_receptor_matched": "matched_q",
    })
    sensitivity = sensitivity.rename(columns={
        "significant_within_receptor_original": "significant_original",
        "significant_within_receptor_matched": "significant_matched",
    })
    out_cols = [c.replace("median_paired_delta_original", "original_delta")
                 .replace("median_paired_delta_matched", "matched_delta")
                 .replace("q_within_receptor_original", "original_q")
                 .replace("q_within_receptor_matched", "matched_q")
                 .replace("significant_within_receptor_original", "significant_original")
                 .replace("significant_within_receptor_matched", "significant_matched")
                for c in out_cols]

    out = TABLES / "functional_site_background_sensitivity.csv"
    sensitivity[out_cols].to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(sensitivity)} rows, "
          f"{len(original)} original tests -> {len(sensitivity)} have a matched-background counterpart)")

    n_dir_preserved = sensitivity["effect_direction_preserved"].sum()
    n_sig_preserved = sensitivity["significance_preserved"].sum()
    n_sig_lost = sensitivity["significance_lost"].sum()
    n_sig_gained = sensitivity["significance_gained"].sum()
    print(f"\nDirection preserved: {n_dir_preserved}/{len(sensitivity)}")
    print(f"Significance preserved (sig in both): {n_sig_preserved}")
    print(f"Significance LOST (sig original, not sig matched): {n_sig_lost}")
    print(f"Significance GAINED (not sig original, sig matched): {n_sig_gained}")

    if n_sig_lost:
        print("\nResults where matching structural context explains away the original effect:")
        print(sensitivity.loc[sensitivity["significance_lost"],
              ["nr_code", "functional_site", "matched_tier_used", "original_delta", "matched_delta",
               "original_q", "matched_q"]].to_string(index=False))


if __name__ == "__main__":
    main()
