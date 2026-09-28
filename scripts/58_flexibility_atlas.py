"""
Part 8-9 of the subfamily flexibility stage: superfamily structural-segment
flexibility atlas. For every realized structural segment (available-case,
receptor-level, all 41 receptors with at least one resolved segment --
never restricted to the 29-receptor primary clustering cohort), compute
descriptive statistics and a bootstrap 95% CI on the median.

Z=0 means "the average residue within that receptor's own instance(s)."
Positive segment Z = relatively flexible WITHIN that receptor's own
structures; negative = relatively rigid WITHIN that receptor's own
structures. This is never an absolute cross-protein mobility scale --
the atlas ranks segments by their typical WITHIN-RECEPTOR relative
flexibility, aggregated across receptors.

Part 9: two complementary descriptive rankings (by median receptor Z, and
top/bottom 10% where the segment count makes that meaningful) -- never a
biological flexible/rigid cutoff imposed on H8 or any unresolved segment.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 0


def bootstrap_median_ci(values: np.ndarray, seed=BOOTSTRAP_SEED, n_boot=BOOTSTRAP_N):
    if len(values) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot_medians = np.array([
        np.median(rng.choice(values, size=len(values), replace=True)) for _ in range(n_boot)
    ])
    return float(np.percentile(boot_medians, 2.5)), float(np.percentile(boot_medians, 97.5))


def build_atlas(profiles: pd.DataFrame, n_coord_receptors: int, value_col: str) -> pd.DataFrame:
    rows = []
    for seg, seg_df in profiles.groupby("structural_segment"):
        vals = seg_df[value_col].dropna().values
        n = len(vals)
        if n == 0:
            continue
        ci_low, ci_high = bootstrap_median_ci(vals)
        rows.append({
            "structural_segment": seg,
            "structural_class": seg_df["structural_class"].mode().iloc[0],
            "n_receptors": n,
            "median_receptor_z": float(np.median(vals)),
            "mean_receptor_z": float(np.mean(vals)),
            "q25": float(np.percentile(vals, 25)),
            "q75": float(np.percentile(vals, 75)),
            "iqr": float(np.percentile(vals, 75) - np.percentile(vals, 25)),
            "bootstrap_95ci_median_low": ci_low,
            "bootstrap_95ci_median_high": ci_high,
            "fraction_receptors_above_0": float((vals > 0).mean()),
            "fraction_receptors_below_0": float((vals < 0).mean()),
            "mapping_support_fraction": n / n_coord_receptors,
        })
    atlas = pd.DataFrame(rows).sort_values("median_receptor_z", ascending=False).reset_index(drop=True)
    return atlas


WELL_SUPPORTED_MIN_N = 10  # a documented, explicit threshold (not a silent filter) below
                            # which a single-digit n makes a segment's rank/CI unreliable
                            # (bootstrap CI is literally undefined at n=1)


def add_rankings(atlas: pd.DataFrame) -> pd.DataFrame:
    atlas = atlas.copy()
    atlas["rank_by_median_z_most_flexible_first"] = atlas["median_receptor_z"].rank(ascending=False, method="min").astype(int)
    atlas["well_supported"] = atlas["n_receptors"] >= WELL_SUPPORTED_MIN_N

    def _top_bottom_10pct(df: pd.DataFrame) -> tuple[set, set]:
        n = len(df)
        if n < 10:
            return set(), set()
        top_n = max(1, round(n * 0.10))
        df_sorted = df.sort_values("median_receptor_z", ascending=False)
        return set(df_sorted.head(top_n)["structural_segment"]), set(df_sorted.tail(top_n)["structural_segment"])

    # unfiltered ranking (every realized segment, however thin its support)
    top_all, bottom_all = _top_bottom_10pct(atlas)
    atlas["top_10pct_most_flexible"] = atlas["structural_segment"].isin(top_all)
    atlas["bottom_10pct_most_rigid"] = atlas["structural_segment"].isin(bottom_all)

    # well-supported-only ranking (n_receptors >= WELL_SUPPORTED_MIN_N) -- the more
    # defensible headline ranking, since a segment observed in only 1-2 receptors has
    # an undefined or near-meaningless bootstrap CI
    well = atlas[atlas["well_supported"]]
    top_well, bottom_well = _top_bottom_10pct(well)
    atlas["top_10pct_most_flexible_well_supported"] = atlas["structural_segment"].isin(top_well)
    atlas["bottom_10pct_most_rigid_well_supported"] = atlas["structural_segment"].isin(bottom_well)

    n_total, n_well = len(atlas), len(well)
    if n_total < 10:
        print(f"NOTE: only {n_total} segments in the atlas -- unfiltered top/bottom 10% would not "
              "be meaningful; flags left False for every row.")
    if n_well < 10:
        print(f"NOTE: only {n_well} segments have n_receptors >= {WELL_SUPPORTED_MIN_N} -- "
              "well-supported top/bottom 10% would not be meaningful; flags left False for every row.")
    return atlas


def main():
    profiles = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    rfs = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")
    n_coord = int((rfs["n_coordinate_candidates"] > 0).sum())

    assert "H8" not in profiles["structural_segment"].unique(), (
        "H8 was never confidently resolved -- it must never receive an atlas value"
    )

    print(f"Building superfamily flexibility atlas across {n_coord} coordinate receptors "
          f"(available-case per segment, no restriction to the primary clustering cohort)")

    atlas = build_atlas(profiles, n_coord, "median_of_entry_median_z")
    atlas = add_rankings(atlas)

    out = ANALYSIS / "superfamily_segment_flexibility_atlas.csv"
    atlas.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(atlas)} segments)")

    print("\nMost flexible (highest median receptor Z):")
    print(atlas.head(5)[["structural_segment", "n_receptors", "median_receptor_z",
                          "bootstrap_95ci_median_low", "bootstrap_95ci_median_high"]].to_string(index=False))
    print("\nMost rigid (lowest median receptor Z):")
    print(atlas.tail(5)[["structural_segment", "n_receptors", "median_receptor_z",
                          "bootstrap_95ci_median_low", "bootstrap_95ci_median_high"]].to_string(index=False))

    print(f"\nTop 10% most flexible (unfiltered, includes thin-n segments): "
          f"{atlas.loc[atlas['top_10pct_most_flexible'], 'structural_segment'].tolist()}")
    print(f"Bottom 10% most rigid (unfiltered): "
          f"{atlas.loc[atlas['bottom_10pct_most_rigid'], 'structural_segment'].tolist()}")
    print(f"\nTop 10% most flexible (well-supported, n>={WELL_SUPPORTED_MIN_N} receptors): "
          f"{atlas.loc[atlas['top_10pct_most_flexible_well_supported'], 'structural_segment'].tolist()}")
    print(f"Bottom 10% most rigid (well-supported, n>={WELL_SUPPORTED_MIN_N} receptors): "
          f"{atlas.loc[atlas['bottom_10pct_most_rigid_well_supported'], 'structural_segment'].tolist()}")

    # secondary: mean-feature version (Part 13 groundwork, reused later for the
    # sensitivity comparison) -- mean_of_entry_mean_z is the exact aggregation
    # underlying receptor_segment_feature_matrix_mean.parquet (entry-level MEAN,
    # then MEAN across entries), the true "MEAN summary" sibling of the approved
    # median pipeline -- not mean_of_entry_median_z, which still uses entry-level
    # medians.
    atlas_mean = build_atlas(profiles, n_coord, "mean_of_entry_mean_z")
    atlas_mean = add_rankings(atlas_mean)
    out_mean = ANALYSIS / "superfamily_segment_flexibility_atlas_mean_sensitivity.csv"
    atlas_mean.to_csv(out_mean, index=False)
    print(f"\nWrote {out_mean.relative_to(PROJECT_ROOT)} (mean-of-entry-median sensitivity view)")


if __name__ == "__main__":
    main()
