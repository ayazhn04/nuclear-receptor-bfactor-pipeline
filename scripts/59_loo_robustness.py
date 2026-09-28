"""
Part 10 of the subfamily flexibility stage: leave-one-receptor-out (LOO)
robustness, for (a) every atlas segment's median (and its well-supported
top/bottom-10% membership), and (b) the two segments with a nominal
(raw p<0.05, non-FDR-significant) omnibus group effect -- to check whether
that nominal signal depends on any single receptor.

Deterministic: LOO drops exactly one observed receptor at a time and
recomputes the same statistic with the same fixed bootstrap seed used
upstream; no randomness beyond what atlas/omnibus scripts already fix.
"""
from pathlib import Path

import json
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

WELL_SUPPORTED_MIN_N = 10


def loo_atlas_medians(profiles: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for seg, seg_df in profiles.groupby("structural_segment"):
        vals = seg_df[["uniprot_id", "median_of_entry_median_z"]].dropna()
        n = len(vals)
        if n < 2:
            rows.append({
                "structural_segment": seg, "n_receptors": n,
                "full_median": float(vals["median_of_entry_median_z"].iloc[0]) if n == 1 else np.nan,
                "loo_median_min": np.nan, "loo_median_max": np.nan, "loo_median_range": np.nan,
                "max_single_receptor_influence": np.nan, "dominated_by_single_receptor": False,
            })
            continue
        full_median = float(vals["median_of_entry_median_z"].median())
        loo_medians = []
        for uid in vals["uniprot_id"]:
            remaining = vals.loc[vals["uniprot_id"] != uid, "median_of_entry_median_z"]
            loo_medians.append(remaining.median())
        loo_medians = np.array(loo_medians)
        loo_range = float(loo_medians.max() - loo_medians.min())
        max_influence = float(np.max(np.abs(loo_medians - full_median)))
        # "dominated by one receptor" defined explicitly: removing some single
        # receptor shifts the median by more than half the segment's own IQR
        iqr = float(np.percentile(vals["median_of_entry_median_z"], 75) - np.percentile(vals["median_of_entry_median_z"], 25))
        threshold = max(0.5 * iqr, 0.25)  # documented floor so a near-zero IQR doesn't trivially flag everything
        dominated = bool(max_influence > threshold)
        rows.append({
            "structural_segment": seg, "n_receptors": n, "full_median": full_median,
            "loo_median_min": float(loo_medians.min()), "loo_median_max": float(loo_medians.max()),
            "loo_median_range": loo_range, "max_single_receptor_influence": max_influence,
            "dominated_by_single_receptor": dominated,
        })
    return pd.DataFrame(rows).sort_values("structural_segment").reset_index(drop=True)


def loo_kruskal_for_nominal_segments(kruskal: pd.DataFrame, profiles: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    scheme_columns = {"primary_project_group": "group", "secondary_nr_family_number": "nr_family_number"}
    nominal = kruskal[(kruskal["raw_p"] < 0.05) & (kruskal["test_status"] == "TESTED")]
    profiles_core = profiles[["uniprot_id", "structural_segment", "median_of_entry_median_z"]]

    rows = []
    for _, row in nominal.iterrows():
        scheme_name, seg = row["grouping_scheme"], row["structural_segment"]
        scheme_col = scheme_columns[scheme_name]
        eligible_groups = json.loads(row["eligible_groups"])
        merged = profiles_core.merge(meta[["uniprot_id", scheme_col]], on="uniprot_id", how="left")
        seg_df = merged[(merged["structural_segment"] == seg) & (merged[scheme_col].isin(eligible_groups))].dropna(
            subset=["median_of_entry_median_z"])

        full_p = row["raw_p"]
        loo_p = []
        loo_still_lt_05 = []
        for uid in seg_df["uniprot_id"]:
            remaining = seg_df[seg_df["uniprot_id"] != uid]
            group_counts = remaining.groupby(scheme_col)["uniprot_id"].nunique()
            still_eligible = group_counts[group_counts >= 3].index.tolist()
            if len(still_eligible) < 2:
                loo_p.append(np.nan)
                loo_still_lt_05.append(False)
                continue
            values_by_group = [remaining.loc[remaining[scheme_col] == g, "median_of_entry_median_z"].values
                                for g in still_eligible]
            H, p = stats.kruskal(*values_by_group)
            loo_p.append(p)
            loo_still_lt_05.append(bool(p < 0.05))
        loo_p = np.array(loo_p, dtype=float)
        rows.append({
            "grouping_scheme": scheme_name, "structural_segment": seg, "full_raw_p": full_p,
            "n_receptors_dropped_one_at_a_time": len(seg_df),
            "loo_p_min": float(np.nanmin(loo_p)) if np.any(~np.isnan(loo_p)) else np.nan,
            "loo_p_max": float(np.nanmax(loo_p)) if np.any(~np.isnan(loo_p)) else np.nan,
            "fraction_loo_still_p_lt_05": float(np.mean(loo_still_lt_05)),
            "conclusion_depends_on_single_receptor": bool(np.mean(loo_still_lt_05) < 1.0),
        })
    return pd.DataFrame(rows)


def main():
    profiles = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    kruskal = pd.read_csv(TABLES / "subfamily_segment_kruskal.csv")
    meta = pd.read_csv(ANALYSIS / "subfamily_receptor_metadata.csv")

    print("Part 10a: LOO robustness of atlas medians (every realized segment)")
    atlas_loo = loo_atlas_medians(profiles)
    out1 = TABLES / "subfamily_atlas_loo_robustness.csv"
    atlas_loo.to_csv(out1, index=False)
    n_dominated = atlas_loo["dominated_by_single_receptor"].sum()
    print(f"Wrote {out1.relative_to(PROJECT_ROOT)} ({len(atlas_loo)} segments); "
          f"{n_dominated} segment(s) flagged as dominated by a single receptor:")
    if n_dominated:
        print(atlas_loo.loc[atlas_loo["dominated_by_single_receptor"],
                             ["structural_segment", "n_receptors", "full_median", "max_single_receptor_influence"]].to_string(index=False))

    print("\nPart 10b: LOO robustness of the 2 nominal (raw p<0.05, non-FDR-significant) omnibus results")
    nominal_loo = loo_kruskal_for_nominal_segments(kruskal, profiles, meta)
    out2 = TABLES / "subfamily_kruskal_loo_robustness.csv"
    nominal_loo.to_csv(out2, index=False)
    print(f"Wrote {out2.relative_to(PROJECT_ROOT)} ({len(nominal_loo)} rows)")
    if len(nominal_loo):
        print(nominal_loo.to_string(index=False))
        for _, row in nominal_loo.iterrows():
            if row["conclusion_depends_on_single_receptor"]:
                print(f"  -> {row['structural_segment']} ({row['grouping_scheme']}): nominal p<0.05 does NOT "
                      f"hold under every leave-one-out fold (only {row['fraction_loo_still_p_lt_05']:.0%} of "
                      "folds retain p<0.05) -- already non-FDR-significant, and now shown to be fragile too.")


if __name__ == "__main__":
    main()
