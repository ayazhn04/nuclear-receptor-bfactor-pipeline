"""
Part 7 of the subfamily flexibility stage: pairwise post-hoc comparisons,
run ONLY for structural segments whose omnibus Kruskal-Wallis test survived
BH-FDR correction (significant_fdr_05 == True) in
reports/tables/subfamily_segment_kruskal.csv. Per instruction, post-hoc
tests are never run or highlighted when the omnibus test does not reach
FDR significance -- this avoids inflating the appearance of evidence with
dozens of uncorrected pairwise tests.

If (as observed) zero segments survive FDR in either grouping scheme, this
script still runs and still writes the output file (with the correct
schema, zero rows), so downstream tooling and tests always find it, and
prints why no post-hoc tests were performed.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

FDR_ALPHA = 0.05

SCHEME_COLUMNS = {
    "primary_project_group": "group",
    "secondary_nr_family_number": "nr_family_number",
}

POSTHOC_COLUMNS = [
    "grouping_scheme", "structural_segment", "group_A", "group_B",
    "n_A", "n_B", "median_A", "median_B", "median_difference",
    "cliffs_delta", "raw_p", "fdr_q",
]


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


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta: (n_greater - n_less) / (n_a * n_b), computed exactly
    (O(n*m), fine at this sample size) rather than via a normal
    approximation. Ranges -1 (all of b > all of a) to +1 (all of a > all
    of b)."""
    n_greater = sum((ai > b).sum() for ai in a)
    n_less = sum((ai < b).sum() for ai in a)
    return (n_greater - n_less) / (len(a) * len(b))


def run_posthoc_for_segment(seg_df: pd.DataFrame, scheme_col: str, eligible_groups: list[str], scheme_name: str, seg: str) -> list[dict]:
    rows = []
    pair_pvals = []
    pair_meta = []
    for i in range(len(eligible_groups)):
        for j in range(i + 1, len(eligible_groups)):
            gA, gB = eligible_groups[i], eligible_groups[j]
            a = seg_df.loc[seg_df[scheme_col] == gA, "median_of_entry_median_z"].dropna().values
            b = seg_df.loc[seg_df[scheme_col] == gB, "median_of_entry_median_z"].dropna().values
            if len(a) < 2 or len(b) < 2:
                continue
            _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            delta = cliffs_delta(a, b)
            pair_pvals.append(p)
            pair_meta.append({
                "grouping_scheme": scheme_name, "structural_segment": seg,
                "group_A": gA, "group_B": gB, "n_A": len(a), "n_B": len(b),
                "median_A": float(np.median(a)), "median_B": float(np.median(b)),
                "median_difference": float(np.median(a) - np.median(b)),
                "cliffs_delta": delta, "raw_p": p,
            })
    if pair_pvals:
        q = benjamini_hochberg(np.array(pair_pvals))
        for m, qi in zip(pair_meta, q):
            m["fdr_q"] = float(qi)
            rows.append(m)
    return rows


def main():
    kruskal = pd.read_csv(TABLES / "subfamily_segment_kruskal.csv")
    profiles = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    meta = pd.read_csv(ANALYSIS / "subfamily_receptor_metadata.csv")
    profiles_core = profiles[["uniprot_id", "structural_segment", "median_of_entry_median_z"]]

    all_rows = []
    for scheme_name, scheme_col in SCHEME_COLUMNS.items():
        sig = kruskal[(kruskal["grouping_scheme"] == scheme_name) & (kruskal["significant_fdr_05"] == True)]  # noqa: E712
        if len(sig) == 0:
            print(f"{scheme_name}: 0 segments survived omnibus BH-FDR -- no post-hoc tests performed "
                  "(per instruction: never run/highlight post-hoc when the omnibus test is not FDR-significant).")
            continue
        merged = profiles_core.merge(meta[["uniprot_id", scheme_col]], on="uniprot_id", how="left")
        for _, row in sig.iterrows():
            seg = row["structural_segment"]
            eligible_groups = json.loads(row["eligible_groups"])
            seg_df = merged[merged["structural_segment"] == seg]
            rows = run_posthoc_for_segment(seg_df, scheme_col, eligible_groups, scheme_name, seg)
            all_rows.extend(rows)
            print(f"{scheme_name}/{seg}: ran {len(rows)} post-hoc pairwise comparisons "
                  f"(omnibus FDR q={row['fdr_q']:.4f})")

    posthoc = pd.DataFrame(all_rows, columns=POSTHOC_COLUMNS)
    out = TABLES / "subfamily_segment_posthoc.csv"
    posthoc.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(posthoc)} rows)")
    if len(posthoc) == 0:
        print("0 post-hoc comparisons in total -- consistent with 0 segments surviving omnibus "
              "BH-FDR correction in either grouping scheme.")


if __name__ == "__main__":
    main()
