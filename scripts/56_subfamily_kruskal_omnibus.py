"""
Part 5-6 of the subfamily flexibility stage: primary Kruskal-Wallis omnibus
test per structural segment (available-case, receptor-level, PRIMARY =
project functional `group`; SECONDARY = formal NR family-number), with
Benjamini-Hochberg FDR correction and a rank-based effect size
(epsilon-squared) with a bootstrap 95% CI.

Outcome value: median_of_entry_median_z from
receptor_segment_flexibility_profiles.parquet -- the approved
instance-Z -> PDB-entry-median -> segment-summary -> receptor-median
hierarchy. No imputation. A segment/group cell is included only if it has
>= 3 observed receptors (Part 4's eligibility table); a segment is tested
only if it has >= 2 eligible groups.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

MIN_N_PER_GROUP = 3
MIN_ELIGIBLE_GROUPS = 2
BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 0
FDR_ALPHA = 0.05

SCHEME_COLUMNS = {
    "primary_project_group": "group",
    "secondary_nr_family_number": "nr_family_number",
}


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Standard BH step-up FDR correction. Returns q-values aligned to input order."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p-value down
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def epsilon_squared(H: float, n: int, k: int) -> float:
    """Rank-based omnibus effect size for Kruskal-Wallis (epsilon-squared),
    documented formula (Tomczak & Tomczak 2014; King & Eckersley 2019):

        epsilon^2 = H / ((n^2 - 1) / (n + 1)) = H * (n + 1) / (n^2 - 1)

    where H is the Kruskal-Wallis statistic, n is the total sample size
    across all eligible groups in the test. Ranges ~0 (no effect) to 1
    (maximal rank separation). k (number of groups) is not used in this
    particular formula but is recorded alongside for documentation.
    """
    if n <= 1:
        return np.nan
    return H * (n + 1) / (n ** 2 - 1)


def bootstrap_epsilon_squared_ci(values_by_group: list[np.ndarray], seed=BOOTSTRAP_SEED, n_boot=BOOTSTRAP_N):
    """Resample receptors WITH replacement within each group (stratified
    bootstrap, preserving group sizes), recompute H and epsilon-squared each
    time; return the 2.5/97.5 percentile CI. Deterministic given the fixed
    seed."""
    rng = np.random.default_rng(seed)
    n_total = sum(len(v) for v in values_by_group)
    boot_eps = []
    for _ in range(n_boot):
        resampled = [rng.choice(v, size=len(v), replace=True) for v in values_by_group]
        try:
            H, _ = stats.kruskal(*resampled)
        except ValueError:
            continue
        boot_eps.append(epsilon_squared(H, n_total, len(values_by_group)))
    boot_eps = np.array(boot_eps)
    if len(boot_eps) == 0:
        return np.nan, np.nan
    return float(np.percentile(boot_eps, 2.5)), float(np.percentile(boot_eps, 97.5))


def run_omnibus_for_scheme(profiles: pd.DataFrame, meta: pd.DataFrame, scheme_col: str, scheme_name: str) -> pd.DataFrame:
    profiles_core = profiles[["uniprot_id", "structural_segment", "median_of_entry_median_z",
                               "mean_of_entry_median_z"]]
    merged = profiles_core.merge(meta[["uniprot_id", scheme_col]], on="uniprot_id", how="left")
    rows = []
    for seg, seg_df in merged.groupby("structural_segment"):
        group_counts = seg_df.groupby(scheme_col)["uniprot_id"].nunique()
        eligible_groups = group_counts[group_counts >= MIN_N_PER_GROUP].index.tolist()
        eligible_groups = sorted(eligible_groups)  # deterministic order
        if len(eligible_groups) < MIN_ELIGIBLE_GROUPS:
            rows.append({
                "grouping_scheme": scheme_name, "structural_segment": seg,
                "n_receptors": int(seg_df["uniprot_id"].nunique()), "n_groups": len(eligible_groups),
                "eligible_groups": json.dumps(eligible_groups),
                "H_statistic": np.nan, "raw_p": np.nan,
                "epsilon_squared": np.nan, "epsilon_squared_ci_low": np.nan, "epsilon_squared_ci_high": np.nan,
                "test_status": "NOT_TESTABLE_FEWER_THAN_2_ELIGIBLE_GROUPS",
            })
            continue

        values_by_group = []
        for g in eligible_groups:
            vals = seg_df.loc[seg_df[scheme_col] == g, "median_of_entry_median_z"].dropna().values
            values_by_group.append(vals)
        n_total = sum(len(v) for v in values_by_group)

        H, p = stats.kruskal(*values_by_group)
        eps = epsilon_squared(H, n_total, len(values_by_group))
        ci_low, ci_high = bootstrap_epsilon_squared_ci(values_by_group)

        rows.append({
            "grouping_scheme": scheme_name, "structural_segment": seg,
            "n_receptors": n_total, "n_groups": len(eligible_groups),
            "eligible_groups": json.dumps(eligible_groups),
            "H_statistic": H, "raw_p": p,
            "epsilon_squared": eps, "epsilon_squared_ci_low": ci_low, "epsilon_squared_ci_high": ci_high,
            "test_status": "TESTED",
        })
    return pd.DataFrame(rows)


def apply_fdr(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["fdr_q"] = np.nan
    df["significant_fdr_05"] = False
    tested = df["test_status"] == "TESTED"
    if tested.sum() > 0:
        q = benjamini_hochberg(df.loc[tested, "raw_p"].values)
        df.loc[tested, "fdr_q"] = q
        df.loc[tested, "significant_fdr_05"] = q <= FDR_ALPHA
    return df


def main():
    profiles = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    meta = pd.read_csv(ANALYSIS / "subfamily_receptor_metadata.csv")

    all_results = []
    for scheme_name, scheme_col in SCHEME_COLUMNS.items():
        print(f"\n--- Omnibus Kruskal-Wallis: {scheme_name} (column '{scheme_col}') ---")
        result = run_omnibus_for_scheme(profiles, meta, scheme_col, scheme_name)
        result = apply_fdr(result)
        all_results.append(result)
        n_tested = (result["test_status"] == "TESTED").sum()
        n_sig_raw = (result["raw_p"] < 0.05).sum()
        n_sig_fdr = result["significant_fdr_05"].sum()
        print(f"  {n_tested} segments tested; {n_sig_raw} with raw p<0.05; "
              f"{n_sig_fdr} significant after BH-FDR (q<=0.05)")
        if n_sig_raw > n_sig_fdr:
            print(f"  NOTE: {n_sig_raw - n_sig_fdr} segment(s) with raw p<0.05 do NOT survive "
                  "FDR correction and must not be described as significant.")

    combined = pd.concat(all_results, ignore_index=True)
    combined = combined.sort_values(["grouping_scheme", "raw_p"], na_position="last").reset_index(drop=True)

    out = TABLES / "subfamily_segment_kruskal.csv"
    combined.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(combined)} rows)")

    primary_sig = combined[(combined["grouping_scheme"] == "primary_project_group") & (combined["significant_fdr_05"])]
    print(f"\nPRIMARY grouping FDR-significant segments: {primary_sig['structural_segment'].tolist()}")


if __name__ == "__main__":
    main()
