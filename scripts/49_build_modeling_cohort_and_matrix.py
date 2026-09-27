"""
Part 3-5 of the segment-level PCA/clustering modeling stage.

Part 3: define the PRIMARY modeling matrix. Starting candidate is the
objectively-derived >=50%-receptor-support feature set from
reports/tables/modeling_matrix_candidates.csv (candidate "C") -- NOT
hard-coded here. The result is reported explicitly, not silently accepted:
if the true >=50% set differs from the 7-core-segment expectation (it does
-- one loop column, loop_H3_H5, also clears 50% support), that is reported
rather than forced. A receptor enters the PRIMARY analysis cohort only if
it has observed (non-imputed) values for >=60% of the selected features;
this is an analysis-specific cohort rule recorded in
data/analysis/modeling_receptor_cohort.csv and does NOT touch
structures.csv or any master/normalization file.

Part 4: missing-value handling. PRIMARY = per-feature median imputation
using only the median across INCLUDED cohort receptors (never 0, never
supervised information). Sensitivity B = KNN imputation. Sensitivity C =
complete-case analysis (receptors with zero missing selected features),
run only if the resulting n is adequate for PCA/clustering.

Part 5: feature scaling. Per-instance Z-normalization (raw B -> z, already
done far upstream) and this cross-receptor feature standardization are
different operations -- after imputation, each segment feature is
standardized (mean 0, SD 1) across the primary cohort so segments with
larger cross-receptor variance don't dominate PCA/clustering purely from
scale. The biological segment Z summaries themselves are never overwritten
-- scaling is only ever applied to a copy used for modeling.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

PRIMARY_SUPPORT_THRESHOLD = 0.50  # the ">=50%" candidate set, chosen as the
                                   # initial primary candidate per Part 3
COHORT_OBSERVED_FRACTION_RULE = 0.60  # a receptor needs >=60% of selected
                                       # primary features observed to enter
                                       # the primary clustering cohort
KNN_SENSITIVITY_MAX_K = 5
COMPLETE_CASE_MIN_N = 10  # minimum complete-case n considered "adequate"


def select_primary_features() -> list[str]:
    candidates = pd.read_csv(TABLES / "modeling_matrix_candidates.csv")
    row = candidates[np.isclose(candidates["support_threshold"], PRIMARY_SUPPORT_THRESHOLD)]
    assert len(row) == 1, "expected exactly one >=50% candidate row"
    features = sorted(json.loads(row.iloc[0]["feature_names_json"]))
    print(f"Part 3: primary candidate features (>= {PRIMARY_SUPPORT_THRESHOLD:.0%} of 42 "
          f"coordinate receptors), objectively derived, {len(features)} features:")
    print(f"  {features}")
    core_helix_beta = [f for f in features if not f.startswith("loop_")]
    loops = [f for f in features if f.startswith("loop_")]
    print(f"  ({len(core_helix_beta)} core helix/beta-sheet segments, {len(loops)} loop segments "
          f"-- reported explicitly rather than silently forced to a pre-conceived count)")
    return features


def build_cohort(primary_features: list[str]) -> pd.DataFrame:
    rfs = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")
    coord = rfs[rfs["n_coordinate_candidates"] > 0][
        ["uniprot_id", "nr_code", "group"]
    ].sort_values("uniprot_id").reset_index(drop=True)

    median = pd.read_parquet(ANALYSIS / "receptor_segment_feature_matrix_median.parquet")
    median_full = median.set_index("uniprot_id").reindex(coord["uniprot_id"])

    n_selected = len(primary_features)
    observed = median_full[primary_features].notna()
    n_observed = observed.sum(axis=1).reindex(coord["uniprot_id"]).values
    observed_fraction = n_observed / n_selected

    cohort = coord.copy()
    cohort["n_selected_features"] = n_selected
    cohort["n_observed_features"] = n_observed
    cohort["observed_fraction"] = observed_fraction
    cohort["included_primary"] = observed_fraction >= COHORT_OBSERVED_FRACTION_RULE

    def reason(row):
        if row["included_primary"]:
            return ""
        if row["n_observed_features"] == 0:
            return (
                "No observed value for any of the primary selected segments "
                "(e.g. no usable structural alignment to any anchor template) "
                "-- analysis-cohort exclusion only, not removed from structures.csv."
            )
        return (
            f"Observed only {int(row['n_observed_features'])}/{n_selected} "
            f"({row['observed_fraction']:.0%}) of the primary selected segments, "
            f"below the {COHORT_OBSERVED_FRACTION_RULE:.0%} inclusion rule -- "
            "analysis-cohort exclusion only, not removed from structures.csv."
        )

    cohort["exclusion_reason"] = cohort.apply(reason, axis=1)
    cohort = cohort.sort_values("uniprot_id").reset_index(drop=True)

    out = ANALYSIS / "modeling_receptor_cohort.csv"
    cohort.to_csv(out, index=False)
    n_incl = int(cohort["included_primary"].sum())
    print(f"\nPart 3: wrote {out.relative_to(PROJECT_ROOT)}")
    print(f"  {n_incl}/{len(cohort)} coordinate receptors included in the primary cohort "
          f"(>= {COHORT_OBSERVED_FRACTION_RULE:.0%} of {n_selected} selected features observed)")
    excluded = cohort[~cohort["included_primary"]]
    if len(excluded):
        print(f"  Excluded: {excluded['uniprot_id'].tolist()}")
    return cohort


def median_impute(matrix: pd.DataFrame, included_ids: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-feature median imputation using only the median across included cohort
    receptors' OWN observed values. Never 0. Never uses group/cluster/ligand info."""
    included = matrix.loc[included_ids]
    feature_medians = included.median(axis=0, skipna=True)
    mask = included.isna()
    imputed = included.fillna(feature_medians)
    return imputed, mask


def knn_impute(matrix: pd.DataFrame, included_ids: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    from sklearn.impute import KNNImputer

    included = matrix.loc[included_ids]
    mask = included.isna()
    k = max(1, min(KNN_SENSITIVITY_MAX_K, len(included_ids) - 1))
    imputer = KNNImputer(n_neighbors=k, weights="distance")
    values = imputer.fit_transform(included.values)
    imputed = pd.DataFrame(values, index=included.index, columns=included.columns)
    return imputed, mask


def complete_case(matrix: pd.DataFrame, included_ids: list[str]) -> pd.DataFrame:
    included = matrix.loc[included_ids]
    complete = included.dropna(axis=0, how="any")
    return complete


def scale_features(matrix: pd.DataFrame) -> pd.DataFrame:
    """Cross-receptor feature standardization (mean 0, SD 1 per column), distinct
    from the upstream per-instance B-factor Z-normalization. Never overwrites the
    biological segment Z summaries -- operates on a modeling-only copy."""
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    values = scaler.fit_transform(matrix.values)
    return pd.DataFrame(values, index=matrix.index, columns=matrix.columns)


def main():
    primary_features = select_primary_features()
    cohort = build_cohort(primary_features)
    included_ids = cohort.loc[cohort["included_primary"], "uniprot_id"].tolist()

    median_mat = pd.read_parquet(ANALYSIS / "receptor_segment_feature_matrix_median.parquet")
    mean_mat = pd.read_parquet(ANALYSIS / "receptor_segment_feature_matrix_mean.parquet")

    median_full = median_mat.set_index("uniprot_id").reindex(cohort["uniprot_id"])[primary_features]
    mean_full = mean_mat.set_index("uniprot_id").reindex(cohort["uniprot_id"])[primary_features]

    # --- missingness before imputation, explicitly reported (Part 3) ---
    included_sub = median_full.loc[included_ids]
    obs = included_sub.notna()
    total_missing_before = 1.0 - obs.values.sum() / obs.size
    print(f"\nPart 3: PRIMARY cohort matrix before imputation: "
          f"{len(included_ids)} receptors x {len(primary_features)} features, "
          f"missingness = {total_missing_before:.1%}")

    # --- Part 4: PRIMARY median imputation (median matrix, included cohort) ---
    imputed_median, mask_median = median_impute(median_full, included_ids)
    unscaled_out = ANALYSIS / "receptor_segment_matrix_primary_unscaled.parquet"
    imputed_median.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(unscaled_out, index=False)
    mask_out = ANALYSIS / "receptor_segment_matrix_primary_imputation_mask.parquet"
    mask_median.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(mask_out, index=False)
    print(f"Part 4: wrote {unscaled_out.relative_to(PROJECT_ROOT)} "
          f"({imputed_median.shape[0]} x {imputed_median.shape[1]}, median-imputed) "
          f"and its imputation mask {mask_out.name}")

    also_out = ANALYSIS / "receptor_segment_matrix_primary.parquet"
    imputed_median.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(also_out, index=False)
    print(f"Part 4: wrote {also_out.relative_to(PROJECT_ROOT)} (alias of the primary median-imputed matrix)")

    # --- Part 5: scaling (after imputation) ---
    scaled_median = scale_features(imputed_median)
    scaled_out = ANALYSIS / "receptor_segment_matrix_primary_scaled.parquet"
    scaled_median.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(scaled_out, index=False)
    print(f"Part 5: wrote {scaled_out.relative_to(PROJECT_ROOT)} "
          f"(standardized mean=0/SD=1 per feature across the {len(included_ids)}-receptor primary cohort)")

    # --- Sensitivity A: mean-of-entry-median feature matrix, median-imputed ---
    imputed_mean, mask_mean = median_impute(mean_full, included_ids)
    mean_out = ANALYSIS / "receptor_segment_matrix_sensitivity_meanmatrix_imputed.parquet"
    imputed_mean.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(mean_out, index=False)
    scaled_mean = scale_features(imputed_mean)
    scaled_mean_out = ANALYSIS / "receptor_segment_matrix_sensitivity_meanmatrix_scaled.parquet"
    scaled_mean.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(scaled_mean_out, index=False)
    print(f"\nSensitivity (median-feature vs mean-feature matrix): wrote {mean_out.name} / {scaled_mean_out.name}")

    # --- Sensitivity B: KNN imputation ---
    imputed_knn, mask_knn = knn_impute(median_full, included_ids)
    knn_out = ANALYSIS / "receptor_segment_matrix_sensitivity_knn_imputed.parquet"
    imputed_knn.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(knn_out, index=False)
    knn_mask_out = ANALYSIS / "receptor_segment_matrix_sensitivity_knn_imputation_mask.parquet"
    mask_knn.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(knn_mask_out, index=False)
    scaled_knn = scale_features(imputed_knn)
    scaled_knn_out = ANALYSIS / "receptor_segment_matrix_sensitivity_knn_scaled.parquet"
    scaled_knn.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(scaled_knn_out, index=False)
    print(f"Sensitivity (median vs KNN imputation): wrote {knn_out.name} / {scaled_knn_out.name}")

    # --- Sensitivity C: complete-case analysis ---
    complete = complete_case(median_full, included_ids)
    n_complete = len(complete)
    adequate = n_complete >= COMPLETE_CASE_MIN_N
    complete_out = ANALYSIS / "receptor_segment_matrix_sensitivity_complete_case.parquet"
    complete.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(complete_out, index=False)
    status = "ADEQUATE" if adequate else "INADEQUATE_SAMPLE_SIZE"
    print(f"Sensitivity (complete-case): {n_complete} fully-observed receptors "
          f"(threshold for adequacy = {COMPLETE_CASE_MIN_N}) -> {status}; wrote {complete_out.name}")
    if adequate:
        scaled_complete = scale_features(complete)
        scaled_complete_out = ANALYSIS / "receptor_segment_matrix_sensitivity_complete_case_scaled.parquet"
        scaled_complete.reset_index().rename(columns={"index": "uniprot_id"}).to_parquet(scaled_complete_out, index=False)
        print(f"  wrote {scaled_complete_out.name}")

    summary = {
        "primary_features": primary_features,
        "n_primary_features": len(primary_features),
        "n_cohort_total": len(cohort),
        "n_cohort_included": len(included_ids),
        "cohort_observed_fraction_rule": COHORT_OBSERVED_FRACTION_RULE,
        "primary_support_threshold": PRIMARY_SUPPORT_THRESHOLD,
        "missingness_before_imputation": total_missing_before,
        "n_complete_case_receptors": n_complete,
        "complete_case_status": status,
        "knn_k": max(1, min(KNN_SENSITIVITY_MAX_K, len(included_ids) - 1)),
    }
    summary_out = ANALYSIS / "modeling_cohort_matrix_summary.json"
    summary_out.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nWrote {summary_out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
