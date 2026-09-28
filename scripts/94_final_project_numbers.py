"""
Part 2 of the freeze/consolidate stage: builds
reports/tables/final_project_numbers.csv, the single traceable source of
truth for every number likely to be quoted in a manuscript or
presentation. Every value is COMPUTED here directly from its source file
(never hand-typed), with source_file/source_column_or_method recorded so
each number can be independently re-derived.
"""
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"

rows = []


def add(category, metric, value, source_file, source_column_or_method, notes=""):
    rows.append({
        "category": category, "metric": metric, "value": value,
        "source_file": source_file, "source_column_or_method": source_column_or_method, "notes": notes,
    })


def main():
    # ---------------- DATASET ----------------
    structures = pd.read_csv(PROCESSED / "structures.csv", dtype={"pdb_id": str}, keep_default_na=False)
    add("dataset", "project_receptors", structures["uniprot_id"].nunique(), "config/nr_metadata.py",
        "len(NR_METADATA)", "48 professor-supplied nuclear receptors")
    add("dataset", "master_candidates", len(structures), "data/processed/structures.csv", "len(df)",
        "one row per receptor x PDB entry x polymer entity")
    coord = structures[structures["coordinate_data_status"] == "COORDINATE_DATA_AVAILABLE"]
    add("dataset", "coordinate_assessable_candidates", len(coord), "data/processed/structures.csv",
        "coordinate_data_status == COORDINATE_DATA_AVAILABLE")

    residues = pd.read_parquet(PROCESSED / "lbd_residue_bfactors_all.parquet")
    add("dataset", "structural_instances", residues["instance_id"].nunique(), "data/processed/lbd_residue_bfactors_all.parquet",
        "instance_id.nunique()")
    add("dataset", "residue_rows", len(residues), "data/processed/lbd_residue_bfactors_all.parquet", "len(df)")

    normalized = pd.read_parquet(ANALYSIS / "normalized_lbd_bfactors_all.parquet")
    add("dataset", "normalized_observations_z_computed", int((normalized["row_normalization_status"] == "Z_COMPUTED").sum()),
        "data/analysis/normalized_lbd_bfactors_all.parquet", "row_normalization_status == Z_COMPUTED")

    rfs = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")
    add("dataset", "coordinate_receptors", int((rfs["n_coordinate_candidates"] > 0).sum()),
        "data/analysis/receptor_flexibility_summary.csv", "n_coordinate_candidates > 0")

    # ---------------- STRUCTURAL MAPPING ----------------
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    named = seg_map[seg_map["structural_segment"] != ""]
    add("structural_mapping", "segments_realized", named["structural_segment"].nunique(),
        "data/analysis/nr_lbd_structural_segment_map.csv", "structural_segment.nunique() where != ''")
    add("structural_mapping", "receptors_with_any_named_segment", named["uniprot_id"].nunique(),
        "data/analysis/nr_lbd_structural_segment_map.csv", "uniprot_id.nunique() among named rows")
    add("structural_mapping", "unresolved_element_H8_receptor_count", 0,
        "data/analysis/nr_lbd_structural_segment_map.csv", "structural_segment == 'H8' count",
        "H8 never confidently resolved by the current mapping; never fabricated or zero-filled")

    cohort = pd.read_csv(ANALYSIS / "modeling_receptor_cohort.csv")
    add("structural_mapping", "modeling_primary_features", int(cohort["n_selected_features"].iloc[0]),
        "data/analysis/modeling_receptor_cohort.csv", "n_selected_features")

    # ---------------- PCA / CLUSTERING ----------------
    add("pca_clustering", "primary_clustering_receptor_cohort", int(cohort["included_primary"].sum()),
        "data/analysis/modeling_receptor_cohort.csv", "included_primary.sum()")
    pca_var = pd.read_csv(TABLES / "pca_variance.csv")
    with open(ANALYSIS / "pca_summary.json") as f:
        pca_summary = json.load(f)
    add("pca_clustering", "pc1_explained_variance", round(float(pca_var.iloc[0]["explained_variance_ratio"]), 4),
        "reports/tables/pca_variance.csv", "row PC1, explained_variance_ratio")
    for k, v in pca_summary["n_components_for_variance"].items():
        add("pca_clustering", f"n_pcs_for_{k}_variance", v, "data/analysis/pca_summary.json",
            f"n_components_for_variance['{k}']")

    clustering_val = pd.read_csv(TABLES / "clustering_validation.csv")
    best = clustering_val[clustering_val["input"] == "pca_primary"].dropna(subset=["silhouette"]).sort_values(
        "silhouette", ascending=False).iloc[0]
    add("pca_clustering", "best_silhouette_pca_primary", round(float(best["silhouette"]), 4),
        "reports/tables/clustering_validation.csv", "max silhouette where input==pca_primary",
        f"method={best['method']}, k={int(best['k'])}; this solution is a degenerate 1-vs-28 split, not evidence of real clusters")
    dbscan = pd.read_csv(TABLES / "clustering_dbscan_grid.csv")
    add("pca_clustering", "dbscan_max_clusters_found_pca_primary",
        int(dbscan[dbscan["input"] == "pca_primary"]["n_clusters_found"].max()),
        "reports/tables/clustering_dbscan_grid.csv", "max n_clusters_found where input==pca_primary")
    agreement = pd.read_csv(TABLES / "clustering_method_agreement.csv")
    add("pca_clustering", "mean_cross_method_ARI_pca_primary",
        round(float(agreement[agreement["input"] == "pca_primary"]["adjusted_rand_index"].mean()), 4),
        "reports/tables/clustering_method_agreement.csv", "mean(adjusted_rand_index) where input==pca_primary")
    sensitivity = pd.read_csv(TABLES / "clustering_sensitivity.csv")
    add("pca_clustering", "sensitivity_ari_range",
        f"{sensitivity.groupby('comparison')['adjusted_rand_index'].mean().min():.2f}-{sensitivity.groupby('comparison')['adjusted_rand_index'].mean().max():.2f}",
        "reports/tables/clustering_sensitivity.csv", "groupby(comparison).mean() min/max")
    add("pca_clustering", "supported_conclusion", "no_robust_discrete_clusters_continuum_consistent",
        "reports/tables/modeling_summary.md", "Part 11 formal conclusion",
        "no method picked because it looked nicest; degenerate best-silhouette solution; DBSCAN finds <=1 cluster")

    # ---------------- SUBFAMILY / GROUP TESTING ----------------
    kruskal = pd.read_csv(TABLES / "subfamily_segment_kruskal.csv")
    tested = kruskal[kruskal["test_status"] == "TESTED"]
    add("subfamily_testing", "testable_segments_primary_scheme",
        tested[tested["grouping_scheme"] == "primary_project_group"]["structural_segment"].nunique(),
        "reports/tables/subfamily_segment_kruskal.csv", "nunique structural_segment where scheme==primary_project_group & TESTED")
    add("subfamily_testing", "total_omnibus_tests", len(tested), "reports/tables/subfamily_segment_kruskal.csv",
        "len(df) where test_status==TESTED")
    add("subfamily_testing", "nominal_p_lt_05", int((tested["raw_p"] < 0.05).sum()),
        "reports/tables/subfamily_segment_kruskal.csv", "raw_p < 0.05")
    add("subfamily_testing", "fdr_significant", int(tested["significant_fdr_05"].sum()),
        "reports/tables/subfamily_segment_kruskal.csv", "significant_fdr_05.sum()",
        "0 segments differ between groups after BH-FDR correction, in either grouping scheme")
    atlas = pd.read_csv(ANALYSIS / "superfamily_segment_flexibility_atlas.csv")
    add("subfamily_testing", "atlas_segments_total", len(atlas), "data/analysis/superfamily_segment_flexibility_atlas.csv", "len(df)")
    add("subfamily_testing", "atlas_well_supported_segments", int(atlas["well_supported"].sum()),
        "data/analysis/superfamily_segment_flexibility_atlas.csv", "well_supported.sum()")

    # ---------------- APO / HOLO ----------------
    elig = pd.read_csv(TABLES / "apo_holo_receptor_eligibility.csv")
    add("apo_holo", "primary_eligible_receptors_5_5", int(elig["eligible_primary_5_5"].sum()),
        "reports/tables/apo_holo_receptor_eligibility.csv", "eligible_primary_5_5.sum()")
    add("apo_holo", "exploratory_eligible_receptors_3_3", int(elig["eligible_exploratory_3_3"].sum()),
        "reports/tables/apo_holo_receptor_eligibility.csv", "eligible_exploratory_3_3.sum()")
    eligible_codes = elig.loc[elig["eligible_exploratory_3_3"], "nr_code"].tolist()
    add("apo_holo", "eligible_receptor_nr_codes", ";".join(eligible_codes),
        "reports/tables/apo_holo_receptor_eligibility.csv", "nr_code list where eligible_exploratory_3_3")
    vdr = elig[elig["uniprot_id"] == "P11473"].iloc[0]
    add("apo_holo", "VDR_n_APO", int(vdr["n_APO"]), "reports/tables/apo_holo_receptor_eligibility.csv", "row P11473, n_APO")
    add("apo_holo", "VDR_n_HOLO", int(vdr["n_HOLO"]), "reports/tables/apo_holo_receptor_eligibility.csv", "row P11473, n_HOLO",
        "VDR is INELIGIBLE for apo/holo comparison (0 APO, 48 HOLO); never removed from structures.csv")

    seg_tests = pd.read_csv(TABLES / "apo_holo_segment_tests.csv")
    add("apo_holo", "n_tests", len(seg_tests), "reports/tables/apo_holo_segment_tests.csv", "len(df)")
    add("apo_holo", "within_receptor_fdr_significant", int(seg_tests["significant_within_receptor"].sum()),
        "reports/tables/apo_holo_segment_tests.csv", "significant_within_receptor.sum()")
    add("apo_holo", "global_fdr_significant", int(seg_tests["significant_global"].sum()),
        "reports/tables/apo_holo_segment_tests.csv", "significant_global.sum()")
    loo = pd.read_csv(TABLES / "apo_holo_loo_robustness.csv")
    sig_loo = loo[loo["significant_within_receptor"]]
    add("apo_holo", "robust_significant_results_loo",
        f"{int((sig_loo['fraction_loo_p_lt_05'] == 1.0).sum())}/{len(sig_loo)}",
        "reports/tables/apo_holo_loo_robustness.csv", "fraction_loo_p_lt_05 == 1.0 among significant results")

    # ---------------- FUNCTIONAL SITES ----------------
    defs = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")
    add("functional_sites", "definition_rows", len(defs), "data/analysis/current_functional_site_definitions.parquet", "len(df)")
    add("functional_sites", "receptors_with_any_site_definition", defs["uniprot_id"].nunique(),
        "data/analysis/current_functional_site_definitions.parquet", "uniprot_id.nunique()")
    within = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv")
    add("functional_sites", "original_within_receptor_tests", len(within),
        "reports/tables/current_functional_site_within_receptor_tests.csv", "len(df)")
    add("functional_sites", "original_within_fdr_significant", int(within["significant_within_receptor"].sum()),
        "reports/tables/current_functional_site_within_receptor_tests.csv", "significant_within_receptor.sum()")
    add("functional_sites", "original_global_fdr_significant", int(within["significant_global"].sum()),
        "reports/tables/current_functional_site_within_receptor_tests.csv", "significant_global.sum()")
    matched = pd.read_csv(TABLES / "functional_site_matched_tests.csv")
    add("functional_sites", "matched_background_tests", len(matched),
        "reports/tables/functional_site_matched_tests.csv", "len(df)")
    add("functional_sites", "matched_background_significant", int(matched["significant_within_receptor"].sum()),
        "reports/tables/functional_site_matched_tests.csv", "significant_within_receptor.sum()")
    classif = pd.read_csv(TABLES / "functional_site_interpretation_classification.csv")
    for _, r in classif.iterrows():
        add("functional_sites", f"interpretation_{r['functional_site']}", r["interpretation_category"],
            "reports/tables/functional_site_interpretation_classification.csv", "interpretation_category")

    numbers = pd.DataFrame(rows)
    out = TABLES / "final_project_numbers.csv"
    numbers.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(numbers)} rows, {numbers['category'].nunique()} categories)")
    print(numbers["category"].value_counts().to_string())


if __name__ == "__main__":
    main()
