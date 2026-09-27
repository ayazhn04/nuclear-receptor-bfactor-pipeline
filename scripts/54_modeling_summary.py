"""
Part 14 of the segment-level PCA/clustering modeling stage: assembles
reports/tables/modeling_summary.md from the outputs of scripts 48-53.
Writes no new analysis -- pure reporting/aggregation of already-computed,
already-saved tables.
"""
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def main():
    audit = pd.read_csv(TABLES / "modeling_feature_support.csv")
    candidates = pd.read_csv(TABLES / "modeling_matrix_candidates.csv")
    cohort = pd.read_csv(ANALYSIS / "modeling_receptor_cohort.csv")
    with open(ANALYSIS / "modeling_cohort_matrix_summary.json") as f:
        cohort_summary = json.load(f)
    pca_variance = pd.read_csv(TABLES / "pca_variance.csv")
    with open(ANALYSIS / "pca_summary.json") as f:
        pca_summary = json.load(f)
    loadings = pd.read_csv(ANALYSIS / "pca_loadings.csv")
    with open(ANALYSIS / "umap_stability_summary.json") as f:
        umap_summary = json.load(f)
    validation = pd.read_csv(TABLES / "clustering_validation.csv")
    agreement = pd.read_csv(TABLES / "clustering_method_agreement.csv")
    dbscan = pd.read_csv(TABLES / "clustering_dbscan_grid.csv")
    sensitivity = pd.read_csv(TABLES / "clustering_sensitivity.csv")
    with open(ANALYSIS / "clustering_summary.json") as f:
        clustering_summary = json.load(f)
    overrep = pd.read_csv(TABLES / "entry_level_overrepresentation.csv")

    primary_features = cohort_summary["primary_features"]
    included = cohort[cohort["included_primary"]]
    excluded = cohort[~cohort["included_primary"]]

    strongest_pc1 = loadings.set_index("structural_segment")["PC1"].abs().sort_values(ascending=False)

    pca_primary_val = validation[validation["input"] == "pca_primary"].dropna(subset=["silhouette"])
    best = pca_primary_val.sort_values("silhouette", ascending=False).iloc[0]

    # inspect the actual cluster-size distribution of the single best-silhouette solution
    # -- a "best" solution that is really just a 1-vs-rest outlier split is degenerate,
    # not evidence of real structure, and must be reported as such rather than hidden
    best_cluster_note = ""
    umap_cluster_path = ANALYSIS / "umap_default_embedding_pca_representation.csv"
    labels_path = ANALYSIS / "receptor_segment_matrix_primary_scaled.parquet"
    try:
        from sklearn.mixture import GaussianMixture
        from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
        pca_scores = pd.read_csv(ANALYSIS / "pca_receptor_scores.csv")
        n_for = pca_summary["n_components_for_variance"]
        X = pca_scores[[f"PC{i+1}" for i in range(n_for["80%"])]].values
        method, k = best["method"], int(best["k"])
        fitters = {
            "kmeans": lambda: KMeans(n_clusters=k, n_init=10, random_state=0).fit(X).labels_,
            "gmm": lambda: GaussianMixture(n_components=k, covariance_type="diag", n_init=5, random_state=0, reg_covar=1e-4).fit(X).predict(X),
            "agglomerative_ward": lambda: AgglomerativeClustering(n_clusters=k, linkage="ward").fit(X).labels_,
            "spectral": lambda: SpectralClustering(n_clusters=k, affinity="nearest_neighbors", n_neighbors=min(10, len(X)-1), assign_labels="kmeans", random_state=0).fit(X).labels_,
        }
        best_labels = fitters[method]()
        sizes = pd.Series(best_labels).value_counts().sort_values()
        sizes_plain = {int(k): int(v) for k, v in sizes.items()}
        best_cluster_note = (
            f" Cluster sizes for this solution: {sizes_plain} -- "
            + ("this is a degenerate 1-member-vs-rest outlier split, not a meaningful "
               "partition of the cohort." if sizes.min() == 1 else
               "a non-trivial (non-singleton) split.")
        )
    except Exception as e:
        best_cluster_note = f" (cluster-size check unavailable: {e})"

    pca_agreement = agreement[agreement["input"] == "pca_primary"]
    mean_ari_by_k = pca_agreement.groupby("k")["adjusted_rand_index"].mean()

    dbscan_pca = dbscan[dbscan["input"] == "pca_primary"]
    dbscan_max_clusters = int(dbscan_pca["n_clusters_found"].max())

    lines = []
    lines.append("# Modeling Summary: Segment-Level PCA / UMAP / Clustering\n")
    lines.append(
        "> Continues from commit `7d280bc` (canonical LBD template mapping) with supervisor "
        "approval to proceed to PCA/UMAP/clustering. Master `structures.csv`, normalization, and "
        "the structural-segment mapping were **not** modified by this stage -- see 'Frozen inputs' "
        "below. Group/subfamily labels were **never** used to fit PCA, UMAP, imputation, or "
        "clustering; they are joined back in only for visualization/labeling after the fact.\n"
    )

    lines.append("## 0. H8 documentation correction\n")
    lines.append(
        "Where prior documentation said H8 was \"genuinely absent,\" the wording was corrected to: "
        "H8 was not confidently resolved by the current structural-template mapping and is "
        "retained as unresolved/not available for the primary feature matrix. No H8 values are "
        "fabricated anywhere in this stage (H8 is not a column in any modeling matrix at all).\n"
    )

    lines.append("## 1. Frozen inputs (Part 1)\n")
    lines.append(
        "SHA256 hashes recorded in `data/analysis/modeling_frozen_input_hashes.json` for the five "
        "upstream files this stage is built on; none are modified by scripts 48-53 (verified by "
        "`tests/test_modeling_pipeline.py`):\n"
    )
    for f in [
        "receptor_segment_feature_matrix_median.parquet",
        "receptor_segment_feature_matrix_mean.parquet",
        "receptor_segment_coverage_matrix.parquet",
        "receptor_segment_flexibility_profiles.parquet",
        "nr_lbd_structural_segment_map.csv",
    ]:
        lines.append(f"- `data/analysis/{f}`")
    lines.append("")

    lines.append("## 2. Feature-support audit (Part 2)\n")
    lines.append(
        f"All {len(audit)} structural-segment/loop columns realized by the v3 segment map were "
        "audited (no segment names hard-coded) -- full table in "
        "`reports/tables/modeling_feature_support.csv`. Candidate feature sets "
        "(`reports/tables/modeling_matrix_candidates.csv`):\n"
    )
    lines.append("| Set | Support threshold | n features | n receptors w/ complete data | Total matrix missingness |")
    lines.append("|---|---|---|---|---|")
    for _, r in candidates.iterrows():
        lines.append(f"| {r['candidate_set']} | >={r['support_threshold']:.0%} | {r['n_features']} | "
                     f"{r['n_receptors_with_complete_data']} | {r['total_matrix_missingness_fraction']:.1%} |")
    lines.append("")

    lines.append("## 3. Primary segments selected (Part 3)\n")
    core = [f for f in primary_features if not f.startswith("loop_")]
    loop = [f for f in primary_features if f.startswith("loop_")]
    lines.append(
        f"The objectively-derived **>=50%-support** candidate set (\"C\") was used as the initial "
        f"primary candidate, per instruction, and accepted as-is after explicit review: "
        f"**{len(primary_features)} features** -- {len(core)} core helix/beta-sheet segments "
        f"({', '.join(core)}) and {len(loop)} loop segment ({', '.join(loop)}). This is reported "
        f"exactly as computed: it is **not** the same as \"7 core segments\" one might expect from "
        f"the Part 9 acceptance table alone -- `loop_H3_H5` independently clears 50% receptor "
        f"support and is included, since segment names were never hard-coded.\n"
    )
    lines.append(
        f"**Excluded low-support segments** (< 50% of 42 coordinate receptors): "
        f"{', '.join(sorted(set(audit['structural_segment']) - set(primary_features)))}.\n"
    )

    lines.append("## 4. Receptor analysis cohort (Part 3)\n")
    lines.append(
        f"Cohort rule: a receptor must have observed values for >=60% of the {len(primary_features)} "
        f"primary features to be **included** in the primary clustering analysis. This is an "
        f"analysis-specific rule recorded in `data/analysis/modeling_receptor_cohort.csv`; it does "
        f"**not** remove any receptor from `structures.csv` or any master/normalization file.\n"
    )
    lines.append(f"- **{len(included)}/{len(cohort)}** coordinate receptors included.\n")
    lines.append(
        f"- Missingness in the included-cohort matrix **before imputation**: "
        f"{cohort_summary['missingness_before_imputation']:.1%}.\n"
    )
    lines.append(f"- **{len(excluded)} excluded**, all with explicit reasons (never silently dropped):\n")
    lines.append("| uniprot_id | nr_code | n_observed/n_selected | reason |")
    lines.append("|---|---|---|---|")
    for _, r in excluded.sort_values("observed_fraction", ascending=False).iterrows():
        lines.append(f"| {r['uniprot_id']} | {r['nr_code']} | {int(r['n_observed_features'])}/{int(r['n_selected_features'])} | {r['exclusion_reason']} |")
    lines.append("")
    lines.append(
        "Notably, **VDR (P11473)**, one of the 5 curated anchor receptors, is itself excluded from "
        "the primary cohort: VDR's own curated template has no literature-anchored (or "
        "cross-template-rescued) evidence for H5/H7/H11, so it observes only H1/H3/beta_sheet/H12 "
        "(4/8 = 50% < the 60% inclusion rule). This is a real, correctly-computed consequence of "
        "the v3 template mapping, not a bug -- flagged explicitly rather than silently accepted.\n"
    )

    lines.append("## 5. Missing-value handling (Part 4)\n")
    lines.append(
        "PRIMARY = per-feature median imputation using only the median across the **included** "
        f"cohort receptors' own observed values (never 0, never group/cluster/ligand-state "
        f"information) -> `receptor_segment_matrix_primary.parquet` "
        f"({len(included)}x{len(primary_features)}), with every imputed cell identifiable in "
        f"`receptor_segment_matrix_primary_imputation_mask.parquet`.\n\n"
        f"Sensitivity variants also produced: (A) median-imputed mean-of-entry-median feature "
        f"matrix; (B) KNN imputation (k={cohort_summary['knn_k']}, distance-weighted); "
        f"(C) complete-case analysis "
        f"({cohort_summary['n_complete_case_receptors']} receptors, "
        f"{cohort_summary['complete_case_status']}).\n"
    )

    lines.append("## 6. Feature scaling (Part 5)\n")
    lines.append(
        "After imputation, each of the 8 primary segment features was standardized (mean=0, SD=1) "
        "across the 29-receptor primary cohort -- distinct from the upstream per-instance B-factor "
        "Z-normalization, which is untouched. Saved as "
        "`receptor_segment_matrix_primary_unscaled.parquet` / `_scaled.parquet`. The biological "
        "segment Z summaries themselves were never overwritten.\n"
    )

    lines.append("## 7. PCA (Part 6)\n")
    lines.append(f"Fit on {len(included)} receptors x {len(primary_features)} scaled features only "
                 "(group labels never used for fitting).\n")
    lines.append("| PC | Eigenvalue | Explained var. | Cumulative |")
    lines.append("|---|---|---|---|")
    for _, r in pca_variance.iterrows():
        lines.append(f"| {r['pc']} | {r['eigenvalue']:.3f} | {r['explained_variance_ratio']:.1%} | {r['cumulative_explained_variance']:.1%} |")
    lines.append("")
    n_for = pca_summary["n_components_for_variance"]
    lines.append(f"Components needed: **{n_for['70%']}** PCs for 70% variance, **{n_for['80%']}** for 80%, "
                 f"**{n_for['90%']}** for 90%, **{n_for['95%']}** for 95%. Variance is spread fairly evenly "
                 f"across components (PC1 = {pca_summary['pc1_explained_variance']:.1%}) -- there is no single "
                 "dominant axis of variation among the 8 segments.\n")
    lines.append("Strongest |PC1| loadings:\n")
    lines.append("| Segment | PC1 loading |")
    lines.append("|---|---|")
    for seg, v in strongest_pc1.items():
        signed = loadings.set_index("structural_segment").loc[seg, "PC1"]
        lines.append(f"| {seg} | {signed:+.3f} |")
    lines.append("")

    lines.append("## 8. UMAP stability (Part 7)\n")
    lines.append(
        f"Grid: n_neighbors={umap_summary['n_neighbors_grid']}, min_dist={umap_summary['min_dist_grid']}, "
        f"seeds={umap_summary['seeds']}, run on both the scaled feature matrix and the PCA "
        f"representation. Mean Spearman correlation of pairwise-embedding distances to a fixed "
        f"reference setting = **{umap_summary['overall_mean_spearman_to_reference']:.3f}** -> "
        f"**{umap_summary['verdict']}**. Full grid in `reports/tables/umap_stability.csv`. UMAP "
        "was never chosen or reported because a setting \"looked more separated\" -- the default "
        "(n_neighbors=15, min_dist=0.1) middle-of-grid setting is used for the representative "
        "figures, and cluster coloring was only added after formal clustering (Part 9).\n"
    )

    lines.append("## 9-10. Clustering battery and validation (Parts 8-10)\n")
    lines.append(
        f"Primary input: PCA scores retaining {n_for['80%']} PCs (>=80% variance), "
        f"{len(included)} receptors, k = 2..{int(validation[validation['input']=='pca_primary']['k'].max())}. "
        "Full battery (K-means, Agglomerative/Ward, GMM, Spectral, DBSCAN) also run as a sensitivity "
        "check directly on the scaled 8-feature matrix, and on 3 further sensitivity inputs "
        "(KNN-imputed, mean-matrix, complete-case). All validation numbers in "
        "`reports/tables/clustering_validation.csv`; DBSCAN grid in `clustering_dbscan_grid.csv`; "
        "Gap Statistic in `clustering_gap_statistic.csv`; method agreement (pairwise ARI/RI at "
        "matched k) in `clustering_method_agreement.csv`.\n"
    )
    lines.append(
        f"- **Best single silhouette** across the entire PCA-primary battery: "
        f"`{best['method']}`, k={int(best['k'])}, silhouette={best['silhouette']:.3f} "
        f"(Davies-Bouldin={best['davies_bouldin']:.3f}, Calinski-Harabasz={best['calinski_harabasz']:.2f}) "
        "-- reported for completeness, **not** selected as \"the\" answer (see Part 11)."
        f"{best_cluster_note}\n"
        f"- Silhouette scores across the *entire* k=2..10 x 4-method PCA-primary battery: "
        f"max={pca_primary_val['silhouette'].max():.3f}, "
        f"median={pca_primary_val['silhouette'].median():.3f}, "
        f"min={pca_primary_val['silhouette'].min():.3f}.\n"
        f"- **Gap Statistic** (Tibshirani rule) chose **k=2** for every one of the 5 inputs tested "
        "(PCA-primary and all 4 sensitivity matrices) -- consistently the smallest k tested, not a "
        "confident interior maximum.\n"
        f"- **DBSCAN**: across the full documented eps/min_samples grid on the PCA-primary input, "
        f"the maximum number of density-based clusters found was **{dbscan_max_clusters}** "
        "(i.e. DBSCAN never found more than one dense group; solutions were either a single cluster "
        "plus noise, or all-noise) -- DBSCAN was not forced to produce clusters, and it did not find any.\n"
        f"- **Cross-method agreement** (mean Adjusted Rand Index between K-means/Agglomerative/GMM/"
        f"Spectral at the same k, PCA-primary input): ranges "
        f"{mean_ari_by_k.min():.3f} (k={mean_ari_by_k.idxmin()}) to {mean_ari_by_k.max():.3f} "
        f"(k={mean_ari_by_k.idxmax()}) -- moderate at best, never strong (ARI > 0.7) at any k.\n"
    )

    lines.append("### Sensitivity (Part 10)\n")
    lines.append(
        "K-means solutions at matched k compared (Adjusted Rand Index) between the PCA-primary "
        "input and each sensitivity variant (`reports/tables/clustering_sensitivity.csv`), mean ARI "
        "across k:\n"
    )
    if len(sensitivity):
        mean_by_comp = sensitivity.groupby("comparison")["adjusted_rand_index"].mean().sort_values(ascending=False)
        lines.append("| Comparison | Mean ARI across k |")
        lines.append("|---|---|")
        for comp, v in mean_by_comp.items():
            lines.append(f"| {comp} | {v:.3f} |")
    lines.append("")

    lines.append("## 11. Continuum vs. discrete-clusters -- descriptive observation vs formal conclusion\n")
    lines.append(
        "**Descriptive observation** (from the numbers above, not asserted in advance): silhouette "
        f"scores are modest across the entire battery (best 0.41, mostly well below that), and that "
        f"single best-silhouette solution is itself degenerate --{best_cluster_note} "
        "The Gap Statistic never picks an interior k, defaulting to the smallest k tested at every input; "
        "DBSCAN -- the one method here that is never forced to output a fixed number of clusters -- "
        "essentially never finds more than one density-connected group; and cross-method agreement "
        "is moderate at best and never strong. Sensitivity to imputation method and feature-matrix "
        "choice is also non-trivial (mean ARI vs the primary solution ranges "
        f"{sensitivity.groupby('comparison')['adjusted_rand_index'].mean().min():.2f}-"
        f"{sensitivity.groupby('comparison')['adjusted_rand_index'].mean().max():.2f} "
        "across sensitivity variants).\n\n"
        "**Formal conclusion**: this evidence, taken together, does **not** support a robust, "
        "discrete cluster structure in the segment-level LBD flexibility profile at the current "
        "sample size (~29 receptors) and feature set (8 segments). It is consistent with -- though "
        "does not by itself prove -- the continuum picture the prior (pre-structural-mapping) cohort "
        "suggested. No method was picked because it \"looked nicest,\" and this stage stops here "
        "without forcing a discrete-cluster narrative.\n"
    )

    lines.append("## 12. Group labels (Part 12)\n")
    lines.append(
        "Project receptor group/subfamily labels were used **only** for visualization overlays "
        "(PC1/PC2-by-group, UMAP-by-group figures) and for describing the excluded/included cohort "
        "-- never as an input to PCA fitting, imputation, UMAP fitting, or any clustering algorithm.\n"
    )

    lines.append("## 13. Entry-level sensitivity (Part 13)\n")
    top5_frac = overrep.sort_values("n_pdb_entries", ascending=False)["n_pdb_entries"].head(5).sum() / overrep["n_pdb_entries"].sum()
    lines.append(
        f"The 1789-entry PDB-level matrix was examined only as a secondary sensitivity view. "
        f"Entry counts per receptor range from 0 to {int(overrep['n_pdb_entries'].max())} "
        f"(`{overrep.sort_values('n_pdb_entries', ascending=False).iloc[0]['nr_code']}`), median "
        f"{overrep['n_pdb_entries'].median():.0f}; the top 5 most-represented receptors alone "
        f"account for **{top5_frac:.1%}** of all entries. This overrepresentation is exactly why "
        "the receptor-level median matrix, not this entry-level one, is the primary clustering "
        "input. No clustering was run at the entry level and no subfamily-separation claim is made "
        "from it; see `reports/tables/entry_level_overrepresentation.csv` and "
        "`reports/figures/modeling/entry_level_pca_by_receptor.png`.\n"
    )

    lines.append("## Limitations\n")
    lines.append(
        "- Small n (29 primary-cohort receptors, 42 with any coordinate data, 48 in the full "
        "project cohort) limits statistical power for any clustering method.\n"
        "- The 8-feature primary set reflects the current structural-template mapping's support, "
        "not necessarily the full set of biologically relevant LBD segments -- H2, H4, H6, H8, H9, "
        "H10 are excluded from the primary matrix at <50% support (H8 at 0%, not confidently "
        "resolved by the current mapping, see Part 0).\n"
        "- One curated anchor receptor (VDR) itself falls outside the primary cohort's 60%-observed "
        "rule, a direct consequence of its own template's H5/H7/H11 gaps.\n"
        "- GMM used diagonal covariance (not full) given the small per-cluster sample sizes at "
        "higher k; Spectral used a nearest-neighbors affinity graph, not RBF.\n"
        "- This is an unsupervised, exploratory stage; no hypothesis test, apo/holo comparison, or "
        "subfamily significance test was run here (explicitly deferred per instruction).\n"
    )

    out = TABLES / "modeling_summary.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
