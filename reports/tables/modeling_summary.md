# Modeling Summary: Segment-Level PCA / UMAP / Clustering

> Continues from commit `7d280bc` (canonical LBD template mapping) with supervisor approval to proceed to PCA/UMAP/clustering. Master `structures.csv`, normalization, and the structural-segment mapping were **not** modified by this stage -- see 'Frozen inputs' below. Group/subfamily labels were **never** used to fit PCA, UMAP, imputation, or clustering; they are joined back in only for visualization/labeling after the fact.

## 0. H8 documentation correction

Where prior documentation said H8 was "genuinely absent," the wording was corrected to: H8 was not confidently resolved by the current structural-template mapping and is retained as unresolved/not available for the primary feature matrix. No H8 values are fabricated anywhere in this stage (H8 is not a column in any modeling matrix at all).

## 1. Frozen inputs (Part 1)

SHA256 hashes recorded in `data/analysis/modeling_frozen_input_hashes.json` for the five upstream files this stage is built on; none are modified by scripts 48-53 (verified by `tests/test_modeling_pipeline.py`):

- `data/analysis/receptor_segment_feature_matrix_median.parquet`
- `data/analysis/receptor_segment_feature_matrix_mean.parquet`
- `data/analysis/receptor_segment_coverage_matrix.parquet`
- `data/analysis/receptor_segment_flexibility_profiles.parquet`
- `data/analysis/nr_lbd_structural_segment_map.csv`

## 2. Feature-support audit (Part 2)

All 40 structural-segment/loop columns realized by the v3 segment map were audited (no segment names hard-coded) -- full table in `reports/tables/modeling_feature_support.csv`. Candidate feature sets (`reports/tables/modeling_matrix_candidates.csv`):

| Set | Support threshold | n features | n receptors w/ complete data | Total matrix missingness |
|---|---|---|---|---|
| A | >=70% | 2 | 27 | 20.2% |
| B | >=60% | 6 | 19 | 30.2% |
| C | >=50% | 8 | 11 | 33.9% |
| D | >=40% | 14 | 5 | 43.5% |

## 3. Primary segments selected (Part 3)

The objectively-derived **>=50%-support** candidate set ("C") was used as the initial primary candidate, per instruction, and accepted as-is after explicit review: **8 features** -- 7 core helix/beta-sheet segments (H1, H11, H12, H3, H5, H7, beta_sheet) and 1 loop segment (loop_H3_H5). This is reported exactly as computed: it is **not** the same as "7 core segments" one might expect from the Part 9 acceptance table alone -- `loop_H3_H5` independently clears 50% receptor support and is included, since segment names were never hard-coded.

**Excluded low-support segments** (< 50% of 42 coordinate receptors): H10, H2, H4, H6, H8, H9, loop_H10_H12, loop_H11_H12, loop_H1_H2, loop_H1_H3, loop_H1_H4, loop_H1_H5, loop_H2_H3, loop_H3_H3, loop_H3_H4, loop_H3_beta_sheet, loop_H4_H11, loop_H4_H9, loop_H4_beta_sheet, loop_H5_H6, loop_H5_H7, loop_H5_beta_sheet, loop_H6_H10, loop_H6_H11, loop_H6_H7, loop_H7_H10, loop_H7_H11, loop_H7_H12, loop_H9_H10, loop_beta_sheet_H11, loop_beta_sheet_H6, loop_beta_sheet_H9.

## 4. Receptor analysis cohort (Part 3)

Cohort rule: a receptor must have observed values for >=60% of the 8 primary features to be **included** in the primary clustering analysis. This is an analysis-specific rule recorded in `data/analysis/modeling_receptor_cohort.csv`; it does **not** remove any receptor from `structures.csv` or any master/normalization file.

- **29/42** coordinate receptors included.

- Missingness in the included-cohort matrix **before imputation**: 12.5%.

- **13 excluded**, all with explicit reasons (never silently dropped):

| uniprot_id | nr_code | n_observed/n_selected | reason |
|---|---|---|---|
| P11473 | NR1I1 | 4/8 | Observed only 4/8 (50%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| P10827 | NR1A1 | 3/8 | Observed only 3/8 (38%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| Q92731 | NR3A2 | 3/8 | Observed only 3/8 (38%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| P48443 | NR2B3 | 2/8 | Observed only 2/8 (25%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| P06401 | NR3C3 | 1/8 | Observed only 1/8 (12%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| P08235 | NR3C2 | 1/8 | Observed only 1/8 (12%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| P28702 | NR2B2 | 1/8 | Observed only 1/8 (12%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| P41235 | NR2A1 | 1/8 | Observed only 1/8 (12%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| P43354 | NR4A2 | 1/8 | Observed only 1/8 (12%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| Q07869 | NR1C1 | 1/8 | Observed only 1/8 (12%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| Q9Y5X4 | NR2E3 | 1/8 | Observed only 1/8 (12%) of the primary selected segments, below the 60% inclusion rule -- analysis-cohort exclusion only, not removed from structures.csv. |
| P49116 | NR2C2 | 0/8 | No observed value for any of the primary selected segments (e.g. no usable structural alignment to any anchor template) -- analysis-cohort exclusion only, not removed from structures.csv. |
| Q15466 | NR0B2 | 0/8 | No observed value for any of the primary selected segments (e.g. no usable structural alignment to any anchor template) -- analysis-cohort exclusion only, not removed from structures.csv. |

Notably, **VDR (P11473)**, one of the 5 curated anchor receptors, is itself excluded from the primary cohort: VDR's own curated template has no literature-anchored (or cross-template-rescued) evidence for H5/H7/H11, so it observes only H1/H3/beta_sheet/H12 (4/8 = 50% < the 60% inclusion rule). This is a real, correctly-computed consequence of the v3 template mapping, not a bug -- flagged explicitly rather than silently accepted.

## 5. Missing-value handling (Part 4)

PRIMARY = per-feature median imputation using only the median across the **included** cohort receptors' own observed values (never 0, never group/cluster/ligand-state information) -> `receptor_segment_matrix_primary.parquet` (29x8), with every imputed cell identifiable in `receptor_segment_matrix_primary_imputation_mask.parquet`.

Sensitivity variants also produced: (A) median-imputed mean-of-entry-median feature matrix; (B) KNN imputation (k=5, distance-weighted); (C) complete-case analysis (11 receptors, ADEQUATE).

## 6. Feature scaling (Part 5)

After imputation, each of the 8 primary segment features was standardized (mean=0, SD=1) across the 29-receptor primary cohort -- distinct from the upstream per-instance B-factor Z-normalization, which is untouched. Saved as `receptor_segment_matrix_primary_unscaled.parquet` / `_scaled.parquet`. The biological segment Z summaries themselves were never overwritten.

## 7. PCA (Part 6)

Fit on 29 receptors x 8 scaled features only (group labels never used for fitting).

| PC | Eigenvalue | Explained var. | Cumulative |
|---|---|---|---|
| PC1 | 2.350 | 28.4% | 28.4% |
| PC2 | 1.822 | 22.0% | 50.4% |
| PC3 | 1.312 | 15.8% | 66.2% |
| PC4 | 1.034 | 12.5% | 78.7% |
| PC5 | 0.679 | 8.2% | 86.9% |
| PC6 | 0.481 | 5.8% | 92.7% |
| PC7 | 0.390 | 4.7% | 97.4% |
| PC8 | 0.217 | 2.6% | 100.0% |

Components needed: **4** PCs for 70% variance, **5** for 80%, **6** for 90%, **7** for 95%. Variance is spread fairly evenly across components (PC1 = 28.4%) -- there is no single dominant axis of variation among the 8 segments.

Strongest |PC1| loadings:

| Segment | PC1 loading |
|---|---|
| beta_sheet | +0.486 |
| H7 | +0.475 |
| H3 | +0.387 |
| H1 | -0.379 |
| loop_H3_H5 | -0.370 |
| H11 | -0.231 |
| H5 | -0.176 |
| H12 | +0.155 |

## 8. UMAP stability (Part 7)

Grid: n_neighbors=[5, 10, 15], min_dist=[0.1, 0.3, 0.5], seeds=[0, 1, 2], run on both the scaled feature matrix and the PCA representation. Mean Spearman correlation of pairwise-embedding distances to a fixed reference setting = **0.748** -> **STABLE: global pattern correlates strongly (mean rho >= 0.7) across parameter settings**. Full grid in `reports/tables/umap_stability.csv`. UMAP was never chosen or reported because a setting "looked more separated" -- the default (n_neighbors=15, min_dist=0.1) middle-of-grid setting is used for the representative figures, and cluster coloring was only added after formal clustering (Part 9).

## 9-10. Clustering battery and validation (Parts 8-10)

Primary input: PCA scores retaining 5 PCs (>=80% variance), 29 receptors, k = 2..10. Full battery (K-means, Agglomerative/Ward, GMM, Spectral, DBSCAN) also run as a sensitivity check directly on the scaled 8-feature matrix, and on 3 further sensitivity inputs (KNN-imputed, mean-matrix, complete-case). All validation numbers in `reports/tables/clustering_validation.csv`; DBSCAN grid in `clustering_dbscan_grid.csv`; Gap Statistic in `clustering_gap_statistic.csv`; method agreement (pairwise ARI/RI at matched k) in `clustering_method_agreement.csv`.

- **Best single silhouette** across the entire PCA-primary battery: `gmm`, k=2, silhouette=0.410 (Davies-Bouldin=0.423, Calinski-Harabasz=4.37) -- reported for completeness, **not** selected as "the" answer (see Part 11). Cluster sizes for this solution: {0: 1, 1: 28} -- this is a degenerate 1-member-vs-rest outlier split, not a meaningful partition of the cohort.
- Silhouette scores across the *entire* k=2..10 x 4-method PCA-primary battery: max=0.410, median=0.163, min=0.079.
- **Gap Statistic** (Tibshirani rule) chose **k=2** for every one of the 5 inputs tested (PCA-primary and all 4 sensitivity matrices) -- consistently the smallest k tested, not a confident interior maximum.
- **DBSCAN**: across the full documented eps/min_samples grid on the PCA-primary input, the maximum number of density-based clusters found was **1** (i.e. DBSCAN never found more than one dense group; solutions were either a single cluster plus noise, or all-noise) -- DBSCAN was not forced to produce clusters, and it did not find any.
- **Cross-method agreement** (mean Adjusted Rand Index between K-means/Agglomerative/GMM/Spectral at the same k, PCA-primary input): ranges 0.229 (k=2) to 0.517 (k=7) -- moderate at best, never strong (ARI > 0.7) at any k.

### Sensitivity (Part 10)

K-means solutions at matched k compared (Adjusted Rand Index) between the PCA-primary input and each sensitivity variant (`reports/tables/clustering_sensitivity.csv`), mean ARI across k:

| Comparison | Mean ARI across k |
|---|---|
| pca_primary_vs_scaled_feature_matrix | 0.683 |
| pca_primary_vs_complete_case_scaled | 0.527 |
| pca_primary_vs_knn_imputed_scaled | 0.466 |
| pca_primary_vs_meanmatrix_scaled | 0.219 |

## 11. Continuum vs. discrete-clusters -- descriptive observation vs formal conclusion

**Descriptive observation** (from the numbers above, not asserted in advance): silhouette scores are modest across the entire battery (best 0.41, mostly well below that), and that single best-silhouette solution is itself degenerate -- Cluster sizes for this solution: {0: 1, 1: 28} -- this is a degenerate 1-member-vs-rest outlier split, not a meaningful partition of the cohort. The Gap Statistic never picks an interior k, defaulting to the smallest k tested at every input; DBSCAN -- the one method here that is never forced to output a fixed number of clusters -- essentially never finds more than one density-connected group; and cross-method agreement is moderate at best and never strong. Sensitivity to imputation method and feature-matrix choice is also non-trivial (mean ARI vs the primary solution ranges 0.22-0.68 across sensitivity variants).

**Formal conclusion**: this evidence, taken together, does **not** support a robust, discrete cluster structure in the segment-level LBD flexibility profile at the current sample size (~29 receptors) and feature set (8 segments). It is consistent with -- though does not by itself prove -- the continuum picture the prior (pre-structural-mapping) cohort suggested. No method was picked because it "looked nicest," and this stage stops here without forcing a discrete-cluster narrative.

## 12. Group labels (Part 12)

Project receptor group/subfamily labels were used **only** for visualization overlays (PC1/PC2-by-group, UMAP-by-group figures) and for describing the excluded/included cohort -- never as an input to PCA fitting, imputation, UMAP fitting, or any clustering algorithm.

## 13. Entry-level sensitivity (Part 13)

The 1789-entry PDB-level matrix was examined only as a secondary sensitivity view. Entry counts per receptor range from 0 to 379 (`NR1C3`), median 10; the top 5 most-represented receptors alone account for **60.3%** of all entries. This overrepresentation is exactly why the receptor-level median matrix, not this entry-level one, is the primary clustering input. No clustering was run at the entry level and no subfamily-separation claim is made from it; see `reports/tables/entry_level_overrepresentation.csv` and `reports/figures/modeling/entry_level_pca_by_receptor.png`.

## Limitations

- Small n (29 primary-cohort receptors, 42 with any coordinate data, 48 in the full project cohort) limits statistical power for any clustering method.
- The 8-feature primary set reflects the current structural-template mapping's support, not necessarily the full set of biologically relevant LBD segments -- H2, H4, H6, H8, H9, H10 are excluded from the primary matrix at <50% support (H8 at 0%, not confidently resolved by the current mapping, see Part 0).
- One curated anchor receptor (VDR) itself falls outside the primary cohort's 60%-observed rule, a direct consequence of its own template's H5/H7/H11 gaps.
- GMM used diagonal covariance (not full) given the small per-cluster sample sizes at higher k; Spectral used a nearest-neighbors affinity graph, not RBF.
- This is an unsupervised, exploratory stage; no hypothesis test, apo/holo comparison, or subfamily significance test was run here (explicitly deferred per instruction).

