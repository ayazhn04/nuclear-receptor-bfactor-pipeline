# Figure Inventory

> Every figure currently under `reports/figures/`, organized by analysis stage. **Status** legend: `PRIMARY` (candidate for the main manuscript text), `SUPPLEMENTARY` (supports a primary figure or documents a sensitivity/robustness check), `QC_ONLY` (pipeline/data-quality diagnostic, not a scientific result), `DO_NOT_USE` (superseded by a later, corrected version -- kept for provenance, never presented as current). No figure has been deleted; superseded QC figures remain on disk for audit trail. No composite publication figure has been created yet (per instruction).

## Recommended manuscript-level figure groups
1. Dataset/pipeline overview -- *no dedicated figure yet; would need to be newly composed from `reports/tables/final_project_numbers.csv` and the methods description.*
2. Segment-level flexibility atlas / structural mapping -- `subfamily_stats/segment_atlas_median_ci.png`, `alignment_qc/v3_schematic_*.png`
3. PCA/UMAP + clustering validation -- `modeling/pca_*.png`, `modeling/umap_*.png`, `modeling/clustering_validation_metrics.png`
4. Subfamily/group statistical result -- `subfamily_stats/effect_size_by_segment.png`, `subfamily_stats/segment_by_group_violin.png`
5. Apo/holo receptor-specific effects -- `apo_holo/effect_forest_plot.png`, `apo_holo/residue_profile_NR1C3.png`
6. Functional-site matched-background result -- `functional_sites_validation/superfamily_matched_forest_plot.png`, `functional_sites_validation/original_vs_matched_effect_scatter.png`

---

## `reports/figures/alignment_qc/` -- structural-segment mapping (v1 -> v2 -> v3)

| Figure | Stage | Source data | Script | Scientific question | Main result | Status |
|---|---|---|---|---|---|---|
| `01_segment_coverage_heatmap.png` | v1/v2 segment mapping | `nr_lbd_structural_segment_map.csv` (v2-era) | `scripts/38_build_segment_qc_figures.py` | Which segments are covered, per receptor? | Early coverage heatmap, superseded by v3 mapping | DO_NOT_USE |
| `02_segment_pdb_support_by_receptor.png` | v1/v2 | same | `scripts/38_build_segment_qc_figures.py` | PDB-entry support per segment/receptor | Superseded by v3 | DO_NOT_USE |
| `03_segment_missingness_summary.png` | v1/v2 | same | `scripts/38_build_segment_qc_figures.py` | Segment missingness overview | Superseded by v3 (`segment_feature_missingness.md`) | DO_NOT_USE |
| `04_example_lbd_schematic_{VDR,PPARGAMMA,ESR1,GR,NR4A1}.png` | v1 (heuristic topology) | v1 segment map | `scripts/38_build_segment_qc_figures.py` | What does the v1 ordinal-counting mapping look like? | v1 diagnostic baseline; H3/H5/H10 mostly unresolved | DO_NOT_USE (historical baseline only, preserved for provenance at commit `9f15062`) |
| `05_segment_median_z_distributions.png` | v1/v2 | v2 feature matrix | `scripts/38_build_segment_qc_figures.py` | Per-segment Z distribution | Superseded by the v3 atlas figures | DO_NOT_USE |
| `v2_schematic_{VDR,PPARgamma,ESR1,GR,NR4A1}.png` | v2 (structural alignment, no template-order fix) | v2 segment map | `scripts/43_validate_and_rebuild_features_v2.py` | v2 structural-alignment-based mapping, per anchor | H3 recovered broadly, H5-H10 still unresolved (template-construction gap) | DO_NOT_USE (superseded by v3, preserved at commit `abd2257`) |
| `v3_schematic_{VDR,PPARgamma,ESR1,GR,NR4A1}.png` | **v3 (current, curated templates + canonical-order enforcement)** | `canonical_lbd_template_segments_v2.csv`, `nr_lbd_structural_segment_map.csv` | `scripts/47_validate_v3.py` | Final curated H1-H12/beta-sheet map per anchor, visually validated for N->C order and functional-site placement | Correct, complete anchor templates; basis for all downstream analysis | **PRIMARY** (structural-mapping methods figure) |

## `reports/figures/modeling/` -- PCA / UMAP / clustering

| Figure | Source data | Script | Scientific question | Main result | Status |
|---|---|---|---|---|---|
| `pca_scree.png` | `pca_variance.csv` | `scripts/50_pca_analysis.py` | How is variance distributed across PCs? | No single dominant axis (PC1=28.4%) | **PRIMARY** |
| `pca_cumulative_variance.png` | same | same | How many PCs needed for X% variance? | 5 PCs for 80% | SUPPLEMENTARY |
| `pca_pc1_pc2_by_receptor.png` | `pca_receptor_scores.csv` | same | Where do receptors fall in PC space? | Labeled scatter, no obvious separation | SUPPLEMENTARY |
| `pca_pc1_pc2_by_group.png` | same + `group` metadata (not used in fitting) | same | Does PC1/PC2 separate by project group? | No visible group separation | **PRIMARY** |
| `pca_loadings_heatmap.png` | `pca_loadings.csv` | same | Which segments load on which PCs? | beta_sheet/H7/H3 dominate PC1 | SUPPLEMENTARY |
| `umap_feature_matrix_by_{receptor,group}.png` | scaled feature matrix | `scripts/51_umap_analysis.py` | 2D visualization, feature space | Descriptive only, stable across grid (mean rho=0.75) | QC_ONLY |
| `umap_pca_representation_by_{receptor,group}.png` | PCA-reduced representation | same | 2D visualization, PCA space | Same conclusion as feature-space UMAP | QC_ONLY |
| `umap_pca_representation_by_cluster.png` | + best-silhouette cluster labels (post-hoc coloring only) | `scripts/52_clustering_analysis.py` | What does the "best" clustering solution look like? | Reveals the solution is a degenerate 1-vs-28 split | **PRIMARY** (shows why the best result is not real evidence) |
| `clustering_validation_metrics.png` | `clustering_validation.csv` | `scripts/52_clustering_analysis.py` | Do validation metrics support a particular k? | No consistent interior-k signal across methods | **PRIMARY** |
| `entry_level_pca_by_receptor.png` | `pdb_segment_feature_matrix_median.parquet` (entry-level, secondary) | `scripts/53_entry_level_sensitivity.py` | Does PDB-entry-level structure show receptor overrepresentation? | Top-5 receptors = 60% of entries -- justifies receptor-level as primary unit | SUPPLEMENTARY |

## `reports/figures/subfamily_stats/` -- subfamily testing & flexibility atlas

| Figure | Source data | Script | Scientific question | Main result | Status |
|---|---|---|---|---|---|
| `segment_atlas_median_ci.png` | `superfamily_segment_flexibility_atlas.csv` | `scripts/61_subfamily_figures.py` | Which segments are most/least flexible superfamily-wide? | H11/loop regions most flexible; H5/H11/H3 most rigid (well-supported only) | **PRIMARY** |
| `receptor_segment_heatmap.png` | `receptor_segment_flexibility_profiles.parquet` | same | Full receptor x segment picture, missingness explicit | No obvious block structure | SUPPLEMENTARY |
| `segment_by_group_violin.png` | profiles + `group` | same | Per-segment distribution by project group | No segment reaches FDR significance | **PRIMARY** |
| `effect_size_by_segment.png` | `subfamily_segment_kruskal.csv` | same | Effect size (epsilon^2) per segment, PRIMARY grouping | All effect sizes small-to-moderate, wide CIs | SUPPLEMENTARY |
| `pca_annotated_with_subfamily_results.png` | modeling-stage PCA scores (reused, not refit) + subfamily test outcome | same | Does PCA structure relate to subfamily test outcomes? | No segment reached significance; annotation only | SUPPLEMENTARY |

## `reports/figures/apo_holo/` -- apo vs holo comparison

| Figure | Source data | Script | Scientific question | Main result | Status |
|---|---|---|---|---|---|
| `apo_holo_counts_by_receptor.png` | `apo_holo_receptor_eligibility.csv` | `scripts/70_apo_holo_figures.py` | Which receptors have enough APO/HOLO coverage? | 9/42 eligible (7 PRIMARY); VDR 0 APO | **PRIMARY** |
| `receptor_segment_delta_heatmap.png` | `apo_holo_segment_effects.parquet` | same | HOLO-APO delta across receptors/segments | Sparse, receptor-specific pattern | **PRIMARY** |
| `effect_forest_plot.png` | same | same | All tested effects + bootstrap CI | 3/97 significant within-receptor, localized | **PRIMARY** |
| `per_receptor_{NR_CODE}.png` (9 files) | `pdb_segment_flexibility_features.parquet` | same | Per-receptor per-segment APO vs HOLO distributions | Receptor-specific detail behind the forest plot | SUPPLEMENTARY |
| `residue_profile_{NR_CODE}.png` (9 files) | `apo_holo_residue_effects.parquet` | same | Residue-level delta profile per eligible receptor | ESR1 H1/H3 and PPARgamma H10/pocket clusters visible | SUPPLEMENTARY (residue_profile_NR1C3 and _NR3A1 are **PRIMARY** candidates) |

## `reports/figures/functional_sites_current/` -- functional-site flexibility map

| Figure | Source data | Script | Scientific question | Main result | Status |
|---|---|---|---|---|---|
| `receptor_site_median_z_heatmap.png` | `receptor_functional_site_flexibility.parquet` | `scripts/81_functional_site_figures.py` | Flexibility per receptor x site class | AF2/H12 positive (flexible), pocket/groove negative (rigid) | **PRIMARY** |
| `receptor_site_background_delta_heatmap.png` | same | same | Site vs background per receptor | Consistent sign pattern across most receptors | SUPPLEMENTARY |
| `site_class_forest_plot.png` | `current_functional_site_superfamily_tests.csv` | same | Superfamily-level site-class effect | AF2/H12, groove, pocket all FDR-significant; dimer interface not | **PRIMARY** |
| `af2_h12_distribution.png`, `lbp_distribution.png` | `receptor_functional_site_flexibility.parquet` | same | Distribution shape across receptors | Right-/left-shifted distributions consistent with the mean effect | SUPPLEMENTARY |
| `paired_site_vs_background_NR1C3_COACTIVATOR_GROOVE.png` | `pdb_functional_site_flexibility.parquet` | same | Worked single-receptor paired example | Illustrates the paired-PDB-entry test design | SUPPLEMENTARY |
| `apo_holo_site_overlap_schematic.png` | `apo_holo_functional_site_overlap.csv` | same | Where do significant apo/holo effects sit relative to sites? | 3/3 significant apo/holo effects overlap a site (descriptive) | SUPPLEMENTARY |
| `anchor_schematic_{VDR,PPARgamma,ESR1,GR,NR4A1}.png` | `current_functional_site_definitions.parquet` | same | Where do the 4 site classes sit along the LBD, per anchor? | Visual QC of site placement (e.g. ESR1 charge clamp, H10/H11 dimer patch) | QC_ONLY |

## `reports/figures/functional_sites_validation/` -- specificity/confounding validation

| Figure | Source data | Script | Scientific question | Main result | Status |
|---|---|---|---|---|---|
| `original_vs_matched_effect_scatter.png` | `functional_site_background_sensitivity.csv` | `scripts/91_specificity_figures.py` | Do effects survive structural-context matching? | Strong diagonal correlation; 4 lost, 5 gained out of 55 | **PRIMARY** |
| `superfamily_matched_forest_plot.png` | `functional_site_superfamily_matched_reanalysis.csv` | same | Confirmatory superfamily result, matched background | AF2/H12 and LBP remain significant; groove does not | **PRIMARY** |
| `interpretation_category_counts.png` | `functional_site_interpretation_classification.csv` | same | How many site classes fall in each evidence category? | 2 SITE_SPECIFIC_SUPPORTED, 1 STRUCTURAL_CONTEXT_EXPLAINED, 1 INSUFFICIENT | **PRIMARY** |
| `effect_by_definition_confidence.png` | `current_functional_site_within_receptor_tests.csv` + tier lookup | same | Does effect size depend on definition-quality tier? | Signal present in all tiers, strongest in tier A | SUPPLEMENTARY |
| `paired_matched_example_NR1C3_COACTIVATOR_GROOVE.png` | `pdb_functional_site_matched_background.parquet` | same | Worked example of a matched-background paired test | Illustrates the matched-control design | SUPPLEMENTARY |

---

**No figure listed here has been deleted.** QC-only and DO_NOT_USE figures remain on disk under their original stage directories for full audit trail; "DO_NOT_USE" means superseded by a later corrected analysis (v1/v2 structural mapping), never that the figure is wrong for its own historical stage.
