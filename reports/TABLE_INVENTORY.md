# Table Inventory

> Every table under `reports/tables/` (76 files), categorized `MAIN` (final scientific result, manuscript-citable), `SUPPLEMENTARY` (robustness/sensitivity detail behind a MAIN table), `QC` (pipeline/data-quality diagnostic, not a scientific result), or `HISTORICAL_ONLY` (superseded pipeline version, or the Role D 1,382-structure strict-QC analysis -- preserved unmodified, never presented as current). No table has been deleted.

## MAIN -- final scientific results

| Table | Stage | One-line description |
|---|---|---|
| `final_project_numbers.csv` | **This freeze stage** | Single traceable source of truth for every manuscript number |
| `final_data_dictionary.md` | Cumulative | Column-by-column documentation of every current data table |
| `structural_segment_validation_v3.csv` | Structural mapping (v3) | Final H1-H12/beta-sheet boundaries + evidence for the 5 anchor receptors |
| `h12_mapping_audit.csv`, `h5_h10_h11_mapping_audit.csv` | Structural mapping (v3) | A-E category audit of every receptor missing a core segment |
| `pca_variance.csv` | Modeling | PCA eigenvalues/explained variance per PC |
| `clustering_validation.csv`, `clustering_dbscan_grid.csv`, `clustering_method_agreement.csv`, `clustering_sensitivity.csv` | Modeling | Full clustering battery, validation metrics, method agreement, sensitivity |
| `modeling_summary.md` | Modeling | Formal PCA/clustering conclusion (no robust discrete clusters) |
| `subfamily_segment_kruskal.csv`, `subfamily_segment_sample_sizes.csv` | Subfamily testing | Omnibus Kruskal-Wallis results and per-segment sample sizes |
| `subfamily_flexibility_summary.md` | Subfamily testing | Formal conclusion (0 FDR-significant segments) |
| `apo_holo_receptor_eligibility.csv`, `apo_holo_segment_tests.csv`, `apo_holo_receptor_summary.csv` | Apo/holo | Eligibility, segment-level test results, per-receptor summary |
| `apo_holo_summary.md` | Apo/holo | Formal conclusion (3 significant, receptor-specific, localized) |
| `current_functional_site_availability.csv`, `current_functional_site_within_receptor_tests.csv`, `current_functional_site_superfamily_tests.csv` | Functional sites (current) | Site availability, within-receptor and superfamily test results |
| `current_functional_site_summary.md` | Functional sites (current) | Formal conclusion (AF2/H12 flexible, pocket/groove rigid) |
| `functional_site_fdr_audit.csv` | Specificity validation | Independent FDR re-derivation (confirms no bug) |
| `functional_site_matched_tests.csv`, `functional_site_superfamily_matched_reanalysis.csv` | Specificity validation | Structural-context-matched-background confirmatory tests |
| `functional_site_interpretation_classification.csv` | Specificity validation | Final per-site-class evidence classification |
| `functional_site_specificity_summary.md` | Specificity validation | Formal conclusion (AF2/H12 & LBP site-specific; groove context-explained) |
| `functional_site_background_sensitivity.csv` | Specificity validation | Original vs. matched-background effect comparison, per test |

## SUPPLEMENTARY -- robustness, sensitivity, and detail tables

| Table | Supports |
|---|---|
| `apo_holo_loo_robustness.csv`, `apo_holo_overlap_confidence_recheck.csv`, `apo_holo_structural_interpretation.csv`, `apo_holo_redundancy_audit.md`, `apo_holo_pymol_ready_table.csv` | `apo_holo_summary.md` |
| `subfamily_segment_posthoc.csv`, `subfamily_atlas_loo_robustness.csv`, `subfamily_kruskal_loo_robustness.csv`, `subfamily_median_vs_mean_agreement.csv`, `subfamily_segment_kruskal_mean_sensitivity.csv` | `subfamily_flexibility_summary.md` |
| `current_functional_site_class_comparison.csv`, `atlas_functional_site_enrichment.csv`, `atlas_functional_site_enrichment_test.csv`, `apo_holo_functional_site_overlap.csv` | `current_functional_site_summary.md` |
| `functional_site_definition_sensitivity.csv`, `functional_site_loo_pdb_robustness.csv`, `functional_site_loo_receptor_robustness.csv`, `functional_site_matched_loo_pdb_robustness.csv`, `functional_site_matched_loo_receptor_robustness.csv` | `functional_site_specificity_summary.md` |

## QC -- pipeline/data-quality diagnostics (not scientific results)

| Table | Stage |
|---|---|
| `raw_bfactor_resolution_diagnostic.csv`, `stage3a_precoordinate_qc_report.md`, `stage3b_coordinate_completeness_report.md`, `stage3b_data_recovery_note.md`, `stage3b_representative_instance_design_note.md`, `stage2_schema_hygiene_note.md`, `nr_metadata_uniprot_audit.md`, `rcsb_discovery_report.md` | Stage 2/3 data acquisition |
| `analysis_foundation_missingness_report.md`, `common_lbd_alignment_readiness.md` | Early (v1-era) analysis foundation |
| `segment_feature_missingness.md` | Structural mapping (v3 missingness policy) |
| `modeling_feature_support.csv`, `modeling_matrix_candidates.csv`, `umap_stability.csv`, `clustering_gap_statistic.csv`, `entry_level_overrepresentation.csv` | Modeling stage inputs/sensitivity |
| `apo_holo_qc_balance.csv`, `apo_holo_qc_resolution_sensitivity.csv` | Apo/holo QC balance audit |
| `functional_site_confidence_coverage_flags.csv`, `functional_site_matched_control_tiers.csv` | Functional-site definition QC |
| `final_receptor_identity_audit.csv`, `final_statistical_consistency_audit.md` | **This freeze stage** |

## HISTORICAL_ONLY -- superseded pipeline versions / historical Role D analysis

| Table | Why historical |
|---|---|
| `functional_site_definitions.csv`, `functional_site_entry_scores.csv`, `functional_site_normalization_qc.csv`, `functional_site_receptor_scores.csv` | **Historical Role D analysis** -- hand-curated, 5-receptor, 1,382-structure strict-QC view; byte-identical since commit `e714bbd`, never modified or overwritten by the current full-dataset functional-site analysis |
| `structural_segment_validation.csv` | v1 diagnostic baseline (ordinal helix counting); preserved at commit `9f15062` |
| `structural_segment_validation_v2.csv` | v2 structural-alignment mapping (template-construction gap on H5-H10); preserved at commit `abd2257`, superseded by v3 |

---

**Nothing in this inventory has been deleted.** HISTORICAL_ONLY tables remain fully intact and byte-identical to their originating commit; the "historical" label means *do not cite as the current result*, not *unavailable*.
