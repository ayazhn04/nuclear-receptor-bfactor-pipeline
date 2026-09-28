# Final Methods Audit

> Documents the final method exactly as implemented (verified against the committed scripts, not reconstructed from memory). For every method: script, input, output, critical parameter, random seed where relevant. Superseded/historical methods (v1/v2 structural mapping, Role D) are marked as such and are not the current method.

## 1. PDB / UniProt retrieval
- **Script**: `scripts/02_*` through `scripts/09_*` (Stage 1-2 discovery), `scripts/utils/rcsb_search_client.py`, `rcsb_data_client.py`, `rcsb_instance_client.py`, `rcsb_sequence_coordinates_client.py`, `uniprot_client.py`.
- **Input**: RCSB Search/Data/GraphQL APIs, UniProt REST API, professor-supplied `config/nr_metadata.py` (48 receptors).
- **Output**: `data/raw/api/`, `data/manifests/final_counts_by_receptor.csv`, master `data/processed/structures.csv`.
- **Critical parameter**: 48 project receptors; no additional inclusion filter at discovery (`data/raw/mmcif/`: 1777 locally cached mmCIF files).

## 2. Master-retention philosophy
- **Script**: `scripts/31_build_master_inventory.py`.
- **Principle**: QC criteria (resolution, R-free, completeness, identity, fusion status) are **metadata/analysis flags, never deletion criteria**. `structures.csv` retains all 2072 discovery candidates; nothing is filtered from the master by any downstream stage.
- **Output**: `data/processed/structures.csv` (2072 rows), `data/processed/lbd_residue_bfactors_all.parquet` (771,542 rows, all standardized-LBD positions x all 3159 instance/models, no coverage filtering).

## 3. Canonical LBD definition
- **Script**: Stage 3A (`scripts/1[4-9]_*`), `scripts/utils/residue_mapping.py`.
- **Method**: standardized LBD boundaries defined per receptor from UniProt domain annotation, cross-validated against construct mapping; `standardized_lbd_start/_end/_length` in `structures.csv`.
- **Critical parameter**: construct-vs-LBD mapping categories audited individually (`lbd_mapping_category`), not assumed from sequence length alone.

## 4. Positive-occupancy C-alpha definition
- **Script**: Stage 3B (`scripts/2[0-3]_*`), `scripts/utils/mmcif_parse.py::get_heavy_atoms`.
- **Definition**: a C-alpha is "observed" only if `occupancy > 0` in the deposited mmCIF; a zero-occupancy placeholder atom record does **not** count as observed (`ca_positive_occupancy`, `zero_occupancy_only_ca` flags in `lbd_residue_bfactors_all.parquet`).
- **Altloc handling**: `scripts/utils/altloc_selection.py` selects one deterministic alternate-location conformer per position (documented tie-break, never averaged).

## 5. Missing residues
- **Method**: never imputed, never silently dropped. Every standardized-LBD position has an explicit row in `lbd_residue_bfactors_all.parquet` for every instance, with `observation_status`/`mapping_status` flags explaining absence (unresolved in the deposited structure, construct boundary, etc.).

## 6. Per-instance normalization
- **Script**: `scripts/32_build_normalization_foundation.py`.
- **Formula**: `z = (B - instance_mean) / instance_population_SD (ddof=0)`, computed independently per coordinate instance (chain/model), over rows with a usable finite `selected_ca_b_iso_or_equiv`.
- **Output**: `data/analysis/normalized_lbd_bfactors_all.parquet` (771,542 rows); `normalization_status` per instance (`NORMALIZED` / `NOT_NORMALIZABLE_FEWER_THAN_2_USABLE_B` / `NOT_NORMALIZABLE_ZERO_SD`).
- **Critical parameter**: minimum 2 usable B-factor rows per instance to normalize; population (not sample) SD.

## 7. Chain -> PDB median aggregation
- **Script**: `scripts/32_build_normalization_foundation.py` (position profiles), `scripts/36_build_segment_features.py` (segment level).
- **Method**: per-position, median Z across contributing representative-instance chains within one PDB entry -> `pdb_position_bfactor_profiles.parquet`; per-segment, median of the position-level PDB-entry profile -> `pdb_segment_flexibility_features.parquet`.
- **This is the approved hierarchy used everywhere downstream**: instance Z -> PDB-entry median -> segment summary -> receptor summary (median primary, mean sensitivity). Never changed after supervisor approval.

## 8. Structure-based segment mapping (current: v3)
- **Scripts**: `scripts/44_curate_templates_v3.py` -> `scripts/47_validate_v3.py`.
- **Method**: 5 anchor receptors (VDR/P11473, PPARgamma/P37231, ESR1/P03372, GR/P04150, NR4A1/P22736) get literature-curated H1-H12/beta-sheet templates (PRIMARY + up to 2 COMPLEMENTARY reference structures each); every other receptor is CE-aligned (Combinatorial Extension, `Bio.PDB.ccealign`) against all 5 templates; boundaries refined against the target's own secondary-structure evidence (mmCIF `_struct_conf`/`_struct_sheet_range`); a boundary is **rejected, never fabricated**, when no overlapping target secondary-structure run of the expected class exists. `enforce_canonical_order()` resolves cross-segment N->C ordering conflicts.
- **Deterministic tie-breaks**: ordered anchor list (not a Python `set`), explicit `sort_values` on `[z_score, n_aligned, template_uniprot_id]`.
- **Superseded methods**: v1 (`9f15062`, ordinal helix-counting, preserved not used), v2 (`abd2257`, structural alignment with incomplete anchor templates, preserved not used).
- **Output**: `data/analysis/nr_lbd_structural_segment_map.csv` (current, live).

## 9. Segment feature construction
- **Script**: `scripts/36_build_segment_features.py`, `scripts/37_build_segment_missingness_and_alignment.py` (reused unchanged by every later stage via `importlib`).
- **Method**: per receptor x PDB x segment, `n_expected_residues`/`n_observed_residues`/`segment_coverage` (metadata, never a flexibility feature) + `mean_z`/`median_z`/etc. from the PDB-entry position profile.
- **Aggregation to receptor**: median (primary) / mean (sensitivity) across PDB entries -> `receptor_segment_feature_matrix_median.parquet` / `_mean.parquet`.

## 10. Missingness handling
- **Principle**: missing segment values are `NaN`, **never 0** (0 is a real, meaningful mean-flexibility value); no imputation for descriptive/reporting tables. For the modeling stage's clustering input specifically, PRIMARY = per-feature median imputation using only the median across the included analysis cohort's own observed values (never a global constant, never supervised information); sensitivity variants (KNN, complete-case) also computed and compared (`scripts/49_build_modeling_cohort_and_matrix.py`).

## 11. PCA
- **Script**: `scripts/50_pca_analysis.py`.
- **Input**: 29-receptor primary clustering cohort x 8 scaled segment features (feature-support >=50% threshold, >=60%-observed-per-receptor cohort rule).
- **Method**: `sklearn.decomposition.PCA`, full SVD, `n_components = min(n_samples, n_features)` (never hard-coded); group labels never used for fitting.
- **Random seed**: `random_state=0` (deterministic full-SVD solver; seed kept for solver-fallback reproducibility).

## 12. UMAP
- **Script**: `scripts/51_umap_analysis.py`.
- **Method**: `umap.UMAP`, visualization only, never used as a clustering input. Grid: `n_neighbors in {5,10,15}`, `min_dist in {0.1,0.3,0.5}`, `seeds in {0,1,2}` (27 runs per input), on both the scaled feature matrix and the PCA representation.
- **Stability metric**: Spearman correlation of pairwise-embedding distances vs. a fixed reference setting (`n_neighbors=15, min_dist=0.1, seed=0`); mean rho=0.75 (STABLE, >=0.7 threshold).

## 13. Clustering methods
- **Script**: `scripts/52_clustering_analysis.py`.
- **Methods**: K-means (`n_init=10, random_state=0`), Agglomerative/Ward (deterministic, no randomness), Gaussian Mixture (`covariance_type="diag", n_init=5, random_state=0, reg_covar=1e-4`), Spectral (`affinity="nearest_neighbors", random_state=0`), DBSCAN (deterministic, documented eps/min_samples grid, eps from pairwise-distance percentiles, `min_samples in {2,3,4,5}`).
- **k range**: `2..min(10, n_receptors-1)` for the 4 k-based methods, never hard-coded to a fixed value.
- **Primary input**: PCA scores retaining components for >=80% variance (5 PCs). Sensitivity: scaled feature matrix directly, KNN-imputed matrix, mean-feature matrix, complete-case matrix.

## 14. Validation metrics
- Silhouette, Davies-Bouldin, Calinski-Harabasz (`sklearn.metrics`) for every k-based solution with >=2 clusters; Gap Statistic (Tibshirani, custom implementation, `n_refs=20`, `seed=0`, stratified-bootstrap reference uniform datasets) for K-means; AIC/BIC for GMM; pairwise Adjusted Rand Index / Rand Index across methods at matched k.

## 15. Group statistical testing (subfamily)
- **Scripts**: `scripts/55_subfamily_group_definition.py` -> `scripts/62_subfamily_summary.py`.
- **Grouping**: PRIMARY = project-supplied functional `group` (largest, most balanced: 25/12/5); SECONDARY = formal NR family-number; formal NR subfamily-letter reported but not used (inadequate sample sizes).
- **Test**: Kruskal-Wallis (`scipy.stats.kruskal`) per structural segment, available-case (every segment tested on its own receptor set, never restricted to the 29-receptor clustering cohort); eligibility = >=3 receptors per group, >=2 eligible groups.
- **Effect size**: epsilon-squared = `H * (n+1) / (n^2 - 1)`, bootstrap 95% CI (`n_boot=2000, seed=0`, stratified by group).
- **Multiple testing**: BH-FDR across all tested segments within the grouping scheme.

## 16. BH-FDR (used identically in every stage: subfamily, apo/holo, functional-site, specificity validation)
- **Implementation**: standard Benjamini-Hochberg step-up procedure, re-derived independently in each stage's script (and a third time in the Part-1 FDR audit) -- `q_(i) = min_{j>=i} (p_(j) * n / j)`, monotonicity enforced from the largest p-value down, clipped to `[0,1]`.
- **Two-family convention** (apo/holo and functional-site stages): within-receptor BH-FDR is the PRIMARY significance criterion; a global BH-FDR across all tests in the study is reported as SECONDARY, never substituted for the primary criterion, never used to claim a result the primary criterion rejects.

## 17. Apo/holo analysis
- **Scripts**: `scripts/63_apo_holo_eligibility.py` -> `scripts/71_apo_holo_summary.py`.
- **Eligibility**: PRIMARY tier >=5 APO and >=5 HOLO PDB entries; EXPLORATORY tier >=3/>=3; both fixed before testing. Eligibility is receptor **and** segment specific (a segment can fail its own tier even if the receptor passes overall).
- **Unit**: PDB entry (never chain, never individual residue for the primary test).
- **Test**: two-sided Mann-Whitney U (`scipy.stats.mannwhitneyu`) on entry-level segment median Z, HOLO vs. APO; Cliff's delta (exact O(n*m) computation, HOLO-vs-APO direction); bootstrap 95% CI on delta_median (`n_boot=2000, seed=0`).
- **Secondary**: residue-level Mann-Whitney, BH-FDR strictly within receptor, positions never pooled across receptors by UniProt integer.

## 18. Functional-site construction (current, new analysis layer)
- **Scripts**: `scripts/72_ligand_pocket_geometry.py` -> `scripts/82_functional_site_summary.py`.
- **4 site classes**: AF2/H12 (from the current v3 H12 mapping); ligand-binding pocket (geometric, see #19); coactivator groove and dimerization interface (evidence hierarchy: tier A = curated Role-D residues read-only, tier C = H3/H4/H5 or H9/H10/H11 segment correspondence; tier B, bound-peptide/protein-protein contact, not computed -- untracked partner-chain identity outside the historical strict-QC subset).
- **Availability status**: every receptor x site is explicitly AVAILABLE / PARTIAL / UNRESOLVED; unresolved sites never receive a placeholder value.

## 19. 5 Angstrom LBP (ligand-binding pocket) definition
- **Script**: `scripts/72_ligand_pocket_geometry.py`.
- **Method**: for every holo PDB entry with an already-classified `FUNCTIONAL_LBD_LIGAND` (reused read-only from `scripts/26_nonpolymer_inventory_and_ligand_state.py`'s classification, which already excludes crystallization additives/ions/buffer molecules and water), every receptor heavy atom with occupancy>0 within **5.0 Angstrom** of any ligand heavy atom is identified (re-parsing cached mmCIF, `scipy.spatial.distance.cdist`).
- **Consensus**: receptor-level `contact_fraction` = (holo PDBs where a position contacts the ligand) / (holo PDBs where the position is observed); PRIMARY threshold >=25%, sensitivity thresholds >=10% and >=50%; minimum 3 holo PDB entries per receptor to attempt a consensus.

## 20. Matched-background validation
- **Script**: `scripts/85_matched_background_tests.py`.
- **Control hierarchy** (tried in order, minimum 5 control residues): (1) same canonical structural segment, non-site residues; (2) same structural class (HELIX/SHEET/LOOP), non-site residues, any segment; (3) original whole-LBD background (documented fallback). RESIDUE-level site definitions only; SEGMENT-level proxies never matched (would be circular).
- **Test**: paired Wilcoxon signed-rank (`scipy.stats.wilcoxon`, `zero_method="wilcox"`), PDB entry as the paired unit; rank-biserial effect size; BH-FDR within receptor + global.

## 21. Bootstrap / leave-one-out procedures (used identically throughout)
- **Bootstrap**: percentile method, `n_boot=2000`, `seed=0` (`numpy.random.default_rng(0)`), resampling with replacement within each group/receptor as appropriate; 95% CI = [2.5th, 97.5th] percentile.
- **Leave-one-PDB-out**: drop one PDB entry at a time from a receptor's paired site/background data, recompute the test; report fraction of folds retaining p<0.05 and fraction preserving effect direction.
- **Leave-one-receptor-out**: drop one receptor at a time from a superfamily-level test, recompute; report the same two fractions plus whether significance ever flips.
- All bootstrap and LOO procedures are deterministic given the fixed seed and are verified deterministic end-to-end by the test suite's pipeline-rerun tests.
