# Final Analysis Results

> **CORE SCIENTIFIC ANALYSIS COMPLETE.** This is the concise, evidence-based master record of the final scientific results for the *Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families* project, as of commit `f76f352` on `analysis/normalization-v1`. Every number below is traceable to `reports/tables/final_project_numbers.csv` and cross-checked against source tables in `reports/tables/final_statistical_consistency_audit.md`. Labels used throughout: **FACT** (a directly observed/counted quantity), **DESCRIPTIVE RESULT** (a computed summary with no formal hypothesis test), **STATISTICALLY SUPPORTED RESULT** (a formal test with a stated correction and significance criterion), **LIMITATION** (a boundary condition on interpretation). No new interpretation beyond what earlier stage summaries already state is introduced here.

---

## DATASET

- **FACT**: 48 project receptors (`config/nr_metadata.py`, professor-supplied, frozen).
- **FACT**: 2072 master receptor/PDB/polymer-entity candidates (`data/processed/structures.csv`); QC criteria are metadata flags, never deletion criteria -- nothing has been removed from this master inventory at any stage.
- **FACT**: 1840 coordinate-assessable candidates (X-ray + LBD sequence overlap; `coordinate_data_status == COORDINATE_DATA_AVAILABLE`).
- **FACT**: 42/48 receptors have any Stage 3B coordinate data ("coordinate receptors" -- the population every downstream structural analysis draws from).
- **FACT**: 3159 structural instances (chain/model observations) measured.
- **FACT**: 771,542 residue rows in the full, unfiltered residue table (`lbd_residue_bfactors_all.parquet`) -- every standardized-LBD position x every instance, no coverage filtering.
- **FACT**: 720,845 rows reached `Z_COMPUTED` under per-instance normalization (`normalized_lbd_bfactors_all.parquet`); the remainder are rows with missing B-factors or instances that failed the >=2-usable-B / nonzero-SD normalizability criteria, explicitly flagged, never dropped.

## STRUCTURAL MAPPING

- **FACT**: 39 structural segments (H1-H7, H9-H12, beta-sheet, plus named loops) are realized in the current (v3) structural-segment map across 41 receptors with at least one named position.
- **LIMITATION**: H8 has 0% mapping support anywhere in the 42-receptor coordinate cohort. This is a mapping-support statement (no anchor, reference, or complementary structure provided evidence anywhere in the cohort), **not a claim that H8 is biologically absent** -- it is never fabricated or zero-filled.
- **DESCRIPTIVE RESULT**: the primary modeling feature set (objectively derived, >=50% receptor support, no hard-coded segment list) is 8 features: H1, H3, H5, H7, H11, H12, beta_sheet, and loop_H3_H5.
- **FACT**: all 5 anchor receptors (VDR, PPARgamma, ESR1, GR, NR4A1) have complete, literature-cited H1-H12/beta-sheet boundaries (`reports/tables/structural_segment_validation_v3.csv`), visually confirmed for correct N->C order.

## PCA / CLUSTERING

- **FACT**: primary clustering cohort = 29/42 coordinate receptors (>=60% of the 8 primary features observed per receptor; an analysis-specific cohort rule, not a master-data exclusion). VDR itself is excluded from this cohort -- its own curated template lacks confident H5/H7/H11.
- **DESCRIPTIVE RESULT**: PC1 explains 28.4% of variance; no single dominant axis (5 PCs needed for 80% cumulative variance, out of 8 total).
- **STATISTICALLY SUPPORTED RESULT**: full clustering battery (K-means, Agglomerative/Ward, GMM, Spectral, DBSCAN) run on the PCA representation (primary) and 4 sensitivity matrices. Best single silhouette across the entire battery = 0.410 (GMM, k=2) -- but this solution is a **degenerate 1-vs-28 outlier split**, not a meaningful partition. DBSCAN, across a full documented eps/min_samples grid, never found more than 1 density-connected cluster. Mean cross-method Adjusted Rand Index (K-means/Ward/GMM/Spectral at matched k) = 0.35, moderate at best.
- **DESCRIPTIVE RESULT**: UMAP grid (27 settings x 2 input representations) is stable (mean Spearman rho=0.75 to a reference setting) -- used for visualization only, never as a clustering input.
- **SUPPORTED CONCLUSION**: *"Under the structurally aligned segment-level representation, there is no robust evidence for discrete nuclear-receptor flexibility clusters; the results are more consistent with continuous variation."* No method was chosen because it "looked nicest"; this conclusion was not decided before running the statistics.

## SUBFAMILY / GROUP TESTING

- **FACT**: PRIMARY grouping = project-supplied functional `group` (Nutritional sensors n=25, Steroid hormone receptors n=12, Orphan receptors n=5); SECONDARY = formal NR family-number; formal NR subfamily-letter reported but unusable (every subfamily n<=4).
- **FACT**: 19 segments testable under the PRIMARY grouping (>=2 eligible groups, each >=3 receptors), 16 under SECONDARY -- 35 total omnibus tests, available-case per segment (a receptor missing one segment still contributes to every other segment's test).
- **STATISTICALLY SUPPORTED RESULT**: 2/35 nominal (raw p<0.05); **0/35 significant after BH-FDR correction**, in either grouping scheme, under either the median or mean receptor summary (100% agreement between the two).
- **DESCRIPTIVE RESULT (atlas)**: among the 16 well-supported segments (n_receptors>=10), loop_H11_H12 and loop_H1_H3 are the most consistently flexible; H5, H11, and H3 are the most consistently rigid (bootstrap 95% CI reported for every value; unfiltered rankings dominated by thin-n segments are explicitly not used as the headline).
- **RELATIONSHIP TO CLUSTERING**: Outcome A (no segment differs after FDR) -- strengthens, does not contradict, the PCA/clustering continuum conclusion.

## APO / HOLO

- **FACT**: PRIMARY eligibility (>=5 APO and >=5 HOLO PDB entries) = 7/42 receptors; EXPLORATORY (>=3/>=3) = 9/42 (a superset). Both thresholds were fixed before any test was run.
- **FACT**: exact eligible receptors (nr_code): NR3A1 (ESR1), NR1C3 (PPARgamma), NR2B1 (RXRalpha), NR1C1 (PPARalpha), NR1I2 (PXR), NR3B3 (ERRgamma), NR4A1 (NUR77) [PRIMARY]; NR1F3 (RORgamma), NR2B2 (RXRbeta) [EXPLORATORY-only].
- **FACT**: **VDR (P11473) is explicitly INELIGIBLE for apo/holo comparison: 0 APO PDB entries, 48 HOLO PDB entries.** VDR is never removed from `structures.csv` or any master dataset, and is never plotted or tested as though it had an APO comparison.
- **FACT**: 97 receptor x segment tests run (segment-specific eligibility tiers, never a global receptor-level pass-through).
- **STATISTICALLY SUPPORTED RESULT**: 17/97 nominal; **3/97 significant within-receptor FDR** (ESR1 H1 delta=+0.656, ESR1 H3 delta=-0.145, ERRgamma loop_H6_H10 delta=+0.138); **0/97 significant under the secondary global FDR** -- these are receptor-specific, localized findings, not superfamily-wide effects.
- **ROBUSTNESS**: all 3 within-receptor-significant results retain effect direction and raw p<0.05 in **100%** of leave-one-PDB-out folds.

## FUNCTIONAL SITES

- **FACT**: 3684 site-definition rows across 41 receptors, 4 site classes (AF2/H12, ligand-binding pocket, coactivator groove, dimerization interface); every receptor x site is explicitly AVAILABLE/PARTIAL/UNRESOLVED, never zero-filled when unresolved.
- **STATISTICALLY SUPPORTED RESULT (original, whole-LBD background)**: 107 within-receptor paired Wilcoxon tests; 73/107 significant within-receptor FDR, 76/107 significant global FDR (independently re-derived and confirmed exact by a second and third BH implementation; the count discrepancy is a documented, legitimate BH family-size/rank-structure property, not a bug).
- **STATISTICALLY SUPPORTED RESULT (matched-background validation)**: 55 RESIDUE-level tests re-run against a structural-context-matched background (same-segment or same-structural-class, never the circular same-segment-proxy match); 43/55 remain significant; 49/55 preserve effect direction.
- **AF2/H12 final interpretation**: **SITE_SPECIFIC_SUPPORTED**. Significantly more flexible than background both originally and against the matched background (confirmatory superfamily reanalysis, n=20 confident receptors, q<0.0001); 93% of individually-significant original results preserve significance under matching; 100% leave-one-receptor-out stable.
- **LBP (ligand-binding pocket) final interpretation**: **SITE_SPECIFIC_SUPPORTED**. Significantly more rigid than background both originally and against the matched background (n=29, q<0.0001); 87% preservation; 100% leave-one-receptor-out stable.
- **Coactivator-groove corrected interpretation**: **STRUCTURAL_CONTEXT_EXPLAINED**. The original apparent rigidity (significant vs. whole-LBD background) does **not** survive matching against same-segment/same-class background (n=4 confident receptors, not significant) -- better attributed to surrounding helical context than to the groove residues specifically. This is a genuine correction to the face-value original-stage interpretation.
- **Dimer-interface limitation**: **INSUFFICIENT_RESIDUE_LEVEL_MAPPING**. Only 2 receptors (ESR1, GR) have HIGH/MODERATE-confidence residue-level dimerization-interface definitions -- too few to confirm or refute specificity this way; the original (non-significant) whole-LBD-background result stands undisputed but underpowered.
- **Apo/holo overlap (descriptive only, not recomputed)**: all 3 significant apo/holo segment effects overlap a functional site, and all 3 overlaps rest on RESIDUE-level, HIGH-confidence site definitions.

---

*Distinguishing fact from interpretation throughout is intentional: every STATISTICALLY SUPPORTED RESULT above states its correction and criterion; every DESCRIPTIVE RESULT is explicitly not a hypothesis test; every LIMITATION states its consequence for interpretation, not just its existence. See `reports/FINAL_LIMITATIONS.md` for the complete limitations register and `reports/FINAL_METHODS_AUDIT.md` for exact methods.*
