# Subfamily Flexibility Testing & Superfamily Segment Atlas

> Continues from commit `06f7b7e` (segment-level PCA/UMAP/clustering). Prior unsupervised conclusion (unchanged, not re-litigated here): *"Under the structurally aligned segment-level representation, there is no robust evidence for discrete nuclear-receptor flexibility clusters; the results are more consistent with continuous variation."* Structural mapping, normalization, and modeling-stage inputs are untouched by this stage. No apo/holo, ClinVar, or Role D analysis was run here.

## Grouping scheme used

Three grouping fields actually exist in the project metadata (never invented): (1) formal NR subfamily-letter (from `nr_code`, e.g. `NR3C4` -> `NR3C`) -- inadequate sample sizes everywhere (max n=4) and **not used as a formal test grouping**, reported for transparency only; (2) formal NR family-number (`NR3C4` -> `NR3`) -- used as **SECONDARY**; (3) the professor-supplied project functional `group` field (`config/nr_metadata.py`) -- the largest, most balanced groups, used as **PRIMARY**.

| Scheme | Groups (n) |
|---|---|
| PRIMARY: project `group` | {'Nutritional sensors': 25, 'Steroid hormone receptors': 12, 'Orphan receptors': 5} |
| SECONDARY: NR family-number | {'NR0': 1, 'NR1': 18, 'NR2': 10, 'NR3': 9, 'NR4': 2, 'NR5': 2} |
| Not usable: NR subfamily-letter | {'NR0B': 1, 'NR1A': 2, 'NR1B': 3, 'NR1C': 3, 'NR1D': 2, 'NR1F': 2, 'NR1H': 3, 'NR1I': 3, 'NR2A': 2, 'NR2B': 3, 'NR2C': 1, 'NR2E': 2, 'NR2F': 2, 'NR3A': 2, 'NR3B': 3, 'NR3C': 4, 'NR4A': 2, 'NR5A': 2} |

## Segment-specific sample sizes (Part 4)

39 realized structural segments (available-case per segment, never restricted to the 29-receptor primary clustering cohort -- a receptor missing one segment still contributes to every other segment's test). **19/39** segments have >=2 eligible groups (n>=3 receptors each) under the PRIMARY scheme; **16/39** under SECONDARY. Full per-segment-per-group counts in `reports/tables/subfamily_segment_sample_sizes.csv`.

## Omnibus results (Part 5-6)

Kruskal-Wallis across eligible groups, outcome = receptor-level `median_of_entry_median_z`, BH-FDR correction applied within each scheme's family of tests. **FDR-significant segments: 0 (PRIMARY), 0 (SECONDARY) -- zero in both.** 2 segment/scheme pair(s) had raw p<0.05 but did **not** survive BH-FDR and are **not** described as significant:

| Scheme | Segment | n | H | raw p | epsilon^2 [95% CI] | FDR q |
|---|---|---|---|---|---|---|
| primary_project_group | loop_H1_H3 | 17 | 4.45 | 0.0348 | 0.278 [0.023, 0.653] | 0.430 |
| secondary_nr_family_number | loop_H5_beta_sheet | 16 | 6.78 | 0.0092 | 0.452 [0.115, 0.723] | 0.148 |

Effect size = epsilon-squared = `H * (n+1) / (n^2-1)` (Tomczak & Tomczak 2014), with a stratified bootstrap 95% CI (seed=0, 2000 resamples). Full table: `reports/tables/subfamily_segment_kruskal.csv`.

## Post-hoc (Part 7)

**0 segments** survived omnibus BH-FDR in either scheme, so, per instruction, **no post-hoc pairwise tests were run or highlighted** (0 rows in `reports/tables/subfamily_segment_posthoc.csv`, by design).

## Flexibility-atlas ranking (Part 8-9)

39 segments in the full atlas (`data/analysis/superfamily_segment_flexibility_atlas.csv`); **16** are "well-supported" (n_receptors >= 10, an explicit documented threshold below which a bootstrap CI is unreliable or undefined). Z=0 means the average residue within a receptor's own instance(s) -- positive = relatively flexible, negative = relatively rigid, **within that receptor's own structures**, never an absolute cross-protein scale.

**Most consistently flexible (well-supported, highest median Z):**

| Segment | n | median Z | 95% bootstrap CI |
|---|---|---|---|
| loop_H11_H12 | 17 | 1.035 | [0.472, 1.398] |
| loop_H1_H3 | 18 | 0.751 | [0.008, 1.505] |
| loop_H5_beta_sheet | 20 | 0.312 | [-0.428, 0.778] |
| H9 | 11 | 0.202 | [-0.101, 0.614] |
| loop_H6_H7 | 17 | 0.169 | [-0.056, 0.916] |

**Most consistently rigid (well-supported, lowest median Z):**

| Segment | n | median Z | 95% bootstrap CI |
|---|---|---|---|
| H5 | 28 | -0.685 | [-0.754, -0.590] |
| H11 | 21 | -0.386 | [-0.608, -0.200] |
| H3 | 34 | -0.354 | [-0.566, -0.255] |
| loop_H7_H11 | 18 | -0.310 | [-0.437, -0.121] |
| H4 | 12 | -0.285 | [-0.485, -0.166] |

An **unfiltered** top/bottom-10% ranking (all realized segments, including thin-n ones) gives top=['loop_H9_H10', 'loop_H5_H7', 'loop_beta_sheet_H6', 'loop_H11_H12'], bottom=['loop_H6_H10', 'loop_H4_beta_sheet', 'H5', 'loop_H3_H3'] -- these are dominated by segments observed in only 1-3 receptors (undefined or extremely wide bootstrap CIs) and are **not** used as the headline ranking; the well-supported ranking above is preferred. No unresolved segment (H8, absent from this table entirely) is ever assigned a flexible/rigid label.

## Leave-one-receptor-out robustness (Part 10)

Atlas medians: 2/39 segments flagged as dominated by a single receptor (removing one receptor shifts the median by more than half the segment's own IQR, or 0.25 Z, whichever is larger) -- both are very-low-n segments already excluded from the well-supported ranking above (`reports/tables/subfamily_atlas_loo_robustness.csv`).

Nominal (raw p<0.05) omnibus results, LOO-checked (`reports/tables/subfamily_kruskal_loo_robustness.csv`):

- **loop_H1_H3** (primary_project_group): retains p<0.05 in 47% of leave-one-out folds (p range [0.0170, 0.0791]) -- fragile, depends on a specific receptor.
- **loop_H5_beta_sheet** (secondary_nr_family_number): retains p<0.05 in 100% of leave-one-out folds (p range [0.0022, 0.0200]) -- stable across folds, but still not FDR-significant.

## Sensitivity: mean vs median segment summaries (Part 13)

Repeating the identical omnibus procedure with the receptor MEAN segment summary (`mean_of_entry_mean_z`, the aggregation behind `receptor_segment_feature_matrix_mean.parquet`) instead of the median pipeline: **0 FDR-significant segments** (same as median-based: 0). 100% agreement on FDR-significance calls across all 35 segment/scheme tests; Spearman correlation of raw p-values between the two summaries = 0.74. This is the same conclusion from a different summary statistic, not a separate narrative (`reports/tables/subfamily_segment_kruskal_mean_sensitivity.csv`, `subfamily_median_vs_mean_agreement.csv`).

## Relationship to the unsupervised (clustering) conclusion (Part 11)

**Outcome A** (of the three possible outcomes named in the instructions): **no segment differs between groups after FDR correction, in either grouping scheme, under either the median or mean summary.** This directly strengthens the prior unsupervised finding -- flexibility at the structural-segment level appears broadly continuous/conserved across the project's functional groups and NR families, rather than concentrated in group-specific segment differences that clustering simply failed to detect as discrete groups. The two nominal (raw p<0.05) hits are exactly the kind of noise expected among ~35 independent tests at alpha=0.05 by chance, and one of the two does not even survive leave-one-out resampling.

## Limitations

- Small, unbalanced group sizes (Orphan receptors n=5) limit statistical power; a true small-to-moderate group difference could be missed (Type II error), not just correctly rejected.
- Available-case analysis means different segments are tested on different, only partially overlapping receptor subsets -- segment-to-segment comparisons of "who is more flexible" are not perfectly apples-to-apples.
- The formal NR subfamily-letter grouping (the most biologically precise grouping) could not be tested at all due to sample size; a coarser project grouping and NR family-number grouping were used instead.
- H8 remains unresolved by the current structural-template mapping and never receives an atlas value, omnibus test, or flexible/rigid label.
- This is a descriptive/inferential statistical stage; no causal or mechanistic claim is made about *why* any segment is more or less flexible.

