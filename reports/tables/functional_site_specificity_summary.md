# Functional-Site Specificity & Confounding Validation

> Continues from commit `e420454`. A sensitivity/validation layer -- the existing functional-site results are **not modified**. Language throughout uses "consistent with" / "associated with"; this analysis never claims to confirm a mechanism or that ligand binding *causes* a flexibility change.

## FDR audit (Part 1)

Independently re-derived BH-FDR (a second implementation, not a call to the same code) **exactly reproduces** every stored q-value and significance call (PASS). **3** tests are global-significant but not within-receptor-significant; this is explained by BH's dependence on each family's own p-value distribution, not family size alone -- a small within-receptor family whose few members are only moderately significant produces a harsher correction than the same p-value's much better relative rank in the full 107-test global pool. No implementation bug was found; no result was changed. Full trace: `reports/tables/functional_site_fdr_audit.csv`.

## Effect of definition confidence (Part 2)

| Tier | n receptors | n tests | Sig. within-FDR | Sig. global-FDR | Median effect size |
|---|---|---|---|---|---|
| A_residue_high_confidence | 29 | 54 | 41 | 42 | -0.393 |
| B_partial_moderate | 27 | 46 | 27 | 28 | -0.473 |
| C_segment_level_proxy | 6 | 7 | 5 | 6 | -0.795 |

Tier-A (residue-level, high confidence) and tier-C (segment-level proxy) results are never combined into one interpretation without this breakdown -- see `reports/tables/functional_site_definition_sensitivity.csv` for the full per-site-class detail.

## Original vs. structurally matched background (Parts 3-6)

55 receptor x site tests have a valid structural-context-matched background (RESIDUE-level definitions only; SEGMENT-level proxies are never given a same-segment match -- that would be circular). Control hierarchy used: 35 same-segment, 20 same-structural-class. **43/55** remain significant within-receptor FDR against the matched background.

Of the 55 tests with both an original and a matched result: **49** preserve effect direction, **38** preserve significance, **4** LOSE significance under matching (explained by structural context), **5** GAIN significance (the whole-LBD background was diluting a real, more localized effect). Full table: `reports/tables/functional_site_background_sensitivity.csv`.

## Confirmatory superfamily reanalysis (Part 7)

| Site class | n receptors | median matched delta | raw p | BH q | Status |
|---|---|---|---|---|---|
| AF2_H12 | 20 | +0.438 | 3.624e-05 | 7.051e-05 | **SIGNIFICANT** |
| COACTIVATOR_GROOVE | 4 | +0.248 | 0.25 | 0.25 | not significant |
| DIMERIZATION_INTERFACE | 2 | -- | -- | -- | INSUFFICIENT_N |
| LIGAND_BINDING_POCKET | 29 | -0.258 | 4.701e-05 | 7.051e-05 | **SIGNIFICANT** |

## Interpretation classification (Part 8)

| Site class | Category | Rationale |
|---|---|---|
| AF2_H12 | **SITE_SPECIFIC_SUPPORTED** | Matched-background superfamily reanalysis remains FDR-significant (93% of individually-significant original results preserve significance). |
| COACTIVATOR_GROOVE | **STRUCTURAL_CONTEXT_EXPLAINED** | Original whole-LBD-background superfamily effect was significant, but disappears once compared against a structurally matched (same-segment/same-class) background. |
| DIMERIZATION_INTERFACE | **INSUFFICIENT_RESIDUE_LEVEL_MAPPING** | Too few HIGH/MODERATE-confidence RESIDUE-level definitions to test the matched-background superfamily reanalysis at all. |
| LIGAND_BINDING_POCKET | **SITE_SPECIFIC_SUPPORTED** | Matched-background superfamily reanalysis remains FDR-significant (87% of individually-significant original results preserve significance). |

**AF2/H12 and the ligand-binding pocket are SITE-SPECIFIC SUPPORTED**: their flexibility differences are consistent with (not merely an artifact of) the site itself, robust across structural-context matching. **The coactivator groove is STRUCTURAL-CONTEXT EXPLAINED**: its original apparent rigidity is better attributed to the surrounding helical context than to the coactivator-groove residues specifically -- an important correction to the prior stage's face-value interpretation. **The dimerization interface has INSUFFICIENT RESIDUE-LEVEL MAPPING** to be confirmed or refuted this way.

## Apo/holo overlap confidence recheck (Part 9, descriptive only)

| Receptor | Segment | Site | Definition level | Confidence | Category |
|---|---|---|---|---|---|
| NR3A1 | H3 | COACTIVATOR_GROOVE | RESIDUE | HIGH | RESIDUE_LEVEL_HIGH_CONFIDENCE |
| NR3A1 | H3 | LIGAND_BINDING_POCKET | RESIDUE | HIGH | RESIDUE_LEVEL_HIGH_CONFIDENCE |
| NR3B3 | loop_H6_H10 | LIGAND_BINDING_POCKET | RESIDUE | HIGH | RESIDUE_LEVEL_HIGH_CONFIDENCE |

All 3 significant apo/holo x functional-site overlaps rest on RESIDUE-level, HIGH-confidence site definitions -- the apo/holo link is not being overstated on top of a coarse segment-level proxy. The apo/holo statistics themselves were not recomputed.

## Robustness (Part 10)

Leave-one-PDB-out (matched-background, within receptor): **39/43** significant results retain p<0.05 in every fold. Leave-one-receptor-out (matched-background superfamily reanalysis): 
- AF2_H12: fraction of LOO folds with p<0.05 = 100%, direction preserved in 100% of folds, significance flips = False.
- COACTIVATOR_GROOVE: fraction of LOO folds with p<0.05 = 0%, direction preserved in 100% of folds, significance flips = False.
- LIGAND_BINDING_POCKET: fraction of LOO folds with p<0.05 = 100%, direction preserved in 100% of folds, significance flips = False.

## Limitations

- The confirmatory matched-background analysis is only possible for RESIDUE-level definitions; the coactivator groove and dimerization interface are mostly SEGMENT-level proxies, which is exactly why they could not be (coactivator groove) or could not usefully be (dimerization interface, n=2) confirmed this way -- this is a genuine evidentiary gap, not a negative biological result about those sites.
- Matched-background control sizes are sometimes small (minimum 5 control residues, minimum 5 paired PDB entries); receptor-level results built on thin controls should be read with that in mind even where the superfamily-level result is robust.
- This remains a descriptive structural-association analysis. Even the SITE-SPECIFIC-SUPPORTED classification for AF2/H12 and the ligand-binding pocket means the flexibility difference is **consistent with** the site being specifically involved, not that a causal or mechanistic link to ligand binding has been established.

