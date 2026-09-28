# APO vs HOLO Structural-Segment Flexibility Comparison

> Continues from commit `f0c085c` (subfamily flexibility testing & atlas). Prior results (unchanged, not re-litigated): no robust discrete clusters; no segment differed between receptor groups after FDR; the flexibility atlas is descriptive only. `structures.csv` membership, per-instance normalization, structural-segment mapping, clustering outputs, and subfamily statistics are all untouched by this stage. No ClinVar/disease analysis and no Role D changes were made.

## Eligibility (Part 1)

**PRIMARY (>=5 APO and >=5 HOLO PDB entries): 7/42 receptors eligible.** **EXPLORATORY (>=3/>=3): 9/42 eligible** (a strict superset of PRIMARY). Both thresholds were fixed before any test was run and never changed after seeing results. Full counts in `reports/tables/apo_holo_receptor_eligibility.csv`.

**Eligible receptors:**

| Receptor | n_APO | n_HOLO | Tier |
|---|---|---|---|
| NR3A1 (P03372) | 8 | 343 | PRIMARY |
| NR1C3 (P37231) | 88 | 291 | PRIMARY |
| NR1F3 (P51449) | 4 | 155 | EXPLORATORY |
| NR2B1 (P19793) | 12 | 82 | PRIMARY |
| NR1C1 (Q07869) | 5 | 74 | PRIMARY |
| NR1I2 (O75469) | 11 | 69 | PRIMARY |
| NR3B3 (P62508) | 6 | 27 | PRIMARY |
| NR4A1 (P22736) | 12 | 9 | PRIMARY |
| NR2B2 (P28702) | 4 | 3 | EXPLORATORY |

**VDR (P11473): 0 APO, 48 HOLO -> INELIGIBLE**, reported explicitly and never plotted or tested as though it had an APO comparison. VDR remains fully present in `structures.csv` and every master dataset.

**33 receptors excluded** from both tiers (insufficient APO or HOLO coverage); see the eligibility table for exact counts per receptor.

## Segment-level tests (Parts 3-6)

**97** receptor x segment tests run (segment-specific eligibility, same >=5/5 PRIMARY / >=3/3 EXPLORATORY tiers, applied per segment -- never a global receptor-level pass-through). **17** nominal (raw p<0.05); **3** significant within-receptor FDR (the PRIMARY significance criterion); **0** significant under the SECONDARY global FDR correction across all tests. Nominal results that do not survive within-receptor FDR are never described as significant.

**Within-receptor-FDR-significant segment results:**

| Receptor | Segment | apo_n | holo_n | delta (HOLO-APO) | Cliff's delta | raw p | q (within) | q (global) |
|---|---|---|---|---|---|---|---|---|
| NR3A1 | H1 | 8 | 343 | +0.656 | +0.660 | 0.0007 | 0.009 | 0.070 |
| NR3A1 | H3 | 8 | 343 | -0.145 | -0.558 | 0.0054 | 0.032 | 0.105 |
| NR3B3 | loop_H6_H10 | 5 | 27 | +0.138 | +0.793 | 0.0032 | 0.039 | 0.105 |

Note none of these 3 survive the SECONDARY global FDR correction (q_global > 0.05) -- they are receptor-specific, localized findings, not superfamily-wide effects. Negative delta = relatively more rigid in HOLO; positive = relatively more flexible in HOLO; both are relative, within-instance normalized comparisons, never an absolute physical claim.

## Residue-level follow-up (Part 9)

2112 position tests across the 9 segment-eligible receptors; **47** significant within-receptor FDR (BH-FDR always computed strictly within one receptor; positions are never pooled across receptors by UniProt integer). Two notable, well-supported localized clusters (never causally interpreted):

- **NR3A1**: 18 significant positions, concentrated in H1 and H3 (delta range [-1.54, 1.96]).
- **NR1C3**: 21 significant positions, concentrated in H10 (near the ligand pocket) and loop_H7_H10 (delta range [-0.39, 0.16]).

## Structural interpretation (Part 10)

Category breakdown of significant results (existing canonical nomenclature only, no new functional mapping invented):

```
level    structural_category                                                       
residue  major loop                                                                    24
         H9/H10/H11 dimerization region                                                 9
         other helix (H1/H2/H6/H7/H8, not in the named functional categories above)     7
         H3/H4/H5 coactivator/pocket region                                             6
         beta-sheet                                                                     1
segment  other helix (H1/H2/H6/H7/H8, not in the named functional categories above)     1
         H3/H4/H5 coactivator/pocket region                                             1
         major loop                                                                     1
```

Notably, 3 of PPARgamma's 9 significant H10 positions fall within 3 residues of the known ligand-pocket residue His477 (`functional_site_definitions.csv`, read-only) -- structurally sensible (H10 lines the ligand pocket), but this is an **observational association, not a causal claim** that ligand binding rigidifies H10.

**Confounders explicitly not adjusted for**: ligand identity, agonist/antagonist state, mutations, construct length, crystallization conditions, crystal packing, resolution/refinement, and coactivator/corepressor presence. Full list and per-receptor breakdown in `reports/tables/apo_holo_structural_interpretation.csv`.

## QC balance (Part 11)

**3/9** eligible receptors flagged for APO-vs-HOLO resolution/R-free/coverage imbalance (Mann-Whitney p<0.05 on at least one metric): ['RXRα/NR2B1', 'PPARγ/NR1C3', 'PPARα/NR1C1']. ERα/NR3A1 (source of the strongest significant segment-level results) is **not** flagged. PPARγ/NR1C3 (source of the H10/pocket residue-level cluster) **is** flagged for coverage imbalance, but its `loop_H7_H10` effect survives a resolution-overlap-restricted sensitivity re-run (p=0.0095 in-band vs p=0.0027 unrestricted) -- see `reports/tables/apo_holo_qc_resolution_sensitivity.csv`. No structures were deleted from master data for this check.

## Redundancy (Part 12)

No eligible receptor's HOLO group is dominated by a single ligand (>50% threshold; full breakdown in `reports/tables/apo_holo_redundancy_audit.md`). Publication/deposition-series and mutation dominance could not be assessed -- this metadata is not tracked anywhere in the pipeline.

## Leave-one-PDB-out robustness (Part 13)

All **3** within-receptor-FDR-significant results retain their effect direction and raw p<0.05 in **100%** of leave-one-PDB-out folds -- none are driven by a single deposition. Several nominal (non-FDR-significant) results are less robust (least robust: RORγ/NR1F3 H12, retaining p<0.05 in only 41% of folds); see `reports/tables/apo_holo_loo_robustness.csv` for every tested result.

## Limitations

- Only 9 receptors (7 PRIMARY, 2 EXPLORATORY-only) have enough APO+HOLO coverage to test at all -- the great majority of the 42-receptor coordinate cohort (including VDR) cannot be assessed for an apo/holo effect with current PDB coverage.
- 0 results survive the secondary global FDR correction -- the 3 within-receptor-significant findings are best described as receptor-specific/localized, not superfamily-wide.
- Observational PDB comparisons cannot separate a true ligand-binding effect from the confounders listed in Part 10; no causal claim is made.
- 3 receptors (RXRα/NR2B1, PPARγ/NR1C3, PPARα/NR1C1) show APO/HOLO QC imbalance; PPARγ/NR1C3's headline residue-level finding was checked and survives a resolution-restricted sensitivity re-run, but this was not done exhaustively for every nominal result.
- Residue-level tests (Part 9) have higher missingness than segment-level tests and are explicitly secondary.

## Relationship to prior results

This is consistent with, and adds nuance to, the prior stages: no global/superfamily-wide discrete structure was found in unsupervised clustering or in subfamily-level FDR testing, and here too, **0 apo/holo effects survive a global correction**. However, unlike the subfamily stage (which found literally zero FDR-significant segment differences in either grouping scheme), the apo/holo comparison **does** surface 3 robust, LOO-stable, receptor-specific within-receptor-significant effects (ERα/NR3A1 H1/H3, ERRγ/NR3B3 loop_H6_H10 -- two unrelated receptors, not adjacent regions of one receptor) plus a well-supported, structurally sensible residue-level cluster in PPARγ/NR1C3 near its own ligand pocket. This is most consistent with **Outcome D (receptor-specific effects)** combined with elements of **Outcome B (localized, not global, effects)** -- ligand-state flexibility responses appear to differ by receptor and by specific region, rather than being either a uniform superfamily-wide rigidification (Outcome A) or entirely absent (Outcome C). This outcome was not decided in advance of running the statistics.

