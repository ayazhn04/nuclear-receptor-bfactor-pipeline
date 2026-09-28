# Current Full-Dataset Functional-Site Flexibility Map

> **NEW analysis layer**, continuing from commit `fb27a86`, built on the current full master dataset, approved per-instance Z normalization, chain -> PDB-entry-median aggregation, and the validated v3 structural-segment mapping. The historical Role D analysis (`reports/tables/functional_site_definitions.csv` and its downstream `functional_site_*_scores.csv`/`functional_site_normalization_qc.csv`, the old 1,382-structure strict-QC view) is **not modified, not overwritten, byte-identical** -- verified by test. No ClinVar/disease analysis was run.

## Site availability (Part 1)

| Site class | AVAILABLE | PARTIAL | UNRESOLVED |
|---|---|---|---|
| AF2_H12 | 28 | 0 | 14 |
| COACTIVATOR_GROOVE | 4 | 37 | 1 |
| DIMERIZATION_INTERFACE | 2 | 37 | 3 |
| LIGAND_BINDING_POCKET | 30 | 0 | 12 |

Unresolved sites are never assigned a placeholder position or a zero flexibility value anywhere downstream. Full receptor x site status: `reports/tables/current_functional_site_availability.csv`.

## Site definitions and provenance (Parts 2, 4, 5, 6)

`data/analysis/current_functional_site_definitions.parquet` -- 3684 rows, 41 receptors, 4 site classes:

- **AF2_H12**: taken directly from the current v3 canonical H12 mapping (never the old heuristic definition) -- 28 receptors.
- **LIGAND_BINDING_POCKET**: receptor-level consensus from real 5.0 A heavy-atom contact geometry to the already-classified FUNCTIONAL_LBD_LIGAND (never crystallization additives, ions, or buffer molecules -- that exclusion is inherited read-only from the master ligand-state classification), >=25% of eligible holo PDBs (PRIMARY) -- 30 receptors; sensitivity thresholds (10%/50%) also computed and reported (`data/analysis/ligand_pocket_consensus_residues.parquet`).
- **COACTIVATOR_GROOVE**: 8 tier-A (curated residue) rows, 1620 tier-C (H3/H4/H5 segment correspondence) rows. Tier B (bound coactivator-peptide contact) was not computed -- untracked polymer-partner identity outside the historical strict-QC subset -- reported as unavailable, not approximated.
- **DIMERIZATION_INTERFACE**: 5 tier-A rows, 810 tier-C (H9/H10/H11 segment correspondence) rows. Tier B (observed receptor-receptor contacts) likewise not computed for the same reason.

## Within-receptor site-vs-background tests (Part 10)

107 receptor x site paired Wilcoxon tests (>=5 paired PDB entries, PDB entry as the paired unit -- never individual residues). **78** nominal; **73** significant within-receptor FDR (PRIMARY); **76** significant under a secondary global FDR too. Full table: `reports/tables/current_functional_site_within_receptor_tests.csv`.

## Superfamily site-class results (Part 11)

| Site class | n receptors | median delta (site-background) | raw p | FDR q | Significant |
|---|---|---|---|---|---|
| AF2_H12 | 28 | +0.264 | 3.663e-05 | 4.885e-05 | **YES** |
| COACTIVATOR_GROOVE | 41 | -0.340 | 9.268e-06 | 3.707e-05 | **YES** |
| DIMERIZATION_INTERFACE | 39 | +0.091 | 0.5387 | 0.5387 | no |
| LIGAND_BINDING_POCKET | 30 | -0.353 | 3.05e-05 | 4.885e-05 | **YES** |

**AF2/H12 is significantly MORE flexible than background superfamily-wide; the coactivator groove and ligand-binding pocket are both significantly MORE RIGID than background; the dimerization interface shows no significant net effect.** This pattern is directly consistent with the canonical NR mechanism (H12 as a mobile activation-helix "lid"; ligand/coactivator contacts locally rigidifying the pocket and groove).

### Pairwise site-class comparison (Part 12)

| Site A | Site B | n paired | median A-B | effect (rank-biserial) | raw p | FDR q | Sig. |
|---|---|---|---|---|---|---|---|
| AF2_H12 | COACTIVATOR_GROOVE | 28 | +0.684 | +0.980 | 5.215e-08 | 3.129e-07 | **YES** |
| AF2_H12 | DIMERIZATION_INTERFACE | 26 | +0.409 | +0.447 | 0.04633 | 0.0556 | no |
| AF2_H12 | LIGAND_BINDING_POCKET | 22 | +0.644 | +0.968 | 3.338e-06 | 1.001e-05 | **YES** |
| COACTIVATOR_GROOVE | DIMERIZATION_INTERFACE | 39 | -0.234 | -0.472 | 0.009348 | 0.01402 | **YES** |
| COACTIVATOR_GROOVE | LIGAND_BINDING_POCKET | 30 | -0.008 | +0.015 | 0.9515 | 0.9515 | no |
| DIMERIZATION_INTERFACE | LIGAND_BINDING_POCKET | 28 | +0.367 | +0.616 | 0.003439 | 0.006878 | **YES** |

## Effect sizes

Within-receptor tests use the matched-pairs rank-biserial correlation (`r = (W+ - W-)/(W+ + W-)` on Wilcoxon signed ranks); superfamily/site-class comparisons use the same measure on receptor-level paired deltas. Values are stored alongside every test in the corresponding CSV (never reported as p-values alone).

## Relationship to the flexibility atlas (Part 14)

Universe = 16 well-supported atlas segments (n_receptors >= 10, the same threshold already fixed in the prior atlas stage); 14 of them are functional-site-associated. Median-split enrichment test (threshold fixed by the data's own median, not chosen after seeing the result): odds ratio = 1.00, Fisher's exact p = 1.000 -- **no significant enrichment of functional-site association in the more-flexible half**, because functional sites already cover 14/16 (the great majority) of well-supported segments, leaving little room for a discriminating signal. Honestly reported null result, not forced.

## Overlap with apo/holo effects (Part 13, descriptive only)

All 4 within-receptor-FDR-significant apo/holo segment results were cross-referenced (never recomputed) against the current site definitions; **3** overlap a functional site:

| Receptor | Segment | Site | Overlap fraction | apo/holo delta |
|---|---|---|---|---|
| NR3A1 | H3 | COACTIVATOR_GROOVE | 4% | -0.145 |
| NR3A1 | H3 | LIGAND_BINDING_POCKET | 25% | -0.145 |
| NR3B3 | loop_H6_H10 | LIGAND_BINDING_POCKET | 5% | +0.138 |

No causal inference is drawn from this overlap. Full table: `reports/tables/apo_holo_functional_site_overlap.csv`.

## Robustness (Part 15)

Leave-one-PDB-out (within receptor): **72/73** within-receptor-significant results retain p<0.05 in every fold. Leave-one-receptor-out (superfamily site classes): 
**0/4** site classes have significance that flips when a single receptor is removed (none). Full tables: `reports/tables/functional_site_loo_pdb_robustness.csv`, `functional_site_loo_receptor_robustness.csv`.

12/138 receptor x site rows rest on a LOW/MODERATE-confidence definition; 74/138 are SEGMENT-level-only (imprecise) definitions; both are flagged explicitly in `reports/tables/functional_site_confidence_coverage_flags.csv` and should be weighted accordingly when interpreting any single receptor's result.

## Limitations (site-definition confidence)

- Tier B evidence (bound coactivator/coregulator peptide contacts; observed receptor-receptor protein contacts) was not computed for either the coactivator groove or the dimerization interface -- this pipeline does not track polymer-partner chain identity outside the historical 1,382-structure strict-QC subset. Most receptors' coactivator-groove and dimerization-interface definitions are therefore SEGMENT-level (H3/H4/H5 or H9/H10/H11 correspondence), not exact residue-level annotations, and are labeled as such everywhere.
- AF2/H12 and dimerization-interface definitions are only AVAILABLE for the subset of receptors with a confident v3 H12/H9-H10-H11 structural mapping; the remainder are UNRESOLVED, never zero-filled.
- The ligand-binding-pocket consensus requires >=3 holo PDB entries; receptors with fewer (or none) are UNRESOLVED for that site, not assigned a homologous pocket.
- This is a descriptive structural association analysis. The strong, consistent, and leave-one-out-robust site-vs-background differences reported here (especially AF2/H12's flexibility and the pocket/groove's rigidity) are **structural associations, not a causal or mechanistic claim** -- they do not by themselves establish *why* a region is more or less flexible, only that it consistently is, across a large fraction of the superfamily.

