# Phase 2 methodology pre-flight (corrected v2 foundation)

> **Superseded (2026-10-21).** Pre-flight record kept for history. The approved decisions differ from the proposals below: S2 uses the independent HELIX/SHEET/COIL consensus (not a helix/inter-helix split), and S3 grants eligibility through provenance without relabelling `review_required` as HIGH. The GR peptide statement in S1 was corrected. Current method: `README.md`.

Nothing in Phase 2 has been run. This note re-derives each old analysis from the v2
data using the previously declared rules, and lists **three methodology decisions
(S1–S3)** that must be approved before any test is run.

## Inputs common to all analyses

Phase 2 re-runs the legacy scripts unchanged in the v2 sandbox. They read **only two
things** from the superseded v3 stage, both now replaced by professor-element
equivalents:

- the segment map: `nr_lbd_structural_segment_map.csv`, which in v2 holds professor
  elements
- the element features from script 36 (`receptor_segment_*`, `pdb_segment_*`)

The literature-based `reports/tables/functional_site_definitions.csv` (Role D, UniProt
numbering) is unaffected.

## Per analysis

| Analysis (scripts) | Old input | Policy kept (declared rules) | Must be re-derived (old v3 consequences) | Test still appropriate? | Corrected n (derived, not hard-coded) | Missing data | BH-FDR family | Robustness |
|---|---|---|---|---|---|---|---|---|
| **Feature support + matrix** (48, 49) | v3 segment matrix | features = elements with ≥50% receptor support; cohort = receptors with ≥60% of features observed; median imputation (primary), KNN and complete-case sensitivity (complete-case only if n ≥ 10); standardize after imputation | feature set and cohort | yes | **10 features** (H1, H3, H4, H5, H6, H7, H8, H9, H10/H11, H12; v1: 8); **cohort 41** (v1: 29; only SHP excluded, as it has no professor map); complete-case n = 31 | as before | n/a | KNN, complete-case, mean matrix |
| **PCA / UMAP** (50, 51) | primary scaled matrix | PCA on scaled features; PCs to ≥80% variance; UMAP visualization only | n, features | yes | 41 × 10 | as above | n/a | UMAP stability |
| **Clustering** (52, 53) | PCA scores | full battery, k = 2..min(10, n−1); every metric reported; fixed seeds; entry-level sensitivity | input dimensions | yes | 41 receptors | as above | n/a | input/imputation/threshold/median-mean |
| **Group tests** (55–57, 59, 60) | receptor × segment medians | Kruskal–Wallis per element; ≥3 receptors per group, ≥2 eligible groups; ε² with bootstrap CI; median primary, mean sensitivity | element list and eligibility | yes | **11 eligible elements** under the primary grouping (the test count follows from both schemes and the rules) | available-case | each scheme × statistic family, as before | LOO receptor, mean sensitivity |
| **Atlas** (58) | receptor segment profiles | descriptive median, bootstrap CI; well-supported = ≥10 receptors | element list | yes | 16 elements; 11 with ≥10 receptors (H2 has 15) | available-case | none (descriptive) | mean sensitivity |
| **Apo/holo** (63–71) | entry-level segment features | 5/5 primary, 3/3 exploratory (fixed before testing); Mann–Whitney per receptor × element; within-receptor FDR primary, global FDR secondary; LOO-PDB | eligibility, test set | yes, **but see S1** | eligibility derived from v2 ligand states | no imputation | within-receptor; global (secondary) | LOO-PDB, QC balance, redundancy |
| **Functional sites** (72–92) | v3 H12 and segment proxies; geometric pocket | pocket = ≥25% contact fraction across ≥3 holo entries (10%/50% sensitivity); site vs whole-LBD background (paired Wilcoxon, ≥5 pairs); superfamily one-sample Wilcoxon; matched background; interpretation rules (70% preservation) | AF-2 = professor H12; groove/dimer proxies from professor elements; all site n | yes, **but see S2 and S3** | derived after definitions | no imputation | within-receptor (primary) and global; across the 4 site classes at superfamily level | LOO-PDB/receptor, threshold sensitivity, definition tiers |
| **Final numbers / audits** (93–95) | all of the above | recomputed from v2 outputs | all values | n/a | — | — | — | — |

## Technical adaptations (vocabulary only; no change of method)

These would be applied as logged patches to **sandbox copies** of the legacy scripts.
The worktree scripts stay byte-identical.

- **A1 (48):** the "core vocabulary" list H1–H12 + β-sheet → the professor element
  vocabulary (`ELEMENT_ORDER`). It only drives a reporting table.
- **A2 (58):** remove the assertion that H8 never appears. That was a v3-specific
  invariant; the professor map defines H8.
- **A3 (66, 71):** add interpretation names for the professor labels H2′, H3n, H10/H11
  and H11′.
- **A4 (73):** dimerization-interface segment proxy `["H9","H10","H11"]` →
  `["H9","H10","H11","H10/H11"]`. The professor keeps H10/H11 combined except in RORα.
  The coactivator proxy H3/H4/H5 is unchanged.
- **A5 (94):** the "unresolved H8" count becomes an observed count; script 95 compares
  against v2 numbers, not the v1 report text.
- **A6:** RORα has separate professor H10 and H11 and no combined "H10/H11". Under the
  unchanged rules its H10/H11 feature is therefore missing (median-imputed) in the
  modeling matrix. No label is altered. Reported, not "fixed".

## Decisions needed before Phase 2 is run

**S1. Apo/holo eligibility counts peptide fragments as APO structures.**
- Every APO entry with TRACE LBD coverage is a short receptor-derived peptide bound to
  another protein, not an LBD structure:
  - GR: 6YMO, 6YO8, 6YOS (all 3 of GR's APO entries)
  - ERRγ: 6Y58
  - RORγ: 4WQP, 5APK, 7QP4
- In v1 they mostly had no usable segment rows. Under the professor map their few
  residues fall inside professor elements, so **GR would become apo/holo-eligible
  (3 APO / 30 HOLO) purely through peptides**. That contradicts the known fact that no
  apo GR LBD crystal structure exists.
- **Smallest defensible change:** count and test apo/holo only for entries whose
  representative chain is an LBD construct. That means excluding the project's
  existing `lbd_mapping_category == TRACE_LBD`, applied symmetrically to APO and HOLO.
  No new threshold.
- **Effect:** GR is not testable (as the professor stated). ERRγ keeps 5 APO (primary
  tier). RORγ keeps 4 APO (exploratory tier).

**S2. Matched-background structural classes (script 85, tier 2).**
- Tier 2 matched site residues to background residues of the same structural class
  (HELIX/SHEET/LOOP/LOOP_UNCERTAIN). Those classes came from the v3/DSSP mapping.
- The professor table defines only helices. Recreating SHEET/LOOP would need a
  non-professor structural assignment (v3/DSSP), which the ground-truth decision
  forbids mixing in.
- **Smallest defensible change:** tier 2 uses the two classes the professor map
  supports, professor helix vs derived inter-helix. Tier 1 (same professor element) and
  tier 3 (whole LBD) are unchanged.
- Consequence: β-sheet and loop residues form one "inter-helix" background, a coarser
  match for the pocket and groove than before.

**S3. Confidence tier of professor-derived site definitions (scripts 73 → 84, 85, 87).**
- The old pipeline graded segment-derived site definitions by the v3 mapping
  confidence (HIGH/MODERATE/LOW). Script 87's confirmatory reanalysis keeps only
  HIGH/MODERATE.
- In the professor map the corresponding field carries the professor status. Left as
  is, **every professor-defined AF-2/H12 definition would silently drop out** of the
  confirmatory analysis.
- **Proposed:** because the professor table is ground truth, professor-element-derived
  definitions get confidence **HIGH**. This includes `review_required`, per the
  instruction to use those labels. The professor status is kept in a separate column,
  and a sensitivity run excludes `review_required`-derived definitions.
- The definition **level** is unchanged: AF-2 = RESIDUE-level professor H12; groove
  and dimer = SEGMENT-level proxies.
