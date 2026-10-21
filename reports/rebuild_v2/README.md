# Corrected v2 analysis — method note (Phase 1 foundation + Phase 2 rerun)

**Status (2026-10-21):** Phase 1 (corrected foundation) and Phase 2 (corrected rerun of
every earlier analysis) are complete and validated. The frozen v1 release (`main` @
18b8fc9) is untouched: no tracked file is modified, and `data/raw/mmcif` keeps its
1,777 release files. v2 results are **baseline** results; see "Future work".

## What v2 is

One internally consistent dataset built from three corrections, defined in
[`reports/canonical_remap/README.md`](../canonical_remap/README.md):

1. **Professor canonical helix map = final ground truth** (decision of 2026-10-20).
   All 434 intervals and all labels, including `review_required`, are used exactly as
   supplied; professor statuses are kept verbatim; no v3 label is used. A position
   inside two intervals belongs to both elements (PXR 261: H3+H4; CAR 344: H11′+H12).
2. **Corrected PDB → UniProt numbering** from each chain's own sequence (0 numbering
   disagreements; the VDR +51 shift is gone).
3. **Final analysis span** = smallest continuous span covering the UniProt/PROSITE
   "NR LBD" ∪ every professor interval (mapped receptors); unmapped receptors keep
   PROSITE. **TRα is numbered on P10827-2** (TRα1), the professor's accession.

## How it was built (no frozen file can be touched)

Everything runs in a gitignored sandbox copy of the project (`data/interim/rebuild_v2/`)
with network access blocked; results are exported to `data/rebuild_v2/` and
`reports/rebuild_v2/`.

- **Phase 1** (`scripts/101`–`103`): the unmodified legacy pipeline (09–33) on the
  corrected inputs; the professor map in the legacy segment-map schema; element
  features from the unmodified script 36. Proof of mechanism: with the frozen inputs
  the sandbox reproduces 204/211 foundation files byte-for-byte (7 differ only in
  provenance fields); a second v2 build is byte-identical. External data: only the 3
  approved mmCIF files (1KB2, 1KB4, 1KB6), kept in the separate cache
  `data/interim/raw_mmcif_v2_additions/`.
- **Phase 2** (`scripts/104`–`107`): the legacy analyses 48–95 with every approved
  adaptation applied as a logged patch to the **sandbox copy** of the script
  (`sandbox_patches_phase2.json`; each patch must match its target exactly once).
  Before the run all v1 results are removed from the sandbox, and every file each
  script opens is logged; a read of anything that is not v2 foundation, a Phase 2
  output or a listed static input stops the run.
  - Fidelity: the *unpatched* legacy 48–95 on the v1 reproduction sandbox reproduce
    142/143 frozen outputs byte-for-byte (the identity audit differs only in the number
    of files it can scan in the sandbox; all receptors CLEAN).
  - The v1 report generators 54/62/71/82/92 contained v1 numbers and verdicts in their
    code; `scripts/105_phase2_v2_summaries.py` writes the same Markdown files from the
    v2 tables only, with every verdict produced by a stated rule. Script 95 checks them
    unchanged.
  - `scripts/106_validate_phase2_v2.py`: no v3 labels, no v1 conclusion text, missing ≠ 0,
    dual-label handling, count consistency, old-vs-v2 comparison.

## Phase 2 decisions (approved)

- **S1 — apo/holo scope.** A receptor × PDB entry whose every coordinate candidate is
  `TRACE_LBD` (a receptor fragment bound to another protein) stays in the master data
  but is excluded from APO **and** HOLO counts and tests (`scripts/utils/apo_holo_scope.py`).
  GR was already ineligible under the usable-data rule; S1 additionally keeps peptide
  fragments out of any LBD comparison. GR and VDR are not testable.
- **S2 — matched background** (`scripts/utils/matching_context.py`). Tier 1: same
  professor element (dual positions belong to both), or for residues outside professor
  helices the same derived inter-element stretch (`INTER:<before>-<after>`, a matching
  device only, never a biological label). Tier 2: same independent local class from the
  script-34 consensus (HELIX/SHEET/COIL; unobserved positions never match). Tier 3:
  whole-LBD fallback. ≥5 control positions, counted once each.
- **S3 — site-definition eligibility.** Professor-derived definitions are eligible for
  the primary confirmatory analyses through their provenance
  (`definition_source = PROFESSOR_GROUND_TRUTH`), regardless of the legacy HIGH/MODERATE
  filter. `review_required` stays `review_required`; excluding it is a labelled
  sensitivity only.
- Vocabulary/technical adaptations A01–A28 (professor element vocabulary, H10/H11 in the
  dimerization proxy, dual-label lookups, unique-position counting, no hard-coded v1
  conclusions); full table in `PHASE2_ADAPTATIONS_PROPOSAL.md`.
- The derived "H10/H11 region" (A29) is **not** used; only RORα has separate H10 and H11.

## Key results

| Quantity | v1 | v2 |
|---|---|---|
| receptors (project list) | 48 | 48 |
| coordinate-assessable candidates | 1,840 | 1,844 |
| chain/model instances | 3,159 | 3,166 |
| residue rows | 771,542 | 786,887 |
| normalized z values | 720,845 | 733,074 |
| modeling features (≥50% support) | 8 (v3) | 10 professor helices |
| modeling cohort | 29 | 41 (all coordinate receptors except SHP, which has no professor map) |
| imputed matrix cells | — | 13/410, all elements the professor does not define for that receptor |

Headline v2 results (details and every number in `tables/*_summary.md`):

- **Clustering:** 0/4 discrete-cluster criteria met (best silhouette 0.219); no robust
  discrete clusters, consistent with a continuum.
- **Subfamily tests:** 0 FDR-significant elements in either grouping scheme.
- **Atlas:** most flexible well-supported element H2, most rigid H8.
- **APO/HOLO** (associations, not causal): 7 PRIMARY + 2 EXPLORATORY eligible
  receptors; 11 within-receptor FDR-significant element results (all leave-one-PDB-out
  robust), 0 under the global FDR.
- **Functional sites** (baseline): AF-2/H12 more flexible than the rest of the LBD;
  coactivator groove, pocket and dimerization proxy more rigid. Specificity
  classifications: AF-2/H12 and pocket site-specific; groove explained by structural
  context; dimerization not testable at residue level.
- **TRACE_LBD site sensitivity** (`trace_lbd_site_sensitivity.md`): excluding TRACE_LBD
  entries from the functional-site analyses changes 3 receptor × site tests (ERRγ 6Y58)
  and GR's pocket denominator, and no significance call, superfamily effect or
  classification. TRACE_LBD fragments do not drive the conclusions; the primary
  analysis keeps them.
- Old vs v2, metric by metric: `old_vs_v2_comparison.csv`.

## Outputs

- `data/rebuild_v2/{inputs,manifests,processed,analysis}/`: corrected foundation and all
  Phase 2 analysis tables (same file names as v1).
- `reports/rebuild_v2/tables/`: result tables and the five v2 summaries
  (`modeling_summary.md`, `subfamily_flexibility_summary.md`, `apo_holo_summary.md`,
  `current_functional_site_summary.md`, `functional_site_specificity_summary.md`),
  `modeling_imputed_cells.csv`, `final_project_numbers.csv`,
  `final_statistical_consistency_audit.md`; S3 sensitivity files end in
  `_sensitivity_excl_review_required`; the TRACE_LBD sensitivity tables are in
  `tables/trace_lbd_site_sensitivity/`.
- `reports/rebuild_v2/figures/`: regenerated figures.
- Baseline per-receptor inputs for later work: `receptor_position_professor_profile_v2.csv`,
  `pdb_position_bfactor_profiles.parquet`, `normalized_lbd_bfactors_all.parquet`,
  `receptor_segment_flexibility_profiles.parquet` (median, IQR, coverage),
  `receptor_functional_site_flexibility.parquet`, `receptor_dssp_consensus.parquet`.

`final_project_numbers.csv` note: `receptors_with_any_discovery_candidate` (45) counts
receptors with at least one RCSB discovery candidate; the project defines 48 receptors
(`project_receptors_defined`). The v1 file called the 45 `project_receptors`.

## Future work (not done; no conclusions drawn)

- position-level IQR biology and excess-IQR analysis;
- SASA / solvent-burial control (functional-site results stay baseline until then);
- outlier-PDB investigation and apo/holo confounder dissection;
- ClinVar / disease-mutation analysis; manuscript novelty claims.

## Tests

`tests/test_rebuild_v2_foundation.py`, `tests/test_phase2_v2.py`,
`tests/test_v2_dual_label_positions.py`, `tests/test_canonical_helix_map.py`,
`tests/test_residue_numbering_crosswalk.py`; the full legacy suite also runs unchanged.
