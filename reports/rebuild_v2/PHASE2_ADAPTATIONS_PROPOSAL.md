# Phase 2: S2 audit, complete adaptation table, report-script audit

> **Applied (2026-10-21).** G34, A01–A28, S1–S3 were approved and applied as logged sandbox patches in Phase 2 (`sandbox_patches_phase2.json`); A29 was not used. Sections saying "not yet applied" describe the state at proposal time. Current method and results: `README.md`.

**Status:** proposal for approval. **No Phase 2 analysis has been run.** No adaptation is
applied, apart from the two logged script-34 count guards needed for the S2 audit.
Nothing is committed.

All adaptations would be applied as **logged text patches to sandbox copies only**
(`scripts/101_rebuild_foundation_sandbox.py`, `SANDBOX_PATCHES_V2`, log
`logs/sandbox_patches_v2.json`). The working-tree legacy scripts stay byte-identical.

## 1. S2: independent local secondary-structure class (script 34 on v2)

**Guard patches:** sandbox only, logged, applied.

| Guard | v1 literal | Data-derived replacement | v1 result | v2 result |
|---|---|---|---|---|
| G34a | `assert len(coord) == 1840` | coordinate candidates == selected coordinate-pool rows | 1840 == 1840 ✓ (same as old) | 1844 == 1844 ✓ |
| G34b | `assert rep["instance_id"].nunique() == len(rep_ids) == 1840` | `… == len(rep_ids) == len(coord)` | 1840 ✓ (same as old) | 1844 ✓ |

**Freshness check:** `receptor_dssp_consensus.parquet` was regenerated in the v2 sandbox
on 2026-10-21 14:10. The sandbox previously held a stale v1 copy dated Oct 17. The new
file contains v2-only data:
- PPARα 202–238: 37 positions, median 79 voting entries; absent in v1.
- GR 759–767: 9 positions.
- TRα1 371–406: 36 positions.
- VDR 377–423: median 48 voting entries (v1: 3, because 45 entries were mis-numbered).
- 9,550 consensus positions (v1: 9,434); 276 positions exist only in v2.

**Inputs to script 34:**
- `data/processed/structures.csv`: the coordinate candidates and their representative
  instance (v2).
- `data/analysis/normalized_lbd_bfactors_all.parquet`: **only** the identity and
  position columns (`instance_id`, `pdb_id`, `label_asym_id`, `entity_label_seq_id`,
  `canonical_uniprot_position`, `is_construct_mapped`, `ca_record_present`). No B or z
  values are read.
- `data/raw/mmcif/<PDB>.cif.gz`: the deposited `_struct_conf` (helix records, any
  `HELX*` type, so 3₁₀ helices included) and `_struct_sheet_range` (strands), read with
  `scripts/utils/mmcif_parse.py:get_secondary_structure_ranges`.

**Rule:**
1. **Per entry first.** For each representative instance (one per coordinate
   candidate) and each observed residue (construct-mapped and Cα present):
   - HELIX if the residue's `label_seq_id` lies in a helix record of that chain;
   - otherwise SHEET if it lies in a strand record;
   - otherwise COIL (the deposited files have no explicit coil record).
2. **Then consensus** per receptor × UniProt position. Count one vote per
   representative instance, and take `consensus_dssp` = the state with the most votes.
   Exact ties (119 positions in v2) resolve deterministically in the order
   HELIX > SHEET > COIL (`idxmax`). `consensus_fraction` records the winning share.
   1,834 positions are non-unanimous.
3. Positions with no observed residue in any representative instance get no consensus
   row: **NO_OBSERVATION** (missing, never a class).

**Counts over the 11,653 v2 analysis-span positions:** HELIX 7,368; COIL 1,875;
SHEET 307; NO_OBSERVATION 2,103. Of the NO_OBSERVATION positions, 1,424 belong to
receptors without coordinate data or without a professor map.

| Professor class (v2 map) | COIL | HELIX | SHEET | NO_OBSERVATION |
|---|---|---|---|---|
| professor element (HELIX) | 86 | 7,037 | 3 | 174 |
| inter-helix (derived) | 1,782 | 314 | 304 | 505 |
| no professor mapping | 7 | 17 | 0 | 1,424 |

**Local reconstruction and independence:**
- All inputs are local: the mmCIF cache, plus the 3 approved files in the separate
  additions cache.
- No v3 field is used. A grep of script 34 and the mmCIF utilities finds no segment
  map, segment name, structural class or template.
- 98.8% of observed professor-element positions are locally HELIX. This is an
  implementation consistency check, not a test of the professor map.

**Finding about the old method:** the old matched-background "structural class" was
**not** an independent class. It came from the v3 map, where 2,765 positions labelled
LOOP have a deposited consensus of HELIX (largely v3's unnamed helices), and its
LOOP_UNCERTAIN class was defined by v3 bracketing.

**Verdict: yes.** The matched-background analysis can be reproduced scientifically with
this independent variable. The tier hierarchy (same element → same local class →
whole-LBD fallback), the ≥5 control residues per tier and every statistic stay as they
are. The local class vocabulary becomes HELIX/SHEET/COIL; LOOP_UNCERTAIN has no
independent analogue. The required code changes are A18 below. The matched-background
test has **not** been run.

## 2. Complete Phase 2 adaptation table

"v1 preserved" = on the frozen v1 inputs the patched script gives the same result as
the unpatched one.

| ID | Script(s) | Analysis | Old behaviour | Why wrong / incompatible in v2 | Proposed corrected behaviour | Logic change? | Scientific / statistical consequence | v1 preserved? |
|---|---|---|---|---|---|---|---|---|
| G34 | 34 | secondary-structure consensus | asserts 1840 coordinate candidates | v2 has 1,844 (approved span) | data-derived identities above (**applied**, sandbox) | NO | none | yes (verified) |
| A01 | 36 | element features (Phase 1E, already run) | one segment per position assumed in `load_segment_lookup` (only the unused `pos_to_seg` dict); features built by merges | dual positions PXR 261 (H3+H4) and CAR 344 (H11′+H12) | **no change**: the merges give each label its own row, so a dual position counts in both elements, per the professor table | NO | PXR H3 and H4 each include 261; CAR H11′ and H12 each include 344 | yes |
| A02 | 48 | feature-support audit | `CANONICAL_CORE_SEGMENTS = H1..H12 + beta_sheet` appended at 0 support; reports v3 confidence/method fractions (HIGH/MODERATE/LOW, own_canonical_template, …) | v3 vocabulary; β-sheet is not a professor element; v3 columns would be empty or misleading | vocabulary = professor `ELEMENT_ORDER`; report professor assignment-status fractions (published / homology / review_required) instead of v3 method/confidence; v3 behaviour kept when the map is v3 | NO for candidate selection (support-based); YES for the reporting columns | candidate feature sets unchanged by this patch (≥50% support rule) | yes (switch on map type) |
| A03 | 49–53 | cohort, matrix, PCA, UMAP, clustering, entry-level | rule-based; docstrings mention "29 receptors × 8 segments", "5 PCs" | docstrings only | **no change** (docstrings ignored) | NO | cohort and features derived: 10 features, 41 receptors | yes |
| A04 | 49 | matrix | features named by element | RORα has professor H10 and H11 separately, no combined H10/H11 | **no change**: RORα's H10/H11 feature is missing and median-imputed under the unchanged rule; reported. No label merged or renamed | NO | RORα contributes imputed H10/H11; H10 and H11 (1 receptor each) fail the ≥50% support rule | yes |
| A05 | 55–57, 59, 60, 61 | group tests, atlas robustness, figures | rule-based per element | none | **no change** | NO | test count derived | yes |
| A06 | 58 | atlas | `assert "H8" not in profiles` | v3 invariant; the professor map defines H8 in 41 receptors | remove the assertion | NO | H8 receives atlas values | yes (assertion true on v1) |
| A07 | 63, 64, 65, 67, 68, 69, 70 | apo/holo (**S1**) | APO/HOLO sets = `ligand_state ∈ {APO, HOLO}` | trace-coverage peptide chains are not LBD structures | one shared helper (new `scripts/utils/apo_holo_scope.py`) applied after each table load: entries with master `lbd_mapping_category == TRACE_LBD` are excluded from **both** APO and HOLO. Master and provenance untouched; 63 reports an extra column `n_excluded_TRACE_LBD` | **YES** (approved S1) | removes ERRγ 6Y58, GR 6YMO, GR 6YOS from the comparison; eligibility = 7 primary + 2 exploratory (same receptors as v1, freshly derived); GR stays untestable | no in general: on v1 inputs it removes the same kinds of entries |
| A08 | 65 | apo/holo residue tests | `seg_lookup = seg_map.set_index(position)` then `.get()` | non-unique index at PXR 261 / CAR 344: `.get` returns a Series and `if seg:` raises (verified) | annotation from `canonical_helix.position_memberships`; one test row per position, `structural_segment` = `display_label` (`"H3|H4"`), display only | NO (annotation only) | tests, p-values and FDR unchanged; one row per position (tested) | yes (v1 index unique) |
| A09 | 79 | atlas × site enrichment | `pos_to_seg.get()` on a position index | same failure at dual positions | membership lookup; **each** element of a site position is added to the site-associated set; the joined string is never added | NO | a site residue at 344 marks both H11′ and H12 (tested) | yes |
| A10 | 66 | apo/holo interpretation | `CATEGORY_MAP` knows H3/H4/H5, H9/H10/H11, H12, beta_sheet; `loop_*` → major loop; other H → "other helix" | `H10/H11` would fall to "other helix"; H2′/H3n/H11′ unknown; a joined residue label (`H3|H4`) would fall to "other helix" | add `H10/H11` → dimerization region; `H2′`, `H3n`, `H11′` → "other helix"; a joined label is split and its categories joined (`H3|H4` → coactivator/pocket) | NO statistical; YES descriptive | interpretation labels only | yes |
| A11 | 73 | site definitions: dimer proxy | `DIMERIZATION_SEGMENTS = ["H9","H10","H11"]` | the professor table combines H10/H11 in 41 receptors; unpatched, the dimer proxy would be H9 only (plus RORα's H10 and H11) | `["H9","H10","H11","H10/H11"]` | **YES** (definition content) | dimerization proxy covers the same helices as intended in v1; segment level unchanged | yes (v1 has no `H10/H11`) |
| A12 | 73 | site definitions: text | notes say "taken directly from the current v3 canonical H12 mapping" | wrong provenance text | "professor canonical H12 (final table)" | NO | none | n/a (text) |
| A13 | 73 (+74, 75 pass-through) | site definitions (**S3**) | `assignment_confidence` copied from the map; no provenance flag | professor statuses are not HIGH/MODERATE | add `definition_source = PROFESSOR_GROUND_TRUTH` for rows derived from professor elements (H12 = AF-2; H3/H4/H5 groove proxy; dimer proxy); **`assignment_confidence` keeps the professor status verbatim**; 74/75 carry the new column through (`any()` per receptor × site) | YES (new eligibility input) | enables S3 without falsifying confidence | yes (column absent on v1) |
| A14 | 84 | definition-quality tiers (**S3**) | tier A = RESIDUE & HIGH; tier B = RESIDUE, or SEGMENT & HIGH | professor H12 would land in tier B | tier A = RESIDUE & (HIGH or PROFESSOR_GROUND_TRUTH); tier B = RESIDUE, or SEGMENT & (HIGH or PROFESSOR_GROUND_TRUTH); `review_required` provenance kept in the output | **YES** | professor AF-2/H12 counted as residue-level ground truth | yes |
| A15 | 87, 89 | confirmatory matched reanalysis; matched robustness (**S3**) | keep `site_definition_confidence ∈ {HIGH, MODERATE}` | professor-derived AF-2 would be silently dropped | keep if confidence ∈ {HIGH, MODERATE} **or** `definition_source == PROFESSOR_GROUND_TRUTH`; plus a separately labelled **sensitivity** output excluding definitions whose professor status is `review_required` (never replaces the primary) | **YES** | AF-2/H12 receptors enter the primary confirmatory test | yes |
| A16 | 88 | apo/holo × site confidence labels | RESIDUE & HIGH → `RESIDUE_LEVEL_HIGH_CONFIDENCE` | professor definitions mislabelled "moderate or lower" | new category `RESIDUE_LEVEL_PROFESSOR_GROUND_TRUTH` (status not relabelled) | YES (labels only) | descriptive | yes |
| A17 | 80 | site robustness flags | `low_confidence_definition = confidence ∈ {LOW, MODERATE}` | professor statuses fall outside the list only by accident | explicit: professor ground truth is not low confidence; separate `professor_review_required` flag | YES (descriptive flag) | descriptive | yes |
| A18 | 85 (+86 reads its output) | matched background (**S2**) | candidate positions = named v3 segments only; tier-1 = same v3 segment (incl. named `loop_*`); tier-2 = same v3 `structural_class` | v3 class is label-derived; professor map leaves inter-helix positions unnamed, so pocket/groove residues in sheet/loop regions would vanish | candidate positions = all v2 analysis-span positions with an independent local class (script-34 v2 consensus: HELIX/SHEET/COIL; NO_OBSERVATION excluded); **tier-1** = same professor element (dual positions belong to both); for inter-helix residues the tier-1 unit is the derived stretch between the flanking professor elements (`INTER:H3–H4`, matching context only, never a label anywhere else, the analogue of v3's flank-named `loop_*`); **tier-2** = same independent local class; tier-3 and ≥5-control rule unchanged | **YES** (approved S2 direction) | background matched on independent structural context; site n derived | no (v1 used v3 class); a v1 rerun with this patch would use the v1 consensus |
| A19 | 90 | interpretation classification | rule-based; rationale text says "HIGH/MODERATE-confidence" | wording no longer exact | wording "eligible (HIGH/MODERATE or professor ground truth)" | NO | none | n/a (text) |
| A20 | 70 | apo/holo figures | `assert VDR n_APO == 0` | still true in v2 (DBD entries are AMBIGUOUS) | **no change** | NO | none | yes |
| A21 | 81, 91 | figures | anchors and file names by receptor / NR code; class colours | none (no element names in file names; "H10/H11" never used in a path) | **no change** | NO | none | yes |
| A22 | 93 | identity audit | scans the whole sandbox | runs in sandbox, not the frozen report | **no change** | NO | none | yes |
| A23 | 94 | final numbers | `unresolved_element_H8_receptor_count = 0` literal; `supported_conclusion = "no_robust_discrete_clusters_continuum_consistent"` literal | hard-coded v1 facts | H8: computed count of receptors with H8; conclusion field **removed** (numbers only); interpretation lives in the new table-driven summary | YES (removes a hard-coded conclusion) | no v1 verdict carried into v2 | no (by design) |
| A24 | 54, 62, 71, 82, 92 | narrative summaries | Markdown with v1 numbers and verdicts written into the code | would publish v1 conclusions as v2 results | **not run**; replaced by C-refactors (section 3) | YES (reporting) | v2 text derived from v2 tables only | n/a |
| A25 | 95 | claim-vs-table consistency audit | recomputes claims from tables and checks they appear in the Markdown summaries | sound design, but expects the old summary phrasing | **no change** to 95; the refactored summaries must emit the same claim lines for each checked quantity, so 95 then validates the v2 text against v2 tables | NO | independent check of the new summaries | yes |
| A26 | 73 (+74, 78, 82, 94) | site definitions at dual positions | segment-level proxy rows built per segment-map row | PXR 261 enters the coactivator-groove proxy through H3 **and** H4 (2 rows, one per membership, intended) | **no change** to the rows; every statistical consumer (74, 78, 85) uses position **sets**, so 261 counts once per site; summaries report definition rows **and** unique positions | NO | no double counting in any statistic | yes |
| A27 | 85 | matched background tier counts (S2 detail) | `len(tier1) >= 5` counts rows | dual positions would count twice | count unique positions | NO | at most 1 position difference at the threshold; identical on v1 | yes (v1 has no duplicates) |
| A28 | 101 export | script-34 v2 output | v2 local consensus exists only in the sandbox | 85 (A18) needs it as an input | export `receptor_dssp_consensus.parquet` to `data/rebuild_v2/analysis/` | NO | none | n/a |
| A29 (optional, **not proposed for primary**) | 36/49 | derived comparative "H10/H11 region" | RORα has separate H10 and H11 | 1 nomenclature-driven imputed cell (RORα H10/H11) | only if wanted as a sensitivity: union of professor H10+H11 for RORα in a derived comparative column; professor labels unchanged | YES (new derived feature) | support 40/42 → 41/42; imputed cells 13 → 12 of 410; features and cohort unchanged | n/a |

## 3. Report-generation scripts: hard-coded conclusion audit

| Script | Output | Hard-coded content found | Class |
|---|---|---|---|
| 54 modeling_summary | `modeling_summary.md` | "29-receptor primary cohort", "8 segments", "`loop_H3_H5` clears 50%", continuum narrative | **C** |
| 62 subfamily_summary | `subfamily_flexibility_summary.md` | "no robust evidence…", "**0 segments** survived", "all 35 segment/scheme tests", "~35 tests", "29-receptor" | **C** |
| 71 apo_holo_summary | `apo_holo_summary.md` | ESR1 "H1 and H3", PPARγ "H10 … loop_H7_H10", ERRγ "loop_H6_H10", "3 within-receptor", "0 global" | **C** |
| 82 functional_site_summary | `current_functional_site_summary.md` | "AF2/H12 significantly MORE flexible … groove and pocket MORE RIGID … dimerization no net effect" | **C** |
| 92 specificity_summary | `functional_site_specificity_summary.md` | "AF2/H12 and pocket SITE-SPECIFIC SUPPORTED … groove STRUCTURAL-CONTEXT EXPLAINED … dimerization INSUFFICIENT…" | **C** |
| 94 final_project_numbers | `final_project_numbers.csv` | H8 = 0 literal; `supported_conclusion` literal | **C** (A23) |
| 95 consistency audit | `final_statistical_consistency_audit.md` | none: every expected value recomputed from tables | **A** (works once the C-refactors keep the checked claim lines) |
| 90 interpretation_classification | classification CSV | rule-based categories; wording only | **B** (A19) |
| 66 apo_holo_structural_interpretation | interpretation CSV | category vocabulary | **B** (A10) |
| 48 feature_support_audit | support CSV / candidates | v3 vocabulary and columns | **B** (A02) |
| 50, 51, 52, 53, 61, 70, 81, 91 figures | PNGs | none in titles; only docstrings | **A** |
| 93 identity audit | identity CSV | static aliases (valid) | **A** |
| 30, 33 (Phase 1, already run) | release JSON / foundation reports | data-driven | **A** |

No script is classified D. All C-class summaries would be rewritten to **generate every
number and every named finding from the current v2 result tables**. Sentences whose
truth depends on results are emitted only when the table supports them, worded with
the project's approved cautious phrasing, e.g. "no robust evidence for discrete
clusters" only if the metrics say so.

## 4. S1 status (approved rule, counts re-derived; not yet applied)

Rule: TRACE_LBD entries stay in the master and provenance data but are excluded from
both APO and HOLO eligibility and tests.

- **Excluded entries** (normalizable APO/HOLO with TRACE_LBD): ERRγ 6Y58, GR 6YMO,
  GR 6YOS (all APO).
- **Before → after:** GR 2 → 0 APO (stays untestable); ERRγ 6 → 5 APO (stays
  PRIMARY); RORγ 4/155 unchanged (EXPLORATORY). Every other receptor is unchanged.
- **Eligible:** PRIMARY = PPARγ, NUR77, RXRα, PXR, ERα, ERRγ, PPARα; EXPLORATORY =
  RXRβ, RORγ.

**Correction to the earlier report.** The pre-flight said GR would become eligible
(3 APO / 30 HOLO) through peptides. That was wrong: it counted raw master ligand
states. Script 63 counts only entries that enter the element-level analysis, and just
2 of GR's 3 peptide "APO" entries do (6YMO, 6YOS; 6YO8 is not normalizable). **GR was
not eligible even before S1.** S1 still matters, because those peptide fragments must
not take part in an LBD comparison at all, for example in residue-level tests or
robustness runs.

## 5. S3 status (approved policy; implementation proposed, not applied)

Eligibility = legacy HIGH/MODERATE **or** `definition_source == PROFESSOR_GROUND_TRUTH`
(provenance from script 73). The professor `assignment_status`, including
`review_required`, is preserved verbatim in `assignment_confidence`; nothing is
relabelled HIGH.

The scripts involved are 73 (flag), 74/75 (pass-through), 84, 87, 88, 89, 80 and 90:
A13–A17 and A19. A sensitivity output excluding `review_required`-derived definitions
is clearly labelled and never replaces the primary analysis.

## 6. Outdated documentation (to update only after approval)

- `reports/canonical_remap/README.md` still describes the pre-decision policy
  (withheld `review_required` labels, the old TRα crosswalk, plot input on frozen
  profiles).
- `reports/rebuild_v2/PHASE2_METHODOLOGY_PREFLIGHT.md` has the GR claim corrected in
  section 4 here; it will be superseded by the final methods document.
- `reports/rebuild_v2/README.md` needs the S1/S2/S3 sections once applied.
