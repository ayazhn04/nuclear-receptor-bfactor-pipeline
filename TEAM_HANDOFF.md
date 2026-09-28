# Team Handoff — Nuclear Receptor Structural Dataset


> **QC criteria are metadata/analysis flags. Structures are retained in the
> master inventory rather than deleted based on coverage thresholds.**
> Coverage, resolution, R-free, method, taxonomy, fusion status etc. are
> columns you filter on yourself — nothing was removed from the master
> dataset because of them (revised per professor's requirement).

> **CORE SCIENTIFIC ANALYSIS COMPLETE** (as of commit `f76f352`). Data &
> Infrastructure, full-dataset B-factor normalization, structure-based
> canonical LBD segment mapping (v3), segment-level PCA/UMAP/clustering,
> subfamily/group statistical testing + a superfamily flexibility atlas,
> an apo-vs-holo flexibility comparison, a current full-dataset
> functional-site flexibility map, and a specificity/confounding
> validation layer for that map are all complete, tested, and frozen.
> **Not yet performed**: ClinVar/disease mutation analysis. Start from
> [`reports/FINAL_ANALYSIS_RESULTS.md`](reports/FINAL_ANALYSIS_RESULTS.md)
> for the master evidence-based record of every final result,
> [`reports/FINAL_METHODS_AUDIT.md`](reports/FINAL_METHODS_AUDIT.md) for
> exact methods, and [`reports/FINAL_LIMITATIONS.md`](reports/FINAL_LIMITATIONS.md)
> before citing anything. This document (`TEAM_HANDOFF.md`) remains the
> guided tour of every data file and analysis layer, oldest to newest.

This document distinguishes five layers, from most to least current —
**never delete or silently overwrite an earlier one**:

1. **MASTER DATA** — `structures.csv` / `lbd_residue_bfactors_all.parquet` (below).
2. **CURRENT FINAL ANALYSES** — normalization foundation, v3 segment mapping,
   PCA/UMAP/clustering, subfamily testing + atlas, apo/holo comparison,
   current functional-site map + specificity validation (each documented in
   its own section further down this file).
3. **HISTORICAL STRICT-QC VIEW** — the 1,382-structure subset (below).
4. **HISTORICAL ROLE D ANALYSIS** — `functional_site_definitions.csv` and its
   downstream scores/QC files; byte-identical since `e714bbd`, superseded
   (never overwritten) by the current functional-site analysis.
5. **DEPRECATED DIAGNOSTIC v1/v2 STRUCTURAL MAPPING** — preserved at commits
   `9f15062` (v1) and `abd2257` (v2); superseded by the current v3 mapping.

## Conceptual model — which file is which

| File | Role |
|---|---|
| [`data/processed/structures.csv`](data/processed/structures.csv) | **MASTER inventory.** All 2072 Stage 2 receptor x polymer-entity candidates, one row each, with every QC value as metadata/status flags. |
| [`data/processed/lbd_residue_bfactors_all.parquet`](data/processed/lbd_residue_bfactors_all.parquet) | **Full residue-level analysis source.** All 771,542 corrected Stage 3B residue observations (3159 measured instance/models x every standardized-LBD position), unfiltered, with raw Ca B-factors and one recommended altloc record. |
| [`data/processed/structures_primary_qc_subset.csv`](data/processed/structures_primary_qc_subset.csv) | Historical strict QC subset (1382 rows). **NOT the master dataset.** |
| [`data/processed/final_lbd_residue_map_primary_qc_subset.parquet`](data/processed/final_lbd_residue_map_primary_qc_subset.parquet) | Historical residue map for that strict subset only. |
| [`data/processed/excluded_structures.csv`](data/processed/excluded_structures.csv) | Why each candidate was outside the strict subset (multi-reason lineage; historical view). |

**Ligand state covers the whole coordinate set.** All 1840 candidates with
Stage 3B coordinate data have an APO/HOLO/AMBIGUOUS assessment (HOLO 1640,
APO 199, AMBIGUOUS 1), regardless of coverage, resolution, R-free or strict-subset
membership. The 232 candidates without coordinate data stay in the master
with `ligand_state = NOT_ASSESSABLE_NO_COORDINATE_DATA` (never labelled APO).
QC thresholds are flags, not deletion rules.

**Nested analysis views (counts are views, not deletions):**

```
2072  Stage 2 discovery candidates          = structures.csv (MASTER, all retained)
 -> 1840  X-ray + LBD-overlap candidates with Stage 3B coordinate data (has_stage3b_coordinate_data)
     -> 1382  previous strict primary-QC subset (previous_primary_qc_member)
```

**1382 is NOT the master dataset.** For example VDR (P11473) has 52
master rows / 48 with coordinate data / 0 in the strict subset (best
positive-occupancy LBD coverage 0.83) — it is fully present in
`structures.csv` and `lbd_residue_bfactors_all.parquet` and you choose your
own coverage cutoff.

To reproduce the old strict view from the master: filter
`previous_primary_qc_member == True` (or apply your own thresholds using
`coverage_ge_0XX`, `resolution_project_range_status`, `r_free_status`,
`method_is_xray`, `identity_class`, `fusion_or_chimera`).

Supporting files:

| Need | File |
|---|---|
| Ligand identity and evidence for all 1840 coordinate candidates (`ligand_state`, `ligand_annotation_status` are also columns of the master) | [`data/manifests/final_ligand_annotations.csv`](data/manifests/final_ligand_annotations.csv) |
| Every nonpolymer component in each of the 1777 coordinate-candidate PDB entries | [`data/manifests/final_nonpolymer_inventory.csv`](data/manifests/final_nonpolymer_inventory.csv) |
| Coregulator/partner annotations (strict subset) | [`data/manifests/final_polymer_partner_annotations.csv`](data/manifests/final_polymer_partner_annotations.csv) |
| Representative-instance selection for every candidate | [`data/manifests/final_representative_instance_selection.csv`](data/manifests/final_representative_instance_selection.csv) |
| Standardized LBD definitions | [`data/manifests/nr_lbd_reference.csv`](data/manifests/nr_lbd_reference.csv) |
| Column meanings | [`reports/tables/final_data_dictionary.md`](reports/tables/final_data_dictionary.md) |
| Per-receptor counts (strict subset) | [`data/manifests/final_counts_by_receptor.csv`](data/manifests/final_counts_by_receptor.csv) |
| Machine-readable release summary + checksums | [`data/manifests/final_data_infrastructure_release.json`](data/manifests/final_data_infrastructure_release.json) |

**Start with `structures.csv` (choose rows) and `lbd_residue_bfactors_all.parquet`
(get residue-level data; join on `instance_id`, use `is_representative_instance`
to get one instance per receptor x PDB x entity).**

## Analysis Foundation (normalized data, descriptive only)

> **Raw B-factor averages are descriptive QC only. Cross-structure biological
> comparisons must use within-instance normalized B-factors. QC thresholds
> remain metadata filters: no structure is removed from the master by this
> stage. The common cross-receptor alignment is still pending professor
> confirmation — do not pivot UniProt positions across receptors.**

Three layers — keep them distinct:

| Layer | Files |
|---|---|
| **RAW / MASTER DATA** | `data/processed/structures.csv` (2072 candidates), `data/processed/lbd_residue_bfactors_all.parquet` (771,542 residue rows) |
| **NORMALIZED ANALYSIS FOUNDATION** (`data/analysis/`) | `normalized_lbd_bfactors_all.parquet`, `instance_bfactor_summary.csv`, `pdb_position_bfactor_profiles.parquet`, `receptor_position_flexibility_profiles.parquet`, `receptor_flexibility_summary.csv` |
| **HISTORICAL STRICT-QC VIEW** (untouched) | `data/processed/structures_primary_qc_subset.csv`, `data/processed/final_lbd_residue_map_primary_qc_subset.parquet` |

How it is built (`scripts/32_build_normalization_foundation.py`, reports by `scripts/33_build_analysis_foundation_reports.py`):

1. **Normalization unit = one coordinate instance.** `z = (B - mean) / SD` with population SD (`ddof=0`), using the already-selected raw Ca B-factor (`selected_ca_b_iso_or_equiv`) on rows where `ca_bfactor_usable` is True. Applied to **all 3159 instances** — no filtering by coverage, resolution, R-free, ligand state or strict-subset membership. Missing LBD residues stay as explicit rows with missing B and missing z (`row_normalization_status`). An instance with <2 usable B or zero SD is kept and marked `NOT_NORMALIZABLE_*` (currently one: `1YNW.D`, RXR-alpha DBD-only entity, 0 usable B).
2. **`instance_bfactor_summary.csv`** — one row per instance with usable counts, raw-B summaries (**QC only, not a flexibility measure**), z summaries and structure metadata.
3. **`structures.csv`** gained eight descriptive `representative_*` columns (usable count, raw B mean/median/SD/IQR for QC, fraction z > 1.5 / < -1.5, normalization status) taken from each candidate's deterministic representative instance. Blank for the 232 candidates without coordinate data. Mean z is deliberately not added (0 by construction). Master membership is unchanged (2072 rows).
4. **`pdb_position_bfactor_profiles.parquet`** — representative instances only; one observation per receptor x PDB x canonical UniProt position. If several representative candidates of one PDB contribute (e.g. two entities of the same receptor) their z-values are reduced by **median** (`n_contributing_representative_candidates`), so a PDB never counts twice. PDB-level `ligand_state` is the shared state; disagreeing candidates give `AMBIGUOUS` with `ligand_state_conflict=True` and the states in `contributing_ligand_states_json`.
5. **`receptor_position_flexibility_profiles.parquet`** — one row per receptor x canonical position over the standardized LBD (all positions kept, unobserved ones with n = 0): distribution of PDB-level z (SD across PDB entries uses `ddof=1`, blank when n < 2) plus separate APO and HOLO descriptive fields and `apo_holo_*_delta` (only when both states are observed). AMBIGUOUS/not-assessable states are never assigned to APO or HOLO. **No significance tests were run.**
6. **`receptor_flexibility_summary.csv`** — all 48 receptors; receptors with no structures remain explicit rows (`profile_status` = `NO_STRUCTURES_IN_MASTER` / `NO_COORDINATE_DATA`) with blank numeric fields. Low-N receptors are flagged (`PROFILE_AVAILABLE_LOW_N`), never dropped.

Reports: `reports/tables/analysis_foundation_missingness_report.md` (coverage/missingness/APO-HOLO balance), `reports/tables/raw_bfactor_resolution_diagnostic.csv` (raw B vs resolution; descriptive only, never used to alter z), `reports/tables/common_lbd_alignment_readiness.md` (why identical UniProt integers across receptors are not homologous, what exists, and the decision requested from the professor).

**Not done by this normalization-foundation stage itself (by design)** --
PCA/UMAP/t-SNE, clustering, subfamily tests, and any cross-receptor
common-position model were all completed in later stages, documented in
their own sections further down this file (Structural-Segment Alignment,
Segment-Level Modeling, Subfamily/Group Testing, Apo vs Holo, Current
Functional-Site Flexibility Map). Role D's functional-site files were not
touched by this or any later stage.

## How residues are numbered — read this before joining anything

```
PDB atom record
    -> label_asym_id + label_seq_id      (mmCIF construct-relative numbering)
    -> RCSB Sequence Coordinates API alignment
    -> canonical_uniprot_position         (UniProt canonical, 1-based, THE join key)
```

**Never use author numbering (`auth_seq_id`) as a join key or as "the"
residue number.** Author numbering is construct-specific, inconsistent
across depositors, and sometimes offset or discontinuous. It is preserved
in `lbd_residue_bfactors_all.parquet` purely for cross-referencing against the
original PDB file — never for joining across structures. Every table in
this handoff is keyed, or can be joined, on `canonical_uniprot_position`.

## Structural-Segment Alignment (canonical H1-H12/beta-sheet/loop segments)

> **Supervisor-approved hierarchy: raw B -> per-instance Z-score -> chain/instance
> -> PDB entry (MEDIAN) -> canonical structural segment -> receptor (MEDIAN,
> primary; MEAN, sensitivity) -> feature matrix. UniProt residue numbers are
> receptor-specific and are never pivoted directly across receptors for
> comparison -- canonical structural segments (H1-H12, beta-sheet, loops) are
> the cross-receptor coordinate. Segment-level features are PRIMARY;
> residue-level aligned profiles are SECONDARY. Coverage is metadata, not a
> flexibility feature. Missing structural segments are `NaN`, never `0`
> (`0` is a real, meaningful mean-flexibility value). PCA/UMAP/t-SNE/clustering
> have NOT been run.**

### Method history: v1 (superseded) -> v2 (superseded) -> v3 (current)

**v1 (commit `9f15062`, kept as an intact diagnostic baseline, not deleted or
rewritten)** used a topology algorithm (`scripts/utils/segment_topology.py`)
that anchored H1/H12/beta-sheet at the domain termini/unique sheet, then tried
to resolve interior helices H2-H11 by checking whether the *number* of
observed helices between two anchors matched the canonical count. This failed
broadly: H2-H11 were each resolved in only 1 receptor (**MR/NR3C2 alone** --
an audit before the v2 rebuild found the original report's claim of
"GR/MR/PR" was incorrect; GR and PR's interior helices were AMBIGUOUS in v1,
not assigned). Failure of helix counting does not mean H3/H5/H10/H11 are
biologically absent or non-homologous -- DSSP/depositor annotations can split,
merge or omit pieces of the same canonical helix, and counting cannot tell the
difference. v1 is preserved for comparison; it must not be read as a
biological claim about which receptors "have" which helices.

**v2 (superseded, commit `abd2257`)** replaced ordinal counting with a real
reference-based 3D structural mapping:
1. One reference LBD structure is selected per receptor
   (`data/analysis/lbd_structural_reference_structures.csv`, Part 2).
2. The 5 anchor receptors (VDR, PPARgamma, ESR1, GR, NR4A1) get an explicit
   H1-H12/beta-sheet **template** built from literature charge-clamp/
   dimerization landmarks (Nolte 1998, Shiau 1998, Bledsoe 2002, Rochel 2000,
   UniProt) cross-referenced against that receptor's own secondary-structure
   evidence -- never from counting
   (`data/analysis/canonical_lbd_template_segments.csv`, Part 3).
3. Every other receptor's reference structure is structurally aligned against
   all 5 templates with **Combinatorial Extension** (CE; Shindyalov & Bourne
   1998) -- already shipped inside Biopython (`Bio.PDB.ccealign`), so **no new
   software was installed**. CE returns actual residue-residue correspondence
   (not just a global score), which is used to transfer the template's
   segment labels, then the boundary is refined against the target's own
   secondary structure (`data/analysis/lbd_structural_alignment_scores.csv`,
   `lbd_structural_label_transfer.csv`, Parts 4-6).

v2's validated CE transfer machinery proved H3 mapped to 39/42 receptors, but
H5-H10 stayed at 0 -- **not because those helices are structurally absent**,
but because v2's 5 anchor templates only had literature landmarks for a
handful of helices each; a segment that has no anchored template definition
in *any* anchor can never be transferred, regardless of how well CE aligns.

**v3 (current, this section, scripts `44_curate_templates_v3.py` ->
`47_validate_v3.py`)** closes that template-construction gap without
touching normalization, aggregation, or the alignment algorithm itself:
1. **Complete reference-template curation** (Part 1-2): each of the 5
   anchors gets a PRIMARY structure plus up to 2 COMPLEMENTARY structures
   (`data/analysis/canonical_lbd_reference_structures_v2.csv`), and H1-H12 +
   beta_sheet boundaries are curated from published structural descriptions,
   the reference structures' own secondary-structure evidence, and known
   functional landmarks -- in that priority order, never from ordinal helix
   counting (`data/analysis/canonical_lbd_template_segments_v2.csv`). Where
   an anchor's own structures still leave a core segment (H5/H7/H10/H11/H12)
   un-anchored, **cross-template rescue** CE-aligns that anchor against the
   *other* anchors that do have it and transfers the boundary -- only when
   the transferred span overlaps a real secondary-structure run in the
   target anchor and doesn't collide with an already-anchored segment. This
   raised anchor-level core-segment coverage from ~3-4/12 to 9-11/12 per
   anchor (H8 remains unresolved in all 5 anchors -- no literature landmark
   or reference-structure evidence was found for it anywhere).
2. **Multi-structure receptor mapping** (Part 4): for non-anchor receptors,
   if the chosen primary structure doesn't resolve a core segment,
   *complementary* candidate structures for that same receptor are searched
   and CE-aligned; a successful rescue is recorded with
   `assignment_method = complementary_structure_rescue` and traced in
   `data/analysis/receptor_segment_mapping_sources.csv`. Structures are never
   averaged together for B-factor normalization -- rescue only ever
   transfers a segment *boundary*, never a flexibility value.
3. **Boundary refinement never fabricates a span** (Part 5): every transfer
   (main CE transfer, cross-template rescue, and complementary rescue) must
   find a real overlapping secondary-structure run of the expected class in
   the *target's own* structure. If no such run exists, the assignment is
   **rejected outright** -- earlier drafts of this stage twice fell back to
   a raw min/max span across the aligned core positions when no overlapping
   run was found, which silently bridged genuine evidence gaps (unobserved/
   zero-occupancy/coil positions) and mislabeled them HELIX/SHEET; all three
   occurrences of this bug class (script 41 in v2, the v3 cross-template
   rescue, and the v3 main transfer loop) were found and fixed the same way:
   reject, never fabricate.
4. **Canonical N->C order enforcement** (Part 7, new in v3): because
   segments for one receptor can now come from different templates/rescue
   paths, local per-alignment monotonicity is not enough -- a mandatory
   visual inspection of the five anchor schematics caught VDR's
   cross-template-rescued "H5" landing *before* VDR's own literature-anchored
   "H3". `enforce_canonical_order()` in `46_rebuild_segment_map_v3.py` now
   greedily resolves any cross-segment order inversion per receptor, keeping
   the higher-confidence/more-direct-evidence block and dropping (never
   silently discarding -- every drop is logged) the conflicting one; 868 of
   6157 named position-rows were dropped this way across ~30 segment
   instances in 16+ receptors.
5. **H12 audit** (Part 6): every coordinate-containing receptor without a
   mapped H12 is classified A (no structure observes it) / B (primary misses,
   a complementary structure observes it -- rescued automatically where
   possible) / C (alignment fails despite observed H12) / D (no anchor
   template curates the segment at all) / E (true receptor-specific
   exception) in `reports/tables/h12_mapping_audit.csv` (and more briefly for
   H5/H10/H11 in `h5_h10_h11_mapping_audit.csv`). All 14 v3 H12-missing
   receptors land in category C: each one has real helix evidence
   (reference or complementary) and a usable CE alignment (z >= 5.86), but
   the transfer's own boundary either found no overlapping target run or
   didn't clear the confidence threshold -- i.e. genuine per-receptor
   alignment/evidence shortfalls, not a missing-template (D) or
   absent-helix (A) problem.

### Acceptance-criteria result (Part 9) -- PASSED, reported not hidden

| Segment | v1 (48 = whole cohort) | v2 (42 w/ coords) | v3 (42 w/ coords) |
|---|---|---|---|
| H1 | 37 | 34 (81%) | 33 (79%) |
| H2 | -- | -- | 9 (21%) |
| H3 | **1** | **39 (93%)** | 34 (81%) |
| H4 | 1 | 12 (29%) | 12 (29%) |
| H5 | 0-1 | **0** | **28 (67%)** |
| H6 | 0-1 | **0** | **20 (48%)** |
| H7 | 0-1 | **0** | **25 (60%)** |
| H8 | 0-1 | **0** | **0 (0%)** |
| H9 | 0-1 | **0** | **11 (26%)** |
| H10 | 0-1 | **0** | **8 (19%)** |
| H11 | 0 | 9 (21%) | 21 (50%) |
| H12 | 39 | 26 (62%) | 28 (67%) |
| beta_sheet | 23 | 11 (26%) | 27 (64%) |

H3, H5, H7, H10/H11, H12 all have broad coverage where coordinates exist, and
none of H5/H7/H10 are near zero -- **acceptance criteria PASSED**. H8 is at
`0%`: no anchor structure, reference structure, or complementary structure
anywhere in the 42-receptor cohort was found with curatable secondary-
structure evidence for it. **This is a mapping-support statement, not a
biological claim** -- H8 was not confidently resolved by the current
structural-template mapping and is therefore retained as unresolved/not
available for the primary feature matrix, never fabricated to fill the gap.
Some v3 numbers (H1, H3, H4, H12) moved slightly *down* from a pre-fix
in-session draft of the v3 numbers, even though template curation only grew
-- this is the boundary-fabrication fix (item 3 above) correctly rejecting
assignments that a raw min/max fallback had over-reached on (most visibly for
beta_sheet, a short element where a fabricated span is especially likely to
bridge real gaps: the pre-fix draft over-reported 97.6%, corrected here to
64% once every fallback-span assignment without real target secondary-
structure evidence was rejected). Reported v3 numbers throughout this handoff
are the final, post-fix numbers. See
`reports/tables/structural_segment_validation_v3.csv` and
`reports/tables/segment_feature_missingness.md` for the full picture, and
`data/analysis/nr_lbd_structural_segment_map.csv`'s `evidence_notes` for every
receptor's specific reasoning. v1 and v2 are preserved in git history at
commits `9f15062` and `abd2257` respectively and are never rewritten.

### Files (`data/analysis/`, plus `reports/tables/` reports and QC figures)

| File | Contents |
|---|---|
| `receptor_dssp_consensus.parquet` | Per-receptor consensus secondary structure (HELIX/SHEET/COIL) at every observed canonical LBD position, from mmCIF `_struct_conf`/`_struct_sheet_range` (no DSSP binary is installed; this is the depositor-supplied secondary-structure evidence already in every cached mmCIF file) |
| `lbd_structural_reference_structures.csv` | v2: the one reference LBD structure selected per receptor (Part 2), with selection criteria recorded |
| `canonical_lbd_reference_structures_v2.csv` | **v3 current**: PRIMARY + up to 2 COMPLEMENTARY reference structures per anchor receptor (Part 1), with role/resolution/coverage/segments-supported/publication |
| `canonical_lbd_template_segments.csv` | v2: the 5 anchor receptors' literature-anchored H1-H12/beta-sheet templates (Part 3) |
| `canonical_lbd_template_segments_v2.csv` | **v3 current**: complete H1-H12/beta-sheet anchor templates, curated from published structural description > reference-structure secondary structure > functional landmarks > cross-template correspondence (Part 2), each row citing `evidence_type`/`publication`/`confidence` |
| `lbd_structural_alignment_scores.csv` | Every receptor-vs-template CE alignment: n_aligned, z_score, RMSD, coverage, sequence identity (v3 reruns this against the richer v3 templates; same schema as v2) |
| `lbd_structural_label_transfer.csv` | **v3 current**: per-target transferred + boundary-refined segment calls, now including `complementary_structure_rescue` rows (Part 5-6); boundary refinement rejects (does not write) any assignment with no overlapping target secondary-structure run -- never a fabricated min/max span |
| `receptor_segment_mapping_sources.csv` | **v3 new**: which source structure (primary vs. complementary) and method supplied each receptor's each segment (Part 4) |
| `reports/tables/h12_mapping_audit.csv` / `h5_h10_h11_mapping_audit.csv` | **v3 new**: A/B/C/D/E category audit for every coordinate receptor missing H12 (and, more briefly, H5/H10/H11) (Part 6) |
| `nr_lbd_structural_segment_map.csv` | **Current/live**: canonical UniProt position -> structural segment, with `reference_template`, `structural_alignment_score`, `assignment_confidence`/`assignment_status`/`evidence_notes` per position (v3 method, including `enforce_canonical_order()` cross-segment N->C conflict resolution; v1's version is preserved only in git history at commit `9f15062`, v2's at `abd2257`) |
| `pdb_segment_coverage.parquet` | receptor x PDB x segment: expected/observed residues, coverage (metadata) |
| `pdb_segment_flexibility_features.parquet` | receptor x PDB x segment: z-score summaries (from the already-approved PDB-median position profiles) + coverage + ligand state |
| `receptor_segment_flexibility_profiles.parquet` | receptor x segment: median/mean-of-entry-median z, APO/HOLO descriptive splits (no significance tests) |
| `receptor_segment_feature_matrix_median.parquet` / `_mean.parquet` | receptor x segment PRIMARY / sensitivity feature matrices (missing = NaN); v3 shape is 41 receptors x 39 segment columns |
| `pdb_segment_feature_matrix_median.parquet` / `_mean.parquet` | entry-level equivalents; v3 shape is 1789 PDB entries x 39 segment columns |
| `receptor_segment_coverage_matrix.parquet` | receptor x segment coverage (kept separate from flexibility matrices) |
| `nr_lbd_common_alignment_map.parquet` | SECONDARY per-residue common-position map, built only from the same conservative segment assignments (gaps explicit, never forced) |
| `reports/tables/structural_segment_validation_v2.csv` / `_v3.csv` | Anchor-receptor segment-by-segment validation: alignment support, secondary-structure support, functional-site consistency |
| `reports/tables/segment_feature_missingness.md` | Missingness policy audit (candidate imputation strategies reported, none chosen); v3: 0 segments >=90%, 0 >=80%, 1 >=70%, 7 >=50% of the 48-receptor cohort |
| `reports/figures/alignment_qc/*.png` | Coverage heatmap, PDB-support heatmap, missingness summary, 5 anchor-receptor schematics (v1 `04_example_...`, v2 `v2_schematic_...`, v3 `v3_schematic_...` -- all visually inspected for correct N->C order and functional-site-marker placement), segment-Z distributions |

## Segment-Level Modeling (PCA / UMAP / clustering, scripts `48`-`54`)

> Continues from commit `7d280bc` (v3 structural-segment mapping) with
> supervisor approval to proceed to unsupervised modeling. **This stage
> never modifies** `structures.csv`, normalization, or the structural-
> segment mapping -- the five upstream files it reads are hash-frozen in
> `data/analysis/modeling_frozen_input_hashes.json` and checked unchanged
> by `tests/test_modeling_pipeline.py`. Receptor group/subfamily labels are
> **never** used to fit PCA, UMAP, imputation, or any clustering algorithm
> -- only joined back in afterward for visualization.

An objective, non-hard-coded feature-support audit
(`reports/tables/modeling_feature_support.csv`) selected an 8-feature
primary set (7 core segments + `loop_H3_H5`, all with >=50% support among
the 42 coordinate receptors); a 60%-observed cohort rule
(`data/analysis/modeling_receptor_cohort.csv`) then included 29 of those
42 receptors (13 excluded with explicit reasons, including VDR itself --
its own template lacks H5/H7/H11). PCA (8 features -> 5 PCs for 80%
variance) and a full clustering battery (K-means, Ward, GMM, Spectral,
DBSCAN; validated with Silhouette/Davies-Bouldin/Calinski-Harabasz/Gap
Statistic/AIC-BIC/ARI) were run on the resulting matrix, with sensitivity
checks across imputation method, feature matrix, and PCA- vs feature-space
clustering. **Result: the evidence does not support a robust discrete
cluster structure** -- silhouette scores are modest, the single
best-silhouette solution is a degenerate 1-vs-28 outlier split, DBSCAN
never finds more than one density group, and cross-method agreement is
only moderate -- consistent with (not proof of) a continuum. Full details,
every number, and all limitations are in
`reports/tables/modeling_summary.md`; per-stage tables/figures are under
`reports/tables/modeling_*.csv`, `reports/tables/clustering_*.csv`,
`reports/tables/pca_variance.csv`, `reports/tables/umap_stability.csv`, and
`reports/figures/modeling/`. No hypothesis test (apo/holo, subfamily,
functional-site, disease) was run in this stage -- explicitly deferred.

## Subfamily/Group Testing & Flexibility Atlas (scripts `55`-`62`)

RECEPTOR is the inferential unit (`median_of_entry_median_z` from
`receptor_segment_flexibility_profiles.parquet`); every one of the 39
realized structural segments is tested independently on its own
available-case receptor set (a receptor missing one segment still
contributes to every other segment's test -- never restricted to the
29-receptor primary clustering cohort). PRIMARY grouping = the
professor-supplied project `group` field (largest, most balanced: 25/12/5);
SECONDARY = formal NR family-number; formal NR subfamily-letter reported
but unusable (every subfamily n<=4). Kruskal-Wallis omnibus per eligible
segment (>=3 receptors/group, >=2 eligible groups), epsilon-squared effect
size with bootstrap CI, BH-FDR within each grouping scheme.

**Result: 0/35 segment tests significant after FDR**, in either grouping
scheme, under either the median or mean receptor summary -- strengthens
the PCA/clustering continuum conclusion rather than contradicting it
(Outcome A). `data/analysis/superfamily_segment_flexibility_atlas.csv`
separately ranks every segment by median receptor Z (positive =
relatively flexible, negative = relatively rigid, **within that
receptor's own structures**, never an absolute cross-protein scale) --
descriptive, not a group-difference claim. Full detail:
`reports/tables/subfamily_flexibility_summary.md`.

## Apo vs Holo Flexibility (scripts `63`-`71`)

PDB entry is the analysis unit (`pdb_segment_flexibility_features.parquet`,
itself already chain -> per-instance Z -> entry-median). Eligibility is
receptor **and** segment specific (>=5 APO and >=5 HOLO PDB entries
PRIMARY, >=3/>=3 EXPLORATORY, both fixed before testing) --
`reports/tables/apo_holo_receptor_eligibility.csv`. **VDR is explicitly
INELIGIBLE: 0 APO, 48 HOLO PDB entries** -- never removed from any master
file, never plotted as though it had an APO comparison. Only 9/42
receptors are eligible at all.

Two-sided Mann-Whitney U + Cliff's delta per receptor x segment, BH-FDR
within receptor (PRIMARY) and globally (SECONDARY, never substituted).
**Result: 3/97 tests significant within-receptor FDR (ESR1 H1/H3, ERRgamma
loop_H6_H10), 0/97 significant globally** -- receptor-specific, localized
findings, not superfamily-wide, and all 3 are 100% leave-one-PDB-out
robust. Full detail: `reports/tables/apo_holo_summary.md`.

## Current Functional-Site Flexibility Map + Specificity Validation (scripts `72`-`92`)

A NEW analysis layer, entirely separate from the historical Role D
analysis below (never merged, never overwrites it). Four site classes:
AF2/H12 (from the current v3 H12 mapping), the ligand-binding pocket
(receptor-level consensus from real 5.0 Angstrom heavy-atom ligand-contact
geometry, reusing the already-validated FUNCTIONAL_LBD_LIGAND
classification), the coactivator groove, and the dimerization interface
(both via an evidence hierarchy: curated residues where available,
H3/H4/H5 or H9/H10/H11 segment correspondence otherwise, explicitly
labeled SEGMENT-level when so). Every receptor x site is explicitly
AVAILABLE/PARTIAL/UNRESOLVED -- never zero-filled.

Within-receptor paired Wilcoxon (site vs. whole-LBD background, PDB entry
as the paired unit) found 73/107 significant within-receptor FDR. A
follow-up specificity-validation stage then re-tested every RESIDUE-level
site against a **structurally matched** background (same segment, or same
structural class) instead of the whole LBD -- SEGMENT-level proxies are
never matched this way (would be circular). **Result: AF2/H12 and the
ligand-binding pocket are SITE_SPECIFIC_SUPPORTED** (confirmed flexible /
confirmed rigid against the matched background, 100% leave-one-receptor-
out stable); **the coactivator groove is STRUCTURAL_CONTEXT_EXPLAINED**
(its original apparent rigidity is better attributed to surrounding
helical context than to the groove residues specifically -- a genuine
correction to the face-value original interpretation); **the dimerization
interface has INSUFFICIENT_RESIDUE_LEVEL_MAPPING** to confirm or refute.
Full detail: `reports/tables/current_functional_site_summary.md`,
`reports/tables/functional_site_specificity_summary.md`.

## Historical strict QC subset — how `structures_primary_qc_subset.csv` was defined

A structure was in the strict subset (`previous_primary_qc_member == True`) only if **all** of:
- **X-ray crystallography**
- **>= 90% of the standardized LBD has an observed (occupancy > 0) Ca** —
  this is the Stage 3B.1-corrected "observed" definition; a zero-occupancy
  placeholder atom record does NOT count as observed.
- **Resolution between 1.8 A and 3.5 A inclusive**
- **R-free <= 0.30**
- **Human or human-derived identity** (confirmed by an explicit,
  individually-documented taxonomy audit, not just a taxonomy-ID filter —
  see `identity_class` in the data dictionary)

Within each (receptor x PDB entry x receptor polymer entity) group with
more than one candidate instance, exactly one **representative instance**
is chosen by a predeclared, fully deterministic tie-break (highest
observed coverage, then B-factor-usable coverage, then fewest internal
missing residues, then fewest zero-occupancy-only residues, then fewest
total missing residues, then lexicographic chain ID) — **never** by mean
B-factor, B-factor shape, ligand state, or any downstream result. Full
detail and every candidate's numbers are in
`data/manifests/final_representative_instance_selection.csv`.

**Ultrahigh-resolution structures (< 1.8 A) are never deleted** — the
project's historical/prior-cohort protocol specifies the 1.8 A lower
boundary for the primary set, but every ultrahigh structure that would
otherwise qualify is preserved, flagged, and fully explained in
`excluded_structures.csv` (reason code `RESOLUTION_LT_1_8`) and in
`final_structure_selection_audit.csv`'s `in_ultrahigh_resolution_sensitivity`
column, for robustness analysis.

Multi-UniProt fusion/chimera constructs (e.g. a receptor LBD fused to a
coactivator peptide, or expressed as an MBP fusion) are **not**
automatically excluded — every one of the 46 project fusion entities has a
sequence-region-disjoint fusion partner that never overlaps the
receptor's own mapped sequence (confirmed programmatically, not assumed),
so a fusion structure stays in the primary set as long as it passes every
other criterion. the master's `fusion_or_chimera` column flags these,
and `in_no_fusion_sensitivity` lets you exclude them for a robustness
check without re-deriving anything.

## Ligand state (APO / HOLO / AMBIGUOUS) — all coordinate candidates

- **HOLO**: at least one nonpolymer component is classified
  `FUNCTIONAL_LBD_LIGAND` (multi-evidence: RCSB's own "subject of
  investigation" annotation, deposited `_struct_site` binding-site
  records, geometric contact with the LBD within 4.5 A, and heavy-atom
  size) AND contacts the LBD.
- **APO**: no functional ligand, and no unresolved ambiguity.
- **AMBIGUOUS**: something contacts the LBD that could not be confidently
  classified functional vs. incidental.

**Water can never by itself cause a HOLO call. A simple monatomic ion
(Na+, Cl-, Zn2+, ...) can never by itself cause a HOLO call. A generic
crystallization additive (glycerol, PEG, sulfate, HEPES, ...) can never by
itself cause a HOLO call merely because it happens to be near the
protein.** See `scripts/utils/nonpolymer_classification.py` for the full,
documented, multi-evidence classification cascade — it is not a small
hand-picked exclusion list assumed complete.

**No pharmacology is inferred.** This dataset never labels a ligand
agonist, antagonist, partial agonist, or inverse agonist — that would
require an authoritative structured source this project does not have,
and guessing from chemical identity was explicitly out of scope.

## B-factor handoff — what's already done, what's NOT

Already done, in `lbd_residue_bfactors_all.parquet` (and its strict-subset counterpart):
- Every standardized LBD position, for every primary structure, whether
  observed or not (`is_construct_mapped`, `ca_positive_occupancy`, etc. —
  never silently dropped).
- Raw `B_iso_or_equiv` values for every deposited Ca altloc
  (`b_iso_values_json`, parallel to `altloc_ids_json`/`occupancies_json`).
- ONE recommended Ca record per residue (`selected_ca_*` columns), chosen
  purely by occupancy and altloc identity — **B-factor magnitude was never
  used to pick a structure, chain, altloc, or representative instance,
  anywhere in this pipeline.**

**NOT done, and intentionally left for the downstream analysis:**
- B-factor normalization (e.g. per-structure Z-scoring) — B-factor scale
  is refinement-program- and resolution-dependent and is not comparable
  raw across structures.
- Any clustering, dimensionality reduction, or statistical testing.
- Any ligand-state effect analysis.

## Important warnings — read before you build on this

- The standardized LBD interval is a **UniProt/PROSITE-based definition**
  ("Option A", PS51843-backed), not the previous cohort's ad hoc VDR
  interval — see `data/manifests/nr_lbd_reference.csv` and
  `reports/tables/stage3a_lbd_boundary_audit.md` for the reconciliation
  against the prior cohort's repository.
- **Sequence-mapping coverage, observed-coordinate coverage, and
  B-factor-usable coverage are three different numbers** — never conflate
  `construct_lbd_mapping_coverage`, `positive_occupancy_lbd_coverage`, and
  `bfactor_usable_lbd_coverage`.
- Zero-occupancy atom records are deposited placeholders, **not**
  observations — this project's "observed" always means
  occupancy > 0.
- Fusion constructs are flagged, not hidden — check `fusion_or_chimera`
  before assuming a chain is "just" the receptor.
- Sensitivity flags (`coverage_ge_0XX`, `in_completeness_*_sensitivity`,
  `in_no_fusion_sensitivity`, `in_ultrahigh_resolution_sensitivity`) exist specifically so you can test how
  sensitive your downstream conclusions are to this project's threshold
  choices — use them before assuming the primary set is the only valid
  cut.
- Missing residues (unmapped construct positions, or mapped-but-unobserved
  positions) are never silently removed from `lbd_residue_bfactors_all.parquet`
  — every standardized LBD position has a row for every primary structure,
  with explicit boolean flags for why it may lack usable coordinates.
- The master inventory has rows for 45 of the 48 receptors; RORβ, TR2 and
  GCNF have no rows because RCSB returned zero experimental candidates
  (`final_counts_by_receptor.csv` still lists all 48). 42 receptors have
  Stage 3B coordinate data; COUP-TFI, NOR-1 and DAX-1 have inventory rows
  only (no X-ray + LBD-overlap coordinate data). Missing values are
  explicit (blank) with `coordinate_data_status`, never silently defaulted.
- `ligand_state` is assessed for all 1840 coordinate candidates using the
  same conservative rules. Contact geometry uses only *observed*
  (occupancy > 0) LBD atoms of the representative instance, so in
  low-coverage structures an unresolved pocket could hide a ligand: treat APO
  calls at low `positive_occupancy_lbd_coverage` with care (filter on coverage
  yourself). One candidate (1YNW_4, RXRα, a DBD-only entity with no observed
  LBD atoms) is AMBIGUOUS / REVIEW_NEEDED because contact cannot be measured.
  Polymer-partner/coregulator annotations (`final_polymer_partner_annotations.csv`)
  still cover the strict subset only.

## Regenerating the handoff tables locally

The processed tables are tracked in Git (small enough — see `.gitignore`'s explicit exceptions
under `data/processed/`), but if you ever need to regenerate it from
scratch (all inputs are already-cached, offline, no network access
required once `data/raw/mmcif/` and `data/raw/api/` are populated):

```bash
conda activate biol363_nr
python scripts/24_final_structure_selection.py
python scripts/25_fetch_nonpolymer_metadata.py
python scripts/26_nonpolymer_inventory_and_ligand_state.py
python scripts/27_polymer_partner_annotation.py
python scripts/28_build_final_residue_map.py
python scripts/29_build_structures_csv.py   # historical strict subset
python scripts/31_build_master_inventory.py # MASTER structures.csv + lbd_residue_bfactors_all.parquet
```

Each script is idempotent and was verified to reproduce byte-identical
output (aside from timestamp fields) across repeated reruns from the
cached inputs.
