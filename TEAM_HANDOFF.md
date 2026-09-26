# Team Handoff — Nuclear Receptor Structural Dataset


> **QC criteria are metadata/analysis flags. Structures are retained in the
> master inventory rather than deleted based on coverage thresholds.**
> Coverage, resolution, R-free, method, taxonomy, fusion status etc. are
> columns you filter on yourself — nothing was removed from the master
> dataset because of them (revised per professor's requirement).

**Data & Infrastructure responsibility is complete as of this document.
No B-factor normalization, clustering, PCA/t-SNE/UMAP, statistical
testing, or biological-conclusion analysis has been performed anywhere in
this repository.** That work starts from the files below.

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
