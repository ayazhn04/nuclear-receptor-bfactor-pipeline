# Team Handoff — Nuclear Receptor Structural Dataset

This document is for teammates who did **not** build this pipeline and
need to start the downstream B-factor/flexibility analysis. It tells you
exactly which file to open, what the numbers mean, and what has and has
NOT already been done.

**Data & Infrastructure responsibility is complete as of this document.
No B-factor normalization, clustering, PCA/t-SNE/UMAP, statistical
testing, or biological-conclusion analysis has been performed anywhere in
this repository.** That work starts from the files below.

## What to use

| Need | File |
|---|---|
| Which structures are in the dataset, their QC, ligand state, resolution, etc. | [`data/processed/structures.csv`](data/processed/structures.csv) |
| Why a specific PDB/receptor is NOT in the dataset | [`data/processed/excluded_structures.csv`](data/processed/excluded_structures.csv) |
| Per-residue canonical numbering + raw Ca coordinatable data for B-factor work | [`data/processed/final_lbd_residue_map.parquet`](data/processed/final_lbd_residue_map.parquet) |
| Ligand identity and evidence behind each APO/HOLO/AMBIGUOUS call | [`data/manifests/final_ligand_annotations.csv`](data/manifests/final_ligand_annotations.csv) |
| Every nonpolymer component (ligand, ion, water, buffer...) in each structure | [`data/manifests/final_nonpolymer_inventory.csv`](data/manifests/final_nonpolymer_inventory.csv) |
| Coactivator/corepressor peptides or other protein partners bound near the LBD | [`data/manifests/final_polymer_partner_annotations.csv`](data/manifests/final_polymer_partner_annotations.csv) |
| The standardized LBD interval definition used throughout | [`data/manifests/nr_lbd_reference.csv`](data/manifests/nr_lbd_reference.csv) |
| Column-by-column meaning of every field above | [`reports/tables/final_data_dictionary.md`](reports/tables/final_data_dictionary.md) |
| Full receptor-by-receptor breakdown (primary/apo/holo/sensitivity counts, all 48 receptors) | [`data/manifests/final_counts_by_receptor.csv`](data/manifests/final_counts_by_receptor.csv) |
| One-shot machine-readable summary of the whole release | [`data/manifests/final_data_infrastructure_release.json`](data/manifests/final_data_infrastructure_release.json) |

**Start with `structures.csv` and `final_lbd_residue_map.parquet`.** Those
two tables are what a typical B-factor/flexibility analysis needs.

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
in `final_lbd_residue_map.parquet` purely for cross-referencing against the
original PDB file — never for joining across structures. Every table in
this handoff is keyed, or can be joined, on `canonical_uniprot_position`.

## Primary structure selection — what's in `structures.csv` and why

A structure is in the primary set only if **all** of:
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
other criterion. `structures.csv`'s `fusion_or_chimera` column flags these,
and `passes_no_fusion_sensitivity` lets you exclude them for a robustness
check without re-deriving anything.

## Ligand state (APO / HOLO / AMBIGUOUS)

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

Already done, in `final_lbd_residue_map.parquet`:
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
- Sensitivity flags (`passes_95pct_completeness`, `passes_99pct_completeness`,
  `passes_no_fusion_sensitivity`, and the ultrahigh-resolution records in
  `excluded_structures.csv`) exist specifically so you can test how
  sensitive your downstream conclusions are to this project's threshold
  choices — use them before assuming the primary set is the only valid
  cut.
- Missing residues (unmapped construct positions, or mapped-but-unobserved
  positions) are never silently removed from `final_lbd_residue_map.parquet`
  — every standardized LBD position has a row for every primary structure,
  with explicit boolean flags for why it may lack usable coordinates.
- 33 of the 48 project receptors have >= 1 primary structure; 15 do not
  (either zero RCSB search results, no coordinate candidate ever entered
  the Stage 3B pool, or none of their candidates simultaneously satisfied
  every primary criterion — `final_counts_by_receptor.csv` tells you
  exactly which, for every one of the 48).

## Regenerating `final_lbd_residue_map.parquet` locally

It is tracked in Git (small enough — see `.gitignore`'s explicit exceptions
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
python scripts/29_build_structures_csv.py
```

Each script is idempotent and was verified to reproduce byte-identical
output (aside from timestamp fields) across repeated reruns from the
cached inputs.
