# Canonical remap stage — method note

**Status (2026-10-21):** final. This stage defines the two corrections that the corrected
v2 analysis is built on. Both are applied in the v2 rebuild described in
[`reports/rebuild_v2/README.md`](../rebuild_v2/README.md), which also reports the v2
results. The frozen v1 release (`main` @ 18b8fc9) is not modified by anything here.

1. **Professor canonical helix map = final ground truth** (decision of 2026-10-20).
2. **Corrected PDB → UniProt residue numbering**, rebuilt from sequence evidence.

## 1. Where everything lives

| Content | Location |
|---|---|
| Professor sources, byte-for-byte (+ SHA-256 provenance) | `data/sources/professor_canonical_helix_table/` |
| Stage data | `data/canonical_remap/` |
| Stage audits and reports | `reports/canonical_remap/` |
| Stage scripts | `scripts/96_…`–`100_…`, `scripts/utils/{sequence_numbering,canonical_helix}.py` |
| Tests | `tests/test_residue_numbering_crosswalk.py`, `tests/test_canonical_helix_map.py` |

The historical identity audit (`scripts/93`) scans `data/analysis`, `data/manifests`,
`data/processed` and `reports/tables`, and its frozen report
`reports/tables/final_receptor_identity_audit.csv` must stay byte-identical. This
stage (and v2) therefore writes **nothing** into those folders. Its own identity audit,
`scripts/99_audit_remap_stage.py`, applies the unmodified checks of `scripts/93` to the
stage and v2 folders: all 48 receptors CLEAN
(`stage_receptor_identity_audit_v1.csv`).

## 2. Professor canonical helix map (final policy)

Source: `NR_canonical_helix_boundaries.csv` (434 intervals, 42 receptors; 27
`published_boundaries`, 354 `structural_homology`, 53 `review_required`).

- The table is the **ground truth**. Every interval is used exactly as supplied:
  boundaries are not modified, labels are not renamed, and nothing is filled or mixed
  from the superseded v3 map. v3 is kept only as audit history
  (`canonical_helix_v3_comparison_*`).
- **All** labels are used, including `review_required`. The professor status is kept
  verbatim on every row (it is never relabelled); excluding `review_required` is done
  only as a labelled sensitivity analysis.
- A position inside two professor intervals belongs to **both** elements (PXR 261: H3
  and H4; CAR 344: H11′ and H12). Joined strings such as `H3|H4` are display labels
  for one-row lookups only, never an element (`tests/test_v2_dual_label_positions.py`).
- An element may consist of more than one interval (RORγ H10/H11: 459–466 and
  470–490). The positions between its intervals are not part of the element.
- Positions in the analysis span outside every professor interval are
  `INTER_HELIX_UNASSIGNED_DERIVED`: derived, never called loops.
- Six receptors have no professor map (RORβ, TR2, COUP-TFI, GCNF, DAX-1, SHP):
  `CANONICAL_MAPPING_UNAVAILABLE`, not biological absence.
- **TRα** is numbered on **P10827-2** (TRα1), the professor's accession; every other
  receptor uses its canonical UniProt entry. P10827 and P10827-2 are identical over
  1–370 only.

Outputs: `data/canonical_remap/nr_lbd_canonical_helix_intervals_professor_2026.csv`,
`nr_lbd_canonical_helix_map_professor_2026.csv` (one row per position and label;
11,653 positions; summary in `canonical_helix_map_summary.json`). All 434 intervals lie
inside the final analysis span (§4); none is truncated.

## 3. Corrected residue numbering

**Defect (found 2026-10-17).** The frozen residue tables took entity → UniProt regions
from the RCSB sequence-coordinates service without checking residue identity. For
constructs with internal deletions the service could return one contiguous region (e.g.
VDR 1DB1: 118–376 ↔ entity 1–259, true construct 118–164 + 216–427), so every residue
after the deletion was misplaced (VDR +51 after residue 165).

**Ground truth** (`scripts/utils/sequence_numbering.py`), in order of trust:

1. the chain's polymer sequence from mmCIF `_pdbx_poly_seq_scheme`;
2. cross-checked against the residues actually observed (0 disagreements);
3. local alignment (BLOSUM62, gap −10/−0.5) against the reference UniProt sequence.

Author numbering, `_struct_ref_seq(_dif)` and the old RCSB regions are corroboration
only; B-factors are never used. A residue gets a position only if every co-optimal
alignment places it identically, the alignment is locally supported (identical
residue: ≥0.8 identity over 10 pairs on one side or ≥0.6 on both; substitution: ≥0.6 on
both sides), the depositor does not declare it a tag/linker, and the chain has enough
evidence (≥20 identical residues at ≥50% identity, or an exact unique match of ≥8
residues for short peptides). Ambiguous residues are withheld, never forced.

Results (`residue_numbering_summary_v1.json`, `residue_numbering_*_v1.csv`):
3,159 frozen instances (3,043 CLEAN, 53 deletion/insertion offset, 4 isoform offset,
2 register offset, 57 with ≤4 withheld junction residues); 10,599 of the 12,975
previously non-identical observed residues were numbering errors. Receptors with
shifted residues include VDR (45/48 entries), Rev-erbα, CAR, LXRβ, PXR, ERα, RXRα, ERβ
and FXR. In v2 every observed residue sits at its sequence-evidence position with 0
numbering disagreements (VDR identity 99.96%).
`entity_uniprot_regions_corrected_v1.csv` re-expresses this per polymer entity in the
legacy region schema used by the rebuild.

## 4. Final analysis span

For each **mapped** receptor the analysis span is the smallest continuous span that
covers the UniProt/PROSITE "NR LBD" domain (PS51843) **and** every professor interval.
Unmapped receptors keep the PROSITE span. 23 receptors are extended (e.g. PPARα and
PPARδ by 37 N-terminal residues for H1/H2; GR, MR, PR by 9 C-terminal residues for
H12); the rest are unchanged. File: `data/rebuild_v2/inputs/analysis_span_v2.csv`
(rule `PROSITE_NR_LBD_UNION_PROFESSOR_HELICES`). The z-score normalization of each
instance uses this span.

## 5. What happened next

- v2 foundation rebuild and the corrected rerun of all earlier analyses (Phase 2):
  [`reports/rebuild_v2/README.md`](../rebuild_v2/README.md).

## 6. History (superseded; kept for audit)

- `lbd_window_*_v1.csv` (`scripts/100`): the pre-build audit showing that the PROSITE
  window left 27 of the 434 professor intervals truncated or outside it; resolved by the final span rule (§4).
- `canonical_helix_endpoint_concordance*.csv` (`scripts/98`): professor endpoints vs
  deposited helix records (748/753 endpoints match under independent alignment).
- `canonical_helix_v3_comparison_*`, `numbering_impact_by_*`, `vdr_*`: the v1/v3
  impact assessment that motivated the rebuild.
- An earlier draft policy that withheld `review_required` labels and joined TRα
  through P10827 was replaced by the final policy in §2.

Open notes: two 1-residue deletions are not declared by their depositors (PXR 6HJ2 at
192, FXR 5Y44 at 471), and LXRβ 1UPV/1UPW depositor numbering appears to follow an
older UniProt version; sequence evidence is used in all three cases.
