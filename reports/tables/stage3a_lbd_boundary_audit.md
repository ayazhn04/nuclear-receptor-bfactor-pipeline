# Stage 3A.1 — LBD Boundary Definition Reconciliation

## A. Source definitions

**UniProt/PROSITE profile interval (`profile_lbd_*`)**
Generated from the UniProt feature `type: "Domain"`, `description: "NR LBD"`,
evidenced by `PROSITE-ProRule:PRUxxxxx`. Live queries to the official EBI
InterPro API confirm this ProRule wraps the PROSITE **profile** `PS51843`
("Nuclear receptor (NR) ligand-binding (LBD) domain profile", short name
`NR_LBD`, `type: domain`), which is integrated into InterPro entry
**IPR000536** ("Nuclear hormone receptor, ligand-binding domain") alongside
two independent domain-family resources that agree on the same boundary
family: **Pfam PF00104** and **SMART SM00430**. This is a computational
**domain-profile match**, not a short sequence motif — three independently
curated domain databases converge on it. InterPro's own description of the
domain explicitly states it "also contains a ligand-dependent AF-2 region,"
i.e. the profile is intended to include AF-2/H12, not stop short of it.
We did not independently verify that UniProt's ProRule numeric output is
byte-identical to a raw PS51843 profile scan (that would require
querying UniProt's internal ProRule database directly), so this is
reported as strong, but not absolutely proven, correspondence.

**Previous-cohort interval (`previous_cohort_lbd_*`)**
A **hard-coded script default** (`--start 118 --end 427`), found in 6
VDR-only analysis scripts in the previous cohort's repository. It is not
stored in any data file, is not explained or justified anywhere in that
repository (no README, no comments beyond "VDR LBD"), and — critically —
**427 is VDR's full canonical sequence length**: the window's C-terminal
boundary is simply "the end of the protein," not a specific structural or
domain annotation. It most likely reflects the common structural-biology
convention of using a "near-full-length" VDR LBD construct (since VDR's
LBD is unusually close to its own C-terminus already, with almost no
tail past H12), rather than an independently derived domain boundary.

**No other authoritative interval was retrieved** — PROSITE/InterPro
(covered above) is the only additional official source consulted, per
Section 8; no literature-curated boundary was manually typed for any
receptor.

## B. Previous cohort data availability

- **Source repository:** `https://github.com/sumeshi3648/BioInformatics-Bfactors.git`
- **Commit inspected:** `69899f60d7ab0eed19d3029687fe31f03b835ac2` (HEAD, a merge commit)
- **Files inspected:** `proteins.csv`, `data/meta/structures.csv`, `data/processed/ca_bfactors.csv`, all 14 files under `scripts/`
- **Receptors with an explicit LBD range: 1 of 48 (VDR only)**

`proteins.csv` (48 rows, matching this project's 48 receptors by UniProt
accession) contains only `subfamily, symbol, uniprot` — **no coordinate
columns at all**. No metadata sheet with "UniProt IDs, chains, and LBD
ranges" (as described in the project background) actually exists in this
repository. This is reported explicitly rather than assumed or worked
around. Relevant source files and their SHA256 checksums are snapshotted
under `data/raw/external/previous_cohort_repo/` (git-ignored, per project
convention for raw external data) with full provenance in
`data/manifests/nr_lbd_boundary_audit_provenance.json`.

## C. Discrepancy summary

- Exact matches: **0**
- Minor differences (≤2 residues at both termini): **0**
- Material differences: **1** (VDR)
- Previous range unavailable: **47**

Given only one receptor has both definitions, no aggregate "narrower/wider
at N-/C-terminus" pattern can be established across the family — this
question cannot be answered systematically for the other 47 receptors
from the available data (see Section G).

## D. Full 48-receptor boundary table

| common_name | uniprot_id | profile interval | previous-project interval | start Δ | end Δ | status |
|---|---|---|---|---|---|---|
| AR | P10275 | 669-900 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| ERα | P03372 | 311-547 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| ERβ | Q92731 | 264-498 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| GR | P04150 | 524-758 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| MR | P08235 | 726-964 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| PR | P06401 | 679-913 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RARα | P10276 | 183-417 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RARβ | P10826 | 183-417 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RARγ | P13631 | 185-419 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| TRα | P10827 | 163-407 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| TRβ | P10828 | 217-461 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| **VDR** | **P11473** | **127-423** | **118-427** | **-9** | **+4** | **MATERIAL_DIFFERENCE** |
| CAR | Q14994 | 109-352 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| ERRα | P11474 | 193-421 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| ERRβ | O95718 | 208-432 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| ERRγ | P62508 | 233-457 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| FXR | Q96RI1 | 262-486 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| HNF4α | P41235 | 147-377 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| HNF4γ | Q14541 | 99-328 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| LRH-1 | O00482 | 300-539 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| LXRα | Q13133 | 209-447 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| LXRβ | P55055 | 222-460 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| NUR77 | P22736 | 360-595 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| NURR1 | P43354 | 360-595 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| PPARα | Q07869 | 239-466 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| PPARδ | Q03181 | 211-439 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| PPARγ | P37231 | 238-503 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| PXR | O75469 | 146-433 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| Rev-erbα | P20393 | 284-614 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| Rev-erbβ | Q14995 | 369-579 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RORα | P35398 | 272-510 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RORβ | Q92753 | 222-460 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RORγ | P51449 | 269-508 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RXRα | P19793 | 227-458 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RXRβ | P28702 | 296-529 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| RXRγ | P48443 | 231-459 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| SF-1 | Q13285 | 222-459 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| TR2 | P13056 | 348-590 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| TR4 | P49116 | 341-583 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| COUP-TFI | P10589 | 184-410 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| COUP-TFII | P24468 | 177-403 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| DAX-1 | P51843 | 205-469 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| EAR-2 | P10588 | 165-393 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| GCNF | Q15406 | 249-480 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| NOR-1 | Q92570 | 394-623 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| PNR | Q9Y5X4 | 169-410 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| SHP | Q15466 | 16-257 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |
| TLX | Q9Y466 | 155-383 | N/A | — | — | PREVIOUS_RANGE_UNAVAILABLE |

## E. VDR deep check

| Quantity | Value |
|---|---|
| Full canonical UniProt sequence length | **427** |
| Profile (PROSITE PS51843 / ProRule) interval | **127-423** (length 297) |
| Previous-cohort script-default interval | **118-427** (length 310) |
| Last annotated Helix feature (UniProt) | **416-422** (this is H12/AF-2) |
| Motif feature "9aaTAD" (AF-2-associated transactivation motif) | **416-424** |
| Residues 425-427 | Unannotated C-terminal tail (`E-I-S`), beyond any Helix or Motif feature |

**Residue-by-residue answer for 424/425/426/427** (VDR has no residue 428 —
the canonical sequence is only 427 residues long):
- All four (424-427) belong to the **canonical sequence** (1-427).
- **424** is the last residue of the 9aaTAD/AF-2-associated motif (416-424); it is 1 residue past the last annotated Helix (416-422).
- **425-427** are **not part of any annotated Helix or Motif** — they are the literal unstructured C-terminal tail of the full-length protein.
- All four are within the **previous-cohort window** (118-427, since it runs to the literal C-terminus).
- **None** of the four are within the **profile interval** (127-423) — the profile interval ends exactly 1 residue before the 9aaTAD motif's C-terminal edge (424) and 4 residues before the true C-terminus.

**Conclusion:** the previous cohort's window is not derived from a specific
domain/motif boundary — it is "everything from a chosen N-terminal cutoff
(118) to the literal end of the protein (427)." The profile interval, by
contrast, is a computed domain-family match that happens to end 1 residue
short of the annotated 9aaTAD/AF-2 motif's C-terminal edge and captures the
entire annotated H12 helix (416-422). At the N-terminus, the previous
cohort's window is 9 residues wider (118 vs 127) than the profile
interval, with no comparable annotation-based justification found for
either specific number.

## F. Sensitivity impact (VDR only — the only receptor with both definitions)

Recomputed LBD-mapping coverage for all 52 VDR Stage 2 candidates under
both boundary definitions (`rcsb_lbd_boundary_sensitivity.csv`), without
touching the existing Stage 3A result:

- **Category changed: 4 / 52** — all four are the **DBD-only** VDR
  structures (1KB2_3, 1KB4_3, 1KB6_3, 1YNW_3; VDR DBD bound to DNA response
  elements). Under the profile interval they correctly show
  `NO_LBD_OVERLAP` (0% coverage); under the wider previous-cohort window
  they show a spurious `TRACE_LBD` (2.6% coverage) — a **false-positive
  LBD signal** created purely because the previous window's N-terminal
  boundary (118) reaches back far enough to graze the tail of the DBD
  construct (which extends to ~residue 125), even though these structures
  contain no LBD sequence at all.
- **Metadata-prefilter consequence changed: 4 / 52** — the same four
  records: under the profile definition they correctly fail the prefilter
  (`NO_LBD_OVERLAP`); under the previous-cohort definition they would
  spuriously pass the LBD-overlap component of the prefilter.
- **Receptor most affected:** VDR is the only receptor with a boundary
  definition to compare, so it is trivially "most affected"; within VDR,
  the effect is concentrated entirely in DBD-only constructs — no
  genuine LBD-containing VDR structure changes category under either
  definition.

**This is not a negligible difference.** The boundary choice determines
whether 4 real structures are correctly recognized as containing zero LBD
sequence, or incorrectly flagged as marginal LBD candidates.

## G. Recommendation options (no default chosen)

**Option A — Use the UniProt/PROSITE profile interval as the standardized superfamily analysis interval.**
Consequence: reproducible, uniform, and independently defined across all
48 receptors from a single authoritative, versioned source (auditable via
UniProt entry/ProRule version, cross-validated by Pfam/SMART through
InterPro). Cost: not identical to the previous cohort's VDR-only window,
so VDR results will not be a literal reproduction of last year's exact
numbers (they will differ in which residues are counted, particularly for
constructs whose boundary sits between 118-127 or 423-427).

**Option B — Use the previous-cohort biological LBD window for continuity with the 2025 pipeline.**
Consequence: only defined for VDR; no systematic, reproducible source
exists for the other 47 receptors, so this option cannot be applied
project-wide without either (a) restricting the superfamily analysis to
VDR only, or (b) inventing 47 more windows by some other means (explicitly
prohibited). Also, per Section F, the previous window measurably
introduces false-positive LBD signal on DBD-only VDR structures.

**Option C — Maintain both, distinctly labeled.**
Keep `profile_lbd_*` (all 48 receptors, the standardized/reproducible
superfamily-wide interval) as the primary analysis definition, and retain
`previous_cohort_lbd_*` (VDR only) purely as a documented historical/
continuity cross-check, never conflated with the primary definition.
Consequence: most transparent option; requires maintaining the
distinction in terminology (`profile_lbd_start`/`NEAR_COMPLETE_PROFILE_LBD`,
etc., per Section 9 of the instructions) throughout Stage 3B and beyond,
and accepting that VDR-specific continuity comparisons are the only ones
possible with current data.

No option has been silently selected. `nr_lbd_reference.csv` (the Stage 3A
profile-derived table) remains unchanged and is not being treated as
final pending your decision.
