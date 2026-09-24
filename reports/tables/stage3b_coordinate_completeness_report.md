# Stage 3B — mmCIF Coordinate Acquisition and Observed LBD Cα Completeness Audit (corrected at Stage 3B.1)

## Coverage definitions (Stage 3B.1 correction)

- **Atom-record-present coverage**: at least one Cα `_atom_site` record exists for the mapped position — includes zero-occupancy placeholder atoms. This is the historical `ca_modeled` concept; it is NOT equivalent to observation.
- **Positive-occupancy observed coverage**: at least one Cα record has `occupancy > 0`. **From Stage 3B.1 onward, this — not mere record presence — is what "observed Cα coverage" means throughout this project.**
- **B-factor-usable coverage**: at least one Cα record has `occupancy > 0` AND a finite numeric `B_iso_or_equiv`.

For a residue with multiple altloc Cα records, presence/positive-occupancy/usability are each true if ANY altloc satisfies the condition — occupancies are never summed across altlocs, and a residue is never counted more than once.

## Occupancy-correction impact

- Zero-occupancy-only LBD residues found: 281
- Instance/models with at least one zero-occupancy-only residue: 36
- Instance/models whose coverage bin or exact-complete status changed: 36
- Full detail: `data/manifests/stage3b_occupancy_correction_audit.csv`

## Three coverage distributions (full-LBD, all measured instance-models)


**Atom-record-present:**
- [0.95,0.99): 1212
- [0.90,0.95): 907
- 1.000: 360
- [0.80,0.90): 318
- [0.99,1.00): 292
- <0.80: 70

**Positive-occupancy observed:**
- [0.95,0.99): 1233
- [0.90,0.95): 909
- 1.000: 341
- [0.80,0.90): 318
- [0.99,1.00): 288
- <0.80: 70

**B-factor-usable:**
- [0.95,0.99): 1233
- [0.90,0.95): 909
- 1.000: 341
- [0.80,0.90): 318
- [0.99,1.00): 288
- <0.80: 70

## Deposited-unobserved-residue crosscheck (4 categories)


- A_DEPOSITED_UNOBSERVED_NO_CA_RECORD: 32756
- B_DEPOSITED_UNOBSERVED_ZERO_OCCUPANCY_CA: 270
- D_NOT_DEPOSITED_UNOBSERVED_NO_CA_RECORD: 13

CROSSCHECK_REVIEW_NEEDED count: 0

## Threshold sensitivity preview (positive-occupancy coverage; no threshold chosen)


| threshold | instances_retained | candidates_retained | unique_pdb_entries_retained | receptors_represented | receptors_lost_count |
|---|---|---|---|---|---|
| 1.0 | 341 | 236 | 233 | 15 | 27 |
| 0.99 | 629 | 464 | 459 | 22 | 20 |
| 0.95 | 1862 | 1262 | 1243 | 29 | 13 |
| 0.9 | 2771 | 1680 | 1636 | 34 | 8 |
| 0.8 | 3089 | 1811 | 1755 | 39 | 3 |

## Design note for future representative-chain selection (Stage 3C)

Representative-chain selection must operate at the level of **project receptor × PDB entry × receptor polymer entity** — NOT globally one chain per PDB entry. A single PDB entry can contain separate polymer entities corresponding to different project receptors (e.g. an RAR/RXR heterodimer entry); selecting one chain per entry would incorrectly discard one receptor's data. Within repeated instances of the SAME receptor polymer entity in one entry, Stage 3C may choose one deterministic representative instance — that choice is deferred, not made here.

## Multi-model observation

No multi-model structures were observed among the 1777 selected X-ray entries in this Stage 3B coordinate-audit pool.

## Important statements


No final LBD-coordinate completeness threshold was selected at Stage 3B or Stage 3B.1.

No representative chain was selected at Stage 3B or Stage 3B.1.

No B-factor normalization or downstream statistical analysis was performed.
