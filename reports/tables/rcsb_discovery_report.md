# RCSB Candidate Structure Discovery — Stage 2 (corrected at Stage 2.1)

## Methods

For each of the 48 UniProt-validated human nuclear receptors (`config/nr_metadata.py`), the official RCSB Search API (`https://search.rcsb.org/rcsbsearch/v2/query`) was queried for current, EXPERIMENTAL polymer entities whose `reference_sequence_identifiers` mapping includes the receptor's validated UniProt accession. Computed structure models were excluded via `rcsb_entry_info.structure_determination_methodology == "experimental"`; no restriction to a specific experimental method (X-ray/cryo-EM/NMR) was applied at this stage. Discovered polymer entities and their parent PDB entries were enriched via the official RCSB Data API (GraphQL, batched). GraphQL introspection against `RcsbPolymerEntityAlign` confirmed the Data API exposes no precomputed coverage-fraction field — only `aligned_regions` (entity_beg_seq_id, ref_beg_seq_id, length). `reference_sequence_coverage` and `entity_sequence_coverage` are therefore computed deterministically from these raw region lengths and the already-retrieved entity/reference sequence lengths (Stage 1 UniProt data for the project's 48 accessions); reference coverage is left null for non-project (e.g. fusion partner) accessions rather than guessed. No mmCIF or PDB coordinate files were downloaded at any point.

## Overall counts

- number_of_validated_receptors_queried: 48
- number_of_successful_receptor_queries: 48
- number_of_failed_receptor_queries: 0
- total_receptor_polymer_entity_relationships: 2072
- unique_polymer_entities: 2072
- unique_pdb_entries: 1996
- receptors_with_zero_candidates: 3
- receptors_with_1_to_2_candidates: 10
- receptors_with_3_or_more_candidates: 35

## All 48 receptors

| common_name | nr_code | uniprot_id | polymer_entity_count | unique_pdb_entry_count | query_status |
|---|---|---|---|---|---|
| AR | NR3C4 | P10275 | 102 | 95 | SUCCEEDED |
| ERα | NR3A1 | P03372 | 517 | 492 | SUCCEEDED |
| ERβ | NR3A2 | Q92731 | 39 | 39 | SUCCEEDED |
| GR | NR3C1 | P04150 | 57 | 57 | SUCCEEDED |
| MR | NR3C2 | P08235 | 30 | 30 | SUCCEEDED |
| PR | NR3C3 | P06401 | 21 | 21 | SUCCEEDED |
| RARα | NR1B1 | P10276 | 14 | 14 | SUCCEEDED |
| RARβ | NR1B2 | P10826 | 9 | 9 | SUCCEEDED |
| RARγ | NR1B3 | P13631 | 11 | 11 | SUCCEEDED |
| TRα | NR1A1 | P10827 | 10 | 10 | SUCCEEDED |
| TRβ | NR1A2 | P10828 | 34 | 34 | SUCCEEDED |
| VDR | NR1I1 | P11473 | 52 | 52 | SUCCEEDED |
| CAR | NR1I3 | Q14994 | 2 | 2 | SUCCEEDED |
| ERRα | NR3B1 | P11474 | 6 | 6 | SUCCEEDED |
| ERRβ | NR3B2 | O95718 | 3 | 3 | SUCCEEDED |
| ERRγ | NR3B3 | P62508 | 43 | 43 | SUCCEEDED |
| FXR | NR1H4 | Q96RI1 | 92 | 92 | SUCCEEDED |
| HNF4α | NR2A1 | P41235 | 8 | 8 | SUCCEEDED |
| HNF4γ | NR2A2 | Q14541 | 1 | 1 | SUCCEEDED |
| LRH-1 | NR5A2 | O00482 | 28 | 28 | SUCCEEDED |
| LXRα | NR1H3 | Q13133 | 7 | 7 | SUCCEEDED |
| LXRβ | NR1H2 | P55055 | 25 | 25 | SUCCEEDED |
| NUR77 | NR4A1 | P22736 | 22 | 22 | SUCCEEDED |
| NURR1 | NR4A2 | P43354 | 10 | 8 | SUCCEEDED |
| PPARα | NR1C1 | Q07869 | 79 | 79 | SUCCEEDED |
| PPARδ | NR1C2 | Q03181 | 56 | 56 | SUCCEEDED |
| PPARγ | NR1C3 | P37231 | 383 | 380 | SUCCEEDED |
| PXR | NR1I2 | O75469 | 80 | 80 | SUCCEEDED |
| Rev-erbα | NR1D1 | P20393 | 5 | 5 | SUCCEEDED |
| Rev-erbβ | NR1D2 | Q14995 | 6 | 6 | SUCCEEDED |
| RORα | NR1F1 | P35398 | 3 | 3 | SUCCEEDED |
| RORβ | NR1F2 | Q92753 | 0 | 0 | SUCCEEDED_ZERO_RESULTS |
| RORγ | NR1F3 | P51449 | 165 | 162 | SUCCEEDED |
| RXRα | NR2B1 | P19793 | 111 | 110 | SUCCEEDED |
| RXRβ | NR2B2 | P28702 | 7 | 7 | SUCCEEDED |
| RXRγ | NR2B3 | P48443 | 2 | 2 | SUCCEEDED |
| SF-1 | NR5A1 | Q13285 | 6 | 6 | SUCCEEDED |
| TR2 | NR2C1 | P13056 | 0 | 0 | SUCCEEDED_ZERO_RESULTS |
| TR4 | NR2C2 | P49116 | 5 | 5 | SUCCEEDED |
| COUP-TFI | NR2F1 | P10589 | 1 | 1 | SUCCEEDED |
| COUP-TFII | NR2F2 | P24468 | 1 | 1 | SUCCEEDED |
| DAX-1 | NR0B1 | P51843 | 1 | 1 | SUCCEEDED |
| EAR-2 | NR2F6 | P10588 | 1 | 1 | SUCCEEDED |
| GCNF | NR6A1 | Q15406 | 0 | 0 | SUCCEEDED_ZERO_RESULTS |
| NOR-1 | NR4A3 | Q92570 | 1 | 1 | SUCCEEDED |
| PNR | NR2E3 | Q9Y5X4 | 1 | 1 | SUCCEEDED |
| SHP | NR0B2 | Q15466 | 14 | 14 | SUCCEEDED |
| TLX | NR2E1 | Q9Y466 | 1 | 1 | SUCCEEDED |

## Receptors with zero candidates

RORβ, TR2, GCNF — all query_status=SUCCEEDED_ZERO_RESULTS (Search API succeeded, 0 matches), not API failures.

## Receptors with 1-2 candidates

EAR-2, COUP-TFI, COUP-TFII, RXRγ, DAX-1, HNF4γ, CAR, NOR-1, TLX, PNR

## Experimental-method breakdown (candidate rows)

- X-ray: 2055
- cryo-EM: 6
- NMR: 11
- Other experimental: 0

## Source taxonomy audit (candidates lacking taxon 9606)

9 candidates audited. EXPLAINED: 4, REVIEW_NEEDED: 5.

| common_name | pdb_id | entity | source | host | status |
|---|---|---|---|---|---|
| AR | 3ZQT | 3ZQT_1 | ["ESCHERICHIA COLI"] | ["ESCHERICHIA COLI"] | EXPLAINED |
| AR | 5CJ6 | 5CJ6_2 | ["synthetic construct"] | [] | EXPLAINED |
| ERα | 3CBM | 3CBM_2 | ["synthetic construct"] | [] | EXPLAINED |
| ERα | 8ARW | 8ARW_2 | ["Escherichia coli"] | [] | REVIEW_NEEDED |
| GR | 5CBX | 5CBX_1 | ["unclassified"] | ["Escherichia coli"] | REVIEW_NEEDED |
| GR | 5CBY | 5CBY_1 | ["UNCLASSIFIED"] | ["Escherichia coli"] | REVIEW_NEEDED |
| GR | 5CBZ | 5CBZ_1 | ["unidentified"] | ["Escherichia coli"] | REVIEW_NEEDED |
| PR | 5CC0 | 5CC0_1 | ["synthetic construct"] | ["Escherichia coli"] | EXPLAINED |
| ERRγ | 1VJB | 1VJB_1 | ["Mus musculus"] | ["Escherichia coli"] | REVIEW_NEEDED |

## Multi-UniProt entity audit

46 multi-mapped entities. By classification: {'FUSION_OR_CHIMERA': 46}

### AMBIGUOUS_REVIEW_NEEDED items

None.

## Cross-receptor PDB entries

35 PDB entries contain polymer entities mapped to more than one project receptor. This is ENTRY-LEVEL co-occurrence of distinct polymer entities, which is different from one POLYMER ENTITY mapping to multiple UniProt accessions (see `ENTITY_MULTI_MAPPED_UNIPROT` / the multi-UniProt audit above). No physical interaction/heterodimerization is claimed unless assembly-level contacts have been examined, which this stage does not do.

| pdb_id | receptors (uniprot_id) |
|---|---|
| 1DSZ | ["P10276","P19793"] |
| 1FM6 | ["P19793","P37231"] |
| 1FM9 | ["P19793","P37231"] |
| 1K74 | ["P19793","P37231"] |
| 1RDT | ["P19793","P37231"] |
| 1UHL | ["P28702","Q13133"] |
| 1XV9 | ["P19793","Q14994"] |
| 1XVP | ["P19793","Q14994"] |
| 1YNW | ["P11473","P19793"] |
| 1YUC | ["O00482","Q15466"] |
| 2NLL | ["P10828","P19793"] |
| 2Z4J | ["P10275","Q15466"] |
| 3DZU | ["P19793","P37231"] |
| 3DZY | ["P19793","P37231"] |
| 3E00 | ["P19793","P37231"] |
| 3H0A | ["P19793","P37231"] |
| 4DOR | ["O00482","Q15466"] |
| 4J5W | ["O75469","P19793"] |
| 4J5X | ["O75469","P19793"] |
| 4NQA | ["P19793","P55055"] |
| 4ONI | ["O00482","Q15466"] |
| 4RWV | ["O00482","P51843"] |
| 4ZO1 | ["P10828","P19793"] |
| 5HJP | ["P28702","P55055"] |
| 5I4V | ["P28702","P55055"] |
| 5JI0 | ["P19793","P37231"] |
| 5KYA | ["P28702","P55055"] |
| 5KYJ | ["P28702","P55055"] |
| 5UAN | ["P10826","P19793"] |
| 5Z12 | ["P19793","Q96RI1"] |
| 6A5Y | ["P19793","Q96RI1"] |
| 6A5Z | ["P19793","Q96RI1"] |
| 6A60 | ["P19793","Q96RI1"] |
| 6XWG | ["P10276","P19793"] |
| 8HBM | ["P19793","Q96RI1"] |

## Important statements

No final structural QC filters were applied at Stage 2 or Stage 2.1.

No mmCIF coordinate files were downloaded at Stage 2 or Stage 2.1.

## Stage 2.2 — Official Sequence Coordinates Validation

Three official, distinct RCSB services were used across Stage 2 / 2.1 / 2.2 / 2.3, each for a different purpose:
- **Search API** (`search.rcsb.org/rcsbsearch/v2/query`) = candidate discovery (Stage 2)
- **Data API** (`data.rcsb.org/graphql`) = entity/entry metadata enrichment (Stage 2, 2.1)
- **Sequence Coordinates API** (`sequence-coordinates.rcsb.org/graphql`) = independent sequence-alignment/coverage validation (Stage 2.2)

All 48 project UniProt accessions were queried once each (UNIPROT -> PDB_ENTITY), rather than once per polymer entity. 2075 experimental-PDB-entity target alignments were returned (computed structure models such as AlphaFold were identified by target ID and excluded); 2072 of the 2072 Stage 2 candidates were matched, 0 were missing, and 3 target(s) fell outside the Stage 2 Search API candidate scope. Stage 2.3 inspected these against authoritative RCSB Data API metadata (not titles) and found they are legitimately out of scope: their `structure_determination_methodology` is `integrative`, not `experimental` — Stage 2's Search query explicitly restricts to `experimental`, while the Sequence Coordinates API applies no such filter. They were NOT recently deposited (2018/2019) and their presence is not evidence of any indexing-timing discrepancy between RCSB services. See `rcsb_sequence_coordinates_extra_targets.csv` for the full per-target classification; they were not added to the Stage 2 candidate inventory.

Coverage agreement (official RCSB `query_coverage`/`target_coverage` vs. this project's independently computed union-interval coverage): MATCH=2072, MINOR_NUMERIC_DIFFERENCE=0, REVIEW_NEEDED=0 (max reference diff=4.98e-07, max entity diff=5.00e-07). Reference-interval overlaps found: 0; entity-interval overlaps found: 0 (out of 2118 normalized mapping rows checked) — where zero, the interval-union computation reduces exactly to the prior simple sum, i.e. the Stage 2.1 naive-sum coverage was already numerically correct for every such row; the union logic is retained as the robust general method going forward.

176 candidates were flagged as short mappings (reference coverage < 0.05 and/or aligned length < 30 residues) in `rcsb_short_mapping_audit.csv` — informational only, for Stage 3 to identify peptide/fragment constructs that are not full LBD structures.

`rcsb_candidate_inventory.csv` now carries BOTH official (`rcsb_reference_sequence_coverage`, `rcsb_entity_sequence_coverage`) and independently computed (`computed_reference_sequence_coverage`, `computed_entity_sequence_coverage`) coverage fields side by side; the official value is null (never guessed) where the Sequence Coordinates API did not return a match.

The Stage 2.1 taxonomy audit's `EXPLAINED` status describes only that the unusual source annotation has been characterized, NOT that the structure is suitable for the final human NR B-factor dataset — all 9 taxonomy-audit rows now carry `qc_suitability = NOT_EVALUATED`; no inclusion/exclusion decision has been made.

No final structural QC filters were applied at Stage 2, 2.1, 2.2, or 2.3.

No mmCIF coordinate files were downloaded at Stage 2, 2.1, 2.2, or 2.3.
