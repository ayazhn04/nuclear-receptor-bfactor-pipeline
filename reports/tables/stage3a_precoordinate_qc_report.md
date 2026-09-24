# Stage 3A — LBD Reference, Instance Mapping, and Pre-Coordinate QC

## 1. Purpose

Establish authoritative canonical-UniProt LBD coordinates for each of the 48 validated receptors, map every Stage 2 candidate's sequence-level coverage against that LBD (and, diagnostically, the DBD), resolve exact RCSB polymer-instance/chain identifiers, and compute descriptive pre-coordinate quality flags. No mmCIF coordinates were downloaded and no final structure inclusion/exclusion decision was made.

## 2. Authoritative LBD-coordinate method

Primary source: the current reviewed UniProt entry's `features[]` list, accepting only a feature with `type == "Domain"` and `description == "NR LBD"` (PROSITE-ProRule evidenced). All 48 receptors had exactly one such feature in their cached Stage 1 raw UniProt response — no PROSITE/InterPro fallback request was required. DBD interval (diagnostic only) uses `type == "DNA binding"`, `description == "Nuclear receptor"`; absent for DAX-1 and SHP, which are documented atypical orphan receptors lacking a classical DBD — no interval was inferred for them. Coordinates are 1-based, inclusive.

## 2a. Standardized LBD policy (Stage 3A.2)

The LBD interval used throughout this report is a **standardized project-wide canonical-sequence interval** derived from the current UniProt NR-LBD domain annotation. It is used because the same systematic annotation is available for all 48 receptors. It is **not asserted to be the only possible biological definition** of an LBD.

The previous cohort's VDR-only 118-427 range is preserved as a historical sensitivity definition (`data/manifests/nr_lbd_boundary_audit.csv`, `rcsb_lbd_boundary_sensitivity.csv`) but is **not used for primary filtering** anywhere in this report.

Functional motifs/sites (e.g. the VDR 9aaTAD/AF-2-associated motif) are represented separately in `data/manifests/nr_functional_features.csv` and may extend outside the standardized profile-defined interval. See `data/manifests/nr_lbd_analysis_policy.json` for the full policy record.

## 3. LBD reference table (all 48 receptors)

| common_name | nr_code | uniprot_id | LBD start-end | length | source | status |
|---|---|---|---|---|---|---|
| AR | NR3C4 | P10275 | 669-900 | 232 | UniProt | PASS |
| ERα | NR3A1 | P03372 | 311-547 | 237 | UniProt | PASS |
| ERβ | NR3A2 | Q92731 | 264-498 | 235 | UniProt | PASS |
| GR | NR3C1 | P04150 | 524-758 | 235 | UniProt | PASS |
| MR | NR3C2 | P08235 | 726-964 | 239 | UniProt | PASS |
| PR | NR3C3 | P06401 | 679-913 | 235 | UniProt | PASS |
| RARα | NR1B1 | P10276 | 183-417 | 235 | UniProt | PASS |
| RARβ | NR1B2 | P10826 | 183-417 | 235 | UniProt | PASS |
| RARγ | NR1B3 | P13631 | 185-419 | 235 | UniProt | PASS |
| TRα | NR1A1 | P10827 | 163-407 | 245 | UniProt | PASS |
| TRβ | NR1A2 | P10828 | 217-461 | 245 | UniProt | PASS |
| VDR | NR1I1 | P11473 | 127-423 | 297 | UniProt | PASS |
| CAR | NR1I3 | Q14994 | 109-352 | 244 | UniProt | PASS |
| ERRα | NR3B1 | P11474 | 193-421 | 229 | UniProt | PASS |
| ERRβ | NR3B2 | O95718 | 208-432 | 225 | UniProt | PASS |
| ERRγ | NR3B3 | P62508 | 233-457 | 225 | UniProt | PASS |
| FXR | NR1H4 | Q96RI1 | 262-486 | 225 | UniProt | PASS |
| HNF4α | NR2A1 | P41235 | 147-377 | 231 | UniProt | PASS |
| HNF4γ | NR2A2 | Q14541 | 99-328 | 230 | UniProt | PASS |
| LRH-1 | NR5A2 | O00482 | 300-539 | 240 | UniProt | PASS |
| LXRα | NR1H3 | Q13133 | 209-447 | 239 | UniProt | PASS |
| LXRβ | NR1H2 | P55055 | 222-460 | 239 | UniProt | PASS |
| NUR77 | NR4A1 | P22736 | 360-595 | 236 | UniProt | PASS |
| NURR1 | NR4A2 | P43354 | 360-595 | 236 | UniProt | PASS |
| PPARα | NR1C1 | Q07869 | 239-466 | 228 | UniProt | PASS |
| PPARδ | NR1C2 | Q03181 | 211-439 | 229 | UniProt | PASS |
| PPARγ | NR1C3 | P37231 | 238-503 | 266 | UniProt | PASS |
| PXR | NR1I2 | O75469 | 146-433 | 288 | UniProt | PASS |
| Rev-erbα | NR1D1 | P20393 | 284-614 | 331 | UniProt | PASS |
| Rev-erbβ | NR1D2 | Q14995 | 369-579 | 211 | UniProt | PASS |
| RORα | NR1F1 | P35398 | 272-510 | 239 | UniProt | PASS |
| RORβ | NR1F2 | Q92753 | 222-460 | 239 | UniProt | PASS |
| RORγ | NR1F3 | P51449 | 269-508 | 240 | UniProt | PASS |
| RXRα | NR2B1 | P19793 | 227-458 | 232 | UniProt | PASS |
| RXRβ | NR2B2 | P28702 | 296-529 | 234 | UniProt | PASS |
| RXRγ | NR2B3 | P48443 | 231-459 | 229 | UniProt | PASS |
| SF-1 | NR5A1 | Q13285 | 222-459 | 238 | UniProt | PASS |
| TR2 | NR2C1 | P13056 | 348-590 | 243 | UniProt | PASS |
| TR4 | NR2C2 | P49116 | 341-583 | 243 | UniProt | PASS |
| COUP-TFI | NR2F1 | P10589 | 184-410 | 227 | UniProt | PASS |
| COUP-TFII | NR2F2 | P24468 | 177-403 | 227 | UniProt | PASS |
| DAX-1 | NR0B1 | P51843 | 205-469 | 265 | UniProt | PASS |
| EAR-2 | NR2F6 | P10588 | 165-393 | 229 | UniProt | PASS |
| GCNF | NR6A1 | Q15406 | 249-480 | 232 | UniProt | PASS |
| NOR-1 | NR4A3 | Q92570 | 394-623 | 230 | UniProt | PASS |
| PNR | NR2E3 | Q9Y5X4 | 169-410 | 242 | UniProt | PASS |
| SHP | NR0B2 | Q15466 | 16-257 | 242 | UniProt | PASS |
| TLX | NR2E1 | Q9Y466 | 155-383 | 229 | UniProt | PASS |

## 4. Mapping-category counts

- NEAR_COMPLETE_LBD: 1731
- NO_LBD_OVERLAP: 225
- SUBSTANTIAL_LBD: 91
- TRACE_LBD: 23
- PARTIAL_LBD: 2

## 5. Experimental method counts

- xray: 2055
- cryo_em: 6
- nmr: 11
- other: 0

## 6. Resolution categories

- PASS_PROJECT_RANGE: 1690
- ULTRAHIGH_LT_1_8: 360
- TOO_LOW_RESOLUTION_GT_3_5: 11
- MISSING: 11

## 7. R-free categories

- PASS_ALL_LE_0_30: 1963
- FAIL_ALL_GT_0_30: 83
- MISSING: 26

## 8. Project metadata-prefilter count

- PASS: 1516
- NOT PASS: 556

## 9. Receptors with zero or very few LBD-mapped candidates

| common_name | near_complete_lbd | substantial_lbd | stage2_entities |
|---|---|---|---|
| CAR | 2 | 0 | 2 |
| ERRβ | 2 | 0 | 3 |
| HNF4γ | 1 | 0 | 1 |
| Rev-erbα | 0 | 0 | 5 |
| RORβ | 0 | 0 | 0 |
| RXRγ | 2 | 0 | 2 |
| TR2 | 0 | 0 | 0 |
| TR4 | 2 | 0 | 5 |
| COUP-TFI | 0 | 0 | 1 |
| COUP-TFII | 1 | 0 | 1 |
| DAX-1 | 0 | 0 | 1 |
| EAR-2 | 0 | 1 | 1 |
| GCNF | 0 | 0 | 0 |
| NOR-1 | 1 | 0 | 1 |
| PNR | 0 | 1 | 1 |
| SHP | 0 | 0 | 14 |
| TLX | 0 | 1 | 1 |

## 10. Short-fragment audit summary

22 of the 176 Stage 2 short mappings overlap the LBD; 154 do not (they lie elsewhere in the receptor sequence, e.g. AF-2/coactivator-groove peptides outside the canonical domain boundary, or non-LBD fragments). See `rcsb_short_mapping_lbd_audit.csv`.

## 11. Instance/chain mapping summary

- Total polymer instances: 3486
- Entities with one instance: 994
- Entities with multiple instances: 1078
- Maximum instances for one entity: 16
- label_asym_id <-> auth_asym_id pairing was resolved from true instance-level RCSB records (`CorePolymerEntityInstance`), never by pairing the Stage 2 entity-level array positions.

## 12. Taxonomy / multi-UniProt caveats

- 46 candidates are multi-UniProt-mapped entities (fusion/chimera constructs, carried forward — not excluded).
- 9 candidates carry a taxonomy-audit status from Stage 2.1/2.2 (not used as an automatic exclusion criterion).

## 13. Important statement

Sequence-mapping coverage is not equivalent to observed crystallographic coordinate coverage.

## 14. Important statement

No final structure inclusion/exclusion decision was made at Stage 3A.

## 15. Important statement

No mmCIF coordinate files were downloaded at Stage 3A.
