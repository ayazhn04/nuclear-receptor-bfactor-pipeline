# Receptor catalogue (v2)

One page per receptor showing how much data exists, its ligand/partner context, and where
flexibility varies between PDB entries. Built from the corrected v2 outputs and the v2
IQR work (branch `analysis/role-d-v2-iqr`, 37c1c15: professor canonical helix map,
corrected numbering, TRα on P10827-2). No analysis table is changed. The page shows a
visual selection; the summary table keeps every field.

## Outputs

| File | Content |
|---|---|
| `receptor_catalogue_v2.pdf` | 2 dataset-at-a-glance pages, then 48 receptor pages (NR-code order) |
| `figures/receptor_<NRcode>.png` | one 1920×1080 page per receptor |
| `figures/overview_glance_{1,2}.png` | dataset at a glance (24 receptors each) |
| `figures/overview_heterogeneity.png` | median vs 90th-percentile IQR per receptor (not in the deck) |
| `tables/receptor_catalogue_summary_v2.csv` | one row per receptor: every count, ligand, partner, construct, sequence-difference, crystal-form and IQR field, the summary line and the caveats |
| `tables/catalogue_entry_context_v2.csv` | one row per receptor × PDB entry |
| `tables/catalogue_partner_contacts_v2.csv` | every polymer chain within 5 Å of a receptor's LBD, with its class |
| `tables/catalogue_sequence_differences_v2.csv` | every residue that differs from the reference sequence, with its audited class and evidence |
| `tables/catalogue_mutation_audit_v2.csv` | engineered-mutant audit, one row per usable entry |
| `tables/catalogue_ligand_partner_contexts_v2.csv` | every ligand × partner and state × partner context per receptor, with counts |
| `tables/catalogue_partner_validation_vs_script27.csv` | agreement with the pipeline's partner table |
| `tables/catalogue_rules_v2.json` | every threshold used for labels, classes, the construct line and caveats |

Rebuild (needs the local mmCIF caches, see `scripts/108_raw_cache_requirements.py`):

```
python scripts/112_catalogue_entry_context.py
python scripts/113_build_receptor_catalogue.py
```

## Receptor page

1. **Header**: receptor, NR code, number of usable PDB entries, and one summary line
   (composition · ligands · partners · IQR heterogeneity).
2. **Structures**: one bar of HOLO / APO / ambiguous / TRACE_LBD (TRACE_LBD-only entries
   are receptor fragments, kept in the profiles but never counted as APO or HOLO); entries
   found and with coordinates; median fraction of the LBD each entry observes. "no APO"
   is flagged when a receptor has no APO entry.
3. **Context**: number of distinct functional ligands with the top three; one partner
   bar (coregulator peptide, NR / heterodimer partner, other partner, standalone; one
   class per entry, NR partner first), or "Standalone: 100%".
4. **Flexibility** (main plot): median per-entry z per canonical position (line), IQR
   across entries where ≥5 entries observe the position (band), dashed line where 1–4
   entries observe it, grey bands where no structure observes the sequence, professor
   helices (hatched = review_required, orange = H12 / AF-2), pocket ticks (geometric
   contact consensus) and groove ticks (curated charge-clamp residues), and the coverage
   track (red line = 5 structures). The z axis is fixed (−2.5 to 5) on every page.
   Beside the axis title: "z-scored / total entries z-scored (%)" (see "Z-score coverage"
   below); "incl. N TRACE_LBD" or "all TRACE_LBD" when TRACE_LBD-only fragments are among
   the z-scored entries.
5. **Between-PDB heterogeneity**: class (tertile of median position IQR among the 30
   receptors with ≥20 supported positions: low ≤ 0.459 < moderate ≤ 0.562 < high),
   median and 90th-percentile IQR, supported positions; "not assessable" with the reason
   otherwise. Same IQR definition as `reports/role_d/contextual_iqr` (v2 profile, ≥5 PDB
   entries).
6. **Constructs & crystals**: construct boundary patterns (full UniProt alignment of each
   construct), crystal forms, resolution; engineered mutants, fusions and internal
   deletions only when notable.
7. **Caveats**: at most three receptor-specific points that could change interpretation.

## Added context

- **Engineered-mutant audit** (`tables/catalogue_mutation_audit_v2.csv`, one row per usable
  entry; residue detail in `tables/catalogue_sequence_differences_v2.csv`). Every residue of
  the representative chain that differs from the reference sequence (corrected v2 numbering,
  so position errors are already removed) gets one class:
  - `ENGINEERED`: the depositor marks it "engineered mutation" in `_struct_ref_seq_dif` and
    its reference residue equals ours, and/or the entity's `_entity.pdbx_mutation` lists the
    same residues at that position;
  - `NATURAL_VARIANT_OR_ISOFORM`: depositor "variant", or a UniProt natural variant or
    single-residue isoform difference;
  - `MODIFIED_RESIDUE`: non-standard component (MSE, CSO, ...);
  - `MAPPING_REFERENCE_DIFFERENCE`: "engineered" but with a different reference residue
    (none found);
  - `UNCERTAIN`: "conflict", other annotations, or none.

  Only `ENGINEERED` counts. Per entry: `ENGINEERED_MUTANT`, `WT_NO_ENGINEERED_SUBSTITUTION`,
  or `UNCERTAIN` (only uncertain differences, or an entity-declared substitution that no
  observed residue confirms). Of 1,418 differences: 980 engineered (716 confirmed by both
  sources, 264 by `_struct_ref_seq_dif` alone), 18 natural variant/isoform, 223 modified,
  197 uncertain (192 depositor "conflict"). Entries: 537 engineered mutant, 1,158 without
  an engineered substitution, 110 uncertain. Several "conflict" differences look like
  engineered surface or cysteine mutations (FXR E291A/E364A, MR C808S/C910S); they are not
  counted, so engineered fractions are lower bounds where uncertainty is high (FXR 49%).
- **Partner identities** (`scripts/utils/partner_identity.py`): the deposited description
  and UniProt accession stay in `catalogue_partner_contacts_v2.csv`; `partner_identity` is
  a separate normalized gene-level name. The rules are applied in this order: explicit
  accession table, name synonyms (e.g. NCOA2 = SRC-2 = TIF2 = GRIP1), NR-box sequence motif,
  generic labels counted as one "unidentified peptide", accessions learned from name-matched
  rows, else the description itself (uncertain identities are not merged). DNA and other
  entities of the same receptor are not partners. 180 deposited names → 57 identities.
- **Ligand × partner contexts** (`tables/catalogue_ligand_partner_contexts_v2.csv`): per
  entry with an LBD context (TRACE_LBD-only excluded), ligand state × partner category and
  functional-ligand CCD id(s) × partner identity(ies); several ligands or partners are
  joined in sorted order, so combinations are deterministic and synonyms cannot inflate
  counts.
- **Crystallographic metadata** from each entry's own mmCIF: space group, unit cell,
  deposition date, polymer chains and receptor copies in the asymmetric unit. A crystal
  form = same space group and unit cell within 2% (lengths) / 2° (angles).
- **Partner context** for all 1,810 coordinate entries with the pipeline's own 5 Å
  contact rule and classifier (`scripts/27`). The pipeline table covers 1,395 primary-QC
  entries and skips other project receptors by design; the catalogue reproduces all 726
  of its chains and adds heterodimer partners and non-primary entries.

## Dataset counts: reconciliation

Five counts appear across the project. They count different units, so they are not
inconsistencies. Each is nested inside the one before it:

| Count | Unit | Source |
|---|---|---|
| 2,072 | receptor × PDB entry × polymer entity: discovery candidates (one row each, 45 receptors, 1,996 distinct PDB IDs) | rows of `data/rebuild_v2/processed/structures.csv` |
| 2,031 | unique receptor × PDB entries found (entities of one receptor in the same entry collapsed) | unique (`uniprot_id`, `pdb_id`) in `structures.csv` = rows of `catalogue_entry_context_v2.csv` = sum of `n_entries_total` |
| 1,844 | receptor × PDB entry × polymer entity with coordinate data (entity-level coordinate candidates) | `structures.csv` rows with `has_stage3b_coordinate_data` |
| 1,810 | unique receptor × PDB entries with coordinate data | unique (`uniprot_id`, `pdb_id`) of those rows = sum of `n_entries_with_coordinates` |
| 1,805 | unique receptor × PDB entries with ≥ 1 finite B-factor z-score | unique (`uniprot_id`, `pdb_id`) with finite `b_factor_z` in `data/rebuild_v2/analysis/pdb_position_bfactor_profiles.parquet` = sum of `n_entries_with_zscores` |

The two differences that change the unit (entities collapsed into one entry) are
2,072 → 2,031 (41 extra entity rows) and 1,844 → 1,810 (34). The differences between
entry-level counts are the excluded entries: 2,031 → 1,810 = 221 without coordinate data
(204 receptor chain outside the LBD, 17 not X-ray), and 1,810 → 1,805 = 5 with
coordinates but 0 usable CA atoms. Catalogue pages and tables use only the three
entry-level counts (2,031 / 1,810 / 1,805).

## Z-score coverage

The label beside each z plot and the first glance column use one denominator:

- **total** (`n_entries_total`, = `n_entries_discovered`): every receptor × PDB entry found
  at discovery (2,031), including NMR/EM entries and entries whose receptor chain lies
  outside the LBD. The page's "N found · M with coordinates" line uses the same total.
- **with coordinates** (`n_entries_with_coordinates`): entries whose LBD coordinates entered
  the structural dataset (1,810).
- **z-scored** (`n_entries_with_zscores`, = `n_entries_contributing`): entries with at least
  one finite per-chain B-factor z-score in the v2 profiles (1,805). These are the "usable
  PDBs" in the header and the denominator of every composition value.

Every entry without z-scores (226) has exactly one reason column:
`n_not_zscored_receptor_chain_outside_lbd` (204), `n_not_zscored_not_xray` (17),
`n_not_zscored_no_usable_ca` (5: VDR 1KB2/1KB4/1KB6/1YNW and RXRα 1YNW; 0 usable CA in
the LBD, all TRACE_LBD-only), `n_not_zscored_too_few_values` (< 2 usable B; 0),
`n_not_zscored_zero_sd` (0), `n_not_zscored_no_coordinates_other` and
`n_not_zscored_other` (0). `fraction_entries_zscored` is empty when no entry was found
(RORβ, TR2, GCNF); the label shows ">99%", never "100%",
when an entry is missing. TRACE_LBD: `n_trace_lbd_only` (24 found) =
`n_trace_lbd_zscored` (19, kept in the profiles; e.g. GR 4, SHP 14, ERRγ 1) +
`n_trace_lbd_not_zscored` (5); `n_zscored_lbd_context` excludes them. Normalization and
TRACE_LBD handling are unchanged.

## Definitions and caveats

- Denominator of every composition value: the receptor's entries that contribute z-scores.
- Entry = receptor × PDB entry.
- The per-position profiles pool all ligand states and TRACE_LBD fragments; they describe
  the deposited data, not an apo/holo or burial-controlled comparison.
- IQR is between-PDB crystallographic heterogeneity, not motion in solution, and is not
  yet controlled for SASA/burial. Ligand and partner context are associations only.

## Known pre-existing test failure

`tests/test_final_freeze_consolidation.py::test_freeze_pipeline_deterministic` fails on a
clean checkout of `37c1c15` (the base of this branch), before any catalogue file exists.
The identity audit (`scripts/93`) flags ERα and PPARγ wording in the committed Role D
literature markdown and rewrites `reports/tables/final_receptor_identity_audit.csv`. With the
catalogue files present the failure and the diff are identical; they add no new failure. Restore that file after running the test
(`git restore -- reports/tables/final_receptor_identity_audit.csv`).
