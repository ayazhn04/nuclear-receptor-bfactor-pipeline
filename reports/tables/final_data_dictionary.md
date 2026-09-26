# Final Data Release — Data Dictionary

> **QC criteria are metadata/analysis flags. Structures are retained in the
> master inventory rather than deleted based on coverage thresholds.**

Nested views (not deletions): 2072 discovery candidates (master
`structures.csv`) -> 1840 with Stage 3B coordinate data -> 1382 previous
strict primary-QC subset. **1382 is not the master dataset.**

Tables: master `structures.csv` and `lbd_residue_bfactors_all.parquet`
(documented first), then the historical strict-subset tables
(`structures_primary_qc_subset.csv` — same columns as the section titled
"structures.csv — ONE ROW PER PRIMARY SELECTED STRUCTURE" below, which now
describes that file; `final_lbd_residue_map_primary_qc_subset.parquet` —
same columns as the section titled "final_lbd_residue_map.parquet" below),
`excluded_structures.csv`, and `final_ligand_annotations.csv`. See
[`TEAM_HANDOFF.md`](../../TEAM_HANDOFF.md) for which file to use for what.

Missing-value convention throughout: an empty string / `NaN` means the
value was not available or not applicable for that row — never a silently
substituted default. A JSON-array column with `[]` means "computed, and
the list is genuinely empty" (not missing).

## MASTER `data/processed/structures.csv` — ONE ROW PER STAGE 2 CANDIDATE (2072)

One row per (project receptor x PDB entry x polymer entity). Nothing is
removed for coverage, resolution, R-free, method, taxonomy or fusion status.
Where a value cannot be computed (no Stage 3B coordinate data) it is left
blank/NA and `coordinate_data_status` says why.

| Column | Type | Meaning / missing semantics |
|---|---|---|
| `uniprot_id`, `nr_code`, `common_name`, `group` | str | Receptor identity |
| `pdb_id`, `polymer_entity_id`, `entity_id` | str | Candidate identity (unique key: uniprot_id + pdb_id + polymer_entity_id) |
| `experimental_method` | str | Deposited method(s) |
| `has_stage3b_coordinate_data` | bool | True for the 1840 X-ray + LBD-overlap candidates whose coordinates were audited |
| `coordinate_data_status` | str | `COORDINATE_DATA_AVAILABLE` / `NO_STAGE3B_COORDINATE_DATA` |
| `method_is_xray`, `mapping_overlaps_lbd` | bool | The two conditions for entering Stage 3B |
| `lbd_mapping_category` | str | Stage 3A construct-vs-LBD mapping category |
| `standardized_lbd_start/_end/_length` | int | Standardized LBD (UniProt canonical, 1-based inclusive) |
| `construct_lbd_mapping_coverage` | float | Sequence-mapping coverage of the LBD (not observed coordinates) |
| `positive_occupancy_lbd_coverage` | float | Observed (occupancy > 0) Ca fraction of the standardized LBD, at the representative instance; **blank if no coordinate data** |
| `bfactor_usable_lbd_coverage` | float | Same with finite B-factor; blank if no coordinate data |
| `coverage_ge_080/090/095/099`, `coverage_eq_100` | nullable bool | Threshold flags on `positive_occupancy_lbd_coverage`; **blank (NA), not False,** when no coordinate data |
| `candidate_instance_count` | int | Measured instances for this candidate; blank if none |
| `resolution` / `resolution_project_range_status` | float / str | Angstrom; status `PASS_PROJECT_RANGE`, `ULTRAHIGH_LT_1_8`, `TOO_LOW_RESOLUTION_GT_3_5`, `MISSING` |
| `r_free` / `r_free_status` | float / str | Status `PASS_ALL_LE_0_30`, `FAIL_ALL_GT_0_30`, `MISSING` |
| `identity_class` | str | `HUMAN_OR_HUMAN_DERIVED`, `NON_HUMAN_ORTHOLOG`, `ANCESTRAL_RECONSTRUCTION`, `NON_HUMAN_CONSTRUCT`, `AMBIGUOUS_IDENTITY` |
| `fusion_or_chimera` | bool | Multi-UniProt fusion/chimera entity |
| `taxonomy_audit_status` | str | Stage 2.1 audit status for unusual-taxonomy entities; blank otherwise |
| `previous_primary_qc_member` | bool | Member of the historical strict 1382-row subset |
| `previous_primary_exclusion_reasons_json` | JSON list | All reasons it was outside that subset (`[]` for members) — informational, NOT a deletion reason |
| `representative_instance_id`, `label_asym_id`, `auth_asym_id` | str | Deterministic representative instance (coverage-first rule); blank if no coordinate data |
| `ligand_state` | str | `APO`/`HOLO`/`AMBIGUOUS` for all 1840 coordinate-assessable candidates (independent of coverage/resolution/R-free); `NOT_ASSESSABLE_NO_COORDINATE_DATA` for the 232 others (never APO) |
| `functional_ligand_component_ids` | JSON list | CCD IDs classified `FUNCTIONAL_LBD_LIGAND`; blank when not assessable |
| `ligand_annotation_status` | str | `RESOLVED`/`REVIEW_NEEDED`, or `NOT_ASSESSABLE_NO_COORDINATE_DATA` |
| `mmcif_sha256` | str | SHA256 of the cached mmCIF; blank if no coordinate data |
| `in_completeness_95/99/80_sensitivity`, `in_ultrahigh_resolution_sensitivity`, `in_no_fusion_sensitivity` | bool | Sensitivity-view membership from the earlier selection audit |

## `data/processed/lbd_residue_bfactors_all.parquet` — ALL corrected Stage 3B residue observations

771,542 rows = 3159 measured instance/models x every standardized-LBD
position. No filtering by coverage, resolution, R-free, or subset
membership; unmapped and unobserved positions are explicit rows. Raw
values, **no normalization**.

Columns: everything in the Stage 3B observation table (`uniprot_id, nr_code,
common_name, pdb_id, polymer_entity_id, entity_id, instance_id,
label_asym_id, auth_asym_id, model_num, canonical_uniprot_position,
lbd_relative_position, entity_label_seq_id, is_construct_mapped,
ca_atom_record_count, *_json altloc/occupancy/B/auth-numbering arrays,
mapping_status, observation_status, ca_record_present, ca_positive_occupancy,
ca_bfactor_usable, zero_occupancy_only_ca, positive_occupancy_ca_altloc_count`;
Stage 3B.1 positive-occupancy semantics unchanged — `ca_positive_occupancy`
is "observed"; zero-occupancy placeholders are not) plus:

| Column | Meaning |
|---|---|
| `canonical_residue_name` | 3-letter code from the UniProt canonical sequence |
| `is_representative_instance` | True for the deterministic representative instance of a candidate (use to get one instance per receptor x PDB x entity) |
| `selected_ca_altloc` / `selected_ca_occupancy` / `selected_ca_b_iso_or_equiv` | One recommended raw Ca record (highest occupancy; ties: blank altloc, then A, then lexicographic; B magnitude never used); null if no occupancy>0 finite-B record |

The legacy alias column `ca_modeled` (== `ca_record_present`) is dropped to
avoid confusion with the "observed" definition. Join to the master on
`instance_id` = `representative_instance_id` or on (uniprot_id, pdb_id,
polymer_entity_id).

## `data/processed/structures_primary_qc_subset.csv` (historical strict subset) — "structures.csv — ONE ROW PER PRIMARY SELECTED STRUCTURE"

One row per (project receptor x PDB entry x receptor polymer entity) that
passed every primary-selection rule, at its single deterministically
selected representative instance. 1382 rows.

| Column | Type | Meaning | Units / allowed values |
|---|---|---|---|
| `uniprot_id` | str | Project receptor's UniProt accession | e.g. `P11473` |
| `nr_code` | str | Nuclear receptor nomenclature code | e.g. `NR1I1` |
| `common_name` | str | Receptor common name | e.g. `VDR` |
| `group` | str | Receptor subfamily group (from professor-supplied metadata) | |
| `pdb_id` | str | PDB entry ID | 4-character PDB code |
| `polymer_entity_id` | str | RCSB polymer entity ID | `{pdb_id}_{entity_id}` |
| `entity_id` | str/int | mmCIF entity number within the entry | |
| `instance_id` | str | Selected representative instance ID | `{pdb_id}.{label_asym_id}` |
| `label_asym_id` | str | mmCIF chain identifier (canonical, use this for programmatic joins) | |
| `auth_asym_id` | str | Author/PDB-deposited chain identifier (provenance/display only — never join on this) | |
| `experimental_method` | str | Deposited experimental method | always `X-RAY DIFFRACTION` in the strict subset by primary policy |
| `resolution` | float | Reported crystallographic resolution | Angstrom; always in [1.8, 3.5] in this table |
| `r_free` | float | Reported R-free | fraction; always <= 0.30 in this table |
| `standardized_lbd_start` / `_end` / `_length` | int | Standardized LBD interval (UniProt canonical, 1-based inclusive; PROSITE PS51843-backed, "Option A") | residue positions / residue count |
| `construct_lbd_mapping_coverage` | float | Fraction of standardized LBD covered by the deposited construct's sequence (regardless of observed coordinates) | 0-1 |
| `positive_occupancy_lbd_coverage` | float | Fraction of standardized LBD with an observed (occupancy > 0) Ca | 0-1; always >= 0.90 in this table |
| `bfactor_usable_lbd_coverage` | float | Fraction of standardized LBD with an occupancy > 0 AND finite-B-factor Ca | 0-1 |
| `mapped_lbd_residues` | int | Count of standardized LBD residues the construct's sequence mapping reaches | |
| `positive_occupancy_lbd_residues` | int | Count of standardized LBD residues with an observed Ca | |
| `mapped_but_unobserved_lbd_residues` | int | Mapped positions with no positive-occupancy Ca (record-absent + zero-occupancy-only, combined) | |
| `zero_occupancy_only_lbd_residues` | int | Mapped positions whose only Ca record(s) have occupancy == 0 (deposited placeholder, not an observation) | |
| `internal_missing_run_count` / `longest_internal_missing_run` | int | Count / length of contiguous mapped-but-unobserved runs that do NOT touch the standardized LBD's own N-/C-terminal edge (i.e. genuinely internal disorder, not ordinary terminal truncation) | residues |
| `anisotropic_category_present` | bool | Whether the deposited entry includes `_atom_site_anisotrop` (anisotropic B-factor refinement) | |
| `identity_class` | str | `HUMAN_OR_HUMAN_DERIVED` always, in this table (only class allowed in primary) | |
| `fusion_or_chimera` | bool | Whether this polymer entity is a multi-UniProt fusion/chimera construct (Stage 2.1 audit) | |
| `source_organism` | str (JSON list) | Deposited source organism name(s) | |
| `taxonomy_audit_status` | str | Stage 2.1 taxonomy-audit status if this entity was ever flagged unusual; blank otherwise | `EXPLAINED` / `REVIEW_NEEDED` / blank |
| `ligand_state` | str | APO / HOLO / AMBIGUOUS classification for this representative instance | `APO`, `HOLO`, `AMBIGUOUS` |
| `functional_ligand_component_ids` | str (JSON list) | CCD component ID(s) classified `FUNCTIONAL_LBD_LIGAND` for this instance | |
| `functional_ligand_names` | str (JSON list) | Chemical name(s) of the above | |
| `coregulator_or_partner_present` | bool | Whether any polymer partner (coregulator peptide, other receptor, other protein) contacts this instance's LBD (see `final_polymer_partner_annotations.csv`) | |
| `polymer_partner_summary` | str | Human-readable one-line summary of contacting polymer partner(s), if any | |
| `passes_95pct_completeness` / `passes_99pct_completeness` | bool | Whether this structure would ALSO pass a stricter completeness threshold (still requires all other primary criteria) | |
| `passes_no_fusion_sensitivity` | bool | True for every row in this table except fusion/chimera constructs (i.e. `in_primary_set AND NOT fusion_or_chimera`) | |
| `is_ultrahigh_resolution` | bool | Always `False` in this table by construction — resolution < 1.8 A structures are never primary; see excluded_structures.csv / `in_ultrahigh_resolution_sensitivity` in final_structure_selection_audit.csv for the preserved ultrahigh-resolution sensitivity records | |
| `mmcif_sha256` | str | SHA256 of the exact downloaded `.cif.gz` this row derives from | |
| `selection_policy_version` | str | Policy version tag (`final-structure-selection-v1`) | |

## `data/processed/excluded_structures.csv` — EVERY NON-PRIMARY STAGE 2 CANDIDATE

690 rows = 2072 Stage 2 candidates minus the 1382 primary rows. Answers
"why isn't PDB X in structures.csv?" without rerunning the pipeline.

Shares most QC columns with `structures.csv` (same meaning). Distinct columns:

| Column | Type | Meaning |
|---|---|---|
| `stage3b_coordinate_selected` | bool | Whether this candidate ever entered the Stage 3B coordinate-completeness audit (X-ray + LBD-overlap) |
| `representative_instance_id` | str | Best-coverage instance found during representative selection, if any (blank if `stage3b_coordinate_selected` is False) |
| `completeness_ge_090` .. `completeness_eq_100` | bool | Whether the best/representative instance reaches each threshold |
| `primary_method_pass` / `primary_completeness_pass` / `primary_resolution_pass` / `primary_rfree_pass` / `primary_identity_pass` | bool | Per-criterion pass/fail |
| `primary_selection_status` | str | `NO_COORDINATE_CANDIDATE` (never entered Stage 3B pool) or `EXCLUDED_PRIMARY` (entered but failed >=1 criterion) |
| `primary_exclusion_reasons_json` | str (JSON list) | ALL applicable reason codes (never just the first) — see `reports/tables/final_data_dictionary.md`'s reason vocabulary below |

**Exclusion reason vocabulary** (multiple may apply to one row):
`NON_XRAY`, `NO_LBD_OVERLAP`, `LBD_COVERAGE_LT_0_90`, `RESOLUTION_LT_1_8`,
`RESOLUTION_GT_3_5`, `RESOLUTION_MISSING`, `RFREE_GT_0_30`,
`RFREE_MISSING`, `NON_HUMAN_ORTHOLOG`, `ANCESTRAL_RECONSTRUCTION`,
`AMBIGUOUS_IDENTITY`, `NON_HUMAN_CONSTRUCT`.

## `data/manifests/final_ligand_annotations.csv` — ONE ROW PER COORDINATE CANDIDATE

1840 rows, one per master `structures.csv` row with `has_stage3b_coordinate_data`
(joins on `pdb_id` + `uniprot_id` + `polymer_entity_id`); the 1382 strict-subset
rows are unchanged from the earlier release. Contact geometry uses observed
(occupancy > 0) LBD atoms of the representative instance; if none exist the row is
`AMBIGUOUS`/`REVIEW_NEEDED` rather than APO.

| Column | Type | Meaning |
|---|---|---|
| `ligand_state` | str | `APO` / `HOLO` / `AMBIGUOUS` |
| `functional_ligand_component_ids_json` / `functional_ligand_names_json` | JSON list | Ligand(s) driving a HOLO call |
| `all_lbd_contacting_nonpolymer_ids_json` | JSON list | Every nonpolymer instance ID (any relevance class) within 4.5 A of the LBD, for full transparency |
| `ligand_annotation_status` | str | `RESOLVED` or `REVIEW_NEEDED` |
| `ligand_annotation_evidence` | str | Machine-generated one-line tally of the evidence behind the call |
| `ligand_annotation_notes` | str | Human-readable caveat when `REVIEW_NEEDED` |

## `data/manifests/final_nonpolymer_inventory.csv` — EVERY NONPOLYMER COMPONENT/INSTANCE

262,783 rows for all 1777 coordinate-candidate PDB entries (257,731 of them water, enumerated once per PDB entry —
water's classification never depends on geometry or receptor context, see
`scripts/utils/nonpolymer_classification.py`). All other nonpolymer
instances are evaluated once per primary representative-receptor context
in that PDB (relevant for PDB entries with two primary
representative receptors, e.g. RXR heterodimer partners, where the same
physical molecule can be functionally relevant to one receptor's LBD and
irrelevant to the other's).

| Column | Type | Meaning |
|---|---|---|
| `nonpolymer_relevance_class` | str | `WATER`, `ION`, `CRYSTALLIZATION_ADDITIVE_OR_BUFFER`, `SOLVENT`, `COVALENT_OR_MODIFIED_COMPONENT`, `FUNCTIONAL_LBD_LIGAND`, `POSSIBLE_FUNCTIONAL_LIGAND`, `OTHER_NONPOLYMER`, or `REVIEW_NEEDED` |
| `nonpolymer_relevance_evidence` | str | The specific evidence combination that produced the classification |
| `min_heavy_atom_distance_to_lbd` | float | Minimum heavy-atom (positive-occupancy, non-H) distance to the representative receptor's standardized-LBD heavy atoms | Angstrom; blank for water |
| `contact_within_4_5A` / `contact_within_5_0A` | bool | Contact flags at each cutoff |
| `contact_lbd_canonical_positions_json` | JSON list | Canonical UniProt LBD positions within 4.5 A |
| `rcsb_subject_of_investigation` | bool | RCSB's own per-instance "ligand of interest" annotation (present for some, not all, true ligands — an incomplete but authoritative positive signal, never used as a negative signal) |
| `struct_site_record` | bool | Deposited `_struct_site` binding-site record names this exact component |
| `covalently_linked` | bool | Deposited `_struct_conn` covalent linkage involves this instance |

## `data/manifests/final_polymer_partner_annotations.csv` — CONTACTING POLYMER PARTNERS

724 rows — only polymer chains (other than the receptor itself) within 5.0
A of a primary representative receptor's LBD are recorded (informational
only; never alters primary selection).

| Column | Type | Meaning |
|---|---|---|
| `partner_category` | str | `COREGULATOR_CANDIDATE`, `OTHER_RECEPTOR_PARTNER`, `OTHER_PROTEIN_PARTNER`, or `UNCLASSIFIED_PARTNER` |
| `partner_category_evidence` | str | Why (UniProt match to a project receptor, coregulator-keyword description match, or neither) |
| `partner_mapped_uniprot_ids_json` | JSON list | UniProt accession(s) this partner entity maps to, if any |
| `minimum_distance_to_receptor_lbd` | float | Angstrom |
| `contacting_receptor_canonical_positions_json` | JSON list | Canonical LBD positions contacted |

## `data/processed/final_lbd_residue_map_primary_qc_subset.parquet` (historical strict subset) — "final_lbd_residue_map.parquet — CANONICAL RESIDUE HANDOFF"

338,012 rows: one row per (primary representative structure x standardized
LBD position), including every position with no observed coordinate — the
canonical residue-numbering and raw-Ca handoff table.

| Column | Type | Meaning |
|---|---|---|
| `canonical_uniprot_position` | int | 1-based UniProt canonical position — **the only numbering downstream code should join on** |
| `lbd_relative_position` | int | 1-based position within the standardized LBD interval (`canonical_uniprot_position - standardized_lbd_start + 1`) |
| `entity_label_seq_id` | float (nullable int) | mmCIF `label_seq_id` for this residue in this entity, if construct-mapped |
| `auth_seq_ids_json` / `insertion_codes_json` | JSON list | Deposited author numbering, for provenance/display only — never a join key |
| `canonical_residue_name` | str | 3-letter amino acid code from the UniProt canonical sequence at this position |
| `observed_residue_names_json` | JSON list | 3-letter code(s) actually observed at this position (per altloc) — may differ from canonical (engineered mutations, etc.) |
| `is_construct_mapped` | bool | Whether this position is covered by this instance's deposited construct |
| `ca_record_present` / `ca_positive_occupancy` / `ca_bfactor_usable` | bool | The three Stage 3B.1 Ca-observation tiers (record-present includes zero-occupancy placeholders; positive-occupancy is "observed"; bfactor-usable additionally requires a finite B) |
| `zero_occupancy_only_ca` | bool | Every Ca record at this position has occupancy == 0 |
| `altloc_ids_json` / `occupancies_json` / `b_iso_values_json` | JSON list, parallel | ALL deposited Ca altloc records at this position, unfiltered — the full raw data |
| `selected_ca_altloc` / `selected_ca_occupancy` / `selected_ca_b_iso_or_equiv` | str/float/float, nullable | The ONE recommended Ca record, chosen deterministically by occupancy then altloc identity only (see `scripts/utils/altloc_selection.py`) — **B-factor magnitude never participates in this choice**; all three are null if no occupancy>0-and-finite-B record exists at this position |
