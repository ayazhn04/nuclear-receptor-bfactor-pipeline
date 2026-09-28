# Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families

**Course:** BIOL363, Nazarbayev University

## Current repository purpose

**CORE ANALYSIS COMPLETE.** The Data & Infrastructure pipeline, full-dataset
B-factor normalization, structure-based canonical LBD segment mapping,
segment-level PCA/UMAP/clustering, subfamily/group statistical testing, a
superfamily flexibility atlas, an apo-vs-holo flexibility comparison, a
current full-dataset functional-site flexibility map, and a specificity/
confounding validation layer for that map are all complete, tested, and
frozen. See [`reports/FINAL_ANALYSIS_RESULTS.md`](reports/FINAL_ANALYSIS_RESULTS.md)
for the master evidence-based record of every final scientific result,
[`reports/FINAL_METHODS_AUDIT.md`](reports/FINAL_METHODS_AUDIT.md) for the
exact method behind each result, and
[`reports/FINAL_LIMITATIONS.md`](reports/FINAL_LIMITATIONS.md) for what those
results do and do not support. **Not yet performed**: ClinVar/disease
mutation analysis. See [`TEAM_HANDOFF.md`](TEAM_HANDOFF.md) for the current
master dataset and a guided tour of every analysis layer.

This repository distinguishes five layers, from most to least current —
never delete or silently overwrite an earlier one:

1. **MASTER DATA** — `data/processed/structures.csv` (2072 rows) and
   `lbd_residue_bfactors_all.parquet` (771,542 rows). QC criteria are
   metadata/analysis flags here, never deletion criteria.
2. **CURRENT FINAL ANALYSES** — everything built on the master data with the
   approved per-instance normalization, chain→PDB-median aggregation, and
   the v3 structural-segment mapping: normalization foundation, v3 segment
   mapping, segment-level PCA/UMAP/clustering, subfamily testing +
   flexibility atlas, apo/holo comparison, and the current full-dataset
   functional-site flexibility map + its specificity validation. This is
   what `reports/FINAL_ANALYSIS_RESULTS.md` reports.
3. **HISTORICAL STRICT-QC VIEW** — the 1,382-structure
   `structures_primary_qc_subset.csv` / `final_lbd_residue_map_primary_qc_subset.parquet`,
   preserved for reproducibility, not the master dataset for new analyses.
4. **HISTORICAL ROLE D ANALYSIS** — `reports/tables/functional_site_definitions.csv`
   and its downstream scores/QC files: a hand-curated, 5-receptor descriptive
   analysis on the strict-QC subset, byte-identical since commit `e714bbd`,
   superseded (not replaced or deleted) by the current full-dataset
   functional-site analysis.
5. **DEPRECATED DIAGNOSTIC v1/v2 STRUCTURAL MAPPING** — the ordinal-
   helix-counting v1 mapping (commit `9f15062`) and the incomplete-template
   v2 structural-alignment mapping (commit `abd2257`), both preserved in git
   history and in `structural_segment_validation.csv`/`_v2.csv`, superseded
   by the current v3 mapping.

Stages completed:

- **Stage 0** — repository/environment bootstrap
- **Stage 1 / 1.5** — professor metadata ingestion, UniProt identity validation, frozen receptor master
- **Stage 2 / 2.1 / 2.2 / 2.3 / 2.4** — RCSB Search API candidate discovery, taxonomy/fusion audits, Sequence Coordinates cross-validation, frozen candidate inventory
- **Stage 3A / 3A.1 / 3A.2** — standardized LBD definition (UniProt/PROSITE PS51843), canonical residue mapping, frozen LBD reference
- **Stage 3B / 3B.1 / 3B.2** — mmCIF acquisition (1777 files), observed-coordinate completeness audit (positive-occupancy Ca model), threshold-attrition and residual-missingness audit
- **Final Data Release** — structure selection, deterministic representative-instance selection, ligand-state annotation, coregulator/partner annotation, canonical residue handoff

> **QC criteria are metadata/analysis flags. Structures are retained in the
> master inventory rather than deleted based on coverage thresholds.**

Nested analysis views (counts are views, **not** deletions):

```
2072 Stage 2 discovery candidates  = data/processed/structures.csv (MASTER)
 -> 1840 X-ray + LBD-overlap coordinate-audit candidates (has_stage3b_coordinate_data)
     -> 1382 previous strict primary-QC subset (structures_primary_qc_subset.csv)
```

**1382 is not the master dataset.** Residue-level data for all 3159 measured
instance/models (771,542 residue rows, no coverage/resolution/R-free
filtering, e.g. VDR) is in `data/processed/lbd_residue_bfactors_all.parquet`.

Current release snapshot (see
[`data/manifests/final_data_infrastructure_release.json`](data/manifests/final_data_infrastructure_release.json)
for the authoritative, checksummed machine-readable version):

- 48 validated human nuclear receptor identities (`config/nr_metadata.py`).
- 2072 Stage 2 candidate receptor→polymer-entity relationships fully
  triaged, with lineage preserved for every one (selected or not).
- 1777 cached mmCIF structure files, SHA256-verified.
- **`structures.csv`: master inventory of all 2072 candidates** with
  coverage, resolution, R-free, method, identity, fusion and taxonomy as
  metadata/status flags (1840 with coordinate data).
- **1382-structure historical strict QC subset** (X-ray, >=90%
  positive-occupancy LBD coverage, 1.8-3.5 A, R-free <= 0.30,
  human/human-derived) preserved in `structures_primary_qc_subset.csv`;
  690 non-members with multi-reason lineage in `excluded_structures.csv`.
- Ligand state (APO/HOLO/AMBIGUOUS) and nonpolymer inventory for all 1840
  coordinate-assessable candidates (QC is metadata only); the 232 candidates
  without coordinate data are retained as `NOT_ASSESSABLE_NO_COORDINATE_DATA`.
  Coregulator/partner annotations cover the strict subset.
- `data/processed/lbd_residue_bfactors_all.parquet`: all 771,542 corrected
  residue observations — every standardized LBD position, observed or not —
  with a deterministic (never B-factor-based) raw Ca handoff record.
  `final_lbd_residue_map_primary_qc_subset.parquet` is the historical
  338,012-row strict-subset map.

## Pipeline (as implemented)

```
Professor-supplied NR metadata
        ↓
current database/nomenclature validation
        ↓
validated receptor master table
        ↓
RCSB structure discovery
        ↓
structure metadata retrieval
        ↓
standardized LBD definition + canonical residue mapping
        ↓
mmCIF download
        ↓
observed-coordinate completeness QC (positive-occupancy Ca model)
        ↓
master structure inventory (QC as flags) + deterministic representative-instance selection
        ↓
ligand-state + coregulator/partner annotation
        ↓
canonical residue-level handoff (data/processed/lbd_residue_bfactors_all.parquet)
        ↓
per-instance B-factor Z-normalization + chain→PDB-median aggregation
        ↓
structure-based canonical LBD segment mapping (v3: curated anchor templates
+ CE structural alignment + secondary-structure-refined transfer)
        ↓
segment-level feature construction (coverage/missingness policy)
        ↓
   ┌────────────────────┬─────────────────────┬───────────────────────┐
   ↓                    ↓                     ↓                       ↓
PCA/UMAP/clustering   subfamily/group      apo vs holo            functional-site
(no robust discrete   Kruskal-Wallis       flexibility            flexibility map
 clusters found)       + BH-FDR             comparison             + specificity
                       (0 significant)      (3 significant,        validation
                                             localized)             (matched-
                                                                     background
                                                                     confirmed)
   └────────────────────┴─────────────────────┴───────────────────────┘
        ↓
>>> CORE ANALYSIS COMPLETE -- see reports/FINAL_ANALYSIS_RESULTS.md <<<
        ↓
>>> ClinVar/disease mutation analysis (not yet performed) <<<
```

Each arrow above corresponds to a distinct, auditable stage with its own
scripts, tests, and manifests under `scripts/` and `data/manifests/`.

## Repository structure

```
biol363_nr_bfactor/
├── README.md                  This file
├── environment.yml            Human-readable conda environment spec
├── environment-lock.yml       Exact resolved-version reproducibility snapshot
├── .gitignore
│
├── config/
│   ├── README.md              Provenance rules for receptor metadata
│   └── nr_metadata.py         NR_METADATA placeholder (populated in Stage 1)
│
├── data/
│   ├── raw/
│   │   ├── mmcif/             Downloaded structure files (gitignored)
│   │   └── api/                Raw API responses (gitignored)
│   ├── interim/
│   │   ├── mappings/           Intermediate ID/chain/residue mappings (gitignored)
│   │   ├── qc/                 Quality-control intermediates (gitignored)
│   │   └── ligand_annotations/ Intermediate ligand-state annotations (gitignored)
│   ├── processed/               Analysis-ready datasets (gitignored)
│   └── manifests/               Provenance manifests (VERSION-CONTROLLED)
│
├── scripts/
│   ├── 00_validate_environment.py
│   └── utils/
│
├── notebooks/
├── logs/                        Run logs (gitignored)
├── reports/
│   ├── figures/
│   └── tables/
└── tests/
```

### Data directory semantics

| Directory              | Contents                                  | In Git? |
|-------------------------|--------------------------------------------|---------|
| `data/raw/`             | Unmodified downloads (mmCIF, API responses) | No |
| `data/interim/`         | Intermediate, regenerable working data      | No |
| `data/processed/`       | Final, analysis-ready team-handoff tables   | **The five handoff tables only** (`structures.csv`, `lbd_residue_bfactors_all.parquet`, `structures_primary_qc_subset.csv`, `final_lbd_residue_map_primary_qc_subset.parquet`, `excluded_structures.csv`) — everything else in this directory stays ignored |
| `data/manifests/`       | Provenance records (what was fetched, when, from where, with what parameters/hashes) | **Yes** |

Raw and most derived biological datasets are never committed to Git — they
are large, regenerable, and third-party-licensed. Manifest files are
committed because they document *provenance* (what was retrieved, from
where, when, and how) and are small, human-readable, and essential for
reproducibility and auditing. The three final team-handoff tables under
`data/processed/` are the one exception — small enough for normal Git and
needed as the actual deliverable (see [`TEAM_HANDOFF.md`](TEAM_HANDOFF.md)).

**Biological identifiers (UniProt accessions, PDB IDs, NR nomenclature
codes, etc.) must never be silently modified, inferred from filenames, or
"corrected" without an explicit, logged justification.** See
[`config/README.md`](config/README.md) for the full provenance policy.

## Environment setup

### Create the environment

```bash
conda env create -f environment.yml
```

### Activate the environment

```bash
conda activate biol363_nr
```

### Validate the environment

```bash
python scripts/00_validate_environment.py
```

This script checks Python version, required package availability and
versions, Parquet round-trip I/O, required directory structure, and that
`config/nr_metadata.py` imports successfully with the expected 48
UniProt-validated records (unique `uniprot_id` and `nr_code` values). It
performs **no network access**.

### Reproducing the exact environment

`environment.yml` is the human-readable specification used to create the
environment. `environment-lock.yml` is a full export of the *exact*
resolved package versions from a successful build, generated with:

```bash
conda env export -n biol363_nr --no-builds > environment-lock.yml
```

Use `environment.yml` to (re)create the environment; use
`environment-lock.yml` as the reproducibility record of exactly what was
installed.

## Testing

```bash
pytest -q
```

## Functional-site analysis

Two functional-site analyses exist. **They are separate layers and are
never merged**:

- **HISTORICAL ROLE D** (`scripts/31_analyze_functional_site_bfactors.py`):
  a hand-curated, 5-receptor descriptive analysis on the historical
  strict-QC subset (`structures_primary_qc_subset.csv`,
  `final_lbd_residue_map_primary_qc_subset.parquet`,
  [`functional_site_definitions.csv`](reports/tables/functional_site_definitions.csv)).
  Byte-identical since commit `e714bbd`; preserved for reproducibility, not
  re-run or modified by any later stage.
- **CURRENT (full-dataset)**, `scripts/72_ligand_pocket_geometry.py` through
  `scripts/92_specificity_summary.py`: a new analysis layer across the full
  42-receptor coordinate cohort, using the approved normalization/
  aggregation hierarchy and the current v3 structural-segment mapping, with
  a geometry-derived ligand-binding pocket (5.0 Å heavy-atom contact),
  literature/segment-correspondence coactivator-groove and dimerization-
  interface definitions, and a subsequent structural-context-matched-
  background specificity validation. See
  [`reports/tables/current_functional_site_summary.md`](reports/tables/current_functional_site_summary.md)
  and [`reports/tables/functional_site_specificity_summary.md`](reports/tables/functional_site_specificity_summary.md).

New functional-site work should build on the CURRENT layer, not the
historical Role D one.

## Status

**CORE SCIENTIFIC ANALYSIS COMPLETE**, through: Data & Infrastructure,
full-dataset normalization, v3 structural-segment mapping, segment-level
PCA/UMAP/clustering, subfamily/group statistical testing + flexibility
atlas, apo-vs-holo flexibility comparison, current full-dataset
functional-site flexibility map, and its specificity/confounding
validation. Every stage is tested (`pytest -q`) and deterministic.

Milestone git tags (`stage2-candidate-inventory-v1`,
`stage3a-lbd-prefilter-v1`, `stage3b-coordinate-audit-v1`) were cut in the
original development repository's history; a clean checkout of this
release repository does not carry them (a fresh clone's commit history is
separate from a different repository's tags). The frozen milestones those
tags represent are instead verified by content, directly against each
stage's release-summary JSON, in `tests/test_final_data_handoff.py`.

**Not yet performed**: ClinVar/disease mutation analysis. See
[`reports/FINAL_ANALYSIS_RESULTS.md`](reports/FINAL_ANALYSIS_RESULTS.md) for
the complete final numbers, [`reports/FINAL_REPRODUCIBILITY_CHECK.md`](reports/FINAL_REPRODUCIBILITY_CHECK.md)
for the authoritative script sequence, and
[`reports/FINAL_LIMITATIONS.md`](reports/FINAL_LIMITATIONS.md) before citing
any result.
