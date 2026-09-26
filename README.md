# Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families

**Course:** BIOL363, Nazarbayev University

## Current repository purpose

**The Data & Infrastructure pipeline is COMPLETE.** This repository
provides a fully reproducible, auditable, version-controlled scaffold,
receptor identity table, structure discovery/selection pipeline, and final
team-handoff dataset for downstream B-factor/flexibility analysis.
**Downstream B-factor normalization, clustering, and statistical analysis
have NOT been performed anywhere in this repository** — see
[`TEAM_HANDOFF.md`](TEAM_HANDOFF.md) to start that work.

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
>>> handoff to downstream B-factor normalization/analysis (not yet performed) <<<
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

Role D's descriptive functional-site analysis uses the primary handoff
tables and the literature-curated definitions in
[`reports/tables/functional_site_definitions.csv`](reports/tables/functional_site_definitions.csv).
It normalizes usable C-alpha B-factors to Z-scores within each selected
receptor instance, reduces multiple selected entities to a median profile per
PDB entry, then reports a median complete-site score per receptor. PDB entries
with any missing predefined site residue retain their coverage QC but are not
included in the receptor-level site median.

```bash
python scripts/31_analyze_functional_site_bfactors.py
```

The command writes normalization QC, PDB-entry site scores, and
receptor-level site summaries to `reports/tables/`. These are descriptive
results only. APO/HOLO inference, pharmacology labels, and cross-receptor
statistical tests require separately specified comparison cohorts and are not
performed by this script.

## Status

**Data & Infrastructure pipeline COMPLETE.** All stages listed above,
through the final structure selection, representative-instance selection,
ligand-state annotation, and canonical-residue handoff, are done, tested,
and frozen at annotated Git tags (`stage2-candidate-inventory-v1`,
`stage3a-lbd-prefilter-v1`, `stage3b-coordinate-audit-v1`,
`data-infrastructure-v1`; the master-inventory revision is a later commit).

**Downstream B-factor normalization, clustering, PCA/t-SNE/UMAP,
statistical testing, functional-site correlation analysis, and any
biological conclusions have NOT been performed anywhere in this
repository.** Start that work from [`TEAM_HANDOFF.md`](TEAM_HANDOFF.md).
