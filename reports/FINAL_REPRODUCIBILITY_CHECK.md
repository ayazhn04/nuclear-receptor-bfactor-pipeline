# Final Reproducibility Check

> Verifies the complete CURRENT-line supported pipeline reproduces deterministically from the already-cached local data (`data/raw/mmcif/`, `data/raw/api/`; no new network access anywhere in this project). Historical/deprecated v1/v2 structural-mapping scripts are documented but explicitly **not** rerun over current outputs, per instruction -- they remain frozen in git history at their own commits and are never re-executed to "refresh" a v3-superseded file.

## Authoritative script sequence (current line only)

Scripts not listed here (Stage 0-2 discovery/acquisition, Stage 3A/3B coordinate QC, `scripts/24`-`30` final release, `scripts/31` Role D, `scripts/38`-`43` v1/v2 structural mapping) are earlier/historical stages, already frozen, and are not part of this reproducibility check -- see `reports/TABLE_INVENTORY.md` for their HISTORICAL_ONLY status.

| # | Script | Stage |
|---|---|---|
| 32 | `32_build_normalization_foundation.py` | Per-instance Z-normalization + chain→PDB aggregation |
| 33 | `33_build_analysis_foundation_reports.py` | Normalization QC reports |
| 34 | `34_build_dssp_consensus.py` | Secondary-structure consensus (input to v3 mapping) |
| 44-47 | `44_curate_templates_v3.py` → `47_validate_v3.py` | **Current (v3)** structural-segment mapping |
| 48-54 | `48_feature_support_audit.py` → `54_modeling_summary.py` | PCA/UMAP/clustering |
| 55-62 | `55_subfamily_group_definition.py` → `62_subfamily_summary.py` | Subfamily testing + flexibility atlas |
| 63-71 | `63_apo_holo_eligibility.py` → `71_apo_holo_summary.py` | Apo vs holo comparison |
| 72-82 | `72_ligand_pocket_geometry.py` → `82_functional_site_summary.py` | Current functional-site flexibility map |
| 83-92 | `83_functional_site_fdr_audit.py` → `92_specificity_summary.py` | Specificity/confounding validation |
| 93-95 | `93_final_receptor_identity_audit.py` → `95_final_statistical_consistency_audit.py` | **This freeze stage** |
| 30 | `30_finalize_release.py` | Release manifest (rerun last, after any documentation update, to keep checksums current) |

## Verification performed

1. **Environment validator**: `python scripts/00_validate_environment.py` → **PASS** (package versions, Parquet round-trip, 15 required directories, 48-record `NR_METADATA` import). No network access used.
2. **Deterministic outputs**: every stage listed above (32-95) has its own `test_..._pipeline_deterministic_and_idempotent` (or equivalently named) test in `tests/`, which reruns that stage's scripts fresh and asserts byte-identical SHA256 hashes on every key output file. All such tests are part of the full `pytest -q` run below and all pass.
3. **All final hashes stable**: confirmed both via the per-stage determinism tests (item 2) and via this freeze stage's own hash-frozen upstream-input checks (`data/analysis/modeling_frozen_input_hashes.json`, re-verified by every stage from the modeling stage onward).
4. **Tests pass**: `pytest -q` → see the exact result in the commit message / final report for this stage (0 failures after the two fixes in Part 9: the stale `release.json` checksum was refreshed through its own intended generator, `scripts/30_finalize_release.py`; the milestone-tag test was adapted to this clean release repository's actual provenance model, matching the precedent already established for other development-only assumptions in `tests/test_master_inventory.py`).
5. **No source/master membership changes**: `data/processed/structures.csv` remains 2072 rows throughout (checked by a dedicated test in nearly every stage's test file, and independently reconfirmed here). No candidate has ever been removed from the master inventory by any stage in this project.

## What "reproducibility" means in this project

A literal from-scratch rerun starting from zero (deleting `data/raw/` and re-downloading from RCSB/UniProt) is **not** how this project's reproducibility is defined or verified, and was not performed here: `data/raw/` is gitignored by design (large, regenerable, third-party-licensed), and every stage's scripts and documentation consistently state "no new network access" as an explicit design constraint. Reproducibility here means: **given the already-cached raw inputs, every current-line script produces byte-identical output across repeated runs**, which is exactly what the per-stage determinism tests verify, and what was reconfirmed for this freeze stage's own new scripts (93-95) during this same check (a genuine non-determinism was found and fixed here -- see the freeze-stage commit message).
