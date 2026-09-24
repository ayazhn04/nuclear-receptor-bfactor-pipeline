# Stage 3B Data Recovery Note

## What happened

During Stage 3B.1, an intermediate draft of the metrics-recomputation script
(`scripts/22_recompute_stage3b_metrics.py`) accidentally overwrote
`data/manifests/mmcif_unobserved_residue_crosscheck.csv` — a DERIVATIVE
manifest (computed from mmCIF content, not raw source data) — with a buggy,
empty (0-row) recomputation, before the bug was caught by re-running the
project's test suite.

## What was NOT affected

The **1777 authoritative cached mmCIF files** (`data/raw/mmcif/*.cif.gz`)
were **never touched** by this bug. They remained exactly as originally
downloaded from `https://files.rcsb.org/`, verified by SHA256 checksum
before and after the incident.

## How it was recovered

The overwritten derivative was rebuilt **locally, directly from the intact
mmCIF cache** (`scripts/utils/rebuild_unobserved_ground_truth.py`), by
re-extracting the `_pdbx_unobs_or_zero_occ_residues` mmCIF category for all
1777 selected entries. **No network request was made** — this is a local
decompression + parse of files already on disk.

## Verification

The regenerated crosscheck counts reproduced the previously observed Stage 3B
results exactly:

| Category | Count |
|---|---|
| A (deposited-unobserved, no Cα record) | 32,756 |
| B (deposited-unobserved, zero-occupancy-only Cα) | 270 |
| C (deposited-unobserved, positive-occupancy Cα) | 0 |
| D (not deposited-unobserved, no Cα record) | 13 |

Stage 3B.2 additionally extracted 73545 total deposited-unobserved
records from the same cache (saved this time to
`data/manifests/mmcif_unobserved_residue_ground_truth.csv` for durable
provenance) and used them to complete the residual-missingness audit
(category E: zero-occupancy-only, not listed as unobserved — 11 cases).

Subsequent offline tests and two full deterministic reruns of the
occupancy-correction and metrics pipeline (Stage 3B.1) produced
byte-identical outputs, confirming the recovery was complete and correct.
