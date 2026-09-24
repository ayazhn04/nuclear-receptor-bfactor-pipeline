# Stage 2 Schema Hygiene Note — `rcsb_candidate_inventory.csv`

## What was found

While building Stage 3A, we found that the frozen, tagged Stage 2 artifact

```
data/manifests/rcsb_candidate_inventory.csv
```

(as committed at `ffd9472` / tag `stage2-candidate-inventory-v1`) contains
**redundant, `_x`/`_y`-suffixed duplicate columns** alongside the canonical
ones:

- `computed_reference_sequence_coverage` (+ `_x`, `_y`)
- `computed_entity_sequence_coverage` (+ `_x`, `_y`)
- `rcsb_reference_sequence_coverage` (+ `_x`, `_y`)
- `rcsb_entity_sequence_coverage` (+ `_x`, `_y`)

## Root cause

During Stage 2.2/2.3, `scripts/09_finalize_stage2_coverage.py` was rerun
several times while an idempotency bug was being fixed (a rerun would try
to re-merge coverage columns it had already written, without first
dropping its own prior output). One of those pre-fix reruns wrote a
CSV containing `_x`/`_y`-suffixed columns from a pandas merge collision.
After the bug was fixed, a subsequent run correctly regenerated the
clean, unsuffixed columns — but the fix's drop-list only named the
canonical column names, not the already-present `_x`/`_y` leftovers, so
those stale duplicates were never removed before the file was committed.

## Verified: values are identical, no scientific content is affected

For all 2072 rows, `computed_reference_sequence_coverage_x`,
`computed_reference_sequence_coverage_y`, and
`computed_reference_sequence_coverage` are **byte-for-byte identical**
(and likewise for the other three column families). This was verified
programmatically before writing this note. No coverage value, mapping
result, audit status, or any other scientific conclusion from Stage 2 is
affected by this defect — it is a redundant-column hygiene issue only.

## Which columns downstream scripts use

All Stage 3A scripts (and this note) read only the **canonical, unsuffixed**
column names:

- `computed_reference_sequence_coverage`
- `computed_entity_sequence_coverage`
- `rcsb_reference_sequence_coverage`
- `rcsb_entity_sequence_coverage`

The `_x`/`_y` variants are ignored entirely and should be treated as dead
weight, not as a second data source.

## Why the tagged artifact was not rewritten

`stage2-candidate-inventory-v1` is the frozen, reviewed, approved snapshot
of the Stage 2 discovery layer. Rewriting a tagged commit's file contents
after the fact — even for a purely cosmetic fix — would undermine the
reproducibility guarantee that tag is meant to provide (anyone checking out
that tag should see exactly what was reviewed). Per explicit instruction,
this defect is **documented here rather than silently corrected or
retroactively rewritten**. A follow-up hygiene commit (dropping the
duplicate columns, values unchanged) may be applied on top of the current
`main` branch state whenever convenient, without touching history.
