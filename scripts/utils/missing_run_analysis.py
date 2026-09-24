"""Decompose per-residue standardized-LBD observation status into
contiguous "missing runs", distinguishing:

- unmapped runs: standardized LBD positions the construct's sequence
  mapping never reached (is_construct_mapped == False) — a construct
  truncation / expression-boundary effect.
- mapped-but-unmodeled runs: positions the construct covers, but with no
  positive-occupancy Cα (is_construct_mapped == True, ca_positive_occupancy
  == False) — an unresolved-coordinate effect.

Each mapped-but-unmodeled run is further classified terminal/internal:
terminal if it touches lbd_relative_position 1 or lbd_length (the standard
LBD's own N-/C-terminal edge), internal otherwise (flanked by observed
residues on both sides). This distinguishes ordinary N-/C-terminal
disorder from a genuinely internal disordered loop.

Operates purely on an already-built per-residue observation group (e.g. one
(pdb_id, polymer_entity_id, instance_id, model_num) slice of
data/interim/qc/stage3b_lbd_ca_observations.parquet) — no mmCIF re-parse.
"""

from __future__ import annotations

import pandas as pd


def _find_runs(flags: list[bool]) -> list[tuple[int, int]]:
    """Return (start_idx, end_idx) inclusive 0-based index pairs for each
    maximal contiguous run where flags[i] is True."""
    runs = []
    start = None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        elif not f and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(flags) - 1))
    return runs


def analyze_missing_runs(group: pd.DataFrame, lbd_length: int) -> dict:
    """`group` must be one instance/model's rows, one per
    lbd_relative_position from 1..lbd_length (as stored in the Stage 3B
    observations parquet)."""
    g = group.sort_values("lbd_relative_position")
    assert len(g) == lbd_length, f"expected {lbd_length} rows, got {len(g)}"
    assert list(g["lbd_relative_position"]) == list(range(1, lbd_length + 1))

    is_mapped = list(g["is_construct_mapped"])
    is_observed = list(g["ca_positive_occupancy"])

    unmapped_flags = [not m for m in is_mapped]
    unmapped_runs = _find_runs(unmapped_flags)

    mapped_unmodeled_flags = [m and not o for m, o in zip(is_mapped, is_observed)]
    mapped_unmodeled_runs = _find_runs(mapped_unmodeled_flags)

    internal_runs = []
    terminal_runs = []
    for start, end in mapped_unmodeled_runs:
        touches_edge = (start == 0) or (end == lbd_length - 1)
        (terminal_runs if touches_edge else internal_runs).append((start, end))

    internal_residue_count = sum(end - start + 1 for start, end in internal_runs)
    terminal_residue_count = sum(end - start + 1 for start, end in terminal_runs)
    longest_internal_run = max((end - start + 1 for start, end in internal_runs), default=0)
    longest_mapped_unmodeled_run = max((end - start + 1 for start, end in mapped_unmodeled_runs), default=0)

    return {
        "unmapped_run_count": len(unmapped_runs),
        "unmapped_residue_count": sum(end - start + 1 for start, end in unmapped_runs),
        "mapped_but_unmodeled_run_count": len(mapped_unmodeled_runs),
        "mapped_but_unmodeled_residue_count": sum(end - start + 1 for start, end in mapped_unmodeled_runs),
        "internal_missing_run_count": len(internal_runs),
        "internal_missing_residue_count": internal_residue_count,
        "terminal_missing_run_count": len(terminal_runs),
        "terminal_missing_residue_count": terminal_residue_count,
        "longest_internal_missing_run": longest_internal_run,
        "longest_mapped_unmodeled_run": longest_mapped_unmodeled_run,
    }
