"""Shared, tested interval arithmetic for 1-based inclusive residue ranges.

Used by both the Stage 2 RCSB coverage computation (rcsb_extract.py) and
Stage 3A LBD/DBD mapping (scripts/12_map_candidates_to_lbd.py), so there is
exactly one implementation of "merge overlapping intervals" and "intersect
against a reference interval" in this project.

Convention: every interval is an inclusive (begin, end) pair, begin <= end,
both 1-based (i.e. length = end - begin + 1).
"""

from __future__ import annotations


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent (begin, end) intervals into the minimal
    sorted set of disjoint intervals covering the same union of positions."""
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for beg, end in ordered[1:]:
        last_beg, last_end = merged[-1]
        if beg <= last_end + 1:
            merged[-1] = (last_beg, max(last_end, end))
        else:
            merged.append((beg, end))
    return merged


def total_length(intervals: list[tuple[int, int]]) -> int:
    """Sum of lengths of a set of (possibly overlapping) intervals, counting
    each covered position once (i.e. length of the merged union)."""
    return sum(end - beg + 1 for beg, end in merge_intervals(intervals))


def intersect_with_reference(
    intervals: list[tuple[int, int]], ref_start: int, ref_end: int
) -> list[tuple[int, int]]:
    """Intersect a (merged) set of intervals against a single fixed
    reference interval [ref_start, ref_end], returning the overlapping
    sub-intervals (still expressed in the original coordinate system)."""
    merged = merge_intervals(intervals)
    overlaps = []
    for beg, end in merged:
        lo, hi = max(beg, ref_start), min(end, ref_end)
        if hi >= lo:
            overlaps.append((lo, hi))
    return overlaps


def overlap_fraction(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> float:
    """Fraction of the SHORTER interval set's total span covered by overlap
    with the other (used for the Stage 2 multi-UniProt classifier)."""
    a_merged, b_merged = merge_intervals(a), merge_intervals(b)
    if not a_merged or not b_merged:
        return 0.0
    overlap = 0
    for a_beg, a_end in a_merged:
        for b_beg, b_end in b_merged:
            lo, hi = max(a_beg, b_beg), min(a_end, b_end)
            if hi >= lo:
                overlap += hi - lo + 1
    span_a = total_length(a_merged)
    span_b = total_length(b_merged)
    shorter_span = min(span_a, span_b)
    return overlap / shorter_span if shorter_span else 0.0
