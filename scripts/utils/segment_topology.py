"""
Structural-segment alignment, Part 4 core logic: turn one receptor's
per-residue secondary-structure consensus into a canonical H1-H12 / beta-sheet
/ loop segment assignment.

This is deliberately conservative. The nuclear-receptor LBD fold (a
three-layered alpha-helical sandwich of up to ~12 helices plus a small
beta-hairpin, with H12 as the C-terminal AF-2 activation helix immediately
preceding the domain's own C-terminus, and H1 as the first helix after the
domain's own N-terminus) is well-established, superfamily-wide structural
biology (Wurtz et al. 1996; Renaud & Moras 2000; Bourguet et al. 1995) -- not
something derived from a new alignment or new software. It is used only to
place two INDEPENDENT anchors (H1, H12) plus the beta-sheet, and to attempt
sequential interior labels ONLY where the observed interior helix count
between two already-anchored landmarks exactly matches the canonical
expectation. Everywhere else, the assignment is left AMBIGUOUS/UNASSIGNED
rather than guessed -- the nth detected helix is never assumed to be Hn.
"""

from __future__ import annotations

from dataclasses import dataclass, field

CANONICAL_INTERIOR_BEFORE_SHEET = ["H2", "H3", "H4", "H5"]   # between H1 and beta_sheet
CANONICAL_INTERIOR_AFTER_SHEET = ["H6", "H7", "H8", "H9", "H10", "H11"]  # between beta_sheet and H12
CANONICAL_INTERIOR_NO_SHEET = ["H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11"]

TERMINUS_TOLERANCE_HIGH = 5
TERMINUS_TOLERANCE_MODERATE = 15
SHEET_MERGE_GAP = 5


@dataclass
class Run:
    state: str  # HELIX / SHEET
    start: int
    end: int


@dataclass
class SegmentCall:
    start: int
    end: int
    structural_segment: str  # e.g. 'H1', 'beta_sheet', '' if unresolved
    structural_class: str  # HELIX / SHEET / LOOP / UNRESOLVED_SSE
    assignment_method: str
    assignment_confidence: str  # HIGH / MODERATE / LOW / NONE
    assignment_status: str
    evidence_notes: str


def build_runs(positions: list[int], states: list[str]) -> list[Run]:
    runs: list[Run] = []
    for pos, st in sorted(zip(positions, states)):
        if runs and runs[-1].state == st and pos == runs[-1].end + 1:
            runs[-1].end = pos
        else:
            runs.append(Run(st, pos, pos))
    return runs


def merge_sheet_runs(runs: list[Run]) -> list[Run]:
    """Merge SHEET runs separated only by COIL gaps <= SHEET_MERGE_GAP (a beta
    hairpin's two strands are typically joined by a short turn), never across
    an intervening HELIX run."""
    out: list[Run] = []
    i = 0
    while i < len(runs):
        r = runs[i]
        if r.state != "SHEET":
            out.append(r)
            i += 1
            continue
        j = i
        merged = Run("SHEET", r.start, r.end)
        while j + 1 < len(runs) and runs[j + 1].state == "SHEET" and (runs[j + 1].start - merged.end - 1) <= SHEET_MERGE_GAP:
            merged.end = runs[j + 1].end
            j += 1
        out.append(merged)
        i = j + 1
    return out


def assign_receptor_segments(
    positions: list[int], states: list[str], lbd_start: int, lbd_end: int,
    functional_terminal_helix_ranges: list[tuple[int, int]] | None = None,
) -> list[SegmentCall]:
    """positions/states: one entry per residue with secondary-structure
    evidence (COIL residues included) for one receptor's standardized LBD.
    Returns SegmentCall entries covering every HELIX/SHEET run (loops are
    added separately by the caller, which also knows the full LBD position
    range including residues with no secondary-structure evidence at all)."""
    functional_terminal_helix_ranges = functional_terminal_helix_ranges or []
    raw_runs = build_runs(positions, states)
    runs = [r for r in raw_runs if r.state in ("HELIX", "SHEET")]
    runs = merge_sheet_runs(runs)
    helix_runs = [r for r in runs if r.state == "HELIX"]
    sheet_runs = [r for r in runs if r.state == "SHEET"]

    calls: list[SegmentCall] = []
    if not helix_runs:
        return calls  # caller marks the whole receptor UNASSIGNED (no resolvable helices)

    # ---- H1 anchor: first helix, near the domain's own N-terminus ----
    h1 = helix_runs[0]
    offset = h1.start - lbd_start
    if 0 <= offset <= TERMINUS_TOLERANCE_HIGH:
        conf, status, note = "HIGH", "ASSIGNED_HIGH_CONFIDENCE", f"First detected helix starts {offset} residue(s) after lbd_start ({lbd_start})."
    elif 0 <= offset <= TERMINUS_TOLERANCE_MODERATE:
        conf, status, note = "MODERATE", "ASSIGNED_MODERATE_CONFIDENCE", f"First detected helix starts {offset} residues after lbd_start ({lbd_start})."
    else:
        conf, status, note = "NONE", "AMBIGUOUS", f"First detected helix starts {offset} residues after lbd_start ({lbd_start}) -- too far to confidently anchor as H1."
    calls.append(SegmentCall(h1.start, h1.end, "H1" if status != "AMBIGUOUS" else "", "HELIX", "topology_terminus_anchor", conf, status, note))

    # ---- H12 anchor: last helix, near the domain's own C-terminus ----
    h12 = helix_runs[-1]
    offset = lbd_end - h12.end
    cross_check = any(a <= h12.start <= b or a <= h12.end <= b for a, b in functional_terminal_helix_ranges)
    if len(helix_runs) > 1 and 0 <= offset <= TERMINUS_TOLERANCE_HIGH:
        conf, status = "HIGH", "ASSIGNED_HIGH_CONFIDENCE"
        note = f"Last detected helix ends {offset} residue(s) before lbd_end ({lbd_end})."
        if cross_check:
            note += " Overlaps a repository-documented terminal-helix/AF-2 functional feature."
    elif len(helix_runs) > 1 and 0 <= offset <= TERMINUS_TOLERANCE_MODERATE:
        conf, status = "MODERATE", "ASSIGNED_MODERATE_CONFIDENCE"
        note = f"Last detected helix ends {offset} residues before lbd_end ({lbd_end})."
        if cross_check:
            conf, status = "HIGH", "ASSIGNED_HIGH_CONFIDENCE"
            note += " Upgraded to high confidence: overlaps a repository-documented terminal-helix/AF-2 functional feature."
    else:
        conf, status = "NONE", "AMBIGUOUS"
        note = f"Last detected helix ends {offset} residues before lbd_end ({lbd_end}) -- too far to confidently anchor as H12."
    calls.append(SegmentCall(h12.start, h12.end, "H12" if status != "AMBIGUOUS" else "", "HELIX", "topology_terminus_anchor", conf, status, note))

    interior_helices = helix_runs[1:-1] if len(helix_runs) > 2 else []

    # ---- beta_sheet ----
    sheet_before_idx = None
    if len(sheet_runs) == 1:
        sh = sheet_runs[0]
        n_before = sum(1 for h in helix_runs if h.end < sh.start)
        n_after = sum(1 for h in helix_runs if h.start > sh.end)
        if n_before > 0 and n_after > 0:
            calls.append(SegmentCall(sh.start, sh.end, "beta_sheet", "SHEET", "topology_unique_sheet",
                                      "HIGH", "ASSIGNED_HIGH_CONFIDENCE",
                                      f"Single consensus beta-sheet region, {n_before} helix run(s) before it and {n_after} after."))
            sheet_before_idx = n_before
    elif len(sheet_runs) > 1:
        for sh in sheet_runs:
            calls.append(SegmentCall(sh.start, sh.end, "", "SHEET", "topology_multiple_sheets",
                                      "NONE", "AMBIGUOUS",
                                      f"{len(sheet_runs)} disjoint consensus beta-sheet regions detected (not merged: separated by a helix or a gap > {SHEET_MERGE_GAP} residues); "
                                      "cannot confidently identify the canonical beta-hairpin without further evidence."))

    def label_block(block: list[Run], expected_names: list[str], context: str) -> None:
        if len(block) == len(expected_names):
            for helix, name in zip(block, expected_names):
                calls.append(SegmentCall(
                    helix.start, helix.end, name, "HELIX", "topology_ordinal_bracketed",
                    "MODERATE", "ASSIGNED_MODERATE_CONFIDENCE",
                    f"{context}: {len(block)} interior helices detected, matching the canonical count -> assigned by order.",
                ))
        else:
            for helix in block:
                calls.append(SegmentCall(
                    helix.start, helix.end, "", "HELIX", "topology_ordinal_bracketed",
                    "NONE", "AMBIGUOUS",
                    f"{context}: {len(block)} interior helices detected, expected {len(expected_names)} -- "
                    "count mismatch (extra/insertion or missing helix); not assigned to a specific H-number.",
                ))

    if sheet_before_idx is not None:
        before = interior_helices[:sheet_before_idx - 1] if sheet_before_idx >= 1 else []
        # interior_helices excludes h1 already; sheet_before_idx counts helices before sheet INCLUDING h1
        before = interior_helices[: max(sheet_before_idx - 1, 0)]
        after = interior_helices[max(sheet_before_idx - 1, 0):]
        label_block(before, CANONICAL_INTERIOR_BEFORE_SHEET, "Helices between H1 and beta_sheet")
        label_block(after, CANONICAL_INTERIOR_AFTER_SHEET, "Helices between beta_sheet and H12")
    elif interior_helices:
        label_block(interior_helices, CANONICAL_INTERIOR_NO_SHEET, "Helices between H1 and H12 (no confidently placed beta-sheet)")

    return calls
