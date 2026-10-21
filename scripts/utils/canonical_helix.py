"""
Professor-approved conventional NR LBD helix map: explicit policy, identity
crosswalk and interval expansion helpers (used by
scripts/97_build_canonical_helix_map.py and tests/test_canonical_helix_map.py).

TWO POLICIES live here, deliberately separate:

  FINAL (professor decision of 2026-10-20; used for every final analysis):
    the professor table is ground truth. Every professor label is used,
    including `review_required` (status kept as metadata only). A position
    inside two professor intervals carries both labels. TRalpha is numbered on
    the professor's accession P10827-2 (TRalpha1). See the FINAL section below.

  PRE-DECISION (audit history only; first remap pass, 2026-10-17): the rules
  described next. Kept so the audit-history outputs stay reproducible; never
  used for final labels.

Rules encoded here (see reports/CANONICAL_HELIX_MAP_README.md):
  * Every source interval keeps its professor-assigned `canonical_helix`
    name verbatim (conventional NR nomenclature, incl. H2', H3n, H11',
    H10/H11). Names are never renumbered.
  * `assignment_status` drives an explicit label policy. `review_required`
    intervals stay visibly uncertain: their name is kept as a *proposed*
    label but the analysis label is withheld.
  * Positions covered by no interval are NOT_ASSIGNED_IN_SOURCE -- never
    "loop", never biological absence. Positions covered by two different
    elements are CONFLICT and keep both names; nothing is resolved silently.
  * Source boundaries are inclusive UniProt positions in the stated
    `sequence_accession`. Joining to project tables (keyed by the project
    UniProt accession) goes through an explicit crosswalk; positions with no
    verified counterpart are never shifted or forced in.
"""

from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Source schema
# ---------------------------------------------------------------------------

SOURCE_COLUMNS = [
    "receptor", "nr_code", "sequence_accession", "canonical_helix", "start_uniprot", "end_uniprot",
    "reference_pdb", "reference_chain", "assignment_status", "boundary_basis", "source_url", "note",
    "aligned_support_residues", "mean_aligned_distance_A",
]
SOURCE_ENCODING = "utf-8-sig"  # the supplied file carries a UTF-8 byte-order mark

# ---------------------------------------------------------------------------
# Assignment-status policy
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StatusPolicy:
    evidence_category: str
    label_policy: str
    label_usable_for_analysis: bool
    review_flag: bool


STATUS_POLICY: dict[str, StatusPolicy] = {
    "published_boundaries": StatusPolicy(
        "PUBLISHED_RECEPTOR_SPECIFIC_BOUNDARY", "LABEL_ASSIGNED", True, False),
    "structural_homology": StatusPolicy(
        "STRUCTURAL_HOMOLOGY_REFERENCE_ENDPOINTS", "LABEL_ASSIGNED", True, False),
    "review_required": StatusPolicy(
        "PROPOSAL_REVIEW_REQUIRED", "LABEL_WITHHELD_PENDING_REVIEW", False, True),
}


def status_policy(assignment_status: str) -> StatusPolicy:
    """Policy for one source status. Unknown statuses are an error, never a
    silent default (a new status must be added here deliberately)."""
    try:
        return STATUS_POLICY[assignment_status]
    except KeyError:
        raise ValueError(f"Unknown assignment_status {assignment_status!r}; "
                         f"known: {sorted(STATUS_POLICY)}") from None


# ---------------------------------------------------------------------------
# Element vocabulary (conventional NR nomenclature, N->C order)
# ---------------------------------------------------------------------------

ELEMENT_ORDER = ["H1", "H2", "H2'", "H3n", "H3", "H4", "H5", "H6", "H7", "H8", "H9",
                 "H10", "H10/H11", "H11", "H11'", "H12"]
ELEMENT_RANK = {e: i for i, e in enumerate(ELEMENT_ORDER)}
RECEPTOR_SPECIFIC_ELEMENTS = {"H2'", "H3n", "H11'"}


def element_rank(name: str) -> int:
    if name not in ELEMENT_RANK:
        raise ValueError(f"Unknown canonical element {name!r}; known: {ELEMENT_ORDER}")
    return ELEMENT_RANK[name]


# ---------------------------------------------------------------------------
# Position status vocabulary for the residue-level map
# ---------------------------------------------------------------------------

POS_ASSIGNED = "ASSIGNED"
POS_REVIEW = "REVIEW_REQUIRED"
POS_CONFLICT = "CONFLICT_MULTIPLE_ELEMENTS"
POS_GAP = "NOT_ASSIGNED_IN_SOURCE"
POS_NO_RECEPTOR = "RECEPTOR_NOT_IN_SOURCE"
POS_ISOFORM = "PROJECT_ISOFORM_REGION_NOT_COVERED_BY_SOURCE"
POSITION_STATUSES = [POS_ASSIGNED, POS_REVIEW, POS_CONFLICT, POS_GAP, POS_NO_RECEPTOR, POS_ISOFORM]

# ---------------------------------------------------------------------------
# Identity crosswalk: source sequence_accession -> project uniprot_id
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Crosswalk:
    source_accession: str
    project_uniprot_id: str
    # Inclusive source-position range with a verified 1:1 identical-numbering
    # counterpart in the project accession; None = whole sequence identical
    # accession (no isoform difference).
    identical_range: tuple[int, int] | None
    basis: str


# TRalpha: the professor table uses TRalpha1 (P10827-2, 410 aa); the project
# uses the canonical UniProt entry P10827 (TRalpha2, 490 aa). The two isoform
# sequences are identical over residues 1-370 and diverge from 371 (verified
# against the UniProt FASTA files stored in
# data/sources/professor_canonical_helix_table/uniprot_isoform_sequences/ and
# re-checked by the tests). Positions 371-410 of P10827-2 therefore have NO
# counterpart in project numbering and are never joined.
TRA_SOURCE = "P10827-2"
TRA_PROJECT = "P10827"
TRA_IDENTICAL_END = 370

ISOFORM_CROSSWALKS = {
    TRA_SOURCE: Crosswalk(
        TRA_SOURCE, TRA_PROJECT, (1, TRA_IDENTICAL_END),
        "TRalpha1 (P10827-2) and canonical TRalpha2 (P10827) share identical sequence and numbering "
        "for residues 1-370; P10827-2 residues 371-410 (TRalpha1-specific C-terminus) have no counterpart "
        "in P10827 numbering and are not joined."),
}


def crosswalk_for(source_accession: str, project_uniprot_id: str) -> Crosswalk:
    """Explicit crosswalk. Same accession -> identity. A known isoform ->
    its documented crosswalk. Anything else is an error, never a guess."""
    if source_accession == project_uniprot_id:
        return Crosswalk(source_accession, project_uniprot_id, None, "identical accession")
    cw = ISOFORM_CROSSWALKS.get(source_accession)
    if cw is not None and cw.project_uniprot_id == project_uniprot_id:
        return cw
    raise ValueError(f"No documented crosswalk from {source_accession!r} to project accession "
                     f"{project_uniprot_id!r}")


def to_project_position(cw: Crosswalk, source_position: int) -> int | None:
    """Project position for one source position, or None when it has no
    verified counterpart. Numbering is never shifted."""
    if cw.identical_range is None:
        return source_position
    lo, hi = cw.identical_range
    return source_position if lo <= source_position <= hi else None


def first_sequence_difference(a: str, b: str) -> int | None:
    """1-based position of the first residue that differs (or where the
    shorter sequence ends); None if identical."""
    for i, (x, y) in enumerate(zip(a, b), start=1):
        if x != y:
            return i
    return None if len(a) == len(b) else min(len(a), len(b)) + 1


def read_fasta_sequence(path) -> str:
    with open(path, encoding="utf-8") as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))


# ---------------------------------------------------------------------------
# Interval expansion / conflicts
# ---------------------------------------------------------------------------


def expand_interval(start: int, end: int) -> list[int]:
    """Inclusive expansion. start > end is an error (reported upstream)."""
    if start > end:
        raise ValueError(f"start {start} > end {end}")
    return list(range(start, end + 1))


def position_assignments(intervals: list[dict]) -> dict[int, list[dict]]:
    """position -> list of covering intervals (each dict needs start, end).
    Kept as a list so overlaps are visible, never collapsed."""
    out: dict[int, list[dict]] = {}
    for iv in intervals:
        for p in expand_interval(iv["start"], iv["end"]):
            out.setdefault(p, []).append(iv)
    return out


def classify_position(covering: list[dict]) -> str:
    """Position status from the intervals covering it."""
    if not covering:
        return POS_GAP
    names = {iv["canonical_helix"] for iv in covering}
    if len(names) > 1:
        return POS_CONFLICT
    if all(iv["assignment_status"] == "review_required" for iv in covering):
        return POS_REVIEW
    if any(iv["assignment_status"] == "review_required" for iv in covering):
        # Same element covered by both a confirmed and a review interval:
        # stays visibly uncertain rather than being promoted.
        return POS_REVIEW
    return POS_ASSIGNED


# ===========================================================================
# FINAL POLICY (professor decision 2026-10-20) -- used for all final analyses
# ===========================================================================

FINAL_POLICY_VERSION = "professor_final_2026-10-20"
FINAL_LABEL_POLICY = "PROFESSOR_LABEL_FINAL"

# Residue-numbering reference sequence per project receptor. The project keeps
# its receptor identity key (uniprot_id); only the numbering reference differs.
# TRalpha: the professor table is on TRalpha1 (P10827-2), so TRalpha structures
# (TRalpha1 constructs) are numbered on P10827-2.
REFERENCE_ACCESSION = {TRA_PROJECT: TRA_SOURCE}


def reference_accession(project_uniprot_id: str) -> str:
    return REFERENCE_ACCESSION.get(project_uniprot_id, project_uniprot_id)


def final_crosswalk(source_accession: str, project_uniprot_id: str) -> Crosswalk:
    """Under the final policy the project numbering reference IS the professor's
    accession, so the crosswalk is the identity. Anything else is an error."""
    ref = reference_accession(project_uniprot_id)
    if source_accession != ref:
        raise ValueError(f"{project_uniprot_id}: professor accession {source_accession!r} != final "
                         f"numbering reference {ref!r}")
    return Crosswalk(source_accession, project_uniprot_id, None,
                     "identity: project numbering reference is the professor accession"
                     + (" (TRalpha1 isoform P10827-2)" if ref != project_uniprot_id else ""))


FPOS_ELEMENT = "PROFESSOR_ELEMENT"
FPOS_MULTI = "PROFESSOR_MULTIPLE_ELEMENTS"        # inside two professor intervals: both labels apply
FPOS_INTER = "INTER_HELIX_UNASSIGNED_DERIVED"     # inside the analysis span, outside every professor interval
FPOS_NO_MAPPING = "CANONICAL_MAPPING_UNAVAILABLE"  # receptor has no professor mapping
FINAL_POSITION_STATUSES = [FPOS_ELEMENT, FPOS_MULTI, FPOS_INTER, FPOS_NO_MAPPING]


def classify_position_final(covering: list[dict], receptor_mapped: bool = True) -> str:
    if not receptor_mapped:
        return FPOS_NO_MAPPING
    if not covering:
        return FPOS_INTER
    return FPOS_MULTI if len({iv["canonical_helix"] for iv in covering}) > 1 else FPOS_ELEMENT


# --- dual-label positions (Phase 2 lookups) -----------------------------------
# A position inside two professor intervals (PXR 261: H3 and H4; CAR 344: H11' and
# H12) is a member of BOTH elements. Statistics and aggregation use the membership
# tuple (one entry per element). The joined string is for display and one-row
# lookups only; it is never an element of its own.
DISPLAY_JOIN = "|"


def position_memberships(rows) -> dict[tuple, tuple[str, ...]]:
    """(uniprot_id, position) -> professor elements containing it, in element order.

    `rows` yields (uniprot_id, position, structural_segment); blank segments are
    ignored, so inter-helix and unmapped positions are simply absent."""
    out: dict[tuple, set] = {}
    for uid, pos, seg in rows:
        if seg:
            out.setdefault((uid, int(pos)), set()).add(seg)
    return {k: tuple(sorted(v, key=element_rank)) for k, v in out.items()}


def display_label(labels) -> str:
    """Joined label for display / one-row lookup only (e.g. "H3|H4")."""
    return DISPLAY_JOIN.join(labels)
