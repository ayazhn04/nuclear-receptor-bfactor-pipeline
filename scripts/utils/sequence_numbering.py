"""
Sequence-evidence residue numbering: map a chain's own polymer sequence onto
a UniProt sequence by explicit alignment, deletion/insertion aware, without
ever forcing a mapping (used by scripts/96_build_residue_numbering_crosswalk.py).

Evidence hierarchy (strongest first):
  1. the chain's polymer sequence from mmCIF _pdbx_poly_seq_scheme
  2. the residues actually observed in that chain (cross-check)
  3. explicit alignment of (1) against the target UniProt sequence
  4. author numbering -- supporting metadata only
  5. the previous RCSB-service mapping -- comparison only
Depositor-declared deletions/mutations (_struct_ref_seq, _struct_ref_seq_dif)
are recorded as corroboration, never as the mapping itself.

A residue gets a position only if (a) every co-optimal alignment places it
at the same UniProt position and (b) the alignment around it is locally
well-supported. Otherwise it is withheld. B-factors are never used.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio import Align
from Bio.Align import substitution_matrices

from scripts.utils.mmcif_parse import as_list

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
# Common modified residues -> parent (used only when the file does not name the parent itself).
MODIFIED_PARENT = {"MSE": "MET", "CSO": "CYS", "CME": "CYS", "CSD": "CYS", "OCS": "CYS", "CAS": "CYS", "CSS": "CYS",
                   "SEP": "SER", "TPO": "THR", "PTR": "TYR", "KCX": "LYS", "MLY": "LYS", "HYP": "PRO", "ALY": "LYS"}

GAP_OPEN, GAP_EXTEND = -10.0, -0.5
WINDOW = 10                 # one-sided window of aligned pairs used for local support
MIN_WINDOW_IDENTITY = 0.6   # substitutions: support needed on BOTH sides
MIN_EDGE_IDENTITY = 0.8     # identical residues: strong support on one side, or MIN_WINDOW_IDENTITY on both
MIN_ALIGNED_IDENTICAL = 20  # whole-chain evidence floor
MIN_ALIGNMENT_IDENTITY = 0.5
MIN_PEPTIDE_IDENTICAL = 8   # short peptide chains: exact, unique match of >= 8 residues
MAX_OPTIMAL_ALIGNMENTS = 200

_ALIGNER = None


def aligner() -> Align.PairwiseAligner:
    global _ALIGNER
    if _ALIGNER is None:
        _ALIGNER = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                         open_gap_score=GAP_OPEN, extend_gap_score=GAP_EXTEND)
    return _ALIGNER


def one_letter(mon_id: str, parents: dict[str, str] | None = None) -> str:
    mon = mon_id.upper()
    if mon in THREE_TO_ONE:
        return THREE_TO_ONE[mon]
    parent = (parents or {}).get(mon) or MODIFIED_PARENT.get(mon)
    return THREE_TO_ONE.get(parent, "X") if parent else "X"


# ---------------------------------------------------------------------------
# mmCIF evidence
# ---------------------------------------------------------------------------


def modified_residue_parents(d: dict) -> dict[str, str]:
    if "_pdbx_struct_mod_residue.label_comp_id" not in d:
        return {}
    return dict(zip(as_list(d["_pdbx_struct_mod_residue.label_comp_id"]),
                    as_list(d["_pdbx_struct_mod_residue.parent_comp_id"])))


def chain_scheme(d: dict, label_asym: str) -> list[dict]:
    """Polymer residues of one chain in label_seq_id order (first residue at a
    microheterogeneous position; heterogeneity flagged)."""
    asym = as_list(d["_pdbx_poly_seq_scheme.asym_id"])
    cols = {c: as_list(d[f"_pdbx_poly_seq_scheme.{c}"]) for c in
            ["seq_id", "mon_id", "auth_seq_num", "pdb_mon_id", "pdb_ins_code"]}
    out: dict[int, dict] = {}
    for i, a in enumerate(asym):
        if a != label_asym:
            continue
        sid = int(cols["seq_id"][i])
        if sid in out:
            out[sid]["heterogeneous"] = True
            continue
        out[sid] = {"label_seq_id": sid, "mon_id": cols["mon_id"][i], "auth_seq_num": cols["auth_seq_num"][i],
                    "pdb_ins_code": cols["pdb_ins_code"][i], "modelled": cols["pdb_mon_id"][i] not in ("?", "."),
                    "heterogeneous": False}
    return [out[k] for k in sorted(out)]


def seq_dif_for_strand(d: dict, auth_strand: str, accession: str) -> tuple[dict[int, str], set[int]]:
    """Depositor-declared differences for one strand: {label_seq_id: details}
    for residues present in the chain, and the set of reference (UniProt)
    positions declared deleted."""
    if "_struct_ref_seq_dif.details" not in d:
        return {}, set()
    n = len(as_list(d["_struct_ref_seq_dif.details"]))
    get = lambda c: as_list(d.get(f"_struct_ref_seq_dif.{c}", ["?"] * n))  # noqa: E731
    strand, seq_num, db_num, details, acc = (get("pdbx_pdb_strand_id"), get("seq_num"), get("pdbx_seq_db_seq_num"),
                                             get("details"), get("pdbx_seq_db_accession_code"))
    by_label, deleted = {}, set()
    for s, sn, dn, det, a in zip(strand, seq_num, db_num, details, acc):
        if s != auth_strand:
            continue
        det = (det or "").strip().lower()
        if det == "deletion":
            if a in (accession, "?", "") and str(dn).lstrip("-").isdigit():
                deleted.add(int(dn))
        elif str(sn).isdigit():
            by_label[int(sn)] = det
    return by_label, deleted


def depositor_segments(d: dict, auth_strand: str, accession: str) -> list[tuple[int, int, int, int]]:
    """(seq_align_beg, seq_align_end, db_align_beg, db_align_end) for this strand/accession."""
    if "_struct_ref_seq.pdbx_strand_id" not in d:
        return []
    rows = zip(*(as_list(d[f"_struct_ref_seq.{c}"]) for c in
                 ["pdbx_strand_id", "seq_align_beg", "seq_align_end", "db_align_beg", "db_align_end",
                  "pdbx_db_accession"]))
    out = []
    for s, sb, se, db, de, a in rows:
        if s == auth_strand and a.split("-")[0] == accession.split("-")[0] and all(
                str(x).lstrip("-").isdigit() for x in (sb, se, db, de)):
            out.append((int(sb), int(se), int(db), int(de)))
    return sorted(out)


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------


@dataclass
class AlignmentResult:
    mapping: dict[int, int]                  # query index -> target index (0-based), consensus only
    withheld_alternative: set[int]           # query indices placed differently by co-optimal alignments
    n_optimal: int
    blocks: list[tuple[tuple[int, int], tuple[int, int]]]  # gapless blocks of the first optimal alignment
    score: float
    window_left: dict[int, float] = field(default_factory=dict)
    window_right: dict[int, float] = field(default_factory=dict)


def _pairs(aln) -> list[tuple[int, int]]:
    out = []
    for (qs, qe), (ts, te) in zip(*aln.aligned):
        out.extend(zip(range(qs, qe), range(ts, te)))
    return out


def align_unique(query: str, target: str) -> AlignmentResult:
    alns = aligner().align(query, target)
    n_opt = len(alns)
    first = alns[0]
    first_pairs = dict(_pairs(first))
    withheld: set[int] = set()
    if n_opt > 1:
        for k, aln in enumerate(alns):
            if k >= MAX_OPTIMAL_ALIGNMENTS:
                break
            other = dict(_pairs(aln))
            for q in set(first_pairs) | set(other):
                if first_pairs.get(q) != other.get(q):
                    withheld.add(q)
        if n_opt > MAX_OPTIMAL_ALIGNMENTS:  # cannot enumerate all: be conservative
            withheld = set(first_pairs)
    mapping = {q: t for q, t in first_pairs.items() if q not in withheld}
    blocks = [((qs, qe), (ts, te)) for (qs, qe), (ts, te) in zip(*first.aligned)]
    res = AlignmentResult(mapping, withheld, n_opt, blocks, float(first.score))
    ordered = sorted(first_pairs.items())
    ident = [query[q] == target[t] for q, t in ordered]
    for i, (q, _) in enumerate(ordered):
        left = ident[max(0, i - WINDOW + 1): i + 1]
        right = ident[i: i + WINDOW]
        res.window_left[q] = sum(left) / len(left)
        res.window_right[q] = sum(right) / len(right)
    return res


def gaps_in_alignment(blocks) -> list[dict]:
    """Deletions (target positions skipped) and insertions (query residues
    skipped) between consecutive gapless blocks, 0-based."""
    out = []
    for ((q1s, q1e), (t1s, t1e)), ((q2s, q2e), (t2s, t2e)) in zip(blocks, blocks[1:]):
        if t2s > t1e:
            out.append({"type": "DELETION_IN_CHAIN", "after_query": q1e - 1, "target_from": t1e, "target_to": t2s - 1})
        if q2s > q1e:
            out.append({"type": "INSERTION_IN_CHAIN", "after_query": q1e - 1, "query_from": q1e, "query_to": q2s - 1})
    return out


@dataclass
class ResidueCall:
    target_index: int | None   # 0-based UniProt index, None = withheld/unaligned
    reason: str                # MAPPED / WITHHELD_* / NOT_ALIGNED
    window_left: float | None
    window_right: float | None


NON_REFERENCE_SEQ_DIF = ("tag", "linker", "cloning", "insertion", "initiating", "leader", "expression")


def is_non_reference_residue(seq_dif_details: str) -> bool:
    """Depositor declares the residue is not part of the reference sequence
    (expression tag, linker, cloning artifact, ...): never map it."""
    d = (seq_dif_details or "").lower()
    return any(k in d for k in NON_REFERENCE_SEQ_DIF)


def isoform_missing_match(lo: int, hi: int, missing_segments: list[tuple[int, int]]) -> bool:
    """A chain deletion is an isoform difference only if it has the same length
    as a UniProt isoform 'Missing' segment and starts within +-1 of it (the
    +-1 covers co-optimal placement of a repeated boundary residue)."""
    return any(hi - lo == b - a and abs(lo - a) <= 1 for a, b in missing_segments)


def call_residues(query: str, target: str, res: AlignmentResult,
                  non_reference: set[int] | None = None) -> tuple[list[ResidueCall], dict]:
    """Per-query-residue decision plus whole-chain evidence summary.
    `non_reference`: query indices the depositor declares as tag/linker/etc."""
    non_reference = non_reference or set()
    n_ident = sum(1 for q, t in res.mapping.items() if query[q] == target[t])
    aligned = len(res.mapping) + len(res.withheld_alternative)
    identity = n_ident / len(res.mapping) if res.mapping else 0.0
    full_chain = n_ident >= MIN_ALIGNED_IDENTICAL and identity >= MIN_ALIGNMENT_IDENTITY
    # A short peptide (e.g. a receptor-derived peptide bound to another receptor)
    # cannot reach MIN_ALIGNED_IDENTICAL, but an exact match with a single
    # optimal placement is unambiguous sequence evidence.
    short_peptide = (not full_chain and n_ident >= MIN_PEPTIDE_IDENTICAL and n_ident == len(res.mapping)
                     and res.n_optimal == 1 and not res.withheld_alternative)
    sufficient = full_chain or short_peptide
    calls = []
    for q in range(len(query)):
        wl, wr = res.window_left.get(q), res.window_right.get(q)
        if q in non_reference:
            calls.append(ResidueCall(None, "WITHHELD_DEPOSITOR_NON_REFERENCE_RESIDUE", wl, wr))
        elif q in res.withheld_alternative:
            calls.append(ResidueCall(None, "WITHHELD_ALTERNATIVE_OPTIMAL_ALIGNMENT", wl, wr))
        elif q not in res.mapping:
            calls.append(ResidueCall(None, "NOT_ALIGNED", wl, wr))
        elif not sufficient:
            calls.append(ResidueCall(None, "WITHHELD_INSUFFICIENT_CHAIN_EVIDENCE", wl, wr))
        elif query[q] == target[res.mapping[q]] and max(wl, wr) < MIN_EDGE_IDENTITY and min(wl, wr) < MIN_WINDOW_IDENTITY:
            calls.append(ResidueCall(None, "WITHHELD_LOW_LOCAL_IDENTITY", wl, wr))
        elif query[q] != target[res.mapping[q]] and min(wl, wr) < MIN_WINDOW_IDENTITY:
            # a non-identical pair is only a point substitution when the alignment
            # is supported on both sides; at an alignment edge it may be divergence
            calls.append(ResidueCall(None, "WITHHELD_UNSUPPORTED_SUBSTITUTION", wl, wr))
        else:
            calls.append(ResidueCall(res.mapping[q], "MAPPED", wl, wr))
    summary = {"n_aligned_pairs": aligned, "n_identical": n_ident, "alignment_identity": round(identity, 4),
               "n_optimal_alignments": res.n_optimal, "sufficient_evidence": sufficient,
               "evidence_basis": ("FULL_CHAIN_ALIGNMENT" if full_chain else
                                  "SHORT_PEPTIDE_EXACT_UNIQUE" if short_peptide else "INSUFFICIENT"),
               "alignment_score": res.score}
    return calls, summary
