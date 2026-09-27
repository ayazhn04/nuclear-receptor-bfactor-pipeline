"""
Real 3D protein structural alignment: Combinatorial Extension (CE;
Shindyalov & Bourne 1998, "Protein structure alignment by incremental
combinatorial extension of the optimal path", Protein Engineering 11(9):
739-747). No new software was installed for this -- CE is already shipped
inside Biopython (`Bio.PDB.ccealign`, used by `Bio.PDB.CEAligner`, originally
adapted from open-source PyMOL). This module calls the same underlying
`run_cealign` C routine directly on raw Ca coordinate lists (bypassing
Biopython's Structure-object wrapper, which we don't need since our
coordinates already come from this project's own mmCIF/label_seq_id
pipeline), and additionally computes the Kabsch-superposition RMSD and a
sequence-order-respecting residue-residue correspondence -- the explicit
output this project's structural-segment mapping requires, not just a
global score.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from Bio.PDB.ccealign import run_cealign

CE_WINDOW_SIZE = 8
CE_MAX_GAP = 30


@dataclass
class StructuralAlignmentResult:
    n_aligned: int
    z_score: float
    rmsd: float
    # index pairs into the ORIGINAL coordinate/position arrays passed in
    idx_a: list[int]
    idx_b: list[int]


def _kabsch_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    a = coords_a - coords_a.mean(axis=0)
    b = coords_b - coords_b.mean(axis=0)
    h = a.T @ b
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    corr = np.diag([1, 1, d])
    r = vt.T @ corr @ u.T
    a_rot = (r @ a.T).T
    return float(np.sqrt(np.mean(np.sum((a_rot - b) ** 2, axis=1))))


def align_structures(
    coords_a: list[list[float]], coords_b: list[list[float]],
    window_size: int = CE_WINDOW_SIZE, max_gap: int = CE_MAX_GAP,
) -> StructuralAlignmentResult | None:
    """coords_a/coords_b: ordered lists of [x, y, z] guide (Ca) coordinates
    (structure A = reference/template, B = mobile/target), already restricted
    to the region of interest (here: the standardized LBD's observed,
    positive-occupancy Ca atoms) and sorted by residue order."""
    if len(coords_a) < window_size * 2 or len(coords_b) < window_size * 2:
        return None

    alignments = run_cealign(coords_a, coords_b, window_size, max_gap)
    if not alignments:
        return None
    longest = [a for a in alignments if a.length == alignments[0].length]

    best = None
    best_rmsd = float("inf")
    for aln in longest:
        idx_a, idx_b = aln.path
        ca = np.array([coords_a[i] for i in idx_a])
        cb = np.array([coords_b[i] for i in idx_b])
        rmsd = _kabsch_rmsd(ca, cb)
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best = aln

    idx_a, idx_b = best.path
    # CE paths are already monotonic in both index arrays by construction
    # (it aligns ordered fragments); assert this explicitly since downstream
    # label transfer relies on residue-order continuity.
    assert all(x < y for x, y in zip(idx_a, idx_a[1:])), "CE path not monotonic in A"
    assert all(x < y for x, y in zip(idx_b, idx_b[1:])), "CE path not monotonic in B"

    return StructuralAlignmentResult(
        n_aligned=len(idx_a), z_score=float(best.z_score), rmsd=best_rmsd,
        idx_a=list(idx_a), idx_b=list(idx_b),
    )
