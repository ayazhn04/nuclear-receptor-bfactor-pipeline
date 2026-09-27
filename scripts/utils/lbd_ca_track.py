"""Extract an ordered (Ca coordinate, canonical_uniprot_position) track for
one representative coordinate instance's standardized LBD, for use as CE
structural-alignment guide atoms. Positive-occupancy, construct-mapped
residues only (the same "observed" definition used throughout this
project since Stage 3B.1)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.utils.mmcif_parse import get_heavy_atoms, parse_mmcif_gz

RAW_MMCIF_DIR = None  # set by caller via set_mmcif_dir


def build_ca_track(
    pdb_id: str, label_asym_id: str, residues: pd.DataFrame, mmcif_dir: Path,
) -> tuple[list[list[float]], list[int]]:
    """residues: rows of the normalized/observation table for exactly this
    (pdb_id, label_asym_id) instance, already filtered to
    is_construct_mapped & ca_positive_occupancy, with entity_label_seq_id and
    canonical_uniprot_position columns, sorted by entity_label_seq_id."""
    d = parse_mmcif_gz(mmcif_dir / f"{pdb_id}.cif.gz")
    atoms = get_heavy_atoms(d)
    by_seq: dict[int, dict] = {}
    for a in atoms:
        if a["label_asym_id"] != label_asym_id or a["label_atom_id"] != "CA":
            continue
        if not ((a["occupancy"] or 0) > 0):
            continue
        sid = a["label_seq_id"]
        if sid is None:
            continue
        # first positive-occupancy CA record wins ties (mirrors the
        # project's altloc-agnostic "observed" definition -- any positive-
        # occupancy record satisfies observation).
        by_seq.setdefault(sid, a)

    coords, positions = [], []
    for _, r in residues.sort_values("entity_label_seq_id").iterrows():
        sid = int(r["entity_label_seq_id"])
        a = by_seq.get(sid)
        if a is not None:
            coords.append([a["x"], a["y"], a["z"]])
            positions.append(int(r["canonical_uniprot_position"]))
    return coords, positions
