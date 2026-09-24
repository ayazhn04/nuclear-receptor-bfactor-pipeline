"""One-off recovery helper: re-extract the deposited unobserved-residue
ground truth (`_pdbx_unobs_or_zero_occ_residues`) from the already-cached
mmCIF files (no download, no network) for every selected PDB entry.

This exists because an earlier Stage 3B.1 script accidentally overwrote
data/manifests/mmcif_unobserved_residue_crosscheck.csv (the Stage 3B
ground-truth table) with an empty, buggy recomputation before the bug was
caught. Rather than trust a partially-overwritten file, this rebuilds the
ground truth directly and deterministically from the local .cif.gz cache.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.mmcif_parse import get_unobserved_residues, parse_mmcif_gz  # noqa: E402

RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"


def rebuild(pdb_ids: list[str]) -> pd.DataFrame:
    rows = []
    for pdb_id in pdb_ids:
        path = RAW_MMCIF_DIR / f"{pdb_id.upper()}.cif.gz"
        if not path.exists():
            continue
        d = parse_mmcif_gz(path)
        for r in get_unobserved_residues(d):
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "label_asym_id": r["label_asym_id"],
                    "model_num": r["model_num"],
                    "label_seq_id": r["label_seq_id"],
                }
            )
    return pd.DataFrame(rows)
