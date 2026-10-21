"""
Apo/holo analysis scope (approved decision S1, Phase 2 of the v2 rebuild).

A receptor x PDB entry whose every coordinate candidate for that receptor has
lbd_mapping_category == TRACE_LBD is a short receptor-derived fragment (peptide)
bound to another protein, not an LBD structure. Such entries stay in the master
and provenance tables, but they cannot count toward APO or HOLO eligibility or
enter any apo/holo test. The rule is symmetric: APO and HOLO alike.

Entries with at least one non-TRACE_LBD candidate are kept unchanged (their
ligand state is whatever the foundation assigned, e.g. AMBIGUOUS when entity
states disagree).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXCLUDED_STATE = "EXCLUDED_TRACE_LBD"
TRACE_CATEGORY = "TRACE_LBD"


def trace_lbd_entries(project_root: Path = PROJECT_ROOT) -> set[tuple[str, str]]:
    """(uniprot_id, pdb_id) pairs whose every coordinate candidate is TRACE_LBD."""
    s = pd.read_csv(project_root / "data" / "processed" / "structures.csv", dtype={"pdb_id": str},
                    keep_default_na=False, low_memory=False)
    coord = s[s["coordinate_data_status"] == "COORDINATE_DATA_AVAILABLE"]
    only_trace = coord.groupby(["uniprot_id", "pdb_id"])["lbd_mapping_category"].agg(
        lambda c: set(c) == {TRACE_CATEGORY})
    return set(only_trace[only_trace].index)


def apply(df: pd.DataFrame, project_root: Path = PROJECT_ROOT, state_col: str = "ligand_state") -> pd.DataFrame:
    """Copy of df with ligand_state set to EXCLUDED_TRACE_LBD for out-of-scope entries,
    so every APO/HOLO filter downstream drops them on both sides."""
    excluded = trace_lbd_entries(project_root)
    out = df.copy()
    mask = pd.Series([(u, p) in excluded for u, p in zip(out["uniprot_id"], out["pdb_id"])], index=out.index)
    out.loc[mask, state_col] = EXCLUDED_STATE
    return out


def drop(df: pd.DataFrame, project_root: Path = PROJECT_ROOT) -> pd.DataFrame:
    """df without the out-of-scope (all-TRACE_LBD) receptor x PDB entries. Used only by the
    labelled TRACE_LBD functional-site sensitivity run; the primary site analysis keeps them."""
    excluded = trace_lbd_entries(project_root)
    keep = [(u, p) not in excluded for u, p in zip(df["uniprot_id"], df["pdb_id"])]
    return df[keep].copy()
