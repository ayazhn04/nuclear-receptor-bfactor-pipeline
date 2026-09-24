"""Low-level mmCIF parsing helpers built on Bio.PDB.MMCIF2Dict, giving
direct access to raw `_atom_site` parallel columns (label_seq_id,
label_asym_id, label_atom_id, occupancy, B_iso_or_equiv, label_alt_id,
group_PDB, auth_*) — the exact fields Stage 3B's per-residue observation
model requires. Schema fields confirmed by live inspection of a real
downloaded file (see reports/tables/stage3b_coordinate_completeness_report.md
"Methods").
"""

from __future__ import annotations

import gzip
import os
import tempfile
from pathlib import Path
from typing import Any

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def parse_mmcif_gz(path: Path) -> dict[str, Any]:
    """Decompress and parse one .cif.gz file into an MMCIF2Dict. MMCIF2Dict
    requires a real file path, so the decompressed text is written to a
    temporary file first (removed immediately after parsing)."""
    with gzip.open(path, "rt", encoding="utf-8") as f:
        text = f.read()
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as tmp:
            tmp.write(text)
            tmp_path = tmp.name
        return MMCIF2Dict(tmp_path)
    finally:
        if tmp_path is not None:
            os.unlink(tmp_path)


def as_list(value) -> list:
    """MMCIF2Dict returns a plain string (not a list) for single-row loops.
    Normalize to a list either way."""
    if isinstance(value, list):
        return value
    return [value]


def to_float_or_none(value: str) -> float | None:
    if value in (None, "?", ".", ""):
        return None
    try:
        f = float(value)
    except ValueError:
        return None
    if f != f or f in (float("inf"), float("-inf")):  # NaN / inf guard
        return None
    return f


def to_int_or_none(value: str) -> int | None:
    if value in (None, "?", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


class AtomSiteIndex:
    """Pre-indexed view of `_atom_site` for O(1) lookups by
    (label_asym_id, model_num, label_seq_id) -> list of CA atom rows.

    Only rows with label_atom_id == 'CA' and a valid integer label_seq_id
    are indexed (all other atoms are irrelevant to Stage 3B's per-residue
    Cα observation model).
    """

    def __init__(self, mmcif_dict: dict[str, Any]):
        if "_atom_site.id" not in mmcif_dict:
            self.available = False
            self.index: dict[tuple[str, str], dict[int, list[dict]]] = {}
            self.asym_to_entity: dict[str, str] = {}
            self.model_nums: set[str] = set()
            return

        self.available = True
        label_asym = as_list(mmcif_dict["_atom_site.label_asym_id"])
        label_atom = as_list(mmcif_dict["_atom_site.label_atom_id"])
        label_seq = as_list(mmcif_dict["_atom_site.label_seq_id"])
        label_alt = as_list(mmcif_dict["_atom_site.label_alt_id"])
        label_entity = as_list(mmcif_dict["_atom_site.label_entity_id"])
        label_comp = as_list(mmcif_dict["_atom_site.label_comp_id"])
        auth_comp = as_list(mmcif_dict["_atom_site.auth_comp_id"])
        auth_seq = as_list(mmcif_dict["_atom_site.auth_seq_id"])
        auth_asym = as_list(mmcif_dict["_atom_site.auth_asym_id"])
        ins_code = as_list(mmcif_dict.get("_atom_site.pdbx_PDB_ins_code", []))
        group_pdb = as_list(mmcif_dict["_atom_site.group_PDB"])
        occupancy = as_list(mmcif_dict["_atom_site.occupancy"])
        b_iso = as_list(mmcif_dict["_atom_site.B_iso_or_equiv"])
        model_num = as_list(mmcif_dict["_atom_site.pdbx_PDB_model_num"])
        atom_id = as_list(mmcif_dict["_atom_site.id"])

        n = len(atom_id)
        has_ins_code = len(ins_code) == n

        self.index = {}
        self.asym_to_entity = {}
        self.model_nums = set()

        for i in range(n):
            if label_atom[i] != "CA":
                continue
            seq_id = to_int_or_none(label_seq[i])
            if seq_id is None:
                continue
            asym = label_asym[i]
            model = model_num[i]
            self.model_nums.add(model)
            self.asym_to_entity[asym] = label_entity[i]

            key = (asym, model)
            self.index.setdefault(key, {}).setdefault(seq_id, []).append(
                {
                    "atom_id": atom_id[i],
                    "label_alt_id": label_alt[i],
                    "label_comp_id": label_comp[i],
                    "auth_comp_id": auth_comp[i],
                    "auth_seq_id": auth_seq[i],
                    "auth_asym_id": auth_asym[i],
                    "pdbx_PDB_ins_code": ins_code[i] if has_ins_code else "",
                    "group_PDB": group_pdb[i],
                    "occupancy": to_float_or_none(occupancy[i]),
                    "B_iso_or_equiv": to_float_or_none(b_iso[i]),
                }
            )

    def get_residue_ca_rows(self, label_asym_id: str, model_num: str, label_seq_id: int) -> list[dict]:
        return self.index.get((label_asym_id, model_num), {}).get(label_seq_id, [])


def get_unobserved_residues(mmcif_dict: dict[str, Any]) -> list[dict]:
    """Returns [{label_asym_id, label_seq_id, model_num, occupancy_flag}]
    from `_pdbx_unobs_or_zero_occ_residues`, or [] if the category is
    absent from this entry."""
    if "_pdbx_unobs_or_zero_occ_residues.label_asym_id" not in mmcif_dict:
        return []
    asym = as_list(mmcif_dict["_pdbx_unobs_or_zero_occ_residues.label_asym_id"])
    seq = as_list(mmcif_dict["_pdbx_unobs_or_zero_occ_residues.label_seq_id"])
    model = as_list(mmcif_dict.get("_pdbx_unobs_or_zero_occ_residues.PDB_model_num", ["1"] * len(asym)))
    occ_flag = as_list(mmcif_dict.get("_pdbx_unobs_or_zero_occ_residues.occupancy_flag", [""] * len(asym)))
    rows = []
    for a, s, m, o in zip(asym, seq, model, occ_flag):
        seq_id = to_int_or_none(s)
        if seq_id is None:
            continue
        rows.append({"label_asym_id": a, "label_seq_id": seq_id, "model_num": m, "occupancy_flag": o})
    return rows


def get_entity_types(mmcif_dict: dict[str, Any]) -> dict[str, str]:
    """{entity_id: type} from `_entity`, e.g. {'1': 'polymer', '2':
    'non-polymer', '3': 'water'}."""
    if "_entity.id" not in mmcif_dict:
        return {}
    ids = as_list(mmcif_dict["_entity.id"])
    types = as_list(mmcif_dict["_entity.type"])
    return dict(zip(ids, types))


def get_entity_descriptions(mmcif_dict: dict[str, Any]) -> dict[str, str]:
    if "_entity.id" not in mmcif_dict:
        return {}
    ids = as_list(mmcif_dict["_entity.id"])
    descs = as_list(mmcif_dict.get("_entity.pdbx_description", [""] * len(ids)))
    return dict(zip(ids, descs))


def get_struct_site_records(mmcif_dict: dict[str, Any]) -> list[dict]:
    """Deposited (often software-detected, occasionally author-curated)
    binding-site records from `_struct_site`, keyed by the auth-numbering
    identity of the site's subject component — the exact fields
    `_struct_site` uses, not `_struct_site_gen` (which lists the
    surrounding residues, not the ligand itself)."""
    if "_struct_site.id" not in mmcif_dict:
        return []
    site_id = as_list(mmcif_dict["_struct_site.id"])
    n = len(site_id)
    auth_comp = as_list(mmcif_dict.get("_struct_site.pdbx_auth_comp_id", [""] * n))
    auth_asym = as_list(mmcif_dict.get("_struct_site.pdbx_auth_asym_id", [""] * n))
    auth_seq = as_list(mmcif_dict.get("_struct_site.pdbx_auth_seq_id", [""] * n))
    evidence = as_list(mmcif_dict.get("_struct_site.pdbx_evidence_code", [""] * n))
    details = as_list(mmcif_dict.get("_struct_site.details", [""] * n))
    return [
        {
            "site_id": site_id[i], "auth_comp_id": auth_comp[i], "auth_asym_id": auth_asym[i],
            "auth_seq_id": auth_seq[i], "evidence_code": evidence[i], "details": details[i],
        }
        for i in range(n)
    ]


def get_covalent_linkages(mmcif_dict: dict[str, Any]) -> list[dict]:
    """Covalent connections from `_struct_conn` where conn_type_id contains
    'covale' (covalent linkage between a nonpolymer component and the
    polymer, or between two nonpolymer components) — standard mmCIF
    dictionary field names (ptnr1/ptnr2 label_asym_id/label_comp_id/
    label_seq_id)."""
    if "_struct_conn.id" not in mmcif_dict:
        return []
    conn_type = as_list(mmcif_dict.get("_struct_conn.conn_type_id", []))
    n = len(conn_type)
    if n == 0:
        return []
    p1_asym = as_list(mmcif_dict.get("_struct_conn.ptnr1_label_asym_id", [""] * n))
    p1_comp = as_list(mmcif_dict.get("_struct_conn.ptnr1_label_comp_id", [""] * n))
    p1_seq = as_list(mmcif_dict.get("_struct_conn.ptnr1_label_seq_id", [""] * n))
    p2_asym = as_list(mmcif_dict.get("_struct_conn.ptnr2_label_asym_id", [""] * n))
    p2_comp = as_list(mmcif_dict.get("_struct_conn.ptnr2_label_comp_id", [""] * n))
    p2_seq = as_list(mmcif_dict.get("_struct_conn.ptnr2_label_seq_id", [""] * n))
    rows = []
    for i in range(n):
        if "coval" not in conn_type[i].lower():
            continue
        rows.append(
            {
                "conn_type": conn_type[i],
                "ptnr1_label_asym_id": p1_asym[i], "ptnr1_label_comp_id": p1_comp[i], "ptnr1_label_seq_id": p1_seq[i],
                "ptnr2_label_asym_id": p2_asym[i], "ptnr2_label_comp_id": p2_comp[i], "ptnr2_label_seq_id": p2_seq[i],
            }
        )
    return rows


def get_heavy_atoms(mmcif_dict: dict[str, Any]) -> list[dict]:
    """Every non-hydrogen `_atom_site` row (all entities, all models), with
    Cartesian coordinates. Used for nonpolymer-inventory enumeration and
    LBD-contact geometry — a separate, fuller pass than AtomSiteIndex
    (which only indexes Cα rows for the Stage 3B coverage model)."""
    if "_atom_site.id" not in mmcif_dict:
        return []
    label_asym = as_list(mmcif_dict["_atom_site.label_asym_id"])
    n = len(label_asym)
    label_atom = as_list(mmcif_dict["_atom_site.label_atom_id"])
    label_seq = as_list(mmcif_dict["_atom_site.label_seq_id"])
    label_alt = as_list(mmcif_dict["_atom_site.label_alt_id"])
    label_entity = as_list(mmcif_dict["_atom_site.label_entity_id"])
    label_comp = as_list(mmcif_dict["_atom_site.label_comp_id"])
    auth_seq = as_list(mmcif_dict["_atom_site.auth_seq_id"])
    auth_asym = as_list(mmcif_dict["_atom_site.auth_asym_id"])
    ins_code = as_list(mmcif_dict.get("_atom_site.pdbx_PDB_ins_code", []))
    has_ins_code = len(ins_code) == n
    group_pdb = as_list(mmcif_dict["_atom_site.group_PDB"])
    occupancy = as_list(mmcif_dict["_atom_site.occupancy"])
    model_num = as_list(mmcif_dict["_atom_site.pdbx_PDB_model_num"])
    type_symbol = as_list(mmcif_dict["_atom_site.type_symbol"])
    x = as_list(mmcif_dict["_atom_site.Cartn_x"])
    y = as_list(mmcif_dict["_atom_site.Cartn_y"])
    z = as_list(mmcif_dict["_atom_site.Cartn_z"])

    rows = []
    for i in range(n):
        element = type_symbol[i].strip().upper()
        if element in ("H", "D"):
            continue
        rows.append(
            {
                "label_asym_id": label_asym[i], "label_atom_id": label_atom[i],
                "label_seq_id": to_int_or_none(label_seq[i]), "label_alt_id": label_alt[i],
                "label_entity_id": label_entity[i], "label_comp_id": label_comp[i],
                "auth_seq_id": auth_seq[i], "auth_asym_id": auth_asym[i],
                "pdbx_PDB_ins_code": ins_code[i] if has_ins_code else "",
                "group_PDB": group_pdb[i], "occupancy": to_float_or_none(occupancy[i]),
                "model_num": model_num[i], "element": element,
                "x": to_float_or_none(x[i]), "y": to_float_or_none(y[i]), "z": to_float_or_none(z[i]),
            }
        )
    return rows


def get_entry_metadata(mmcif_dict: dict[str, Any]) -> dict[str, Any]:
    entry_id = as_list(mmcif_dict.get("_entry.id", [""]))[0]
    deposition_date = as_list(mmcif_dict.get("_pdbx_database_status.recvd_initial_deposition_date", [""]))[0]
    revision_dates = as_list(mmcif_dict.get("_pdbx_audit_revision_history.revision_date", []))
    latest_revision = max(revision_dates) if revision_dates else ""
    anisotrop_present = "_atom_site_anisotrop.id" in mmcif_dict
    anisotrop_count = len(as_list(mmcif_dict["_atom_site_anisotrop.id"])) if anisotrop_present else 0
    anisotrop_atom_ids = set(as_list(mmcif_dict["_atom_site_anisotrop.id"])) if anisotrop_present else set()
    return {
        "entry_id": entry_id,
        "deposition_date": deposition_date,
        "revision_dates": revision_dates,
        "latest_revision_date": latest_revision,
        "anisotrop_category_present": anisotrop_present,
        "anisotrop_atom_count": anisotrop_count,
        "anisotrop_atom_ids": anisotrop_atom_ids,
    }
