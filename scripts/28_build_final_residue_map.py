"""
Final Data Release, Steps 11-12: canonical residue handoff.

Builds one row per (selected primary representative structure x
standardized-LBD canonical residue), covering EVERY standardized LBD
position (including ones with no observed coordinate, so downstream users
never have to guess whether a position was silently dropped), plus a
deterministic single recommended raw Cα record per residue (occupancy and
altloc-identity only — B-factor magnitude never participates in the
choice; see scripts/utils/altloc_selection.py).

Reads:
    data/manifests/final_structure_selection_audit.csv
    data/interim/qc/stage3b_lbd_ca_observations.parquet
    data/raw/api/uniprot/{UNIPROT_ID}.json   (already-cached canonical sequence)

Writes:
    data/processed/final_lbd_residue_map.parquet
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.altloc_selection import select_representative_ca_record  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
OBSERVATIONS_PARQUET = PROJECT_ROOT / "data" / "interim" / "qc" / "stage3b_lbd_ca_observations.parquet"
UNIPROT_RAW_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"

SELECTION_AUDIT_CSV = MANIFESTS_DIR / "final_structure_selection_audit.csv"
RESIDUE_MAP_PARQUET = PROCESSED_DIR / "final_lbd_residue_map.parquet"

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL", "U": "SEC", "O": "PYL",
}


def load_canonical_sequences(uniprot_ids: list[str]) -> dict[str, str]:
    sequences = {}
    for uid in uniprot_ids:
        path = UNIPROT_RAW_DIR / f"{uid}.json"
        if not path.exists():
            raise RuntimeError(f"STOPPING: no cached UniProt entry for {uid} at {path}")
        raw = json.loads(path.read_text())
        seq = raw.get("sequence", {}).get("value")
        if not seq:
            raise RuntimeError(f"STOPPING: cached UniProt entry for {uid} has no sequence.value")
        sequences[uid] = seq
    return sequences


def canonical_residue_name(sequence: str, canonical_position: int) -> str:
    if canonical_position < 1 or canonical_position > len(sequence):
        return ""
    one_letter = sequence[canonical_position - 1]
    return ONE_TO_THREE.get(one_letter, one_letter)


def main() -> None:
    audit = pd.read_csv(SELECTION_AUDIT_CSV, dtype={"pdb_id": str})
    primary = audit[audit["in_primary_set"]].copy()
    print(f"Primary representative structures: {len(primary)}")

    obs = pd.read_parquet(OBSERVATIONS_PARQUET)
    rep_instance_ids = set(primary["representative_instance_id"])
    obs_primary = obs[obs["instance_id"].isin(rep_instance_ids)].copy()
    print(f"Observation rows for primary representative instances: {len(obs_primary)} "
          f"(expect one row per LBD position per primary structure)")

    canonical_sequences = load_canonical_sequences(sorted(primary["uniprot_id"].unique()))

    rows = []
    for _, r in obs_primary.iterrows():
        seq = canonical_sequences[r["uniprot_id"]]
        occupancies = json.loads(r["occupancies_json"])
        b_isos = json.loads(r["b_iso_values_json"])
        altlocs = json.loads(r["altloc_ids_json"])
        selected = select_representative_ca_record(occupancies, b_isos, altlocs)

        rows.append(
            {
                "uniprot_id": r["uniprot_id"], "nr_code": r["nr_code"], "common_name": r["common_name"],
                "pdb_id": r["pdb_id"], "polymer_entity_id": r["polymer_entity_id"], "entity_id": r["entity_id"],
                "instance_id": r["instance_id"], "label_asym_id": r["label_asym_id"], "auth_asym_id": r["auth_asym_id"],
                "canonical_uniprot_position": r["canonical_uniprot_position"],
                "lbd_relative_position": r["lbd_relative_position"],
                "entity_label_seq_id": r["entity_label_seq_id"],
                "auth_seq_ids_json": r["auth_seq_ids_json"], "insertion_codes_json": r["insertion_codes_json"],
                "canonical_residue_name": canonical_residue_name(seq, r["canonical_uniprot_position"]),
                "observed_residue_names_json": r["label_comp_ids_json"],
                "is_construct_mapped": r["is_construct_mapped"],
                "ca_record_present": r["ca_record_present"], "ca_positive_occupancy": r["ca_positive_occupancy"],
                "ca_bfactor_usable": r["ca_bfactor_usable"], "zero_occupancy_only_ca": r["zero_occupancy_only_ca"],
                "altloc_ids_json": r["altloc_ids_json"], "occupancies_json": r["occupancies_json"],
                "b_iso_values_json": r["b_iso_values_json"],
                "selected_ca_altloc": selected["altloc"], "selected_ca_occupancy": selected["occupancy"],
                "selected_ca_b_iso_or_equiv": selected["b_iso_or_equiv"],
            }
        )

    residue_map = pd.DataFrame(rows)

    expected_rows = sum(int(r["standardized_lbd_length"]) for _, r in primary.iterrows())
    assert len(residue_map) == expected_rows, (
        f"expected {expected_rows} rows (sum of standardized_lbd_length across {len(primary)} primary structures), "
        f"got {len(residue_map)}"
    )

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    residue_map.to_parquet(RESIDUE_MAP_PARQUET, index=False)
    print(f"\nWrote {RESIDUE_MAP_PARQUET.relative_to(PROJECT_ROOT)} ({len(residue_map)} rows)")
    print(f"  usable selected Ca rows (selected_ca_occupancy notna): {residue_map['selected_ca_occupancy'].notna().sum()}")
    print(f"  rows with no usable Ca (selected fields null): {residue_map['selected_ca_occupancy'].isna().sum()}")


if __name__ == "__main__":
    main()
