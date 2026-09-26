"""
Final Data Release, Steps 5-9: enumerate every nonpolymer component and
instance in each PRIMARY-selected PDB entry, compute geometric contact to
each primary representative receptor's standardized-LBD atoms, classify
relevance, and derive APO/HOLO/AMBIGUOUS ligand state.

Water is enumerated once per PDB (not duplicated per receptor context,
since its classification never depends on geometry or receptor context —
see scripts/utils/nonpolymer_classification.py). All other nonpolymer
instances are evaluated once per PRIMARY representative-receptor context
in that PDB (38 primary-set PDB entries have two primary representative
receptors — e.g. RXR heterodimer partners — and the same physical ligand
can be the functional ligand of one receptor's LBD while being irrelevant
to the other's, so contact/classification is never collapsed across
receptor contexts).

Reads (no new network access — API responses already cached by
scripts/25_fetch_nonpolymer_metadata.py; mmCIF already cached locally):
    data/manifests/final_structure_selection_audit.csv
    data/interim/qc/stage3b_lbd_ca_observations.parquet
    data/raw/api/rcsb/data/nonpolymer_entities_batch_*.json
    data/raw/mmcif/{PDB_ID}.cif.gz

Writes:
    data/manifests/final_nonpolymer_inventory.csv
    data/manifests/final_ligand_annotations.csv
"""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.mmcif_parse import (  # noqa: E402
    get_covalent_linkages,
    get_entity_types,
    get_heavy_atoms,
    get_struct_site_records,
    parse_mmcif_gz,
)
from scripts.utils.nonpolymer_classification import (  # noqa: E402
    CONTACT_THRESHOLD_A,
    classify_nonpolymer_instance,
)

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"
RAW_API_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "data"
OBSERVATIONS_PARQUET = PROJECT_ROOT / "data" / "interim" / "qc" / "stage3b_lbd_ca_observations.parquet"

SELECTION_AUDIT_CSV = MANIFESTS_DIR / "final_structure_selection_audit.csv"
NONPOLYMER_INVENTORY_CSV = MANIFESTS_DIR / "final_nonpolymer_inventory.csv"
LIGAND_ANNOTATIONS_CSV = MANIFESTS_DIR / "final_ligand_annotations.csv"

WATER_COMP_IDS = {"HOH", "DOD"}
CONTACT_5_0_A = 5.0


def load_nonpolymer_entity_metadata() -> dict[str, dict]:
    """{'PDBID_N': {comp_id, name, formula, atom_count_heavy, subject_of_investigation, instances:[...]}}"""
    by_entity: dict[str, dict] = {}
    for path in sorted(RAW_API_DIR.glob("nonpolymer_entities_batch_*.json")):
        wrapper = json.loads(Path(path).read_text())
        for ent in wrapper["response"]["data"]["nonpolymer_entities"]:
            ids = ent["rcsb_nonpolymer_entity_container_identifiers"]
            cc = ent["nonpolymer_comp"]["chem_comp"] if ent["nonpolymer_comp"] else {}
            info = (ent["nonpolymer_comp"] or {}).get("rcsb_chem_comp_info") or {}
            annotations = ent["rcsb_nonpolymer_entity_annotation"] or []
            by_entity[ent["rcsb_id"]] = {
                "entry_id": ids["entry_id"], "entity_id": ids["entity_id"],
                "comp_id": ids["nonpolymer_comp_id"],
                "name": cc.get("name") or "", "formula": cc.get("formula") or "",
                "formula_weight": cc.get("formula_weight"), "chem_comp_type": cc.get("type") or "",
                "atom_count_heavy": info.get("atom_count_heavy"),
                "subject_of_investigation": any(a["type"] == "SUBJECT_OF_INVESTIGATION" for a in annotations),
                "annotation_types_json": json.dumps(sorted({a["type"] for a in annotations})),
                "instances": [
                    inst["rcsb_nonpolymer_entity_instance_container_identifiers"]
                    for inst in ent["nonpolymer_entity_instances"]
                ],
            }
    return by_entity


def build_lbd_atom_arrays(obs_rep: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, set[int]]:
    """From the representative instance's obs-parquet rows (already
    filtered to is_construct_mapped), return (coords[N,3], canonical_pos[N],
    entity_label_seq_id set) — placeholder; actual coords come from mmCIF,
    this just returns the label_seq_id -> canonical_position lookup."""
    mapped = obs_rep[obs_rep["is_construct_mapped"]]
    label_seq_ids = set(int(v) for v in mapped["entity_label_seq_id"].dropna())
    pos_by_label_seq = dict(zip(mapped["entity_label_seq_id"].astype("Int64"), mapped["canonical_uniprot_position"]))
    return label_seq_ids, pos_by_label_seq


def main() -> None:
    audit = pd.read_csv(SELECTION_AUDIT_CSV, dtype={"pdb_id": str})
    # Target population: ALL 1840 coordinate-assessable candidates; QC (coverage,
    # resolution, R-free, ...) is metadata only. Same classification rules as before.
    primary = audit[audit["stage3b_coordinate_selected"]].copy()
    print(f"Coordinate candidates: {len(primary)} across {primary['pdb_id'].nunique()} unique PDB entries")

    obs = pd.read_parquet(OBSERVATIONS_PARQUET)
    obs_by_instance = {iid: g for iid, g in obs.groupby("instance_id", sort=False)}

    nonpolymer_meta = load_nonpolymer_entity_metadata()
    print(f"Loaded RCSB metadata for {len(nonpolymer_meta)} non-polymer entities")

    inventory_rows: list[dict] = []
    ligand_rows: list[dict] = []

    by_pdb = primary.groupby("pdb_id")
    for pdb_num, (pdb_id, pdb_rows) in enumerate(by_pdb, start=1):
        mmcif_path = RAW_MMCIF_DIR / f"{pdb_id}.cif.gz"
        d = parse_mmcif_gz(mmcif_path)
        entity_types = get_entity_types(d)
        heavy_atoms = get_heavy_atoms(d)
        struct_sites = get_struct_site_records(d)
        covalent = get_covalent_linkages(d)

        # Group heavy atoms by label_asym_id, splitting nonpolymer vs water vs polymer.
        atoms_by_asym: dict[str, list[dict]] = {}
        for a in heavy_atoms:
            atoms_by_asym.setdefault(a["label_asym_id"], []).append(a)

        nonpolymer_asyms = []
        water_asyms = []
        for asym, atoms in atoms_by_asym.items():
            etype = entity_types.get(atoms[0]["label_entity_id"], "")
            if etype == "non-polymer":
                nonpolymer_asyms.append(asym)
            elif etype == "water":
                water_asyms.append(asym)

        # --- water: one row per PDB, no geometry ---
        water_by_seq: dict[tuple[str, str], list[dict]] = {}
        for asym in water_asyms:
            for a in atoms_by_asym[asym]:
                key = (asym, a["auth_seq_id"])
                water_by_seq.setdefault(key, []).append(a)
        for (asym, auth_seq), atoms in water_by_seq.items():
            inventory_rows.append(
                {
                    "pdb_id": pdb_id, "component_id": atoms[0]["label_comp_id"], "chemical_name": "water",
                    "formula": "H2 O", "formula_weight": 18.015, "chem_comp_type": "water",
                    "nonpolymer_entity_id": "", "nonpolymer_instance_id": f"{pdb_id}.{asym}",
                    "label_asym_id": asym, "auth_asym_id": atoms[0]["auth_asym_id"], "auth_seq_id": auth_seq,
                    "atom_count_heavy": 1, "representative_instance_id": "", "uniprot_id": "", "nr_code": "",
                    "min_heavy_atom_distance_to_lbd": "", "contact_within_4_5A": "", "contact_within_5_0A": "",
                    "number_of_lbd_residues_within_4_5A": "", "contact_lbd_canonical_positions_json": "[]",
                    "rcsb_subject_of_investigation": False, "struct_site_record": False, "covalently_linked": False,
                    "nonpolymer_relevance_class": "WATER", "nonpolymer_relevance_evidence": "Water molecule (HOH/DOD); excluded from geometric analysis by policy.",
                }
            )

        # --- nonpolymer (non-water): once per primary representative receptor context ---
        for _, prow in pdb_rows.iterrows():
            rep_instance_id = prow["representative_instance_id"]
            if rep_instance_id not in obs_by_instance:
                continue
            obs_rep = obs_by_instance[rep_instance_id]
            lbd_label_seq_ids, pos_by_label_seq = build_lbd_atom_arrays(obs_rep)

            rep_label_asym = prow["label_asym_id"]
            lbd_atoms = [
                a for a in atoms_by_asym.get(rep_label_asym, [])
                if a["label_seq_id"] in lbd_label_seq_ids and (a["occupancy"] or 0) > 0
            ]
            if lbd_atoms:
                lbd_coords = np.array([[a["x"], a["y"], a["z"]] for a in lbd_atoms])
                lbd_pos = np.array([pos_by_label_seq.get(a["label_seq_id"], -1) for a in lbd_atoms])
            else:
                lbd_coords = np.zeros((0, 3))
                lbd_pos = np.array([])

            no_lbd_atoms_note = ""
            this_structure_instances: list[dict] = []
            for asym in nonpolymer_asyms:
                atoms = atoms_by_asym[asym]
                comp_id = atoms[0]["label_comp_id"]
                entity_id = atoms[0]["label_entity_id"]
                nonpolymer_entity_id = f"{pdb_id}_{entity_id}"
                meta = nonpolymer_meta.get(nonpolymer_entity_id, {})

                pos_atoms = [a for a in atoms if (a["occupancy"] or 0) > 0]
                calc_atoms = pos_atoms if pos_atoms else atoms

                if len(lbd_coords) and calc_atoms:
                    inst_coords = np.array([[a["x"], a["y"], a["z"]] for a in calc_atoms])
                    dmat = cdist(inst_coords, lbd_coords)
                    min_dist = float(dmat.min())
                    within_4_5 = bool((dmat <= CONTACT_THRESHOLD_A).any())
                    within_5_0 = bool((dmat <= CONTACT_5_0_A).any())
                    contact_mask_per_lbd_atom = (dmat <= CONTACT_THRESHOLD_A).any(axis=0)
                    contact_positions = sorted({int(p) for p in lbd_pos[contact_mask_per_lbd_atom] if p >= 0})
                else:
                    min_dist, within_4_5, within_5_0, contact_positions = None, False, False, []

                auth_seq = atoms[0]["auth_seq_id"]
                auth_asym = atoms[0]["auth_asym_id"]
                has_struct_site = any(
                    s["auth_comp_id"] == comp_id and s["auth_asym_id"] == auth_asym and s["auth_seq_id"] == auth_seq
                    for s in struct_sites
                )
                is_covalent = any(
                    (c["ptnr1_label_asym_id"] == asym) or (c["ptnr2_label_asym_id"] == asym) for c in covalent
                )

                classification, evidence = classify_nonpolymer_instance(
                    comp_id=comp_id,
                    atom_count_heavy=meta.get("atom_count_heavy") or len({a["label_atom_id"] for a in atoms}),
                    contacts_lbd=within_4_5,
                    min_distance_to_lbd=min_dist,
                    has_subject_of_investigation_annotation=meta.get("subject_of_investigation", False),
                    has_struct_site_record=has_struct_site,
                    is_covalently_linked=is_covalent,
                )

                inventory_row = {
                    "pdb_id": pdb_id, "component_id": comp_id,
                    "chemical_name": meta.get("name", ""), "formula": meta.get("formula", ""),
                    "formula_weight": meta.get("formula_weight"), "chem_comp_type": meta.get("chem_comp_type", ""),
                    "nonpolymer_entity_id": nonpolymer_entity_id, "nonpolymer_instance_id": f"{pdb_id}.{asym}",
                    "label_asym_id": asym, "auth_asym_id": auth_asym, "auth_seq_id": auth_seq,
                    "atom_count_heavy": meta.get("atom_count_heavy"),
                    "representative_instance_id": rep_instance_id,
                    "uniprot_id": prow["uniprot_id"], "nr_code": prow["nr_code"],
                    "min_heavy_atom_distance_to_lbd": min_dist,
                    "contact_within_4_5A": within_4_5, "contact_within_5_0A": within_5_0,
                    "number_of_lbd_residues_within_4_5A": len(contact_positions),
                    "contact_lbd_canonical_positions_json": json.dumps(contact_positions),
                    "rcsb_subject_of_investigation": meta.get("subject_of_investigation", False),
                    "struct_site_record": has_struct_site, "covalently_linked": is_covalent,
                    "nonpolymer_relevance_class": classification, "nonpolymer_relevance_evidence": evidence,
                }
                inventory_rows.append(inventory_row)
                this_structure_instances.append(inventory_row)

            # --- aggregate to this representative structure's ligand state ---
            functional = [r for r in this_structure_instances if r["nonpolymer_relevance_class"] == "FUNCTIONAL_LBD_LIGAND"]
            possible_or_review = [
                r for r in this_structure_instances
                if r["nonpolymer_relevance_class"] in ("POSSIBLE_FUNCTIONAL_LIGAND", "REVIEW_NEEDED")
                and r["contact_within_4_5A"]
            ]
            all_lbd_contacting = [r for r in this_structure_instances if r["contact_within_4_5A"]]

            if not len(lbd_coords):
                # No observed LBD atoms to measure contact against: cannot safely say APO.
                ligand_state = "AMBIGUOUS"
                status = "REVIEW_NEEDED"
                no_lbd_atoms_note = "No positive-occupancy LBD heavy atoms available for contact geometry."
            elif functional:
                ligand_state = "HOLO"
                status = "RESOLVED"
            elif possible_or_review:
                ligand_state = "AMBIGUOUS"
                status = "REVIEW_NEEDED"
            else:
                ligand_state = "APO"
                status = "RESOLVED"

            ligand_rows.append(
                {
                    "pdb_id": pdb_id, "uniprot_id": prow["uniprot_id"], "nr_code": prow["nr_code"],
                    "common_name": prow["common_name"], "polymer_entity_id": prow["polymer_entity_id"],
                    "label_asym_id": rep_label_asym, "auth_asym_id": prow["auth_asym_id"],
                    "ligand_state": ligand_state,
                    "functional_ligand_component_ids_json": json.dumps(sorted({r["component_id"] for r in functional})),
                    "functional_ligand_names_json": json.dumps(sorted({r["chemical_name"] for r in functional if r["chemical_name"]})),
                    "all_lbd_contacting_nonpolymer_ids_json": json.dumps(sorted({r["nonpolymer_instance_id"] for r in all_lbd_contacting})),
                    "ligand_annotation_status": status,
                    "ligand_annotation_evidence": (
                        f"{len(functional)} FUNCTIONAL_LBD_LIGAND, {len(possible_or_review)} POSSIBLE/REVIEW "
                        f"LBD-contacting, {len(all_lbd_contacting)} total LBD-contacting nonpolymer instance(s)."
                    ),
                    "ligand_annotation_notes": (
                        no_lbd_atoms_note if not len(lbd_coords) else
                        "" if status == "RESOLVED" else
                        "One or more LBD-contacting nonpolymer components could not be confidently classified as "
                        "functional vs. incidental; see final_nonpolymer_inventory.csv for the specific instance(s) "
                        "and evidence."
                    ),
                }
            )

        if pdb_num % 200 == 0:
            print(f"  processed {pdb_num}/{len(by_pdb)} primary PDB entries")

    inventory_df = pd.DataFrame(inventory_rows)
    inventory_df.to_csv(NONPOLYMER_INVENTORY_CSV, index=False)
    print(f"\nWrote {NONPOLYMER_INVENTORY_CSV.relative_to(PROJECT_ROOT)} ({len(inventory_df)} rows)")
    print(inventory_df["nonpolymer_relevance_class"].value_counts().to_string())

    ligand_df = pd.DataFrame(ligand_rows)
    ligand_df.to_csv(LIGAND_ANNOTATIONS_CSV, index=False)
    print(f"\nWrote {LIGAND_ANNOTATIONS_CSV.relative_to(PROJECT_ROOT)} ({len(ligand_df)} rows)")
    print(ligand_df["ligand_state"].value_counts().to_string())


if __name__ == "__main__":
    main()
