"""
Part 3 of the current functional-site flexibility stage (continuing from
commit fb27a86 on analysis/normalization-v1): compute the ligand-binding
pocket (LBP) as a receptor-level CONSENSUS across holo structures, from
real atomic geometry -- never a single universal UniProt residue list
carried over from one structure or from literature.

For every receptor with >= MIN_HOLO_PDBS_FOR_POCKET holo PDB entries:
1. the biologically relevant bound ligand per entry is taken from the
   ALREADY-COMPUTED, already-validated nonpolymer relevance classification
   (data/manifests/final_nonpolymer_inventory.csv,
   nonpolymer_relevance_class == "FUNCTIONAL_LBD_LIGAND" -- this
   classification already excludes crystallization additives, ions,
   buffer molecules, and other irrelevant HETATM records; see
   scripts/utils/nonpolymer_classification.py, unmodified, reused
   read-only here);
2. for that entry, every receptor heavy atom (occupancy > 0) within 5.0 A
   of any heavy atom of that functional ligand is identified by re-parsing
   the cached mmCIF file (the master ligand-state pipeline only stored
   exact contact positions at 4.5 A plus a boolean at 5.0 A -- this script
   recomputes the position list at exactly the 5.0 A the current
   instruction specifies, reusing the same parsing utilities);
3. contacts are mapped to canonical UniProt positions and aggregated
   across all eligible holo entries into a per-residue contact_fraction.

This is a genuinely new, compute-heavy pass over cached mmCIF files (no
network access) -- separated into its own script because of that cost.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.mmcif_parse import get_entity_types, get_heavy_atoms, parse_mmcif_gz  # noqa: E402

MANIFESTS = PROJECT_ROOT / "data" / "manifests"
PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"

CONTACT_5_0_A = 5.0
MIN_HOLO_PDBS_FOR_POCKET = 3  # documented minimum for a "consensus" -- fixed before computing


def build_lbd_position_lookup() -> tuple[dict, dict]:
    """instance_id -> {label_seq_id: canonical_uniprot_position}, and
    instance_id -> label_asym_id, both from the protected, read-only master
    residue table (never modified)."""
    residues = pd.read_parquet(PROCESSED / "lbd_residue_bfactors_all.parquet")
    rep_all = residues[residues["is_representative_instance"]]
    asym_lookup = rep_all.groupby("instance_id")["label_asym_id"].first().to_dict()
    rep = rep_all[rep_all["is_construct_mapped"]]
    pos_lookup = {}
    for iid, g in rep.groupby("instance_id"):
        pos_lookup[iid] = dict(zip(g["entity_label_seq_id"].astype("Int64"), g["canonical_uniprot_position"]))
    return pos_lookup, asym_lookup


def main():
    inventory = pd.read_csv(MANIFESTS / "final_nonpolymer_inventory.csv", dtype={"pdb_id": str}, keep_default_na=False)
    functional = inventory[inventory["nonpolymer_relevance_class"] == "FUNCTIONAL_LBD_LIGAND"].copy()
    print(f"{len(functional)} FUNCTIONAL_LBD_LIGAND instances across "
          f"{functional['pdb_id'].nunique()} PDB entries, {functional['uniprot_id'].nunique()} receptors "
          "(classification reused read-only from the master ligand-state pipeline)")

    rfs = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")
    eligible_receptors = []
    per_receptor_holo_pdbs = functional.groupby("uniprot_id")["pdb_id"].nunique()
    for uid, n in per_receptor_holo_pdbs.items():
        if n >= MIN_HOLO_PDBS_FOR_POCKET:
            eligible_receptors.append(uid)
    print(f"\n{len(eligible_receptors)} receptors have >= {MIN_HOLO_PDBS_FOR_POCKET} holo PDB entries "
          "with a classified functional ligand -- eligible for a consensus LBP")

    lbd_lookup, asym_lookup = build_lbd_position_lookup()

    contact_rows = []  # per (pdb_id, uniprot_id, canonical_position): contacted True/False
    processed_pdbs = 0
    for pdb_id, g in functional[functional["uniprot_id"].isin(eligible_receptors)].groupby("pdb_id"):
        mmcif_path = RAW_MMCIF_DIR / f"{pdb_id}.cif.gz"
        if not mmcif_path.exists():
            print(f"  WARNING: {pdb_id}.cif.gz not cached -- skipping (no network access used)")
            continue
        d = parse_mmcif_gz(mmcif_path)
        heavy_atoms = get_heavy_atoms(d)
        atoms_by_asym: dict[str, list[dict]] = {}
        for a in heavy_atoms:
            atoms_by_asym.setdefault(a["label_asym_id"], []).append(a)

        for _, row in g.iterrows():
            uid, rep_instance_id, ligand_asym = row["uniprot_id"], row["representative_instance_id"], row["label_asym_id"]
            if rep_instance_id not in lbd_lookup:
                continue
            pos_by_label_seq = lbd_lookup[rep_instance_id]
            recv_label_asym = asym_lookup.get(rep_instance_id)
            if recv_label_asym is None:
                continue

            lbd_atoms = [
                a for a in atoms_by_asym.get(recv_label_asym, [])
                if a["label_seq_id"] in pos_by_label_seq and (a["occupancy"] or 0) > 0
            ]
            ligand_atoms = [a for a in atoms_by_asym.get(ligand_asym, []) if (a["occupancy"] or 0) > 0]
            if not lbd_atoms or not ligand_atoms:
                continue

            lbd_coords = np.array([[a["x"], a["y"], a["z"]] for a in lbd_atoms])
            lbd_pos = np.array([pos_by_label_seq.get(a["label_seq_id"], -1) for a in lbd_atoms])
            lig_coords = np.array([[a["x"], a["y"], a["z"]] for a in ligand_atoms])

            dmat = cdist(lig_coords, lbd_coords)
            contact_mask = (dmat <= CONTACT_5_0_A).any(axis=0)
            contact_positions = set(int(p) for p in lbd_pos[contact_mask] if p >= 0)
            all_positions = set(int(p) for p in np.unique(lbd_pos) if p >= 0)

            for pos in all_positions:
                contact_rows.append({
                    "uniprot_id": uid, "pdb_id": pdb_id,
                    "canonical_uniprot_position": pos, "contacts_ligand_5_0A": pos in contact_positions,
                })
        processed_pdbs += 1
        if processed_pdbs % 100 == 0:
            print(f"  processed {processed_pdbs} PDB entries...")

    contacts = pd.DataFrame(contact_rows)
    raw_out = ANALYSIS / "ligand_pocket_raw_contacts.parquet"
    contacts.to_parquet(raw_out, index=False)
    print(f"\nWrote {raw_out.relative_to(PROJECT_ROOT)} ({len(contacts)} rows, "
          f"{contacts['pdb_id'].nunique()} PDB entries processed)")

    # aggregate to receptor x position consensus
    agg_rows = []
    for (uid, pos), g in contacts.groupby(["uniprot_id", "canonical_uniprot_position"]):
        n_observed = g["pdb_id"].nunique()
        n_contacting = g.loc[g["contacts_ligand_5_0A"], "pdb_id"].nunique()
        agg_rows.append({
            "uniprot_id": uid, "canonical_uniprot_position": pos,
            "n_holo_pdbs_observed": n_observed, "n_holo_pdbs_contacting": n_contacting,
            "contact_fraction": n_contacting / n_observed if n_observed else np.nan,
        })
    consensus = pd.DataFrame(agg_rows).sort_values(["uniprot_id", "canonical_uniprot_position"]).reset_index(drop=True)
    consensus["in_pocket_ge_10pct"] = consensus["contact_fraction"] >= 0.10
    consensus["in_pocket_ge_25pct_primary"] = consensus["contact_fraction"] >= 0.25
    consensus["in_pocket_ge_50pct"] = consensus["contact_fraction"] >= 0.50

    out = ANALYSIS / "ligand_pocket_consensus_residues.parquet"
    consensus.to_parquet(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(consensus)} candidate residues)")

    for col, label in [("in_pocket_ge_10pct", ">=10%"), ("in_pocket_ge_25pct_primary", ">=25% (PRIMARY)"),
                        ("in_pocket_ge_50pct", ">=50%")]:
        n_res = consensus[col].sum()
        n_rec = consensus.loc[consensus[col], "uniprot_id"].nunique()
        print(f"  {label}: {n_res} pocket residues across {n_rec} receptors")

    n_apo_only = len(rfs) - len(eligible_receptors)  # descriptive only, not used downstream
    print(f"\n{len(eligible_receptors)} receptors got a consensus LBP; the remainder "
          "(including apo-only or holo-poor receptors) are marked pocket-UNRESOLVED downstream "
          "-- no homologous pocket is invented for them.")


if __name__ == "__main__":
    main()
