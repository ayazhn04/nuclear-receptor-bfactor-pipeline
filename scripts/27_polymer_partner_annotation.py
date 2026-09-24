"""
Final Data Release, Step 10: identify polymeric receptor-interacting
partners (coregulator peptides, heterodimer receptor partners, other
protein chains) in each PRIMARY-selected PDB entry. Informational only —
does not alter primary structure selection.

Fetches (official RCSB Data API, cached/idempotent) metadata for every
polymer entity in each primary PDB entry OTHER than the project receptor's
own entity, then computes geometric contact to the representative
receptor's standardized-LBD atoms from the already-cached mmCIF files.

Reads:
    data/manifests/final_structure_selection_audit.csv
    data/interim/qc/stage3b_lbd_ca_observations.parquet
    data/raw/api/rcsb/data/entry_container_identifiers_batch_*.json  (from script 25)
    data/raw/mmcif/{PDB_ID}.cif.gz

Writes (raw, cached):
    data/raw/api/rcsb/data/other_polymer_entities_batch_XXXX.json

Writes (tracked):
    data/manifests/final_polymer_partner_annotations.csv
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402
from scripts.utils.mmcif_parse import get_entity_types, get_heavy_atoms, parse_mmcif_gz  # noqa: E402
from scripts.utils.rcsb_data_client import DEFAULT_BATCH_SIZE, fetch_batches  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
RAW_API_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "data"
RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"
OBSERVATIONS_PARQUET = PROJECT_ROOT / "data" / "interim" / "qc" / "stage3b_lbd_ca_observations.parquet"

SELECTION_AUDIT_CSV = MANIFESTS_DIR / "final_structure_selection_audit.csv"
PARTNER_ANNOTATIONS_CSV = MANIFESTS_DIR / "final_polymer_partner_annotations.csv"

CONTACT_THRESHOLD_A = 5.0
COREGULATOR_KEYWORDS = [
    "coactivator", "corepressor", "coregulator", "nuclear receptor box", "lxxll",
    "src-1", "src-2", "src-3", "src1", "src2", "src3", "ncoa", "ncor", "smrt",
    "pgc-1", "pgc1", "steroid receptor coactivator", "tif2", "tif-2", "grip1", "grip-1",
    "actr", "med1", "mediator", "peroxisome proliferator-activated receptor gamma coactivator",
    "silencing mediator",
]
COREGULATOR_MAX_LENGTH = 40

PROJECT_RECEPTOR_UNIPROT_IDS = {r["uniprot_id"] for r in NR_METADATA}
PROJECT_RECEPTOR_NAMES = {r["uniprot_id"]: r["common_name"] for r in NR_METADATA}


def cached_ids(prefix: str) -> set[str]:
    covered: set[str] = set()
    for path in sorted(RAW_API_DIR.glob(f"{prefix}_batch_*.json")):
        wrapper = json.loads(Path(path).read_text())
        covered.update(wrapper.get("ids_requested", []))
    return covered


def fetch_and_cache(prefix: str, batch_type: str, ids: list[str]) -> None:
    already = cached_ids(prefix)
    missing = [i for i in ids if i not in already]
    if not missing:
        print(f"  {prefix}: all {len(ids)} IDs already cached — reusing.")
        return
    print(f"  {prefix}: fetching {len(missing)} of {len(ids)} IDs...")
    start_number = len(sorted(RAW_API_DIR.glob(f"{prefix}_batch_*.json")))
    results = fetch_batches(batch_type, missing, batch_size=DEFAULT_BATCH_SIZE)
    RAW_API_DIR.mkdir(parents=True, exist_ok=True)
    for r in results:
        batch_number = start_number + r.batch_number
        raw_path = RAW_API_DIR / f"{prefix}_batch_{batch_number:04d}.json"
        wrapper = {
            "batch_type": batch_type, "batch_number": batch_number, "ids_requested": r.ids_requested,
            "retrieved_at_utc": r.retrieved_at_utc, "http_status": r.http_status,
            "success": r.success, "error_message": r.error_message, "response": r.raw_json,
        }
        raw_path.write_bytes(json.dumps(wrapper, indent=2, ensure_ascii=False).encode("utf-8"))
        print(f"    batch {batch_number:4d} ({len(r.ids_requested):4d} ids) -> {'OK' if r.success else 'FAILED'}")
        if not r.success:
            raise RuntimeError(f"STOPPING: {prefix} batch {batch_number} failed: {r.error_message}")


def load_polymer_entity_ids_by_pdb(pdb_ids: list[str]) -> dict[str, list[str]]:
    by_pdb: dict[str, list[str]] = {}
    for path in sorted(RAW_API_DIR.glob("entry_container_identifiers_batch_*.json")):
        wrapper = json.loads(Path(path).read_text())
        for entry in wrapper["response"]["data"]["entries"]:
            pdb_id = entry["rcsb_id"]
            ids = entry["rcsb_entry_container_identifiers"].get("polymer_entity_ids") or []
            by_pdb[pdb_id] = [f"{pdb_id}_{eid}" for eid in ids]
    return {pdb_id: by_pdb.get(pdb_id, []) for pdb_id in pdb_ids}


def load_other_polymer_entity_metadata() -> dict[str, dict]:
    meta: dict[str, dict] = {}
    for path in sorted(RAW_API_DIR.glob("other_polymer_entities_batch_*.json")):
        wrapper = json.loads(Path(path).read_text())
        for ent in wrapper["response"]["data"]["polymer_entities"]:
            ids = ent["rcsb_polymer_entity_container_identifiers"]
            refs = ids.get("reference_sequence_identifiers") or []
            uniprot_accessions = [r["database_accession"] for r in refs if r["database_name"] == "UniProt"]
            meta[ent["rcsb_id"]] = {
                "entry_id": ids["entry_id"], "entity_id": ids["entity_id"],
                "asym_ids": ids.get("asym_ids") or [], "auth_asym_ids": ids.get("auth_asym_ids") or [],
                "description": (ent.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "",
                "sequence_length": (ent.get("entity_poly") or {}).get("rcsb_sample_sequence_length"),
                "uniprot_accessions": uniprot_accessions,
            }
    return meta


def classify_partner(
    description: str, sequence_length: int | None, uniprot_accessions: list[str], own_receptor_uniprot_id: str
) -> tuple[str, str]:
    matched_receptors = [u for u in uniprot_accessions if u in PROJECT_RECEPTOR_UNIPROT_IDS]
    if matched_receptors:
        names = ", ".join(f"{u} ({PROJECT_RECEPTOR_NAMES[u]})" for u in matched_receptors)
        if matched_receptors == [own_receptor_uniprot_id]:
            relation = "same receptor — a second copy/homodimer partner chain, not a heterodimer partner"
        else:
            relation = "a different project-receptor heterodimer partner"
        return "OTHER_RECEPTOR_PARTNER", f"UniProt mapping matches project receptor(s): {names} ({relation})."

    desc_lower = (description or "").lower()
    matched_kw = [kw for kw in COREGULATOR_KEYWORDS if kw in desc_lower]
    if matched_kw and sequence_length is not None and sequence_length <= COREGULATOR_MAX_LENGTH:
        return "COREGULATOR_CANDIDATE", f"Short polymer ({sequence_length} aa) with coregulator-keyword description match: {matched_kw}."
    if matched_kw:
        return "COREGULATOR_CANDIDATE", f"Description keyword match ({matched_kw}) despite length {sequence_length} aa."

    if sequence_length is not None and sequence_length <= COREGULATOR_MAX_LENGTH:
        return "UNCLASSIFIED_PARTNER", f"Short polymer ({sequence_length} aa) but no coregulator-keyword evidence in description ('{description}') — not assumed to be a coregulator without supporting metadata."

    return "OTHER_PROTEIN_PARTNER", f"Longer polymer ({sequence_length} aa), no project-receptor UniProt match, no coregulator keyword: '{description}'."


def main() -> None:
    audit = pd.read_csv(SELECTION_AUDIT_CSV, dtype={"pdb_id": str})
    primary = audit[audit["in_primary_set"]].copy()
    primary_pdb_ids = sorted(primary["pdb_id"].unique().tolist())

    print("Enumerating all polymer entities per primary PDB (via cached entry_container_identifiers)...")
    all_entities_by_pdb = load_polymer_entity_ids_by_pdb(primary_pdb_ids)

    own_entity_ids_by_pdb: dict[str, set[str]] = {}
    for pdb_id, rows in primary.groupby("pdb_id"):
        own_entity_ids_by_pdb[pdb_id] = {f"{pdb_id}_{eid}" for eid in rows["entity_id"].astype(str)}

    other_entity_ids = sorted(
        {
            eid
            for pdb_id, ids in all_entities_by_pdb.items()
            for eid in ids
            if eid not in own_entity_ids_by_pdb.get(pdb_id, set())
        }
    )
    print(f"Non-project-receptor polymer entities across primary set: {len(other_entity_ids)}")

    if other_entity_ids:
        fetch_and_cache("other_polymer_entities", "polymer_entity", other_entity_ids)
    other_meta = load_other_polymer_entity_metadata()
    print(f"Loaded metadata for {len(other_meta)} other polymer entities")

    obs = pd.read_parquet(OBSERVATIONS_PARQUET)
    obs_by_instance = {iid: g for iid, g in obs.groupby("instance_id", sort=False)}

    rows_out: list[dict] = []
    for pdb_num, (pdb_id, pdb_rows) in enumerate(primary.groupby("pdb_id"), start=1):
        other_ids_here = [eid for eid in all_entities_by_pdb.get(pdb_id, []) if eid not in own_entity_ids_by_pdb.get(pdb_id, set())]
        if not other_ids_here:
            continue

        d = parse_mmcif_gz(RAW_MMCIF_DIR / f"{pdb_id}.cif.gz")
        entity_types = get_entity_types(d)
        heavy_atoms = get_heavy_atoms(d)
        atoms_by_asym: dict[str, list[dict]] = {}
        for a in heavy_atoms:
            atoms_by_asym.setdefault(a["label_asym_id"], []).append(a)

        for _, prow in pdb_rows.iterrows():
            rep_instance_id = prow["representative_instance_id"]
            if rep_instance_id not in obs_by_instance:
                continue
            obs_rep = obs_by_instance[rep_instance_id]
            mapped = obs_rep[obs_rep["is_construct_mapped"]]
            lbd_label_seq_ids = set(int(v) for v in mapped["entity_label_seq_id"].dropna())
            pos_by_label_seq = dict(zip(mapped["entity_label_seq_id"].astype("Int64"), mapped["canonical_uniprot_position"]))

            rep_label_asym = prow["label_asym_id"]
            lbd_atoms = [
                a for a in atoms_by_asym.get(rep_label_asym, [])
                if a["label_seq_id"] in lbd_label_seq_ids and (a["occupancy"] or 0) > 0
            ]
            if not lbd_atoms:
                continue
            lbd_coords = np.array([[a["x"], a["y"], a["z"]] for a in lbd_atoms])
            lbd_pos = np.array([pos_by_label_seq.get(a["label_seq_id"], -1) for a in lbd_atoms])

            for entity_full_id in other_ids_here:
                entity_id = entity_full_id.split("_", 1)[1]
                meta = other_meta.get(entity_full_id, {})
                if entity_types.get(entity_id) != "polymer":
                    continue

                partner_asyms = [asym for asym, atoms in atoms_by_asym.items() if atoms[0]["label_entity_id"] == entity_id]
                for asym in partner_asyms:
                    atoms = [a for a in atoms_by_asym[asym] if (a["occupancy"] or 0) > 0]
                    if not atoms:
                        continue
                    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms])
                    dmat = cdist(coords, lbd_coords)
                    min_dist = float(dmat.min())
                    contact_mask = (dmat <= CONTACT_THRESHOLD_A).any(axis=0)
                    contact_positions = sorted({int(p) for p in lbd_pos[contact_mask] if p >= 0})

                    if min_dist > CONTACT_THRESHOLD_A:
                        continue  # not contacting the LBD — not recorded (informational table, LBD-relevant only)

                    classification, reason = classify_partner(
                        meta.get("description", ""), meta.get("sequence_length"), meta.get("uniprot_accessions", []),
                        prow["uniprot_id"],
                    )
                    rows_out.append(
                        {
                            "pdb_id": pdb_id, "receptor_uniprot_id": prow["uniprot_id"], "receptor_nr_code": prow["nr_code"],
                            "receptor_common_name": prow["common_name"], "receptor_representative_instance_id": rep_instance_id,
                            "partner_polymer_entity_id": entity_full_id, "partner_chain": asym,
                            "partner_auth_asym_id": atoms[0]["auth_asym_id"],
                            "partner_description": meta.get("description", ""),
                            "partner_sequence_length": meta.get("sequence_length"),
                            "partner_mapped_uniprot_ids_json": json.dumps(meta.get("uniprot_accessions", [])),
                            "minimum_distance_to_receptor_lbd": min_dist,
                            "contacting_receptor_canonical_positions_json": json.dumps(contact_positions),
                            "partner_category": classification, "partner_category_evidence": reason,
                        }
                    )

        if pdb_num % 300 == 0:
            print(f"  processed {pdb_num} primary PDB entries")

    out_df = pd.DataFrame(rows_out)
    out_df.to_csv(PARTNER_ANNOTATIONS_CSV, index=False)
    print(f"\nWrote {PARTNER_ANNOTATIONS_CSV.relative_to(PROJECT_ROOT)} ({len(out_df)} rows)")
    if len(out_df):
        print(out_df["partner_category"].value_counts().to_string())


if __name__ == "__main__":
    main()
