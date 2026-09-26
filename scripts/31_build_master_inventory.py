"""
Final handoff revision: QC criteria are metadata/analysis flags, never
deletion criteria.

Builds, purely from already-frozen/derived local files (no network, no mmCIF
re-parse, no recomputation of Stage 2/3A/3B results):

  data/processed/structures.csv
      MASTER inventory: exactly one row per Stage 2 receptor x polymer-entity
      candidate (2072), with every QC value carried as a metadata/status flag.
  data/processed/lbd_residue_bfactors_all.parquet
      ALL corrected Stage 3B per-residue observations (every measured
      instance/model x every standardized-LBD position), with the
      deterministic raw Ca altloc handoff record. No filtering, no
      normalization.

The previous strict subset is preserved separately (structures_primary_qc_subset.csv,
final_lbd_residue_map_primary_qc_subset.parquet) and is used here only to
populate previous_primary_* lineage columns.
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

AUDIT_CSV = MANIFESTS_DIR / "final_structure_selection_audit.csv"
PRECOORD_CSV = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"
INSTANCE_QC_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_instance_qc.csv"
REP_CSV = MANIFESTS_DIR / "final_representative_instance_selection.csv"
LIGAND_CSV = MANIFESTS_DIR / "final_ligand_annotations.csv"
DOWNLOAD_MANIFEST_CSV = MANIFESTS_DIR / "mmcif_download_manifest.csv"
SUBSET_CSV = PROCESSED_DIR / "structures_primary_qc_subset.csv"

MASTER_CSV = PROCESSED_DIR / "structures.csv"
RESIDUE_ALL_PARQUET = PROCESSED_DIR / "lbd_residue_bfactors_all.parquet"

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL", "U": "SEC", "O": "PYL",
}

POLICY_NOTE = "previous strict primary-QC subset (final-structure-selection-v1)"


def build_master() -> pd.DataFrame:
    audit = pd.read_csv(AUDIT_CSV, dtype={"pdb_id": str})
    pq = pd.read_csv(PRECOORD_CSV, dtype={"pdb_id": str}).set_index(["uniprot_id", "polymer_entity_id"]).sort_index()
    rep = pd.read_csv(REP_CSV, dtype={"pdb_id": str}).set_index(["uniprot_id", "pdb_id", "polymer_entity_id"]).sort_index()
    iqc = pd.read_csv(INSTANCE_QC_CSV, dtype={"pdb_id": str}).set_index("instance_id")
    lig = pd.read_csv(LIGAND_CSV, dtype={"pdb_id": str}).set_index(["pdb_id", "uniprot_id", "polymer_entity_id"]).sort_index()
    dl = pd.read_csv(DOWNLOAD_MANIFEST_CSV, dtype={"pdb_id": str}).set_index("pdb_id")["sha256_gzip"].to_dict()

    rows = []
    for _, a in audit.iterrows():
        p = pq.loc[(a["uniprot_id"], a["polymer_entity_id"])]
        has_coord = bool(a["stage3b_coordinate_selected"])
        rep_key = (a["uniprot_id"], a["pdb_id"], a["polymer_entity_id"])
        r = rep.loc[rep_key] if has_coord else None
        inst = iqc.loc[r["selected_instance_id"]] if has_coord else None
        lig_key = (a["pdb_id"], a["uniprot_id"], a["polymer_entity_id"])
        has_lig = lig_key in lig.index
        assert has_lig == has_coord, f"ligand assessment must exist for exactly the coordinate candidates: {lig_key}"

        def cov_flag(v):
            return pd.NA if not has_coord else bool(v)

        reasons = json.loads(a["primary_exclusion_reasons_json"])
        prev_member = a["primary_selection_status"] == "SELECTED_PRIMARY"

        rows.append(
            {
                "uniprot_id": a["uniprot_id"], "nr_code": a["nr_code"], "common_name": a["common_name"], "group": a["group"],
                "pdb_id": a["pdb_id"], "polymer_entity_id": a["polymer_entity_id"], "entity_id": a["entity_id"],
                "experimental_method": a["experimental_methods"],
                "has_stage3b_coordinate_data": has_coord,
                "method_is_xray": bool(p["method_is_xray"]),
                "mapping_overlaps_lbd": bool(p["mapping_overlaps_lbd"]),
                "lbd_mapping_category": p["lbd_mapping_category"],
                "standardized_lbd_start": p["lbd_start"], "standardized_lbd_end": p["lbd_end"],
                "standardized_lbd_length": p["lbd_length"],
                "construct_lbd_mapping_coverage": p["lbd_mapping_coverage"],
                "positive_occupancy_lbd_coverage": a["positive_occupancy_lbd_coverage"] if has_coord else pd.NA,
                "bfactor_usable_lbd_coverage": a["bfactor_usable_lbd_coverage"] if has_coord else pd.NA,
                "coverage_ge_080": cov_flag(a["completeness_ge_080"]),
                "coverage_ge_090": cov_flag(a["completeness_ge_090"]),
                "coverage_ge_095": cov_flag(a["completeness_ge_095"]),
                "coverage_ge_099": cov_flag(a["completeness_ge_099"]),
                "coverage_eq_100": cov_flag(a["completeness_eq_100"]),
                "coordinate_data_status": "COORDINATE_DATA_AVAILABLE" if has_coord else "NO_STAGE3B_COORDINATE_DATA",
                "candidate_instance_count": int(r["candidate_instance_count"]) if has_coord else pd.NA,
                "resolution": a["resolution"], "resolution_project_range_status": p["resolution_project_range_status"],
                "r_free": a["r_free"], "r_free_status": p["r_free_status"],
                "identity_class": a["identity_class"], "fusion_or_chimera": bool(a["fusion_or_chimera"]),
                "taxonomy_audit_status": p["taxonomy_audit_status"] if isinstance(p["taxonomy_audit_status"], str) else "",
                "previous_primary_qc_member": prev_member,
                "previous_primary_exclusion_reasons_json": json.dumps(reasons),
                "representative_instance_id": r["selected_instance_id"] if has_coord else "",
                "label_asym_id": r["selected_label_asym_id"] if has_coord else "",
                "auth_asym_id": r["selected_auth_asym_id"] if has_coord else "",
                "ligand_state": lig.loc[lig_key, "ligand_state"] if has_lig else "NOT_ASSESSABLE_NO_COORDINATE_DATA",
                "ligand_annotation_status": (
                    lig.loc[lig_key, "ligand_annotation_status"] if has_lig else "NOT_ASSESSABLE_NO_COORDINATE_DATA"
                ),
                "functional_ligand_component_ids": lig.loc[lig_key, "functional_ligand_component_ids_json"] if has_lig else "",
                "mmcif_sha256": dl.get(a["pdb_id"], "") if has_coord else "",
                "in_completeness_95_sensitivity": bool(a["in_completeness_95_sensitivity"]),
                "in_completeness_99_sensitivity": bool(a["in_completeness_99_sensitivity"]),
                "in_completeness_80_sensitivity": bool(a["in_completeness_80_sensitivity"]),
                "in_ultrahigh_resolution_sensitivity": bool(a["in_ultrahigh_resolution_sensitivity"]),
                "in_no_fusion_sensitivity": bool(a["in_no_fusion_sensitivity"]),
            }
        )
    master = pd.DataFrame(rows)
    for c in ["coverage_ge_080", "coverage_ge_090", "coverage_ge_095", "coverage_ge_099", "coverage_eq_100"]:
        master[c] = master[c].astype("boolean")
    return master


def build_residue_all() -> pd.DataFrame:
    obs = pd.read_parquet(OBSERVATIONS_PARQUET)
    audit = pd.read_csv(AUDIT_CSV, dtype={"pdb_id": str})
    rep_ids = set(audit.loc[audit["representative_instance_id"].fillna("") != "", "representative_instance_id"])

    seqs = {}
    for uid in obs["uniprot_id"].unique():
        seqs[uid] = json.loads((UNIPROT_RAW_DIR / f"{uid}.json").read_text())["sequence"]["value"]

    sel_alt, sel_occ, sel_b = [], [], []
    for occ_j, b_j, alt_j in zip(obs["occupancies_json"], obs["b_iso_values_json"], obs["altloc_ids_json"]):
        s = select_representative_ca_record(json.loads(occ_j), json.loads(b_j), json.loads(alt_j))
        sel_alt.append(s["altloc"]); sel_occ.append(s["occupancy"]); sel_b.append(s["b_iso_or_equiv"])

    out = obs.drop(columns=["ca_modeled"]).copy()  # legacy alias of ca_record_present
    out["canonical_residue_name"] = [
        ONE_TO_THREE.get(seqs[u][p - 1], seqs[u][p - 1]) if 1 <= p <= len(seqs[u]) else ""
        for u, p in zip(out["uniprot_id"], out["canonical_uniprot_position"])
    ]
    out["is_representative_instance"] = out["instance_id"].isin(rep_ids)
    out["selected_ca_altloc"] = sel_alt
    out["selected_ca_occupancy"] = sel_occ
    out["selected_ca_b_iso_or_equiv"] = sel_b
    assert len(out) == len(obs)
    return out


def main() -> None:
    master = build_master()
    assert len(master) == 2072
    assert not master.duplicated(subset=["uniprot_id", "pdb_id", "polymer_entity_id"]).any()
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    master.to_csv(MASTER_CSV, index=False)
    print(f"Wrote {MASTER_CSV.relative_to(PROJECT_ROOT)} ({len(master)} rows)")
    print(f"  with coordinate data: {int(master['has_stage3b_coordinate_data'].sum())}")
    print(f"  previous primary QC members: {int(master['previous_primary_qc_member'].sum())}")

    subset = pd.read_csv(SUBSET_CSV, dtype={"pdb_id": str})
    m = master[master["previous_primary_qc_member"]]
    assert set(zip(subset["uniprot_id"], subset["pdb_id"], subset["polymer_entity_id"])) == set(
        zip(m["uniprot_id"], m["pdb_id"], m["polymer_entity_id"])
    ), "previous_primary_qc_member disagrees with preserved strict subset"

    res = build_residue_all()
    res.to_parquet(RESIDUE_ALL_PARQUET, index=False)
    print(f"Wrote {RESIDUE_ALL_PARQUET.relative_to(PROJECT_ROOT)} ({len(res)} rows, "
          f"{res['instance_id'].nunique()} instance/models)")

    v = master[master["uniprot_id"] == "P11473"]
    vr = res[res["uniprot_id"] == "P11473"]
    cov = v["positive_occupancy_lbd_coverage"].dropna().astype(float)
    print("\nVDR (P11473) verification:")
    print(f"  master candidates: {len(v)}; PDB entries: {v['pdb_id'].nunique()}; with coordinate data: {int(v['has_stage3b_coordinate_data'].sum())}")
    print(f"  residue rows: {len(vr)} across {vr['instance_id'].nunique()} instances / {vr['pdb_id'].nunique()} PDBs")
    print(f"  representative positive-occupancy coverage: best={cov.max():.4f} min={cov.min():.4f} median={cov.median():.4f}")


if __name__ == "__main__":
    main()
