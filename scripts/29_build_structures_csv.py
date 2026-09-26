"""
Final Data Release, Steps 13-16: build the definitive team handoff tables
structures.csv (primary only) and excluded_structures.csv (full Stage 2
lineage minus primary), run primary-dataset validation assertions, and
build the 48-row final_counts_by_receptor.csv.

Reads:
    data/manifests/final_structure_selection_audit.csv
    data/manifests/final_representative_instance_selection.csv
    data/manifests/final_ligand_annotations.csv
    data/manifests/final_polymer_partner_annotations.csv
    data/manifests/mmcif_download_manifest.csv
    data/manifests/rcsb_candidate_inventory.csv
    data/manifests/rcsb_precoordinate_qc.csv
    data/interim/qc/stage3b_lbd_ca_observations.parquet
    data/raw/mmcif/{PDB_ID}.cif.gz   (anisotropic-category flag only)

Writes:
    data/processed/structures_primary_qc_subset.csv   (historical strict subset; master structures.csv is built by script 31)
    data/processed/excluded_structures.csv
    data/manifests/final_counts_by_receptor.csv
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402
from scripts.utils.missing_run_analysis import analyze_missing_runs  # noqa: E402
from scripts.utils.mmcif_parse import get_entry_metadata, parse_mmcif_gz  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"
OBSERVATIONS_PARQUET = PROJECT_ROOT / "data" / "interim" / "qc" / "stage3b_lbd_ca_observations.parquet"

SELECTION_AUDIT_CSV = MANIFESTS_DIR / "final_structure_selection_audit.csv"
REP_SELECTION_CSV = MANIFESTS_DIR / "final_representative_instance_selection.csv"
LIGAND_ANNOTATIONS_CSV = MANIFESTS_DIR / "final_ligand_annotations.csv"
PARTNER_ANNOTATIONS_CSV = MANIFESTS_DIR / "final_polymer_partner_annotations.csv"
DOWNLOAD_MANIFEST_CSV = MANIFESTS_DIR / "mmcif_download_manifest.csv"
CANDIDATE_INVENTORY_CSV = MANIFESTS_DIR / "rcsb_candidate_inventory.csv"
PRECOORDINATE_QC_CSV = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"
INSTANCE_QC_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_instance_qc.csv"
LBD_REFERENCE_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"

STRUCTURES_CSV = PROCESSED_DIR / "structures_primary_qc_subset.csv"
EXCLUDED_STRUCTURES_CSV = PROCESSED_DIR / "excluded_structures.csv"
FINAL_COUNTS_BY_RECEPTOR_CSV = MANIFESTS_DIR / "final_counts_by_receptor.csv"

SELECTION_POLICY_VERSION = "final-structure-selection-v1"


def load_anisotropic_flags(pdb_ids: list[str]) -> dict[str, bool]:
    flags = {}
    for pdb_id in pdb_ids:
        d = parse_mmcif_gz(RAW_MMCIF_DIR / f"{pdb_id}.cif.gz")
        meta = get_entry_metadata(d)
        flags[pdb_id] = meta["anisotrop_category_present"]
    return flags


def main() -> None:
    audit = pd.read_csv(SELECTION_AUDIT_CSV, dtype={"pdb_id": str})
    primary = audit[audit["in_primary_set"]].copy()
    excluded = audit[~audit["in_primary_set"]].copy()

    rep_selection = pd.read_csv(REP_SELECTION_CSV, dtype={"pdb_id": str}).set_index(
        ["uniprot_id", "pdb_id", "polymer_entity_id"]
    ).sort_index()
    ligand = pd.read_csv(LIGAND_ANNOTATIONS_CSV, dtype={"pdb_id": str}).set_index(
        ["pdb_id", "uniprot_id", "polymer_entity_id"]
    ).sort_index()
    partners = pd.read_csv(PARTNER_ANNOTATIONS_CSV, dtype={"pdb_id": str}) if PARTNER_ANNOTATIONS_CSV.exists() else pd.DataFrame(
        columns=["pdb_id", "receptor_uniprot_id", "partner_category", "partner_description"]
    )
    downloads = pd.read_csv(DOWNLOAD_MANIFEST_CSV, dtype={"pdb_id": str}).set_index("pdb_id")["sha256_gzip"].to_dict()
    inventory = pd.read_csv(CANDIDATE_INVENTORY_CSV, dtype={"pdb_id": str})[
        ["pdb_id", "polymer_entity_id", "source_organism_names"]
    ].drop_duplicates(subset=["pdb_id", "polymer_entity_id"]).set_index(["pdb_id", "polymer_entity_id"])
    precoordinate_qc = pd.read_csv(PRECOORDINATE_QC_CSV, dtype={"pdb_id": str}).set_index(["uniprot_id", "polymer_entity_id"]).sort_index()
    instance_qc = pd.read_csv(INSTANCE_QC_CSV, dtype={"pdb_id": str}).set_index("instance_id")
    lbd_reference = pd.read_csv(LBD_REFERENCE_CSV)
    lbd_lengths = dict(zip(lbd_reference["uniprot_id"], lbd_reference["lbd_length"]))

    obs = pd.read_parquet(OBSERVATIONS_PARQUET)

    print("Computing internal-vs-terminal missing-run stats for primary representative instances...")
    rep_instance_ids = set(primary["representative_instance_id"])
    run_stats: dict[str, dict] = {}
    for instance_id, group in obs[obs["instance_id"].isin(rep_instance_ids)].groupby("instance_id", sort=False):
        uniprot_id = group["uniprot_id"].iloc[0]
        run_stats[instance_id] = analyze_missing_runs(group, lbd_lengths[uniprot_id])

    print("Loading anisotropic-refinement flags for all primary PDB entries...")
    anisotropic_flags = load_anisotropic_flags(sorted(primary["pdb_id"].unique()))

    print("\nBuilding structures.csv (primary only)...")
    partner_by_key = partners.groupby(["pdb_id", "receptor_uniprot_id"]) if len(partners) else None

    structure_rows = []
    for _, r in primary.iterrows():
        rep_key = (r["uniprot_id"], r["pdb_id"], r["polymer_entity_id"])
        rep = rep_selection.loc[rep_key]
        lig = ligand.loc[(r["pdb_id"], r["uniprot_id"], r["polymer_entity_id"])]
        pq = precoordinate_qc.loc[(r["uniprot_id"], r["polymer_entity_id"])]
        org = inventory.loc[(r["pdb_id"], r["polymer_entity_id"]), "source_organism_names"] if (r["pdb_id"], r["polymer_entity_id"]) in inventory.index else "[]"
        run = run_stats.get(r["representative_instance_id"], {})

        if partner_by_key is not None and (r["pdb_id"], r["uniprot_id"]) in partner_by_key.groups:
            p_rows = partner_by_key.get_group((r["pdb_id"], r["uniprot_id"]))
            coregulator_or_partner_present = bool(len(p_rows))
            polymer_partner_summary = "; ".join(
                f"{cat} ({desc})" for cat, desc in zip(p_rows["partner_category"], p_rows["partner_description"])
            )
        else:
            coregulator_or_partner_present = False
            polymer_partner_summary = ""

        structure_rows.append(
            {
                "uniprot_id": r["uniprot_id"], "nr_code": r["nr_code"], "common_name": r["common_name"], "group": r["group"],
                "pdb_id": r["pdb_id"], "polymer_entity_id": r["polymer_entity_id"], "entity_id": r["entity_id"],
                "instance_id": r["representative_instance_id"], "label_asym_id": r["label_asym_id"], "auth_asym_id": r["auth_asym_id"],
                "experimental_method": r["experimental_methods"], "resolution": r["resolution"], "r_free": r["r_free"],
                "standardized_lbd_start": r["standardized_lbd_start"], "standardized_lbd_end": r["standardized_lbd_end"],
                "standardized_lbd_length": r["standardized_lbd_length"],
                "construct_lbd_mapping_coverage": pq["lbd_mapping_coverage"],
                "positive_occupancy_lbd_coverage": r["positive_occupancy_lbd_coverage"],
                "bfactor_usable_lbd_coverage": r["bfactor_usable_lbd_coverage"],
                "mapped_lbd_residues": pq["mapped_lbd_residues"],
                "positive_occupancy_lbd_residues": round(r["positive_occupancy_lbd_coverage"] * r["standardized_lbd_length"]),
                "mapped_but_unobserved_lbd_residues": run.get("mapped_but_unmodeled_residue_count"),
                "zero_occupancy_only_lbd_residues": instance_qc.loc[r["representative_instance_id"], "mapped_with_zero_occupancy_only_ca_residues"],
                "internal_missing_run_count": run.get("internal_missing_run_count"),
                "longest_internal_missing_run": run.get("longest_internal_missing_run"),
                "anisotropic_category_present": anisotropic_flags.get(r["pdb_id"], False),
                "identity_class": r["identity_class"], "fusion_or_chimera": r["fusion_or_chimera"],
                "source_organism": org, "taxonomy_audit_status": pq["taxonomy_audit_status"],
                "ligand_state": lig["ligand_state"],
                "functional_ligand_component_ids": lig["functional_ligand_component_ids_json"],
                "functional_ligand_names": lig["functional_ligand_names_json"],
                "coregulator_or_partner_present": coregulator_or_partner_present,
                "polymer_partner_summary": polymer_partner_summary,
                "passes_95pct_completeness": r["in_completeness_95_sensitivity"],
                "passes_99pct_completeness": r["in_completeness_99_sensitivity"],
                "passes_no_fusion_sensitivity": r["in_no_fusion_sensitivity"],
                "is_ultrahigh_resolution": r["in_ultrahigh_resolution_sensitivity"],
                "mmcif_sha256": downloads.get(r["pdb_id"], ""),
                "selection_policy_version": SELECTION_POLICY_VERSION,
            }
        )

    structures_df = pd.DataFrame(structure_rows)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    structures_df.to_csv(STRUCTURES_CSV, index=False)
    print(f"Wrote {STRUCTURES_CSV.relative_to(PROJECT_ROOT)} ({len(structures_df)} rows)")

    print("\nBuilding excluded_structures.csv...")
    excluded_cols = [
        "uniprot_id", "nr_code", "common_name", "group", "pdb_id", "polymer_entity_id", "entity_id",
        "stage3b_coordinate_selected", "representative_instance_id", "label_asym_id", "auth_asym_id",
        "experimental_methods", "resolution", "r_free",
        "standardized_lbd_start", "standardized_lbd_end", "standardized_lbd_length",
        "positive_occupancy_lbd_coverage", "bfactor_usable_lbd_coverage",
        "completeness_ge_090", "completeness_ge_095", "completeness_ge_099", "completeness_ge_080", "completeness_eq_100",
        "identity_class", "fusion_or_chimera",
        "primary_method_pass", "primary_completeness_pass", "primary_resolution_pass", "primary_rfree_pass", "primary_identity_pass",
        "primary_selection_status", "primary_exclusion_reasons_json",
    ]
    excluded_df = excluded[excluded_cols].copy()
    excluded_df.to_csv(EXCLUDED_STRUCTURES_CSV, index=False)
    print(f"Wrote {EXCLUDED_STRUCTURES_CSV.relative_to(PROJECT_ROOT)} ({len(excluded_df)} rows)")
    assert len(structures_df) + len(excluded_df) == 2072

    print("\nRunning primary dataset validation (Step 15)...")
    errors = []
    if not (structures_df["experimental_method"].str.contains("X-RAY", case=False, na=False)).all():
        errors.append("non-X-ray row(s) found in structures.csv")
    if not (structures_df["positive_occupancy_lbd_coverage"] >= 0.90).all():
        errors.append("coverage < 0.90 row(s) found in structures.csv")
    if not (structures_df["resolution"] >= 1.8).all():
        errors.append("resolution < 1.8 row(s) found in structures.csv")
    if not (structures_df["resolution"] <= 3.5).all():
        errors.append("resolution > 3.5 row(s) found in structures.csv")
    if not (structures_df["r_free"] <= 0.30).all():
        errors.append("r_free > 0.30 row(s) found in structures.csv")
    if not (structures_df["identity_class"] == "HUMAN_OR_HUMAN_DERIVED").all():
        errors.append("non-human-derived identity_class row(s) found in structures.csv")
    if structures_df["instance_id"].isna().any() or (structures_df["instance_id"] == "").any():
        errors.append("missing representative instance_id")
    if structures_df["label_asym_id"].isna().any() or structures_df["auth_asym_id"].isna().any():
        errors.append("missing label/auth chain identifiers")
    if structures_df["mmcif_sha256"].isna().any() or (structures_df["mmcif_sha256"] == "").any():
        errors.append("missing mmCIF checksum")
    if not structures_df["ligand_state"].isin(["APO", "HOLO", "AMBIGUOUS"]).all():
        errors.append("ligand_state outside {APO,HOLO,AMBIGUOUS}")
    dup_key = structures_df.duplicated(subset=["uniprot_id", "pdb_id", "polymer_entity_id"])
    if dup_key.any():
        errors.append("duplicate (uniprot_id, pdb_id, polymer_entity_id) logical primary key")

    if errors:
        raise RuntimeError("STOPPING: primary dataset validation failed:\n  - " + "\n  - ".join(errors))
    print("  all primary dataset validation checks PASSED")

    print("\nBuilding final_counts_by_receptor.csv (48 rows)...")
    stage2_pool = pd.read_csv(MANIFESTS_DIR / "stage3b_coordinate_pool.csv", dtype={"pdb_id": str})
    stage2_counts = stage2_pool.groupby("uniprot_id").size().to_dict()
    stage3b_counts = stage2_pool[stage2_pool["coordinate_audit_selected"]].groupby("uniprot_id").size().to_dict()

    receptor_rows = []
    for rec in NR_METADATA:
        uid = rec["uniprot_id"]
        prim = structures_df[structures_df["uniprot_id"] == uid]
        receptor_rows.append(
            {
                "uniprot_id": uid, "nr_code": rec["nr_code"], "common_name": rec["common_name"],
                "stage2_candidates": stage2_counts.get(uid, 0),
                "stage3b_coordinate_candidates": stage3b_counts.get(uid, 0),
                "primary_selected_structures": len(prim),
                "primary_unique_pdb_entries": prim["pdb_id"].nunique(),
                "apo_count": int((prim["ligand_state"] == "APO").sum()),
                "holo_count": int((prim["ligand_state"] == "HOLO").sum()),
                "ambiguous_ligand_count": int((prim["ligand_state"] == "AMBIGUOUS").sum()),
                "95pct_sensitivity_count": int(prim["passes_95pct_completeness"].sum()),
                "99pct_sensitivity_count": int(prim["passes_99pct_completeness"].sum()),
                "fusion_count": int(prim["fusion_or_chimera"].sum()),
                "no_fusion_sensitivity_count": int(prim["passes_no_fusion_sensitivity"].sum()),
            }
        )
    counts_df = pd.DataFrame(receptor_rows)
    assert len(counts_df) == 48
    counts_df.to_csv(FINAL_COUNTS_BY_RECEPTOR_CSV, index=False)
    print(f"Wrote {FINAL_COUNTS_BY_RECEPTOR_CSV.relative_to(PROJECT_ROOT)} ({len(counts_df)} rows)")
    print(f"  receptors with primary_selected_structures > 0: {(counts_df['primary_selected_structures'] > 0).sum()}")
    print(f"  receptors with primary_selected_structures == 0: {(counts_df['primary_selected_structures'] == 0).sum()}")


if __name__ == "__main__":
    main()
