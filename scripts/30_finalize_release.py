"""
Final Data Release, Step 19: build the machine-readable release summary
data/manifests/final_data_infrastructure_release.json, including SHA256
checksums of every team-handoff deliverable and the parent Git commit.

Reads all already-built final manifests/processed tables — performs no
new computation, network access, or mmCIF parsing.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
REPORTS_TABLES_DIR = PROJECT_ROOT / "reports" / "tables"

RELEASE_JSON = MANIFESTS_DIR / "final_data_infrastructure_release.json"

CHECKSUM_FILES = {
    "structures_csv": PROCESSED_DIR / "structures.csv",
    "lbd_residue_bfactors_all_parquet": PROCESSED_DIR / "lbd_residue_bfactors_all.parquet",
    "structures_primary_qc_subset_csv": PROCESSED_DIR / "structures_primary_qc_subset.csv",
    "final_lbd_residue_map_primary_qc_subset_parquet": PROCESSED_DIR / "final_lbd_residue_map_primary_qc_subset.parquet",
    "excluded_structures_csv": PROCESSED_DIR / "excluded_structures.csv",
    "final_ligand_annotations_csv": MANIFESTS_DIR / "final_ligand_annotations.csv",
    "final_polymer_partner_annotations_csv": MANIFESTS_DIR / "final_polymer_partner_annotations.csv",
    "team_handoff_md": PROJECT_ROOT / "TEAM_HANDOFF.md",
    "final_data_dictionary_md": REPORTS_TABLES_DIR / "final_data_dictionary.md",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_head_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:
        return ""


def main() -> None:
    structures = pd.read_csv(PROCESSED_DIR / "structures_primary_qc_subset.csv", dtype={"pdb_id": str})
    master = pd.read_csv(PROCESSED_DIR / "structures.csv", dtype={"pdb_id": str})
    excluded = pd.read_csv(PROCESSED_DIR / "excluded_structures.csv", dtype={"pdb_id": str})
    rep_selection = pd.read_csv(MANIFESTS_DIR / "final_representative_instance_selection.csv", dtype={"pdb_id": str})
    residue_map = pd.read_parquet(PROCESSED_DIR / "final_lbd_residue_map_primary_qc_subset.parquet")
    residue_all = pd.read_parquet(PROCESSED_DIR / "lbd_residue_bfactors_all.parquet")
    ligand_all = pd.read_csv(MANIFESTS_DIR / "final_ligand_annotations.csv", dtype={"pdb_id": str})
    strict_keys = set(zip(structures["pdb_id"], structures["uniprot_id"], structures["polymer_entity_id"]))
    ligand = ligand_all[[k in strict_keys for k in zip(ligand_all["pdb_id"], ligand_all["uniprot_id"], ligand_all["polymer_entity_id"])]]
    download_manifest = pd.read_csv(MANIFESTS_DIR / "mmcif_download_manifest.csv", dtype={"pdb_id": str})
    counts_by_receptor = pd.read_csv(MANIFESTS_DIR / "final_counts_by_receptor.csv")

    import json as _json
    exclusion_reason_counts: dict[str, int] = {}
    for reasons_json in excluded["primary_exclusion_reasons_json"].dropna():
        for reason in _json.loads(reasons_json):
            exclusion_reason_counts[reason] = exclusion_reason_counts.get(reason, 0) + 1

    checksums = {}
    for key, path in CHECKSUM_FILES.items():
        if path.exists():
            checksums[key] = {"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        else:
            checksums[key] = {"path": str(path.relative_to(PROJECT_ROOT)), "sha256": None, "size_bytes": None}

    release = {
        "project": "Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families",
        "release_name": "final-data-infrastructure-v2-master-inventory",
        "validated_receptors": 48,
        "dataset_model": "structures.csv is the MASTER inventory (QC criteria are metadata/analysis flags, never deletion criteria); primary_qc_subset files are a historical strict analysis view",
        "nested_views": {
            "stage2_discovery_candidates_master_structures_csv": len(master),
            "xray_lbd_overlap_coordinate_audit_candidates": int(master["has_stage3b_coordinate_data"].sum()),
            "previous_strict_primary_qc_subset": int(master["previous_primary_qc_member"].sum()),
        },
        "master_receptors_represented": int(master["uniprot_id"].nunique()),
        "lbd_residue_bfactors_all_rows": len(residue_all),
        "lbd_residue_bfactors_all_instance_models": int(residue_all["instance_id"].nunique()),
        "primary_structure_count": len(structures),
        "primary_pdb_count": int(structures["pdb_id"].nunique()),
        "receptors_represented_in_primary": int(structures["uniprot_id"].nunique()),
        "receptors_with_zero_primary_structures": int((counts_by_receptor["primary_selected_structures"] == 0).sum()),
        "ligand_state_counts": {
            "APO": int((ligand["ligand_state"] == "APO").sum()),
            "HOLO": int((ligand["ligand_state"] == "HOLO").sum()),
            "AMBIGUOUS": int((ligand["ligand_state"] == "AMBIGUOUS").sum()),
        },
        "ligand_state_counts_all_coordinate_candidates": {
            "APO": int((ligand_all["ligand_state"] == "APO").sum()),
            "HOLO": int((ligand_all["ligand_state"] == "HOLO").sum()),
            "AMBIGUOUS": int((ligand_all["ligand_state"] == "AMBIGUOUS").sum()),
            "NOT_ASSESSABLE_NO_COORDINATE_DATA": int((~master["has_stage3b_coordinate_data"]).sum()),
        },
        "ligand_annotation_review_needed_count_all_coordinate_candidates": int((ligand_all["ligand_annotation_status"] == "REVIEW_NEEDED").sum()),
        "excluded_candidate_count": len(excluded),
        "exclusion_reason_counts": exclusion_reason_counts,
        "sensitivity_set_counts": {
            "passes_95pct_completeness": int(structures["passes_95pct_completeness"].sum()),
            "passes_99pct_completeness": int(structures["passes_99pct_completeness"].sum()),
            "passes_no_fusion_sensitivity": int(structures["passes_no_fusion_sensitivity"].sum()),
            "ultrahigh_resolution_preserved_records": int(
                (excluded["primary_exclusion_reasons_json"] == '["RESOLUTION_LT_1_8"]').sum()
            ),
        },
        "representative_instances_selected": len(rep_selection),
        "residue_map_rows": len(residue_map),
        "residue_map_usable_ca_rows": int(residue_map["selected_ca_occupancy"].notna().sum()),
        "mmcif_files_locally_available": int(download_manifest["download_status"].isin(["SUCCESS", "CACHED_VALID"]).sum()),
        "ligand_annotation_complete_count": int((ligand["ligand_annotation_status"] == "RESOLVED").sum()),
        "ligand_annotation_review_needed_count": int((ligand["ligand_annotation_status"] == "REVIEW_NEEDED").sum()),
        "final_selection_policy": {
            "version": "final-structure-selection-v1",
            "method": "X-RAY DIFFRACTION only",
            "completeness_threshold": "positive-occupancy full standardized-LBD Ca coverage >= 0.90",
            "resolution_range_angstrom": [1.8, 3.5],
            "r_free_max": 0.30,
            "identity": "HUMAN_OR_HUMAN_DERIVED only (individually audited taxonomy)",
            "representative_instance_granularity": "project receptor x PDB entry x receptor polymer entity",
        },
        "checksums": checksums,
        "parent_git_commit": git_head_commit(),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }

    RELEASE_JSON.write_text(json.dumps(release, indent=2, default=str) + "\n")
    print(f"Wrote {RELEASE_JSON.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
