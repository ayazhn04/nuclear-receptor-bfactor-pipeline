"""
Stage 3B, step 3: aggregate the instance/model-level coordinate-QC table
into candidate-level and receptor-level summaries, the machine-readable
Stage 3B summary, the human-readable report, and provenance.

Reads:
    data/manifests/stage3b_coordinate_pool.csv
    data/manifests/mmcif_download_manifest.csv
    data/manifests/rcsb_lbd_coordinate_instance_qc.csv
    data/manifests/mmcif_unobserved_residue_crosscheck.csv
    data/manifests/rcsb_precoordinate_qc.csv
    data/manifests/rcsb_source_taxonomy_audit.csv
    data/manifests/rcsb_multi_uniprot_entity_audit.csv

Writes:
    data/manifests/rcsb_lbd_coordinate_candidate_summary.csv
    data/manifests/stage3b_counts_by_receptor.csv
    data/manifests/stage3b_coordinate_summary.json
    reports/tables/stage3b_coordinate_completeness_report.md
    data/manifests/stage3b_provenance.json

No final keep/exclude decision, no representative-chain selection, no
completeness threshold is chosen here.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
REPORTS_TABLES_DIR = PROJECT_ROOT / "reports" / "tables"

POOL_CSV = MANIFESTS_DIR / "stage3b_coordinate_pool.csv"
DOWNLOAD_MANIFEST_CSV = MANIFESTS_DIR / "mmcif_download_manifest.csv"
INSTANCE_QC_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_instance_qc.csv"
CROSSCHECK_CSV = MANIFESTS_DIR / "mmcif_unobserved_residue_crosscheck.csv"
PRECOORDINATE_QC_CSV = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"
TAXONOMY_AUDIT_CSV = MANIFESTS_DIR / "rcsb_source_taxonomy_audit.csv"
MULTI_UNIPROT_CSV = MANIFESTS_DIR / "rcsb_multi_uniprot_entity_audit.csv"
STAGE3A_RELEASE_SUMMARY = MANIFESTS_DIR / "stage3a_release_summary.json"
LBD_ANALYSIS_POLICY = MANIFESTS_DIR / "nr_lbd_analysis_policy.json"

CANDIDATE_SUMMARY_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_candidate_summary.csv"
COUNTS_BY_RECEPTOR_CSV = MANIFESTS_DIR / "stage3b_counts_by_receptor.csv"
SUMMARY_JSON = MANIFESTS_DIR / "stage3b_coordinate_summary.json"
REPORT_MD = REPORTS_TABLES_DIR / "stage3b_coordinate_completeness_report.md"
PROVENANCE_JSON = MANIFESTS_DIR / "stage3b_provenance.json"


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def coverage_bin(value) -> str:
    if pd.isna(value) or value == "":
        return "MISSING"
    v = float(value)
    if v == 1.0:
        return "1.000"
    if v >= 0.99:
        return "[0.99,1.00)"
    if v >= 0.95:
        return "[0.95,0.99)"
    if v >= 0.90:
        return "[0.90,0.95)"
    if v >= 0.80:
        return "[0.80,0.90)"
    return "<0.80"


def count_bin(value) -> str:
    if pd.isna(value) or value == "":
        return "MISSING"
    v = int(value)
    if v == 0:
        return "0"
    if v <= 2:
        return "1-2"
    if v <= 5:
        return "3-5"
    if v <= 10:
        return "6-10"
    return ">10"


def build_candidate_summary(instance_qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    measured = instance_qc[instance_qc["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]
    grouped = instance_qc.groupby(["uniprot_id", "nr_code", "common_name", "group", "pdb_id", "polymer_entity_id", "entity_id"], dropna=False)

    for key, g in grouped:
        uid, nr_code, common_name, grp, pdb_id, entity_id_full, entity_id = key
        g_measured = g[g["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]

        instance_count = g["instance_id"].nunique()
        single_model = g.groupby("instance_id")["model_num"].nunique().eq(1).sum()
        multi_model = g.groupby("instance_id")["model_num"].nunique().gt(1).sum()

        full_cov = pd.to_numeric(g_measured["observed_full_lbd_ca_coverage"], errors="coerce")
        within_cov = pd.to_numeric(g_measured["observed_within_construct_ca_coverage"], errors="coerce")
        bfactor_full_cov = pd.to_numeric(g_measured["bfactor_usable_full_lbd_coverage"], errors="coerce")
        lbd_length = pd.to_numeric(g_measured["lbd_length"], errors="coerce")
        mapped = pd.to_numeric(g_measured["mapped_lbd_residues"], errors="coerce")
        modeled = pd.to_numeric(g_measured["modeled_lbd_ca_residues"], errors="coerce")

        construct_full = (mapped == lbd_length) if len(mapped) else pd.Series([], dtype=bool)
        ca_full_complete = (modeled == lbd_length) if len(modeled) else pd.Series([], dtype=bool)

        review_needed = (g["coordinate_qc_status"] != "MEASURED").any()

        rows.append(
            {
                "uniprot_id": uid, "nr_code": nr_code, "common_name": common_name, "group": grp,
                "pdb_id": pdb_id, "polymer_entity_id": entity_id_full, "entity_id": entity_id,
                "instance_count": instance_count,
                "single_model_instance_count": int(single_model),
                "multi_model_instance_count": int(multi_model),
                "min_full_lbd_ca_coverage": full_cov.min() if len(full_cov) else "",
                "max_full_lbd_ca_coverage": full_cov.max() if len(full_cov) else "",
                "median_full_lbd_ca_coverage": full_cov.median() if len(full_cov) else "",
                "min_within_construct_ca_coverage": within_cov.min() if within_cov.notna().any() else "",
                "max_within_construct_ca_coverage": within_cov.max() if within_cov.notna().any() else "",
                "max_bfactor_usable_full_lbd_coverage": bfactor_full_cov.max() if len(bfactor_full_cov) else "",
                "any_instance_full_lbd_ca_complete": bool(ca_full_complete.any()) if len(ca_full_complete) else False,
                "all_instances_full_lbd_ca_complete": bool(ca_full_complete.all()) if len(ca_full_complete) else False,
                "construct_maps_entire_standardized_lbd": bool(construct_full.any()) if len(construct_full) else False,
                "coordinate_review_needed": bool(review_needed),
            }
        )
    return pd.DataFrame(rows)


def build_counts_by_receptor(pool: pd.DataFrame, candidate_summary: pd.DataFrame, instance_qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in NR_METADATA:
        uid = r["uniprot_id"]
        stage2_candidates = int((pool["uniprot_id"] == uid).sum())
        sel = pool[(pool["uniprot_id"] == uid) & (pool["coordinate_audit_selected"])]
        stage3b_selected = len(sel)
        selected_pdbs = sel["pdb_id"].nunique()

        cs = candidate_summary[candidate_summary["uniprot_id"] == uid]
        selected_instances = int(cs["instance_count"].sum()) if len(cs) else 0
        iq = instance_qc[instance_qc["uniprot_id"] == uid]
        selected_instance_models = len(iq[iq["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])])

        construct_full_count = int(cs["construct_maps_entire_standardized_lbd"].sum()) if len(cs) else 0
        full_ca_complete_count = int(cs["any_instance_full_lbd_ca_complete"].sum()) if len(cs) else 0

        full_cov = pd.to_numeric(iq[iq["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]["observed_full_lbd_ca_coverage"], errors="coerce")
        bfactor_cov = pd.to_numeric(iq[iq["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]["bfactor_usable_full_lbd_coverage"], errors="coerce")
        lbd_len = pd.to_numeric(iq["lbd_length"], errors="coerce")
        modeled_ca = pd.to_numeric(iq["modeled_lbd_ca_residues"], errors="coerce")
        bfactor_complete_count = int((pd.to_numeric(iq["bfactor_usable_lbd_ca_residues"], errors="coerce") == lbd_len).sum())

        review_count = int((iq["coordinate_qc_status"] != "MEASURED").sum())

        rows.append(
            {
                "common_name": r["common_name"], "nr_code": r["nr_code"], "uniprot_id": uid,
                "stage2_candidates": stage2_candidates,
                "stage3b_selected_candidates": stage3b_selected,
                "selected_unique_pdb_entries": selected_pdbs,
                "selected_instances": selected_instances,
                "selected_instance_models": selected_instance_models,
                "construct_full_lbd_count": construct_full_count,
                "full_lbd_ca_complete_count": full_ca_complete_count,
                "full_lbd_bfactor_complete_count": bfactor_complete_count,
                "coverage_ge_0_99_count": int((full_cov >= 0.99).sum()),
                "coverage_ge_0_95_count": int((full_cov >= 0.95).sum()),
                "coverage_ge_0_90_count": int((full_cov >= 0.90).sum()),
                "coordinate_review_needed_count": review_count,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    pool = pd.read_csv(POOL_CSV, dtype={"pdb_id": str})
    downloads = pd.read_csv(DOWNLOAD_MANIFEST_CSV, dtype={"pdb_id": str})
    instance_qc = pd.read_csv(INSTANCE_QC_CSV, dtype={"pdb_id": str})
    crosscheck = pd.read_csv(CROSSCHECK_CSV, dtype={"pdb_id": str}) if CROSSCHECK_CSV.exists() else pd.DataFrame()
    precoordinate_qc = pd.read_csv(PRECOORDINATE_QC_CSV, dtype={"pdb_id": str})
    taxonomy_audit = pd.read_csv(TAXONOMY_AUDIT_CSV, dtype={"pdb_id": str})
    multi_uniprot = pd.read_csv(MULTI_UNIPROT_CSV) if MULTI_UNIPROT_CSV.exists() else pd.DataFrame(columns=["polymer_entity_id"])

    print("Building candidate-level coordinate summary...")
    candidate_summary = build_candidate_summary(instance_qc)
    candidate_summary.to_csv(CANDIDATE_SUMMARY_CSV, index=False)
    print(f"Wrote {CANDIDATE_SUMMARY_CSV.relative_to(PROJECT_ROOT)} ({len(candidate_summary)} rows)")

    print("\nBuilding per-receptor counts...")
    counts_by_receptor = build_counts_by_receptor(pool, candidate_summary, instance_qc)
    counts_by_receptor.to_csv(COUNTS_BY_RECEPTOR_CSV, index=False)
    print(f"Wrote {COUNTS_BY_RECEPTOR_CSV.relative_to(PROJECT_ROOT)} ({len(counts_by_receptor)} rows)")

    # --- Stage 3B summary JSON ---
    selected = pool[pool["coordinate_audit_selected"]]
    measured = instance_qc[instance_qc["coordinate_qc_status"].isin(["MEASURED", "MULTI_MODEL_REVIEW"])]
    full_cov = pd.to_numeric(measured["observed_full_lbd_ca_coverage"], errors="coerce")

    dl_status_counts = downloads["download_status"].value_counts().to_dict()
    multi_model_entries = downloads.merge(
        instance_qc[instance_qc["coordinate_qc_status"] == "MULTI_MODEL_REVIEW"][["pdb_id"]].drop_duplicates(),
        on="pdb_id",
    )

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection": {
            "total_stage2_candidates": len(pool),
            "selected_coordinate_audit_candidates": len(selected),
            "selected_unique_entities": int(selected["polymer_entity_id"].nunique()),
            "selected_unique_pdb_entries": int(selected["pdb_id"].nunique()),
        },
        "downloads": {
            "requested": int(selected["pdb_id"].nunique()),
            "successful": int(dl_status_counts.get("SUCCESS", 0)),
            "cached_valid": int(dl_status_counts.get("CACHED_VALID", 0)),
            "failed": int(dl_status_counts.get("FAILED", 0)),
            "failed_pdb_ids": downloads.loc[downloads["download_status"] == "FAILED", "pdb_id"].tolist(),
        },
        "instances_models": {
            "instance_model_rows": len(instance_qc),
            "measured_rows": len(measured),
            "multi_model_entries": int(multi_model_entries["pdb_id"].nunique()),
        },
        "coordinate_completeness": {
            "exact_full_lbd_ca_complete": int((pd.to_numeric(measured["modeled_lbd_ca_residues"], errors="coerce") == pd.to_numeric(measured["lbd_length"], errors="coerce")).sum()),
            "exact_mapped_construct_ca_complete": int((pd.to_numeric(measured["modeled_lbd_ca_residues"], errors="coerce") == pd.to_numeric(measured["mapped_lbd_residues"], errors="coerce")).sum()),
            "exact_full_lbd_bfactor_complete": int((pd.to_numeric(measured["bfactor_usable_lbd_ca_residues"], errors="coerce") == pd.to_numeric(measured["lbd_length"], errors="coerce")).sum()),
            "coverage_bins": full_cov.apply(coverage_bin).value_counts().to_dict(),
        },
        "missingness": {
            "mapped_but_unmodeled_bins": pd.to_numeric(measured["mapped_but_unmodeled_lbd_residues"], errors="coerce").apply(count_bin).value_counts().to_dict(),
            "unmapped_lbd_bins": pd.to_numeric(measured["unmapped_lbd_residues"], errors="coerce").apply(count_bin).value_counts().to_dict(),
        },
        "bfactor_usability": {
            "modeled_but_bfactor_unusable_total": int(pd.to_numeric(measured["modeled_but_bfactor_unusable_residues"], errors="coerce").sum()),
            "zero_occupancy_ca_residue_total": int(pd.to_numeric(measured["zero_occupancy_ca_residue_count"], errors="coerce").sum()),
            "alternate_ca_residue_total": int(pd.to_numeric(measured["alternate_ca_residue_count"], errors="coerce").sum()),
        },
        "anisotropy": {
            "files_with_anisotrop_category": int(measured[measured["anisotrop_category_present"] == True]["pdb_id"].nunique()),  # noqa: E712
        },
        "review": {
            "mapping_review_needed": int((instance_qc["coordinate_qc_status"] == "MAPPING_REVIEW_NEEDED").sum()),
            "parse_review_needed": int((instance_qc["coordinate_qc_status"] == "PARSE_REVIEW_NEEDED").sum()),
            "multi_model_review": int((instance_qc["coordinate_qc_status"] == "MULTI_MODEL_REVIEW").sum()),
        },
        # Total distinct mmCIF files successfully acquired and validated
        # (fresh downloads this run + previously-cached files re-validated
        # this run) — NOT just this run's fresh-download count, which would
        # misleadingly read as 0 on a pure cache-validation rerun even
        # though 1777 valid files are present on disk.
        "coordinate_files_downloaded": int(dl_status_counts.get("SUCCESS", 0)) + int(dl_status_counts.get("CACHED_VALID", 0)),
        "final_structure_qc_applied": False,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"\nWrote {SUMMARY_JSON.relative_to(PROJECT_ROOT)}")
    print(json.dumps(summary["coordinate_completeness"]["coverage_bins"], indent=2))

    write_report(pool, downloads, instance_qc, candidate_summary, counts_by_receptor, summary, precoordinate_qc, taxonomy_audit, multi_uniprot)
    print(f"Wrote {REPORT_MD.relative_to(PROJECT_ROOT)}")

    write_provenance(summary)
    print(f"Wrote {PROVENANCE_JSON.relative_to(PROJECT_ROOT)}")


def write_report(pool, downloads, instance_qc, candidate_summary, counts_by_receptor, summary, precoordinate_qc, taxonomy_audit, multi_uniprot) -> None:
    lines = ["# Stage 3B — mmCIF Coordinate Acquisition and Observed LBD Cα Completeness Audit\n"]
    lines.append(
        "## 1. Scope and selection rule\n\n"
        f"Coordinate-audit pool: all {summary['selection']['total_stage2_candidates']} Stage 2 candidates. "
        f"Selected for mmCIF download/analysis: `method_is_xray AND mapping_overlaps_lbd` (standardized "
        f"Stage 3A LBD interval) — **{summary['selection']['selected_coordinate_audit_candidates']}** "
        f"candidates, {summary['selection']['selected_unique_entities']} unique entities, "
        f"{summary['selection']['selected_unique_pdb_entries']} unique PDB entries. No resolution, R-free, "
        f"taxonomy, fusion, or metadata-prefilter restriction was applied to this selection.\n"
    )
    lines.append(
        "## 2. Download provenance\n\n"
        f"Requested: {summary['downloads']['requested']}, successful: {summary['downloads']['successful']}, "
        f"cached_valid: {summary['downloads']['cached_valid']}, failed: {summary['downloads']['failed']}.\n"
    )
    if summary["downloads"]["failed_pdb_ids"]:
        lines.append(f"Failed PDB IDs: {summary['downloads']['failed_pdb_ids']}\n")
    lines.append(
        "## 3. Why label_seq_id, not author numbering\n\n"
        "Canonical residue mapping is built exclusively from mmCIF `_atom_site.label_seq_id` "
        "(polymer entity sequence position), matched against the official RCSB Sequence Coordinates "
        "alignment's `target_begin`/`target_end` (also entity-sequence-numbered). Author residue "
        "numbers (`auth_seq_id`) are carried through the observation table only as display/provenance "
        "metadata — they are never used to compute canonical UniProt positions, since author numbering "
        "can be discontinuous, offset, or reused across chains/constructs.\n"
    )
    lines.append(
        "## 4. Coverage metric definitions\n\n"
        "- **Construct mapping coverage** (Stage 3A): mapped_lbd_residues / lbd_length — sequence-level "
        "construct coverage only.\n"
        "- **Observed full-LBD Cα coverage**: modeled_lbd_ca_residues / lbd_length.\n"
        "- **Observed within-construct Cα coverage**: modeled_lbd_ca_residues / mapped_lbd_residues.\n"
        "- **B-factor-usable** variants of the above use residues with a numeric, finite `B_iso_or_equiv` "
        "and positive occupancy.\n"
    )
    lines.append("## 5. Download success/failures\n\n" + f"See Section 2.\n")
    lines.append("## 6. Coverage distributions (observed full-LBD Cα coverage, measured instance-models)\n\n")
    for k, v in summary["coordinate_completeness"]["coverage_bins"].items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append(
        f"- Exact full-LBD Cα complete: {summary['coordinate_completeness']['exact_full_lbd_ca_complete']}\n"
        f"- Exact mapped-construct Cα complete: {summary['coordinate_completeness']['exact_mapped_construct_ca_complete']}\n"
        f"- Exact full-LBD B-factor complete: {summary['coordinate_completeness']['exact_full_lbd_bfactor_complete']}\n"
    )
    lines.append("## 7. Missingness topology\n\n")
    lines.append(f"mapped-but-unmodeled bins: {summary['missingness']['mapped_but_unmodeled_bins']}")
    lines.append(f"unmapped-LBD bins: {summary['missingness']['unmapped_lbd_bins']}\n")
    lines.append("## 8. B-factor availability\n\n")
    lines.append(
        f"- modeled-but-B-factor-unusable (total residues): {summary['bfactor_usability']['modeled_but_bfactor_unusable_total']}\n"
        f"- zero-occupancy Cα residues (total): {summary['bfactor_usability']['zero_occupancy_ca_residue_total']}\n"
        f"- alternate-conformer Cα residues (total): {summary['bfactor_usability']['alternate_ca_residue_total']}\n"
    )
    lines.append("## 9. Anisotropic-refinement results\n\n")
    lines.append(f"Files with an anisotropic atom category present: {summary['anisotropy']['files_with_anisotrop_category']}\n")
    lines.append("## 10. Per-receptor summary\n\n")
    lines.append("See `data/manifests/stage3b_counts_by_receptor.csv` (all 48 receptors).\n")
    lines.append("## 11. Important review pools\n\n")
    lines.append(
        f"- Mapping review needed: {summary['review']['mapping_review_needed']}\n"
        f"- Parse review needed: {summary['review']['parse_review_needed']}\n"
        f"- Multi-model review: {summary['review']['multi_model_review']}\n"
    )
    lines.append("## 12. Important statement\n\nNo final LBD-coordinate completeness threshold was selected at Stage 3B.\n")
    lines.append("## 13. Important statement\n\nNo representative chain was selected at Stage 3B.\n")
    lines.append("## 14. Important statement\n\nNo B-factor normalization or downstream statistical analysis was performed at Stage 3B.\n")

    REPORT_MD.write_text("\n".join(lines))


def write_provenance(summary: dict) -> None:
    stage3a_release = json.loads(STAGE3A_RELEASE_SUMMARY.read_text()) if STAGE3A_RELEASE_SUMMARY.exists() else {}
    provenance = {
        "project": "Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families",
        "stage_description": "Stage 3B: selective mmCIF acquisition and observed LBD Cα coordinate completeness audit.",
        "standardized_lbd_policy_sha256": sha256_of_file(LBD_ANALYSIS_POLICY),
        "stage3a_release_summary_sha256": sha256_of_file(STAGE3A_RELEASE_SUMMARY),
        "coordinate_pool_sha256": sha256_of_file(POOL_CSV),
        "mmcif_download_manifest_sha256": sha256_of_file(DOWNLOAD_MANIFEST_CSV),
        "coordinate_qc_table_sha256": sha256_of_file(INSTANCE_QC_CSV),
        "candidate_summary_sha256": sha256_of_file(CANDIDATE_SUMMARY_CSV),
        "per_receptor_summary_sha256": sha256_of_file(COUNTS_BY_RECEPTOR_CSV),
        "report_sha256": sha256_of_file(REPORT_MD),
        "official_rcsb_file_endpoint": "https://files.rcsb.org/download/<PDB_ID>.cif.gz",
        "number_of_files": summary["downloads"]["successful"] + summary["downloads"]["cached_valid"],
        "python_version": sys.version,
        "generating_scripts": [
            "scripts/18_download_mmcif.py",
            "scripts/19_process_lbd_coordinates.py",
            "scripts/20_build_stage3b_summaries.py",
        ],
        "coordinate_convention": "1-based inclusive UniProt canonical positions; mmCIF label_seq_id for entity mapping",
        "atom_presence_definition": "ca_modeled = at least one _atom_site record with label_atom_id=='CA' for the mapped label_seq_id (any altloc)",
        "bfactor_usable_definition": "ca_bfactor_usable = at least one CA record with numeric finite B_iso_or_equiv AND occupancy > 0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    PROVENANCE_JSON.write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
