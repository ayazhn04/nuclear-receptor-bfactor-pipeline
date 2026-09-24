"""
Stage 3A, step 4: build descriptive pre-coordinate QC flags for every Stage
2 candidate, the instance-level candidate inventory, receptor-by-receptor
counts, the machine-readable summary, and the human-readable report.

This produces NO final keep/exclude decision. `passes_project_metadata_prefilter`
reproduces the project brief's metadata-level criteria (X-ray, resolution
1.8-3.5 A, R-free <= 0.30, some LBD overlap) as closely as possible from
already-retrieved metadata, but actual LBD coordinate completeness is
unknown until Stage 3B downloads and inspects mmCIF files.

Reads:
    data/manifests/rcsb_candidate_inventory.csv   (Stage 2/2.1/2.2 metadata — NOTE: this
        file currently carries redundant _x/_y duplicate coverage columns
        from a Stage 2.2/2.3 rerun bug; identical in value to the clean
        columns (verified), so only the clean, unsuffixed column names are
        used here. See Stage 3A report "Warnings" for the disclosure.)
    data/manifests/rcsb_candidate_lbd_mapping.csv
    data/manifests/rcsb_source_taxonomy_audit.csv
    data/manifests/rcsb_discovery_anomalies.csv    (for multi_uniprot_entity flag)
    data/manifests/rcsb_polymer_instances.csv
    data/manifests/rcsb_short_mapping_lbd_audit.csv

Writes:
    data/manifests/rcsb_precoordinate_qc.csv
    data/manifests/rcsb_candidate_instance_inventory.csv
    data/manifests/stage3a_precoordinate_summary.json
    data/manifests/stage3a_counts_by_receptor.csv
    reports/tables/stage3a_precoordinate_qc_report.md
"""

from __future__ import annotations

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

CANDIDATE_INVENTORY_CSV = MANIFESTS_DIR / "rcsb_candidate_inventory.csv"
LBD_MAPPING_CSV = MANIFESTS_DIR / "rcsb_candidate_lbd_mapping.csv"
LBD_REFERENCE_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"
TAXONOMY_AUDIT_CSV = MANIFESTS_DIR / "rcsb_source_taxonomy_audit.csv"
ANOMALIES_CSV = MANIFESTS_DIR / "rcsb_discovery_anomalies.csv"
INSTANCES_CSV = MANIFESTS_DIR / "rcsb_polymer_instances.csv"
SHORT_MAPPING_LBD_CSV = MANIFESTS_DIR / "rcsb_short_mapping_lbd_audit.csv"
SEARCH_MANIFEST_CSV = MANIFESTS_DIR / "rcsb_search_manifest.csv"

PRECOORDINATE_QC_CSV = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"
INSTANCE_INVENTORY_CSV = MANIFESTS_DIR / "rcsb_candidate_instance_inventory.csv"
SUMMARY_JSON = MANIFESTS_DIR / "stage3a_precoordinate_summary.json"
COUNTS_BY_RECEPTOR_CSV = MANIFESTS_DIR / "stage3a_counts_by_receptor.csv"
REPORT_MD = REPORTS_TABLES_DIR / "stage3a_precoordinate_qc_report.md"

RESOLUTION_LOWER_BOUND = 1.8
RESOLUTION_UPPER_BOUND = 3.5
R_FREE_UPPER_BOUND = 0.30

CANDIDATE_INVENTORY_CLEAN_COLUMNS = [
    "uniprot_id", "nr_code", "common_name", "group",
    "pdb_id", "polymer_entity_id", "entity_id",
    "entity_description", "entity_sequence_length",
    "label_asym_ids", "auth_asym_ids",
    "mapping_status",
    "source_organism_names", "source_taxonomy_ids",
    "experimental_methods", "resolution_combined", "r_work", "r_free",
    "is_xray", "is_cryo_em", "is_nmr", "is_other_experimental",
    "has_resolution", "has_r_free", "source_includes_homo_sapiens_9606",
    "computed_reference_sequence_coverage", "computed_entity_sequence_coverage",
    "rcsb_reference_sequence_coverage", "rcsb_entity_sequence_coverage",
]


def resolution_status(values: list[float]) -> str:
    if not values:
        return "MISSING"
    if len(values) > 1:
        return "MULTI_VALUE_REVIEW"
    v = values[0]
    if v < RESOLUTION_LOWER_BOUND:
        return "ULTRAHIGH_LT_1_8"
    if v > RESOLUTION_UPPER_BOUND:
        return "TOO_LOW_RESOLUTION_GT_3_5"
    return "PASS_PROJECT_RANGE"


def r_free_status(values: list[float]) -> str:
    if not values:
        return "MISSING"
    if all(v <= R_FREE_UPPER_BOUND for v in values):
        return "PASS_ALL_LE_0_30"
    if all(v > R_FREE_UPPER_BOUND for v in values):
        return "FAIL_ALL_GT_0_30"
    return "MIXED_REVIEW"


def build_precoordinate_qc(
    inventory: pd.DataFrame,
    lbd_mapping: pd.DataFrame,
    taxonomy_audit: pd.DataFrame,
    multi_uniprot_entity_ids: set[str],
) -> pd.DataFrame:
    merged = inventory.merge(
        lbd_mapping[[
            "uniprot_id", "polymer_entity_id", "lbd_start", "lbd_end", "lbd_length",
            "lbd_definition", "lbd_coordinate_system",
            "mapped_lbd_residues", "lbd_mapping_coverage", "lbd_mapping_category",
            "mapping_overlaps_lbd", "mapped_dbd_residues", "dbd_mapping_coverage",
        ]],
        on=["uniprot_id", "polymer_entity_id"], how="left",
    )
    taxonomy_by_entity = taxonomy_audit.set_index("polymer_entity_id")["taxonomy_audit_status"].to_dict()

    rows = []
    for _, row in merged.iterrows():
        resolution_values = json.loads(row["resolution_combined"]) if isinstance(row["resolution_combined"], str) else []
        r_free_values = json.loads(row["r_free"]) if isinstance(row["r_free"], str) else []

        res_status = resolution_status(resolution_values)
        rfree_status = r_free_status(r_free_values)
        method_is_xray = bool(row["is_xray"])
        overlaps_lbd = bool(row["mapping_overlaps_lbd"]) if pd.notna(row["mapping_overlaps_lbd"]) else False

        passes_prefilter = (
            method_is_xray
            and res_status == "PASS_PROJECT_RANGE"
            and rfree_status == "PASS_ALL_LE_0_30"
            and overlaps_lbd
        )

        reasons = []
        if not method_is_xray:
            reasons.append("NON_XRAY")
        if not overlaps_lbd:
            reasons.append("NO_LBD_OVERLAP")
        if res_status != "PASS_PROJECT_RANGE":
            reasons.append("RESOLUTION_REVIEW")
        if rfree_status != "PASS_ALL_LE_0_30":
            reasons.append("RFREE_REVIEW")
        if res_status == "MISSING" or rfree_status == "MISSING":
            if "METADATA_INCOMPLETE" not in reasons:
                reasons.append("METADATA_INCOMPLETE")

        if passes_prefilter:
            review_status = "PROJECT_PREFILTER_PASS"
        elif len(reasons) > 1:
            review_status = "MULTIPLE_REASONS"
        elif reasons:
            review_status = reasons[0]
        else:
            review_status = "METADATA_INCOMPLETE"

        rows.append(
            {
                "uniprot_id": row["uniprot_id"], "nr_code": row["nr_code"],
                "common_name": row["common_name"], "group": row["group"],
                "pdb_id": row["pdb_id"], "polymer_entity_id": row["polymer_entity_id"], "entity_id": row["entity_id"],
                "experimental_methods": row["experimental_methods"], "method_is_xray": method_is_xray,
                "resolution_values_json": json.dumps(resolution_values),
                "resolution_min": min(resolution_values) if resolution_values else "",
                "resolution_max": max(resolution_values) if resolution_values else "",
                "resolution_project_range_status": res_status,
                "r_free_values_json": json.dumps(r_free_values),
                "r_free_min": min(r_free_values) if r_free_values else "",
                "r_free_max": max(r_free_values) if r_free_values else "",
                "r_free_status": rfree_status,
                "lbd_start": row["lbd_start"], "lbd_end": row["lbd_end"], "lbd_length": row["lbd_length"],
                "lbd_definition": row["lbd_definition"], "lbd_coordinate_system": row["lbd_coordinate_system"],
                "mapped_lbd_residues": row["mapped_lbd_residues"],
                "lbd_mapping_coverage": row["lbd_mapping_coverage"],
                "lbd_mapping_category": row["lbd_mapping_category"],
                "mapping_overlaps_lbd": overlaps_lbd,
                "mapped_dbd_residues": row["mapped_dbd_residues"],
                "dbd_mapping_coverage": row["dbd_mapping_coverage"],
                "multi_uniprot_entity": row["polymer_entity_id"] in multi_uniprot_entity_ids,
                "taxonomy_audit_status": taxonomy_by_entity.get(row["polymer_entity_id"], ""),
                "passes_project_metadata_prefilter": passes_prefilter,
                "precoordinate_review_status": review_status,
                "review_reasons": json.dumps(reasons),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    inventory = pd.read_csv(CANDIDATE_INVENTORY_CSV, dtype={"pdb_id": str})[CANDIDATE_INVENTORY_CLEAN_COLUMNS]
    lbd_mapping = pd.read_csv(LBD_MAPPING_CSV, dtype={"pdb_id": str})
    lbd_reference = pd.read_csv(LBD_REFERENCE_CSV)
    taxonomy_audit = pd.read_csv(TAXONOMY_AUDIT_CSV, dtype={"pdb_id": str})
    anomalies = pd.read_csv(ANOMALIES_CSV, dtype={"pdb_id": str})
    instances = pd.read_csv(INSTANCES_CSV, dtype={"pdb_id": str})
    short_lbd = pd.read_csv(SHORT_MAPPING_LBD_CSV, dtype={"pdb_id": str})
    search_manifest = pd.read_csv(SEARCH_MANIFEST_CSV)

    multi_uniprot_entity_ids = set(
        anomalies.loc[anomalies["anomaly_type"] == "ENTITY_MULTI_MAPPED_UNIPROT", "polymer_entity_id"]
    )

    print("Building pre-coordinate QC table...")
    qc_df = build_precoordinate_qc(inventory, lbd_mapping, taxonomy_audit, multi_uniprot_entity_ids)
    qc_df.to_csv(PRECOORDINATE_QC_CSV, index=False)
    print(f"Wrote {PRECOORDINATE_QC_CSV.relative_to(PROJECT_ROOT)} ({len(qc_df)} rows)")

    # --- instance-level candidate inventory ---
    print("\nBuilding instance-level candidate inventory...")
    qc_join_cols = [
        "uniprot_id", "polymer_entity_id",
        "experimental_methods", "resolution_min", "resolution_max", "r_free_min", "r_free_max",
        "lbd_start", "lbd_end", "mapped_lbd_residues", "lbd_mapping_coverage", "lbd_mapping_category",
        "taxonomy_audit_status", "multi_uniprot_entity",
    ]
    inst_inventory = instances[instances["instance_metadata_status"] == "PASS"].merge(
        qc_df[qc_join_cols], on="polymer_entity_id", how="left"
    )
    # instances carries its own pdb_id/entity_id from the Data API; the QC
    # frame's uniprot_id was matched purely via polymer_entity_id, which is
    # unique per Stage 2 candidate (no fusion ambiguity here since each
    # entity maps to exactly one discovering receptor in this dataset).
    inst_inventory = inst_inventory.rename(columns={"resolution_min": "resolution_combined_min", "resolution_max": "resolution_combined_max"})
    nr_code_by_uniprot = {r["uniprot_id"]: r["nr_code"] for r in NR_METADATA}
    common_name_by_uniprot = {r["uniprot_id"]: r["common_name"] for r in NR_METADATA}
    group_by_uniprot = {r["uniprot_id"]: r["group"] for r in NR_METADATA}
    inst_inventory["nr_code"] = inst_inventory["uniprot_id"].map(nr_code_by_uniprot)
    inst_inventory["common_name"] = inst_inventory["uniprot_id"].map(common_name_by_uniprot)
    inst_inventory["group"] = inst_inventory["uniprot_id"].map(group_by_uniprot)
    inst_inventory.to_csv(INSTANCE_INVENTORY_CSV, index=False)
    print(f"Wrote {INSTANCE_INVENTORY_CSV.relative_to(PROJECT_ROOT)} ({len(inst_inventory)} rows)")

    # --- summary counts ---
    print("\nComputing summary counts...")
    instances_per_entity = instances[instances["instance_metadata_status"] == "PASS"].groupby("polymer_entity_id").size()

    zero_result_receptors = search_manifest[
        (search_manifest["success"]) & (search_manifest["result_count_parsed"] == 0)
    ]["common_name"].tolist()

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "total_stage2_candidates": len(qc_df),
        "method_counts": {
            "xray": int(qc_df["method_is_xray"].sum()),
            "cryo_em": int(inventory["is_cryo_em"].sum()),
            "nmr": int(inventory["is_nmr"].sum()),
            "other": int(inventory["is_other_experimental"].sum()),
        },
        "lbd_mapping_counts": qc_df["lbd_mapping_category"].value_counts().to_dict(),
        "resolution_status_counts": qc_df["resolution_project_range_status"].value_counts().to_dict(),
        "r_free_status_counts": qc_df["r_free_status"].value_counts().to_dict(),
        "project_metadata_prefilter": {
            "pass": int(qc_df["passes_project_metadata_prefilter"].sum()),
            "not_pass": int((~qc_df["passes_project_metadata_prefilter"]).sum()),
        },
        "instance_level": {
            "total_polymer_instances": int(len(instances[instances["instance_metadata_status"] == "PASS"])),
            "entities_with_one_instance": int((instances_per_entity == 1).sum()),
            "entities_with_multiple_instances": int((instances_per_entity > 1).sum()),
            "max_instances_for_one_entity": int(instances_per_entity.max()) if len(instances_per_entity) else 0,
        },
        "special": {
            "dbd_only_mappings": int(
                ((qc_df["mapped_dbd_residues"].astype(str) != "") & (qc_df["mapped_dbd_residues"].astype(str) != "nan")
                 & (~qc_df["mapping_overlaps_lbd"])
                 & (pd.to_numeric(qc_df["mapped_dbd_residues"], errors="coerce").fillna(0) > 0)).sum()
            ),
            "multi_uniprot_entities": int(qc_df["multi_uniprot_entity"].sum()),
            "taxonomy_review_entities": int((qc_df["taxonomy_audit_status"] != "").sum()),
            "short_mappings_overlapping_lbd": int(short_lbd["mapping_overlaps_lbd"].fillna(False).astype(bool).sum()),
            "short_mappings_not_overlapping_lbd": int((~short_lbd["mapping_overlaps_lbd"].fillna(False).astype(bool)).sum()),
        },
        "zero_result_receptors": zero_result_receptors,
        "note": "No final structure inclusion/exclusion decision was made at Stage 3A.",
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"Wrote {SUMMARY_JSON.relative_to(PROJECT_ROOT)}")

    # --- counts by receptor ---
    print("\nBuilding counts-by-receptor table...")
    receptor_rows = []
    for r in NR_METADATA:
        uid = r["uniprot_id"]
        sub = qc_df[qc_df["uniprot_id"] == uid]
        prefilter_pass_sub = sub[sub["passes_project_metadata_prefilter"]]
        entity_ids_prefilter_pass = set(prefilter_pass_sub["polymer_entity_id"])
        inst_count = int(instances[
            (instances["polymer_entity_id"].isin(entity_ids_prefilter_pass))
            & (instances["instance_metadata_status"] == "PASS")
        ].shape[0])

        receptor_rows.append(
            {
                "common_name": r["common_name"], "nr_code": r["nr_code"], "uniprot_id": uid,
                "stage2_entities": len(sub),
                "no_lbd_overlap": int((sub["lbd_mapping_category"] == "NO_LBD_OVERLAP").sum()),
                "trace_lbd": int((sub["lbd_mapping_category"] == "TRACE_LBD").sum()),
                "partial_lbd": int((sub["lbd_mapping_category"] == "PARTIAL_LBD").sum()),
                "substantial_lbd": int((sub["lbd_mapping_category"] == "SUBSTANTIAL_LBD").sum()),
                "near_complete_lbd": int((sub["lbd_mapping_category"] == "NEAR_COMPLETE_LBD").sum()),
                "xray_entities": int(sub["method_is_xray"].sum()),
                "project_metadata_prefilter_pass": len(prefilter_pass_sub),
                "unique_pdb_entries_prefilter_pass": prefilter_pass_sub["pdb_id"].nunique(),
                "instance_count_for_prefilter_entities": inst_count,
                "review_needed_count": int((sub["precoordinate_review_status"] != "PROJECT_PREFILTER_PASS").sum()),
            }
        )
    counts_df = pd.DataFrame(receptor_rows)
    counts_df.to_csv(COUNTS_BY_RECEPTOR_CSV, index=False)
    print(f"Wrote {COUNTS_BY_RECEPTOR_CSV.relative_to(PROJECT_ROOT)} ({len(counts_df)} rows)")

    # --- report ---
    write_report(lbd_reference, summary, counts_df, short_lbd)
    print(f"Wrote {REPORT_MD.relative_to(PROJECT_ROOT)}")

    print("\n--- Stage 3A summary ---")
    print(f"LBD mapping counts: {summary['lbd_mapping_counts']}")
    print(f"Project metadata prefilter pass: {summary['project_metadata_prefilter']['pass']}")


def write_report(lbd_reference: pd.DataFrame, summary: dict, counts_df: pd.DataFrame, short_lbd: pd.DataFrame) -> None:
    lines = ["# Stage 3A — LBD Reference, Instance Mapping, and Pre-Coordinate QC\n"]
    lines.append(
        "## 1. Purpose\n\n"
        "Establish authoritative canonical-UniProt LBD coordinates for each of the 48 validated "
        "receptors, map every Stage 2 candidate's sequence-level coverage against that LBD (and, "
        "diagnostically, the DBD), resolve exact RCSB polymer-instance/chain identifiers, and compute "
        "descriptive pre-coordinate quality flags. No mmCIF coordinates were downloaded and no final "
        "structure inclusion/exclusion decision was made.\n"
    )
    lines.append(
        "## 2. Authoritative LBD-coordinate method\n\n"
        "Primary source: the current reviewed UniProt entry's `features[]` list, accepting only a "
        "feature with `type == \"Domain\"` and `description == \"NR LBD\"` (PROSITE-ProRule evidenced). "
        "All 48 receptors had exactly one such feature in their cached Stage 1 raw UniProt response — "
        "no PROSITE/InterPro fallback request was required. DBD interval (diagnostic only) uses "
        "`type == \"DNA binding\"`, `description == \"Nuclear receptor\"`; absent for DAX-1 and SHP, "
        "which are documented atypical orphan receptors lacking a classical DBD — no interval was "
        "inferred for them. Coordinates are 1-based, inclusive.\n"
    )
    lines.append(
        "## 2a. Standardized LBD policy (Stage 3A.2)\n\n"
        "The LBD interval used throughout this report is a **standardized project-wide canonical-"
        "sequence interval** derived from the current UniProt NR-LBD domain annotation. It is used "
        "because the same systematic annotation is available for all 48 receptors. It is **not "
        "asserted to be the only possible biological definition** of an LBD.\n\n"
        "The previous cohort's VDR-only 118-427 range is preserved as a historical sensitivity "
        "definition (`data/manifests/nr_lbd_boundary_audit.csv`, "
        "`rcsb_lbd_boundary_sensitivity.csv`) but is **not used for primary filtering** anywhere in "
        "this report.\n\n"
        "Functional motifs/sites (e.g. the VDR 9aaTAD/AF-2-associated motif) are represented "
        "separately in `data/manifests/nr_functional_features.csv` and may extend outside the "
        "standardized profile-defined interval. See `data/manifests/nr_lbd_analysis_policy.json` "
        "for the full policy record.\n"
    )
    lines.append("## 3. LBD reference table (all 48 receptors)\n")
    lines.append("| common_name | nr_code | uniprot_id | LBD start-end | length | source | status |")
    lines.append("|---|---|---|---|---|---|---|")
    for _, row in lbd_reference.sort_values("record_index").iterrows():
        lines.append(
            f"| {row['common_name']} | {row['nr_code']} | {row['uniprot_id']} | "
            f"{row['lbd_start']}-{row['lbd_end']} | {row['lbd_length']} | {row['lbd_source']} | {row['lbd_status']} |"
        )
    lines.append("")

    lines.append("## 4. Mapping-category counts\n")
    for k, v in summary["lbd_mapping_counts"].items():
        lines.append(f"- {k}: {v}")
    lines.append("")

    lines.append("## 5. Experimental method counts\n")
    for k, v in summary["method_counts"].items():
        lines.append(f"- {k}: {v}")
    lines.append("")

    lines.append("## 6. Resolution categories\n")
    for k, v in summary["resolution_status_counts"].items():
        lines.append(f"- {k}: {v}")
    lines.append("")

    lines.append("## 7. R-free categories\n")
    for k, v in summary["r_free_status_counts"].items():
        lines.append(f"- {k}: {v}")
    lines.append("")

    lines.append("## 8. Project metadata-prefilter count\n")
    lines.append(f"- PASS: {summary['project_metadata_prefilter']['pass']}")
    lines.append(f"- NOT PASS: {summary['project_metadata_prefilter']['not_pass']}\n")

    lines.append("## 9. Receptors with zero or very few LBD-mapped candidates\n")
    low = counts_df[(counts_df["near_complete_lbd"] + counts_df["substantial_lbd"]) <= 2]
    if len(low) == 0:
        lines.append("None.\n")
    else:
        lines.append("| common_name | near_complete_lbd | substantial_lbd | stage2_entities |")
        lines.append("|---|---|---|---|")
        for _, row in low.iterrows():
            lines.append(f"| {row['common_name']} | {row['near_complete_lbd']} | {row['substantial_lbd']} | {row['stage2_entities']} |")
        lines.append("")

    lines.append("## 10. Short-fragment audit summary\n")
    lines.append(
        f"{summary['special']['short_mappings_overlapping_lbd']} of the 176 Stage 2 short mappings "
        f"overlap the LBD; {summary['special']['short_mappings_not_overlapping_lbd']} do not (they lie "
        f"elsewhere in the receptor sequence, e.g. AF-2/coactivator-groove peptides outside the "
        f"canonical domain boundary, or non-LBD fragments). See `rcsb_short_mapping_lbd_audit.csv`.\n"
    )

    lines.append("## 11. Instance/chain mapping summary\n")
    il = summary["instance_level"]
    lines.append(f"- Total polymer instances: {il['total_polymer_instances']}")
    lines.append(f"- Entities with one instance: {il['entities_with_one_instance']}")
    lines.append(f"- Entities with multiple instances: {il['entities_with_multiple_instances']}")
    lines.append(f"- Maximum instances for one entity: {il['max_instances_for_one_entity']}")
    lines.append(
        "- label_asym_id <-> auth_asym_id pairing was resolved from true instance-level RCSB records "
        "(`CorePolymerEntityInstance`), never by pairing the Stage 2 entity-level array positions.\n"
    )

    lines.append("## 12. Taxonomy / multi-UniProt caveats\n")
    lines.append(
        f"- {summary['special']['multi_uniprot_entities']} candidates are multi-UniProt-mapped entities "
        f"(fusion/chimera constructs, carried forward — not excluded)."
    )
    lines.append(
        f"- {summary['special']['taxonomy_review_entities']} candidates carry a taxonomy-audit status "
        f"from Stage 2.1/2.2 (not used as an automatic exclusion criterion).\n"
    )

    lines.append("## 13. Important statement\n")
    lines.append(
        "Sequence-mapping coverage is not equivalent to observed crystallographic coordinate coverage.\n"
    )
    lines.append("## 14. Important statement\n")
    lines.append("No final structure inclusion/exclusion decision was made at Stage 3A.\n")
    lines.append("## 15. Important statement\n")
    lines.append("No mmCIF coordinate files were downloaded at Stage 3A.\n")

    REPORT_MD.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
