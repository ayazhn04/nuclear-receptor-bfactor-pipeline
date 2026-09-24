"""
Stage 2, step 3: build the master candidate inventory by joining validated
receptor metadata, Search API discovery, and Data API metadata for polymer
entities and their parent PDB entries.

Reads:
    config/nr_metadata.py
    data/manifests/rcsb_search_manifest.csv
    data/manifests/rcsb_discovery_by_receptor.csv
    data/raw/api/rcsb/data/polymer_entities_batch_*.json
    data/raw/api/rcsb/data/entries_batch_*.json
    data/raw/api/rcsb/holdings/current_entry_ids.json  (fetched here if absent)

Writes:
    data/manifests/rcsb_polymer_entities.csv
    data/manifests/rcsb_entries.csv
    data/manifests/rcsb_candidate_inventory.csv
    data/manifests/rcsb_discovery_anomalies.csv
    data/manifests/rcsb_discovery_summary.json
    data/manifests/rcsb_stage2_provenance.json
    reports/tables/rcsb_discovery_report.md

Note: the RCSB Data API silently OMITS unknown/unrecognized IDs from its
response array rather than returning a null placeholder (confirmed via a
pilot test with a deliberately invalid ID). This script therefore matches
returned records by their `rcsb_id` field, not by request-list position,
and explicitly flags any requested ID that never comes back in any batch.

This script performs NO network access beyond the optional one-shot
holdings fetch (reused if already cached) and downloads NO coordinate files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402
from scripts.utils.rcsb_extract import (  # noqa: E402
    compute_mapping_status,
    compute_method_flags,
    extract_entry_fields,
    extract_polymer_entity_fields,
)

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
REPORTS_TABLES_DIR = PROJECT_ROOT / "reports" / "tables"
RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "data"
RAW_HOLDINGS_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "holdings"
HOLDINGS_URL = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"

SEARCH_MANIFEST_PATH = MANIFESTS_DIR / "rcsb_search_manifest.csv"
DISCOVERY_BY_RECEPTOR_PATH = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"

POLYMER_ENTITIES_CSV = MANIFESTS_DIR / "rcsb_polymer_entities.csv"
ENTRIES_CSV = MANIFESTS_DIR / "rcsb_entries.csv"
CANDIDATE_INVENTORY_CSV = MANIFESTS_DIR / "rcsb_candidate_inventory.csv"
ANOMALIES_CSV = MANIFESTS_DIR / "rcsb_discovery_anomalies.csv"
SUMMARY_JSON = MANIFESTS_DIR / "rcsb_discovery_summary.json"
PROVENANCE_JSON = MANIFESTS_DIR / "rcsb_stage2_provenance.json"
REPORT_MD = REPORTS_TABLES_DIR / "rcsb_discovery_report.md"


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_field(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


# --------------------------------------------------------------------------
# Load raw Data API batches, matching by returned rcsb_id (NOT position).
# --------------------------------------------------------------------------

def load_polymer_entity_records() -> tuple[dict[str, dict], set[str]]:
    """Returns (records_by_id, all_requested_ids)."""
    records: dict[str, dict] = {}
    requested: set[str] = set()
    for path in sorted(RAW_DATA_DIR.glob("polymer_entities_batch_*.json")):
        wrapper = json.loads(path.read_text())
        requested.update(wrapper.get("ids_requested", []))
        if not wrapper.get("success"):
            continue
        response = wrapper.get("response") or {}
        entities = (response.get("data") or {}).get("polymer_entities") or []
        for e in entities:
            if e and e.get("rcsb_id"):
                records[e["rcsb_id"]] = e
    return records, requested


def load_entry_records() -> tuple[dict[str, dict], set[str]]:
    records: dict[str, dict] = {}
    requested: set[str] = set()
    for path in sorted(RAW_DATA_DIR.glob("entries_batch_*.json")):
        wrapper = json.loads(path.read_text())
        requested.update(wrapper.get("ids_requested", []))
        if not wrapper.get("success"):
            continue
        response = wrapper.get("response") or {}
        entries = (response.get("data") or {}).get("entries") or []
        for e in entries:
            if e and e.get("rcsb_id"):
                records[e["rcsb_id"]] = e
    return records, requested


# --------------------------------------------------------------------------
# Holdings check (optional; used if reachable, never blocks the pipeline)
# --------------------------------------------------------------------------

def get_current_holdings(refresh: bool) -> set[str] | None:
    RAW_HOLDINGS_DIR.mkdir(parents=True, exist_ok=True)
    raw_path = RAW_HOLDINGS_DIR / "current_entry_ids.json"

    if not refresh and raw_path.exists():
        try:
            wrapper = json.loads(raw_path.read_text())
            return set(wrapper["response"])
        except (json.JSONDecodeError, KeyError):
            pass

    try:
        response = requests.get(HOLDINGS_URL, timeout=60)
        response.raise_for_status()
        ids = response.json()
    except (requests.exceptions.RequestException, ValueError) as exc:
        print(f"  WARNING: could not retrieve current holdings ({exc}); "
              f"is_current_pdb_entry will be left blank for all entries.")
        return None

    wrapper = {
        "request_url": HOLDINGS_URL,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "response": ids,
    }
    raw_path.write_text(json.dumps(wrapper))
    return set(ids)


# --------------------------------------------------------------------------
# Build unique polymer entity table
# --------------------------------------------------------------------------

def build_polymer_entities_table(
    discovery_df: pd.DataFrame,
    entity_records: dict[str, dict],
    requested_entity_ids: set[str],
) -> pd.DataFrame:
    unique_entity_ids = sorted(discovery_df["polymer_entity_id"].unique())
    # A discovery row's query_uniprot_id tells us which accession(s) discovered
    # this entity; for mapping verification we use the FIRST discovering
    # query (an entity can be discovered by more than one query — see the
    # cross-mapping anomaly detection below).
    first_query_by_entity = (
        discovery_df.sort_values("query_retrieved_at_utc")
        .groupby("polymer_entity_id")["query_uniprot_id"]
        .first()
        .to_dict()
    )

    rows = []
    for entity_id in unique_entity_ids:
        raw = entity_records.get(entity_id)
        if raw is None:
            rows.append(
                {
                    "polymer_entity_id": entity_id,
                    "pdb_id": entity_id.split("_", 1)[0],
                    "entity_id": entity_id.split("_", 1)[1] if "_" in entity_id else "",
                    "entity_description": "",
                    "polymer_type": "",
                    "entity_sequence_length": "",
                    "source_organism_names": json_field([]),
                    "source_taxonomy_ids": json_field([]),
                    "label_asym_ids": json_field([]),
                    "auth_asym_ids": json_field([]),
                    "all_mapped_uniprot_accessions": json_field([]),
                    "mapping_status": "REVIEW_NEEDED",
                    "mapping_notes": (
                        "Entity was returned by the Search API but the Data API never "
                        "returned a record for it in any batch (silently omitted, not null) "
                        "— metadata unavailable."
                    ),
                }
            )
            continue

        extracted = extract_polymer_entity_fields(raw)
        query_uniprot_id = first_query_by_entity.get(entity_id, "")
        _, mapping_status, mapping_notes = compute_mapping_status(query_uniprot_id, extracted)

        rows.append(
            {
                "polymer_entity_id": entity_id,
                "pdb_id": extracted["entry_id"],
                "entity_id": extracted["entity_id"],
                "entity_description": extracted["entity_description"],
                "polymer_type": extracted["polymer_type"],
                "entity_sequence_length": extracted["entity_sequence_length"],
                "source_organism_names": json_field(extracted["source_organism_names"]),
                "source_taxonomy_ids": json_field(extracted["source_taxonomy_ids"]),
                "label_asym_ids": json_field(extracted["label_asym_ids"]),
                "auth_asym_ids": json_field(extracted["auth_asym_ids"]),
                "all_mapped_uniprot_accessions": json_field(extracted["all_mapped_uniprot_accessions"]),
                "mapping_status": mapping_status,
                "mapping_notes": mapping_notes,
            }
        )

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Build unique entry table
# --------------------------------------------------------------------------

def build_entries_table(
    discovery_df: pd.DataFrame,
    entry_records: dict[str, dict],
    holdings: set[str] | None,
) -> pd.DataFrame:
    unique_pdb_ids = sorted(discovery_df["pdb_id"].unique())

    rows = []
    for pdb_id in unique_pdb_ids:
        raw = entry_records.get(pdb_id)
        is_current = (pdb_id in holdings) if holdings is not None else ""
        if raw is None:
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "structure_title": "",
                    "structure_determination_methodology": "",
                    "experimental_methods": json_field([]),
                    "resolution_combined": json_field([]),
                    "r_work": json_field([]),
                    "r_free": json_field([]),
                    "deposit_date": "",
                    "initial_release_date": "",
                    "revision_date": "",
                    "space_group": "",
                    "polymer_entity_count": "",
                    "protein_entity_count": "",
                    "nonpolymer_entity_count": "",
                    "deposited_polymer_instance_count": "",
                    "primary_citation_title": "",
                    "primary_citation_doi": "",
                    "primary_citation_year": "",
                    "is_current_pdb_entry": is_current,
                    "data_api_metadata_missing": True,
                }
            )
            continue

        extracted = extract_entry_fields(raw)
        rows.append(
            {
                "pdb_id": pdb_id,
                "structure_title": extracted["structure_title"],
                "structure_determination_methodology": extracted["structure_determination_methodology"],
                "experimental_methods": json_field(extracted["experimental_methods"]),
                "resolution_combined": json_field(extracted["resolution_combined"]),
                "r_work": json_field(extracted["r_work"]),
                "r_free": json_field(extracted["r_free"]),
                "deposit_date": extracted["deposit_date"],
                "initial_release_date": extracted["initial_release_date"],
                "revision_date": extracted["revision_date"],
                "space_group": extracted["space_group"],
                "polymer_entity_count": extracted["polymer_entity_count"],
                "protein_entity_count": extracted["protein_entity_count"],
                "nonpolymer_entity_count": extracted["nonpolymer_entity_count"],
                "deposited_polymer_instance_count": extracted["deposited_polymer_instance_count"],
                "primary_citation_title": extracted["primary_citation_title"],
                "primary_citation_doi": extracted["primary_citation_doi"],
                "primary_citation_year": extracted["primary_citation_year"],
                "is_current_pdb_entry": is_current,
                "data_api_metadata_missing": False,
            }
        )

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Candidate inventory (receptor query x polymer entity, joined)
# --------------------------------------------------------------------------

def build_candidate_inventory(
    discovery_df: pd.DataFrame,
    entities_df: pd.DataFrame,
    entries_df: pd.DataFrame,
) -> pd.DataFrame:
    merged = discovery_df.merge(entities_df, on="polymer_entity_id", how="left", suffixes=("", "_entity"))
    merged = merged.merge(entries_df, on="pdb_id", how="left", suffixes=("", "_entry"))

    def parse_json_col(col):
        return col.apply(lambda v: json.loads(v) if isinstance(v, str) and v else [])

    experimental_methods_parsed = parse_json_col(merged["experimental_methods"])
    method_flags = experimental_methods_parsed.apply(compute_method_flags)

    merged["is_xray"] = method_flags.apply(lambda d: d["is_xray"])
    merged["is_cryo_em"] = method_flags.apply(lambda d: d["is_cryo_em"])
    merged["is_nmr"] = method_flags.apply(lambda d: d["is_nmr"])
    merged["is_other_experimental"] = method_flags.apply(lambda d: d["is_other_experimental"])

    resolution_parsed = parse_json_col(merged["resolution_combined"])
    r_free_parsed = parse_json_col(merged["r_free"])
    merged["has_resolution"] = resolution_parsed.apply(lambda v: len(v) > 0)
    merged["has_r_free"] = r_free_parsed.apply(lambda v: len(v) > 0)

    mapped_accessions_parsed = parse_json_col(merged["all_mapped_uniprot_accessions"])
    merged["mapping_contains_query_uniprot"] = [
        q in mapped for q, mapped in zip(merged["query_uniprot_id"], mapped_accessions_parsed)
    ]

    taxon_ids_parsed = parse_json_col(merged["source_taxonomy_ids"])
    merged["source_includes_homo_sapiens_9606"] = taxon_ids_parsed.apply(lambda v: 9606 in v)

    columns = [
        "uniprot_id", "nr_code", "common_name", "group",
        "pdb_id", "polymer_entity_id", "entity_id",
        "entity_description", "entity_sequence_length",
        "label_asym_ids", "auth_asym_ids",
        "all_mapped_uniprot_accessions", "mapping_contains_query_uniprot", "mapping_status",
        "source_organism_names", "source_taxonomy_ids",
        "experimental_methods", "structure_determination_methodology",
        "resolution_combined", "r_work", "r_free",
        "deposit_date", "initial_release_date",
        "is_xray", "is_cryo_em", "is_nmr", "is_other_experimental",
        "has_resolution", "has_r_free", "source_includes_homo_sapiens_9606",
        "is_current_pdb_entry",
    ]
    return merged[columns]


# --------------------------------------------------------------------------
# Anomaly detection
# --------------------------------------------------------------------------

def build_anomalies_table(
    discovery_df: pd.DataFrame,
    entities_df: pd.DataFrame,
    search_manifest_df: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    # A. same polymer entity returned for more than one receptor query
    entity_to_receptors = discovery_df.groupby("polymer_entity_id")["uniprot_id"].apply(lambda s: sorted(set(s)))
    for entity_id, receptors in entity_to_receptors.items():
        if len(receptors) > 1:
            rows.append(
                {
                    "anomaly_type": "CROSS_MAPPED_ENTITY_MULTIPLE_RECEPTOR_QUERIES",
                    "uniprot_id": json_field(receptors),
                    "nr_code": "",
                    "pdb_id": entity_id.split("_", 1)[0],
                    "polymer_entity_id": entity_id,
                    "details": f"Polymer entity {entity_id} was returned by {len(receptors)} distinct receptor queries: {receptors}",
                    "review_status": "REVIEW_NEEDED",
                }
            )

    # B. polymer entities mapping to more than one UniProt accession
    for _, row in entities_df.iterrows():
        mapped = json.loads(row["all_mapped_uniprot_accessions"]) if row["all_mapped_uniprot_accessions"] else []
        if len(mapped) > 1:
            rows.append(
                {
                    "anomaly_type": "ENTITY_MULTI_MAPPED_UNIPROT",
                    "uniprot_id": json_field(mapped),
                    "nr_code": "",
                    "pdb_id": row["pdb_id"],
                    "polymer_entity_id": row["polymer_entity_id"],
                    "details": f"Entity maps to multiple UniProt accessions: {mapped}",
                    "review_status": "REVIEW_NEEDED",
                }
            )

    # C. same PDB entry containing more than one of our 48 receptor UniProt accessions
    entry_to_receptors = discovery_df.groupby("pdb_id")["uniprot_id"].apply(lambda s: sorted(set(s)))
    for pdb_id, receptors in entry_to_receptors.items():
        if len(receptors) > 1:
            rows.append(
                {
                    "anomaly_type": "PDB_ENTRY_MULTIPLE_PROJECT_RECEPTORS",
                    "uniprot_id": json_field(receptors),
                    "nr_code": "",
                    "pdb_id": pdb_id,
                    "polymer_entity_id": "",
                    "details": f"PDB entry {pdb_id} contains polymer entities mapped to {len(receptors)} distinct project receptors: {receptors}",
                    "review_status": "REVIEW_NEEDED",
                }
            )

    # D & E. zero-result / low-count receptors (from search manifest, successful queries only)
    for _, row in search_manifest_df.iterrows():
        if not row["success"]:
            continue
        count = row["result_count_parsed"]
        if count == 0:
            rows.append(
                {
                    "anomaly_type": "ZERO_CANDIDATE_RECEPTOR",
                    "uniprot_id": row["uniprot_id"],
                    "nr_code": row["nr_code"],
                    "pdb_id": "",
                    "polymer_entity_id": "",
                    "details": f"{row['common_name']} ({row['nr_code']}): 0 experimental polymer entities found for a QUERY_SUCCEEDED_ZERO_RESULTS search.",
                    "review_status": "REVIEW_NEEDED",
                }
            )
        elif count in (1, 2):
            rows.append(
                {
                    "anomaly_type": "LOW_CANDIDATE_COUNT_RECEPTOR",
                    "uniprot_id": row["uniprot_id"],
                    "nr_code": row["nr_code"],
                    "pdb_id": "",
                    "polymer_entity_id": "",
                    "details": f"{row['common_name']} ({row['nr_code']}): only {count} experimental polymer entity/entities found.",
                    "review_status": "REVIEW_NEEDED",
                }
            )

    # F. failed API retrievals (search-level)
    for _, row in search_manifest_df.iterrows():
        if not row["success"]:
            rows.append(
                {
                    "anomaly_type": "QUERY_FAILED",
                    "uniprot_id": row["uniprot_id"],
                    "nr_code": row["nr_code"],
                    "pdb_id": "",
                    "polymer_entity_id": "",
                    "details": f"Search API query failed for {row['common_name']} ({row['nr_code']}): {row['error_message']}",
                    "review_status": "REVIEW_NEEDED",
                }
            )

    # F (cont). mapping-status REVIEW_NEEDED entities from the unique entity table
    for _, row in entities_df.iterrows():
        if row["mapping_status"] == "REVIEW_NEEDED":
            rows.append(
                {
                    "anomaly_type": "UNCONFIRMED_UNIPROT_MAPPING",
                    "uniprot_id": "",
                    "nr_code": "",
                    "pdb_id": row["pdb_id"],
                    "polymer_entity_id": row["polymer_entity_id"],
                    "details": row["mapping_notes"],
                    "review_status": "REVIEW_NEEDED",
                }
            )

    columns = ["anomaly_type", "uniprot_id", "nr_code", "pdb_id", "polymer_entity_id", "details", "review_status"]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows)[columns]


# --------------------------------------------------------------------------
# Summary + report
# --------------------------------------------------------------------------

def build_summary(
    search_manifest_df: pd.DataFrame,
    discovery_df: pd.DataFrame,
    entities_df: pd.DataFrame,
    entries_df: pd.DataFrame,
    inventory_df: pd.DataFrame,
    anomalies_df: pd.DataFrame,
) -> dict:
    counts_by_receptor = (
        discovery_df.groupby(["uniprot_id", "nr_code", "common_name"])
        .agg(polymer_entity_count=("polymer_entity_id", "nunique"), unique_pdb_entry_count=("pdb_id", "nunique"))
        .reset_index()
    )
    zero_result_receptors = search_manifest_df[
        (search_manifest_df["success"]) & (search_manifest_df["result_count_parsed"] == 0)
    ]["common_name"].tolist()
    counts_1_2 = counts_by_receptor[
        counts_by_receptor["polymer_entity_count"].isin([1, 2])
    ]["common_name"].tolist()
    counts_ge3 = counts_by_receptor[counts_by_receptor["polymer_entity_count"] >= 3]["common_name"].tolist()

    method_counts = {
        "xray": int(inventory_df["is_xray"].sum()),
        "cryo_em": int(inventory_df["is_cryo_em"].sum()),
        "nmr": int(inventory_df["is_nmr"].sum()),
        "other_experimental": int(inventory_df["is_other_experimental"].sum()),
    }

    return {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "number_of_validated_receptors_queried": 48,
        "number_of_successful_receptor_queries": int(search_manifest_df["success"].sum()),
        "number_of_failed_receptor_queries": int((~search_manifest_df["success"]).sum()),
        "total_receptor_polymer_entity_relationships": len(discovery_df),
        "unique_polymer_entities": int(discovery_df["polymer_entity_id"].nunique()),
        "unique_pdb_entries": int(discovery_df["pdb_id"].nunique()),
        "receptors_with_zero_candidates": len(zero_result_receptors),
        "receptors_with_zero_candidates_list": zero_result_receptors,
        "receptors_with_1_to_2_candidates": len(counts_1_2),
        "receptors_with_1_to_2_candidates_list": counts_1_2,
        "receptors_with_3_or_more_candidates": len(counts_ge3),
        "candidate_counts_by_receptor": counts_by_receptor.to_dict(orient="records"),
        "experimental_method_counts_across_candidate_rows": method_counts,
        "candidates_with_resolution": int(inventory_df["has_resolution"].sum()),
        "candidates_without_resolution": int((~inventory_df["has_resolution"]).sum()),
        "candidates_with_r_free": int(inventory_df["has_r_free"].sum()),
        "candidates_without_r_free": int((~inventory_df["has_r_free"]).sum()),
        "candidates_source_includes_taxon_9606": int(inventory_df["source_includes_homo_sapiens_9606"].sum()),
        "candidates_uniprot_mapping_confirmed": int((entities_df["mapping_status"].isin(["CONFIRMED", "CONFIRMED_MULTI_MAPPED"])).sum()),
        "candidates_uniprot_mapping_review_needed": int((entities_df["mapping_status"] == "REVIEW_NEEDED").sum()),
        "number_of_cross_mapped_entities": int((anomalies_df["anomaly_type"] == "CROSS_MAPPED_ENTITY_MULTIPLE_RECEPTOR_QUERIES").sum()),
        "number_of_pdb_entries_with_multiple_project_receptors": int((anomalies_df["anomaly_type"] == "PDB_ENTRY_MULTIPLE_PROJECT_RECEPTORS").sum()),
        "number_of_anomalies_total": len(anomalies_df),
        "note": "This is a discovery-stage inventory. No final structural QC filters were applied.",
    }


def write_report(
    search_manifest_df: pd.DataFrame,
    counts_by_receptor: pd.DataFrame,
    summary: dict,
    anomalies_df: pd.DataFrame,
) -> None:
    lines = ["# RCSB Candidate Structure Discovery — Stage 2\n"]
    lines.append(
        "## Methods\n\n"
        "For each of the 48 UniProt-validated human nuclear receptors "
        "(`config/nr_metadata.py`), the official RCSB Search API "
        "(`https://search.rcsb.org/rcsbsearch/v2/query`) was queried for "
        "current, EXPERIMENTAL polymer entities whose "
        "`reference_sequence_identifiers` mapping includes the receptor's "
        "validated UniProt accession. Computed structure models were "
        "excluded via `rcsb_entry_info.structure_determination_methodology "
        "== \"experimental\"`; no restriction to a specific experimental "
        "method (X-ray/cryo-EM/NMR) was applied at this stage. Discovered "
        "polymer entities and their parent PDB entries were then enriched "
        "via the official RCSB Data API (GraphQL, batched). No mmCIF or PDB "
        "coordinate files were downloaded at any point.\n"
    )
    lines.append("## Overall counts\n")
    for key in (
        "number_of_validated_receptors_queried",
        "number_of_successful_receptor_queries",
        "number_of_failed_receptor_queries",
        "total_receptor_polymer_entity_relationships",
        "unique_polymer_entities",
        "unique_pdb_entries",
        "receptors_with_zero_candidates",
        "receptors_with_1_to_2_candidates",
        "receptors_with_3_or_more_candidates",
    ):
        lines.append(f"- {key}: {summary[key]}")
    lines.append("")

    lines.append("## All 48 receptors\n")
    lines.append("| common_name | nr_code | uniprot_id | polymer_entity_count | unique_pdb_entry_count |")
    lines.append("|---|---|---|---|---|")
    merged = search_manifest_df.merge(
        counts_by_receptor, on=["uniprot_id", "nr_code", "common_name"], how="left"
    )
    merged["polymer_entity_count"] = merged["polymer_entity_count"].fillna(0).astype(int)
    merged["unique_pdb_entry_count"] = merged["unique_pdb_entry_count"].fillna(0).astype(int)
    for _, row in merged.sort_values("record_index").iterrows():
        lines.append(
            f"| {row['common_name']} | {row['nr_code']} | {row['uniprot_id']} | "
            f"{row['polymer_entity_count']} | {row['unique_pdb_entry_count']} |"
        )
    lines.append("")

    lines.append("## Receptors with zero candidates\n")
    zero_list = summary["receptors_with_zero_candidates_list"]
    lines.append(", ".join(zero_list) if zero_list else "None.")
    lines.append("")

    lines.append("## Receptors with 1-2 candidates\n")
    low_list = summary["receptors_with_1_to_2_candidates_list"]
    lines.append(", ".join(low_list) if low_list else "None.")
    lines.append("")

    lines.append("## Experimental-method breakdown (candidate rows)\n")
    mc = summary["experimental_method_counts_across_candidate_rows"]
    lines.append(f"- X-ray: {mc['xray']}")
    lines.append(f"- cryo-EM: {mc['cryo_em']}")
    lines.append(f"- NMR: {mc['nmr']}")
    lines.append(f"- Other experimental: {mc['other_experimental']}")
    lines.append("")

    lines.append("## Mapping anomalies\n")
    mapping_anomalies = anomalies_df[anomalies_df["anomaly_type"] == "UNCONFIRMED_UNIPROT_MAPPING"]
    lines.append(f"{len(mapping_anomalies)} entities flagged REVIEW_NEEDED for unconfirmed UniProt mapping.\n")

    lines.append("## Cross-receptor PDB entries\n")
    cross_entries = anomalies_df[anomalies_df["anomaly_type"] == "PDB_ENTRY_MULTIPLE_PROJECT_RECEPTORS"]
    if len(cross_entries) == 0:
        lines.append("None detected.\n")
    else:
        lines.append("| pdb_id | receptors (uniprot_id) |")
        lines.append("|---|---|")
        for _, row in cross_entries.iterrows():
            lines.append(f"| {row['pdb_id']} | {row['uniprot_id']} |")
        lines.append("")

    lines.append("## Important statements\n")
    lines.append("No final structural QC filters were applied at Stage 2.\n")
    lines.append("No mmCIF coordinate files were downloaded at Stage 2.\n")

    REPORT_MD.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-holdings", action="store_true", help="Force re-fetching the current holdings list.")
    args = parser.parse_args()

    search_manifest_df = pd.read_csv(SEARCH_MANIFEST_PATH)
    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR_PATH, dtype={"pdb_id": str, "entity_id": str})

    entity_records, requested_entity_ids = load_polymer_entity_records()
    entry_records, requested_pdb_ids = load_entry_records()

    print(f"Loaded {len(entity_records)} polymer entity records "
          f"({len(requested_entity_ids) - len(entity_records)} requested-but-missing).")
    print(f"Loaded {len(entry_records)} entry records "
          f"({len(requested_pdb_ids) - len(entry_records)} requested-but-missing).")

    print("\nChecking current PDB holdings...")
    holdings = get_current_holdings(refresh=args.refresh_holdings)
    if holdings is not None:
        print(f"  {len(holdings)} current entries in holdings snapshot.")

    entities_df = build_polymer_entities_table(discovery_df, entity_records, requested_entity_ids)
    entities_df.to_csv(POLYMER_ENTITIES_CSV, index=False)
    print(f"\nWrote {POLYMER_ENTITIES_CSV.relative_to(PROJECT_ROOT)} ({len(entities_df)} unique polymer entities)")

    entries_df = build_entries_table(discovery_df, entry_records, holdings)
    entries_df.to_csv(ENTRIES_CSV, index=False)
    print(f"Wrote {ENTRIES_CSV.relative_to(PROJECT_ROOT)} ({len(entries_df)} unique PDB entries)")

    inventory_df = build_candidate_inventory(discovery_df, entities_df, entries_df)
    inventory_df.to_csv(CANDIDATE_INVENTORY_CSV, index=False)
    print(f"Wrote {CANDIDATE_INVENTORY_CSV.relative_to(PROJECT_ROOT)} ({len(inventory_df)} rows)")

    anomalies_df = build_anomalies_table(discovery_df, entities_df, search_manifest_df)
    anomalies_df.to_csv(ANOMALIES_CSV, index=False)
    print(f"Wrote {ANOMALIES_CSV.relative_to(PROJECT_ROOT)} ({len(anomalies_df)} anomaly rows)")

    summary = build_summary(search_manifest_df, discovery_df, entities_df, entries_df, inventory_df, anomalies_df)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"Wrote {SUMMARY_JSON.relative_to(PROJECT_ROOT)}")

    counts_by_receptor = (
        discovery_df.groupby(["uniprot_id", "nr_code", "common_name"])
        .agg(polymer_entity_count=("polymer_entity_id", "nunique"), unique_pdb_entry_count=("pdb_id", "nunique"))
        .reset_index()
    )
    write_report(search_manifest_df, counts_by_receptor, summary, anomalies_df)
    print(f"Wrote {REPORT_MD.relative_to(PROJECT_ROOT)}")

    provenance = {
        "project_name": "Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families",
        "description": (
            "Stage 2: discovery of current experimental RCSB polymer-entity "
            "candidates for the 48 validated nuclear receptors, via the "
            "official RCSB Search API and Data API. Inventory only — no "
            "final structural QC filtering, no coordinate downloads."
        ),
        "validated_receptor_master_sha256": hashlib.sha256(
            (PROJECT_ROOT / "config" / "nr_metadata.py").read_bytes()
        ).hexdigest(),
        "rcsb_search_api_endpoint": "https://search.rcsb.org/rcsbsearch/v2/query",
        "rcsb_data_api_endpoint": "https://data.rcsb.org/graphql",
        "rcsb_holdings_endpoint": HOLDINGS_URL if holdings is not None else None,
        "retrieval_window_utc": {
            "search_manifest_min": str(search_manifest_df["retrieved_at_utc"].min()),
            "search_manifest_max": str(search_manifest_df["retrieved_at_utc"].max()),
        },
        "number_of_search_api_queries": len(search_manifest_df),
        "number_of_data_api_batches": len(list(RAW_DATA_DIR.glob("*_batch_*.json"))),
        "generated_manifests": {
            "rcsb_search_manifest.csv": sha256_of_file(SEARCH_MANIFEST_PATH),
            "rcsb_discovery_by_receptor.csv": sha256_of_file(DISCOVERY_BY_RECEPTOR_PATH),
            "rcsb_polymer_entities.csv": sha256_of_file(POLYMER_ENTITIES_CSV),
            "rcsb_entries.csv": sha256_of_file(ENTRIES_CSV),
            "rcsb_candidate_inventory.csv": sha256_of_file(CANDIDATE_INVENTORY_CSV),
            "rcsb_discovery_anomalies.csv": sha256_of_file(ANOMALIES_CSV),
        },
        "python_version": sys.version,
        "generating_scripts": [
            "scripts/04_discover_rcsb_candidates.py",
            "scripts/05_fetch_rcsb_metadata.py",
            "scripts/06_build_rcsb_inventory.py",
        ],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    PROVENANCE_JSON.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Wrote {PROVENANCE_JSON.relative_to(PROJECT_ROOT)}")

    print("\n--- Stage 2 discovery summary ---")
    for key in (
        "number_of_successful_receptor_queries", "number_of_failed_receptor_queries",
        "total_receptor_polymer_entity_relationships", "unique_polymer_entities", "unique_pdb_entries",
        "receptors_with_zero_candidates", "receptors_with_1_to_2_candidates",
        "candidates_uniprot_mapping_confirmed", "candidates_uniprot_mapping_review_needed",
        "number_of_anomalies_total",
    ):
        print(f"{key}: {summary[key]}")


if __name__ == "__main__":
    main()
