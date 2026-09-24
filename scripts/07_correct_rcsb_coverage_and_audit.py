"""
Stage 2.1: correct the reference/entity sequence coverage representation in
the Stage 2 outputs, and audit the source-taxonomy and multi-UniProt-mapping
anomalies flagged during Stage 2.

Background / correction:
    Stage 2's rcsb_polymer_entities.csv incorrectly stored the SAME raw
    alignment-region JSON in both a `reference_sequence_coverage` and an
    `entity_sequence_coverage` column. These are different quantities, and
    GraphQL introspection against https://data.rcsb.org/graphql (see
    reports/tables/rcsb_discovery_report.md "Methods") confirmed the RCSB
    Data API does NOT expose a precomputed coverage-fraction field at all —
    `RcsbPolymerEntityAlign` only has `aligned_regions` (entity_beg_seq_id,
    ref_beg_seq_id, length). This script computes both coverage fractions
    deterministically from data already retrieved:
      - entity_sequence_coverage = (aligned length) / (entity's own
        sequence length) — always computable from already-fetched data.
      - reference_sequence_coverage = (aligned length) / (reference
        accession's UniProt sequence length) — computable ONLY for the
        project's 48 validated receptors, whose sequence lengths were
        already retrieved in Stage 1 (data/manifests/nr_metadata_uniprot_audit.csv).
        For a fusion/chimera partner accession (e.g. a coactivator peptide)
        outside the 48-receptor master, no sequence length is available in
        this project's data, and NO additional UniProt request is made
        for it — the value is left null rather than guessed or fetched
        outside this stage's declared scope.
    No new RCSB Data API fields were needed for coverage. ONE new, narrowly
    scoped field (`rcsb_entity_host_organism`) was added and fetched ONLY
    for the small set of polymer entities whose source organism does not
    include taxon 9606, to distinguish source organism from expression host
    per this stage's explicit instruction.

Reads:
    config/nr_metadata.py
    data/manifests/nr_metadata_uniprot_audit.csv          (Stage 1: reference sequence lengths)
    data/manifests/rcsb_discovery_by_receptor.csv
    data/manifests/rcsb_search_manifest.csv
    data/raw/api/rcsb/data/polymer_entities_batch_*.json  (already fetched; reused)
    data/raw/api/rcsb/data/entries_batch_*.json           (already fetched; reused)

Writes (corrected):
    data/manifests/rcsb_polymer_entities.csv
    data/manifests/rcsb_candidate_inventory.csv
    data/manifests/rcsb_discovery_anomalies.csv
    data/manifests/rcsb_discovery_summary.json
    data/manifests/rcsb_stage2_provenance.json
    reports/tables/rcsb_discovery_report.md

Writes (new):
    data/manifests/rcsb_polymer_entity_uniprot_mappings.csv
    data/manifests/rcsb_source_taxonomy_audit.csv
    data/manifests/rcsb_multi_uniprot_entity_audit.csv
    data/raw/api/rcsb/data/host_organism_supplement.json  (gitignored; minimal targeted fetch)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402
from scripts.utils.rcsb_data_client import RCSB_DATA_API_URL, fetch_batches  # noqa: E402
from scripts.utils.rcsb_extract import (  # noqa: E402
    classify_multi_uniprot_entity,
    compute_mapping_status,
    compute_method_flags,
    extract_entry_fields,
    extract_polymer_entity_fields,
    extract_uniprot_mappings,
)

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
REPORTS_TABLES_DIR = PROJECT_ROOT / "reports" / "tables"
RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "data"

DISCOVERY_BY_RECEPTOR_PATH = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
SEARCH_MANIFEST_PATH = MANIFESTS_DIR / "rcsb_search_manifest.csv"
UNIPROT_AUDIT_PATH = MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv"

POLYMER_ENTITIES_CSV = MANIFESTS_DIR / "rcsb_polymer_entities.csv"
ENTRIES_CSV = MANIFESTS_DIR / "rcsb_entries.csv"
CANDIDATE_INVENTORY_CSV = MANIFESTS_DIR / "rcsb_candidate_inventory.csv"
ANOMALIES_CSV = MANIFESTS_DIR / "rcsb_discovery_anomalies.csv"
SUMMARY_JSON = MANIFESTS_DIR / "rcsb_discovery_summary.json"
PROVENANCE_JSON = MANIFESTS_DIR / "rcsb_stage2_provenance.json"
REPORT_MD = REPORTS_TABLES_DIR / "rcsb_discovery_report.md"

MAPPINGS_CSV = MANIFESTS_DIR / "rcsb_polymer_entity_uniprot_mappings.csv"
TAXONOMY_AUDIT_CSV = MANIFESTS_DIR / "rcsb_source_taxonomy_audit.csv"
MULTI_UNIPROT_AUDIT_CSV = MANIFESTS_DIR / "rcsb_multi_uniprot_entity_audit.csv"
HOST_SUPPLEMENT_RAW_PATH = RAW_DATA_DIR / "host_organism_supplement.json"

PROJECT_ACCESSIONS = {r["uniprot_id"] for r in NR_METADATA}


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_field(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def load_polymer_entity_records() -> dict[str, dict]:
    records: dict[str, dict] = {}
    for path in sorted(RAW_DATA_DIR.glob("polymer_entities_batch_*.json")):
        wrapper = json.loads(path.read_text())
        if not wrapper.get("success"):
            continue
        entities = ((wrapper.get("response") or {}).get("data") or {}).get("polymer_entities") or []
        for e in entities:
            if e and e.get("rcsb_id"):
                records[e["rcsb_id"]] = e
    return records


def load_entry_records() -> dict[str, dict]:
    records: dict[str, dict] = {}
    for path in sorted(RAW_DATA_DIR.glob("entries_batch_*.json")):
        wrapper = json.loads(path.read_text())
        if not wrapper.get("success"):
            continue
        entries = ((wrapper.get("response") or {}).get("data") or {}).get("entries") or []
        for e in entries:
            if e and e.get("rcsb_id"):
                records[e["rcsb_id"]] = e
    return records


def load_reference_sequence_lengths() -> dict[str, int]:
    """uniprot_id -> UniProt sequence_length, for the project's 48 receptors
    ONLY (from the Stage 1 audit, already retrieved — no new request)."""
    df = pd.read_csv(UNIPROT_AUDIT_PATH)
    lengths = {}
    for _, row in df.iterrows():
        if pd.notna(row.get("sequence_length")) and row.get("resolved_primary_accession"):
            lengths[row["resolved_primary_accession"]] = int(row["sequence_length"])
    return lengths


# --------------------------------------------------------------------------
# Normalized mapping table (Section 3)
# --------------------------------------------------------------------------

def build_mapping_table(
    discovery_df: pd.DataFrame,
    entity_records: dict[str, dict],
    reference_lengths: dict[str, int],
) -> pd.DataFrame:
    query_uniprot_by_entity = (
        discovery_df.sort_values("query_retrieved_at_utc")
        .groupby("polymer_entity_id")["query_uniprot_id"]
        .first()
        .to_dict()
    )

    rows = []
    for entity_id, raw in entity_records.items():
        entity_poly = raw.get("entity_poly") or {}
        entity_length = entity_poly.get("rcsb_sample_sequence_length")
        container = raw.get("rcsb_polymer_entity_container_identifiers") or {}
        pdb_id = container.get("entry_id")
        entity_num = container.get("entity_id")
        query_uniprot = query_uniprot_by_entity.get(entity_id, "")

        mappings = extract_uniprot_mappings(raw, entity_length)
        for m in mappings:
            accession = m["mapped_uniprot_accession"]
            ref_length = reference_lengths.get(accession)
            reference_coverage = (
                round(m["total_aligned_length"] / ref_length, 6)
                if (ref_length and m["total_aligned_length"] is not None)
                else None
            )
            rows.append(
                {
                    "polymer_entity_id": entity_id,
                    "pdb_id": pdb_id,
                    "entity_id": entity_num,
                    "mapped_database_name": m["mapped_database_name"],
                    "mapped_uniprot_accession": accession,
                    "reference_sequence_coverage": reference_coverage,
                    "entity_sequence_coverage": m["entity_sequence_coverage"],
                    "alignment_regions_json": json_field(m["aligned_regions"]),
                    "query_project_uniprot": query_uniprot,
                    "is_project_receptor_accession": accession in PROJECT_ACCESSIONS,
                    "is_query_receptor_accession": accession == query_uniprot,
                }
            )

    columns = [
        "polymer_entity_id", "pdb_id", "entity_id",
        "mapped_database_name", "mapped_uniprot_accession",
        "reference_sequence_coverage", "entity_sequence_coverage",
        "alignment_regions_json",
        "query_project_uniprot", "is_project_receptor_accession", "is_query_receptor_accession",
    ]
    df = pd.DataFrame(rows, columns=columns)

    # Flag duplicate (entity, accession) mapping rows explicitly, rather
    # than silently collapsing them (Section 5 requirement).
    dup_mask = df.duplicated(subset=["polymer_entity_id", "mapped_uniprot_accession"], keep=False)
    df["duplicate_mapping_for_same_accession"] = dup_mask
    return df


# --------------------------------------------------------------------------
# Corrected unique polymer entity table (Section 4)
# --------------------------------------------------------------------------

def build_polymer_entities_table_corrected(
    discovery_df: pd.DataFrame,
    entity_records: dict[str, dict],
    mapping_df: pd.DataFrame,
) -> pd.DataFrame:
    unique_entity_ids = sorted(discovery_df["polymer_entity_id"].unique())
    first_query_by_entity = (
        discovery_df.sort_values("query_retrieved_at_utc")
        .groupby("polymer_entity_id")["query_uniprot_id"]
        .first()
        .to_dict()
    )
    mapping_by_entity = defaultdict(list)
    for _, row in mapping_df.iterrows():
        mapping_by_entity[row["polymer_entity_id"]].append(row)

    rows = []
    for entity_id in unique_entity_ids:
        raw = entity_records.get(entity_id)
        query_uniprot_id = first_query_by_entity.get(entity_id, "")

        if raw is None:
            rows.append(
                {
                    "polymer_entity_id": entity_id,
                    "pdb_id": entity_id.split("_", 1)[0],
                    "entity_id": entity_id.split("_", 1)[1] if "_" in entity_id else "",
                    "entity_description": "", "polymer_type": "", "entity_sequence_length": "",
                    "source_organism_names": json_field([]), "source_taxonomy_ids": json_field([]),
                    "label_asym_ids": json_field([]), "auth_asym_ids": json_field([]),
                    "all_mapped_uniprot_accessions": json_field([]),
                    "query_uniprot_reference_sequence_coverage": "",
                    "query_uniprot_entity_sequence_coverage": "",
                    "uniprot_alignment_regions": json_field([]),
                    "mapping_status": "REVIEW_NEEDED",
                    "mapping_notes": "Entity returned by Search API but never returned by Data API in any batch.",
                }
            )
            continue

        extracted = extract_polymer_entity_fields(raw)
        _, mapping_status, mapping_notes = compute_mapping_status(query_uniprot_id, extracted)

        entity_mappings = mapping_by_entity.get(entity_id, [])
        query_mapping_rows = [m for m in entity_mappings if m["mapped_uniprot_accession"] == query_uniprot_id]
        if len(query_mapping_rows) == 1:
            qm = query_mapping_rows[0]
            query_ref_coverage = qm["reference_sequence_coverage"]
            query_entity_coverage = qm["entity_sequence_coverage"]
            alignment_regions = json.loads(qm["alignment_regions_json"])
        elif len(query_mapping_rows) > 1:
            # Duplicate mapping rows for the query accession — do not
            # silently pick one; leave coverage blank and rely on the
            # normalized mapping table (flagged there) as authoritative.
            query_ref_coverage = ""
            query_entity_coverage = ""
            alignment_regions = []
            mapping_notes = (mapping_notes + " | " if mapping_notes else "") + (
                f"RCSB supplied {len(query_mapping_rows)} duplicate alignment mappings for the "
                f"query accession {query_uniprot_id} — see rcsb_polymer_entity_uniprot_mappings.csv."
            )
            mapping_status = "REVIEW_NEEDED"
        else:
            query_ref_coverage = ""
            query_entity_coverage = ""
            alignment_regions = []

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
                "query_uniprot_reference_sequence_coverage": query_ref_coverage,
                "query_uniprot_entity_sequence_coverage": query_entity_coverage,
                "uniprot_alignment_regions": json_field(alignment_regions),
                "mapping_status": mapping_status,
                "mapping_notes": mapping_notes,
            }
        )

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Source taxonomy audit (Section 6)
# --------------------------------------------------------------------------

def fetch_host_organism_supplement(entity_ids: list[str], refresh: bool) -> dict[str, dict]:
    """Minimal, targeted fetch of rcsb_entity_host_organism for the given
    entity IDs only. Cached in a single wrapper file; reused unless
    --refresh. This is the ONLY new RCSB Data API field/request added by
    Stage 2.1."""
    if not entity_ids:
        return {}

    if not refresh and HOST_SUPPLEMENT_RAW_PATH.exists():
        wrapper = json.loads(HOST_SUPPLEMENT_RAW_PATH.read_text())
        cached_ids = set(wrapper.get("ids_requested", []))
        if cached_ids >= set(entity_ids):
            print(f"  host_organism_supplement: {len(entity_ids)} IDs already covered by cache — reusing.")
            entities = ((wrapper.get("response") or {}).get("data") or {}).get("polymer_entities") or []
            return {e["rcsb_id"]: e for e in entities if e and e.get("rcsb_id")}

    print(f"  host_organism_supplement: fetching host organism for {len(entity_ids)} entities "
          f"lacking taxon 9606 source (minimal targeted request)...")
    results = fetch_batches("host_organism_supplement", sorted(entity_ids), batch_size=len(entity_ids) or 1)
    RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)

    records: dict[str, dict] = {}
    all_ids_requested: list[str] = []
    combined_entities: list[dict] = []
    success = True
    for r in results:
        all_ids_requested.extend(r.ids_requested)
        success = success and r.success
        if r.success and r.raw_json:
            entities = ((r.raw_json.get("data") or {}).get("polymer_entities")) or []
            combined_entities.extend(entities)
            for e in entities:
                if e and e.get("rcsb_id"):
                    records[e["rcsb_id"]] = e

    wrapper = {
        "batch_type": "host_organism_supplement",
        "ids_requested": all_ids_requested,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "success": success,
        "response": {"data": {"polymer_entities": combined_entities}},
    }
    HOST_SUPPLEMENT_RAW_PATH.write_text(json.dumps(wrapper, indent=2, ensure_ascii=False))
    return records


def build_taxonomy_audit(
    inventory_df: pd.DataFrame,
    entity_records: dict[str, dict],
    entities_df: pd.DataFrame,
    refresh: bool,
) -> pd.DataFrame:
    non_human = inventory_df[~inventory_df["source_includes_homo_sapiens_9606"]]
    if len(non_human) == 0:
        return pd.DataFrame(
            columns=[
                "uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id",
                "source_organism_names", "source_taxonomy_ids",
                "host_organism_names", "host_taxonomy_ids",
                "all_mapped_uniprot_accessions", "expected_uniprot_mapping_confirmed",
                "taxonomy_audit_status", "taxonomy_audit_reason",
            ]
        )

    entity_ids_needed = sorted(non_human["polymer_entity_id"].unique())
    host_records = fetch_host_organism_supplement(entity_ids_needed, refresh=refresh)

    entities_by_id = entities_df.set_index("polymer_entity_id")

    rows = []
    for _, row in non_human.iterrows():
        entity_id = row["polymer_entity_id"]
        raw = entity_records.get(entity_id) or {}
        polymer_entity = raw.get("rcsb_polymer_entity") or {}
        entity_description = polymer_entity.get("pdbx_description", "")

        host_raw = (host_records.get(entity_id) or {}).get("rcsb_entity_host_organism") or []
        host_names = sorted({h["scientific_name"] for h in host_raw if h.get("scientific_name")})
        host_taxon_ids = sorted({h["ncbi_taxonomy_id"] for h in host_raw if h.get("ncbi_taxonomy_id") is not None})

        entity_row = entities_by_id.loc[entity_id] if entity_id in entities_by_id.index else None
        mapping_status = entity_row["mapping_status"] if entity_row is not None else ""
        mapping_confirmed = mapping_status in ("CONFIRMED", "CONFIRMED_MULTI_MAPPED")

        source_names = json.loads(row["source_organism_names"]) if isinstance(row["source_organism_names"], str) else []
        source_taxa = json.loads(row["source_taxonomy_ids"]) if isinstance(row["source_taxonomy_ids"], str) else []

        # Conservative, evidence-based classification (never infer "just an
        # expression artifact" without host-organism data actually showing
        # a mismatch between source and host that plausibly explains it).
        if not source_names and not source_taxa:
            status, reason = "REVIEW_NEEDED", "No source organism recorded for this entity at all."
        elif any(t in (9606,) for t in source_taxa):
            status, reason = "EXPLAINED", "Source organism list actually includes taxon 9606 despite the flag (recheck)."
        elif host_taxon_ids and 9606 not in source_taxa and set(source_taxa) == set(host_taxon_ids):
            status, reason = (
                "EXPLAINED",
                f"Source organism and expression host are the same ({source_names}) — this entity's "
                f"deposited source is genuinely non-human, not a human protein misattributed to its "
                f"expression host.",
            )
        elif any(name.lower() in ("synthetic construct",) for name in source_names):
            status, reason = (
                "EXPLAINED",
                "Source organism is 'synthetic construct' (taxid 32630) — a chemically/recombinantly "
                "synthesized peptide or engineered sequence, not an extant-organism source; its "
                "UniProt mapping reflects sequence identity/similarity to the reference accession, "
                "not literal biological provenance from that organism.",
            )
        elif any(name.lower() in ("unclassified", "unidentified") for name in source_names):
            status, reason = (
                "REVIEW_NEEDED",
                "Source organism is 'unclassified'/'unidentified' — consistent with a computationally "
                "reconstructed ancestral sequence (not an extant organism); entity_description should "
                "be manually reviewed to confirm before any downstream use.",
            )
        else:
            status, reason = (
                "REVIEW_NEEDED",
                f"Source organism {source_names} (taxon {source_taxa}) is genuinely non-human and "
                f"does not match a recognized expression-host artifact pattern; requires manual review.",
            )

        rows.append(
            {
                "uniprot_id": row["uniprot_id"],
                "nr_code": row["nr_code"],
                "common_name": row["common_name"],
                "pdb_id": row["pdb_id"],
                "polymer_entity_id": entity_id,
                "entity_description": entity_description,
                "source_organism_names": json_field(source_names),
                "source_taxonomy_ids": json_field(source_taxa),
                "host_organism_names": json_field(host_names),
                "host_taxonomy_ids": json_field(host_taxon_ids),
                "all_mapped_uniprot_accessions": row["all_mapped_uniprot_accessions"],
                "expected_uniprot_mapping_confirmed": mapping_confirmed,
                "taxonomy_audit_status": status,
                "taxonomy_audit_reason": reason,
            }
        )

    columns = [
        "uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id", "entity_description",
        "source_organism_names", "source_taxonomy_ids",
        "host_organism_names", "host_taxonomy_ids",
        "all_mapped_uniprot_accessions", "expected_uniprot_mapping_confirmed",
        "taxonomy_audit_status", "taxonomy_audit_reason",
    ]
    return pd.DataFrame(rows)[columns]


# --------------------------------------------------------------------------
# Multi-UniProt entity audit (Section 7)
# --------------------------------------------------------------------------

def build_multi_uniprot_audit(
    mapping_df: pd.DataFrame,
    entities_df: pd.DataFrame,
    discovery_df: pd.DataFrame,
) -> pd.DataFrame:
    receptor_by_uniprot = {r["uniprot_id"]: r for r in NR_METADATA}
    entities_by_id = entities_df.set_index("polymer_entity_id")

    mapping_counts = mapping_df.groupby("polymer_entity_id")["mapped_uniprot_accession"].nunique()
    multi_entity_ids = sorted(mapping_counts[mapping_counts > 1].index)

    rows = []
    for entity_id in multi_entity_ids:
        entity_mappings = mapping_df[mapping_df["polymer_entity_id"] == entity_id]
        query_uniprot = entity_mappings["query_project_uniprot"].iloc[0]
        receptor = receptor_by_uniprot.get(query_uniprot, {})

        project_mapping_rows = entity_mappings[entity_mappings["mapped_uniprot_accession"] == query_uniprot]
        other_mapping_rows = entity_mappings[entity_mappings["mapped_uniprot_accession"] != query_uniprot]

        project_regions = []
        for _, pm in project_mapping_rows.iterrows():
            project_regions.extend(json.loads(pm["alignment_regions_json"]))

        other_mappings_for_classifier = [
            {
                "mapped_uniprot_accession": om["mapped_uniprot_accession"],
                "aligned_regions": json.loads(om["alignment_regions_json"]),
            }
            for _, om in other_mapping_rows.iterrows()
        ]

        if not project_regions:
            classification, reason = (
                "AMBIGUOUS_REVIEW_NEEDED",
                "No alignment region data available for the project receptor's own mapping.",
            )
        else:
            classification, reason = classify_multi_uniprot_entity(project_regions, other_mappings_for_classifier)

        entity_row = entities_by_id.loc[entity_id] if entity_id in entities_by_id.index else None
        all_accessions = sorted(entity_mappings["mapped_uniprot_accession"].unique())

        project_ref_coverage = project_mapping_rows["reference_sequence_coverage"].iloc[0] if len(project_mapping_rows) else None
        project_entity_coverage = project_mapping_rows["entity_sequence_coverage"].iloc[0] if len(project_mapping_rows) else None

        all_ranges = {
            row["mapped_uniprot_accession"]: json.loads(row["alignment_regions_json"])
            for _, row in entity_mappings.iterrows()
        }

        rows.append(
            {
                "polymer_entity_id": entity_id,
                "pdb_id": entity_mappings["pdb_id"].iloc[0],
                "entity_id": entity_mappings["entity_id"].iloc[0],
                "project_receptor_uniprot": query_uniprot,
                "project_receptor_name": receptor.get("common_name", ""),
                "entity_description": (entity_row["entity_description"] if entity_row is not None else ""),
                "entity_sequence_length": (entity_row["entity_sequence_length"] if entity_row is not None else ""),
                "all_mapped_uniprot_accessions": json_field(all_accessions),
                "project_receptor_reference_coverage": project_ref_coverage,
                "project_receptor_entity_coverage": project_entity_coverage,
                "all_mapping_ranges_json": json_field(all_ranges),
                "source_organism_names": (entity_row["source_organism_names"] if entity_row is not None else json_field([])),
                "source_taxonomy_ids": (entity_row["source_taxonomy_ids"] if entity_row is not None else json_field([])),
                "classification": classification,
                "classification_reason": reason,
            }
        )

    columns = [
        "polymer_entity_id", "pdb_id", "entity_id",
        "project_receptor_uniprot", "project_receptor_name",
        "entity_description", "entity_sequence_length",
        "all_mapped_uniprot_accessions",
        "project_receptor_reference_coverage", "project_receptor_entity_coverage",
        "all_mapping_ranges_json",
        "source_organism_names", "source_taxonomy_ids",
        "classification", "classification_reason",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows)[columns]


# --------------------------------------------------------------------------
# Candidate inventory (corrected, Section 5)
# --------------------------------------------------------------------------

def build_candidate_inventory_corrected(
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

    merged["query_uniprot_reference_sequence_coverage"] = merged["query_uniprot_reference_sequence_coverage"]
    merged["query_uniprot_entity_sequence_coverage"] = merged["query_uniprot_entity_sequence_coverage"]

    columns = [
        "uniprot_id", "nr_code", "common_name", "group",
        "pdb_id", "polymer_entity_id", "entity_id",
        "entity_description", "entity_sequence_length",
        "label_asym_ids", "auth_asym_ids",
        "all_mapped_uniprot_accessions", "mapping_contains_query_uniprot", "mapping_status",
        "query_uniprot_reference_sequence_coverage", "query_uniprot_entity_sequence_coverage",
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
# Anomalies (Section 8 wording correction + coverage sanity check)
# --------------------------------------------------------------------------

def build_anomalies_table_corrected(
    discovery_df: pd.DataFrame,
    entities_df: pd.DataFrame,
    mapping_df: pd.DataFrame,
    search_manifest_df: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    entity_to_receptors = discovery_df.groupby("polymer_entity_id")["uniprot_id"].apply(lambda s: sorted(set(s)))
    for entity_id, receptors in entity_to_receptors.items():
        if len(receptors) > 1:
            rows.append(
                {
                    "anomaly_type": "CROSS_MAPPED_ENTITY_MULTIPLE_RECEPTOR_QUERIES",
                    "uniprot_id": json_field(receptors), "nr_code": "",
                    "pdb_id": entity_id.split("_", 1)[0], "polymer_entity_id": entity_id,
                    "details": f"Polymer entity {entity_id} was returned by {len(receptors)} distinct receptor queries: {receptors}",
                    "review_status": "REVIEW_NEEDED",
                }
            )

    for _, row in entities_df.iterrows():
        mapped = json.loads(row["all_mapped_uniprot_accessions"]) if row["all_mapped_uniprot_accessions"] else []
        if len(mapped) > 1:
            rows.append(
                {
                    "anomaly_type": "ENTITY_MULTI_MAPPED_UNIPROT",
                    "uniprot_id": json_field(mapped), "nr_code": "",
                    "pdb_id": row["pdb_id"], "polymer_entity_id": row["polymer_entity_id"],
                    "details": (
                        f"Entity maps to multiple UniProt accessions: {mapped}. "
                        f"See rcsb_multi_uniprot_entity_audit.csv for the evidence-based classification "
                        f"(fusion/chimera vs. same-protein vs. ambiguous)."
                    ),
                    "review_status": "REVIEW_NEEDED",
                }
            )

    entry_to_receptors = discovery_df.groupby("pdb_id")["uniprot_id"].apply(lambda s: sorted(set(s)))
    for pdb_id, receptors in entry_to_receptors.items():
        if len(receptors) > 1:
            rows.append(
                {
                    "anomaly_type": "PDB_ENTRY_MULTIPLE_PROJECT_RECEPTORS",
                    "uniprot_id": json_field(receptors), "nr_code": "",
                    "pdb_id": pdb_id, "polymer_entity_id": "",
                    "details": (
                        f"PDB entry {pdb_id} contains polymer entities mapped to {len(receptors)} distinct "
                        f"project receptors: {receptors}. This is entry-level co-occurrence of separate "
                        f"polymer entities — NOT the same as one polymer entity mapping to multiple UniProt "
                        f"accessions (see ENTITY_MULTI_MAPPED_UNIPROT rows). No physical "
                        f"interaction/heterodimerization is claimed here beyond co-occurrence in the same "
                        f"deposited entry; that would require examining assembly-level contacts, which this "
                        f"stage does not do."
                    ),
                    "review_status": "REVIEW_NEEDED",
                }
            )

    for _, row in search_manifest_df.iterrows():
        if not row["success"]:
            continue
        count = row["result_count_parsed"]
        if count == 0:
            rows.append(
                {
                    "anomaly_type": "ZERO_CANDIDATE_RECEPTOR",
                    "uniprot_id": row["uniprot_id"], "nr_code": row["nr_code"],
                    "pdb_id": "", "polymer_entity_id": "",
                    "details": (
                        f"{row['common_name']} ({row['nr_code']}): query_status=SUCCEEDED_ZERO_RESULTS — "
                        f"the Search API query succeeded and returned 0 experimental polymer entities. "
                        f"This is NOT an API failure."
                    ),
                    "review_status": "REVIEW_NEEDED",
                }
            )
        elif count in (1, 2):
            rows.append(
                {
                    "anomaly_type": "LOW_CANDIDATE_COUNT_RECEPTOR",
                    "uniprot_id": row["uniprot_id"], "nr_code": row["nr_code"],
                    "pdb_id": "", "polymer_entity_id": "",
                    "details": f"{row['common_name']} ({row['nr_code']}): only {count} experimental polymer entity/entities found.",
                    "review_status": "REVIEW_NEEDED",
                }
            )

    for _, row in search_manifest_df.iterrows():
        if not row["success"]:
            rows.append(
                {
                    "anomaly_type": "QUERY_FAILED",
                    "uniprot_id": row["uniprot_id"], "nr_code": row["nr_code"],
                    "pdb_id": "", "polymer_entity_id": "",
                    "details": f"Search API query FAILED for {row['common_name']} ({row['nr_code']}): {row['error_message']}",
                    "review_status": "REVIEW_NEEDED",
                }
            )

    for _, row in entities_df.iterrows():
        if row["mapping_status"] == "REVIEW_NEEDED":
            rows.append(
                {
                    "anomaly_type": "UNCONFIRMED_UNIPROT_MAPPING",
                    "uniprot_id": "", "nr_code": "",
                    "pdb_id": row["pdb_id"], "polymer_entity_id": row["polymer_entity_id"],
                    "details": row["mapping_notes"],
                    "review_status": "REVIEW_NEEDED",
                }
            )

    coverage_over_1 = mapping_df[
        (mapping_df["reference_sequence_coverage"].notna() & (mapping_df["reference_sequence_coverage"] > 1.0))
        | (mapping_df["entity_sequence_coverage"].notna() & (mapping_df["entity_sequence_coverage"] > 1.0))
    ]
    for _, row in coverage_over_1.iterrows():
        rows.append(
            {
                "anomaly_type": "COVERAGE_EXCEEDS_1",
                "uniprot_id": row["mapped_uniprot_accession"], "nr_code": "",
                "pdb_id": row["pdb_id"], "polymer_entity_id": row["polymer_entity_id"],
                "details": (
                    f"Computed coverage exceeds 1.0 (reference={row['reference_sequence_coverage']}, "
                    f"entity={row['entity_sequence_coverage']}) — likely overlapping aligned regions; "
                    f"not clipped, flagged for manual review."
                ),
                "review_status": "REVIEW_NEEDED",
            }
        )

    columns = ["anomaly_type", "uniprot_id", "nr_code", "pdb_id", "polymer_entity_id", "details", "review_status"]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows)[columns]


# --------------------------------------------------------------------------
# Summary + report (regenerated with corrected data)
# --------------------------------------------------------------------------

def build_summary(
    search_manifest_df: pd.DataFrame,
    discovery_df: pd.DataFrame,
    entities_df: pd.DataFrame,
    inventory_df: pd.DataFrame,
    anomalies_df: pd.DataFrame,
    mapping_df: pd.DataFrame,
    taxonomy_audit_df: pd.DataFrame,
    multi_uniprot_audit_df: pd.DataFrame,
) -> dict:
    counts_by_receptor = (
        discovery_df.groupby(["uniprot_id", "nr_code", "common_name"])
        .agg(polymer_entity_count=("polymer_entity_id", "nunique"), unique_pdb_entry_count=("pdb_id", "nunique"))
        .reset_index()
    )
    zero_result_receptors = search_manifest_df[
        (search_manifest_df["success"]) & (search_manifest_df["result_count_parsed"] == 0)
    ]["common_name"].tolist()
    counts_1_2 = counts_by_receptor[counts_by_receptor["polymer_entity_count"].isin([1, 2])]["common_name"].tolist()
    counts_ge3 = counts_by_receptor[counts_by_receptor["polymer_entity_count"] >= 3]["common_name"].tolist()

    method_counts = {
        "xray": int(inventory_df["is_xray"].sum()),
        "cryo_em": int(inventory_df["is_cryo_em"].sum()),
        "nmr": int(inventory_df["is_nmr"].sum()),
        "other_experimental": int(inventory_df["is_other_experimental"].sum()),
    }

    ref_cov = pd.to_numeric(mapping_df["reference_sequence_coverage"], errors="coerce")
    ent_cov = pd.to_numeric(mapping_df["entity_sequence_coverage"], errors="coerce")

    return {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "Stage 2.1 correction/audit",
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
        "candidates_uniprot_mapping_confirmed": int(entities_df["mapping_status"].isin(["CONFIRMED", "CONFIRMED_MULTI_MAPPED"]).sum()),
        "candidates_uniprot_mapping_review_needed": int((entities_df["mapping_status"] == "REVIEW_NEEDED").sum()),
        "number_of_cross_mapped_entities": int((anomalies_df["anomaly_type"] == "CROSS_MAPPED_ENTITY_MULTIPLE_RECEPTOR_QUERIES").sum()),
        "number_of_pdb_entries_with_multiple_project_receptors": int((anomalies_df["anomaly_type"] == "PDB_ENTRY_MULTIPLE_PROJECT_RECEPTORS").sum()),
        "number_of_anomalies_total": len(anomalies_df),
        "coverage_correction": {
            "total_mapping_rows": len(mapping_df),
            "unique_polymer_entities_in_mapping_table": int(mapping_df["polymer_entity_id"].nunique()),
            "entities_with_exactly_one_mapping": int((mapping_df.groupby("polymer_entity_id").size() == 1).sum()),
            "entities_with_multiple_mappings": int((mapping_df.groupby("polymer_entity_id").size() > 1).sum()),
            "reference_coverage_available": int(ref_cov.notna().sum()),
            "reference_coverage_missing": int(ref_cov.isna().sum()),
            "reference_coverage_min": float(ref_cov.min()) if ref_cov.notna().any() else None,
            "reference_coverage_max": float(ref_cov.max()) if ref_cov.notna().any() else None,
            "entity_coverage_available": int(ent_cov.notna().sum()),
            "entity_coverage_missing": int(ent_cov.isna().sum()),
            "entity_coverage_min": float(ent_cov.min()) if ent_cov.notna().any() else None,
            "entity_coverage_max": float(ent_cov.max()) if ent_cov.notna().any() else None,
        },
        "taxonomy_audit": {
            "candidates_lacking_taxon_9606": len(taxonomy_audit_df),
            "explained": int((taxonomy_audit_df["taxonomy_audit_status"] == "EXPLAINED").sum()) if len(taxonomy_audit_df) else 0,
            "review_needed": int((taxonomy_audit_df["taxonomy_audit_status"] == "REVIEW_NEEDED").sum()) if len(taxonomy_audit_df) else 0,
        },
        "multi_uniprot_audit": {
            "total_multi_mapped_entities": len(multi_uniprot_audit_df),
            "by_classification": multi_uniprot_audit_df["classification"].value_counts().to_dict() if len(multi_uniprot_audit_df) else {},
        },
        "note": "This is a discovery-stage inventory, corrected/audited at Stage 2.1. No final structural QC filters were applied.",
    }


def write_report(
    search_manifest_df: pd.DataFrame,
    counts_by_receptor: pd.DataFrame,
    summary: dict,
    anomalies_df: pd.DataFrame,
    taxonomy_audit_df: pd.DataFrame,
    multi_uniprot_audit_df: pd.DataFrame,
) -> None:
    lines = ["# RCSB Candidate Structure Discovery — Stage 2 (corrected at Stage 2.1)\n"]
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
        "polymer entities and their parent PDB entries were enriched via "
        "the official RCSB Data API (GraphQL, batched). GraphQL "
        "introspection against `RcsbPolymerEntityAlign` confirmed the Data "
        "API exposes no precomputed coverage-fraction field — only "
        "`aligned_regions` (entity_beg_seq_id, ref_beg_seq_id, length). "
        "`reference_sequence_coverage` and `entity_sequence_coverage` are "
        "therefore computed deterministically from these raw region "
        "lengths and the already-retrieved entity/reference sequence "
        "lengths (Stage 1 UniProt data for the project's 48 accessions); "
        "reference coverage is left null for non-project (e.g. fusion "
        "partner) accessions rather than guessed. No mmCIF or PDB "
        "coordinate files were downloaded at any point.\n"
    )
    lines.append("## Overall counts\n")
    for key in (
        "number_of_validated_receptors_queried", "number_of_successful_receptor_queries",
        "number_of_failed_receptor_queries", "total_receptor_polymer_entity_relationships",
        "unique_polymer_entities", "unique_pdb_entries", "receptors_with_zero_candidates",
        "receptors_with_1_to_2_candidates", "receptors_with_3_or_more_candidates",
    ):
        lines.append(f"- {key}: {summary[key]}")
    lines.append("")

    lines.append("## All 48 receptors\n")
    lines.append("| common_name | nr_code | uniprot_id | polymer_entity_count | unique_pdb_entry_count | query_status |")
    lines.append("|---|---|---|---|---|---|")
    merged = search_manifest_df.merge(counts_by_receptor, on=["uniprot_id", "nr_code", "common_name"], how="left")
    merged["polymer_entity_count"] = merged["polymer_entity_count"].fillna(0).astype(int)
    merged["unique_pdb_entry_count"] = merged["unique_pdb_entry_count"].fillna(0).astype(int)
    for _, row in merged.sort_values("record_index").iterrows():
        status = "SUCCEEDED_ZERO_RESULTS" if (row["success"] and row["polymer_entity_count"] == 0) else ("FAILED" if not row["success"] else "SUCCEEDED")
        lines.append(
            f"| {row['common_name']} | {row['nr_code']} | {row['uniprot_id']} | "
            f"{row['polymer_entity_count']} | {row['unique_pdb_entry_count']} | {status} |"
        )
    lines.append("")

    lines.append("## Receptors with zero candidates\n")
    zero_list = summary["receptors_with_zero_candidates_list"]
    lines.append(
        (", ".join(zero_list) if zero_list else "None.")
        + " — all query_status=SUCCEEDED_ZERO_RESULTS (Search API succeeded, 0 matches), not API failures.\n"
    )

    lines.append("## Receptors with 1-2 candidates\n")
    low_list = summary["receptors_with_1_to_2_candidates_list"]
    lines.append((", ".join(low_list) if low_list else "None.") + "\n")

    lines.append("## Experimental-method breakdown (candidate rows)\n")
    mc = summary["experimental_method_counts_across_candidate_rows"]
    lines.append(f"- X-ray: {mc['xray']}")
    lines.append(f"- cryo-EM: {mc['cryo_em']}")
    lines.append(f"- NMR: {mc['nmr']}")
    lines.append(f"- Other experimental: {mc['other_experimental']}")
    lines.append("")

    lines.append("## Source taxonomy audit (candidates lacking taxon 9606)\n")
    lines.append(f"{len(taxonomy_audit_df)} candidates audited. "
                  f"EXPLAINED: {summary['taxonomy_audit']['explained']}, "
                  f"REVIEW_NEEDED: {summary['taxonomy_audit']['review_needed']}.\n")
    if len(taxonomy_audit_df):
        lines.append("| common_name | pdb_id | entity | source | host | status |")
        lines.append("|---|---|---|---|---|---|")
        for _, row in taxonomy_audit_df.iterrows():
            lines.append(
                f"| {row['common_name']} | {row['pdb_id']} | {row['polymer_entity_id']} | "
                f"{row['source_organism_names']} | {row['host_organism_names']} | {row['taxonomy_audit_status']} |"
            )
        lines.append("")

    lines.append("## Multi-UniProt entity audit\n")
    lines.append(f"{len(multi_uniprot_audit_df)} multi-mapped entities. By classification: "
                  f"{summary['multi_uniprot_audit']['by_classification']}\n")
    ambiguous = multi_uniprot_audit_df[multi_uniprot_audit_df["classification"] == "AMBIGUOUS_REVIEW_NEEDED"] if len(multi_uniprot_audit_df) else multi_uniprot_audit_df
    lines.append("### AMBIGUOUS_REVIEW_NEEDED items\n")
    if len(ambiguous) == 0:
        lines.append("None.\n")
    else:
        lines.append("| polymer_entity_id | project_receptor_name | all_mapped_uniprot_accessions | reason |")
        lines.append("|---|---|---|---|")
        for _, row in ambiguous.iterrows():
            lines.append(
                f"| {row['polymer_entity_id']} | {row['project_receptor_name']} | "
                f"{row['all_mapped_uniprot_accessions']} | {row['classification_reason']} |"
            )
        lines.append("")

    lines.append("## Cross-receptor PDB entries\n")
    cross_entries = anomalies_df[anomalies_df["anomaly_type"] == "PDB_ENTRY_MULTIPLE_PROJECT_RECEPTORS"]
    lines.append(
        f"{len(cross_entries)} PDB entries contain polymer entities mapped to more than one project "
        f"receptor. This is ENTRY-LEVEL co-occurrence of distinct polymer entities, which is different "
        f"from one POLYMER ENTITY mapping to multiple UniProt accessions (see "
        f"`ENTITY_MULTI_MAPPED_UNIPROT` / the multi-UniProt audit above). No physical "
        f"interaction/heterodimerization is claimed unless assembly-level contacts have been examined, "
        f"which this stage does not do.\n"
    )
    if len(cross_entries) > 0:
        lines.append("| pdb_id | receptors (uniprot_id) |")
        lines.append("|---|---|")
        for _, row in cross_entries.iterrows():
            lines.append(f"| {row['pdb_id']} | {row['uniprot_id']} |")
        lines.append("")

    lines.append("## Important statements\n")
    lines.append("No final structural QC filters were applied at Stage 2 or Stage 2.1.\n")
    lines.append("No mmCIF coordinate files were downloaded at Stage 2 or Stage 2.1.\n")

    REPORT_MD.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-host-organism", action="store_true",
                         help="Force re-fetching the host-organism supplement even if cached.")
    args = parser.parse_args()

    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR_PATH, dtype={"pdb_id": str, "entity_id": str})
    search_manifest_df = pd.read_csv(SEARCH_MANIFEST_PATH)
    reference_lengths = load_reference_sequence_lengths()

    entity_records = load_polymer_entity_records()
    entry_records = load_entry_records()
    print(f"Loaded {len(entity_records)} polymer entity records, {len(entry_records)} entry records "
          f"(all reused from existing raw batches — no re-fetch needed for coverage correction).")
    print(f"Loaded reference sequence lengths for {len(reference_lengths)} project accessions from Stage 1.")

    print("\nBuilding normalized UniProt mapping table...")
    mapping_df = build_mapping_table(discovery_df, entity_records, reference_lengths)
    mapping_df.to_csv(MAPPINGS_CSV, index=False)
    print(f"Wrote {MAPPINGS_CSV.relative_to(PROJECT_ROOT)} ({len(mapping_df)} rows)")

    print("\nRebuilding corrected unique polymer entity table...")
    entities_df = build_polymer_entities_table_corrected(discovery_df, entity_records, mapping_df)
    entities_df.to_csv(POLYMER_ENTITIES_CSV, index=False)
    print(f"Wrote {POLYMER_ENTITIES_CSV.relative_to(PROJECT_ROOT)} ({len(entities_df)} rows)")

    entries_df = pd.read_csv(ENTRIES_CSV, dtype={"pdb_id": str})

    print("\nRebuilding corrected candidate inventory...")
    inventory_df = build_candidate_inventory_corrected(discovery_df, entities_df, entries_df)
    inventory_df.to_csv(CANDIDATE_INVENTORY_CSV, index=False)
    print(f"Wrote {CANDIDATE_INVENTORY_CSV.relative_to(PROJECT_ROOT)} ({len(inventory_df)} rows)")

    print("\nBuilding source taxonomy audit (fetching host organism for non-9606-source entities)...")
    taxonomy_audit_df = build_taxonomy_audit(inventory_df, entity_records, entities_df, refresh=args.refresh_host_organism)
    taxonomy_audit_df.to_csv(TAXONOMY_AUDIT_CSV, index=False)
    print(f"Wrote {TAXONOMY_AUDIT_CSV.relative_to(PROJECT_ROOT)} ({len(taxonomy_audit_df)} rows)")

    print("\nBuilding multi-UniProt entity audit...")
    multi_uniprot_audit_df = build_multi_uniprot_audit(mapping_df, entities_df, discovery_df)
    multi_uniprot_audit_df.to_csv(MULTI_UNIPROT_AUDIT_CSV, index=False)
    print(f"Wrote {MULTI_UNIPROT_AUDIT_CSV.relative_to(PROJECT_ROOT)} ({len(multi_uniprot_audit_df)} rows)")

    print("\nRebuilding anomalies table (corrected wording + coverage sanity check)...")
    anomalies_df = build_anomalies_table_corrected(discovery_df, entities_df, mapping_df, search_manifest_df)
    anomalies_df.to_csv(ANOMALIES_CSV, index=False)
    print(f"Wrote {ANOMALIES_CSV.relative_to(PROJECT_ROOT)} ({len(anomalies_df)} rows)")

    summary = build_summary(
        search_manifest_df, discovery_df, entities_df, inventory_df, anomalies_df,
        mapping_df, taxonomy_audit_df, multi_uniprot_audit_df,
    )
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"Wrote {SUMMARY_JSON.relative_to(PROJECT_ROOT)}")

    counts_by_receptor = (
        discovery_df.groupby(["uniprot_id", "nr_code", "common_name"])
        .agg(polymer_entity_count=("polymer_entity_id", "nunique"), unique_pdb_entry_count=("pdb_id", "nunique"))
        .reset_index()
    )
    write_report(search_manifest_df, counts_by_receptor, summary, anomalies_df, taxonomy_audit_df, multi_uniprot_audit_df)
    print(f"Wrote {REPORT_MD.relative_to(PROJECT_ROOT)}")

    provenance = {
        "project_name": "Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families",
        "description": (
            "Stage 2.1: correction of the reference/entity sequence coverage representation, plus a "
            "source-taxonomy and multi-UniProt-mapping audit of the Stage 2 RCSB discovery inventory. "
            "No final structural QC filtering, no coordinate downloads."
        ),
        "correction_summary": (
            "rcsb_polymer_entities.csv previously stored identical raw alignment-region JSON in both "
            "reference_sequence_coverage and entity_sequence_coverage columns. GraphQL introspection "
            "confirmed RCSB exposes no precomputed coverage field; both are now computed deterministically "
            "from aligned_regions (entity_beg_seq_id/ref_beg_seq_id/length) and already-retrieved sequence "
            "lengths. One new field (rcsb_entity_host_organism) was fetched, narrowly, only for entities "
            "lacking taxon 9606 source."
        ),
        "validated_receptor_master_sha256": hashlib.sha256((PROJECT_ROOT / "config" / "nr_metadata.py").read_bytes()).hexdigest(),
        "rcsb_search_api_endpoint": "https://search.rcsb.org/rcsbsearch/v2/query",
        "rcsb_data_api_endpoint": RCSB_DATA_API_URL,
        "new_data_api_fields_added_at_stage_2_1": ["rcsb_entity_host_organism (targeted fetch only)"],
        "number_of_search_api_queries": len(search_manifest_df),
        "number_of_data_api_batches_reused_from_stage2": len(list(RAW_DATA_DIR.glob("polymer_entities_batch_*.json"))) + len(list(RAW_DATA_DIR.glob("entries_batch_*.json"))),
        "host_organism_supplement_fetched_for_n_entities": int(len(pd.read_csv(TAXONOMY_AUDIT_CSV))) if TAXONOMY_AUDIT_CSV.exists() else 0,
        "generated_manifests": {
            "rcsb_search_manifest.csv": sha256_of_file(SEARCH_MANIFEST_PATH),
            "rcsb_discovery_by_receptor.csv": sha256_of_file(DISCOVERY_BY_RECEPTOR_PATH),
            "rcsb_polymer_entities.csv": sha256_of_file(POLYMER_ENTITIES_CSV),
            "rcsb_entries.csv": sha256_of_file(ENTRIES_CSV),
            "rcsb_candidate_inventory.csv": sha256_of_file(CANDIDATE_INVENTORY_CSV),
            "rcsb_discovery_anomalies.csv": sha256_of_file(ANOMALIES_CSV),
            "rcsb_polymer_entity_uniprot_mappings.csv": sha256_of_file(MAPPINGS_CSV),
            "rcsb_source_taxonomy_audit.csv": sha256_of_file(TAXONOMY_AUDIT_CSV),
            "rcsb_multi_uniprot_entity_audit.csv": sha256_of_file(MULTI_UNIPROT_AUDIT_CSV),
        },
        "python_version": sys.version,
        "generating_scripts": [
            "scripts/04_discover_rcsb_candidates.py",
            "scripts/05_fetch_rcsb_metadata.py",
            "scripts/06_build_rcsb_inventory.py",
            "scripts/07_correct_rcsb_coverage_and_audit.py",
        ],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    PROVENANCE_JSON.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Wrote {PROVENANCE_JSON.relative_to(PROJECT_ROOT)}")

    print("\n--- Stage 2.1 correction/audit summary ---")
    print(f"mapping rows: {len(mapping_df)}  unique entities in mapping table: {mapping_df['polymer_entity_id'].nunique()}")
    print(f"taxonomy audit rows: {len(taxonomy_audit_df)}  multi-uniprot audit rows: {len(multi_uniprot_audit_df)}")
    print(f"anomalies total: {len(anomalies_df)}")


if __name__ == "__main__":
    main()
