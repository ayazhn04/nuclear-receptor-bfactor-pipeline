"""
Stage 1, step 2: validate every professor-supplied UniProt accession against
the CURRENT official UniProt REST API (https://rest.uniprot.org/).

Reads:
    config/nr_metadata_professor.py

Writes:
    data/raw/api/uniprot/{accession}.json         (raw responses, gitignored)
    data/manifests/uniprot_api_manifest.csv        (fetch audit trail)
    data/manifests/nr_metadata_uniprot_audit.csv   (main audit table)
    data/manifests/nr_metadata_uniprot_audit_summary.json
    reports/tables/nr_metadata_uniprot_audit.md

Idempotent by default: if a valid cached raw response already exists for an
accession, it is reused and NOT re-fetched. Pass --refresh to force
re-fetching every accession from the live API.

This script NEVER modifies config/nr_metadata_professor.py or
config/nr_metadata.py, and never overwrites a professor-supplied value with
a UniProt-resolved value.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata_professor import NR_METADATA  # noqa: E402
from scripts.utils.uniprot_client import UNIPROT_API_BASE, UniProtFetchResult, fetch_uniprot_entry  # noqa: E402
from scripts.utils.uniprot_extract import check_identity, decide_audit_status, extract_fields  # noqa: E402

RAW_UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"
MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
REPORTS_TABLES_DIR = PROJECT_ROOT / "reports" / "tables"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def utc_mtime_iso(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def load_cached_response(accession: str) -> dict | None:
    raw_path = RAW_UNIPROT_DIR / f"{accession}.json"
    if not raw_path.exists():
        return None
    try:
        return json.loads(raw_path.read_text())
    except json.JSONDecodeError:
        return None


def resolve_one(accession: str, refresh: bool) -> dict:
    """Fetch (or reuse cached) UniProt data for one accession. Returns a
    manifest-row dict; always produced, success or failure."""
    raw_path = RAW_UNIPROT_DIR / f"{accession}.json"

    if not refresh:
        cached = load_cached_response(accession)
        if cached is not None:
            raw_bytes = raw_path.read_bytes()
            return {
                "requested_accession": accession,
                "resolved_primary_accession": cached.get("primaryAccession"),
                "request_url": f"{UNIPROT_API_BASE}/{accession}.json",
                "http_status": 200,
                "retrieved_at_utc": utc_mtime_iso(raw_path),
                "raw_file": str(raw_path.relative_to(PROJECT_ROOT)),
                "sha256": sha256_bytes(raw_bytes),
                "success": True,
                "error_message": "",
                "_raw_json": cached,
                "_from_cache": True,
            }

    result: UniProtFetchResult = fetch_uniprot_entry(accession)

    if result.success and result.raw_json is not None:
        RAW_UNIPROT_DIR.mkdir(parents=True, exist_ok=True)
        raw_bytes = json.dumps(result.raw_json, indent=2, ensure_ascii=False).encode("utf-8")
        raw_path.write_bytes(raw_bytes)
        return {
            "requested_accession": accession,
            "resolved_primary_accession": result.raw_json.get("primaryAccession"),
            "request_url": result.request_url,
            "http_status": result.http_status,
            "retrieved_at_utc": utc_mtime_iso(raw_path),
            "raw_file": str(raw_path.relative_to(PROJECT_ROOT)),
            "sha256": sha256_bytes(raw_bytes),
            "success": True,
            "error_message": "",
            "_raw_json": result.raw_json,
            "_from_cache": False,
        }

    return {
        "requested_accession": accession,
        "resolved_primary_accession": "",
        "request_url": result.request_url,
        "http_status": result.http_status if result.http_status is not None else "",
        "retrieved_at_utc": result.retrieved_at_utc,
        "raw_file": "",
        "sha256": "",
        "success": False,
        "error_message": result.error_message,
        "_raw_json": None,
        "_from_cache": False,
    }


def build_audit_row(record_index: int, professor_record: dict, manifest_row: dict) -> dict:
    professor_row = {
        "record_index": record_index,
        "professor_uniprot_id": professor_record["uniprot_id"],
        "professor_nr_code": professor_record["nr_code"],
        "professor_common_name": professor_record["common_name"],
        "professor_group": professor_record["group"],
        "professor_search_terms": json.dumps(
            professor_record["search_terms"], ensure_ascii=False, separators=(",", ":")
        ),
        "requested_accession": manifest_row["requested_accession"],
        "http_status": manifest_row["http_status"],
        "retrieved_at_utc": manifest_row["retrieved_at_utc"],
    }

    if manifest_row["success"] and manifest_row["_raw_json"] is not None:
        extracted = extract_fields(manifest_row["_raw_json"])
        checks = check_identity(professor_record, extracted)
        audit_status, audit_reason = decide_audit_status(True, manifest_row["http_status"], checks)
        professor_row.update(extracted)
        professor_row.update(checks)
    else:
        # Failure: emit an audit row with empty/False fields rather than
        # dropping the record. Never fabricate resolved data.
        empty_extracted = {
            "resolved_primary_accession": "",
            "uniprot_entry_name": "",
            "entry_type": "",
            "reviewed_status": "",
            "organism_scientific_name": "",
            "organism_common_name": "",
            "organism_taxon_id": "",
            "recommended_protein_name": "",
            "primary_gene_name": "",
            "gene_synonyms": "",
            "sequence_length": "",
            "sequence_version": "",
            "entry_version": "",
            "first_public_date": "",
            "last_sequence_update_date": "",
            "last_annotation_update_date": "",
        }
        empty_checks = {
            "accession_resolved": False,
            "primary_accession_matches_professor": False,
            "organism_is_homo_sapiens": False,
            "organism_taxon_is_9606": False,
            "entry_is_reviewed": False,
            "gene_name_available": False,
            "protein_name_available": False,
            "sequence_available": False,
            "sequence_length_positive": False,
            "naming_overlap_found": False,
        }
        audit_status, audit_reason = decide_audit_status(False, manifest_row["http_status"], None)
        audit_reason = f"{audit_reason} error_message={manifest_row['error_message']!r}"
        professor_row.update(empty_extracted)
        professor_row.update(empty_checks)

    professor_row["audit_status"] = audit_status
    professor_row["audit_reason"] = audit_reason
    professor_row["nomenclature_check"] = "NOT_AUTOMATED"
    return professor_row


AUDIT_COLUMNS = [
    "record_index",
    "professor_uniprot_id", "professor_nr_code", "professor_common_name",
    "professor_group", "professor_search_terms",
    "requested_accession", "resolved_primary_accession",
    "uniprot_entry_name", "reviewed_status",
    "organism_scientific_name", "organism_taxon_id",
    "recommended_protein_name", "primary_gene_name", "gene_synonyms",
    "sequence_length", "sequence_version", "entry_version",
    "first_public_date", "last_sequence_update_date", "last_annotation_update_date",
    "http_status", "retrieved_at_utc",
    "primary_accession_matches_professor", "organism_is_homo_sapiens",
    "organism_taxon_is_9606", "entry_is_reviewed",
    "naming_overlap_found", "nomenclature_check",
    "audit_status", "audit_reason",
]

MANIFEST_COLUMNS = [
    "requested_accession", "resolved_primary_accession", "request_url",
    "http_status", "retrieved_at_utc", "raw_file", "sha256", "success", "error_message",
]


def write_markdown_report(audit_df: pd.DataFrame, summary: dict) -> Path:
    lines = []
    lines.append("# Nuclear Receptor Metadata — Current UniProt Audit\n")
    lines.append(
        "## Methods\n\n"
        "Every UniProt accession in the professor-supplied `NR_METADATA` "
        "snapshot (`config/nr_metadata_professor.py`) was queried against "
        f"the official UniProt REST API (`{summary['uniprot_api_base_url']}`) "
        "and compared against the professor-supplied values. Comparisons "
        "were conservative: name differences alone were never treated as "
        "proof of a mismatch, and no professor-supplied value was altered. "
        "See `data/manifests/nr_metadata_uniprot_audit.csv` for the full "
        "per-record table and `data/manifests/uniprot_api_manifest.csv` for "
        "the raw fetch audit trail.\n"
    )
    lines.append("## Overall counts\n")
    lines.append(f"- Professor records: {summary['professor_record_count']}")
    lines.append(f"- Successful API responses: {summary['successful_responses']}")
    lines.append(f"- Failed API responses: {summary['failed_responses']}")
    lines.append(f"- PASS: {summary['pass_count']}")
    lines.append(f"- REVIEW_NEEDED: {summary['review_needed_count']}")
    lines.append(f"- FAIL: {summary['fail_count']}")
    lines.append(f"- Human / taxon 9606: {summary['human_taxon_9606_count']}")
    lines.append(f"- Reviewed (Swiss-Prot) entries: {summary['reviewed_count']}")
    lines.append(
        f"- Primary accession differs from professor accession: "
        f"{summary['primary_accession_differs_count']}\n"
    )

    lines.append("## PASS / REVIEW_NEEDED / FAIL summary\n")
    lines.append("| Status | Count |")
    lines.append("|---|---|")
    lines.append(f"| PASS | {summary['pass_count']} |")
    lines.append(f"| REVIEW_NEEDED | {summary['review_needed_count']} |")
    lines.append(f"| FAIL | {summary['fail_count']} |\n")

    non_pass = audit_df[audit_df["audit_status"] != "PASS"]
    lines.append("## Non-PASS records\n")
    if len(non_pass) == 0:
        lines.append("None.\n")
    else:
        lines.append(
            "| record_index | professor_common_name | professor_nr_code | "
            "professor_uniprot_id | resolved_primary_accession | audit_status | audit_reason |"
        )
        lines.append("|---|---|---|---|---|---|---|")
        for _, row in non_pass.iterrows():
            reason = str(row["audit_reason"]).replace("|", "\\|")
            lines.append(
                f"| {row['record_index']} | {row['professor_common_name']} | "
                f"{row['professor_nr_code']} | {row['professor_uniprot_id']} | "
                f"{row['resolved_primary_accession']} | {row['audit_status']} | {reason} |"
            )
        lines.append("")

    lines.append("## Unexpected findings\n")
    accession_changes = audit_df[
        (audit_df["resolved_primary_accession"] != "")
        & (audit_df["resolved_primary_accession"] != audit_df["professor_uniprot_id"])
    ]
    if len(accession_changes) == 0:
        lines.append("No accessions resolved to a different current primary accession.\n")
    else:
        lines.append(
            f"{len(accession_changes)} accession(s) resolved to a different current "
            "primary accession than supplied by the professor. See the accession "
            "change table below; professor values were preserved unchanged.\n"
        )
        lines.append("| professor_common_name | professor_uniprot_id | resolved_primary_accession |")
        lines.append("|---|---|---|")
        for _, row in accession_changes.iterrows():
            lines.append(
                f"| {row['professor_common_name']} | {row['professor_uniprot_id']} | "
                f"{row['resolved_primary_accession']} |"
            )
        lines.append("")

    lines.append("## Data integrity statement\n")
    lines.append(
        "No professor-supplied value (`uniprot_id`, `nr_code`, `common_name`, "
        "`search_terms`, `group`) was modified, overwritten, or silently "
        "corrected during this audit. `config/nr_metadata_professor.py` and "
        "`config/nr_metadata.py` (still `NR_METADATA = []`) were not touched "
        "by this script beyond being read. All discrepancies are recorded "
        "side-by-side in the audit table for manual review.\n"
    )

    out_path = REPORTS_TABLES_DIR / "nr_metadata_uniprot_audit.md"
    out_path.write_text("\n".join(lines))
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh", action="store_true",
        help="Force re-fetching every accession from the live API, ignoring any cached raw response.",
    )
    args = parser.parse_args()

    RAW_UNIPROT_DIR.mkdir(parents=True, exist_ok=True)
    MANIFESTS_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_TABLES_DIR.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    audit_rows = []

    print(f"Validating {len(NR_METADATA)} professor-supplied accessions against {UNIPROT_API_BASE} ...")
    for idx, record in enumerate(NR_METADATA):
        accession = record["uniprot_id"]
        manifest_row = resolve_one(accession, refresh=args.refresh)
        source = "cache" if manifest_row.get("_from_cache") else ("live" if manifest_row["success"] else "FAILED")
        print(f"  [{idx+1:2d}/{len(NR_METADATA)}] {accession:10s} common_name={record['common_name']:10s} -> {source}")

        manifest_rows.append({k: v for k, v in manifest_row.items() if not k.startswith("_")})
        audit_rows.append(build_audit_row(idx, record, manifest_row))

    manifest_df = pd.DataFrame(manifest_rows)[MANIFEST_COLUMNS]
    manifest_path = MANIFESTS_DIR / "uniprot_api_manifest.csv"
    manifest_df.to_csv(manifest_path, index=False)
    print(f"\nWrote {manifest_path.relative_to(PROJECT_ROOT)}")

    audit_df = pd.DataFrame(audit_rows)[AUDIT_COLUMNS]
    audit_path = MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv"
    audit_df.to_csv(audit_path, index=False)
    print(f"Wrote {audit_path.relative_to(PROJECT_ROOT)}")

    status_counts = Counter(audit_df["audit_status"])
    accession_diff_mask = (
        (audit_df["resolved_primary_accession"] != "")
        & (audit_df["resolved_primary_accession"] != audit_df["professor_uniprot_id"])
    )

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "uniprot_api_base_url": UNIPROT_API_BASE,
        "professor_record_count": len(NR_METADATA),
        "successful_responses": int(manifest_df["success"].sum()),
        "failed_responses": int((~manifest_df["success"]).sum()),
        "pass_count": int(status_counts.get("PASS", 0)),
        "review_needed_count": int(status_counts.get("REVIEW_NEEDED", 0)),
        "fail_count": int(status_counts.get("FAIL", 0)),
        "human_taxon_9606_count": int(audit_df["organism_taxon_is_9606"].sum()),
        "reviewed_count": int(audit_df["entry_is_reviewed"].sum()),
        "primary_accession_differs_count": int(accession_diff_mask.sum()),
        "review_needed_records": audit_df.loc[
            audit_df["audit_status"] == "REVIEW_NEEDED", "professor_common_name"
        ].tolist(),
        "fail_records": audit_df.loc[
            audit_df["audit_status"] == "FAIL", "professor_common_name"
        ].tolist(),
        "group_counts": dict(Counter(r["group"] for r in NR_METADATA)),
    }
    summary_path = MANIFESTS_DIR / "nr_metadata_uniprot_audit_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {summary_path.relative_to(PROJECT_ROOT)}")

    report_path = write_markdown_report(audit_df, summary)
    print(f"Wrote {report_path.relative_to(PROJECT_ROOT)}")

    print("\n--- Audit summary ---")
    for key in (
        "professor_record_count", "successful_responses", "failed_responses",
        "pass_count", "review_needed_count", "fail_count",
        "human_taxon_9606_count", "reviewed_count", "primary_accession_differs_count",
    ):
        print(f"{key}: {summary[key]}")


if __name__ == "__main__":
    main()
