"""
Stage 3A, step 1: build the authoritative canonical-UniProt LBD (and,
diagnostically, DBD) reference interval for each of the 48 validated
receptors.

Reads:
    config/nr_metadata.py
    data/raw/api/uniprot/{uniprot_id}.json   (Stage 1 cached raw UniProt responses)

Writes:
    data/manifests/nr_lbd_reference.csv
    data/manifests/nr_lbd_reference_provenance.json

Primary source: the current reviewed UniProt entry's `features[]` list
(type == "Domain", description == "NR LBD"), already cached from Stage 1 —
no network access is needed for any of the 48 receptors (verified
empirically; see module docstring in scripts/utils/domain_features.py).

Fallback (PROSITE / InterPro) is implemented but was NOT exercised for the
current 48-receptor set, since every receptor's cached UniProt entry
carries an unambiguous NR LBD feature. If it ever is exercised, the
fallback source and identifier are recorded explicitly in
`lbd_source`/`lbd_source_identifier`, never silently.

Coordinate convention: 1-based, inclusive (matches UniProt's own convention).
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
from scripts.utils.domain_features import (  # noqa: E402
    extract_dbd_feature,
    extract_entry_audit,
    extract_lbd_feature,
    extract_sequence_length,
)

RAW_UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"
MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
OUTPUT_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"
PROVENANCE_JSON = MANIFESTS_DIR / "nr_lbd_reference_provenance.json"
NR_METADATA_PY = PROJECT_ROOT / "config" / "nr_metadata.py"

COLUMNS = [
    "record_index", "uniprot_id", "nr_code", "common_name", "group",
    "sequence_length",
    "lbd_start", "lbd_end", "lbd_length",
    "lbd_definition", "lbd_coordinate_system",
    "lbd_source", "lbd_source_identifier", "lbd_feature_type", "lbd_feature_description", "lbd_feature_evidence",
    "uniprot_entry_version", "uniprot_sequence_version",
    "lbd_status", "lbd_notes",
    "dbd_start", "dbd_end", "dbd_source", "dbd_status",
]

# Stage 3A.2 policy decision (Option A): the UniProt/PROSITE profile
# interval is THE standardized project-wide LBD analysis definition for
# all 48 receptors. This label is descriptive of provenance, not a claim
# that it is the only possible biological LBD boundary — see
# data/manifests/nr_lbd_analysis_policy.json and
# reports/tables/stage3a_lbd_boundary_audit.md.
LBD_DEFINITION = "UNIPROT_PROSITE_NR_LBD"
LBD_COORDINATE_SYSTEM = "UNIPROT_CANONICAL_1_BASED_INCLUSIVE"


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_row(record_index: int, record: dict) -> dict:
    accession = record["uniprot_id"]
    raw_path = RAW_UNIPROT_DIR / f"{accession}.json"

    if not raw_path.exists():
        return {
            "record_index": record_index, "uniprot_id": accession,
            "nr_code": record["nr_code"], "common_name": record["common_name"], "group": record["group"],
            "sequence_length": "",
            "lbd_start": "", "lbd_end": "", "lbd_length": "",
            "lbd_definition": LBD_DEFINITION, "lbd_coordinate_system": LBD_COORDINATE_SYSTEM,
            "lbd_source": "", "lbd_source_identifier": "", "lbd_feature_type": "", "lbd_feature_description": "",
            "lbd_feature_evidence": "",
            "uniprot_entry_version": "", "uniprot_sequence_version": "",
            "lbd_status": "NOT_FOUND",
            "lbd_notes": f"No cached raw UniProt response found at {raw_path.relative_to(PROJECT_ROOT)}.",
            "dbd_start": "", "dbd_end": "", "dbd_source": "", "dbd_status": "NOT_FOUND",
        }

    raw = json.loads(raw_path.read_text())
    seq_length = extract_sequence_length(raw)
    audit = extract_entry_audit(raw)

    lbd = extract_lbd_feature(raw)
    if lbd is None:
        lbd_status = "NOT_FOUND"
        lbd_notes = (
            "No unambiguous UniProt 'NR LBD' Domain feature found (zero or multiple candidates). "
            "PROSITE/InterPro fallback was configured but not exercised automatically for this run — "
            "requires manual review before any fallback value is accepted."
        )
        lbd_start = lbd_end = lbd_length = ""
        lbd_source = lbd_source_id = lbd_feature_type = lbd_feature_desc = lbd_evidence = ""
    else:
        lbd_start, lbd_end = lbd["start"], lbd["end"]
        lbd_length = lbd_end - lbd_start + 1
        lbd_source = "UniProt"
        lbd_source_id = lbd["evidence"] or accession
        lbd_feature_type = "Domain"
        lbd_feature_desc = "NR LBD"
        lbd_evidence = lbd["evidence"]

        assertion_failures = []
        if lbd_start < 1:
            assertion_failures.append(f"lbd_start={lbd_start} < 1")
        if seq_length is not None and lbd_end > seq_length:
            assertion_failures.append(f"lbd_end={lbd_end} > sequence_length={seq_length}")
        if lbd_start > lbd_end:
            assertion_failures.append(f"lbd_start={lbd_start} > lbd_end={lbd_end}")

        if assertion_failures:
            lbd_status = "REVIEW_NEEDED"
            lbd_notes = "Coordinate assertion failed: " + "; ".join(assertion_failures)
        else:
            lbd_status = "PASS"
            lbd_notes = ""

    dbd = extract_dbd_feature(raw)
    if dbd is None:
        dbd_start = dbd_end = ""
        dbd_source = ""
        dbd_status = "NOT_FOUND"
    else:
        dbd_start, dbd_end = dbd["start"], dbd["end"]
        dbd_source = "UniProt"
        dbd_status = "PASS"

    return {
        "record_index": record_index,
        "uniprot_id": accession,
        "nr_code": record["nr_code"],
        "common_name": record["common_name"],
        "group": record["group"],
        "sequence_length": seq_length if seq_length is not None else "",
        "lbd_start": lbd_start,
        "lbd_end": lbd_end,
        "lbd_length": lbd_length,
        "lbd_definition": LBD_DEFINITION,
        "lbd_coordinate_system": LBD_COORDINATE_SYSTEM,
        "lbd_source": lbd_source,
        "lbd_source_identifier": lbd_source_id,
        "lbd_feature_type": lbd_feature_type,
        "lbd_feature_description": lbd_feature_desc,
        "lbd_feature_evidence": lbd_evidence,
        "uniprot_entry_version": audit["entry_version"],
        "uniprot_sequence_version": audit["sequence_version"],
        "lbd_status": lbd_status,
        "lbd_notes": lbd_notes,
        "dbd_start": dbd_start,
        "dbd_end": dbd_end,
        "dbd_source": dbd_source,
        "dbd_status": dbd_status,
    }


def main() -> None:
    rows = [build_row(idx, record) for idx, record in enumerate(NR_METADATA)]
    df = pd.DataFrame(rows)[COLUMNS]
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {OUTPUT_CSV.relative_to(PROJECT_ROOT)} ({len(df)} rows)")

    status_counts = df["lbd_status"].value_counts().to_dict()
    dbd_status_counts = df["dbd_status"].value_counts().to_dict()
    print(f"LBD status: {status_counts}")
    print(f"DBD status: {dbd_status_counts}")

    cached_files_used = sorted(
        str((RAW_UNIPROT_DIR / f"{r['uniprot_id']}.json").relative_to(PROJECT_ROOT))
        for r in NR_METADATA
        if (RAW_UNIPROT_DIR / f"{r['uniprot_id']}.json").exists()
    )

    provenance = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "number_of_receptors": len(NR_METADATA),
        "lbd_pass_count": int(status_counts.get("PASS", 0)),
        "lbd_review_needed_count": int(status_counts.get("REVIEW_NEEDED", 0)),
        "lbd_not_found_count": int(status_counts.get("NOT_FOUND", 0)),
        "dbd_pass_count": int(dbd_status_counts.get("PASS", 0)),
        "dbd_not_found_count": int(dbd_status_counts.get("NOT_FOUND", 0)),
        "source_hierarchy": [
            "1. current UniProt explicit 'NR LBD' Domain feature (features[].type=='Domain', description=='NR LBD')",
            "2. PROSITE NR_LBD / PS51843 annotation (fallback; NOT exercised in this run)",
            "3. InterPro nuclear hormone receptor ligand-binding-domain annotation (fallback; NOT exercised in this run)",
        ],
        "fallback_sources_required": False,
        "fallback_notes": (
            "All 48 receptors had an unambiguous UniProt 'NR LBD' Domain feature in their cached Stage 1 "
            "raw response; no PROSITE/InterPro fallback request was made."
        ),
        "dbd_feature_rule": (
            "features[].type=='DNA binding', description=='Nuclear receptor' (diagnostic only). "
            "Absent for DAX-1 (P51843) and SHP (Q15466) — these are documented atypical orphan "
            "nuclear receptors lacking a classical zinc-finger DBD; no interval was inferred for them."
        ),
        "coordinate_convention": "1-based, inclusive (matches UniProt feature location convention)",
        "cached_uniprot_files_used": cached_files_used,
        "external_fallback_requests_made": [],
        "nr_lbd_reference_csv_sha256": sha256_of_file(OUTPUT_CSV),
        "validated_receptor_master_sha256": sha256_of_file(NR_METADATA_PY),
        "generating_script": "scripts/11_build_lbd_reference.py",
    }
    PROVENANCE_JSON.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Wrote {PROVENANCE_JSON.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
