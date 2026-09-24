"""
Stage 3A.1: reconcile the UniProt/PROSITE profile-derived LBD interval
(Stage 3A) against the previous BIOL363 cohort's LBD window definition.

Reads:
    data/manifests/nr_lbd_reference.csv               (Stage 3A profile-derived intervals — UNCHANGED)
    data/raw/external/previous_cohort_repo/proteins.csv
    data/raw/api/uniprot/{uniprot_id}.json

Writes:
    data/manifests/nr_lbd_boundary_audit.csv           (48 rows)
    data/manifests/nr_lbd_boundary_audit_provenance.json

Previous-cohort data availability (see report for full detail): the
previous cohort's public repository
(https://github.com/sumeshi3648/BioInformatics-Bfactors.git,
commit 69899f60d7ab0eed19d3029687fe31f03b835ac2) contains NO 48-receptor
LBD-range metadata sheet. `proteins.csv` lists only subfamily/symbol/
uniprot (no coordinates). The only explicit LBD window found anywhere in
the repository is a HARD-CODED script default of 118-427, used only in
VDR-specific analysis scripts (plot_vdr_final.py, compute_dssp_vdr.py,
consensus_sse_vdr.py, build_vdr_consensus_sse.py,
plot_vdr_violin_sse.py, summarize_vdr_flexibility_by_region.py) — not
stored in any data file, and not defined for any of the other 47
receptors. This script therefore populates `previous_cohort_lbd_*` for
VDR only; all other 47 receptors are correctly marked
PREVIOUS_RANGE_UNAVAILABLE — never fabricated.

This script does NOT modify data/manifests/nr_lbd_reference.csv.
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
LBD_REFERENCE_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"
RAW_UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"
PREV_COHORT_DIR = PROJECT_ROOT / "data" / "raw" / "external" / "previous_cohort_repo"

OUTPUT_CSV = MANIFESTS_DIR / "nr_lbd_boundary_audit.csv"
PROVENANCE_JSON = MANIFESTS_DIR / "nr_lbd_boundary_audit_provenance.json"

# The ONLY previous-cohort LBD window found: VDR, hard-coded in 6 scripts
# (script default arguments `--start 118 --end 427`), not a data file.
PREVIOUS_COHORT_VDR_LBD = {"uniprot_id": "P11473", "start": 118, "end": 427}
PREVIOUS_COHORT_SOURCE = (
    "https://github.com/sumeshi3648/BioInformatics-Bfactors.git "
    "@69899f60d7ab0eed19d3029687fe31f03b835ac2 "
    "(scripts/plot_vdr_final.py, compute_dssp_vdr.py, consensus_sse_vdr.py, "
    "build_vdr_consensus_sse.py, plot_vdr_violin_sse.py, "
    "summarize_vdr_flexibility_by_region.py — hard-coded script defaults, not a data file)"
)

COLUMNS = [
    "uniprot_id", "nr_code", "common_name", "sequence_length",
    "profile_lbd_start", "profile_lbd_end", "profile_lbd_length", "profile_source", "profile_source_identifier",
    "previous_cohort_lbd_start", "previous_cohort_lbd_end", "previous_cohort_lbd_length", "previous_cohort_source",
    "start_difference_previous_vs_profile", "end_difference_previous_vs_profile", "length_difference_previous_vs_profile",
    "historical_range_available",
    "boundary_audit_status", "boundary_audit_notes",
]


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def classify(start_diff: int | None, end_diff: int | None) -> tuple[str, str]:
    if start_diff is None or end_diff is None:
        return "PREVIOUS_RANGE_UNAVAILABLE", "No previous-cohort LBD range found for this receptor."
    if start_diff == 0 and end_diff == 0:
        return "MATCH", "Profile and previous-cohort intervals are identical."
    max_abs_diff = max(abs(start_diff), abs(end_diff))
    # Thresholds set AFTER examining the only available data point (VDR:
    # start diff 9, end diff 4) — a difference reaching into a few residues
    # near a domain boundary (esp. C-terminal, where H12/AF-2 sits) is
    # treated as material rather than assumed negligible, per instruction
    # not to assume small C-terminal differences are irrelevant.
    if max_abs_diff <= 2:
        return "MINOR_DIFFERENCE", f"Small boundary difference (max |diff|={max_abs_diff})."
    return "MATERIAL_DIFFERENCE", f"Boundary difference reaches {max_abs_diff} residue(s) at at least one terminus."


def main() -> None:
    lbd_ref = pd.read_csv(LBD_REFERENCE_CSV).set_index("uniprot_id")

    rows = []
    for record in NR_METADATA:
        uid = record["uniprot_id"]
        ref_row = lbd_ref.loc[uid]

        if uid == PREVIOUS_COHORT_VDR_LBD["uniprot_id"]:
            prev_start = PREVIOUS_COHORT_VDR_LBD["start"]
            prev_end = PREVIOUS_COHORT_VDR_LBD["end"]
            prev_length = prev_end - prev_start + 1
            historical_available = True
            prev_source = PREVIOUS_COHORT_SOURCE
            start_diff = prev_start - int(ref_row["lbd_start"])
            end_diff = prev_end - int(ref_row["lbd_end"])
            length_diff = prev_length - int(ref_row["lbd_length"])
        else:
            prev_start = prev_end = prev_length = ""
            historical_available = False
            prev_source = ""
            start_diff = end_diff = length_diff = None

        status, notes = classify(start_diff, end_diff)

        rows.append(
            {
                "uniprot_id": uid, "nr_code": record["nr_code"], "common_name": record["common_name"],
                "sequence_length": ref_row["sequence_length"],
                "profile_lbd_start": ref_row["lbd_start"], "profile_lbd_end": ref_row["lbd_end"],
                "profile_lbd_length": ref_row["lbd_length"],
                "profile_source": ref_row["lbd_source"], "profile_source_identifier": ref_row["lbd_source_identifier"],
                "previous_cohort_lbd_start": prev_start, "previous_cohort_lbd_end": prev_end,
                "previous_cohort_lbd_length": prev_length, "previous_cohort_source": prev_source,
                "start_difference_previous_vs_profile": start_diff if start_diff is not None else "",
                "end_difference_previous_vs_profile": end_diff if end_diff is not None else "",
                "length_difference_previous_vs_profile": length_diff if length_diff is not None else "",
                "historical_range_available": historical_available,
                "boundary_audit_status": status, "boundary_audit_notes": notes,
            }
        )

    df = pd.DataFrame(rows)[COLUMNS]
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {OUTPUT_CSV.relative_to(PROJECT_ROOT)} ({len(df)} rows)")
    print(df["boundary_audit_status"].value_counts().to_string())

    provenance = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "previous_cohort_repository": "https://github.com/sumeshi3648/BioInformatics-Bfactors.git",
        "previous_cohort_commit": "69899f60d7ab0eed19d3029687fe31f03b835ac2",
        "previous_cohort_retrieval_timestamp_utc": "2026-09-23T18:20:00+00:00",
        "previous_cohort_files_inspected": [
            "proteins.csv", "data/meta/structures.csv", "data/processed/ca_bfactors.csv",
            "scripts/*.py (all 14 scripts)",
        ],
        "previous_cohort_finding": (
            "No 48-receptor LBD-range metadata sheet exists in the repository. proteins.csv has "
            "only subfamily/symbol/uniprot columns (48 rows, cross-matches this project's 48 "
            "receptors by UniProt accession). The only explicit LBD window found anywhere is a "
            "hard-coded VDR-only script default (118-427), present in 6 VDR-specific analysis "
            "scripts, not in any data file, and not defined for any other receptor."
        ),
        "snapshotted_source_files": {
            str(p.relative_to(PROJECT_ROOT)): sha256_of_file(p)
            for p in sorted(PREV_COHORT_DIR.rglob("*")) if p.is_file()
        },
        "profile_lbd_reference_unchanged": sha256_of_file(LBD_REFERENCE_CSV),
        "boundary_audit_csv_sha256": sha256_of_file(OUTPUT_CSV),
        "generating_script": "scripts/15_build_lbd_boundary_audit.py",
    }
    PROVENANCE_JSON.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Wrote {PROVENANCE_JSON.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
