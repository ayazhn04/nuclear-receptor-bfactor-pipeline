"""
Stage 3A.1: sensitivity analysis — recompute LBD mapping coverage under
BOTH the profile-derived and previous-cohort LBD boundary definitions, for
every receptor where both are available (currently: VDR only), and report
whether the classification/prefilter consequence actually changes.

Does NOT overwrite data/manifests/rcsb_candidate_lbd_mapping.csv.

Reads:
    data/manifests/nr_lbd_boundary_audit.csv
    data/manifests/rcsb_sequence_coordinate_alignments.csv
    data/manifests/rcsb_discovery_by_receptor.csv
    data/manifests/rcsb_precoordinate_qc.csv   (for prefilter consequence comparison)

Writes:
    data/manifests/rcsb_lbd_boundary_sensitivity.csv
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.intervals import intersect_with_reference, merge_intervals, total_length  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
BOUNDARY_AUDIT_CSV = MANIFESTS_DIR / "nr_lbd_boundary_audit.csv"
ALIGNMENTS_CSV = MANIFESTS_DIR / "rcsb_sequence_coordinate_alignments.csv"
DISCOVERY_CSV = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
PRECOORDINATE_QC_CSV = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"

OUTPUT_CSV = MANIFESTS_DIR / "rcsb_lbd_boundary_sensitivity.csv"


def categorize(coverage: float) -> str:
    if coverage == 0:
        return "NO_LBD_OVERLAP"
    if coverage < 0.25:
        return "TRACE_LBD"
    if coverage < 0.80:
        return "PARTIAL_LBD"
    if coverage < 0.95:
        return "SUBSTANTIAL_LBD"
    return "NEAR_COMPLETE_LBD"


def main() -> None:
    boundary_audit = pd.read_csv(BOUNDARY_AUDIT_CSV)
    both_available = boundary_audit[boundary_audit["historical_range_available"]]

    if len(both_available) == 0:
        pd.DataFrame(
            columns=[
                "uniprot_id", "common_name", "pdb_id", "polymer_entity_id",
                "profile_lbd_start", "profile_lbd_end", "profile_lbd_mapping_coverage", "profile_lbd_mapping_category",
                "previous_lbd_start", "previous_lbd_end", "previous_lbd_mapping_coverage", "previous_lbd_mapping_category",
                "coverage_difference", "category_changed", "prefilter_consequence_changed",
            ]
        ).to_csv(OUTPUT_CSV, index=False)
        print(f"Wrote {OUTPUT_CSV.relative_to(PROJECT_ROOT)} (0 rows — no receptor has both definitions available)")
        return

    alignments = pd.read_csv(ALIGNMENTS_CSV, dtype={"pdb_id": str, "entity_id": str})
    alignments = alignments[alignments["is_stage2_candidate"]]
    discovery = pd.read_csv(DISCOVERY_CSV, dtype={"pdb_id": str, "entity_id": str})
    qc = pd.read_csv(PRECOORDINATE_QC_CSV, dtype={"pdb_id": str})

    rows = []
    for _, boundary_row in both_available.iterrows():
        uid = boundary_row["uniprot_id"]
        profile_start, profile_end = int(boundary_row["profile_lbd_start"]), int(boundary_row["profile_lbd_end"])
        profile_length = int(boundary_row["profile_lbd_length"])
        prev_start, prev_end = int(boundary_row["previous_cohort_lbd_start"]), int(boundary_row["previous_cohort_lbd_end"])
        prev_length = int(boundary_row["previous_cohort_lbd_length"])

        candidates = discovery[discovery["uniprot_id"] == uid]
        for _, cand in candidates.iterrows():
            key = (uid, cand["polymer_entity_id"])
            align_row = alignments[
                (alignments["uniprot_id"] == uid) & (alignments["polymer_entity_id"] == cand["polymer_entity_id"])
            ]
            if len(align_row) == 0:
                continue
            regions = json.loads(align_row.iloc[0]["aligned_regions_json"])
            query_intervals = merge_intervals([(r["query_begin"], r["query_end"]) for r in regions])

            profile_overlap = total_length(intersect_with_reference(query_intervals, profile_start, profile_end))
            profile_coverage = round(profile_overlap / profile_length, 6) if profile_length else 0.0
            profile_category = categorize(profile_coverage)

            prev_overlap = total_length(intersect_with_reference(query_intervals, prev_start, prev_end))
            prev_coverage = round(prev_overlap / prev_length, 6) if prev_length else 0.0
            prev_category = categorize(prev_coverage)

            category_changed = profile_category != prev_category

            qc_row = qc[(qc["uniprot_id"] == uid) & (qc["polymer_entity_id"] == cand["polymer_entity_id"])]
            prefilter_changed = ""
            if len(qc_row):
                qr = qc_row.iloc[0]
                current_prefilter = bool(qr["passes_project_metadata_prefilter"])
                would_overlap_prev_lbd = prev_overlap > 0
                # Recompute prefilter with the previous-cohort LBD-overlap
                # criterion substituted for the profile-based one, holding
                # method/resolution/R-free fixed.
                hypothetical_prefilter = (
                    bool(qr["method_is_xray"])
                    and qr["resolution_project_range_status"] == "PASS_PROJECT_RANGE"
                    and qr["r_free_status"] == "PASS_ALL_LE_0_30"
                    and would_overlap_prev_lbd
                )
                prefilter_changed = current_prefilter != hypothetical_prefilter

            rows.append(
                {
                    "uniprot_id": uid, "common_name": cand["common_name"],
                    "pdb_id": cand["pdb_id"], "polymer_entity_id": cand["polymer_entity_id"],
                    "profile_lbd_start": profile_start, "profile_lbd_end": profile_end,
                    "profile_lbd_mapping_coverage": profile_coverage, "profile_lbd_mapping_category": profile_category,
                    "previous_lbd_start": prev_start, "previous_lbd_end": prev_end,
                    "previous_lbd_mapping_coverage": prev_coverage, "previous_lbd_mapping_category": prev_category,
                    "coverage_difference": round(prev_coverage - profile_coverage, 6),
                    "category_changed": category_changed,
                    "prefilter_consequence_changed": prefilter_changed,
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {OUTPUT_CSV.relative_to(PROJECT_ROOT)} ({len(df)} rows)")
    if len(df):
        print(f"category_changed: {int(df['category_changed'].sum())} / {len(df)}")
        print(f"prefilter_consequence_changed: {int((df['prefilter_consequence_changed'] == True).sum())} / {len(df)}")


if __name__ == "__main__":
    main()
