"""
Structural-segment alignment v2, Part 8 + Part 9: rebuild the per-position
canonical segment map from the reference-based 3D structural mapping (own
canonical template for the 5 anchors; CE structural-alignment label
transfer, boundary-refined against each target's own secondary structure,
for the other 37 receptors), and check the acceptance criteria before any
feature matrix is rebuilt.

Reads:
    data/manifests/nr_lbd_reference.csv
    data/analysis/canonical_lbd_template_segments.csv
    data/analysis/lbd_structural_label_transfer.csv
    data/analysis/receptor_dssp_consensus.parquet (v1 output, reused: still a
        valid, unchanged evidence source -- consensus secondary structure
        does not depend on the labelling method)

Writes:
    data/analysis/nr_lbd_structural_segment_map.csv   (OVERWRITES the v1 file;
        v1's own commit 9f15062 remains an intact diagnostic baseline in git
        history)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"

LBD_REFERENCE_CSV = MANIFESTS / "nr_lbd_reference.csv"
TEMPLATE_CSV = ANALYSIS / "canonical_lbd_template_segments.csv"
TRANSFER_CSV = ANALYSIS / "lbd_structural_label_transfer.csv"
CONSENSUS_PARQUET = ANALYSIS / "receptor_dssp_consensus.parquet"

SEGMENT_MAP_CSV = ANALYSIS / "nr_lbd_structural_segment_map.csv"

CANONICAL_ORDER = {
    "H1": 1, "loop_H1_H2": 2, "H2": 3, "loop_H2_H3": 4, "H3": 5, "loop_H3_H4": 6, "H4": 7,
    "loop_H4_H5": 8, "H5": 9, "loop_H5_beta_sheet": 10, "beta_sheet": 11, "loop_beta_sheet_H6": 12,
    "H6": 13, "loop_H6_H7": 14, "H7": 15, "loop_H7_H8": 16, "H8": 17, "loop_H8_H9": 18, "H9": 19,
    "loop_H9_H10": 20, "H10": 21, "loop_H10_H11": 22, "H11": 23, "loop_H11_H12": 24, "H12": 25,
}
CORE_SEGMENTS = ["H1", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "beta_sheet"]


def build_named_lookup() -> pd.DataFrame:
    tmpl = pd.read_csv(TEMPLATE_CSV, dtype={"pdb_id": str}, keep_default_na=False)
    tmpl_named = tmpl[tmpl["canonical_segment"] != ""][
        ["uniprot_id", "canonical_uniprot_position", "canonical_segment", "structural_class",
         "assignment_confidence", "assignment_status", "evidence_source"]
    ].rename(columns={"canonical_segment": "structural_segment", "evidence_source": "evidence_notes"})
    tmpl_named["reference_template"] = tmpl_named["uniprot_id"]
    tmpl_named["reference_pdb"] = tmpl.set_index(["uniprot_id", "canonical_uniprot_position"]).loc[
        list(zip(tmpl_named["uniprot_id"], tmpl_named["canonical_uniprot_position"])), "pdb_id"
    ].values
    tmpl_named["structural_alignment_score"] = pd.NA
    tmpl_named["structural_alignment_coverage"] = pd.NA
    tmpl_named["assignment_method"] = "own_canonical_template"

    transfer = pd.read_csv(TRANSFER_CSV, keep_default_na=False)
    transfer["assignment_status"] = transfer["assignment_status_hint"].replace("", "AMBIGUOUS")
    transfer = transfer.drop(columns=["assignment_status_hint"])
    transfer["assignment_status"] = transfer["assignment_confidence"].map(
        {"HIGH": "ASSIGNED_HIGH_CONFIDENCE", "MODERATE": "ASSIGNED_MODERATE_CONFIDENCE", "LOW": "AMBIGUOUS"}
    )

    cols = ["uniprot_id", "canonical_uniprot_position", "structural_segment", "structural_class",
            "reference_template", "reference_pdb", "structural_alignment_score", "structural_alignment_coverage",
            "assignment_method", "assignment_confidence", "assignment_status", "evidence_notes"]
    return pd.concat([tmpl_named[cols], transfer[cols]], ignore_index=True)


def build_rows(uid: str, nr_code: str, lbd_start: int, lbd_end: int,
               named: pd.DataFrame, consensus: pd.DataFrame) -> list[dict]:
    named_here = named[named["uniprot_id"] == uid].sort_values("canonical_uniprot_position")
    if len(named_here) == 0 and len(consensus) == 0:
        return [
            {
                "uniprot_id": uid, "nr_code": nr_code, "canonical_uniprot_position": p, "canonical_residue_name": "",
                "structural_segment": "", "structural_class": "NO_EVIDENCE", "segment_order": None,
                "reference_template": "", "reference_pdb": "", "structural_alignment_score": None,
                "structural_alignment_coverage": None, "assignment_method": "no_structures",
                "assignment_confidence": "NONE", "assignment_status": "UNASSIGNED",
                "evidence_notes": "No experimental structure coordinate data available for this receptor.",
            }
            for p in range(lbd_start, lbd_end + 1)
        ]

    named_by_pos = {int(r["canonical_uniprot_position"]): r for _, r in named_here.iterrows()}
    observed_positions = set(consensus["canonical_uniprot_position"]) if len(consensus) else set(named_by_pos)
    names_lookup = consensus.set_index("canonical_uniprot_position")["canonical_residue_name"].to_dict() if len(consensus) else {}

    calls = []
    for _, r in named_here.iterrows():
        if not calls or calls[-1]["structural_segment"] != r["structural_segment"] or int(r["canonical_uniprot_position"]) != calls[-1]["end"] + 1:
            calls.append({"start": int(r["canonical_uniprot_position"]), "end": int(r["canonical_uniprot_position"]), **r.to_dict()})
        else:
            calls[-1]["end"] = int(r["canonical_uniprot_position"])
    calls.sort(key=lambda c: c["start"])

    rows = []
    for p in range(lbd_start, lbd_end + 1):
        if p not in observed_positions and p not in named_by_pos:
            rows.append({"uniprot_id": uid, "nr_code": nr_code, "canonical_uniprot_position": p, "canonical_residue_name": "",
                         "structural_segment": "", "structural_class": "NO_EVIDENCE", "segment_order": None,
                         "reference_template": "", "reference_pdb": "", "structural_alignment_score": None,
                         "structural_alignment_coverage": None, "assignment_method": "no_secondary_structure_evidence",
                         "assignment_confidence": "NONE", "assignment_status": "UNASSIGNED",
                         "evidence_notes": "Never observed with a positive-occupancy Ca in any representative instance of this receptor."})
            continue
        if p in named_by_pos:
            r = named_by_pos[p]
            rows.append({"uniprot_id": uid, "nr_code": nr_code, "canonical_uniprot_position": p,
                         "canonical_residue_name": names_lookup.get(p, ""),
                         "structural_segment": r["structural_segment"], "structural_class": r["structural_class"],
                         "segment_order": CANONICAL_ORDER.get(r["structural_segment"]),
                         "reference_template": r["reference_template"], "reference_pdb": r["reference_pdb"],
                         "structural_alignment_score": r["structural_alignment_score"],
                         "structural_alignment_coverage": r["structural_alignment_coverage"],
                         "assignment_method": r["assignment_method"], "assignment_confidence": r["assignment_confidence"],
                         "assignment_status": r["assignment_status"], "evidence_notes": r["evidence_notes"]})
            continue
        prev_call = max((c for c in calls if c["end"] < p), key=lambda c: c["end"], default=None)
        next_call = min((c for c in calls if c["start"] > p), key=lambda c: c["start"], default=None)
        prev_named = prev_call["structural_segment"] if prev_call else ""
        next_named = next_call["structural_segment"] if next_call else ""
        if prev_named and next_named:
            seg = f"loop_{prev_named}_{next_named}"
            rows.append({"uniprot_id": uid, "nr_code": nr_code, "canonical_uniprot_position": p,
                         "canonical_residue_name": names_lookup.get(p, ""), "structural_segment": seg,
                         "structural_class": "LOOP", "segment_order": (CANONICAL_ORDER.get(prev_named, 0) + CANONICAL_ORDER.get(next_named, 0)) / 2,
                         "reference_template": prev_call["reference_template"], "reference_pdb": prev_call["reference_pdb"],
                         "structural_alignment_score": prev_call["structural_alignment_score"],
                         "structural_alignment_coverage": prev_call["structural_alignment_coverage"],
                         "assignment_method": "loop_bracketed_by_named_segments", "assignment_confidence": "MODERATE",
                         "assignment_status": "ASSIGNED_MODERATE_CONFIDENCE",
                         "evidence_notes": f"Consensus coil residue between confidently named {prev_named} and {next_named}."})
        else:
            reason = "before the first identifiable segment" if prev_call is None else (
                "after the last identifiable segment" if next_call is None else "between two not-confidently-named segments")
            rows.append({"uniprot_id": uid, "nr_code": nr_code, "canonical_uniprot_position": p,
                         "canonical_residue_name": names_lookup.get(p, ""), "structural_segment": "",
                         "structural_class": "LOOP_UNCERTAIN", "segment_order": None,
                         "reference_template": "", "reference_pdb": "", "structural_alignment_score": None,
                         "structural_alignment_coverage": None, "assignment_method": "loop_unbracketed",
                         "assignment_confidence": "NONE",
                         "assignment_status": "AMBIGUOUS" if (prev_call or next_call) else "UNASSIGNED",
                         "evidence_notes": f"Consensus/observed residue {reason}; not assigned a loop label."})
    return rows


def main() -> None:
    lbd_ref = pd.read_csv(LBD_REFERENCE_CSV)
    consensus_all = pd.read_parquet(CONSENSUS_PARQUET)
    named = build_named_lookup()

    all_rows = []
    for _, r in lbd_ref.iterrows():
        uid = r["uniprot_id"]
        c = consensus_all[consensus_all["uniprot_id"] == uid]
        all_rows.extend(build_rows(uid, r["nr_code"], int(r["lbd_start"]), int(r["lbd_end"]), named, c))

    seg_map = pd.DataFrame(all_rows)
    assert not seg_map.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    seg_map.to_csv(SEGMENT_MAP_CSV, index=False)
    print(f"Wrote {SEGMENT_MAP_CSV.relative_to(PROJECT_ROOT)} ({len(seg_map)} rows)")
    print(seg_map["assignment_status"].value_counts().to_string())

    print("\n--- PART 9 ACCEPTANCE CRITERIA ---")
    named_only = seg_map[seg_map["structural_segment"] != ""]
    counts = named_only.groupby("structural_segment")["uniprot_id"].nunique()
    n_with_coord = consensus_all["uniprot_id"].nunique()
    print(f"Receptors with coordinate/secondary-structure evidence: {n_with_coord}")
    for seg in CORE_SEGMENTS:
        n = int(counts.get(seg, 0))
        print(f"  {seg:12s}: {n:3d} / {n_with_coord} receptors ({100 * n / n_with_coord:5.1f}%)")

    core_helices = ["H3", "H5", "H10", "H11", "H12"]
    failing = [s for s in core_helices if counts.get(s, 0) < n_with_coord * 0.5]
    if failing:
        print(f"\nSTOP-CONDITION: {failing} map for fewer than 50% of receptors with coordinate data.")
        print("This is reported, not silently worked around -- see the structural-alignment score")
        print("distribution (data/analysis/lbd_structural_alignment_scores.csv) and the per-segment")
        print("template coverage (data/analysis/canonical_lbd_template_segments.csv) for why.")
    else:
        print("\nAcceptance criteria PASSED for H3/H5/H10/H11/H12.")


if __name__ == "__main__":
    main()
