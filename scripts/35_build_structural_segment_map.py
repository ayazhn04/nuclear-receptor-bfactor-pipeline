"""
Structural-segment alignment, Part 4 + Part 5: canonical H1-H12 / beta-sheet /
loop segment labelling per receptor, plus validation against anchor receptors
with well-characterized LBD architecture.

Method (scripts/utils/segment_topology.py, fully documented there): H1 is
anchored as the first consensus helix near the standardized LBD's own
N-terminus; H12 as the last consensus helix near its C-terminus (cross-checked
against this repository's own documented terminal-helix/AF-2 functional
features where available); beta_sheet as a single, uniquely identified
consensus sheet region flanked by helices on both sides. Interior helices are
labelled sequentially (H2-H5, H6-H11) ONLY when the observed interior helix
count between two already-anchored landmarks exactly matches the canonical
count -- never by raw ordinal position alone. Everything else is left
AMBIGUOUS or UNASSIGNED and reported as such.

Reads:
    data/analysis/receptor_dssp_consensus.parquet
    data/manifests/nr_lbd_reference.csv
    data/manifests/nr_functional_features.csv (terminal-helix cross-check)

Writes:
    data/analysis/nr_lbd_structural_segment_map.csv
    reports/tables/structural_segment_validation.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.segment_topology import assign_receptor_segments  # noqa: E402

ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"
TABLES = PROJECT_ROOT / "reports" / "tables"

CONSENSUS_PARQUET = ANALYSIS / "receptor_dssp_consensus.parquet"
LBD_REFERENCE_CSV = MANIFESTS / "nr_lbd_reference.csv"
FUNCTIONAL_FEATURES_CSV = MANIFESTS / "nr_functional_features.csv"

SEGMENT_MAP_CSV = ANALYSIS / "nr_lbd_structural_segment_map.csv"
VALIDATION_CSV = TABLES / "structural_segment_validation.csv"

ANCHOR_RECEPTORS = {"P11473": "VDR", "P37231": "PPARgamma", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}

CANONICAL_ORDER = {
    "H1": 1, "loop_H1_H2": 2, "H2": 3, "loop_H2_H3": 4, "H3": 5, "loop_H3_H4": 6, "H4": 7,
    "loop_H4_H5": 8, "H5": 9, "loop_H5_beta_sheet": 10, "beta_sheet": 11, "loop_beta_sheet_H6": 12,
    "H6": 13, "loop_H6_H7": 14, "H7": 15, "loop_H7_H8": 16, "H8": 17, "loop_H8_H9": 18, "H9": 19,
    "loop_H9_H10": 20, "H10": 21, "loop_H10_H11": 22, "H11": 23, "loop_H11_H12": 24, "H12": 25,
}


def terminal_helix_ranges(features: pd.DataFrame, uniprot_id: str) -> list[tuple[int, int]]:
    f = features[(features["uniprot_id"] == uniprot_id) & features["feature_name"].str.contains("H12|terminal helix|AF-2|AF2", case=False, regex=True)]
    return list(zip(f["feature_start"], f["feature_end"]))


def build_segment_rows(uniprot_id: str, nr_code: str, lbd_start: int, lbd_end: int,
                        consensus: pd.DataFrame, features: pd.DataFrame) -> list[dict]:
    if len(consensus) == 0:
        return [
            {
                "uniprot_id": uniprot_id, "nr_code": nr_code, "canonical_uniprot_position": p,
                "canonical_residue_name": "", "structural_segment": "", "structural_class": "NO_EVIDENCE",
                "segment_order": None, "assignment_method": "no_structures", "assignment_confidence": "NONE",
                "assignment_status": "UNASSIGNED", "evidence_notes": "No experimental structure coordinate data available for this receptor.",
            }
            for p in range(lbd_start, lbd_end + 1)
        ]

    names_lookup = consensus.set_index("canonical_uniprot_position")["canonical_residue_name"].to_dict()
    calls = assign_receptor_segments(
        consensus["canonical_uniprot_position"].tolist(), consensus["consensus_dssp"].tolist(),
        lbd_start, lbd_end, terminal_helix_ranges(features, uniprot_id),
    )
    calls = sorted(calls, key=lambda c: c.start)

    rows = []
    observed_positions = set(consensus["canonical_uniprot_position"])
    for p in range(lbd_start, lbd_end + 1):
        if p not in observed_positions:
            rows.append(
                {
                    "uniprot_id": uniprot_id, "nr_code": nr_code, "canonical_uniprot_position": p,
                    "canonical_residue_name": "", "structural_segment": "", "structural_class": "NO_EVIDENCE",
                    "segment_order": None, "assignment_method": "no_secondary_structure_evidence",
                    "assignment_confidence": "NONE", "assignment_status": "UNASSIGNED",
                    "evidence_notes": "Never observed with a positive-occupancy Ca in any representative instance of this receptor.",
                }
            )
            continue
        call = next((c for c in calls if c.start <= p <= c.end), None)
        if call is not None:
            rows.append(
                {
                    "uniprot_id": uniprot_id, "nr_code": nr_code, "canonical_uniprot_position": p,
                    "canonical_residue_name": names_lookup.get(p, ""),
                    "structural_segment": call.structural_segment, "structural_class": call.structural_class,
                    "segment_order": CANONICAL_ORDER.get(call.structural_segment),
                    "assignment_method": call.assignment_method, "assignment_confidence": call.assignment_confidence,
                    "assignment_status": call.assignment_status, "evidence_notes": call.evidence_notes,
                }
            )
            continue
        # not inside any HELIX/SHEET call -> loop/coil residue; bracket by neighbouring calls
        prev_call = max((c for c in calls if c.end < p), key=lambda c: c.end, default=None)
        next_call = min((c for c in calls if c.start > p), key=lambda c: c.start, default=None)
        prev_named = prev_call.structural_segment if prev_call else ""
        next_named = next_call.structural_segment if next_call else ""
        if prev_named and next_named:
            seg = f"loop_{prev_named}_{next_named}"
            rows.append(
                {
                    "uniprot_id": uniprot_id, "nr_code": nr_code, "canonical_uniprot_position": p,
                    "canonical_residue_name": names_lookup.get(p, ""), "structural_segment": seg,
                    "structural_class": "LOOP", "segment_order": (CANONICAL_ORDER.get(prev_named, 0) + CANONICAL_ORDER.get(next_named, 0)) / 2,
                    "assignment_method": "loop_bracketed_by_named_segments", "assignment_confidence": "MODERATE",
                    "assignment_status": "ASSIGNED_MODERATE_CONFIDENCE",
                    "evidence_notes": f"Consensus coil residue between confidently named {prev_named} and {next_named}.",
                }
            )
        else:
            reason = "before the first identifiable segment" if prev_call is None else (
                "after the last identifiable segment" if next_call is None else "between two not-confidently-named segments")
            rows.append(
                {
                    "uniprot_id": uniprot_id, "nr_code": nr_code, "canonical_uniprot_position": p,
                    "canonical_residue_name": names_lookup.get(p, ""), "structural_segment": "",
                    "structural_class": "LOOP_UNCERTAIN", "segment_order": None,
                    "assignment_method": "loop_unbracketed", "assignment_confidence": "NONE",
                    "assignment_status": "AMBIGUOUS" if (prev_call or next_call) else "UNASSIGNED",
                    "evidence_notes": f"Consensus coil residue {reason}; not assigned a loop label.",
                }
            )
    return rows


def build_validation_rows(seg_map: pd.DataFrame, consensus: pd.DataFrame) -> pd.DataFrame:
    named = seg_map[seg_map["structural_segment"] != ""]
    rows = []
    for (uid, seg), g in named.groupby(["uniprot_id", "structural_segment"]):
        c = consensus[(consensus["uniprot_id"] == uid) & consensus["canonical_uniprot_position"].isin(g["canonical_uniprot_position"])]
        rows.append(
            {
                "uniprot_id": uid, "nr_code": g["nr_code"].iloc[0], "structural_segment": seg,
                "structural_class": g["structural_class"].iloc[0],
                "start_uniprot": int(g["canonical_uniprot_position"].min()), "end_uniprot": int(g["canonical_uniprot_position"].max()),
                "n_residues": len(g), "supporting_pdbs": int(c["n_pdb_observed"].max()) if len(c) else 0,
                "median_consensus_fraction": float(c["consensus_fraction"].median()) if len(c) else float("nan"),
                "assignment_confidence": g["assignment_confidence"].iloc[0],
                "validation_notes": g["evidence_notes"].iloc[0],
            }
        )
    out = pd.DataFrame(rows).sort_values(["uniprot_id", "start_uniprot"]).reset_index(drop=True)
    return out


def main() -> None:
    consensus_all = pd.read_parquet(CONSENSUS_PARQUET)
    lbd_ref = pd.read_csv(LBD_REFERENCE_CSV)
    features = pd.read_csv(FUNCTIONAL_FEATURES_CSV)

    all_rows = []
    for _, r in lbd_ref.iterrows():
        uid = r["uniprot_id"]
        c = consensus_all[consensus_all["uniprot_id"] == uid]
        all_rows.extend(build_segment_rows(uid, r["nr_code"], int(r["lbd_start"]), int(r["lbd_end"]), c, features))

    seg_map = pd.DataFrame(all_rows)
    assert not seg_map.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    seg_map.to_csv(SEGMENT_MAP_CSV, index=False)
    print(f"Wrote {SEGMENT_MAP_CSV.relative_to(PROJECT_ROOT)} ({len(seg_map)} rows)")
    named = seg_map[seg_map["structural_segment"] != ""]
    print(f"  {named['uniprot_id'].nunique()} receptors have >=1 named segment; "
          f"{seg_map['uniprot_id'].nunique() - named['uniprot_id'].nunique()} have none")
    print(seg_map["assignment_status"].value_counts().to_string())

    print("\nAnchor-receptor H1/H12/beta_sheet coverage:")
    canonical_helices = {f"H{i}" for i in range(1, 13)} | {"beta_sheet"}
    for uid, name in ANCHOR_RECEPTORS.items():
        got = set(seg_map.loc[seg_map["uniprot_id"] == uid, "structural_segment"]) & canonical_helices
        print(f"  {name} ({uid}): {sorted(got, key=lambda s: CANONICAL_ORDER.get(s, 99))}")

    validation = build_validation_rows(seg_map, consensus_all)
    TABLES.mkdir(parents=True, exist_ok=True)
    validation.to_csv(VALIDATION_CSV, index=False)
    print(f"\nWrote {VALIDATION_CSV.relative_to(PROJECT_ROOT)} ({len(validation)} rows)")

    missing_major = []
    for uid, name in ANCHOR_RECEPTORS.items():
        got = set(seg_map.loc[seg_map["uniprot_id"] == uid, "structural_segment"])
        for major in ["H3", "H5", "H10", "H11", "H12"]:
            if major not in got:
                missing_major.append((name, major))
    if missing_major:
        print("\nSTOP-CONDITION CHECK: major canonical elements not consistently mapped in anchor receptors:")
        for name, major in missing_major:
            print(f"  {name}: {major} not assigned")
    else:
        print("\nSTOP-CONDITION CHECK: all of H3/H5/H10/H11/H12 mapped in all 5 anchor receptors.")


if __name__ == "__main__":
    main()
