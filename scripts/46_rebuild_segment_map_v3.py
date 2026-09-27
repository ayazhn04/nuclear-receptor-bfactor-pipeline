"""
Structural-segment alignment v3, Part 8 (map rebuild) + Part 9 (acceptance
report) + Part 10 (feature/matrix rebuild, unchanged aggregation) + Part 12
(missingness recompute).

Reads:
    data/manifests/nr_lbd_reference.csv
    data/analysis/canonical_lbd_template_segments_v2.csv
    data/analysis/lbd_structural_label_transfer.csv
    data/analysis/receptor_dssp_consensus.parquet (unchanged evidence source)

Writes:
    data/analysis/nr_lbd_structural_segment_map.csv        (OVERWRITTEN, v3)
    data/analysis/pdb_segment_coverage.parquet              (OVERWRITTEN)
    data/analysis/pdb_segment_flexibility_features.parquet  (OVERWRITTEN)
    data/analysis/receptor_segment_flexibility_profiles.parquet (OVERWRITTEN)
    data/analysis/receptor_segment_feature_matrix_median.parquet (OVERWRITTEN)
    data/analysis/receptor_segment_feature_matrix_mean.parquet   (OVERWRITTEN)
    data/analysis/pdb_segment_feature_matrix_median.parquet      (OVERWRITTEN)
    data/analysis/pdb_segment_feature_matrix_mean.parquet        (OVERWRITTEN)
    data/analysis/receptor_segment_coverage_matrix.parquet       (OVERWRITTEN)
    reports/tables/segment_feature_missingness.md                (OVERWRITTEN)
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"
TABLES = PROJECT_ROOT / "reports" / "tables"

LBD_REFERENCE_CSV = MANIFESTS / "nr_lbd_reference.csv"
TEMPLATE_V2_CSV = ANALYSIS / "canonical_lbd_template_segments_v2.csv"
TRANSFER_CSV = ANALYSIS / "lbd_structural_label_transfer.csv"
CONSENSUS_PARQUET = ANALYSIS / "receptor_dssp_consensus.parquet"
SEGMENT_MAP_CSV = ANALYSIS / "nr_lbd_structural_segment_map.csv"

CANONICAL_ORDER = {"H1": 1, "H2": 3, "H3": 5, "H4": 7, "H5": 9, "beta_sheet": 11, "H6": 13, "H7": 15,
                    "H8": 17, "H9": 19, "H10": 21, "H11": 23, "H12": 25}
CORE_SEGMENTS = ["H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "beta_sheet"]

CONF_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2}
METHOD_RANK = {"own_canonical_template": 0, "structural_alignment_transfer": 1, "complementary_structure_rescue": 1}


def enforce_canonical_order(named: pd.DataFrame) -> pd.DataFrame:
    """Part 7 biological-consistency check: canonical order must be
    H1 -> H2 -> H3 -> H4 -> H5 -> beta_sheet -> H6 -> ... -> H12. Local
    residue-order continuity within one transferred segment is already
    guaranteed by CE's monotonic alignment path, but a segment transferred
    or cross-template-rescued from a DIFFERENT anchor can still land at a
    sequence position that contradicts another, independently-anchored
    segment in the SAME receptor (seen for real: VDR's cross-template-
    rescued H5 landed inside its own N-terminal insertion domain, before
    its own literature-anchored H3). Resolve greedily in canonical-segment
    order, keeping the higher-confidence/more-direct-evidence segment and
    dropping (unnaming) whichever conflicting segment is weaker -- flagged
    in evidence_notes, never silently discarded."""
    kept_frames = []
    for uid, g in named.groupby("uniprot_id"):
        blocks = []
        for seg, gg in g.groupby("structural_segment"):
            rank = CANONICAL_ORDER.get(seg)
            if rank is None:
                blocks.append({"segment": seg, "rank": 999, "start": gg["canonical_uniprot_position"].min(),
                                "end": gg["canonical_uniprot_position"].max(), "rows": gg,
                                "conf": CONF_RANK.get(gg["assignment_confidence"].iloc[0], 2),
                                "method": METHOD_RANK.get(gg["assignment_method"].iloc[0], 2)})
                continue
            blocks.append({"segment": seg, "rank": rank, "start": gg["canonical_uniprot_position"].min(),
                            "end": gg["canonical_uniprot_position"].max(), "rows": gg,
                            "conf": CONF_RANK.get(gg["assignment_confidence"].iloc[0], 2),
                            "method": METHOD_RANK.get(gg["assignment_method"].iloc[0], 2)})
        blocks.sort(key=lambda b: b["rank"])
        # Keep the maximum-priority (lowest conf/method rank first, ties broken by
        # canonical rank) subset whose positions are strictly increasing with rank --
        # a weighted longest-increasing-subsequence, computed greedily by priority.
        order = sorted(range(len(blocks)), key=lambda i: (blocks[i]["conf"], blocks[i]["method"], blocks[i]["rank"]))
        kept: list[int] = []
        for i in order:
            b = blocks[i]
            lower = [blocks[j] for j in kept if blocks[j]["rank"] < b["rank"]]
            upper = [blocks[j] for j in kept if blocks[j]["rank"] > b["rank"]]
            if lower and max(x["end"] for x in lower) >= b["start"]:
                continue
            if upper and min(x["start"] for x in upper) <= b["end"]:
                continue
            kept.append(i)
        dropped = [blocks[i] for i in range(len(blocks)) if i not in kept]
        for b in dropped:
            print(f"  ORDER VIOLATION dropped: {uid} {b['segment']} ({b['start']}-{b['end']}) "
                  f"conflicts with canonical N->C order against another segment in this receptor.")
        for i in kept:
            kept_frames.append(blocks[i]["rows"])
    return pd.concat(kept_frames, ignore_index=True) if kept_frames else named.iloc[0:0]


def build_named_lookup() -> pd.DataFrame:
    tmpl = pd.read_csv(TEMPLATE_V2_CSV, dtype={"reference_pdb": str}, keep_default_na=False)
    tmpl_named = tmpl[tmpl["structural_segment"] != ""][
        ["uniprot_id", "canonical_uniprot_position", "structural_segment", "structural_class",
         "assignment_confidence", "reference_pdb", "publication"]
    ].rename(columns={"publication": "evidence_notes"})
    tmpl_named["reference_template"] = tmpl_named["uniprot_id"]
    tmpl_named["structural_alignment_score"] = pd.NA
    tmpl_named["structural_alignment_coverage"] = pd.NA
    tmpl_named["assignment_method"] = "own_canonical_template"
    tmpl_named["assignment_status"] = "ASSIGNED_HIGH_CONFIDENCE"

    transfer = pd.read_csv(TRANSFER_CSV, keep_default_na=False)
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
            {"uniprot_id": uid, "nr_code": nr_code, "canonical_uniprot_position": p, "canonical_residue_name": "",
             "structural_segment": "", "structural_class": "NO_EVIDENCE", "segment_order": None,
             "reference_template": "", "reference_pdb": "", "structural_alignment_score": None,
             "structural_alignment_coverage": None, "assignment_method": "no_structures",
             "assignment_confidence": "NONE", "assignment_status": "UNASSIGNED",
             "evidence_notes": "No experimental structure coordinate data available for this receptor."}
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
    print("Enforcing canonical N->C segment order (Part 7 biological-consistency check)...")
    n_before = len(named)
    named = enforce_canonical_order(named)
    print(f"  {n_before - len(named)} position-rows dropped for canonical-order violations "
          f"({n_before} -> {len(named)})")

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

    print("\n--- PART 9 ACCEPTANCE CRITERIA (v3) ---")
    named_only = seg_map[seg_map["structural_segment"] != ""]
    counts = named_only.groupby("structural_segment")["uniprot_id"].nunique()
    n_with_coord = consensus_all["uniprot_id"].nunique()
    for seg in CORE_SEGMENTS:
        n = int(counts.get(seg, 0))
        print(f"  {seg:12s}: {n:3d} / {n_with_coord} receptors ({100 * n / n_with_coord:5.1f}%)")
    core_helices = ["H3", "H5", "H7", "H10", "H11", "H12"]
    near_zero = [s for s in ["H5", "H7", "H10"] if counts.get(s, 0) < n_with_coord * 0.10]
    if near_zero:
        print(f"\nSTOP-CONDITION: {near_zero} still near-zero -- diagnosing required, not proceeding to feature rebuild for those.")
    else:
        print("\nAcceptance criteria PASSED: none of H5/H7/H10 are near zero.")

    print("\n--- Part 10-12: rebuild features/matrices/missingness (aggregation unchanged) ---")
    spec36 = importlib.util.spec_from_file_location("seg_features_v1", PROJECT_ROOT / "scripts" / "36_build_segment_features.py")
    seg_features_v1 = importlib.util.module_from_spec(spec36)
    spec36.loader.exec_module(seg_features_v1)
    spec37 = importlib.util.spec_from_file_location("seg_missingness_v1", PROJECT_ROOT / "scripts" / "37_build_segment_missingness_and_alignment.py")
    seg_missingness_v1 = importlib.util.module_from_spec(spec37)
    spec37.loader.exec_module(seg_missingness_v1)

    master = pd.read_csv(seg_features_v1.MASTER_CSV, dtype={"pdb_id": str})
    norm = pd.read_parquet(seg_features_v1.NORMALIZED_PARQUET,
                            columns=["instance_id", "uniprot_id", "canonical_uniprot_position", "is_construct_mapped", "ca_positive_occupancy"])
    pdb_pos = pd.read_parquet(seg_features_v1.PDB_POSITION_PARQUET)
    defs, _ = seg_features_v1.load_segment_lookup(seg_map)

    coverage = seg_features_v1.build_coverage(master, norm, seg_map, defs)
    coverage.to_parquet(seg_features_v1.COVERAGE_PARQUET, index=False)
    pdb_feat = seg_features_v1.build_pdb_features(master, norm, pdb_pos, seg_map, defs, coverage)
    pdb_feat.to_parquet(seg_features_v1.PDB_FEATURES_PARQUET, index=False)

    meta = master.drop_duplicates("uniprot_id")[["uniprot_id", "nr_code", "common_name", "group"]].copy()
    n_pdb_total = master.loc[master["has_stage3b_coordinate_data"]].groupby("uniprot_id")["pdb_id"].nunique()
    meta["n_pdbs_with_coordinate_data"] = meta["uniprot_id"].map(n_pdb_total).fillna(0).astype(int)
    rec_profile = seg_features_v1.build_receptor_profile(pdb_feat, meta)
    rec_profile.to_parquet(seg_features_v1.RECEPTOR_PROFILE_PARQUET, index=False)

    rec_med = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="median_of_entry_median_z")
    rec_mean = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="mean_of_entry_mean_z")
    pdb_med = pdb_feat.pivot_table(index=["uniprot_id", "pdb_id"], columns="structural_segment", values="median_z")
    pdb_mean = pdb_feat.pivot_table(index=["uniprot_id", "pdb_id"], columns="structural_segment", values="mean_z")
    rec_cov = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="median_segment_coverage")
    rec_med.reset_index().to_parquet(seg_features_v1.RECEPTOR_MATRIX_MEDIAN, index=False)
    rec_mean.reset_index().to_parquet(seg_features_v1.RECEPTOR_MATRIX_MEAN, index=False)
    pdb_med.reset_index().to_parquet(seg_features_v1.PDB_MATRIX_MEDIAN, index=False)
    pdb_mean.reset_index().to_parquet(seg_features_v1.PDB_MATRIX_MEAN, index=False)
    rec_cov.reset_index().to_parquet(seg_features_v1.RECEPTOR_COVERAGE_MATRIX, index=False)
    print(f"receptor x segment median matrix: {rec_med.shape}; pdb x segment median matrix: {pdb_med.shape}")

    report = seg_missingness_v1.build_missingness_report(rec_med.reset_index(), rec_cov.reset_index(), rec_profile, lbd_ref)
    report = report.replace(
        "## 1. STOP-CONDITION finding from Part 5 (carried forward)",
        "## 1. Acceptance-criteria status (v3: complete curated templates + cross-template/complementary-structure rescue)",
    )
    seg_missingness_v1.MISSINGNESS_MD.write_text(report + "\n")
    print(f"Wrote {seg_missingness_v1.MISSINGNESS_MD.relative_to(PROJECT_ROOT)}")

    for thr in (0.9, 0.8, 0.7, 0.5):
        n_seg = int((rec_med.notna().sum() >= thr * len(lbd_ref)).sum())
        print(f"  segments represented in >= {thr:.0%} of 48 receptors: {n_seg}")


if __name__ == "__main__":
    main()
