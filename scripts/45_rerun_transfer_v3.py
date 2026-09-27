"""
Structural-segment alignment v3, Parts 4-6: rerun the CE structural-transfer
pipeline using the fully-curated v3 templates (script 44), and audit H12
(plus H5/H10-11 briefly) to classify every coordinate-containing receptor
without a mapped H12 (Part 6).

Same mechanism as v2 (scripts/41_structural_align_and_transfer.py): align
every non-anchor receptor's reference structure against all 5 (now richer)
anchor templates with Combinatorial Extension, transfer through the actual
aligned residue pairs, refine against the target's own secondary structure,
require residue-order continuity (monotonic CE path, enforced in
scripts/utils/structural_alignment.py), preserve low-confidence assignments
explicitly.

Reads:
    data/analysis/lbd_structural_reference_structures.csv
    data/analysis/canonical_lbd_template_segments_v2.csv
    data/analysis/normalized_lbd_bfactors_all.parquet
    data/raw/mmcif/{PDB_ID}.cif.gz (already-cached)

Writes:
    data/analysis/lbd_structural_alignment_scores.csv   (OVERWRITTEN, v3)
    data/analysis/lbd_structural_label_transfer.csv     (OVERWRITTEN, v3)
    data/analysis/receptor_segment_mapping_sources.csv
    reports/tables/h12_mapping_audit.csv
    reports/tables/h5_h10_h11_mapping_audit.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.lbd_ca_track import build_ca_track  # noqa: E402
from scripts.utils.mmcif_parse import get_secondary_structure_ranges, parse_mmcif_gz  # noqa: E402
from scripts.utils.segment_topology import Run, build_runs, merge_sheet_runs  # noqa: E402
from scripts.utils.structural_alignment import align_structures  # noqa: E402

ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
RAW_MMCIF = PROJECT_ROOT / "data" / "raw" / "mmcif"

REF_CSV = ANALYSIS / "lbd_structural_reference_structures.csv"
TEMPLATE_V2_CSV = ANALYSIS / "canonical_lbd_template_segments_v2.csv"
NORMALIZED_PARQUET = ANALYSIS / "normalized_lbd_bfactors_all.parquet"

SCORES_CSV = ANALYSIS / "lbd_structural_alignment_scores.csv"
TRANSFER_CSV = ANALYSIS / "lbd_structural_label_transfer.csv"
MAPPING_SOURCES_CSV = ANALYSIS / "receptor_segment_mapping_sources.csv"
H12_AUDIT_CSV = TABLES / "h12_mapping_audit.csv"
H5_H10_H11_AUDIT_CSV = TABLES / "h5_h10_h11_mapping_audit.csv"

ANCHOR_UNIPROT_IDS = ["P11473", "P37231", "P03372", "P04150", "P22736"]
MIN_Z_SCORE_USABLE = 3.5


def get_target_runs(pdb_id: str, label_asym_id: str, residues: pd.DataFrame) -> list[Run]:
    d = parse_mmcif_gz(RAW_MMCIF / f"{pdb_id}.cif.gz")
    ranges = [r for r in get_secondary_structure_ranges(d) if r["label_asym_id"] == label_asym_id]
    pos_by_seq = dict(zip(residues["entity_label_seq_id"].astype("Int64"), residues["canonical_uniprot_position"]))
    positions, states = [], []
    for seq_id, pos in pos_by_seq.items():
        if pd.isna(seq_id):
            continue
        state = "COIL"
        for r in ranges:
            if r["beg_label_seq_id"] <= seq_id <= r["end_label_seq_id"]:
                state = "HELIX" if r["sse_type"] == "HELIX" else "SHEET"
                if state == "HELIX":
                    break
        positions.append(int(pos))
        states.append(state)
    runs = build_runs(positions, states)
    return merge_sheet_runs([r for r in runs if r.state in ("HELIX", "SHEET")])


def refine_boundary(core_positions: set[int], target_runs: list[Run], expected_class: str) -> Run | None:
    candidates = [r for r in target_runs if r.state == expected_class and (set(range(r.start, r.end + 1)) & core_positions)]
    if not candidates:
        return None
    best = max(candidates, key=lambda r: len(set(range(r.start, r.end + 1)) & core_positions))
    return Run(expected_class, best.start, best.end)


def main() -> None:
    ref = pd.read_csv(REF_CSV, dtype={"pdb_id": str}).set_index("uniprot_id")
    templates = pd.read_csv(TEMPLATE_V2_CSV, dtype={"reference_pdb": str}, keep_default_na=False)
    norm = pd.read_parquet(
        NORMALIZED_PARQUET,
        columns=["instance_id", "uniprot_id", "pdb_id", "label_asym_id", "canonical_uniprot_position",
                 "canonical_residue_name", "entity_label_seq_id", "is_construct_mapped", "ca_positive_occupancy"],
    )

    def observed(instance_id: str) -> pd.DataFrame:
        return norm[(norm["instance_id"] == instance_id) & norm["is_construct_mapped"] & norm["ca_positive_occupancy"]]

    print("Extracting Ca tracks for all 42 reference structures...")
    tracks: dict[str, tuple[list, list]] = {}
    for uid, r in ref.iterrows():
        residues = observed(r["instance_id"])
        tracks[uid] = build_ca_track(r["pdb_id"], r["label_asym_id"], residues, RAW_MMCIF)
    for uid, (c, p) in tracks.items():
        assert len(c) == len(p) and len(p) == len(set(p))

    print("Aligning every receptor against every anchor template (Combinatorial Extension)...")
    names_lookup = norm.drop_duplicates(["uniprot_id", "canonical_uniprot_position"]).set_index(
        ["uniprot_id", "canonical_uniprot_position"]
    )["canonical_residue_name"]
    score_rows = []
    alignments: dict[tuple[str, str], object] = {}
    for tmpl_uid in ANCHOR_UNIPROT_IDS:
        tmpl_coords, tmpl_positions = tracks[tmpl_uid]
        for target_uid, (coords, positions) in tracks.items():
            if target_uid == tmpl_uid:
                continue
            result = align_structures(tmpl_coords, coords)
            if result is None:
                score_rows.append({"template_uniprot_id": tmpl_uid, "target_uniprot_id": target_uid,
                                    "n_aligned": 0, "z_score": None, "rmsd": None, "coverage": 0.0,
                                    "sequence_identity": None, "usable": False})
                continue
            n_ident = 0
            for ia, ib in zip(result.idx_a, result.idx_b):
                pa, pb = tmpl_positions[ia], positions[ib]
                na = names_lookup.get((tmpl_uid, pa), "")
                nb = names_lookup.get((target_uid, pb), "")
                if na and na == nb:
                    n_ident += 1
            coverage = result.n_aligned / len(positions)
            score_rows.append({
                "template_uniprot_id": tmpl_uid, "target_uniprot_id": target_uid,
                "n_aligned": result.n_aligned, "z_score": result.z_score, "rmsd": result.rmsd,
                "coverage": coverage, "sequence_identity": n_ident / result.n_aligned,
                "usable": result.z_score >= MIN_Z_SCORE_USABLE,
            })
            alignments[(tmpl_uid, target_uid)] = result

    scores = pd.DataFrame(score_rows).sort_values(["target_uniprot_id", "template_uniprot_id"])
    scores.to_csv(SCORES_CSV, index=False)
    print(f"Wrote {SCORES_CSV.relative_to(PROJECT_ROOT)} ({len(scores)} rows); usable: {int(scores['usable'].sum())}/{len(scores)}")

    print("\nTransferring segment labels (all 5 templates contribute per receptor, not just the single best)...")
    named_templates = templates[templates["structural_segment"] != ""]
    transfer_rows = []

    for uid in tracks:
        if uid in ANCHOR_UNIPROT_IDS:
            own = named_templates[named_templates["uniprot_id"] == uid]
            for _, row in own.iterrows():
                transfer_rows.append({
                    "uniprot_id": uid, "canonical_uniprot_position": row["canonical_uniprot_position"],
                    "structural_segment": row["structural_segment"], "structural_class": row["structural_class"],
                    "reference_template": uid, "reference_pdb": row["reference_pdb"],
                    "structural_alignment_score": None, "structural_alignment_coverage": None,
                    "assignment_method": "own_canonical_template",
                    "assignment_confidence": row["assignment_confidence"],
                    "evidence_notes": row["publication"],
                })
            continue

        residues = observed(ref.loc[uid, "instance_id"])
        target_runs = get_target_runs(ref.loc[uid, "pdb_id"], ref.loc[uid, "label_asym_id"], residues)

        # For EACH segment, use the best-scoring template that actually offers that segment
        # (a receptor's H5 might best come from PPARgamma while its H12 best comes from GR).
        candidates_by_target = scores[(scores["target_uniprot_id"] == uid) & scores["usable"]].sort_values(
            ["z_score", "n_aligned", "template_uniprot_id"], ascending=[False, False, True]
        )
        segment_best: dict[str, dict] = {}
        for _, srow in candidates_by_target.iterrows():
            tmpl_uid = srow["template_uniprot_id"]
            tmpl_named = named_templates[named_templates["uniprot_id"] == tmpl_uid]
            for seg in tmpl_named["structural_segment"].unique():
                if seg not in segment_best:
                    segment_best[seg] = {"template_uniprot_id": tmpl_uid, "score_row": srow}

        for seg, info in segment_best.items():
            tmpl_uid = info["template_uniprot_id"]
            srow = info["score_row"]
            result = alignments[(tmpl_uid, uid)]
            tmpl_coords, tmpl_positions = tracks[tmpl_uid]
            _, positions = tracks[uid]
            pos_map = {tmpl_positions[ia]: positions[ib] for ia, ib in zip(result.idx_a, result.idx_b)}
            g = named_templates[(named_templates["uniprot_id"] == tmpl_uid) & (named_templates["structural_segment"] == seg)]

            core_target_positions = {pos_map[p] for p in g["canonical_uniprot_position"] if p in pos_map}
            if not core_target_positions:
                continue
            expected_class = g["structural_class"].iloc[0]
            refined = refine_boundary(core_target_positions, target_runs, expected_class)
            if refined is None:
                # No overlapping target secondary-structure run of the expected class was
                # found -- never fall back to a raw min/max span of the aligned core, since
                # that can bridge genuine evidence gaps (unobserved/coil positions) and label
                # them HELIX/SHEET with no underlying secondary-structure evidence at all
                # (the exact class of bug already fixed once in refine_boundary itself).
                continue
            span = refined
            n_core = len(core_target_positions)
            n_template_seg = len(g)
            core_fraction = n_core / n_template_seg
            tmpl_conf = g["assignment_confidence"].iloc[0]
            if tmpl_conf == "HIGH" and core_fraction >= 0.5:
                confidence, status = "HIGH", "ASSIGNED_HIGH_CONFIDENCE"
            elif core_fraction >= 0.3:
                confidence, status = "MODERATE", "ASSIGNED_MODERATE_CONFIDENCE"
            else:
                confidence, status = "LOW", "AMBIGUOUS"
            note = (
                f"Transferred from {tmpl_uid} template segment {seg} via CE structural alignment "
                f"(z={srow['z_score']:.2f}, {result.n_aligned} aligned residues, "
                f"{core_fraction:.0%} of template segment's aligned core, template confidence={tmpl_conf}); "
                f"boundary refined to the full overlapping target secondary-structure run."
            )
            for pos in range(span.start, span.end + 1):
                transfer_rows.append({
                    "uniprot_id": uid, "canonical_uniprot_position": pos, "structural_segment": seg,
                    "structural_class": expected_class, "reference_template": tmpl_uid,
                    "reference_pdb": ref.loc[uid, "pdb_id"],
                    "structural_alignment_score": srow["z_score"], "structural_alignment_coverage": srow["coverage"],
                    "assignment_method": "structural_alignment_transfer",
                    "assignment_confidence": confidence, "assignment_status_hint": status,
                    "evidence_notes": note,
                })

    transfer = pd.DataFrame(transfer_rows)
    conf_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2}
    transfer["_rank"] = transfer.get("assignment_confidence", "LOW").map(conf_rank).fillna(2)
    # Resolve overlaps between DIFFERENT segments' spans (can happen when two segments'
    # boundary-refined runs collide): higher confidence wins, deterministic tie-break by
    # canonical order (earlier-numbered segment keeps the position).
    transfer = transfer.sort_values(["uniprot_id", "canonical_uniprot_position", "_rank"]).drop_duplicates(
        ["uniprot_id", "canonical_uniprot_position"], keep="first"
    ).drop(columns="_rank")

    # ---- Part 6 rescue: for receptors still missing H12 (or H5/H10/H11), try their OTHER
    # coordinate candidates (not just the single chosen reference) for one whose own
    # secondary structure actually resolves that region -- "use complementary structures
    # to resolve category B whenever possible". ----
    print("\nComplementary-structure rescue for missing H12/H5/H10/H11 (Part 6)...")
    master = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "structures.csv", dtype={"pdb_id": str})
    rescued_rows = []
    for seg_name in ["H12", "H5", "H10", "H11"]:
        have_seg = set(transfer.loc[transfer["structural_segment"] == seg_name, "uniprot_id"]) | set(
            named_templates.loc[named_templates["structural_segment"] == seg_name, "uniprot_id"]
        )
        missing = [u for u in tracks if u not in have_seg and u not in ANCHOR_UNIPROT_IDS]
        for uid in missing:
            best_tmpl_rows = scores[(scores["target_uniprot_id"] == uid) & scores["usable"]].sort_values(
                ["z_score", "n_aligned", "template_uniprot_id"], ascending=[False, False, True]
            )
            tmpl_uid = next(
                (r["template_uniprot_id"] for _, r in best_tmpl_rows.iterrows()
                 if named_templates[(named_templates["uniprot_id"] == r["template_uniprot_id"]) & (named_templates["structural_segment"] == seg_name)].shape[0] > 0),
                None,
            )
            if tmpl_uid is None:
                continue
            g = named_templates[(named_templates["uniprot_id"] == tmpl_uid) & (named_templates["structural_segment"] == seg_name)]
            tmpl_coords, tmpl_positions = tracks[tmpl_uid]

            candidates = master[(master["uniprot_id"] == uid) & master["has_stage3b_coordinate_data"]].sort_values(
                "positive_occupancy_lbd_coverage", ascending=False
            ).head(10)
            for _, cand in candidates.iterrows():
                if cand["representative_instance_id"] == ref.loc[uid, "instance_id"]:
                    continue  # already tried as the primary reference
                cand_residues = observed(cand["representative_instance_id"])
                if len(cand_residues) == 0:
                    continue
                cand_coords, cand_positions = build_ca_track(cand["pdb_id"], cand["label_asym_id"], cand_residues, RAW_MMCIF)
                result = align_structures(tmpl_coords, cand_coords)
                if result is None or result.z_score < MIN_Z_SCORE_USABLE:
                    continue
                pos_map = {tmpl_positions[ia]: cand_positions[ib] for ia, ib in zip(result.idx_a, result.idx_b)}
                core_target = {pos_map[p] for p in g["canonical_uniprot_position"] if p in pos_map}
                if len(core_target) < max(3, 0.3 * len(g)):
                    continue
                cand_runs = get_target_runs(cand["pdb_id"], cand["label_asym_id"], cand_residues)
                expected_class = g["structural_class"].iloc[0]
                refined = refine_boundary(core_target, cand_runs, expected_class)
                if refined is None:
                    continue
                core_fraction = len(core_target) / len(g)
                confidence = "HIGH" if core_fraction >= 0.5 else "MODERATE"
                for pos in range(refined.start, refined.end + 1):
                    rescued_rows.append({
                        "uniprot_id": uid, "canonical_uniprot_position": pos, "structural_segment": seg_name,
                        "structural_class": expected_class, "reference_template": tmpl_uid,
                        "reference_pdb": cand["pdb_id"],
                        "structural_alignment_score": result.z_score, "structural_alignment_coverage": result.n_aligned / len(cand_positions),
                        "assignment_method": "complementary_structure_rescue",
                        "assignment_confidence": confidence,
                        "assignment_status_hint": "ASSIGNED_HIGH_CONFIDENCE" if confidence == "HIGH" else "ASSIGNED_MODERATE_CONFIDENCE",
                        "evidence_notes": f"H12/H5/H10/H11 rescue: {seg_name} not resolved in {uid}'s primary reference "
                                          f"({ref.loc[uid, 'pdb_id']}) but recovered from complementary coordinate candidate "
                                          f"{cand['pdb_id']} via CE alignment to {tmpl_uid} template (z={result.z_score:.2f}, "
                                          f"{core_fraction:.0%} of template core aligned), boundary refined to the "
                                          f"complementary structure's own secondary-structure run.",
                    })
                break  # one successful rescue per receptor per segment is enough

    if rescued_rows:
        rescued = pd.DataFrame(rescued_rows)
        rescued["_rank"] = 0
        transfer = pd.concat([transfer, rescued], ignore_index=True).sort_values(
            ["uniprot_id", "canonical_uniprot_position", "_rank"]
        ).drop_duplicates(["uniprot_id", "canonical_uniprot_position"], keep="first").drop(columns="_rank")
        print(f"  rescued {rescued['uniprot_id'].nunique()} receptor x segment gaps from complementary structures")

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    transfer.to_csv(TRANSFER_CSV, index=False)
    print(f"Wrote {TRANSFER_CSV.relative_to(PROJECT_ROOT)} ({len(transfer)} rows, {transfer['uniprot_id'].nunique()} receptors)")
    print(transfer.groupby("structural_segment")["uniprot_id"].nunique().sort_values(ascending=False).to_string())

    # ---- Part 4: receptor_segment_mapping_sources.csv ----
    print("\nBuilding receptor_segment_mapping_sources.csv...")
    src_rows = []
    for uid in tracks:
        if uid in ANCHOR_UNIPROT_IDS:
            own = named_templates[named_templates["uniprot_id"] == uid]
            for seg, g in own.groupby("structural_segment"):
                src_rows.append({"uniprot_id": uid, "structural_segment": seg, "source_pdb": g["reference_pdb"].iloc[0],
                                  "source_method": "own_canonical_template"})
        else:
            t = transfer[transfer["uniprot_id"] == uid]
            for seg, g in t.groupby("structural_segment"):
                src_rows.append({"uniprot_id": uid, "structural_segment": seg, "source_pdb": ref.loc[uid, "pdb_id"],
                                  "source_method": f"structural_alignment_transfer <- {g['reference_template'].iloc[0]}"})
    pd.DataFrame(src_rows).to_csv(MAPPING_SOURCES_CSV, index=False)
    print(f"Wrote {MAPPING_SOURCES_CSV.relative_to(PROJECT_ROOT)} ({len(src_rows)} rows)")

    # ---- Part 6: H12 (and H5/H10/H11) audit -- run AFTER the complementary-structure
    # rescue above, so a receptor only appears here if rescue genuinely could not resolve it ----
    print("\nBuilding H12 mapping audit...")

    def audit_segment(seg_name: str) -> pd.DataFrame:
        rows = []
        for uid in tracks:
            has_seg = ((uid in ANCHOR_UNIPROT_IDS and (named_templates[(named_templates["uniprot_id"] == uid) & (named_templates["structural_segment"] == seg_name)].shape[0] > 0))
                       or (uid not in ANCHOR_UNIPROT_IDS and (transfer[(transfer["uniprot_id"] == uid) & (transfer["structural_segment"] == seg_name)].shape[0] > 0)))
            if has_seg:
                continue  # audit only receptors WITHOUT the segment mapped, post-rescue
            n_candidates = int((master["uniprot_id"] == uid).sum())
            candidate_pdbs = master[(master["uniprot_id"] == uid) & master["has_stage3b_coordinate_data"]].sort_values(
                "positive_occupancy_lbd_coverage", ascending=False
            ).head(10)
            any_candidate_has_helix_evidence = False
            for _, cand in candidate_pdbs.iterrows():
                cres = observed(cand["representative_instance_id"])
                if len(cres) == 0:
                    continue
                if any(r.state == "HELIX" for r in get_target_runs(cand["pdb_id"], cand["label_asym_id"], cres)):
                    any_candidate_has_helix_evidence = True
                    break
            if uid in ANCHOR_UNIPROT_IDS:
                category, note = "D", f"{uid} is an anchor: {seg_name} not curated in the v3 template (no literature landmark + cross-template rescue failed)."
            elif not any_candidate_has_helix_evidence:
                category, note = "A", f"No helix secondary-structure evidence at all in any of this receptor's top-{len(candidate_pdbs)} coordinate candidates."
            else:
                best_row = scores[(scores["target_uniprot_id"] == uid)].sort_values("z_score", ascending=False).iloc[0] if len(scores[scores["target_uniprot_id"] == uid]) else None
                tmpl_has_seg = any(
                    named_templates[(named_templates["uniprot_id"] == t) & (named_templates["structural_segment"] == seg_name)].shape[0] > 0
                    for t in ANCHOR_UNIPROT_IDS
                )
                if best_row is None or not best_row["usable"]:
                    category, note = "C", f"No usable (z>=3.5) structural alignment to any anchor template (best z={best_row['z_score'] if best_row is not None else 'n/a'})."
                elif not tmpl_has_seg:
                    category, note = "D", f"No anchor template curates {seg_name} at all; nothing to transfer."
                else:
                    category, note = "C", (
                        f"Helix evidence exists in this receptor's reference or a complementary candidate, and a usable "
                        f"alignment (best z={best_row['z_score']:.2f}) and template both exist, but the transferred "
                        f"core/refinement/complementary-rescue attempt did not survive the confidence threshold."
                    )
            rows.append({
                "uniprot_id": uid, "segment": seg_name, "category": category,
                "n_stage2_candidates": n_candidates, "n_coordinate_candidates_checked": len(candidate_pdbs), "note": note,
            })
        return pd.DataFrame(rows)

    h12 = audit_segment("H12")
    h12.to_csv(H12_AUDIT_CSV, index=False)
    print(f"Wrote {H12_AUDIT_CSV.relative_to(PROJECT_ROOT)} ({len(h12)} rows)")
    print(h12["category"].value_counts().to_string() if len(h12) else "  (all receptors have H12 mapped)")

    other = pd.concat([audit_segment("H5"), audit_segment("H10"), audit_segment("H11")], ignore_index=True)
    other.to_csv(H5_H10_H11_AUDIT_CSV, index=False)
    print(f"Wrote {H5_H10_H11_AUDIT_CSV.relative_to(PROJECT_ROOT)} ({len(other)} rows)")


if __name__ == "__main__":
    main()
