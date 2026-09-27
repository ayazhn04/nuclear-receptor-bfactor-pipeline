"""
Structural-segment alignment v2, Parts 4-7: real 3D structural alignment
(Combinatorial Extension -- see scripts/utils/structural_alignment.py for
why no new software was installed) of every receptor's reference LBD
structure against all 5 anchor templates, selection of the best-supported
template per receptor, transfer of canonical segment labels through the
ACTUAL aligned residue-residue correspondence (never through raw UniProt
integer equality, never through helix ordinal position), and boundary
refinement using the target's own secondary-structure evidence.

Reads:
    data/analysis/lbd_structural_reference_structures.csv
    data/analysis/canonical_lbd_template_segments.csv
    data/analysis/normalized_lbd_bfactors_all.parquet (identity/mapping columns)
    data/raw/mmcif/{PDB_ID}.cif.gz (already-cached)

Writes:
    data/analysis/lbd_structural_alignment_scores.csv
    data/analysis/lbd_structural_label_transfer.csv  (intermediate: per
        target-receptor transferred + refined segment calls, before being
        assembled into the final per-position segment map in script 42)
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
RAW_MMCIF = PROJECT_ROOT / "data" / "raw" / "mmcif"

REF_CSV = ANALYSIS / "lbd_structural_reference_structures.csv"
TEMPLATE_CSV = ANALYSIS / "canonical_lbd_template_segments.csv"
NORMALIZED_PARQUET = ANALYSIS / "normalized_lbd_bfactors_all.parquet"

SCORES_CSV = ANALYSIS / "lbd_structural_alignment_scores.csv"
TRANSFER_CSV = ANALYSIS / "lbd_structural_label_transfer.csv"

ANCHOR_UNIPROT_IDS = ["P11473", "P37231", "P03372", "P04150", "P22736"]  # ordered list, not a set: iteration order must be deterministic across process runs (string set order is hash-seed-randomized)
MIN_Z_SCORE_USABLE = 3.5  # CE's own significance threshold (matches Biopython's optimizer trigger)


def canonical_order_key(seg: str) -> tuple:
    order = ["H1", "H2", "H3", "H4", "H5", "beta_sheet", "H6", "H7", "H8", "H9", "H10", "H11", "H12"]
    return (order.index(seg),) if seg in order else (99,)


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
    """Extend a transferred segment's core aligned positions to the full
    extent of the target's own overlapping secondary-structure run (Part 6:
    3D correspondence establishes identity; the target's own DSSP/mmCIF
    evidence then refines the exact boundary).

    If the aligned core touches more than one DISJOINT target run (e.g. a
    real gap of unobserved/zero-occupancy residues splits what the template
    calls one segment), never bridge across that gap -- that would silently
    label unobserved residues as part of the helix. Pick the single run with
    the greatest overlap with the aligned core instead."""
    candidates = [r for r in target_runs if r.state == expected_class and (set(range(r.start, r.end + 1)) & core_positions)]
    if not candidates:
        return None
    best = max(candidates, key=lambda r: len(set(range(r.start, r.end + 1)) & core_positions))
    return Run(expected_class, best.start, best.end)


def main() -> None:
    ref = pd.read_csv(REF_CSV, dtype={"pdb_id": str}).set_index("uniprot_id")
    templates = pd.read_csv(TEMPLATE_CSV, dtype={"pdb_id": str})
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

    scores = pd.DataFrame(score_rows).sort_values(["target_uniprot_id", "z_score"], ascending=[True, False])
    scores.to_csv(SCORES_CSV, index=False)
    print(f"Wrote {SCORES_CSV.relative_to(PROJECT_ROOT)} ({len(scores)} rows)")
    print(f"  usable (z >= {MIN_Z_SCORE_USABLE}) alignments: {int(scores['usable'].sum())} / {len(scores)}")

    print("\nSelecting best-supported template per receptor and transferring labels...")
    named_templates = templates[templates["canonical_segment"] != ""]
    transfer_rows = []

    for uid in tracks:
        if uid in set(ANCHOR_UNIPROT_IDS):
            own = named_templates[named_templates["uniprot_id"] == uid]
            for _, row in own.iterrows():
                transfer_rows.append({
                    "uniprot_id": uid, "canonical_uniprot_position": row["canonical_uniprot_position"],
                    "structural_segment": row["canonical_segment"], "structural_class": row["structural_class"],
                    "reference_template": uid, "reference_pdb": row["pdb_id"],
                    "structural_alignment_score": None, "structural_alignment_coverage": None,
                    "assignment_method": "own_canonical_template",
                    "assignment_confidence": row["assignment_confidence"],
                    "evidence_notes": row["evidence_source"],
                })
            continue

        candidates = scores[(scores["target_uniprot_id"] == uid) & scores["usable"]].sort_values(
            ["z_score", "n_aligned", "template_uniprot_id"], ascending=[False, False, True]
        )
        if len(candidates) == 0:
            continue
        best = candidates.iloc[0]
        tmpl_uid = best["template_uniprot_id"]
        result = alignments[(tmpl_uid, uid)]
        tmpl_coords, tmpl_positions = tracks[tmpl_uid]
        _, positions = tracks[uid]
        pos_map = {tmpl_positions[ia]: positions[ib] for ia, ib in zip(result.idx_a, result.idx_b)}

        tmpl_named = named_templates[named_templates["uniprot_id"] == tmpl_uid]
        residues = observed(ref.loc[uid, "instance_id"])
        target_runs = get_target_runs(ref.loc[uid, "pdb_id"], ref.loc[uid, "label_asym_id"], residues)

        for seg, g in tmpl_named.groupby("canonical_segment"):
            core_target_positions = {pos_map[p] for p in g["canonical_uniprot_position"] if p in pos_map}
            if not core_target_positions:
                continue
            expected_class = g["structural_class"].iloc[0]
            refined = refine_boundary(core_target_positions, target_runs, expected_class)
            span = refined if refined is not None else Run(expected_class, min(core_target_positions), max(core_target_positions))
            n_core = len(core_target_positions)
            n_template_seg = len(g)
            core_fraction = n_core / n_template_seg
            tmpl_conf = g["assignment_confidence"].iloc[0]
            if tmpl_conf == "HIGH" and core_fraction >= 0.5 and refined is not None:
                confidence, status = "HIGH", "ASSIGNED_HIGH_CONFIDENCE"
            elif core_fraction >= 0.3:
                confidence, status = "MODERATE", "ASSIGNED_MODERATE_CONFIDENCE"
            else:
                confidence, status = "LOW", "AMBIGUOUS"
            note = (
                f"Transferred from {tmpl_uid} template segment {seg} via CE structural alignment "
                f"(z={best['z_score']:.2f}, {result.n_aligned} aligned residues, "
                f"{core_fraction:.0%} of template segment's aligned core, template confidence={tmpl_conf}); "
                f"boundary {'refined to the full overlapping target secondary-structure run' if refined is not None else 'NOT refined -- no overlapping target run found, using the aligned-core span only'}."
            )
            for pos in range(span.start, span.end + 1):
                transfer_rows.append({
                    "uniprot_id": uid, "canonical_uniprot_position": pos, "structural_segment": seg,
                    "structural_class": expected_class, "reference_template": tmpl_uid,
                    "reference_pdb": ref.loc[uid, "pdb_id"],
                    "structural_alignment_score": best["z_score"], "structural_alignment_coverage": best["coverage"],
                    "assignment_method": "structural_alignment_transfer",
                    "assignment_confidence": confidence, "assignment_status_hint": status,
                    "evidence_notes": note,
                })

    transfer = pd.DataFrame(transfer_rows)
    # Later (higher-confidence-first) rows win a position if two transferred
    # segments' refined spans happen to overlap (rare, but must never silently
    # produce >1 label for one position).
    conf_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2}
    transfer["_rank"] = transfer.get("assignment_confidence", "LOW").map(conf_rank).fillna(2)
    transfer = transfer.sort_values(["uniprot_id", "canonical_uniprot_position", "_rank"]).drop_duplicates(
        ["uniprot_id", "canonical_uniprot_position"], keep="first"
    ).drop(columns="_rank")

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    transfer.to_csv(TRANSFER_CSV, index=False)
    print(f"Wrote {TRANSFER_CSV.relative_to(PROJECT_ROOT)} ({len(transfer)} rows, {transfer['uniprot_id'].nunique()} receptors)")
    print(transfer.groupby("structural_segment")["uniprot_id"].nunique().sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
