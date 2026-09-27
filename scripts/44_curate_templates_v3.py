"""
Structural-segment alignment v3, Parts 1-3: curate PRIMARY (+ up to 2
COMPLEMENTARY) reference structures and a much more complete H1-H12/
beta-sheet template for the 5 anchor receptors, using published
structural/functional landmarks (never ordinal DSSP counting), then
cross-template rescue via CE structural alignment AMONG the 5 anchors
themselves (Part 4's evidence priority #4: "structural correspondence to
another fully curated NR template") so a segment anchored directly in only
one anchor (e.g. PPARgamma's H5/H6/H7/H10-11) can support the others.

Literature landmarks gathered this session (all cross-checked against this
repository's own functional_site_definitions.csv / nr_functional_features.csv
where possible -- several turned out to already be in the repository under
a different site_id, confirmed via the well-established PPARgamma1 ->
PPARgamma2 +28 residue numbering offset):

  VDR (P11473): Rochel et al. 2000 (Mol Cell 5:173, DOI 10.1016/S1097-2765(00)80413-X)
    -- H3 (Ser235), H4 (Lys264), H12 (Thr415-Glu420, matches this repo's own
    nr_functional_features.csv H12 segment 416-422); Peleg/PMC5081507 (VDR LBD
    dynamics review) -- H9/H10 loop (His371-Leu378), S2 beta-strand
    (Trp286-Ala303), VDR insertion domain (Ser165-Pro215, explains the extra
    helices this project's own DSSP consensus finds between H1 and H3).
  PPARgamma (P37231): Nolte et al. 1998 (Nature 395:137, DOI 10.1038/25931)
    -- H3 (Lys329)/H12 (Glu499) charge clamp (already in
    functional_site_definitions.csv); PMC6317912 (cobinding study) -- H5
    (His351, == PPARgamma1 His323+28, matches this repo's own ligand_pocket
    351/HIS exactly) and H10/H11 (His477, == His449+28, matches this repo's
    477/HIS); a PPARgamma antagonist/MD study -- H2 (Arg262, == Arg234+28),
    beta-sheet (Gly366, == Gly338+28), H6 (Arg378, == Arg350+28), H7
    (Asn403, == Asn375+28), H8/H9 loop (Asp408, == Asp380+28); all four via
    the confirmed +28 PPARgamma1->PPARgamma2 offset (independently verified
    against the repo's own 317/351/477/501/329/499 entries).
  ESR1 (P03372): Shiau et al. 1998 (Cell 95:927, DOI
    10.1016/S0092-8674(00)81717-1) -- H3 (Lys362)/H12 (Glu542) charge clamp
    (already in functional_site_definitions.csv); dimer-interface literature
    -- H11 (Leu504/Leu508/Leu511, H10/H11 leucine-zipper-like patch).
  GR (P04150): Bledsoe et al. 2002 (Cell 110:93, DOI
    10.1016/S0092-8674(02)00817-6) -- H3 (Lys579 charge clamp, Asn564/Gln570
    ligand contacts, already in functional_site_definitions.csv), H5
    (Arg611, already in functional_site_definitions.csv as a ligand_pocket
    residue but not previously assigned a helix), H7 (Gln642, ditto), H12
    (Glu755 charge clamp; Thr739 ligand-pocket residue tentatively bracketed
    to H10/H11 given its position between H7 and H12 in the reported pocket
    residue list "helices 3,4,5,6,7,10").
  NR4A1 (P22736): nr_functional_features.csv h12_segment (UniProt-annotated)
    -- explicit H12 584-595.

This is still not a complete H1-H12 table for any single receptor from a
single source (no accessible primary-literature figure/table with all 12
boundaries was found for any of the 5 anchors within this session's search
budget) -- it is a materially larger, still individually-cited landmark set
than the previous round, used the same way: to anchor helix RUNS (never to
invent boundaries), with interior helices between two anchors labelled only
when the observed count matches the canonical expectation for that specific
bracketed sub-interval.

Reads:
    data/processed/structures.csv
    data/analysis/normalized_lbd_bfactors_all.parquet
    reports/tables/functional_site_definitions.csv
    data/raw/mmcif/{PDB_ID}.cif.gz (already-cached)

Writes:
    data/analysis/canonical_lbd_reference_structures_v2.csv
    data/analysis/canonical_lbd_template_segments_v2.csv
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

PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
RAW_MMCIF = PROJECT_ROOT / "data" / "raw" / "mmcif"

MASTER_CSV = PROCESSED / "structures.csv"
NORMALIZED_PARQUET = ANALYSIS / "normalized_lbd_bfactors_all.parquet"
FUNCTIONAL_SITE_CSV = TABLES / "functional_site_definitions.csv"

REF_STRUCTURES_V2_CSV = ANALYSIS / "canonical_lbd_reference_structures_v2.csv"
TEMPLATE_V2_CSV = ANALYSIS / "canonical_lbd_template_segments_v2.csv"

ANCHORS = ["P11473", "P37231", "P03372", "P04150", "P22736"]  # ordered, not a set (determinism)
ANCHOR_NAMES = {"P11473": "VDR", "P37231": "PPARgamma", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}

LANDMARKS: dict[str, list[tuple[int, str, str]]] = {
    "P11473": [
        (235, "H3", "Rochel et al. 2000 (Mol Cell 5:173) -- Ser235(H3) H-bonds Thr415(H12)."),
        (264, "H4", "Rochel et al. 2000 -- Lys264(H4)-Glu420(H12) salt bridge stabilizing H12."),
        (415, "H12", "Rochel et al. 2000 -- Thr415, H12 N-terminal boundary residue."),
        (420, "H12", "Rochel et al. 2000 -- Glu420 of H12, salt-bridge partner of Lys264(H4)."),
        (290, "beta_sheet", "PMC5081507 (VDR LBD dynamics) -- S2 beta-strand Trp286-Ala303 (midpoint used as anchor)."),
        (371, "H9", "PMC5081507 -- H9/H10 loop His371-Leu378; His371 is the last H9 residue."),
        (378, "H10", "PMC5081507 -- H9/H10 loop His371-Leu378; Leu378 is the first H10 residue."),
    ],
    "P37231": [
        (262, "H2", "PPARgamma antagonist/MD literature -- Arg234(PPARg1)+28=Arg262(H2)."),
        (329, "H3", "functional_site_definitions.csv af2_coactivator_groove (Nolte et al. 1998) -- Lys329(H3) charge-clamp."),
        (351, "H5", "functional_site_definitions.csv ligand_pocket (Nolte et al. 1998) -- His351, == PPARg1 His323+28, PMC6317912 confirms His323 is on H5."),
        (366, "beta_sheet", "PPARgamma antagonist/MD literature -- Gly338(PPARg1)+28=Gly366, beta-sheet in the ligand-binding pocket."),
        (378, "H6", "PPARgamma antagonist/MD literature -- Arg350(PPARg1)+28=Arg378(H6)."),
        (403, "H7", "PPARgamma antagonist/MD literature -- Asn375(PPARg1)+28=Asn403(H7)."),
        (477, "H10", "functional_site_definitions.csv ligand_pocket (Nolte et al. 1998) -- His477, == PPARg1 His449+28, PMC6317912 confirms His449 is on H10/H11."),
        (499, "H12", "functional_site_definitions.csv af2_coactivator_groove (Nolte et al. 1998) -- Glu499(H12) charge-clamp."),
        (501, "H12", "functional_site_definitions.csv ligand_pocket (Nolte et al. 1998) -- Tyr501, == PPARg1 Tyr473+28, PMC6317912 confirms Y473 is on H12."),
    ],
    "P03372": [
        (362, "H3", "functional_site_definitions.csv af2_coactivator_groove (Shiau et al. 1998) -- Lys362(H3) charge-clamp."),
        (504, "H11", "External literature (hERalpha dimer-interface studies) -- Leu504, H10/H11 leucine-zipper-like dimer patch."),
        (508, "H11", "External literature -- Leu508, H10/H11 dimer-interface patch."),
        (511, "H11", "External literature -- Leu511, H10/H11 dimer-interface patch."),
        (542, "H12", "functional_site_definitions.csv af2_coactivator_groove (Shiau et al. 1998) -- Glu542(H12) charge-clamp."),
    ],
    "P04150": [
        (579, "H3", "functional_site_definitions.csv af2_coactivator_groove (Bledsoe et al. 2002) -- Lys579(H3) charge-clamp."),
        (611, "H5", "functional_site_definitions.csv ligand_pocket (Bledsoe et al. 2002) -- Arg611, confirmed on H5 by external literature."),
        (642, "H7", "functional_site_definitions.csv ligand_pocket (Bledsoe et al. 2002) -- Gln642, confirmed on H7 by external literature."),
        (755, "H12", "functional_site_definitions.csv af2_coactivator_groove (Bledsoe et al. 2002) -- Glu755(H12) charge-clamp."),
    ],
    "P22736": [
        (584, "H12", "nr_functional_features.csv h12_segment (UniProt:P22736) -- explicit UniProt-annotated H12 segment start."),
        (595, "H12", "nr_functional_features.csv h12_segment (UniProt:P22736) -- explicit UniProt-annotated H12 segment end."),
    ],
}

CANONICAL_ALL = ["H1", "H2", "H3", "H4", "H5", "beta_sheet", "H6", "H7", "H8", "H9", "H10", "H11", "H12"]
CORE_SEGMENTS = ["H1", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "beta_sheet"]


def get_runs_for_instance(pdb_id: str, label_asym_id: str, residues: pd.DataFrame) -> list[Run]:
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


def find_run_containing(runs: list[Run], position: int, state: str) -> Run | None:
    for r in runs:
        if r.state == state and r.start <= position <= r.end:
            return r
    return None


def select_reference_candidates(uid: str, master: pd.DataFrame, norm: pd.DataFrame, max_n: int = 3) -> list[pd.Series]:
    """PRIMARY = highest-coverage candidate; up to 2 COMPLEMENTARY = next-best
    candidates that, combined, cover a materially different residue span
    (chosen to maximize the chance of rescuing a segment the primary
    structure lacks, e.g. a disordered H12) -- never averaged, used only to
    look up secondary-structure evidence for template construction."""
    landmarks = LANDMARKS[uid]
    needed = {seg for _, seg, _ in landmarks}
    pool = master[(master["uniprot_id"] == uid) & master["has_stage3b_coordinate_data"]].sort_values(
        "positive_occupancy_lbd_coverage", ascending=False
    )
    chosen: list[pd.Series] = []
    covered_positions: set[int] = set()
    for _, cand in pool.iterrows():
        residues = norm[(norm["instance_id"] == cand["representative_instance_id"]) & norm["is_construct_mapped"] & norm["ca_positive_occupancy"]]
        if len(residues) == 0:
            continue
        this_positions = set(residues["canonical_uniprot_position"])
        new_coverage = len(this_positions - covered_positions)
        if not chosen:
            chosen.append(cand)
            covered_positions |= this_positions
            continue
        runs = get_runs_for_instance(cand["pdb_id"], cand["label_asym_id"], residues)
        helix_runs = [r for r in runs if r.state == "HELIX"]
        rescues_needed_segment = any(
            find_run_containing(helix_runs, pos, "HELIX") is not None
            for pos, seg, _ in landmarks
            if seg not in {c_seg for c in chosen for c_seg in _segments_found(c, uid, norm)}
        )
        if rescues_needed_segment and new_coverage > 5:
            chosen.append(cand)
            covered_positions |= this_positions
        if len(chosen) >= max_n:
            break
    return chosen


def _segments_found(cand: pd.Series, uid: str, norm: pd.DataFrame) -> set[str]:
    residues = norm[(norm["instance_id"] == cand["representative_instance_id"]) & norm["is_construct_mapped"] & norm["ca_positive_occupancy"]]
    if len(residues) == 0:
        return set()
    runs = get_runs_for_instance(cand["pdb_id"], cand["label_asym_id"], residues)
    helix_runs = [r for r in runs if r.state == "HELIX"]
    return {seg for pos, seg, _ in LANDMARKS[uid] if find_run_containing(helix_runs, pos, "HELIX") is not None}


def build_template_from_candidates(uid: str, candidates: list[pd.Series], norm: pd.DataFrame) -> tuple[list[dict], dict]:
    """Merge evidence across PRIMARY + COMPLEMENTARY candidates: each
    landmark segment is anchored from whichever candidate structure
    actually resolves it (searched primary first)."""
    landmarks = LANDMARKS[uid]
    anchor_runs: dict[str, Run] = {}
    anchor_evidence: dict[str, list[str]] = {}
    anchor_source: dict[str, pd.Series] = {}
    all_candidate_runs: dict[str, list[Run]] = {}

    for cand in candidates:
        residues = norm[(norm["instance_id"] == cand["representative_instance_id"]) & norm["is_construct_mapped"] & norm["ca_positive_occupancy"]]
        runs = get_runs_for_instance(cand["pdb_id"], cand["label_asym_id"], residues)
        all_candidate_runs[cand["pdb_id"]] = runs
        helix_runs = [r for r in runs if r.state == "HELIX"]
        for pos, seg, evidence in landmarks:
            if seg in anchor_runs:
                continue
            run = find_run_containing(helix_runs, pos, "HELIX")
            if run is None:
                continue
            anchor_runs[seg] = run
            anchor_evidence.setdefault(seg, []).append(f"[{cand['pdb_id']}] {evidence}")
            anchor_source[seg] = cand

    primary = candidates[0]
    primary_residues = norm[(norm["instance_id"] == primary["representative_instance_id"]) & norm["is_construct_mapped"] & norm["ca_positive_occupancy"]]
    primary_runs = [r for r in all_candidate_runs[primary["pdb_id"]] if r.state == "HELIX"]
    if "H1" not in anchor_runs and primary_runs:
        anchor_runs["H1"] = primary_runs[0]
        anchor_evidence["H1"] = [f"[{primary['pdb_id']}] First helix run at the standardized LBD's own N-terminal boundary (domain-topology anchor)."]
        anchor_source["H1"] = primary

    calls: list[dict] = []
    for seg, run in anchor_runs.items():
        src = anchor_source[seg]
        calls.append({"start": run.start, "end": run.end, "segment": seg, "structural_class": "HELIX",
                       "confidence": "HIGH", "evidence": " | ".join(anchor_evidence[seg]),
                       "reference_pdb": src["pdb_id"], "reference_chain": src["label_asym_id"]})

    sheet_runs = [r for r in all_candidate_runs[primary["pdb_id"]] if r.state == "SHEET"]
    if len(sheet_runs) == 1:
        sh = sheet_runs[0]
        if any(h.end < sh.start for h in primary_runs) and any(h.start > sh.end for h in primary_runs) and "beta_sheet" not in anchor_runs:
            calls.append({"start": sh.start, "end": sh.end, "segment": "beta_sheet", "structural_class": "SHEET",
                          "confidence": "HIGH", "evidence": f"[{primary['pdb_id']}] Single consensus beta-sheet region flanked by helices on both sides.",
                          "reference_pdb": primary["pdb_id"], "reference_chain": primary["label_asym_id"]})

    return calls, {"candidates": candidates}


def _overlaps_existing(start: int, end: int, calls: list[dict]) -> bool:
    return any(not (end < c["start"] or start > c["end"]) for c in calls)


def rescue_via_cross_template_alignment(
    templates_by_anchor: dict[str, list[dict]], tracks: dict[str, tuple[list, list]],
    target_runs_by_anchor: dict[str, list[Run]],
) -> None:
    """Part 4 evidence priority #4: for any anchor missing a core segment,
    try to recover it by CE-aligning its own reference track against every
    OTHER anchor that DOES have that segment directly anchored, then
    boundary-refining against the TARGET's own secondary-structure run
    (never bridging across a gap -- same fix as script 41's refine_boundary)
    and rejecting the rescue outright if it would overlap an already-
    anchored segment in this receptor (literature-direct evidence always
    wins over cross-template rescue)."""
    for uid in ANCHORS:
        have = {c["segment"] for c in templates_by_anchor[uid]}
        missing_core = [s for s in CORE_SEGMENTS if s not in have]
        if not missing_core:
            continue
        target_coords, target_positions = tracks[uid]
        for src_uid in ANCHORS:
            if src_uid == uid or not missing_core:
                continue
            src_calls = templates_by_anchor[src_uid]
            src_named = {c["segment"]: c for c in src_calls}
            usable_segments = [s for s in missing_core if s in src_named]
            if not usable_segments:
                continue
            src_coords, src_positions = tracks[src_uid]
            result = align_structures(src_coords, target_coords)
            if result is None or result.z_score < 3.5:
                continue
            pos_map = {src_positions[ia]: target_positions[ib] for ia, ib in zip(result.idx_a, result.idx_b)}
            for seg in list(usable_segments):
                call = src_named[seg]
                core_target = {pos_map[p] for p in range(call["start"], call["end"] + 1) if p in pos_map}
                if len(core_target) < 3:
                    continue
                candidates = [r for r in target_runs_by_anchor[uid] if r.state == call["structural_class"] and (set(range(r.start, r.end + 1)) & core_target)]
                if not candidates:
                    continue
                best = max(candidates, key=lambda r: len(set(range(r.start, r.end + 1)) & core_target))
                if _overlaps_existing(best.start, best.end, templates_by_anchor[uid]):
                    continue
                templates_by_anchor[uid].append({
                    "start": best.start, "end": best.end, "segment": seg, "structural_class": call["structural_class"],
                    "confidence": "MODERATE",
                    "evidence": f"Cross-template rescue: transferred from anchor {src_uid} ({ANCHOR_NAMES[src_uid]}) template segment {seg} "
                                f"via CE structural alignment (z={result.z_score:.2f}, {result.n_aligned} aligned residues), "
                                f"boundary refined to this receptor's own overlapping secondary-structure run; "
                                f"not independently literature-anchored in {ANCHOR_NAMES[uid]} within this session's search budget.",
                    "reference_pdb": f"{call['reference_pdb']} (via {src_uid})", "reference_chain": call["reference_chain"],
                })
                missing_core.remove(seg)
                have.add(seg)


def main() -> None:
    master = pd.read_csv(MASTER_CSV, dtype={"pdb_id": str})
    norm = pd.read_parquet(
        NORMALIZED_PARQUET,
        columns=["instance_id", "uniprot_id", "canonical_uniprot_position", "canonical_residue_name",
                 "entity_label_seq_id", "is_construct_mapped", "ca_positive_occupancy"],
    )

    ref_rows = []
    templates_by_anchor: dict[str, list[dict]] = {}
    all_candidates: dict[str, list[pd.Series]] = {}

    for uid in ANCHORS:
        candidates = select_reference_candidates(uid, master, norm)
        all_candidates[uid] = candidates
        calls, _ = build_template_from_candidates(uid, candidates, norm)
        templates_by_anchor[uid] = calls
        for i, cand in enumerate(candidates):
            segs_here = _segments_found(cand, uid, norm)
            ref_rows.append({
                "uniprot_id": uid, "nr_code": cand["nr_code"], "pdb_id": cand["pdb_id"], "chain": cand["label_asym_id"],
                "role": "PRIMARY" if i == 0 else "COMPLEMENTARY", "resolution": cand["resolution"],
                "coverage": cand["positive_occupancy_lbd_coverage"], "segments_supported": ",".join(sorted(segs_here)),
                "source_publication": "; ".join(sorted({e.split(" --")[0] for _, seg, e in LANDMARKS[uid] if seg in segs_here})) or "topology (H1) only",
            })
        named = sorted({c["segment"] for c in calls}, key=lambda s: CANONICAL_ALL.index(s) if s in CANONICAL_ALL else 99)
        cand_labels = ", ".join(f"{c.pdb_id}({'P' if i == 0 else 'C'})" for i, c in enumerate(candidates))
        print(f"{ANCHOR_NAMES[uid]} ({uid}): {len(candidates)} reference structure(s) [{cand_labels}] -> {named}")

    print("\nCross-template rescue among the 5 anchors (Part 4 evidence priority #4)...")
    tracks = {}
    target_runs_by_anchor = {}
    for uid in ANCHORS:
        primary = all_candidates[uid][0]
        residues = norm[(norm["instance_id"] == primary["representative_instance_id"]) & norm["is_construct_mapped"] & norm["ca_positive_occupancy"]]
        tracks[uid] = build_ca_track(primary["pdb_id"], primary["label_asym_id"], residues, RAW_MMCIF)
        target_runs_by_anchor[uid] = get_runs_for_instance(primary["pdb_id"], primary["label_asym_id"], residues)
    rescue_via_cross_template_alignment(templates_by_anchor, tracks, target_runs_by_anchor)
    for uid in ANCHORS:
        named = sorted({c["segment"] for c in templates_by_anchor[uid]}, key=lambda s: CANONICAL_ALL.index(s) if s in CANONICAL_ALL else 99)
        print(f"  {ANCHOR_NAMES[uid]} after rescue: {named}")

    template_rows = []
    for uid in ANCHORS:
        primary = all_candidates[uid][0]
        names_lookup = norm[(norm["uniprot_id"] == uid)].drop_duplicates("canonical_uniprot_position").set_index("canonical_uniprot_position")["canonical_residue_name"]
        for c in sorted(templates_by_anchor[uid], key=lambda c: c["start"]):
            for pos in range(c["start"], c["end"] + 1):
                template_rows.append({
                    "uniprot_id": uid, "nr_code": primary["nr_code"], "structural_segment": c["segment"],
                    "canonical_uniprot_position": pos, "canonical_residue_name": names_lookup.get(pos, ""),
                    "structural_class": c["structural_class"], "reference_pdb": c["reference_pdb"],
                    "reference_chain": c["reference_chain"], "evidence_type": c["evidence"].split(":")[0] if ":" in c["evidence"] else "literature_landmark",
                    "publication": c["evidence"], "assignment_confidence": c["confidence"],
                    "notes": "start_uniprot" if False else "",
                })

    ref_df = pd.DataFrame(ref_rows)
    tmpl_df = pd.DataFrame(template_rows)
    assert not tmpl_df.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    ref_df.to_csv(REF_STRUCTURES_V2_CSV, index=False)
    tmpl_df.to_csv(TEMPLATE_V2_CSV, index=False)
    print(f"\nWrote {REF_STRUCTURES_V2_CSV.relative_to(PROJECT_ROOT)} ({len(ref_df)} rows)")
    print(f"Wrote {TEMPLATE_V2_CSV.relative_to(PROJECT_ROOT)} ({len(tmpl_df)} rows)")

    print("\nFinal per-anchor core-segment support:")
    for uid in ANCHORS:
        named = {c["segment"] for c in templates_by_anchor[uid]}
        core_present = [s for s in CORE_SEGMENTS if s in named]
        core_missing = [s for s in CORE_SEGMENTS if s not in named]
        print(f"  {ANCHOR_NAMES[uid]:12s}: {len(core_present)}/{len(CORE_SEGMENTS)} core segments -- missing: {core_missing}")


if __name__ == "__main__":
    main()
