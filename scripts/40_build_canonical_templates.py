"""
Structural-segment alignment v2, Part 3: canonical H1-H12/beta-sheet
templates for the 5 anchor receptors (VDR, PPARgamma, ESR1, GR, NR4A1).

Canonical labels come from structural/literature evidence, NEVER from
numbering DSSP helices sequentially. Evidence used, all with an explicit
citation per assignment:

1. This repository's own literature-sourced functional-site data
   (reports/tables/functional_site_definitions.csv,
   data/manifests/nr_functional_features.csv) -- charge-clamp residues
   (Lys on H3, Glu on H12) and, for NR4A1, an explicit UniProt-annotated
   H12 segment.
2. External literature landmarks gathered for this stage (cited inline
   below) confirming the charge-clamp assignments and general NR-LBD
   topology (H5-H6 flanked by a short beta-hairpin; H7/H9/H10/H11 at the
   dimerization interface) -- Nolte et al. 1998 (PPARgamma, PMID 9813012);
   general NR-LBD architecture review (PMC3637177).
3. This receptor's own reference-structure secondary-structure evidence
   (_struct_conf/_struct_sheet_range, already extracted in v1) -- used
   ONLY to locate which helix RUN contains a literature-anchored residue,
   and, for interior helices strictly bracketed by two independent
   anchors, to sequentially label them ONLY when the observed interior
   count exactly matches the canonical expectation for that specific
   bracketed interval (never across the whole H1-H12 span at once, and
   always with a MODERATE, not HIGH, confidence flag and an explicit note).

Where none of the above applies, the position is marked
PENDING_EXTERNAL_VALIDATION rather than guessed.

Reads:
    data/analysis/lbd_structural_reference_structures.csv
    data/raw/mmcif/{PDB_ID}.cif.gz (already-cached; no download)
    data/analysis/normalized_lbd_bfactors_all.parquet (identity/mapping columns only)
    reports/tables/functional_site_definitions.csv
    data/manifests/nr_functional_features.csv

Writes:
    data/analysis/canonical_lbd_template_segments.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.mmcif_parse import get_secondary_structure_ranges, parse_mmcif_gz  # noqa: E402
from scripts.utils.segment_topology import Run, build_runs, merge_sheet_runs  # noqa: E402

ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"
TABLES = PROJECT_ROOT / "reports" / "tables"
RAW_MMCIF = PROJECT_ROOT / "data" / "raw" / "mmcif"

REF_CSV = ANALYSIS / "lbd_structural_reference_structures.csv"
NORMALIZED_PARQUET = ANALYSIS / "normalized_lbd_bfactors_all.parquet"
FUNCTIONAL_SITE_CSV = TABLES / "functional_site_definitions.csv"
FUNCTIONAL_FEATURES_CSV = MANIFESTS / "nr_functional_features.csv"

OUT_CSV = ANALYSIS / "canonical_lbd_template_segments.csv"

ANCHORS = {"P11473": "VDR", "P37231": "PPARgamma", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}

# Literature/repo-sourced landmark residues (canonical UniProt position -> segment),
# each with an explicit citation. H3/H12 charge-clamp pairs are the primary anchors.
LANDMARKS: dict[str, list[tuple[int, str, str]]] = {
    "P11473": [  # VDR
        (235, "H3", "Rochel et al. 2000 (Mol Cell 5:173) -- Ser235(H3) H-bonds Thr415(H12)."),
        (264, "H4", "Rochel et al. 2000 -- Lys264(H4)-Glu420(H12) salt bridge stabilizing H12."),
        (420, "H12", "Rochel et al. 2000 -- Glu420 of H12, salt-bridge partner of Lys264(H4)."),
        (415, "H12", "Rochel et al. 2000 -- Thr415, H12 N-terminal boundary residue."),
    ],
    "P37231": [  # PPARgamma (repo functional_site_definitions.csv canonical numbering)
        (329, "H3", "functional_site_definitions.csv af2_coactivator_groove (Nolte et al. 1998, PMID 9813012) -- Lys(H3) charge-clamp residue."),
        (499, "H12", "functional_site_definitions.csv af2_coactivator_groove (Nolte et al. 1998, PMID 9813012) -- Glu(H12/AF-2) charge-clamp residue."),
    ],
    "P03372": [  # ESR1 / ERalpha
        (362, "H3", "functional_site_definitions.csv af2_coactivator_groove (Shiau et al. 1998, DOI 10.1016/S0092-8674(00)81717-1) -- Lys362(H3) charge-clamp residue."),
        (542, "H12", "functional_site_definitions.csv af2_coactivator_groove (Shiau et al. 1998) -- Glu542(H12/AF-2) charge-clamp residue."),
        (504, "H11", "External literature (dimer-interface studies of hERalpha LBD, e.g. PMC10760953 and related dimerization studies) -- Leu504 of the H10/H11 leucine-zipper-like dimer-interface hydrophobic patch."),
        (508, "H11", "External literature -- Leu508, H10/H11 dimer-interface hydrophobic patch."),
        (511, "H11", "External literature -- Leu511, H10/H11 dimer-interface hydrophobic patch."),
    ],
    "P04150": [  # GR / NR3C1
        (579, "H3", "functional_site_definitions.csv af2_coactivator_groove (Bledsoe et al. 2002, DOI 10.1016/S0092-8674(02)00817-6) -- Lys579(H3) charge-clamp residue."),
        (755, "H12", "functional_site_definitions.csv af2_coactivator_groove (Bledsoe et al. 2002) -- Glu755(H12/AF-2) charge-clamp residue."),
        # NOTE: functional_site_definitions.csv's GR "lbd_dimerization" residues 625/628 are
        # explicitly the D-LOOP (a loop, by its own name) -- tested against this receptor's own
        # helix runs and confirmed to fall in a coil/loop region, not inside a helix, so they
        # cannot anchor a helix segment under this algorithm (a loop landmark, not a helix one).
        # Left out rather than mis-assigned to the nearest helix.
    ],
    "P22736": [  # NR4A1 / NUR77
        (584, "H12", "nr_functional_features.csv h12_segment (UniProt:P22736) -- explicit UniProt-annotated H12 segment start."),
        (595, "H12", "nr_functional_features.csv h12_segment (UniProt:P22736) -- explicit UniProt-annotated H12 segment end."),
    ],
}

CANONICAL_BEFORE_SHEET = ["H2", "H3", "H4", "H5"]
CANONICAL_AFTER_SHEET = ["H6", "H7", "H8", "H9", "H10", "H11"]


def get_template_runs(pdb_id: str, label_asym_id: str, residues: pd.DataFrame) -> list[Run]:
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
    runs = [r for r in runs if r.state in ("HELIX", "SHEET")]
    return merge_sheet_runs(runs)


def find_run_containing(runs: list[Run], position: int, state: str) -> Run | None:
    for r in runs:
        if r.state == state and r.start <= position <= r.end:
            return r
    return None


def build_template(uid: str, name: str, ref_row: pd.Series, residues: pd.DataFrame) -> list[dict]:
    runs = get_template_runs(ref_row["pdb_id"], ref_row["label_asym_id"], residues)
    helix_runs = sorted([r for r in runs if r.state == "HELIX"], key=lambda r: r.start)
    sheet_runs = [r for r in runs if r.state == "SHEET"]

    anchor_runs: dict[str, Run] = {}
    anchor_evidence: dict[str, list[str]] = {}
    for pos, seg, evidence in LANDMARKS[uid]:
        run = find_run_containing(helix_runs, pos, "HELIX")
        if run is None:
            continue
        if seg in anchor_runs and anchor_runs[seg] is not run:
            # two landmarks for the same segment landed in different runs -- keep the
            # union span (e.g. VDR H12: Thr415 and Glu420 both cited, should be one run).
            lo = Run("HELIX", min(anchor_runs[seg].start, run.start), max(anchor_runs[seg].end, run.end))
            anchor_runs[seg] = lo
        else:
            anchor_runs[seg] = run
        anchor_evidence.setdefault(seg, []).append(evidence)

    # H1: first helix run, always anchored by topology (domain's own N-terminus) --
    # this is the one ordinal-independent, universally valid anchor (see design note).
    if "H1" not in anchor_runs and helix_runs:
        anchor_runs["H1"] = helix_runs[0]
        anchor_evidence["H1"] = ["First helix run at the standardized LBD's own N-terminal boundary (domain-topology anchor, not sequential DSSP counting)."]

    calls: list[dict] = []
    for seg, run in anchor_runs.items():
        calls.append({"start": run.start, "end": run.end, "segment": seg, "structural_class": "HELIX",
                       "confidence": "HIGH", "evidence": " | ".join(anchor_evidence[seg])})

    # beta_sheet: single unique sheet run flanked by helices on both sides.
    if len(sheet_runs) == 1:
        sh = sheet_runs[0]
        if any(h.end < sh.start for h in helix_runs) and any(h.start > sh.end for h in helix_runs):
            calls.append({"start": sh.start, "end": sh.end, "segment": "beta_sheet", "structural_class": "SHEET",
                          "confidence": "HIGH", "evidence": "Single consensus beta-sheet region flanked by helices on both sides (Renaud & Moras 2000 topology: beta-hairpin follows H5)."})

    h3 = anchor_runs.get("H3")
    h12 = anchor_runs.get("H12")
    sheet_call = next((c for c in calls if c["segment"] == "beta_sheet"), None)
    already_assigned_ranges = {(r.start, r.end) for r in anchor_runs.values()}

    def interior(lo_run: Run | None, hi_run: Run | None) -> list[Run]:
        lo_end = lo_run.end if lo_run else -10**9
        hi_start = hi_run.start if hi_run else 10**9
        return [h for h in helix_runs if lo_end < h.start and h.end < hi_start and (h.start, h.end) not in already_assigned_ranges]

    if h3 is not None and sheet_call is not None:
        between = interior(h3, Run("SHEET", sheet_call["start"], sheet_call["end"]))
        expected_before = [seg for seg in CANONICAL_BEFORE_SHEET if seg not in anchor_runs]
        if expected_before and len(between) == len(expected_before):
            for helix, seg in zip(between, expected_before):
                calls.append({"start": helix.start, "end": helix.end, "segment": seg, "structural_class": "HELIX",
                              "confidence": "MODERATE",
                              "evidence": f"Bracketed between literature-anchored H3 and the identified beta_sheet; exactly {len(between)} interior helices observed, matching the canonical {expected_before} count -- labelled by order within this short, doubly-anchored interval (not by counting from H1)."})
        else:
            for helix in between:
                calls.append({"start": helix.start, "end": helix.end, "segment": "", "structural_class": "HELIX",
                              "confidence": "NONE",
                              "evidence": f"{len(between)} helices observed between H3 and beta_sheet, expected {len(expected_before)} ({expected_before}) -- count mismatch, not assigned; PENDING_EXTERNAL_VALIDATION."})

    if sheet_call is not None and h12 is not None:
        between = interior(Run("SHEET", sheet_call["start"], sheet_call["end"]), h12)
        expected_after = [seg for seg in CANONICAL_AFTER_SHEET if seg not in anchor_runs]
        if expected_after and len(between) == len(expected_after):
            for helix, seg in zip(between, expected_after):
                calls.append({"start": helix.start, "end": helix.end, "segment": seg, "structural_class": "HELIX",
                              "confidence": "MODERATE",
                              "evidence": f"Bracketed between the identified beta_sheet and literature-anchored H12; exactly {len(between)} interior helices observed, matching the canonical {expected_after} count -- labelled by order within this doubly-anchored interval."})
        else:
            for helix in between:
                calls.append({"start": helix.start, "end": helix.end, "segment": "", "structural_class": "HELIX",
                              "confidence": "NONE",
                              "evidence": f"{len(between)} helices observed between beta_sheet and H12, expected {len(expected_after)} ({expected_after}) -- count mismatch, not assigned; PENDING_EXTERNAL_VALIDATION."})

    assigned_ranges = {(c["start"], c["end"]) for c in calls}
    for helix in helix_runs:
        if (helix.start, helix.end) not in assigned_ranges and not any(c["start"] <= helix.start and helix.end <= c["end"] for c in calls):
            calls.append({"start": helix.start, "end": helix.end, "segment": "", "structural_class": "HELIX",
                          "confidence": "NONE", "evidence": "No literature anchor within reach (not adjacent to a doubly-anchored bracketed interval); PENDING_EXTERNAL_VALIDATION."})
    for sh in sheet_runs:
        if sheet_call is None or (sh.start, sh.end) != (sheet_call["start"], sheet_call["end"]):
            calls.append({"start": sh.start, "end": sh.end, "segment": "", "structural_class": "SHEET",
                          "confidence": "NONE", "evidence": "Additional/disjoint sheet run, not the uniquely identified beta_sheet; PENDING_EXTERNAL_VALIDATION."})

    rows = []
    positions = residues.set_index("entity_label_seq_id")["canonical_uniprot_position"]
    auth_map = residues.set_index("canonical_uniprot_position")["auth_seq_ids_json"] if "auth_seq_ids_json" in residues.columns else None
    for c in sorted(calls, key=lambda c: c["start"]):
        for pos in range(c["start"], c["end"] + 1):
            rows.append(
                {
                    "uniprot_id": uid, "nr_code": ref_row["nr_code"], "pdb_id": ref_row["pdb_id"],
                    "pdb_residue": pos, "canonical_uniprot_position": pos,
                    "canonical_segment": c["segment"], "structural_class": c["structural_class"],
                    "assignment_confidence": c["confidence"],
                    "assignment_status": ("ASSIGNED_HIGH_CONFIDENCE" if c["confidence"] == "HIGH" else
                                           "ASSIGNED_MODERATE_CONFIDENCE" if c["confidence"] == "MODERATE" else
                                           "PENDING_EXTERNAL_VALIDATION"),
                    "evidence_source": c["evidence"],
                }
            )
    return rows


def select_template_candidate(uid: str, master: pd.DataFrame, norm: pd.DataFrame) -> pd.Series:
    """The Part-2 reference (chosen for overall coverage/completeness) is not
    guaranteed to have the SPECIFIC literature-landmark residues resolved --
    e.g. VDR's best-coverage structure still frequently lacks H12 coordinates
    (see reports/tables/segment_feature_missingness.md, VDR H12 example).
    For template construction specifically, search this receptor's full
    coordinate-candidate pool (best coverage first) for the first structure
    whose own secondary-structure evidence actually anchors every landmark
    segment this receptor needs (H3 and H12 at minimum); fall back to the
    Part-2 reference, documented, if none does."""
    landmarks = LANDMARKS[uid]
    needed_segments = {seg for _, seg, _ in landmarks}
    pool = master[(master["uniprot_id"] == uid) & master["has_stage3b_coordinate_data"]].sort_values(
        "positive_occupancy_lbd_coverage", ascending=False
    )
    for _, cand in pool.iterrows():
        residues = norm[(norm["instance_id"] == cand["representative_instance_id"]) & norm["is_construct_mapped"] & norm["ca_positive_occupancy"]]
        if len(residues) == 0:
            continue
        runs = get_template_runs(cand["pdb_id"], cand["label_asym_id"], residues)
        helix_runs = [r for r in runs if r.state == "HELIX"]
        found = {seg for pos, seg, _ in landmarks if find_run_containing(helix_runs, pos, "HELIX") is not None}
        if needed_segments <= found:
            out = cand.copy()
            out["instance_id"] = cand["representative_instance_id"]
            out["template_selection_note"] = (
                f"Selected by landmark-completeness search over {len(pool)} coordinate candidates "
                f"(coverage-ranked): all {sorted(needed_segments)} landmarks resolved here."
            )
            return out
    fallback = pool.iloc[0].copy()
    fallback["instance_id"] = fallback["representative_instance_id"]
    fallback["template_selection_note"] = (
        f"No candidate among {len(pool)} fully resolved all landmarks {sorted(needed_segments)}; "
        "fell back to the Part-2 overall reference structure (highest coverage)."
    )
    return fallback


def main() -> None:
    master = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "structures.csv", dtype={"pdb_id": str})
    norm = pd.read_parquet(
        NORMALIZED_PARQUET,
        columns=["instance_id", "uniprot_id", "canonical_uniprot_position", "entity_label_seq_id",
                 "is_construct_mapped", "ca_positive_occupancy"],
    )

    template_rows = []
    all_rows = []
    for uid, name in ANCHORS.items():
        r = select_template_candidate(uid, master, norm)
        residues = norm[(norm["instance_id"] == r["instance_id"]) & norm["is_construct_mapped"] & norm["ca_positive_occupancy"]]
        rows = build_template(uid, name, r, residues)
        all_rows.extend(rows)
        named = {row["canonical_segment"] for row in rows if row["canonical_segment"]}
        print(f"{name} ({uid}, template {r['pdb_id']}.{r['label_asym_id']}, coverage={r['positive_occupancy_lbd_coverage']:.3f}): "
              f"{sorted(named, key=lambda s: (s != 'H1', s != 'H12', s))}")
        print(f"    {r['template_selection_note']}")
        template_rows.append({"uniprot_id": uid, "nr_code": r["nr_code"], "pdb_id": r["pdb_id"],
                               "instance_id": r["instance_id"], "label_asym_id": r["label_asym_id"],
                               "positive_occupancy_lbd_coverage": r["positive_occupancy_lbd_coverage"],
                               "selection_note": r["template_selection_note"]})
    pd.DataFrame(template_rows).to_csv(ANALYSIS / "canonical_lbd_template_reference_choice.csv", index=False)

    out = pd.DataFrame(all_rows).sort_values(["uniprot_id", "canonical_uniprot_position"]).reset_index(drop=True)
    assert not out.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV.relative_to(PROJECT_ROOT)} ({len(out)} rows)")
    print(out["assignment_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
