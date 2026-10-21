"""
Professor canonical helix map -- FINAL policy (professor decision 2026-10-20).

The professor's interval table is ground truth. This script transfers it
exactly: every interval keeps its professor start/end and label (including
`review_required`, whose status is kept as metadata only); a position inside
two professor intervals carries both labels; positions inside the analysis
span but outside every professor interval are INTER_HELIX_UNASSIGNED_DERIVED
(never a v3 loop name); receptors without professor mapping get
CANONICAL_MAPPING_UNAVAILABLE. TRalpha is numbered on the professor accession
P10827-2. No v3 label is used anywhere in the map.

The v3 comparison tables are AUDIT HISTORY only (they describe how the old map
differed); they do not influence any label.

Reads:
    data/sources/professor_canonical_helix_table/NR_canonical_helix_boundaries.csv   (ground truth, verbatim)
    data/sources/professor_canonical_helix_table/uniprot_isoform_sequences/*.fasta
    data/manifests/nr_metadata_validated.csv, data/manifests/nr_lbd_reference.csv (v1 PROSITE window, metadata)
    data/rebuild_v2/inputs/analysis_span_v2.csv                    (scripts/102)
    data/processed/structures.csv                                  (reference-structure method only)
    data/analysis/nr_lbd_structural_segment_map.csv                (historical v3; audit-history comparison only)

Writes:
    data/canonical_remap/nr_lbd_canonical_helix_intervals_professor_2026.csv
    data/canonical_remap/nr_lbd_canonical_helix_map_professor_2026.csv
    data/canonical_remap/canonical_helix_map_summary.json
    reports/canonical_remap/canonical_helix_identity_isoform_audit.csv
    reports/canonical_remap/canonical_helix_receptor_coverage_audit.csv
    reports/canonical_remap/canonical_helix_conflicts.csv
    reports/canonical_remap/canonical_helix_v3_comparison_elements.csv     (audit history)
    reports/canonical_remap/canonical_helix_v3_comparison_by_receptor.csv  (audit history)

Usage:
    python scripts/97_build_canonical_helix_map.py [--output-root DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils import canonical_helix as ch  # noqa: E402

SOURCE_DIR = PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table"
SOURCE_CSV = SOURCE_DIR / "NR_canonical_helix_boundaries.csv"
SOURCE_SHA256 = "ae1fb1005e8e4a10d7cf568ac3657e5bf692c62f5eac7572117b558b92e4268c"
ISOFORM_DIR = SOURCE_DIR / "uniprot_isoform_sequences"

METADATA_CSV = PROJECT_ROOT / "data" / "manifests" / "nr_metadata_validated.csv"
LBD_REFERENCE_CSV = PROJECT_ROOT / "data" / "manifests" / "nr_lbd_reference.csv"
SPAN_CSV = PROJECT_ROOT / "data" / "rebuild_v2" / "inputs" / "analysis_span_v2.csv"
MASTER_CSV = PROJECT_ROOT / "data" / "processed" / "structures.csv"
V3_MAP_CSV = PROJECT_ROOT / "data" / "analysis" / "nr_lbd_structural_segment_map.csv"

OUT_ANALYSIS = Path("data") / "canonical_remap"
OUT_TABLES = Path("reports") / "canonical_remap"
INTERVALS_OUT = OUT_ANALYSIS / "nr_lbd_canonical_helix_intervals_professor_2026.csv"
MAP_OUT = OUT_ANALYSIS / "nr_lbd_canonical_helix_map_professor_2026.csv"
SUMMARY_OUT = OUT_ANALYSIS / "canonical_helix_map_summary.json"
IDENTITY_OUT = OUT_TABLES / "canonical_helix_identity_isoform_audit.csv"
COVERAGE_OUT = OUT_TABLES / "canonical_helix_receptor_coverage_audit.csv"
CONFLICTS_OUT = OUT_TABLES / "canonical_helix_conflicts.csv"
V3_ELEMENTS_OUT = OUT_TABLES / "canonical_helix_v3_comparison_elements.csv"
V3_RECEPTOR_OUT = OUT_TABLES / "canonical_helix_v3_comparison_by_receptor.csv"

EXPECTED_SOURCE_ROWS = 434
CORE_ELEMENTS = ["H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10/H11", "H12"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(df: pd.DataFrame, root: Path, rel: Path) -> None:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, lineterminator="\n")


def join_unique(values) -> str:
    seen: list[str] = []
    for v in values:
        if v not in seen:
            seen.append(v)
    return "|".join(seen)


def load_source() -> pd.DataFrame:
    if sha256(SOURCE_CSV) != SOURCE_SHA256:
        raise RuntimeError(f"{SOURCE_CSV} does not match the preserved source hash; refusing to build")
    src = pd.read_csv(SOURCE_CSV, encoding=ch.SOURCE_ENCODING, dtype=str, keep_default_na=False)
    if list(src.columns) != ch.SOURCE_COLUMNS or len(src) != EXPECTED_SOURCE_ROWS:
        raise RuntimeError("Unexpected source shape")
    src.insert(0, "source_row_number", range(1, len(src) + 1))
    return src


def load_project_receptors() -> pd.DataFrame:
    meta = pd.read_csv(METADATA_CSV, dtype=str, keep_default_na=False)[["uniprot_id", "nr_code", "common_name", "group"]]
    lbd = pd.read_csv(LBD_REFERENCE_CSV)[["uniprot_id", "sequence_length", "lbd_start", "lbd_end"]]
    lbd = lbd.rename(columns={"lbd_start": "prosite_lbd_start_v1", "lbd_end": "prosite_lbd_end_v1"})
    span = pd.read_csv(SPAN_CSV)[["uniprot_id", "reference_accession", "analysis_span_start", "analysis_span_end"]]
    out = meta.merge(lbd, on="uniprot_id", validate="one_to_one").merge(span, on="uniprot_id", validate="one_to_one")
    assert len(out) == 48
    return out


def sequences() -> dict[str, str]:
    return {p.stem: ch.read_fasta_sequence(p) for p in sorted(ISOFORM_DIR.glob("*.fasta"))}


# ---------------------------------------------------------------------------


def build_identity_audit(src: pd.DataFrame, rec: pd.DataFrame, seqs: dict[str, str]) -> pd.DataFrame:
    rows = []
    for _, r in rec.sort_values("nr_code").iterrows():
        s = src[src["nr_code"] == r["nr_code"]]
        row = {"nr_code": r["nr_code"], "project_uniprot_id": r["uniprot_id"], "project_common_name": r["common_name"],
               "final_reference_accession": r["reference_accession"], "in_source": len(s) > 0}
        if not len(s):
            row.update({"identity_status": "NOT_IN_SOURCE", "source_receptor_names": "", "source_sequence_accessions": "",
                        "name_match": "", "accession_matches_final_reference": "", "reference_sequence_length": "",
                        "n_intervals": 0, "n_intervals_outside_reference_sequence": "", "isoform_sequence_check": ""})
            rows.append(row)
            continue
        names, accs = sorted(s["receptor"].unique()), sorted(s["sequence_accession"].unique())
        if len(names) != 1 or len(accs) != 1:
            raise RuntimeError(f"{r['nr_code']}: inconsistent source identity {names} {accs}")
        acc = accs[0]
        ch.final_crosswalk(acc, r["uniprot_id"])  # raises unless acc is the final numbering reference
        ref_len = len(seqs[acc]) if acc in seqs and acc != r["uniprot_id"] else int(r["sequence_length"])
        n_out = int(((s["start_uniprot"].astype(int) < 1) | (s["end_uniprot"].astype(int) > ref_len)).sum())
        iso = ""
        if acc != r["uniprot_id"]:
            fd = ch.first_sequence_difference(seqs[acc], seqs[r["uniprot_id"]])
            iso = (f"{acc} ({len(seqs[acc])} aa) vs canonical {r['uniprot_id']} ({len(seqs[r['uniprot_id']])} aa): "
                   f"identical 1-{fd - 1}; final numbering uses {acc} throughout")
        status = "MATCH" if acc == r["uniprot_id"] else "MATCH_PROFESSOR_ISOFORM_IS_FINAL_REFERENCE"
        if names[0] != r["common_name"] or n_out:
            status = "CONFLICT"
        row.update({"identity_status": status, "source_receptor_names": names[0], "source_sequence_accessions": acc,
                    "name_match": names[0] == r["common_name"], "accession_matches_final_reference": True,
                    "reference_sequence_length": ref_len, "n_intervals": len(s),
                    "n_intervals_outside_reference_sequence": n_out, "isoform_sequence_check": iso})
        rows.append(row)
    return pd.DataFrame(rows)


def build_intervals(src: pd.DataFrame, rec: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    recs = rec.set_index("nr_code")
    rows, notes = [], []

    def note(kind, code, detail, positions="", rows_="", handling=""):
        notes.append({"issue_type": kind, "nr_code": code, "project_uniprot_id": recs.loc[code, "uniprot_id"],
                      "common_name": recs.loc[code, "common_name"], "positions": positions,
                      "source_row_numbers": rows_, "detail": detail, "final_policy_handling": handling})

    for _, s in src.iterrows():
        r = recs.loc[s["nr_code"]]
        start, end = int(s["start_uniprot"]), int(s["end_uniprot"])
        pol = ch.status_policy(s["assignment_status"])
        ch.element_rank(s["canonical_helix"])
        ch.final_crosswalk(s["sequence_accession"], r["uniprot_id"])
        lo, hi = int(r["analysis_span_start"]), int(r["analysis_span_end"])
        rows.append({
            **{c: s[c] for c in ["source_row_number"] + ch.SOURCE_COLUMNS},
            "project_uniprot_id": r["uniprot_id"], "project_common_name": r["common_name"],
            "final_reference_accession": r["reference_accession"],
            "start": start, "end": end, "interval_length": end - start + 1,
            "element_rank": ch.ELEMENT_RANK[s["canonical_helix"]],
            "is_receptor_specific_element": s["canonical_helix"] in ch.RECEPTOR_SPECIFIC_ELEMENTS,
            "evidence_category": pol.evidence_category, "label_policy": ch.FINAL_LABEL_POLICY,
            "analysis_label": s["canonical_helix"], "review_required_flag": s["assignment_status"] == "review_required",
            "analysis_span_start": lo, "analysis_span_end": hi,
            "fully_inside_analysis_span": bool(lo <= start and end <= hi),
            "prosite_lbd_start_v1": int(r["prosite_lbd_start_v1"]), "prosite_lbd_end_v1": int(r["prosite_lbd_end_v1"]),
            "n_positions_outside_prosite_lbd_v1": sum(1 for p in range(start, end + 1)
                                                      if not r["prosite_lbd_start_v1"] <= p <= r["prosite_lbd_end_v1"]),
        })
        if start > end:
            note("START_GT_END", s["nr_code"], f"{s['canonical_helix']} {start}-{end}", rows_=str(s["source_row_number"]))
    iv = pd.DataFrame(rows)
    iv["overlaps_with"] = ""
    for code, g in iv.groupby("nr_code", sort=True):
        g = g.sort_values(["start", "end", "source_row_number"])
        rr = list(g.itertuples())
        for i, a in enumerate(rr):
            for b in rr[i + 1:]:
                if b.start > a.end:
                    break
                note("PROFESSOR_INTERVALS_OVERLAP", code,
                     f"{a.canonical_helix} {a.start}-{a.end} overlaps {b.canonical_helix} {b.start}-{b.end}",
                     positions=f"{b.start}-{min(a.end, b.end)}", rows_=f"{a.source_row_number}|{b.source_row_number}",
                     handling="BOTH_PROFESSOR_LABELS_APPLIED_AT_SHARED_POSITIONS")
                for x, y in ((a, b), (b, a)):
                    iv.loc[x.Index, "overlaps_with"] = join_unique(
                        [v for v in iv.loc[x.Index, "overlaps_with"].split("|") if v] + [f"{y.canonical_helix}:row{y.source_row_number}"])
        for name, gg in g.groupby("canonical_helix"):
            if len(gg) > 1:
                note("ELEMENT_GIVEN_AS_SEPARATE_INTERVALS", code,
                     f"{name}: " + ", ".join(f"{a}-{b}" for a, b in zip(gg["start"], gg["end"])),
                     rows_=join_unique(gg["source_row_number"].astype(str)),
                     handling="INTERVALS_USED_EXACTLY; INTERVENING_POSITIONS_INTER_HELIX_UNASSIGNED_DERIVED")
    iv = iv.sort_values("source_row_number").reset_index(drop=True)
    nt = pd.DataFrame(notes, columns=["issue_type", "nr_code", "project_uniprot_id", "common_name", "positions",
                                      "source_row_numbers", "detail", "final_policy_handling"])
    return iv, nt.sort_values(["nr_code", "issue_type", "source_row_numbers"]).reset_index(drop=True)


def build_residue_map(iv: pd.DataFrame, rec: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in rec.sort_values("uniprot_id").iterrows():
        uid = r["uniprot_id"]
        g = iv[iv["project_uniprot_id"] == uid]
        mapped = len(g) > 0
        cov: dict[int, list[dict]] = {}
        for t in g.to_dict("records"):
            for p in ch.expand_interval(t["start"], t["end"]):
                cov.setdefault(p, []).append(t)
        lo, hi = int(r["analysis_span_start"]), int(r["analysis_span_end"])
        frame = sorted(set(range(lo, hi + 1)) | set(cov))
        extents = sorted((t["start"], t["end"], t["canonical_helix"]) for t in g.to_dict("records"))
        for p in frame:
            covering = sorted(cov.get(p, []), key=lambda t: t["source_row_number"])
            status = ch.classify_position_final(covering, mapped)
            names = join_unique(t["canonical_helix"] for t in covering)
            before = after = ""
            if status == ch.FPOS_INTER:
                b = [e for e in extents if e[1] < p]
                a = [e for e in extents if e[0] > p]
                before = max(b, key=lambda e: e[1])[2] if b else ""
                after = min(a, key=lambda e: e[0])[2] if a else ""
            rows.append({
                "uniprot_id": uid, "nr_code": r["nr_code"], "common_name": r["common_name"],
                "reference_accession": r["reference_accession"], "canonical_uniprot_position": p,
                "within_analysis_span": bool(lo <= p <= hi),
                "within_prosite_lbd_v1": bool(r["prosite_lbd_start_v1"] <= p <= r["prosite_lbd_end_v1"])
                if r["reference_accession"] == uid else bool(r["prosite_lbd_start_v1"] <= p <= min(r["prosite_lbd_end_v1"], ch.TRA_IDENTICAL_END)),
                "receptor_has_professor_mapping": mapped, "position_status": status,
                "canonical_elements": names, "n_professor_elements": len({t["canonical_helix"] for t in covering}),
                "element_rank": ch.ELEMENT_RANK[names] if status == ch.FPOS_ELEMENT else pd.NA,
                "assignment_status": join_unique(t["assignment_status"] for t in covering),
                "review_required_flag": any(t["assignment_status"] == "review_required" for t in covering),
                "evidence_category": join_unique(t["evidence_category"] for t in covering),
                "reference_pdb": join_unique(t["reference_pdb"] for t in covering),
                "reference_chain": join_unique(t["reference_chain"] for t in covering),
                "boundary_basis": join_unique(t["boundary_basis"] for t in covering),
                "source_url": join_unique(t["source_url"] for t in covering),
                "note": join_unique(t["note"] for t in covering if t["note"]),
                "source_row_numbers": join_unique(str(t["source_row_number"]) for t in covering),
                "inter_helix_flank_before": before, "inter_helix_flank_after": after,
            })
    m = pd.DataFrame(rows).sort_values(["uniprot_id", "canonical_uniprot_position"]).reset_index(drop=True)
    assert not m.duplicated(["uniprot_id", "canonical_uniprot_position"]).any()
    return m


def build_coverage(iv: pd.DataFrame, m: pd.DataFrame, rec: pd.DataFrame, master: pd.DataFrame) -> pd.DataFrame:
    coord = master[master["has_stage3b_coordinate_data"] == "True"].groupby("uniprot_id")["pdb_id"].nunique()
    rows = []
    for _, r in rec.sort_values("nr_code").iterrows():
        uid = r["uniprot_id"]
        g = iv[iv["project_uniprot_id"] == uid].sort_values(["start", "end"])
        mm = m[m["uniprot_id"] == uid]
        st = mm["position_status"].value_counts()
        ref_pdb = join_unique(g["reference_pdb"]) if len(g) else ""
        method = join_unique(master[(master["uniprot_id"] == uid) & (master["pdb_id"] == ref_pdb)]["experimental_method"]) if ref_pdb else ""
        present = set(g["canonical_helix"]) | ({"H10/H11"} if {"H10", "H11"} <= set(g["canonical_helix"]) else set())
        rows.append({
            "nr_code": r["nr_code"], "uniprot_id": uid, "common_name": r["common_name"], "group": r["group"],
            "final_reference_accession": r["reference_accession"],
            "professor_mapping_exists": len(g) > 0,
            "mapping_status": "MAPPED" if len(g) else "CANONICAL_MAPPING_UNAVAILABLE",
            "reference_pdb": ref_pdb, "reference_chain": join_unique(g["reference_chain"]) if len(g) else "",
            "reference_experimental_method": method,
            "analysis_span": f"{r['analysis_span_start']}-{r['analysis_span_end']}",
            "n_interval_records": len(g), "n_distinct_elements": g["canonical_helix"].nunique(),
            "elements_in_sequence_order": " ".join(f"{n}{'(review)' if s == 'review_required' else ''}"
                                                   for n, s in zip(g["canonical_helix"], g["assignment_status"])),
            "receptor_specific_elements": " ".join(join_unique(g[g["is_receptor_specific_element"]]["canonical_helix"]).split("|")) if len(g) else "",
            "core_elements_not_in_professor_table": (" ".join(e for e in CORE_ELEMENTS if e not in present) if len(g) else "ALL"),
            "n_published_boundaries": int((g["assignment_status"] == "published_boundaries").sum()),
            "n_structural_homology": int((g["assignment_status"] == "structural_homology").sum()),
            "n_review_required": int((g["assignment_status"] == "review_required").sum()),
            "n_map_positions": len(mm), "n_positions_professor_element": int(st.get(ch.FPOS_ELEMENT, 0)),
            "n_positions_multiple_professor_elements": int(st.get(ch.FPOS_MULTI, 0)),
            "n_positions_inter_helix_derived": int(st.get(ch.FPOS_INTER, 0)),
            "n_pdb_entries_with_coordinates_v1": int(coord.get(uid, 0)),
            "interpretation_note": ("Helix identities are the professor's; inter-helix positions are derived, not loop names."
                                    if len(g) else "No professor mapping; canonical labels unavailable (not biological absence)."),
        })
    return pd.DataFrame(rows)


def comparison_family(name: str) -> str:
    return "H10/H11" if name in ("H10", "H11", "H10/H11") else name


def build_v3_comparison(iv: pd.DataFrame, m: pd.DataFrame, rec: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """AUDIT HISTORY ONLY. TRalpha positions > 370 are excluded (v3 used P10827 numbering)."""
    v3 = pd.read_csv(V3_MAP_CSV, keep_default_na=False).rename(columns={"structural_segment": "v3_segment"})
    v3 = v3[~((v3["uniprot_id"] == ch.TRA_PROJECT) & (v3["canonical_uniprot_position"] > ch.TRA_IDENTICAL_END))]
    v3["v3_kind"] = "UNASSIGNED"
    v3.loc[v3["v3_segment"].str.startswith("loop_"), "v3_kind"] = "LOOP"
    v3.loc[(v3["v3_segment"] != "") & ~v3["v3_segment"].str.startswith("loop_"), "v3_kind"] = "ELEMENT"
    v3["v3_family"] = v3["v3_segment"].map(lambda s: comparison_family(s) if s and not s.startswith("loop_") else "")
    v3_pos = v3.set_index(["uniprot_id", "canonical_uniprot_position"])
    rows = []
    for _, r in rec.sort_values("nr_code").iterrows():
        uid = r["uniprot_id"]
        g = iv[iv["project_uniprot_id"] == uid]
        if uid == ch.TRA_PROJECT:
            g = g[g["start"] <= ch.TRA_IDENTICAL_END]
        v3r = v3[(v3["uniprot_id"] == uid) & (v3["v3_kind"] == "ELEMENT")]
        v3_fams = {f: set(gg["canonical_uniprot_position"]) for f, gg in v3r.groupby("v3_family")}
        v3_names = {f: join_unique(sorted(gg["v3_segment"].unique())) for f, gg in v3r.groupby("v3_family")}
        elems: dict[str, set] = {}
        for t in g.to_dict("records"):
            elems.setdefault(comparison_family(t["canonical_helix"]), set()).update(range(t["start"], t["end"] + 1))
        for fam, pos in sorted(elems.items(), key=lambda kv: min(kv[1])):
            labs = [v3_pos.loc[(uid, p)] if (uid, p) in v3_pos.index else None for p in sorted(pos)]
            overl = sorted({l["v3_segment"] for l in labs if l is not None and l["v3_kind"] == "ELEMENT"})
            same = v3_fams.get(fam, set())
            change = ("MATCHED_SAME_ELEMENT" if same & pos else "SAME_NAME_IN_V3_BUT_NO_POSITIONAL_OVERLAP" if same
                      else "RENAMED_FROM_V3" if overl else "RECOVERED_NOT_NAMED_IN_V3")
            rows.append({"nr_code": r["nr_code"], "uniprot_id": uid, "common_name": r["common_name"],
                         "direction": "PROFESSOR_ELEMENT_VS_V3", "element": fam, "change_class": change,
                         "professor_start": min(pos), "professor_end": max(pos),
                         "v3_same_name_start": min(same) if same else pd.NA, "v3_same_name_end": max(same) if same else pd.NA,
                         "v3_elements_overlapping": "|".join(overl),
                         "n_positions_v3_loop": sum(1 for l in labs if l is not None and l["v3_kind"] == "LOOP")})
        for fam, pos in sorted(v3_fams.items(), key=lambda kv: min(kv[1])):
            if fam in elems and elems[fam] & pos:
                continue
            now = sorted({n for n in m[(m["uniprot_id"] == uid) & m["canonical_uniprot_position"].isin(pos)]["canonical_elements"] if n})
            change = ("RECEPTOR_NOT_IN_SOURCE" if not len(g) else "NOT_REPRESENTED_IN_HELIX_SOURCE" if fam == "beta_sheet"
                      else ("RENAMED_TO_PROFESSOR_ELEMENT" if now else "NOT_A_PROFESSOR_ELEMENT"))
            rows.append({"nr_code": r["nr_code"], "uniprot_id": uid, "common_name": r["common_name"],
                         "direction": "V3_ELEMENT_VS_PROFESSOR", "element": v3_names[fam], "change_class": change,
                         "professor_start": pd.NA, "professor_end": pd.NA,
                         "v3_same_name_start": min(pos), "v3_same_name_end": max(pos),
                         "v3_elements_overlapping": "|".join(now), "n_positions_v3_loop": 0})
    el = pd.DataFrame(rows)
    rec_rows = []
    for (code, uid, name), e in el.groupby(["nr_code", "uniprot_id", "common_name"], sort=True):
        p = e[e["direction"] == "PROFESSOR_ELEMENT_VS_V3"]
        rec_rows.append({"nr_code": code, "uniprot_id": uid, "common_name": name,
                         "professor_elements_not_named_in_v3": " ".join(p[p["change_class"] == "RECOVERED_NOT_NAMED_IN_V3"]["element"]),
                         "professor_elements_renamed_from_v3": " ".join(f"{a}<-v3:{b}" for a, b in zip(
                             p[p["change_class"] == "RENAMED_FROM_V3"]["element"], p[p["change_class"] == "RENAMED_FROM_V3"]["v3_elements_overlapping"])),
                         "professor_elements_same_name_elsewhere_in_v3": " ".join(
                             p[p["change_class"] == "SAME_NAME_IN_V3_BUT_NO_POSITIONAL_OVERLAP"]["element"]),
                         "n_v3_loop_positions_inside_professor_elements": int(p["n_positions_v3_loop"].sum()),
                         "note": "AUDIT HISTORY ONLY - does not influence any final label"})
    return el, pd.DataFrame(rec_rows)


def build(output_root: Path) -> dict:
    src = load_source()
    rec = load_project_receptors()
    seqs = sequences()
    master = pd.read_csv(MASTER_CSV, dtype=str, keep_default_na=False,
                         usecols=["uniprot_id", "pdb_id", "experimental_method", "has_stage3b_coordinate_data"])
    identity = build_identity_audit(src, rec, seqs)
    iv, notes = build_intervals(src, rec)
    m = build_residue_map(iv, rec)
    cov = build_coverage(iv, m, rec, master)
    v3_el, v3_rec = build_v3_comparison(iv, m, rec)
    for df, rel in [(iv, INTERVALS_OUT), (m, MAP_OUT), (identity, IDENTITY_OUT), (cov, COVERAGE_OUT),
                    (notes, CONFLICTS_OUT), (v3_el, V3_ELEMENTS_OUT), (v3_rec, V3_RECEPTOR_OUT)]:
        write_csv(df, output_root, rel)
    summary = {
        "policy": ch.FINAL_POLICY_VERSION, "source_file": str(SOURCE_CSV.relative_to(PROJECT_ROOT)),
        "source_sha256": SOURCE_SHA256, "n_source_rows": len(src), "n_source_receptors": int(src["nr_code"].nunique()),
        "assignment_status_counts": {k: int(v) for k, v in src["assignment_status"].value_counts().sort_index().items()},
        "receptors_without_professor_mapping": cov[~cov["professor_mapping_exists"]]["common_name"].tolist(),
        "identity_status_counts": {k: int(v) for k, v in identity["identity_status"].value_counts().sort_index().items()},
        "position_status_counts": {k: int(v) for k, v in m["position_status"].value_counts().sort_index().items()},
        "n_map_rows": len(m), "n_intervals_fully_inside_analysis_span": int(iv["fully_inside_analysis_span"].sum()),
        "inputs_sha256": {str(p.relative_to(PROJECT_ROOT)): sha256(p) for p in
                          [METADATA_CSV, LBD_REFERENCE_CSV, SPAN_CSV, MASTER_CSV, V3_MAP_CSV]},
    }
    (output_root / SUMMARY_OUT).parent.mkdir(parents=True, exist_ok=True)
    (output_root / SUMMARY_OUT).write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    print(json.dumps(build(ap.parse_args().output_root), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
