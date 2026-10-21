"""
Pre-build audit: professor canonical helix intervals vs the frozen standardized
LBD analysis window (UniProt 'NR LBD' Domain / PROSITE PS51843, nr_lbd_reference.csv).
Read-only; changes no boundary and no frozen file.

For every professor interval: residues inside/outside the window, truncation
class, the most likely reason for the difference, whether structures in the
project actually model the outside residues (i.e. whether widening the window
would add B-factor data), analysis impact, and a per-receptor recommendation.

"Data if extended" uses the corrected residue-numbering crosswalk v1: a residue
counts when it is mapped there AND modelled in _pdbx_poly_seq_scheme (CA
presence and B-factor usability are only established by the rebuild itself).
TRalpha intervals are in P10827-2 (TRalpha1) numbering; positions 371-410 have
no P10827 counterpart, so for them the chains are aligned to P10827-2 directly.

Reads:
    data/canonical_remap/nr_lbd_canonical_helix_intervals_professor_2026.csv   (scripts/97)
    data/canonical_remap/residue_numbering_crosswalk_v1.parquet                (scripts/96)
    data/manifests/nr_lbd_reference.csv, data/raw/api/uniprot/<acc>.json, data/raw/mmcif/<PDB>.cif.gz
    data/sources/professor_canonical_helix_table/uniprot_isoform_sequences/P10827-2.fasta

Writes:
    reports/canonical_remap/lbd_window_interval_audit_v1.csv
    reports/canonical_remap/lbd_window_receptor_summary_v1.csv
    reports/canonical_remap/lbd_window_helix_summary_v1.csv

Usage:
    python scripts/100_audit_prebuild_lbd_window.py [--output-root DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils import canonical_helix as ch  # noqa: E402
from scripts.utils import sequence_numbering as sn  # noqa: E402
from scripts.utils.mmcif_parse import parse_mmcif_gz  # noqa: E402

INTERVALS = PROJECT_ROOT / "data" / "canonical_remap" / "nr_lbd_canonical_helix_intervals_professor_2026.csv"
CROSSWALK = PROJECT_ROOT / "data" / "canonical_remap" / "residue_numbering_crosswalk_v1.parquet"
LBD_REFERENCE = PROJECT_ROOT / "data" / "manifests" / "nr_lbd_reference.csv"
UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"
MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"
TRA1_FASTA = PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table" / "uniprot_isoform_sequences" / "P10827-2.fasta"
OUT = Path("reports") / "canonical_remap"

MINOR_MAX_RESIDUES = 2
MINOR_MAX_FRACTION = 0.15
RELEVANT_FRACTION = 0.25
FUNCTIONAL_ELEMENTS = {"H12"}  # AF-2 helix: a professor-requested functional site


def uniprot_lbd_feature(acc: str) -> tuple[int, int] | None:
    d = json.loads((UNIPROT_DIR / f"{acc}.json").read_text())
    f = [x for x in d.get("features", []) if x.get("type") == "Domain" and x.get("description") == "NR LBD"]
    if len(f) != 1:
        return None
    return int(f[0]["location"]["start"]["value"]), int(f[0]["location"]["end"]["value"])


def tra1_tail_coverage(cw: pd.DataFrame) -> dict[int, set[str]]:
    """P10827-2 position -> PDB entries whose TRalpha chain models that residue,
    from a direct alignment to TRalpha1 (positions > 370 only)."""
    tra1 = ch.read_fasta_sequence(TRA1_FASTA)
    out: dict[int, set[str]] = {}
    inst = cw[cw["uniprot_id"] == "P10827"][["pdb_id", "label_asym_id"]].drop_duplicates()
    for pdb, g in inst.groupby("pdb_id"):
        d = parse_mmcif_gz(MMCIF_DIR / f"{pdb}.cif.gz")
        parents = sn.modified_residue_parents(d)
        for asym in g["label_asym_id"]:
            scheme = sn.chain_scheme(d, asym)
            q = "".join(sn.one_letter(r["mon_id"], parents) for r in scheme)
            res = sn.align_unique(q, tra1)
            calls, _ = sn.call_residues(q, tra1, res)
            for r, c in zip(scheme, calls):
                if c.target_index is not None and c.target_index + 1 > ch.TRA_IDENTICAL_END and r["modelled"]:
                    out.setdefault(c.target_index + 1, set()).add(pdb)
    return out


def classify_window(start: int, end: int, lo: int, hi: int) -> str:
    if start >= lo and end <= hi:
        return "FULLY_INSIDE"
    if end < lo or start > hi:
        return "FULLY_OUTSIDE"
    if start < lo and end > hi:
        return "PARTIALLY_TRUNCATED_BOTH_ENDS"
    return "PARTIALLY_TRUNCATED_N_TERMINAL" if start < lo else "PARTIALLY_TRUNCATED_C_TERMINAL"


def build(output_root: Path) -> dict:
    iv = pd.read_csv(INTERVALS, keep_default_na=False)
    ref = pd.read_csv(LBD_REFERENCE).set_index("uniprot_id")
    cw = pd.read_parquet(CROSSWALK, columns=["uniprot_id", "pdb_id", "instance_id", "label_asym_id",
                                             "corrected_uniprot_position", "modelled_in_poly_seq_scheme"])
    modelled = cw[cw["corrected_uniprot_position"].notna() & cw["modelled_in_poly_seq_scheme"]]
    entries_at = modelled.groupby(["uniprot_id", "corrected_uniprot_position"])["pdb_id"].agg(set).to_dict()
    n_entries = cw.groupby("uniprot_id")["pdb_id"].nunique()
    tra_tail = tra1_tail_coverage(cw)

    rows = []
    for t in iv.itertuples():
        uid = t.project_uniprot_id
        lo, hi = int(ref.loc[uid, "lbd_start"]), int(ref.loc[uid, "lbd_end"])
        feat = uniprot_lbd_feature(uid)
        start, end = int(t.start), int(t.end)
        positions = list(range(start, end + 1))
        cw_ = ch.crosswalk_for(t.sequence_accession, uid)
        not_joinable = [p for p in positions if ch.to_project_position(cw_, p) is None]
        joinable = [p for p in positions if ch.to_project_position(cw_, p) is not None]
        inside = [p for p in joinable if lo <= p <= hi]
        outside = [p for p in joinable if not lo <= p <= hi]
        win_class = classify_window(start, end, lo, hi) if not not_joinable else (
            "ISOFORM_POSITIONS_NOT_REPRESENTABLE_IN_PROJECT_ACCESSION")
        # data that a wider window (or isoform-aware numbering) would bring in
        gap_positions = outside + not_joinable
        cov = [(entries_at.get((uid, p), set()) if p in outside else tra_tail.get(p, set())) for p in gap_positions]
        entries_any = set().union(*cov) if cov else set()
        n_pos_with_data = sum(1 for c in cov if c)
        n_out = len(gap_positions)
        frac_out = n_out / len(positions)
        # reason
        reasons = []
        if not_joinable:
            reasons.append("ISOFORM_NUMBERING")
        if outside:
            if feat == (lo, hi):
                reasons.append("UNIPROT_PROSITE_NR_LBD_DOMAIN_NARROWER_THAN_STRUCTURAL_HELIX")
            else:
                reasons.append("OLD_PIPELINE_WINDOW_DIFFERS_FROM_UNIPROT_FEATURE")
            if t.nr_code.startswith("NR1C") and t.canonical_helix in ("H1", "H2"):
                reasons.append("RECEPTOR_SPECIFIC_EXTENSION(PPAR N-terminal H1/H2 before the PROSITE profile start)")
            if t.canonical_helix == "H12" and t.nr_code.startswith(("NR3C", "NR3B", "NR4A")):
                reasons.append("RECEPTOR_SPECIFIC_EXTENSION(AF-2/C-terminal region beyond the PROSITE profile end)")
        if n_out and t.assignment_status == "review_required":
            reasons.append("PROFESSOR_INTERVAL_UNCERTAINTY(review_required)")
        # impact
        if not n_out:
            impact = ""
        elif not entries_any:
            impact = "NO_PRACTICAL_IMPACT"
        elif win_class in ("FULLY_OUTSIDE", "ISOFORM_POSITIONS_NOT_REPRESENTABLE_IN_PROJECT_ACCESSION") and not inside:
            impact = "ANALYSIS_RELEVANT" if t.assignment_status == "review_required" else "BLOCKING"
        elif n_out <= MINOR_MAX_RESIDUES and frac_out <= MINOR_MAX_FRACTION:
            impact = "MINOR_TRUNCATION"
        elif t.canonical_helix in FUNCTIONAL_ELEMENTS or frac_out >= RELEVANT_FRACTION:
            impact = "ANALYSIS_RELEVANT"
        else:
            impact = "MINOR_TRUNCATION"
        rows.append({
            "receptor": t.receptor, "nr_code": t.nr_code, "project_uniprot_id": uid,
            "sequence_accession": t.sequence_accession, "canonical_helix": t.canonical_helix,
            "source_row_number": t.source_row_number, "professor_start": start, "professor_end": end,
            "interval_length": len(positions), "project_lbd_start": lo, "project_lbd_end": hi,
            "uniprot_nr_lbd_feature": f"{feat[0]}-{feat[1]}" if feat else "",
            "n_residues_inside_window": len(inside), "n_residues_outside_window": len(outside),
            "n_residues_not_representable_isoform": len(not_joinable),
            "outside_positions": (f"{min(gap_positions)}-{max(gap_positions)}" if gap_positions else ""),
            "window_class": win_class, "assignment_status": t.assignment_status, "boundary_basis": t.boundary_basis,
            "source_url": t.source_url,
            "n_outside_positions_modelled_in_any_entry": n_pos_with_data,
            "n_entries_modelling_any_outside_position": len(entries_any),
            "n_entries_receptor": int(n_entries.get(uid, 0)),
            "data_if_window_extended": ("NO_PROJECT_STRUCTURE_MODELS_THESE_RESIDUES" if n_out and not entries_any else
                                        "YES" if entries_any else ""),
            "difference_reason": "; ".join(reasons), "impact": impact,
        })
    a = pd.DataFrame(rows)

    rec = []
    for (code, name, uid), g in a.groupby(["nr_code", "receptor", "project_uniprot_id"], sort=True):
        d = g[g["impact"] != ""]
        imps = set(d["impact"])
        if d.empty:
            r = "A_KEEP_CURRENT_LBD"
        elif (d["n_residues_not_representable_isoform"] > 0).any():
            r = "D_MANUAL_REVIEW_REQUIRED"
        elif imps & {"BLOCKING", "ANALYSIS_RELEVANT"}:
            r = "B_EXTEND_LBD_TO_CANONICAL_HELIX_BOUNDARY"
        else:
            r = "C_KEEP_CURRENT_LBD_BUT_FLAG_TRUNCATION"
        rec.append({
            "nr_code": code, "receptor": name, "uniprot_id": uid,
            "project_lbd": f"{g['project_lbd_start'].iloc[0]}-{g['project_lbd_end'].iloc[0]}",
            "canonical_span": f"{g['professor_start'].min()}-{g['professor_end'].max()}",
            "n_intervals": len(g), "n_fully_inside": int((g["window_class"] == "FULLY_INSIDE").sum()),
            "n_partially_truncated": int(g["window_class"].str.startswith("PARTIALLY").sum()),
            "n_fully_outside": int((g["window_class"] == "FULLY_OUTSIDE").sum()),
            "n_isoform_not_representable": int((g["window_class"].str.startswith("ISOFORM")).sum()),
            "discrepant_intervals": " ".join(f"{h}[{s}-{e}:{c.split('_')[0]}{'/' + i if i else ''}]" for h, s, e, c, i in zip(
                d["canonical_helix"], d["professor_start"], d["professor_end"], d["window_class"], d["impact"])),
            "h12_affected": bool((d["canonical_helix"] == "H12").any()),
            "worst_impact": next((x for x in ["BLOCKING", "ANALYSIS_RELEVANT", "MINOR_TRUNCATION", "NO_PRACTICAL_IMPACT"]
                                  if x in imps), ""),
            "window_extension_needed_to_cover_canonical_helices": (
                f"{min(lo_ := g['professor_start'].min(), g['project_lbd_start'].iloc[0])}-"
                f"{max(g['professor_end'].max(), g['project_lbd_end'].iloc[0])}" if not d.empty else ""),
            "recommendation": r,
        })
    r = pd.DataFrame(rec)
    hel = a.groupby("canonical_helix").agg(
        n_intervals=("receptor", "size"),
        n_fully_inside=("window_class", lambda s: int((s == "FULLY_INSIDE").sum())),
        n_partially_truncated=("window_class", lambda s: int(s.str.startswith("PARTIALLY").sum())),
        n_fully_outside=("window_class", lambda s: int((s == "FULLY_OUTSIDE").sum())),
        n_isoform_not_representable=("window_class", lambda s: int(s.str.startswith("ISOFORM").sum())),
        receptors_discrepant=("receptor", lambda s: " ".join(s[a.loc[s.index, "impact"] != ""])),
    ).reset_index()
    hel["_o"] = hel["canonical_helix"].map(ch.ELEMENT_RANK)
    hel = hel.sort_values("_o").drop(columns="_o")
    for df, name in [(a, "lbd_window_interval_audit_v1.csv"), (r, "lbd_window_receptor_summary_v1.csv"),
                     (hel, "lbd_window_helix_summary_v1.csv")]:
        p = output_root / OUT / name
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=False, lineterminator="\n")
    return {"window_class": a["window_class"].value_counts().to_dict(),
            "impact": a[a["impact"] != ""]["impact"].value_counts().to_dict(),
            "recommendation": r["recommendation"].value_counts().to_dict(),
            "windows_equal_uniprot_feature": bool((a["uniprot_nr_lbd_feature"] ==
                                                   a["project_lbd_start"].astype(str) + "-" + a["project_lbd_end"].astype(str)).all())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    print(json.dumps(build(ap.parse_args().output_root), indent=2))


if __name__ == "__main__":
    main()
