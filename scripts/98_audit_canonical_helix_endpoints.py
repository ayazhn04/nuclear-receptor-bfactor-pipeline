"""
Canonical helix remapping, step 3: independent numbering check of the
professor helix table against this project's own residue -> UniProt mapping.

For intervals whose stated basis is "deposited segment endpoints", the
professor took boundaries from the selected reference chain's deposited mmCIF
helix annotations (_struct_conf) and converted them to UniProt numbering by
their own sequence alignment. Here the same deposited helix records are
converted with THIS project's mapping (entity_label_seq_id ->
canonical_uniprot_position in data/processed/lbd_residue_bfactors_all.parquet,
the numbering every B-factor profile uses). If both numberings agree, the
professor's start/end positions fall exactly on mapped deposited helix
begin/end positions; a numbering offset would show up as a systematic shift.

As a second, project-independent numbering, each reference chain's own
polymer sequence (_entity_poly_seq) is aligned to the UniProt sequence of the
source accession (Biopython PairwiseAligner, BLOSUM62, local) -- the method
the professor's document describes -- and the same endpoint check is
repeated. Agreement under the alignment numbering but not under the project
numbering points at the project mapping, not the source (see
scripts/96_build_residue_numbering_crosswalk.py).

Read-only. Needs the locally cached mmCIF files (data/raw/mmcif/) and cached
UniProt entries (data/raw/api/uniprot/), both gitignored.

Reads:
    data/canonical_remap/nr_lbd_canonical_helix_intervals_professor_2026.csv  (scripts/97 output)
    data/processed/lbd_residue_bfactors_all.parquet
    data/raw/mmcif/<PDB>.cif.gz
    data/raw/api/uniprot/<accession>.json
    data/sources/professor_canonical_helix_table/uniprot_isoform_sequences/P10827-2.fasta

Writes:
    reports/canonical_remap/canonical_helix_endpoint_concordance.csv
    reports/canonical_remap/canonical_helix_endpoint_concordance_by_receptor.csv

Usage:
    python scripts/98_audit_canonical_helix_endpoints.py [--mmcif-dir DIR] [--output-root DIR]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json  # noqa: E402

from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402

from scripts.utils import canonical_helix as ch  # noqa: E402
from scripts.utils.mmcif_parse import as_list, get_secondary_structure_ranges, parse_mmcif_gz  # noqa: E402

INTERVALS_CSV = PROJECT_ROOT / "data" / "canonical_remap" / "nr_lbd_canonical_helix_intervals_professor_2026.csv"
RESIDUES_PARQUET = PROJECT_ROOT / "data" / "processed" / "lbd_residue_bfactors_all.parquet"
DEFAULT_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"
DEFAULT_UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"
ISOFORM_DIR = PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table" / "uniprot_isoform_sequences"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def source_sequence(accession: str, uniprot_dir: Path) -> str:
    iso = ISOFORM_DIR / f"{accession}.fasta"
    if iso.exists():
        return ch.read_fasta_sequence(iso)
    return json.loads((uniprot_dir / f"{accession}.json").read_text())["sequence"]["value"]


def alignment_label_to_uniprot(mmcif_dict: dict, label_asym: str, uniprot_seq: str) -> tuple[dict[int, int], float]:
    """label_seq_id -> UniProt position from a local alignment of the chain's
    own polymer sequence to the UniProt sequence (identical residues only)."""
    asym_ids = as_list(mmcif_dict["_struct_asym.id"])
    entity = dict(zip(asym_ids, as_list(mmcif_dict["_struct_asym.entity_id"])))[label_asym]
    ent = as_list(mmcif_dict["_entity_poly_seq.entity_id"])
    nums = as_list(mmcif_dict["_entity_poly_seq.num"])
    mons = as_list(mmcif_dict["_entity_poly_seq.mon_id"])
    seen: dict[int, str] = {}
    for e, n, mon in zip(ent, nums, mons):
        if e == entity and int(n) not in seen:  # first listed residue at a heterogeneous position
            seen[int(n)] = THREE_TO_ONE.get(mon, "X")
    label_ids = sorted(seen)
    chain_seq = "".join(seen[i] for i in label_ids)
    aligner = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                    open_gap_score=-10, extend_gap_score=-0.5)
    aln = aligner.align(chain_seq, uniprot_seq)[0]
    mapping, n_aligned, n_ident = {}, 0, 0
    for (qs, qe), (ts, te) in zip(*aln.aligned):
        for q, t in zip(range(qs, qe), range(ts, te)):
            n_aligned += 1
            if chain_seq[q] == uniprot_seq[t]:
                n_ident += 1
                mapping[label_ids[q]] = t + 1
    return mapping, (round(n_ident / n_aligned, 4) if n_aligned else 0.0)

OUT_INTERVALS = Path("reports") / "canonical_remap" / "canonical_helix_endpoint_concordance.csv"
OUT_RECEPTORS = Path("reports") / "canonical_remap" / "canonical_helix_endpoint_concordance_by_receptor.csv"

DEPOSITED_BASIS_MARKER = "deposited segment endpoints"


def endpoint_result(pos: int, mapped: set[int], frame: set[int]) -> tuple[str, object]:
    """EXACT / MISMATCH(+nearest delta) / NOT_CHECKABLE (position outside the
    mapped frame of the reference chain under that numbering)."""
    if pos not in frame:
        return "NOT_CHECKABLE_OUTSIDE_MAPPED_FRAME", pd.NA
    if pos in mapped:
        return "EXACT", 0
    if not mapped:
        return "NO_MAPPED_DEPOSITED_HELIX", pd.NA
    nearest = min(mapped, key=lambda m: (abs(pos - m), m))
    return "MISMATCH", pos - nearest


def label_asym_for_auth_chain(mmcif_dict: dict, auth_chain: str) -> list[str]:
    asym = as_list(mmcif_dict["_pdbx_poly_seq_scheme.asym_id"])
    strand = as_list(mmcif_dict["_pdbx_poly_seq_scheme.pdb_strand_id"])
    return sorted({a for a, s in zip(asym, strand) if s == auth_chain})


def summarize_checks(exp: pd.DataFrame, start_col: str, end_col: str, sd_col: str, ed_col: str) -> dict:
    checks = pd.concat([exp[start_col], exp[end_col]])
    deltas = [d for d in pd.concat([exp[sd_col], exp[ed_col]]) if pd.notna(d) and d != 0]
    n_chk = int(checks.isin(["EXACT", "MISMATCH"]).sum())
    n_exact = int((checks == "EXACT").sum())
    common = Counter(deltas).most_common(1)
    return {"checkable": n_chk, "exact": n_exact, "mismatch": int((checks == "MISMATCH").sum()),
            "not_checkable": int(len(checks) - n_chk),
            "fraction_exact": round(n_exact / n_chk, 4) if n_chk else pd.NA,
            "common_delta": common[0][0] if common else pd.NA, "common_delta_count": common[0][1] if common else 0}


def build(mmcif_dir: Path, output_root: Path, uniprot_dir: Path = DEFAULT_UNIPROT_DIR) -> pd.DataFrame:
    iv = pd.read_csv(INTERVALS_CSV, keep_default_na=False, dtype={"start": int, "end": int})
    res = pd.read_parquet(RESIDUES_PARQUET, columns=[
        "uniprot_id", "pdb_id", "auth_asym_id", "label_asym_id", "model_num", "canonical_uniprot_position",
        "entity_label_seq_id", "is_construct_mapped"])
    res = res[res["is_construct_mapped"]]

    rows = []
    for (uid, acc, pdb, chain), g in iv.groupby(
            ["project_uniprot_id", "sequence_accession", "reference_pdb", "reference_chain"], sort=True):
        path = mmcif_dir / f"{pdb}.cif.gz"
        if path.exists():
            d = parse_mmcif_gz(path)
            label_asyms = label_asym_for_auth_chain(d, chain)
            if len(label_asyms) != 1:
                raise RuntimeError(f"{pdb}: auth chain {chain} -> label asym {label_asyms}")
            la = label_asyms[0]
            helices = [h for h in get_secondary_structure_ranges(d)
                       if h["sse_type"] == "HELIX" and h["label_asym_id"] == la]
        else:  # e.g. NOR-1's NMR reference: never downloaded by this X-ray B-factor project
            d, la, helices = None, "", []

        # (1) project numbering
        r = res[(res["uniprot_id"] == uid) & (res["pdb_id"] == pdb) & (res["label_asym_id"] == la)]
        if len(r):
            model = sorted(r["model_num"].unique(), key=lambda m: int(m))[0]
            r = r[r["model_num"] == model]
        chain_status = ("OK" if len(r) else "REFERENCE_CHAIN_NOT_IN_PROJECT_RESIDUE_TABLE") if d is not None \
            else "REFERENCE_MMCIF_NOT_IN_LOCAL_CACHE"
        l2u = dict(zip(r["entity_label_seq_id"].astype(int), r["canonical_uniprot_position"].astype(int)))
        if len(set(l2u.values())) != len(l2u):
            raise RuntimeError(f"{pdb}/{chain}: label_seq -> UniProt mapping is not one-to-one")
        p_beg = {l2u[h["beg_label_seq_id"]] for h in helices if h["beg_label_seq_id"] in l2u}
        p_end = {l2u[h["end_label_seq_id"]] for h in helices if h["end_label_seq_id"] in l2u}
        p_frame = set(l2u.values())

        # (2) independent alignment numbering against the source accession
        a2u, aln_identity = (alignment_label_to_uniprot(d, la, source_sequence(acc, uniprot_dir))
                             if d is not None else ({}, pd.NA))
        a_beg = {a2u[h["beg_label_seq_id"]] for h in helices if h["beg_label_seq_id"] in a2u}
        a_end = {a2u[h["end_label_seq_id"]] for h in helices if h["end_label_seq_id"] in a2u}
        a_frame = set(a2u.values())

        for t in g.sort_values("source_row_number").itertuples():
            s_res, s_d = endpoint_result(t.start, p_beg, p_frame)
            e_res, e_d = endpoint_result(t.end, p_end, p_frame)
            sa_res, sa_d = endpoint_result(t.start, a_beg, a_frame)
            ea_res, ea_d = endpoint_result(t.end, a_end, a_frame)
            rows.append({
                "source_row_number": t.source_row_number, "nr_code": t.nr_code, "project_uniprot_id": uid,
                "receptor": t.receptor, "sequence_accession": acc,
                "canonical_helix": t.canonical_helix, "start": t.start, "end": t.end,
                "assignment_status": t.assignment_status, "boundary_basis": t.boundary_basis,
                "basis_expects_deposited_endpoints": DEPOSITED_BASIS_MARKER in t.boundary_basis,
                "reference_pdb": pdb, "reference_chain_auth": chain, "reference_label_asym_id": la,
                "reference_chain_status": chain_status, "n_deposited_helix_records_in_chain": len(helices),
                "start_check": s_res, "start_delta_to_nearest_deposited_begin": s_d,
                "end_check": e_res, "end_delta_to_nearest_deposited_end": e_d,
                "alignment_identity_chain_vs_source_sequence": aln_identity,
                "start_check_alignment": sa_res, "start_delta_alignment": sa_d,
                "end_check_alignment": ea_res, "end_delta_alignment": ea_d,
            })
    out = pd.DataFrame(rows).sort_values("source_row_number").reset_index(drop=True)

    summ = []
    for code, g in out.groupby("nr_code", sort=True):
        exp = g[g["basis_expects_deposited_endpoints"]]
        p = summarize_checks(exp, "start_check", "end_check", "start_delta_to_nearest_deposited_begin",
                             "end_delta_to_nearest_deposited_end")
        a = summarize_checks(exp, "start_check_alignment", "end_check_alignment", "start_delta_alignment",
                             "end_delta_alignment")
        summ.append({
            "nr_code": code, "project_uniprot_id": g["project_uniprot_id"].iloc[0], "receptor": g["receptor"].iloc[0],
            "reference_pdb": g["reference_pdb"].iloc[0], "reference_chain_auth": g["reference_chain_auth"].iloc[0],
            "reference_chain_status": g["reference_chain_status"].iloc[0],
            "n_intervals": len(g), "n_intervals_expecting_deposited_endpoints": len(exp),
            "n_endpoints_checkable": p["checkable"], "n_endpoints_exact": p["exact"],
            "n_endpoints_mismatch": p["mismatch"], "n_endpoints_not_checkable": p["not_checkable"],
            "fraction_exact_of_checkable": p["fraction_exact"],
            "most_common_nonzero_delta": p["common_delta"], "most_common_nonzero_delta_count": p["common_delta_count"],
            "alignment_identity_chain_vs_source_sequence": g["alignment_identity_chain_vs_source_sequence"].iloc[0],
            "n_endpoints_checkable_alignment": a["checkable"], "n_endpoints_exact_alignment": a["exact"],
            "n_endpoints_mismatch_alignment": a["mismatch"],
            "fraction_exact_of_checkable_alignment": a["fraction_exact"],
            "most_common_nonzero_delta_alignment": a["common_delta"],
            "most_common_nonzero_delta_alignment_count": a["common_delta_count"],
        })
    summ = pd.DataFrame(summ)
    for df, rel in [(out, OUT_INTERVALS), (summ, OUT_RECEPTORS)]:
        path = output_root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False, lineterminator="\n")
    return summ


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--mmcif-dir", type=Path, default=DEFAULT_MMCIF_DIR)
    ap.add_argument("--uniprot-dir", type=Path, default=DEFAULT_UNIPROT_DIR)
    ap.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    args = ap.parse_args()
    summ = build(args.mmcif_dir, args.output_root, args.uniprot_dir)
    with pd.option_context("display.width", 250, "display.max_rows", 100):
        print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
