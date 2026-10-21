"""
Build the corrected inputs for the v2 foundation rebuild (scripts/101 --variant v2).
Nothing frozen is modified; outputs go to data/rebuild_v2/inputs/.

1. FINAL ANALYSIS SPAN per receptor (professor decision 2026-10-20):
     mapped receptors:   smallest continuous span covering
                         UniProt/PROSITE 'NR LBD' positions  UNION  every position of every
                         professor canonical helix interval (never truncating a professor helix)
     unmapped receptors: the UniProt/PROSITE 'NR LBD' span, unchanged
   TRalpha is numbered on the professor's accession P10827-2 (TRalpha1). Its
   PROSITE window (P10827 163-407) is shared with P10827-2 only for 163-370;
   P10827 371-407 is TRalpha2-specific sequence that no TRalpha1 construct
   contains, so it has no P10827-2 counterpart. TRalpha span = union(163-370,
   professor 161-406) = 161-406 in P10827-2 numbering.
   This is an analysis-span rule only: it assigns no helix identity.

2. nr_lbd_reference.csv (v1 schema) with lbd_start/lbd_end/lbd_length = span,
   TRalpha sequence_length 410 and coordinate system P10827-2.

3. rcsb_sequence_coordinate_alignments.csv (v1 schema) with aligned regions
   replaced by the sequence-evidence regions (scripts/96) for every entity that
   has coordinates; coverage columns recomputed from the new regions with the
   same formula as the RCSB values (aligned residues / length). Rows without
   coordinate evidence keep their RCSB regions and are listed.

4. TRalpha reference entry: a copy of the cached P10827 UniProt entry whose
   sequence is replaced by P10827-2 (stored FASTA). Features lying entirely in
   the shared residues 1-370 are kept; features touching 371+ (TRalpha2-specific
   numbering) are dropped and listed.

Reads:
    data/sources/professor_canonical_helix_table/NR_canonical_helix_boundaries.csv
    data/sources/professor_canonical_helix_table/uniprot_isoform_sequences/P10827-2.fasta
    data/manifests/nr_lbd_reference.csv, data/manifests/rcsb_sequence_coordinate_alignments.csv
    data/canonical_remap/entity_uniprot_regions_corrected_v1.csv   (scripts/96)
    data/raw/api/uniprot/P10827.json

Writes (data/rebuild_v2/inputs/):
    analysis_span_v2.csv, nr_lbd_reference_v2.csv, rcsb_sequence_coordinate_alignments_v2.csv,
    alignment_region_substitution_log_v2.csv, uniprot_P10827_reference_P10827-2.json,
    v2_inputs_manifest.json
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
from scripts.utils import sequence_numbering as sn  # noqa: E402
from scripts.utils.mmcif_parse import as_list, parse_mmcif_gz  # noqa: E402

SOURCE_CSV = PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table" / "NR_canonical_helix_boundaries.csv"
TRA1_FASTA = PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table" / "uniprot_isoform_sequences" / "P10827-2.fasta"
LBD_REF = PROJECT_ROOT / "data" / "manifests" / "nr_lbd_reference.csv"
ALIGN = PROJECT_ROOT / "data" / "manifests" / "rcsb_sequence_coordinate_alignments.csv"
REGIONS = PROJECT_ROOT / "data" / "canonical_remap" / "entity_uniprot_regions_corrected_v1.csv"
UNIPROT_P10827 = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot" / "P10827.json"
UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"
MMCIF_CACHES = [PROJECT_ROOT / "data" / "raw" / "mmcif", PROJECT_ROOT / "data" / "interim" / "raw_mmcif_v2_additions"]


def mmcif_path(pdb: str) -> Path:
    for c in MMCIF_CACHES:
        if (c / f"{pdb}.cif.gz").exists():
            return c / f"{pdb}.cif.gz"
    return MMCIF_CACHES[0] / f"{pdb}.cif.gz"
DOWNLOAD_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "mmcif_download_manifest.csv"
DOWNLOAD_ADDITIONS = PROJECT_ROOT / "data" / "rebuild_v2" / "inputs" / "mmcif_download_additions_v2.csv"
OUT = Path("data") / "rebuild_v2" / "inputs"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def analysis_spans(ref: pd.DataFrame, src: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in ref.sort_values("uniprot_id").itertuples():
        g = src[src["nr_code"] == r.nr_code]
        acc = ch.reference_accession(r.uniprot_id)
        lo, hi = int(r.lbd_start), int(r.lbd_end)
        note = ""
        if acc != r.uniprot_id:  # TRalpha: only the shared-numbering part of the PROSITE window exists in P10827-2
            hi_shared = min(hi, ch.TRA_IDENTICAL_END)
            note = (f"PROSITE window {lo}-{hi} is on P10827; only {lo}-{hi_shared} is shared with P10827-2 "
                    f"(P10827 {hi_shared + 1}-{hi} is TRalpha2-specific, no TRalpha1 counterpart)")
            hi = hi_shared
        if len(g):
            pmin, pmax = int(g["start_uniprot"].astype(int).min()), int(g["end_uniprot"].astype(int).max())
            s, e = min(lo, pmin), max(hi, pmax)
            rule = "PROSITE_NR_LBD_UNION_PROFESSOR_HELICES"
        else:
            pmin = pmax = None
            s, e = lo, hi
            rule = "PROSITE_NR_LBD_ONLY_NO_PROFESSOR_MAPPING"
        rows.append({
            "uniprot_id": r.uniprot_id, "nr_code": r.nr_code, "common_name": r.common_name,
            "reference_accession": acc, "prosite_lbd_start_v1": int(r.lbd_start), "prosite_lbd_end_v1": int(r.lbd_end),
            "professor_min_start": pmin, "professor_max_end": pmax, "professor_mapping": bool(len(g)),
            "analysis_span_start": s, "analysis_span_end": e, "analysis_span_length": e - s + 1,
            "n_positions_added_n_terminal": max(0, int(r.lbd_start) - s),
            "n_positions_added_c_terminal": max(0, e - int(r.lbd_end)), "span_rule": rule, "note": note,
        })
    return pd.DataFrame(rows)


def reference_sequence(uid: str) -> str:
    ref = ch.reference_accession(uid)
    if ref != uid:
        return ch.read_fasta_sequence(TRA1_FASTA)
    return json.loads((UNIPROT_DIR / f"{uid}.json").read_text())["sequence"]["value"]


def entity_alignment_regions(pdb: str, entity_id: str, uid: str) -> tuple[list[dict], dict]:
    """Sequence-evidence regions for an entity the crosswalk did not cover (same
    rules as scripts/96: consensus of co-optimal alignments, local support,
    depositor tag/linker residues withheld)."""
    d = parse_mmcif_gz(mmcif_path(pdb))
    asym_ent = dict(zip(as_list(d["_struct_asym.id"]), as_list(d["_struct_asym.entity_id"])))
    asym = sorted(a for a, e in asym_ent.items() if e == entity_id)[0]
    strand = dict(zip(as_list(d["_pdbx_poly_seq_scheme.asym_id"]), as_list(d["_pdbx_poly_seq_scheme.pdb_strand_id"])))[asym]
    scheme = sn.chain_scheme(d, asym)
    q = "".join(sn.one_letter(r["mon_id"], sn.modified_residue_parents(d)) for r in scheme)
    useq = reference_sequence(uid)
    res = sn.align_unique(q, useq)
    dif, _ = sn.seq_dif_for_strand(d, strand, uid)
    non_ref = {i for i, r in enumerate(scheme) if sn.is_non_reference_residue(dif.get(r["label_seq_id"], ""))}
    calls, summ = sn.call_residues(q, useq, res, non_ref)
    regions, cur = [], None
    for r, c in zip(scheme, calls):
        if c.target_index is None:
            cur = None
            continue
        lsid, pos = r["label_seq_id"], c.target_index + 1
        if cur and lsid == cur["target_end"] + 1 and pos == cur["query_end"] + 1:
            cur["target_end"], cur["query_end"] = lsid, pos
        else:
            cur = {"query_begin": pos, "query_end": pos, "target_begin": lsid, "target_end": lsid}
            regions.append(cur)
    return regions, summ


def build(output_root: Path) -> dict:
    out = output_root / OUT
    out.mkdir(parents=True, exist_ok=True)
    ref = pd.read_csv(LBD_REF, keep_default_na=False)
    src = pd.read_csv(SOURCE_CSV, encoding=ch.SOURCE_ENCODING, dtype=str, keep_default_na=False)
    tra1 = ch.read_fasta_sequence(TRA1_FASTA)

    span = analysis_spans(ref, src)
    span.to_csv(out / "analysis_span_v2.csv", index=False, lineterminator="\n")

    ref2 = ref.copy()
    sp = span.set_index("uniprot_id")
    for i, r in ref2.iterrows():
        s = sp.loc[r["uniprot_id"]]
        ref2.at[i, "lbd_start"] = int(s["analysis_span_start"])
        ref2.at[i, "lbd_end"] = int(s["analysis_span_end"])
        ref2.at[i, "lbd_length"] = int(s["analysis_span_length"])
        ref2.at[i, "lbd_definition"] = f"ANALYSIS_SPAN_V2_{s['span_rule']}"
        ref2.at[i, "lbd_notes"] = (f"v2 analysis span (professor decision 2026-10-20); v1 UniProt NR LBD "
                                   f"{s['prosite_lbd_start_v1']}-{s['prosite_lbd_end_v1']}" + (f"; {s['note']}" if s["note"] else ""))
        if s["reference_accession"] != r["uniprot_id"]:
            ref2.at[i, "sequence_length"] = len(tra1)
            ref2.at[i, "lbd_coordinate_system"] = f"UNIPROT_ISOFORM_{s['reference_accession']}_1_BASED_INCLUSIVE"
    ref2.to_csv(out / "nr_lbd_reference_v2.csv", index=False, lineterminator="\n")

    al = pd.read_csv(ALIGN, keep_default_na=False, dtype=str)
    reg = pd.read_csv(REGIONS, dtype=str).set_index(["uniprot_id", "polymer_entity_id"])
    log = []
    for i, r in al.iterrows():
        key = (r["uniprot_id"], r["polymer_entity_id"])
        qlen = len(tra1) if ch.reference_accession(r["uniprot_id"]) != r["uniprot_id"] else int(r["query_length"])
        cached = mmcif_path(r["pdb_id"]).exists()
        extra = None
        if key not in reg.index and cached:
            regions, summ = entity_alignment_regions(r["pdb_id"], str(r["entity_id"]), r["uniprot_id"])
            extra = (json.dumps(regions), summ)
        if key in reg.index or (extra and extra[1]["sufficient_evidence"] and json.loads(extra[0])):
            new = reg.loc[key, "aligned_regions_json"] if key in reg.index else extra[0]
            n = sum(x["query_end"] - x["query_begin"] + 1 for x in json.loads(new))
            changed = json.loads(new) != json.loads(r["aligned_regions_json"])
            al.at[i, "aligned_regions_json"] = new
            al.at[i, "rcsb_reference_sequence_coverage"] = repr(n / qlen)
            al.at[i, "rcsb_entity_sequence_coverage"] = repr(n / int(r["target_length"]))
            source = "SEQUENCE_EVIDENCE_V1_CROSSWALK" if key in reg.index else "SEQUENCE_EVIDENCE_ENTITY_ALIGNMENT_V2"
        elif extra:
            changed, source = False, "RCSB_SERVICE_KEPT_ALIGNMENT_INSUFFICIENT_EVIDENCE"
        else:
            changed, source = False, "RCSB_SERVICE_UNVERIFIED_NO_LOCAL_MMCIF"
        al.at[i, "query_length"] = str(qlen)
        log.append({"uniprot_id": r["uniprot_id"], "polymer_entity_id": r["polymer_entity_id"],
                    "region_source": source, "regions_changed": changed})
    al.to_csv(out / "rcsb_sequence_coordinate_alignments_v2.csv", index=False, lineterminator="\n")
    lg = pd.DataFrame(log)
    lg.to_csv(out / "alignment_region_substitution_log_v2.csv", index=False, lineterminator="\n")

    d = json.loads(UNIPROT_P10827.read_text())
    kept, dropped = [], []
    for f in d.get("features", []):
        loc = f.get("location", {})
        try:
            e = int(loc["end"]["value"])
        except (KeyError, TypeError, ValueError):
            e = None
        (kept if e is not None and e <= ch.TRA_IDENTICAL_END else dropped).append(f)
    d["features"] = kept
    d["sequence"] = {"value": tra1, "length": len(tra1)}
    d["v2_reference_note"] = ("Sequence replaced by UniProt isoform P10827-2 (TRalpha1), the accession of the "
                              "professor canonical helix table. Features touching residue >370 (TRalpha2-specific "
                              f"numbering) dropped: {len(dropped)}; kept: {len(kept)}.")
    (out / "uniprot_P10827_reference_P10827-2.json").write_text(json.dumps(d, indent=1, ensure_ascii=False))

    dm = pd.read_csv(DOWNLOAD_MANIFEST, dtype=str, keep_default_na=False)
    add = pd.read_csv(DOWNLOAD_ADDITIONS, dtype=str, keep_default_na=False)[list(dm.columns)]
    assert not set(add["pdb_id"]) & set(dm["pdb_id"])
    dm2 = pd.concat([dm, add], ignore_index=True).sort_values("pdb_id")
    dm2.to_csv(out / "mmcif_download_manifest_v2.csv", index=False, lineterminator="\n")

    subs = {
        "data/manifests/mmcif_download_manifest.csv": str(OUT / "mmcif_download_manifest_v2.csv"),
        "data/manifests/nr_lbd_reference.csv": str(OUT / "nr_lbd_reference_v2.csv"),
        "data/manifests/rcsb_sequence_coordinate_alignments.csv": str(OUT / "rcsb_sequence_coordinate_alignments_v2.csv"),
        "data/raw/api/uniprot/P10827.json": str(OUT / "uniprot_P10827_reference_P10827-2.json"),
    }
    manifest = {"description": "v2 rebuild input substitutions (sandbox path -> versioned source)",
                "substitutions": subs,
                "sha256": {k: sha256(output_root / v) for k, v in subs.items()},
                "inputs_sha256": {str(p.relative_to(PROJECT_ROOT)): sha256(p) for p in
                                  [SOURCE_CSV, TRA1_FASTA, LBD_REF, ALIGN, REGIONS, DOWNLOAD_MANIFEST, DOWNLOAD_ADDITIONS]}}
    (out / "v2_inputs_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return {"span_changed_receptors": int((span["analysis_span_length"] != (span["prosite_lbd_end_v1"] - span["prosite_lbd_start_v1"] + 1)).sum()),
            "regions_substituted": int(lg["region_source"].str.startswith("SEQUENCE_EVIDENCE").sum()),
            "region_source_counts": lg["region_source"].value_counts().to_dict(),
            "regions_changed": int(lg["regions_changed"].sum()),
            "rows_unverified": int((~lg["region_source"].str.startswith("SEQUENCE_EVIDENCE")).sum()),
            "tra_features_dropped": len(dropped)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    print(json.dumps(build(ap.parse_args().output_root), indent=2))


if __name__ == "__main__":
    main()
