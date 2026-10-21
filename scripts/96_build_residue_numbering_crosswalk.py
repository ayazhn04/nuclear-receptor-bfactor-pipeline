"""
Residue-numbering correction layer v1: rebuild the PDB -> UniProt residue
numbering of EVERY measured chain instance from primary sequence evidence and
record it, residue by residue, next to the frozen (RCSB-service) numbering.

The frozen tables are NOT modified. This builds a new, versioned crosswalk:
unchanged residues are recorded as UNCHANGED, corrected ones carry the old
and new position plus the evidence, and anything that cannot be placed
unambiguously is WITHHELD (never forced). See scripts/utils/sequence_numbering.py
for the evidence hierarchy and decision rules. B-factor values are never read.

Reads:
    data/processed/lbd_residue_bfactors_all.parquet          (frozen; old numbering + observed residues)
    data/manifests/rcsb_sequence_coordinate_alignments.csv   (frozen; old RCSB-service regions, comparison only)
    data/raw/mmcif/<PDB>.cif.gz                              (gitignored cache; primary sequence evidence)
    data/raw/api/uniprot/<accession>.json                    (gitignored cache; canonical sequences + isoform features)

Writes:
    data/canonical_remap/residue_numbering_crosswalk_v1.parquet
    reports/canonical_remap/residue_numbering_instance_audit_v1.csv
    reports/canonical_remap/residue_numbering_receptor_summary_v1.csv
    reports/canonical_remap/residue_numbering_summary_v1.json
    data/canonical_remap/entity_uniprot_regions_corrected_v1.csv   (rebuild input: same aligned_regions_json
        schema as data/manifests/rcsb_sequence_coordinate_alignments.csv, query = UniProt, target = entity)

Usage:
    python scripts/96_build_residue_numbering_crosswalk.py [--output-root DIR] [--pdb ID ...] [--workers N]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils import canonical_helix as ch  # noqa: E402
from scripts.utils import sequence_numbering as sn  # noqa: E402
from scripts.utils.mmcif_parse import parse_mmcif_gz  # noqa: E402
from scripts.utils.residue_mapping import build_entity_to_uniprot_map  # noqa: E402

RESIDUES_PARQUET = PROJECT_ROOT / "data" / "processed" / "lbd_residue_bfactors_all.parquet"
RCSB_ALIGNMENTS_CSV = PROJECT_ROOT / "data" / "manifests" / "rcsb_sequence_coordinate_alignments.csv"
MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"
UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"
ISOFORM_DIR = PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table" / "uniprot_isoform_sequences"

CROSSWALK_OUT = Path("data") / "canonical_remap" / "residue_numbering_crosswalk_v1.parquet"
INSTANCE_OUT = Path("reports") / "canonical_remap" / "residue_numbering_instance_audit_v1.csv"
RECEPTOR_OUT = Path("reports") / "canonical_remap" / "residue_numbering_receptor_summary_v1.csv"
SUMMARY_OUT = Path("reports") / "canonical_remap" / "residue_numbering_summary_v1.json"
REGIONS_OUT = Path("data") / "canonical_remap" / "entity_uniprot_regions_corrected_v1.csv"

MAPPING_METHOD = (f"local alignment of _pdbx_poly_seq_scheme chain sequence vs the reference UniProt sequence "
                  f"(canonical; TRalpha: P10827-2) "
                  f"(BLOSUM62, gap {sn.GAP_OPEN}/{sn.GAP_EXTEND}); consensus of all co-optimal alignments; "
                  f"identical residues need >= {sn.MIN_EDGE_IDENTITY} identity over {sn.WINDOW} aligned pairs on one side "
                  f"or >= {sn.MIN_WINDOW_IDENTITY} on both, "
                  f"substitutions >= {sn.MIN_WINDOW_IDENTITY} on both sides; depositor tag/linker residues withheld")
OLD_COLS = ["uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id", "entity_id", "instance_id",
            "label_asym_id", "auth_asym_id", "model_num", "canonical_uniprot_position", "entity_label_seq_id",
            "is_construct_mapped", "ca_record_present", "label_comp_ids_json", "canonical_residue_name"]


def uniprot_entry(acc: str) -> tuple[str, list[tuple[int, int]]]:
    """Reference sequence for one project receptor under the FINAL numbering
    policy (TRalpha -> P10827-2, the professor's accession; stored FASTA, no
    isoform features), else the cached canonical UniProt entry."""
    ref = ch.reference_accession(acc)
    if ref != acc:
        return ch.read_fasta_sequence(ISOFORM_DIR / f"{ref}.fasta"), []
    d = json.loads((UNIPROT_DIR / f"{acc}.json").read_text())
    missing = []  # isoform "Missing" segments (alternative isoforms lacking these residues)
    for f in d.get("features", []):
        if f.get("type") == "Alternative sequence" and not f.get("alternativeSequence", {}).get("alternativeSequences"):
            missing.append((int(f["location"]["start"]["value"]), int(f["location"]["end"]["value"])))
    varseq_starts = sorted(int(f["location"]["start"]["value"]) for f in d.get("features", [])
                           if f.get("type") == "Alternative sequence")
    return d["sequence"]["value"], sorted(set(missing)) + [(s, -1) for s in varseq_starts]


def first_json(s):
    if s is None or s in ("", "[]"):
        return None
    v = json.loads(s)
    return v[0] if v else None


def process_pdb(args) -> tuple[list[dict], list[dict]]:
    pdb, instances, mmcif_dir, uniprot_dir = args
    global UNIPROT_DIR
    UNIPROT_DIR = Path(uniprot_dir)
    d = parse_mmcif_gz(Path(mmcif_dir) / f"{pdb}.cif.gz")
    parents = sn.modified_residue_parents(d)
    rows, inst_rows = [], []
    for inst in instances:
        uid, asym, strand = inst["uniprot_id"], inst["label_asym_id"], inst["auth_asym_id"]
        useq, iso = uniprot_entry(uid)
        iso_missing = [(a, b) for a, b in iso if b != -1]
        varseq_starts = [a for a, b in iso if b == -1]
        scheme = sn.chain_scheme(d, asym)
        query = "".join(sn.one_letter(r["mon_id"], parents) for r in scheme)
        res = sn.align_unique(query, useq)
        dif, declared_del = sn.seq_dif_for_strand(d, strand, uid)
        non_ref = {i for i, r in enumerate(scheme) if sn.is_non_reference_residue(dif.get(r["label_seq_id"], ""))}
        calls, summ = sn.call_residues(query, useq, res, non_ref)
        gaps = sn.gaps_in_alignment(res.blocks)
        segs = sn.depositor_segments(d, strand, uid)
        seg_gaps = {p for (a, b), (c, e) in zip([(s[2], s[3]) for s in segs], [(s[2], s[3]) for s in segs[1:]])
                    for p in range(b + 1, c)}
        old_map = inst["old_map"]
        old_tbl = inst["old_table"]  # label_seq -> (old_pos, ca_present, observed comp)
        gap_desc = []
        for g in gaps:
            if g["type"] == "DELETION_IN_CHAIN":
                lo, hi = g["target_from"] + 1, g["target_to"] + 1
                isoform = sn.isoform_missing_match(lo, hi, iso_missing)
                declared = set(range(lo, hi + 1)) & (declared_del | seg_gaps)
                g["desc"] = (f"DELETION_IN_CHAIN UniProt {lo}-{hi} ({hi - lo + 1} aa)"
                             + ("; matches UniProt isoform 'Missing' segment" if isoform else "")
                             + (f"; depositor-declared {len(declared)}/{hi - lo + 1}" if declared else "; not depositor-declared"))
                g["isoform"] = isoform
            else:
                lo, hi = g["query_from"], g["query_to"]
                g["desc"] = f"INSERTION_IN_CHAIN label_seq {scheme[lo]['label_seq_id']}-{scheme[hi]['label_seq_id']} ({hi - lo + 1} aa)"
                g["isoform"] = False
            gap_desc.append(g["desc"])
        aligned_t = sorted(res.mapping.values())
        tail_isoform = bool(aligned_t) and any(
            s - 1 <= aligned_t[-1] + 1 <= s + 1 for s in varseq_starts) and len(query) - 1 > max(res.mapping)

        counts = Counter()
        offsets = Counter()
        seg_offsets = []
        for i, (r, c) in enumerate(zip(scheme, calls)):
            lsid = r["label_seq_id"]
            old = old_map.get(lsid)
            new = c.target_index + 1 if c.target_index is not None else None
            if old is None and new is None:
                continue
            in_tbl = lsid in old_tbl
            t_old, ca, obs_comp = old_tbl.get(lsid, (None, None, None))
            if in_tbl and t_old != old:
                raise RuntimeError(f"{inst['instance_id']} label_seq {lsid}: frozen table {t_old} != RCSB region {old}")
            if old is not None and new is not None:
                status = "UNCHANGED" if old == new else "CORRECTED"
            elif old is not None:
                status = "WITHHELD_PREVIOUSLY_MAPPED"
            else:
                status = "NEWLY_MAPPED"
            counts[status] += 1
            ures = useq[new - 1] if new else ""
            relation = ""
            if new:
                relation = ("IDENTICAL" if query[i] == ures else
                            "NONSTANDARD_UNRESOLVED" if query[i] == "X" else "SUBSTITUTION")
                if relation == "SUBSTITUTION":
                    counts["substitution_documented" if lsid in dif else "substitution_undocumented"] += 1
            upstream = [g["desc"] for g in gaps if g["after_query"] < i]
            if status == "CORRECTED":
                offsets[new - old] += 1
                reason = (f"offset {new - old:+d}; " + ("; ".join(upstream) if upstream else
                          "previous RCSB-service region mis-registered (no alignment gap upstream)"))
            elif status == "WITHHELD_PREVIOUSLY_MAPPED":
                reason = c.reason
            elif status == "NEWLY_MAPPED":
                reason = "position not covered by previous RCSB-service region"
            else:
                reason = ""
            if old is not None:
                key = (new - old) if new is not None else None
                if not seg_offsets or seg_offsets[-1][0] != key:
                    seg_offsets.append([key, lsid, lsid, old])
                else:
                    seg_offsets[-1][2] = lsid
            rows.append({
                "uniprot_id": uid, "nr_code": inst["nr_code"], "common_name": inst["common_name"],
                "uniprot_sequence_accession": ch.reference_accession(uid), "pdb_id": pdb, "entity_id": inst["entity_id"],
                "polymer_entity_id": inst["polymer_entity_id"], "instance_id": inst["instance_id"],
                "label_asym_id": asym, "auth_asym_id": strand, "model_num": inst["model_num"],
                "label_seq_id": lsid, "auth_seq_num": r["auth_seq_num"], "auth_ins_code": r["pdb_ins_code"],
                "chain_mon_id": r["mon_id"], "chain_residue": query[i], "modelled_in_poly_seq_scheme": r["modelled"],
                "microheterogeneous": r["heterogeneous"], "in_frozen_lbd_table": in_tbl,
                "frozen_ca_record_present": ca, "frozen_observed_comp_id": obs_comp or "",
                "old_canonical_uniprot_position": old, "corrected_uniprot_position": new,
                "uniprot_residue_at_old": useq[old - 1] if old and old <= len(useq) else "",
                "uniprot_residue_at_corrected": ures, "residue_relation": relation,
                "depositor_seq_dif": dif.get(lsid, ""), "mapping_method": MAPPING_METHOD,
                "alignment_window_identity_left": round(c.window_left, 4) if c.window_left is not None else None,
                "alignment_window_identity_right": round(c.window_right, 4) if c.window_right is not None else None,
                "call_reason": c.reason, "confidence": "HIGH" if new else "WITHHELD",
                "correction_status": status, "correction_reason": reason,
            })
        n_old = sum(1 for v in old_map.values() if v is not None)
        isoform_gap = any(g["isoform"] for g in gaps)
        if not summ["sufficient_evidence"]:
            cls = "INSUFFICIENT_EVIDENCE"
        elif counts["CORRECTED"] and isoform_gap:
            cls = "ISOFORM_DIFFERENCE_OFFSET"
        elif counts["CORRECTED"] and gaps:
            cls = "DELETION_INSERTION_OFFSET"
        elif counts["CORRECTED"]:
            cls = "REGISTER_OFFSET_NO_GAP"
        elif counts["WITHHELD_PREVIOUSLY_MAPPED"]:
            cls = "PARTIAL_AMBIGUOUS"
        else:
            cls = "CLEAN"
        tbl_obs = [(t, comp) for t, ca, comp in old_tbl.values() if ca and comp]
        inst_rows.append({
            "uniprot_id": uid, "nr_code": inst["nr_code"], "common_name": inst["common_name"], "pdb_id": pdb,
            "entity_id": inst["entity_id"], "instance_id": inst["instance_id"], "label_asym_id": asym,
            "auth_asym_id": strand, "model_num": inst["model_num"], "chain_polymer_length": len(query),
            "n_old_mapped_residues": n_old, "n_frozen_table_rows_construct_mapped": len(old_tbl),
            "n_frozen_table_observed_ca": len(tbl_obs),
            **summ,
            "aligned_uniprot_start": aligned_t[0] + 1 if aligned_t else None,
            "aligned_uniprot_end": aligned_t[-1] + 1 if aligned_t else None,
            "alignment_gaps": " | ".join(gap_desc),
            "depositor_segments": " ".join(f"{a}-{b}:{c}-{e}" for a, b, c, e in segs),
            "n_depositor_declared_deleted_positions": len(declared_del | seg_gaps),
            "n_unchanged": counts["UNCHANGED"], "n_corrected": counts["CORRECTED"],
            "n_withheld_previously_mapped": counts["WITHHELD_PREVIOUSLY_MAPPED"], "n_newly_mapped": counts["NEWLY_MAPPED"],
            "n_substitutions_documented": counts["substitution_documented"],
            "n_substitutions_undocumented": counts["substitution_undocumented"],
            "offsets_new_minus_old": " ".join(f"{k:+d}x{v}" for k, v in sorted(offsets.items())),
            "offset_segments_label_seq": " ".join(
                f"{a}-{b}:{'withheld' if k is None else f'{k:+d}'}" for k, a, b, _ in seg_offsets),
            "isoform_gap": isoform_gap, "isoform_specific_tail_unmapped": tail_isoform,
            "numbering_class": cls,
        })
    return rows, inst_rows


def load_instances(pdbs: list[str] | None) -> dict[str, list[dict]]:
    old = pd.read_parquet(RESIDUES_PARQUET, columns=OLD_COLS)
    if pdbs:
        old = old[old["pdb_id"].isin(pdbs)]
    regions = pd.read_csv(RCSB_ALIGNMENTS_CSV, usecols=["uniprot_id", "polymer_entity_id", "aligned_regions_json"])
    regions = {(u, p): json.loads(j) for u, p, j in regions.itertuples(index=False)}
    by_pdb: dict[str, list[dict]] = {}
    for (uid, inst_id), g in old.groupby(["uniprot_id", "instance_id"], sort=True):
        f = g.iloc[0]
        reg = regions[(uid, f["polymer_entity_id"])]
        # previous mapping exactly as scripts/19 built it (query = UniProt, target = entity)
        m = build_entity_to_uniprot_map(reg)
        if m["status"] != "OK":
            raise RuntimeError(f"{inst_id}: old RCSB region map {m['status']}")
        cm = g[g["is_construct_mapped"]]
        old_tbl = {int(l): (int(p), bool(ca), first_json(c)) for l, p, ca, c in zip(
            cm["entity_label_seq_id"], cm["canonical_uniprot_position"], cm["ca_record_present"], cm["label_comp_ids_json"])}
        by_pdb.setdefault(f["pdb_id"], []).append({
            "uniprot_id": uid, "nr_code": f["nr_code"], "common_name": f["common_name"], "instance_id": inst_id,
            "label_asym_id": f["label_asym_id"], "auth_asym_id": f["auth_asym_id"], "model_num": f["model_num"],
            "entity_id": f["entity_id"], "polymer_entity_id": f["polymer_entity_id"],
            "old_map": m["mapping"], "old_table": old_tbl})
    return by_pdb


def corrected_regions(cw: pd.DataFrame) -> pd.DataFrame:
    """Per (receptor, polymer entity): maximal runs of consecutive label_seq_id
    with a constant (UniProt - label_seq) offset, using only residues mapped
    in EVERY instance of that entity (conservative intersection)."""
    m = cw[cw["corrected_uniprot_position"].notna()]
    n_inst = cw.groupby(["uniprot_id", "polymer_entity_id"])["instance_id"].nunique()
    k = m.groupby(["uniprot_id", "polymer_entity_id", "label_seq_id"]).agg(
        pos=("corrected_uniprot_position", "first"), n=("instance_id", "nunique"), npos=("corrected_uniprot_position", "nunique"))
    if (k["npos"] > 1).any():
        raise RuntimeError("instances of one entity disagree on a corrected position")
    k = k[k["n"] == n_inst.reindex(k.index.droplevel(2)).to_numpy()].reset_index()
    rows = []
    for (uid, pe), g in k.sort_values(["uniprot_id", "polymer_entity_id", "label_seq_id"]).groupby(["uniprot_id", "polymer_entity_id"]):
        regions, cur = [], None
        for lsid, pos in zip(g["label_seq_id"].astype(int), g["pos"].astype(int)):
            if cur and lsid == cur["target_end"] + 1 and pos == cur["query_end"] + 1:
                cur["target_end"], cur["query_end"] = lsid, pos
            else:
                cur = {"query_begin": pos, "query_end": pos, "target_begin": lsid, "target_end": lsid}
                regions.append(cur)
        rows.append({"uniprot_id": uid, "polymer_entity_id": pe, "aligned_regions_json": json.dumps(regions),
                     "n_regions": len(regions), "n_mapped_residues": len(g),
                     "mapping_source": "scripts/96 sequence-evidence crosswalk v1"})
    return pd.DataFrame(rows)


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False), path, compression="zstd",
                   use_dictionary=True, write_statistics=False)


def build(output_root: Path, pdbs: list[str] | None = None, workers: int = 8,
          mmcif_dir: Path = MMCIF_DIR, uniprot_dir: Path = UNIPROT_DIR) -> dict:
    by_pdb = load_instances(pdbs)
    jobs = [(p, by_pdb[p], str(mmcif_dir), str(uniprot_dir)) for p in sorted(by_pdb)]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(process_pdb, jobs, chunksize=8))
    else:
        results = [process_pdb(j) for j in jobs]
    rows = [r for rr, _ in results for r in rr]
    inst = [r for _, ii in results for r in ii]
    cw = pd.DataFrame(rows).sort_values(["uniprot_id", "instance_id", "label_seq_id"]).reset_index(drop=True)
    for c in ["old_canonical_uniprot_position", "corrected_uniprot_position", "label_seq_id"]:
        cw[c] = cw[c].astype("Int64")
    cw["frozen_ca_record_present"] = cw["frozen_ca_record_present"].astype("boolean")
    ia = pd.DataFrame(inst).sort_values(["uniprot_id", "instance_id"]).reset_index(drop=True)

    rec = []
    for uid, g in ia.groupby("uniprot_id", sort=True):
        bad = g[g["numbering_class"] != "CLEAN"]
        changed = g[(g["n_corrected"] > 0) | (g["n_withheld_previously_mapped"] > 0)]
        c = cw[(cw["uniprot_id"] == uid) & cw["correction_status"].isin(["CORRECTED", "WITHHELD_PREVIOUSLY_MAPPED"])]
        rec.append({
            "uniprot_id": uid, "nr_code": g["nr_code"].iloc[0], "common_name": g["common_name"].iloc[0],
            "n_instances": len(g), "n_pdb_entries": g["pdb_id"].nunique(),
            "n_instances_by_class": " ".join(f"{k}:{v}" for k, v in sorted(Counter(g["numbering_class"]).items())),
            "n_pdb_entries_with_changed_numbering": changed["pdb_id"].nunique(),
            "fraction_pdb_entries_with_changed_numbering": round(changed["pdb_id"].nunique() / g["pdb_id"].nunique(), 4),
            "pdb_entries_with_changed_numbering": " ".join(sorted(changed["pdb_id"].unique())),
            "pdb_entries_non_clean": " ".join(sorted(bad["pdb_id"].unique())),
            "n_residue_rows_corrected": int((c["correction_status"] == "CORRECTED").sum()),
            "n_residue_rows_withheld": int((c["correction_status"] == "WITHHELD_PREVIOUSLY_MAPPED").sum()),
            "old_positions_affected_range": (f"{int(c['old_canonical_uniprot_position'].min())}-"
                                             f"{int(c['old_canonical_uniprot_position'].max())}" if len(c) else ""),
        })
    rec = pd.DataFrame(rec)

    # Task F identity statistics on residues OBSERVED in the frozen LBD table
    obs = cw[cw["in_frozen_lbd_table"] & cw["frozen_ca_record_present"].fillna(False)].copy()
    obs["obs1"] = obs["frozen_observed_comp_id"].map(lambda m: sn.one_letter(m) if m else "X")
    before_match = obs["obs1"] == obs["uniprot_residue_at_old"]
    after_mapped = obs["corrected_uniprot_position"].notna()
    after_match = after_mapped & (obs["obs1"] == obs["uniprot_residue_at_corrected"])
    stats = {
        "n_instances": int(len(ia)), "n_crosswalk_rows": int(len(cw)),
        "numbering_class_counts": dict(sorted(Counter(ia["numbering_class"]).items())),
        "correction_status_counts": dict(sorted(Counter(cw["correction_status"]).items())),
        "observed_frozen_rows": int(len(obs)),
        "before_identical": int(before_match.sum()), "before_nonidentical": int((~before_match).sum()),
        "after_mapped": int(after_mapped.sum()), "after_identical": int(after_match.sum()),
        "after_substitution_documented": int((after_mapped & (obs["residue_relation"] == "SUBSTITUTION")
                                              & (obs["depositor_seq_dif"] != "")).sum()),
        "after_substitution_undocumented": int((after_mapped & (obs["residue_relation"] == "SUBSTITUTION")
                                                & (obs["depositor_seq_dif"] == "")).sum()),
        "after_nonstandard_residue_unresolved": int((after_mapped & (obs["residue_relation"] == "NONSTANDARD_UNRESOLVED")).sum()),
        "after_modified_residue_parent_identical": int((after_mapped & ~after_match & (obs["residue_relation"] == "IDENTICAL")).sum()),
        "after_withheld": int((~after_mapped).sum()),
        "before_nonidentical_explained_by_numbering": int((~before_match & (obs["correction_status"] != "UNCHANGED")).sum()),
        "before_nonidentical_real_substitution_same_position": int((~before_match & (obs["correction_status"] == "UNCHANGED")).sum()),
        "poly_seq_vs_observed_comp_mismatch": int((obs["chain_mon_id"] != obs["frozen_observed_comp_id"]).sum()),
    }
    write_parquet(cw, output_root / CROSSWALK_OUT)
    reg = corrected_regions(cw)
    (output_root / REGIONS_OUT).parent.mkdir(parents=True, exist_ok=True)
    reg.to_csv(output_root / REGIONS_OUT, index=False, lineterminator="\n")
    stats["n_entities_with_corrected_regions"] = int(len(reg))
    stats["n_entities_multi_region"] = int((reg["n_regions"] > 1).sum())
    for df, rel in [(ia, INSTANCE_OUT), (rec, RECEPTOR_OUT)]:
        p = output_root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=False, lineterminator="\n")
    (output_root / SUMMARY_OUT).write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    ap.add_argument("--pdb", nargs="*")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    print(json.dumps(build(args.output_root, args.pdb, args.workers), indent=2))


if __name__ == "__main__":
    main()
