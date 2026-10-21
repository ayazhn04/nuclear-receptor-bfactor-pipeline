"""
Phase 1D-1F of the v2 foundation: apply the professor canonical helix map to the
rebuilt (v2) foundation, build canonical-element profiles with the UNCHANGED
legacy aggregation code, and reconcile counts.

1D  professor map -> legacy segment-map schema (data/rebuild_v2/analysis/
    nr_lbd_structural_segment_map.csv). In v2 this file holds PROFESSOR canonical
    elements only: structural_segment = professor label (both labels at a position
    inside two professor intervals -> one row per label), '' between helices
    (structural_class INTER_HELIX_DERIVED), and NO_PROFESSOR_MAPPING for the six
    receptors without a professor map. No v3 label is read or written.
    Receptor-position profiles (v2) joined to the professor annotation.
1E  canonical-element profiles: scripts/36_build_segment_features.py run UNCHANGED
    inside the v2 sandbox on that map (median/mean/q25/q75/IQR, coverage, PDB
    support, apo/holo descriptive values; missing stays NaN).
1F  count reconciliation v1 -> v2, normalization edge-case audit, and the report
    for the VDR DBD entries that the final span admitted (1KB2, 1KB4, 1KB6, 1YNW).

Reads:  data/canonical_remap/nr_lbd_canonical_helix_map_professor_2026.csv (scripts/97)
        data/rebuild_v2/{processed,analysis,manifests}/...   (scripts/101 export)
        data/rebuild_v2/inputs/analysis_span_v2.csv, frozen v1 counterparts (read-only)
Writes: data/rebuild_v2/analysis/nr_lbd_structural_segment_map.csv
        data/rebuild_v2/analysis/{pdb_segment_*, receptor_segment_*}.parquet     (via scripts/36)
        data/rebuild_v2/analysis/receptor_position_professor_profile_v2.csv
        reports/rebuild_v2/count_reconciliation_v2.csv
        reports/rebuild_v2/normalization_edge_case_audit_v2.csv
        reports/rebuild_v2/vdr_dbd_new_entries_v2.csv
        reports/rebuild_v2/canonical_layer_summary_v2.json
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils import canonical_helix as ch  # noqa: E402

V2 = PROJECT_ROOT / "data" / "rebuild_v2"
V2R = PROJECT_ROOT / "reports" / "rebuild_v2"
SANDBOX = PROJECT_ROOT / "data" / "interim" / "rebuild_v2"
PROF_MAP = PROJECT_ROOT / "data" / "canonical_remap" / "nr_lbd_canonical_helix_map_professor_2026.csv"
SPAN = V2 / "inputs" / "analysis_span_v2.csv"
SEGMENT_SCHEMA = ["uniprot_id", "nr_code", "canonical_uniprot_position", "canonical_residue_name", "structural_segment",
                  "structural_class", "segment_order", "reference_template", "reference_pdb",
                  "structural_alignment_score", "structural_alignment_coverage", "assignment_method",
                  "assignment_confidence", "assignment_status", "evidence_notes"]
FEATURE_FILES = ["pdb_segment_coverage.parquet", "pdb_segment_flexibility_features.parquet",
                 "receptor_segment_flexibility_profiles.parquet", "receptor_segment_feature_matrix_median.parquet",
                 "receptor_segment_feature_matrix_mean.parquet", "pdb_segment_feature_matrix_median.parquet",
                 "pdb_segment_feature_matrix_mean.parquet", "receptor_segment_coverage_matrix.parquet"]
NEW_DBD_ENTRIES = ["1KB2", "1KB4", "1KB6", "1YNW"]
MIN_PDB_ENTRIES_FOR_PLOT = 5
AA3 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS",
       "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
       "Y": "TYR", "V": "VAL"}


def ref_sequence(uid: str) -> str:
    if ch.reference_accession(uid) != uid:
        return ch.read_fasta_sequence(PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table"
                                      / "uniprot_isoform_sequences" / f"{ch.reference_accession(uid)}.fasta")
    return json.loads((PROJECT_ROOT / "data" / "raw" / "api" / "uniprot" / f"{uid}.json").read_text())["sequence"]["value"]


def build_segment_map(pm: pd.DataFrame) -> pd.DataFrame:
    rows = []
    seqs = {u: ref_sequence(u) for u in pm["uniprot_id"].unique()}
    for r in pm.itertuples():
        seq = seqs[r.uniprot_id]
        base = {"uniprot_id": r.uniprot_id, "nr_code": r.nr_code, "canonical_uniprot_position": r.canonical_uniprot_position,
                "canonical_residue_name": AA3.get(seq[r.canonical_uniprot_position - 1], "UNK")
                if r.canonical_uniprot_position <= len(seq) else "",
                "reference_template": "", "structural_alignment_score": "", "structural_alignment_coverage": ""}
        if r.position_status in (ch.FPOS_ELEMENT, ch.FPOS_MULTI):
            statuses = r.assignment_status.split("|")
            for k, el in enumerate(r.canonical_elements.split("|")):
                rows.append({**base, "structural_segment": el, "structural_class": "HELIX",
                             "segment_order": ch.ELEMENT_RANK[el] + 1, "reference_pdb": r.reference_pdb,
                             "assignment_method": "PROFESSOR_CANONICAL_TABLE_FINAL",
                             "assignment_confidence": statuses[min(k, len(statuses) - 1)],
                             "assignment_status": "ASSIGNED_PROFESSOR",
                             "evidence_notes": f"{r.boundary_basis}" + (f" | {r.note}" if isinstance(r.note, str) and r.note else "")
                             + (" | position inside two professor intervals: both labels applied" if r.position_status == ch.FPOS_MULTI else "")})
        else:
            inter = r.position_status == ch.FPOS_INTER
            rows.append({**base, "structural_segment": "",
                         "structural_class": "INTER_HELIX_DERIVED" if inter else "NO_PROFESSOR_MAPPING",
                         "segment_order": "", "reference_pdb": "",
                         "assignment_method": "DERIVED_BETWEEN_PROFESSOR_INTERVALS" if inter else "NO_PROFESSOR_MAPPING",
                         "assignment_confidence": "NONE", "assignment_status": "UNASSIGNED",
                         "evidence_notes": (f"between {r.inter_helix_flank_before or 'span start'} and "
                                            f"{r.inter_helix_flank_after or 'span end'} (derived; not a loop label)") if inter
                         else "receptor has no professor canonical mapping"})
    return pd.DataFrame(rows)[SEGMENT_SCHEMA].sort_values(
        ["uniprot_id", "canonical_uniprot_position", "segment_order"], key=lambda c: c if c.name != "segment_order"
        else pd.to_numeric(c, errors="coerce")).reset_index(drop=True)


def run_feature_builder(seg_map: pd.DataFrame) -> None:
    seg_map.to_csv(SANDBOX / "data" / "analysis" / "nr_lbd_structural_segment_map.csv", index=False, lineterminator="\n")
    env = {**os.environ, "HTTP_PROXY": "http://127.0.0.1:9", "HTTPS_PROXY": "http://127.0.0.1:9", "PYTHONHASHSEED": "0"}
    log = SANDBOX / "logs" / "36_build_segment_features_professor.log"
    with open(log, "w") as fh:
        r = subprocess.run([sys.executable, "scripts/36_build_segment_features.py"], cwd=SANDBOX, env=env,
                           stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"scripts/36 failed in sandbox; see {log}")
    for f in FEATURE_FILES:
        shutil.copy2(SANDBOX / "data" / "analysis" / f, V2 / "analysis" / f)


def build_profile_join(pm: pd.DataFrame) -> pd.DataFrame:
    prof = pd.read_parquet(V2 / "analysis" / "receptor_position_flexibility_profiles.parquet")
    keep = ["uniprot_id", "canonical_uniprot_position", "position_status", "canonical_elements", "n_professor_elements",
            "assignment_status", "review_required_flag", "evidence_category", "reference_pdb", "reference_chain",
            "boundary_basis", "source_url", "note", "inter_helix_flank_before", "inter_helix_flank_after",
            "within_analysis_span", "within_prosite_lbd_v1", "reference_accession"]
    j = prof.merge(pm[keep], on=["uniprot_id", "canonical_uniprot_position"], how="left", validate="one_to_one", indicator=True)
    if (j["_merge"] != "both").any():
        bad = j[j["_merge"] != "both"][["uniprot_id", "canonical_uniprot_position"]]
        raise RuntimeError(f"profile positions without a professor-map frame row: {bad.head().to_dict('records')}")
    j = j.drop(columns="_merge")
    j["has_observation"] = j["n_pdb_entries_observed_at_position"] > 0
    j["plot_eligible_min_5_entries"] = j["n_pdb_entries_total_for_receptor"] >= MIN_PDB_ENTRIES_FOR_PLOT
    return j.sort_values(["uniprot_id", "canonical_uniprot_position"]).reset_index(drop=True)


def count_reconciliation() -> pd.DataFrame:
    def counts(root: Path, analysis: Path) -> dict:
        s = pd.read_csv(root / "processed" / "structures.csv", dtype={"pdb_id": str}, low_memory=False)
        res = pd.read_parquet(root / "processed" / "lbd_residue_bfactors_all.parquet",
                              columns=["instance_id", "is_construct_mapped", "ca_record_present", "ca_bfactor_usable"])
        norm = pd.read_parquet(analysis / "normalized_lbd_bfactors_all.parquet", columns=["row_normalization_status", "b_factor_z"])
        dm_path = root / "manifests" / "mmcif_download_manifest.csv"
        if not dm_path.exists():  # v2: input substitution (v1 manifest + the 3 approved downloads)
            dm_path = V2 / "inputs" / "mmcif_download_manifest_v2.csv"
        dm = pd.read_csv(dm_path, dtype={"pdb_id": str})
        coord = s[s["has_stage3b_coordinate_data"].astype(str) == "True"]
        prof = pd.read_parquet(analysis / "receptor_position_flexibility_profiles.parquet",
                               columns=["uniprot_id", "n_pdb_entries_total_for_receptor"])
        nent = prof.groupby("uniprot_id")["n_pdb_entries_total_for_receptor"].first()
        return {
            "receptors (project list)": len(pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "nr_metadata_validated.csv")),
            "receptors with >=1 discovery candidate": s["uniprot_id"].nunique(),
            "discovery candidates (receptor x PDB x polymer-entity master rows)": len(s),
            "unique PDB entries with a downloaded mmCIF": int((dm["download_status"] != "FAILED").sum()),
            "coordinate-assessable candidates (master rows with coordinate data)": len(coord),
            "receptors with coordinate data": coord["uniprot_id"].nunique(),
            "measured chain/model instances": res["instance_id"].nunique(),
            "residue rows (analysis-span frame x instance)": len(res),
            "construct-mapped residue rows": int(res["is_construct_mapped"].sum()),
            "observed CA residue rows": int(res["ca_record_present"].sum()),
            "B-factor-usable residue rows": int(res["ca_bfactor_usable"].sum()),
            "normalized residue rows (finite z)": int(np.isfinite(norm["b_factor_z"]).sum()),
            "receptors with >=5 PDB entries (per-receptor plots)": int((nent >= 5).sum()),
        }
    v1 = counts(PROJECT_ROOT / "data", PROJECT_ROOT / "data" / "analysis")
    v2 = counts(V2, V2 / "analysis")
    span = pd.read_csv(SPAN)
    extra = {"receptors with professor canonical mapping": (42, int(span["professor_mapping"].sum())),
             "receptors without professor canonical mapping": (6, int((~span["professor_mapping"]).sum()))}
    rows = [{"quantity": k, "unit": k.split(" (")[0], "v1_frozen": v1[k], "v2_corrected": v2[k], "change": v2[k] - v1[k]}
            for k in v1]
    rows += [{"quantity": k, "unit": "receptors", "v1_frozen": a, "v2_corrected": b, "change": b - a} for k, (a, b) in extra.items()]
    return pd.DataFrame(rows)


def normalization_edge_audit() -> tuple[pd.DataFrame, dict]:
    inst = pd.read_csv(V2 / "analysis" / "instance_bfactor_summary.csv")
    n_col = "n_ca_usable" if "n_ca_usable" in inst.columns else "instance_n_usable_b"
    s = {"n_instances": len(inst), "min_usable_residues": int(inst[n_col].min()),
         "instances_0_usable": int((inst[n_col] == 0).sum()), "instances_1_usable": int((inst[n_col] == 1).sum()),
         "instances_2_usable": int((inst[n_col] == 2).sum()),
         "normalization_status_counts": inst["normalization_status"].value_counts().to_dict()}
    if "raw_b_sd" in inst.columns:
        s["instances_zero_sd_with_ge2_usable"] = int(((inst[n_col] >= 2) & ~(inst["raw_b_sd"] > 0)).sum())
    edge = inst[(inst[n_col] <= 2) | (inst["normalization_status"] != "NORMALIZED")]
    return edge.sort_values(["uniprot_id", "instance_id"]).reset_index(drop=True), s


def dbd_report() -> pd.DataFrame:
    res = pd.read_parquet(V2 / "processed" / "lbd_residue_bfactors_all.parquet",
                          columns=["uniprot_id", "pdb_id", "instance_id", "canonical_uniprot_position", "is_construct_mapped",
                                   "ca_record_present", "ca_bfactor_usable"])
    norm = pd.read_parquet(V2 / "analysis" / "normalized_lbd_bfactors_all.parquet",
                           columns=["instance_id", "canonical_uniprot_position", "b_factor_z", "normalization_status",
                                    "row_normalization_status", "instance_n_usable_b"])
    rows = []
    for pdb in NEW_DBD_ENTRIES:
        x = res[(res["pdb_id"] == pdb) & (res["uniprot_id"] == "P11473")]
        for inst, g in x.groupby("instance_id"):
            n = norm[norm["instance_id"] == inst]
            m = g[g["is_construct_mapped"]]
            rows.append({"pdb_id": pdb, "instance_id": inst,
                         "mapped_residues_in_final_vdr_span": len(m),
                         "mapped_positions": " ".join(map(str, sorted(m["canonical_uniprot_position"]))),
                         "observed_ca_residues": int(g["ca_record_present"].sum()),
                         "raw_bfactor_usable_rows": int(g["ca_bfactor_usable"].sum()),
                         "normalized_rows_finite_z": int(np.isfinite(n["b_factor_z"]).sum()),
                         "instance_normalization_status": n["normalization_status"].iloc[0] if len(n) else "",
                         "row_normalization_status": "|".join(sorted(n["row_normalization_status"].unique())) if len(n) else "",
                         "new_in_v2": pdb != "1YNW" or True})
    if not rows:
        rows.append({"pdb_id": "|".join(NEW_DBD_ENTRIES), "instance_id": "", "mapped_residues_in_final_vdr_span": 0})
    return pd.DataFrame(rows)


def build(run_features: bool = True) -> dict:
    pm = pd.read_csv(PROF_MAP, keep_default_na=False, dtype={"canonical_uniprot_position": int})
    seg = build_segment_map(pm)
    seg.to_csv(V2 / "analysis" / "nr_lbd_structural_segment_map.csv", index=False, lineterminator="\n")
    if run_features:
        run_feature_builder(seg)
    pj = build_profile_join(pm)
    pj.to_csv(V2 / "analysis" / "receptor_position_professor_profile_v2.csv", index=False, lineterminator="\n")
    V2R.mkdir(parents=True, exist_ok=True)
    cr = count_reconciliation()
    cr.to_csv(V2R / "count_reconciliation_v2.csv", index=False, lineterminator="\n")
    edge, es = normalization_edge_audit()
    edge.to_csv(V2R / "normalization_edge_case_audit_v2.csv", index=False, lineterminator="\n")
    dbd = dbd_report()
    dbd.to_csv(V2R / "vdr_dbd_new_entries_v2.csv", index=False, lineterminator="\n")
    summ = {"policy": ch.FINAL_POLICY_VERSION, "segment_map_rows": len(seg),
            "segment_map_labels": sorted(set(seg["structural_segment"]) - {""}),
            "profile_join_rows": len(pj), "normalization_edges": es}
    (V2R / "canonical_layer_summary_v2.json").write_text(json.dumps(summ, indent=2, default=int) + "\n")
    return summ


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--no-features", action="store_true")
    print(json.dumps(build(not ap.parse_args().no_features), indent=2, default=int))


if __name__ == "__main__":
    main()
