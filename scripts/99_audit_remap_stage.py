"""
Canonical-remap stage audits (read-only; nothing frozen is modified):

  1. VDR per-entry numbering audit (all VDR coordinate entries).
  2. Downstream-impact assessment of the residue-numbering corrections:
     by receptor, and by every legacy tracked table.
  3. Entity-level consistency of the corrected numbering (all instances of
     one polymer entity must receive the same corrected map).
  4. This stage's OWN receptor-identity audit, reusing the unmodified
     checking functions of scripts/93_final_receptor_identity_audit.py on the
     stage directories only -- the historical frozen report
     (reports/tables/final_receptor_identity_audit.csv) is never rewritten.

Reads:
    data/canonical_remap/residue_numbering_crosswalk_v1.parquet   (scripts/96)
    reports/canonical_remap/residue_numbering_instance_audit_v1.csv (scripts/96)
    data/processed/lbd_residue_bfactors_all.parquet               (frozen)
    data/manifests/nr_lbd_reference.csv
    all tracked legacy tables under data/analysis, data/processed, data/manifests, reports/tables

Writes:
    reports/canonical_remap/vdr_entry_numbering_audit_v1.csv
    reports/canonical_remap/numbering_impact_by_receptor_v1.csv
    reports/canonical_remap/numbering_impact_by_table_v1.csv
    reports/canonical_remap/numbering_entity_consistency_v1.csv
    reports/canonical_remap/vdr_downstream_impact_v1.csv
    reports/canonical_remap/stage_receptor_identity_audit_v1.csv
    reports/canonical_remap/stage_receptor_identity_findings_v1.csv

Usage:
    python scripts/99_audit_remap_stage.py [--output-root DIR]
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

CROSSWALK = PROJECT_ROOT / "data" / "canonical_remap" / "residue_numbering_crosswalk_v1.parquet"
INSTANCE_AUDIT = PROJECT_ROOT / "reports" / "canonical_remap" / "residue_numbering_instance_audit_v1.csv"
FROZEN_RESIDUES = PROJECT_ROOT / "data" / "processed" / "lbd_residue_bfactors_all.parquet"
LBD_REFERENCE = PROJECT_ROOT / "data" / "manifests" / "nr_lbd_reference.csv"

OUT = Path("reports") / "canonical_remap"
STAGE_DIRS = ["data/canonical_remap", "reports/canonical_remap", "data/sources/professor_canonical_helix_table",
              "data/rebuild_v2", "reports/rebuild_v2"]
LEGACY_DIRS = ["data/analysis", "data/processed", "data/manifests", "reports/tables"]

# Legacy files produced BEFORE residue numbering (scripts 01-18, metadata / discovery /
# inventory / download / ligand annotation), i.e. independent of the defect.
UPSTREAM_INDEPENDENT_PREFIXES = (
    "nr_metadata", "uniprot_api_manifest", "nr_lbd_reference", "nr_lbd_boundary_audit", "nr_lbd_analysis_policy",
    "rcsb_search_manifest", "rcsb_discovery", "rcsb_entries", "rcsb_data_api_manifest", "rcsb_polymer_entities",
    "rcsb_polymer_entity_uniprot_mappings", "rcsb_polymer_instances", "rcsb_instance_api_manifest",
    "rcsb_instance_mapping_anomalies", "rcsb_candidate_inventory", "rcsb_candidate_instance_inventory",
    "rcsb_multi_uniprot_entity_audit", "rcsb_source_taxonomy_audit",
    "rcsb_sequence_coordinates", "rcsb_short_mapping_audit", "rcsb_stage2_provenance", "stage2_release_summary",
    "mmcif_download_manifest", "nr_functional_features",
)
# NOTE: final_nonpolymer_inventory / final_ligand_annotations (26) and final_polymer_partner_annotations (27)
# are DOWNSTREAM: both decide "LBD-contacting" from the numbered LBD-window residues and record contact
# positions in UniProt numbering. rcsb_nonpolymer_data_api_manifest (25) depends on the selection (24).
# NOTE: rcsb_coverage_validation (09), rcsb_candidate_lbd_mapping (12), rcsb_precoordinate_qc and
# stage3a_* (14), rcsb_lbd_boundary_sensitivity (16) compute LBD coverage from the same RCSB regions
# and are therefore treated as DOWNSTREAM (conservative).
# Downstream tables whose CONTENT does not depend on residue positions:
NUMBERING_INDEPENDENT_BY_CONTENT = {
    "final_receptor_identity_audit.csv": "receptor-identity text scan of other files",
    "subfamily_receptor_metadata.csv": "receptor group labels",
}
SOURCE_OF_DEFECT = {"rcsb_sequence_coordinate_alignments.csv"}
CHANGED = ("CORRECTED", "WITHHELD_PREVIOUSLY_MAPPED")


def write(df: pd.DataFrame, root: Path, name: str) -> None:
    p = root / OUT / name
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False, lineterminator="\n")


def rng(s: pd.Series) -> str:
    s = s.dropna()
    return f"{int(s.min())}-{int(s.max())}" if len(s) else ""


# ---------------------------------------------------------------------------


def vdr_audit(cw: pd.DataFrame, ia: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for t in ia[ia["uniprot_id"] == "P11473"].sort_values("pdb_id").itertuples():
        x = cw[cw["instance_id"] == t.instance_id]
        new = x["corrected_uniprot_position"]
        present_165_215 = int(((new >= 165) & (new <= 215)).sum())
        old_win = x[(x["old_canonical_uniprot_position"] >= 150) & (x["old_canonical_uniprot_position"] <= 240)]
        segs = []
        for off, g in old_win.groupby((old_win["corrected_uniprot_position"] - old_win["old_canonical_uniprot_position"])
                                      .fillna(-9999).ne((old_win["corrected_uniprot_position"]
                                                         - old_win["old_canonical_uniprot_position"]).fillna(-9999).shift())
                                      .cumsum()):
            o, n = g["old_canonical_uniprot_position"], g["corrected_uniprot_position"]
            segs.append(f"old {rng(o)} -> " + (f"new {rng(n)}" if n.notna().all() else "withheld"))
        before = x[x["correction_status"] == "UNCHANGED"]
        after = x[x["correction_status"] == "CORRECTED"]
        obs_after = after[after["frozen_ca_record_present"].fillna(False)]
        rows.append({
            "pdb_id": t.pdb_id, "instance_id": t.instance_id, "chain_polymer_length": t.chain_polymer_length,
            "uniprot_165_215_in_chain": ("PRESENT" if present_165_215 == 51 else
                                         "DELETED" if present_165_215 == 0 else f"PARTIAL({present_165_215}/51)"),
            "aligned_uniprot_span": f"{t.aligned_uniprot_start}-{t.aligned_uniprot_end}",
            "alignment_gaps": t.alignment_gaps, "depositor_segments": t.depositor_segments,
            "old_mapping_150_240": "; ".join(segs),
            "offsets_new_minus_old": t.offsets_new_minus_old, "offset_segments_label_seq": t.offset_segments_label_seq,
            "breakpoint_old_position": (int(after["old_canonical_uniprot_position"].min()) if len(after) else None),
            "n_residues_unchanged": t.n_unchanged, "n_residues_corrected": t.n_corrected,
            "n_residues_withheld": t.n_withheld_previously_mapped,
            "old_mapping_before_breakpoint_correct": bool(len(before)) and bool(
                (before["residue_relation"].isin(["IDENTICAL", "SUBSTITUTION"])).all()),
            "old_mapping_after_breakpoint_identity": (round(float((obs_after["chain_residue"] == obs_after["uniprot_residue_at_old"]).mean()), 4)
                                                      if len(obs_after) else None),
            "corrected_mapping_after_breakpoint_identity": (round(float((obs_after["chain_residue"] == obs_after["uniprot_residue_at_corrected"]).mean()), 4)
                                                            if len(obs_after) else None),
            "numbering_class": t.numbering_class,
            "status": ("OLD_MAPPING_CORRECT" if t.numbering_class == "CLEAN" else
                       "CORRECTED_HIGH_CONFIDENCE" if t.numbering_class == "DELETION_INSERTION_OFFSET" else t.numbering_class),
        })
    return pd.DataFrame(rows)


def impact_by_receptor(cw: pd.DataFrame, ia: pd.DataFrame) -> pd.DataFrame:
    frozen = pd.read_parquet(FROZEN_RESIDUES, columns=["uniprot_id", "instance_id", "entity_label_seq_id",
                                                       "ca_bfactor_usable", "is_construct_mapped"])
    frozen = frozen[frozen["is_construct_mapped"]].rename(columns={"entity_label_seq_id": "label_seq_id"})
    frozen["label_seq_id"] = frozen["label_seq_id"].astype("Int64")
    lbd = pd.read_csv(LBD_REFERENCE).set_index("uniprot_id")[["lbd_start", "lbd_end"]]
    x = cw.merge(frozen[["instance_id", "label_seq_id", "ca_bfactor_usable"]], on=["instance_id", "label_seq_id"], how="left")
    x = x.join(lbd, on="uniprot_id")
    in_old = x["old_canonical_uniprot_position"].between(x["lbd_start"], x["lbd_end"]).fillna(False)
    in_new = x["corrected_uniprot_position"].between(x["lbd_start"], x["lbd_end"]).fillna(False)
    x["leaves_window"] = in_old & ~in_new
    x["enters_window"] = ~in_old & in_new
    x["window_membership_changes"] = x["leaves_window"] | x["enters_window"]
    rows = []
    for uid, g in x.groupby("uniprot_id", sort=True):
        ii = ia[ia["uniprot_id"] == uid]
        obs = g[g["in_frozen_lbd_table"] & g["frozen_ca_record_present"].fillna(False)]
        chg = g[g["correction_status"].isin(CHANGED)]
        chg_obs = obs[obs["correction_status"].isin(CHANGED)]
        usable_window_change = g[g["window_membership_changes"] & (g["ca_bfactor_usable"].fillna(False) | g["enters_window"])]
        z_inst = usable_window_change["instance_id"].unique()
        changed_inst = chg["instance_id"].unique()
        rows.append({
            "uniprot_id": uid, "nr_code": g["nr_code"].iloc[0], "common_name": g["common_name"].iloc[0],
            "n_pdb_entries": ii["pdb_id"].nunique(), "n_instances": len(ii),
            "n_pdb_entries_changed": chg["pdb_id"].nunique(),
            "fraction_pdb_entries_changed": round(chg["pdb_id"].nunique() / ii["pdb_id"].nunique(), 4),
            "n_instances_with_shifted_residues": int((ii["n_corrected"] > 0).sum()),
            "n_instances_with_withheld_residues_only": int(((ii["n_corrected"] == 0) & (ii["n_withheld_previously_mapped"] > 0)).sum()),
            "n_observed_residue_rows": len(obs), "n_observed_residue_rows_changed": len(chg_obs),
            "fraction_observed_residue_rows_changed": round(len(chg_obs) / len(obs), 6) if len(obs) else 0.0,
            "old_positions_affected": rng(chg["old_canonical_uniprot_position"]),
            "corrected_positions_of_shifted_residues": rng(chg["corrected_uniprot_position"]),
            "n_instances_whose_lbd_residue_set_changes": len(z_inst),
            "n_frozen_rows_in_instances_with_any_change": int(cw[cw["instance_id"].isin(changed_inst) & cw["in_frozen_lbd_table"]].shape[0]),
            "severity": ("SHIFTED_RESIDUES" if (ii["n_corrected"] > 0).any() else
                         "BOUNDARY_RESIDUES_WITHHELD_ONLY" if len(chg) else "UNAFFECTED"),
            "pdb_entries_changed": " ".join(sorted(chg["pdb_id"].unique())),
        })
    return pd.DataFrame(rows)


def legacy_tables() -> list[Path]:
    tracked = subprocess.run(["git", "ls-files", *LEGACY_DIRS], cwd=PROJECT_ROOT, capture_output=True, text=True,
                             check=True).stdout.split()
    return [PROJECT_ROOT / f for f in sorted(tracked) if f.endswith((".csv", ".parquet"))]


def impact_by_table(imp: pd.DataFrame, cw: pd.DataFrame) -> pd.DataFrame:
    shifted_rec = set(imp[imp["severity"] == "SHIFTED_RESIDUES"]["uniprot_id"])
    any_rec = set(imp[imp["severity"] != "UNAFFECTED"]["uniprot_id"])
    chg = cw[cw["correction_status"].isin(CHANGED)]
    chg_pairs = set(zip(chg["uniprot_id"], chg["pdb_id"]))
    chg_inst = set(chg["instance_id"])
    name_to_uid = dict(zip(imp["common_name"], imp["uniprot_id"])) | dict(zip(imp["nr_code"], imp["uniprot_id"]))
    rows = []
    for p in legacy_tables():
        rel = str(p.relative_to(PROJECT_ROOT))
        if p.name in SOURCE_OF_DEFECT:
            dep = "SOURCE_OF_DEFECT"
        elif p.name in NUMBERING_INDEPENDENT_BY_CONTENT:
            dep = "NUMBERING_INDEPENDENT_BY_CONTENT"
        elif p.name.startswith(UPSTREAM_INDEPENDENT_PREFIXES):
            dep = "UPSTREAM_OF_RESIDUE_NUMBERING"
        else:
            dep = "DOWNSTREAM_OF_RESIDUE_NUMBERING"
        cols = (pq.read_schema(p).names if p.suffix == ".parquet"
                else pd.read_csv(p, nrows=0).columns.tolist())
        keys = [c for c in ("instance_id", "pdb_id", "uniprot_id", "nr_code", "common_name") if c in cols]
        df = (pd.read_parquet(p, columns=keys) if p.suffix == ".parquet" and keys else
              pd.read_csv(p, usecols=keys, dtype=str, keep_default_na=False) if keys else None)
        if df is not None and "uniprot_id" not in df.columns:
            src = "nr_code" if "nr_code" in df.columns else ("common_name" if "common_name" in df.columns else None)
            if src:
                df["uniprot_id"] = df[src].map(name_to_uid)
        n_rows = len(df) if df is not None else None
        level, n_aff = "NO_RECEPTOR_KEY", None
        if df is not None and "instance_id" in df.columns:
            level, n_aff = "INSTANCE", int(df["instance_id"].isin(chg_inst).sum())
        elif df is not None and "pdb_id" in df.columns and "uniprot_id" in df.columns:
            level = "PDB_ENTRY"
            n_aff = int(pd.Series(list(zip(df["uniprot_id"], df["pdb_id"]))).isin(chg_pairs).sum())
        elif df is not None and "uniprot_id" in df.columns:
            level, n_aff = "RECEPTOR", int(df["uniprot_id"].isin(any_rec).sum())
        has_pos = "canonical_uniprot_position" in cols
        if dep != "DOWNSTREAM_OF_RESIDUE_NUMBERING":
            verdict = {"SOURCE_OF_DEFECT": "SOURCE_OF_DEFECT", "NUMBERING_INDEPENDENT_BY_CONTENT": "UNAFFECTED_BY_CONTENT"}.get(
                dep, "UNAFFECTED_UPSTREAM")
        elif level == "NO_RECEPTOR_KEY":
            verdict = "POTENTIALLY_AFFECTED_CROSS_RECEPTOR_AGGREGATE"
        elif n_aff and level == "RECEPTOR":
            verdict = "RECEPTOR_LEVEL_ROWS_FROM_CHANGED_RECEPTORS"
        elif n_aff:
            verdict = "ENTRY_OR_INSTANCE_ROWS_WITH_CHANGED_NUMBERING"
        else:
            verdict = "NO_ROWS_FROM_CHANGED_ENTRIES"
        n_shift = None
        if df is not None and "uniprot_id" in df.columns:
            n_shift = int(df["uniprot_id"].isin(shifted_rec).sum())
        rows.append({"file": rel, "numbering_dependency": dep, "row_key_level": level, "n_rows": n_rows,
                     "n_rows_from_changed_entries_or_receptors": n_aff, "n_rows_from_receptors_with_shifted_residues": n_shift,
                     "has_canonical_uniprot_position": has_pos, "verdict": verdict})
    return pd.DataFrame(rows)


def vdr_downstream(cw: pd.DataFrame, ia: pd.DataFrame) -> pd.DataFrame:
    """Quantified VDR-specific inheritance of the numbering defect (no biology)."""
    shifted = set(ia[(ia["uniprot_id"] == "P11473") & (ia["n_corrected"] > 0)]["pdb_id"])
    prof = pd.read_parquet(PROJECT_ROOT / "data" / "analysis" / "pdb_position_bfactor_profiles.parquet",
                           columns=["uniprot_id", "pdb_id", "canonical_uniprot_position"],
                           filters=[("uniprot_id", "==", "P11473")])
    def frac(lo, hi):
        x = prof[prof["canonical_uniprot_position"].between(lo, hi)]
        return len(x), int(x["pdb_id"].isin(shifted).sum())
    n220, s220 = frac(220, 320)
    n166, s166 = frac(166, 423)
    n127, s127 = frac(127, 164)
    src = pd.read_csv(PROJECT_ROOT / "data" / "analysis" / "receptor_segment_mapping_sources.csv", keep_default_na=False)
    via_vdr = src[src["source_method"].str.contains("P11473") | src["source_pdb"].str.contains("P11473")]
    refs = pd.read_csv(PROJECT_ROOT / "data" / "analysis" / "canonical_lbd_reference_structures_v2.csv", keep_default_na=False)
    vref = refs[refs.iloc[:, 0] == "P11473"]
    rows = [
        ("old per-position profile, positions 220-320",
         f"{s220}/{n220} entry-level observations come from the {len(shifted)} shifted entries",
         "in those entries old position p holds UniProt residue p+51 (i.e. 271-371)"),
        ("old per-position profile, positions 166-423", f"{s166}/{n166} entry-level observations from shifted entries", ""),
        ("old per-position profile, positions 127-164", f"{s127}/{n127} from shifted entries",
         "numbering before the deletion is correct in every VDR entry"),
        ("old v3 VDR segment map", "derived from receptor_dssp_consensus over the 48 VDR entries (45 shifted) and VDR templates "
         + "; ".join(f"{r.iloc[2]} ({r.iloc[4]}): {'shifted' if r.iloc[2] in shifted else 'correct'}" for _, r in vref.iterrows()),
         "segment positions >165 in the v3 VDR map inherit the shift"),
        ("old v3 labels transferred from the VDR anchor template",
         f"{len(via_vdr)} receptor x segment labels in {via_vdr['uniprot_id'].nunique()} receptors list P11473 in their transfer chain",
         "transfer used VDR template coordinates whose residue numbering is shifted in the PRIMARY template; "
         "v3 is superseded by the professor map, so these labels are not carried forward"),
        ("canonical helix join (plot input)", "VDR receptor_numbering_integrity = SUBSTANTIAL_FRACTION_OF_ENTRIES_SHIFTED_DO_NOT_INTERPRET",
         "the professor intervals are correct UniProt numbering; the frozen VDR profile under them is not"),
        ("flagged entries 5V39 4G2I 3W0Y 2HAR",
         " ".join(f"{p}:{'+51 shifted' if p in shifted else 'correct'}" for p in ["5V39", "4G2I", "3W0Y", "2HAR"]),
         f"same correction as the other shifted entries; correctly numbered VDR entries: "
         + " ".join(sorted(set(ia[ia['uniprot_id'] == 'P11473']['pdb_id']) - shifted))),
    ]
    return pd.DataFrame(rows, columns=["item", "quantity", "note"])


def entity_consistency(cw: pd.DataFrame) -> pd.DataFrame:
    m = cw[cw["corrected_uniprot_position"].notna()]
    g = m.groupby(["uniprot_id", "polymer_entity_id", "label_seq_id"])["corrected_uniprot_position"].nunique()
    bad = g[g > 1].reset_index()
    per = cw.groupby(["uniprot_id", "polymer_entity_id"]).agg(n_instances=("instance_id", "nunique")).reset_index()
    per["n_label_seq_with_conflicting_corrected_positions"] = per.apply(
        lambda r: int(((bad["uniprot_id"] == r["uniprot_id"]) & (bad["polymer_entity_id"] == r["polymer_entity_id"])).sum()), axis=1)
    return per[per["n_instances"] > 1].reset_index(drop=True)


def stage_identity_audit() -> tuple[pd.DataFrame, pd.DataFrame]:
    spec = importlib.util.spec_from_file_location("a93", PROJECT_ROOT / "scripts" / "93_final_receptor_identity_audit.py")
    a93 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(a93)
    files = sorted(f for d in STAGE_DIRS for ext in ("*.csv", "*.parquet") for f in (PROJECT_ROOT / d).rglob(ext))
    issues = [i for f in files for i in a93.audit_structured_file(f)]
    stage_md = {str(f.relative_to(PROJECT_ROOT)) for d in STAGE_DIRS for f in (PROJECT_ROOT / d).rglob("*.md")}
    issues += [i for i in a93.audit_prose_files() if i.get("file") in stage_md]
    rows = []
    for uid, auth in sorted(a93.AUTHORITATIVE.items()):
        gen = sum(1 for i in issues if i.get("severity") == "GENUINE"
                  and (i.get("uniprot_id") == uid or f"uniprot {uid}" in i.get("issue", "")))
        ben = sum(1 for i in issues if i.get("severity") == "BENIGN_ALIAS" and i.get("uniprot_id") == uid)
        rows.append({"uniprot_id": uid, "authoritative_nr_code": auth["nr_code"],
                     "authoritative_common_name": auth["common_name"], "stage_files_checked": len(files),
                     "stage_markdown_files_checked": len(stage_md), "n_inconsistencies": gen,
                     "n_documented_benign_aliases": ben,
                     "status": "GENUINE_INCONSISTENCY_UNRESOLVED" if gen else (
                         "CLEAN_WITH_DOCUMENTED_BENIGN_ALIASES" if ben else "CLEAN")})
    findings = pd.DataFrame(issues, columns=["file", "line", "uniprot_id", "issue", "line_text", "severity"])
    return pd.DataFrame(rows), findings.sort_values(["file", "line"]).reset_index(drop=True)


def build(output_root: Path) -> dict:
    cw = pd.read_parquet(CROSSWALK)
    ia = pd.read_csv(INSTANCE_AUDIT, keep_default_na=False)
    vdr = vdr_audit(cw, ia)
    imp = impact_by_receptor(cw, ia)
    tab = impact_by_table(imp, cw)
    ent = entity_consistency(cw)
    vdr_ds = vdr_downstream(cw, ia)
    ident, findings = stage_identity_audit()
    for df, name in [(vdr, "vdr_entry_numbering_audit_v1.csv"), (imp, "numbering_impact_by_receptor_v1.csv"),
                     (tab, "numbering_impact_by_table_v1.csv"), (ent, "numbering_entity_consistency_v1.csv"), (vdr_ds, "vdr_downstream_impact_v1.csv"),
                     (ident, "stage_receptor_identity_audit_v1.csv"), (findings, "stage_receptor_identity_findings_v1.csv")]:
        write(df, output_root, name)
    return {"vdr_status": vdr["status"].value_counts().to_dict(),
            "receptor_severity": imp["severity"].value_counts().to_dict(),
            "table_verdicts": tab["verdict"].value_counts().to_dict(),
            "entity_conflicts": int(ent["n_label_seq_with_conflicting_corrected_positions"].sum()),
            "stage_identity_status": ident["status"].value_counts().to_dict()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    print(build(ap.parse_args().output_root))


if __name__ == "__main__":
    main()
