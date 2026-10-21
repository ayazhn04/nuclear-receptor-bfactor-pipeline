"""
Versioned foundation rebuild driver.

Runs the UNMODIFIED legacy pipeline scripts (09 ... 33) inside an isolated
sandbox copy of the project, so no frozen file in the working tree can be
touched, then exports the results to a versioned location.

Sandbox (gitignored): data/interim/rebuild_<variant>/
    scripts/, config/                copied from the working tree (legacy code unchanged)
    data/manifests, data/processed, data/analysis, reports/tables   copied tracked v1 files
    data/raw/mmcif                   symlink to the local mmCIF cache (read-only use)
    data/raw/api/uniprot             copy of the cached UniProt entries
    data/raw/api/rcsb/data           RCSB API caches RECONSTRUCTED locally (see below)

Variants:
    v1_repro  frozen inputs; must reproduce the frozen tracked outputs (proves the
              sandbox + metadata reconstruction are faithful)
    v2        corrected inputs from data/rebuild_v2/inputs/ (built by
              scripts/102_build_v2_inputs.py): corrected entity->UniProt regions,
              final analysis spans, TRalpha reference sequence P10827-2

RCSB API caches (not migrated to this machine) are reconstructed WITHOUT network:
    nonpolymer entity metadata   <- frozen data/manifests/final_nonpolymer_inventory.csv
    entry polymer entity ids     <- local mmCIF _entity_poly
    other polymer entity metadata <- frozen final_polymer_partner_annotations.csv where present,
                                     otherwise local mmCIF (_entity, _entity_poly_seq, _struct_ref)
All network access is blocked (proxy env to a closed local port).

Script 18 is NOT run as a whole: only its build_pool() is executed, then the
selected PDB set is checked against the local cache (checkpoint A). Script 25
(network fetch) is not run; checkpoint B verifies it is not needed.

Usage:
    python scripts/101_rebuild_foundation_sandbox.py --variant v1_repro|v2 [--from-step N] [--export]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.mmcif_parse import as_list, parse_mmcif_gz  # noqa: E402

PY = sys.executable
V2_INPUTS = PROJECT_ROOT / "data" / "rebuild_v2" / "inputs"
EXPORT_DATA = PROJECT_ROOT / "data" / "rebuild_v2"
EXPORT_REPORTS = PROJECT_ROOT / "reports" / "rebuild_v2"

# mmCIF caches: the frozen release cache (exactly the 1777 v1 files) plus the separate,
# gitignored cache of the files downloaded with explicit approval for v2 (1KB2/1KB4/1KB6).
MMCIF_CACHES = [PROJECT_ROOT / "data" / "raw" / "mmcif", PROJECT_ROOT / "data" / "interim" / "raw_mmcif_v2_additions"]


def mmcif_path(pdb: str) -> Path:
    for c in MMCIF_CACHES:
        p = c / f"{pdb}.cif.gz"
        if p.exists():
            return p
    return MMCIF_CACHES[0] / f"{pdb}.cif.gz"


COPY_DIRS = ["scripts", "config", "data/manifests", "data/processed", "data/analysis", "reports/tables"]
# Ordered pipeline (dependency order derived from the scripts' actual reads/writes).
# Only when the analysis span / reference definitions change (v2): boundary audit and
# functional features read nr_lbd_reference.csv.
V2_ONLY_STEPS = {"15", "17"}
STEPS: list[tuple[str, str]] = [
    ("15", "15_build_lbd_boundary_audit.py"),
    ("17", "17_build_functional_features.py"),
    ("09", "09_finalize_stage2_coverage.py"),
    ("12", "12_map_candidates_to_lbd.py"),
    ("14", "14_build_precoordinate_qc.py"),
    ("16", "16_lbd_boundary_sensitivity.py"),
    ("18", "CHECKPOINT_A_POOL"),
    ("19", "19_process_lbd_coordinates.py"),
    ("20", "20_build_stage3b_summaries.py"),
    ("21", "21_apply_occupancy_correction.py"),
    ("22", "22_recompute_stage3b_metrics.py"),
    ("23", "23_stage3b2_threshold_audit.py"),
    ("24", "24_final_structure_selection.py"),
    ("25", "CHECKPOINT_B_SELECTION"),
    ("26", "26_nonpolymer_inventory_and_ligand_state.py"),
    ("27", "27_polymer_partner_annotation.py"),
    ("28", "28_build_final_residue_map.py"),
    ("29", "29_build_structures_csv.py"),
    ("30", "30_finalize_release.py"),
    ("31", "31_build_master_inventory.py"),
    ("32", "32_build_normalization_foundation.py"),
    ("33", "33_build_analysis_foundation_reports.py"),
]
# Sandbox-only adaptations of legacy scripts for the v2 variant. Each replaces a
# hard-coded FROZEN-RELEASE count (a regression guard, not a scientific rule) by the
# equivalent data-derived identity, verified to hold exactly on the frozen v1 data.
# The working-tree scripts are never edited; every applied patch is logged.
SANDBOX_PATCHES_V2: list[tuple[str, str, str, str]] = [
    ("32_build_normalization_foundation.py",
     'assert len(out) == len(res) == 771542, "row count changed"',
     'assert len(out) == len(res), "row count changed"',
     "normalization must preserve the residue-row count (v1 literal 771542 replaced by the identity itself)"),
    ("32_build_normalization_foundation.py",
     'assert len(master) == 2072 and int(master["has_stage3b_coordinate_data"].sum()) == 1840',
     'assert len(master) == 2072 and int(master["has_stage3b_coordinate_data"].sum()) == int('
     'pd.read_csv(MANIFESTS / "stage3b_coordinate_pool.csv")["coordinate_audit_selected"].sum())',
     "coordinate candidates must equal the selected coordinate pool (v1: 1840 == 1840)"),
    ("32_build_normalization_foundation.py",
     "assert len(res) == 771542",
     'assert len(res) == int(res.groupby("uniprot_id")["instance_id"].nunique().mul('
     'pd.read_csv(LBD_REFERENCE_CSV).set_index("uniprot_id")["lbd_length"]).dropna().sum())',
     "residue rows must equal sum(instances x analysis-span length) per receptor (v1: 771542 == 771542)"),
    ("34_build_dssp_consensus.py",
     "assert len(coord) == 1840",
     'assert len(coord) == int(pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "stage3b_coordinate_pool.csv")'
     '["coordinate_audit_selected"].sum())',
     "coordinate candidates must equal the selected coordinate pool (v1: 1840 == 1840; v2: 1844 == 1844)"),
    ("34_build_dssp_consensus.py",
     'assert rep["instance_id"].nunique() == len(rep_ids) == 1840',
     'assert rep["instance_id"].nunique() == len(rep_ids) == len(coord)',
     "one representative instance per coordinate candidate, all present in the normalized table (v1: 1840; v2: 1844)"),
]


def apply_sandbox_patches(w: Path) -> list[dict]:
    log = []
    for script, old, new, reason in SANDBOX_PATCHES_V2:
        f = w / "scripts" / script
        text = f.read_text()
        if new in text:
            status = "ALREADY_APPLIED"
        elif text.count(old) == 1:
            f.write_text(text.replace(old, new))
            status = "APPLIED"
        else:
            raise SystemExit(f"patch target not found exactly once in {script}: {old!r}")
        log.append({"script": script, "replaced": old, "with": new, "reason": reason, "status": status})
    (w / "logs" / "sandbox_patches_v2.json").write_text(json.dumps(log, indent=2))
    return log


NO_NETWORK_ENV = {"HTTP_PROXY": "http://127.0.0.1:9", "HTTPS_PROXY": "http://127.0.0.1:9",
                  "ALL_PROXY": "http://127.0.0.1:9", "http_proxy": "http://127.0.0.1:9",
                  "https_proxy": "http://127.0.0.1:9", "NO_PROXY": "", "no_proxy": ""}


def work_dir(variant: str) -> Path:
    # REBUILD_WORK_SUFFIX lets an independent second run (determinism check) use its own sandbox
    return PROJECT_ROOT / "data" / "interim" / f"rebuild_{variant}{os.environ.get('REBUILD_WORK_SUFFIX', '')}"


# ---------------------------------------------------------------------------
# sandbox construction
# ---------------------------------------------------------------------------


def tracked_files(rel_dir: str) -> list[str]:
    out = subprocess.run(["git", "ls-files", rel_dir], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True)
    return out.stdout.split()


def build_sandbox(variant: str) -> Path:
    w = work_dir(variant)
    if w.exists():
        shutil.rmtree(w)
    for d in COPY_DIRS:
        for f in tracked_files(d):
            dst = w / f
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(PROJECT_ROOT / f, dst)
    for f in ["README.md", "TEAM_HANDOFF.md", "environment.yml", "environment-lock.yml", ".gitattributes"]:
        if (PROJECT_ROOT / f).exists():
            shutil.copy2(PROJECT_ROOT / f, w / f)  # hashed as provenance by scripts/30
    # utils that are untracked but imported by nothing in the legacy chain are copied anyway
    for f in (PROJECT_ROOT / "scripts" / "utils").glob("*.py"):
        shutil.copy2(f, w / "scripts" / "utils" / f.name)
    (w / "data" / "raw").mkdir(parents=True, exist_ok=True)
    union = w / "data" / "raw" / "mmcif"  # union of both caches, read-only symlinks
    union.mkdir(parents=True)
    for cache in MMCIF_CACHES:
        if cache.is_dir():
            for f in sorted(cache.glob("*.cif.gz")):
                if not (union / f.name).exists():
                    (union / f.name).symlink_to(f)
    shutil.copytree(PROJECT_ROOT / "data" / "raw" / "api" / "uniprot", w / "data" / "raw" / "api" / "uniprot")
    for d in ["data/interim/qc", "data/interim/mappings", "data/interim/ligand_annotations", "logs"]:
        (w / d).mkdir(parents=True, exist_ok=True)
    reconstruct_rcsb_caches(w)
    if variant == "v2":
        apply_v2_inputs(w)
    return w


def apply_v2_inputs(w: Path) -> None:
    manifest = json.loads((V2_INPUTS / "v2_inputs_manifest.json").read_text())
    for rel_dst, rel_src in manifest["substitutions"].items():
        shutil.copy2(PROJECT_ROOT / rel_src, w / rel_dst)


def reconstruct_rcsb_caches(w: Path) -> None:
    out = w / "data" / "raw" / "api" / "rcsb" / "data"
    out.mkdir(parents=True, exist_ok=True)
    man = PROJECT_ROOT / "data" / "manifests"

    # --- nonpolymer entities (script 26) from the frozen inventory ---
    inv = pd.read_csv(man / "final_nonpolymer_inventory.csv", dtype={"pdb_id": str}, keep_default_na=False)
    inv = inv[inv["nonpolymer_entity_id"] != ""]
    ents = []
    for eid, g in inv.groupby("nonpolymer_entity_id", sort=True):
        f = g.iloc[0]
        for c in ["chemical_name", "formula", "formula_weight", "chem_comp_type", "atom_count_heavy",
                  "rcsb_subject_of_investigation"]:
            if g[c].nunique() > 1:
                raise RuntimeError(f"{eid}: inconsistent {c} across frozen inventory rows")
        entry, ent = eid.split("_", 1)
        ents.append({
            "rcsb_id": eid,
            "rcsb_nonpolymer_entity_container_identifiers": {"entry_id": entry, "entity_id": ent,
                                                             "nonpolymer_comp_id": f["component_id"]},
            "nonpolymer_comp": {"chem_comp": {"name": f["chemical_name"] or None, "formula": f["formula"] or None,
                                              "formula_weight": float(f["formula_weight"]) if f["formula_weight"] != "" else None,
                                              "type": f["chem_comp_type"] or None},
                                "rcsb_chem_comp_info": {"atom_count_heavy": int(float(f["atom_count_heavy"]))
                                                        if f["atom_count_heavy"] != "" else None}},
            "rcsb_nonpolymer_entity_annotation": ([{"type": "SUBJECT_OF_INVESTIGATION"}]
                                                  if str(f["rcsb_subject_of_investigation"]) == "True" else []),
            "nonpolymer_entity_instances": [],
        })
    frozen_np_ids = {e["rcsb_id"] for e in ents}
    np_prov = [{"nonpolymer_entity_id": e["rcsb_id"], "source": "frozen_nonpolymer_inventory"} for e in ents]

    # --- polymer entities (script 27): entry ids from mmCIF, partner metadata from frozen table or mmCIF ---
    audit = pd.read_csv(man / "final_structure_selection_audit.csv", dtype={"pdb_id": str})
    # only entries with a cached mmCIF can enter scripts 26/27 (coordinate-selected set);
    # discovery-only candidates without coordinates are skipped and listed
    all_pdbs = sorted(audit["pdb_id"].unique())
    pdbs = [p for p in all_pdbs if mmcif_path(p).exists()]
    (w / "logs").mkdir(parents=True, exist_ok=True)
    (w / "logs" / "reconstruction_skipped_uncached_pdbs.txt").write_text("\n".join(sorted(set(all_pdbs) - set(pdbs))) + "\n")
    part = pd.read_csv(man / "final_polymer_partner_annotations.csv", dtype={"pdb_id": str}, keep_default_na=False)
    frozen_meta = {r.partner_polymer_entity_id: r for r in part.drop_duplicates("partner_polymer_entity_id").itertuples()}
    entries, polys, prov = [], [], []
    for pdb in pdbs:
        d = parse_mmcif_gz(mmcif_path(pdb))
        poly_ids = as_list(d.get("_entity_poly.entity_id", []))
        entries.append({"rcsb_id": pdb, "rcsb_entry_container_identifiers": {"polymer_entity_ids": poly_ids}})
        desc = dict(zip(as_list(d["_entity.id"]), as_list(d.get("_entity.pdbx_description", []))))
        seq_len: dict[str, int] = {}
        for e, n in zip(as_list(d.get("_entity_poly_seq.entity_id", [])), as_list(d.get("_entity_poly_seq.num", []))):
            seq_len.setdefault(e, set()).add(n)
        unp: dict[str, list[str]] = {}
        if "_struct_ref.entity_id" in d:
            for e, db, acc in zip(as_list(d["_struct_ref.entity_id"]), as_list(d["_struct_ref.db_name"]),
                                  as_list(d["_struct_ref.pdbx_db_accession"])):
                if db == "UNP" and acc not in unp.get(e, []):
                    unp.setdefault(e, []).append(acc)
        # nonpolymer entities absent from the frozen inventory (entries new to the v2 pool):
        # metadata from the local mmCIF (_pdbx_entity_nonpoly, _chem_comp, _pdbx_entity_instance_feature)
        if "_pdbx_entity_nonpoly.entity_id" in d:
            cc = {}
            if "_chem_comp.id" in d:
                for cid, nm, fm, fw, ty in zip(as_list(d["_chem_comp.id"]), as_list(d.get("_chem_comp.name", [])),
                                               as_list(d.get("_chem_comp.formula", [])), as_list(d.get("_chem_comp.formula_weight", [])),
                                               as_list(d.get("_chem_comp.type", []))):
                    cc[cid] = (nm, fm, fw, ty)
            soi = set()
            if "_pdbx_entity_instance_feature.feature_type" in d:
                for ft, comp in zip(as_list(d["_pdbx_entity_instance_feature.feature_type"]),
                                    as_list(d.get("_pdbx_entity_instance_feature.comp_id", []))):
                    if "SUBJECT OF INVESTIGATION" in ft.upper():
                        soi.add(comp)
            for ent_id, comp in zip(as_list(d["_pdbx_entity_nonpoly.entity_id"]), as_list(d["_pdbx_entity_nonpoly.comp_id"])):
                rid = f"{pdb}_{ent_id}"
                if rid in frozen_np_ids or comp in ("HOH", "DOD"):
                    continue
                nm, fm, fw, ty = cc.get(comp, ("", "", "", ""))
                clean = lambda v: None if v in ("?", ".", "") else v  # noqa: E731
                ents.append({
                    "rcsb_id": rid,
                    "rcsb_nonpolymer_entity_container_identifiers": {"entry_id": pdb, "entity_id": ent_id, "nonpolymer_comp_id": comp},
                    "nonpolymer_comp": {"chem_comp": {"name": clean(nm), "formula": clean(fm),
                                                      "formula_weight": float(fw) if clean(fw) else None, "type": clean(ty)},
                                        "rcsb_chem_comp_info": {"atom_count_heavy": None}},
                    "rcsb_nonpolymer_entity_annotation": [{"type": "SUBJECT_OF_INVESTIGATION"}] if comp in soi else [],
                    "nonpolymer_entity_instances": [],
                })
                np_prov.append({"nonpolymer_entity_id": rid, "source": "local_mmcif",
                                "subject_of_investigation_from_mmcif": comp in soi})
        for e in poly_ids:
            rid = f"{pdb}_{e}"
            mm = {"description": desc.get(e, "") if desc.get(e, "") not in ("?", ".") else "",
                  "sequence_length": len(seq_len.get(e, ())), "uniprot": unp.get(e, [])}
            fz = frozen_meta.get(rid)
            if fz is not None:
                meta = {"description": fz.partner_description,
                        "sequence_length": int(float(fz.partner_sequence_length)) if fz.partner_sequence_length != "" else None,
                        "uniprot": json.loads(fz.partner_mapped_uniprot_ids_json), "source": "frozen_partner_table"}
                prov.append({"polymer_entity_id": rid, "source": "frozen_partner_table",
                             "mmcif_description_equal": mm["description"] == meta["description"],
                             "mmcif_length_equal": mm["sequence_length"] == meta["sequence_length"],
                             "mmcif_uniprot_equal": mm["uniprot"] == meta["uniprot"]})
            else:
                meta = {**mm, "source": "local_mmcif"}
                prov.append({"polymer_entity_id": rid, "source": "local_mmcif", "mmcif_description_equal": None,
                             "mmcif_length_equal": None, "mmcif_uniprot_equal": None})
            polys.append({
                "rcsb_id": rid,
                "rcsb_polymer_entity_container_identifiers": {
                    "entry_id": pdb, "entity_id": e, "asym_ids": [], "auth_asym_ids": [],
                    "reference_sequence_identifiers": [{"database_name": "UniProt", "database_accession": a}
                                                       for a in meta["uniprot"]]},
                "rcsb_polymer_entity": {"pdbx_description": meta["description"]},
                "entity_poly": {"rcsb_sample_sequence_length": meta["sequence_length"]},
            })
    (out / "entry_container_identifiers_batch_9000.json").write_text(json.dumps(
        {"batch_type": "entry", "ids_requested": pdbs, "reconstructed_from": "local mmCIF _entity_poly",
         "response": {"data": {"entries": entries}}}, indent=1))
    (out / "other_polymer_entities_batch_9000.json").write_text(json.dumps(
        {"batch_type": "polymer_entity", "ids_requested": [p["rcsb_id"] for p in polys],
         "reconstructed_from": "frozen final_polymer_partner_annotations.csv, else local mmCIF",
         "response": {"data": {"polymer_entities": polys}}}, indent=1))
    (out / "nonpolymer_entities_batch_9000.json").write_text(json.dumps(
        {"batch_type": "nonpolymer_entity", "ids_requested": [e["rcsb_id"] for e in ents],
         "reconstructed_from": "frozen final_nonpolymer_inventory.csv; entries absent there from local mmCIF",
         "response": {"data": {"nonpolymer_entities": ents}}}, indent=1))
    pd.DataFrame(prov).to_csv(w / "logs" / "reconstructed_polymer_metadata_provenance.csv", index=False)
    pd.DataFrame(np_prov).to_csv(w / "logs" / "reconstructed_nonpolymer_metadata_provenance.csv", index=False)


# ---------------------------------------------------------------------------
# checkpoints
# ---------------------------------------------------------------------------


def checkpoint_a_pool(w: Path) -> dict:
    """Script 18's own build_pool(); no download. Stop if a selected entry is not cached."""
    code = (
        "import sys, importlib.util; from pathlib import Path; import pandas as pd\n"
        "sys.path.insert(0, '.')\n"
        "spec = importlib.util.spec_from_file_location('s18', 'scripts/18_download_mmcif.py')\n"
        "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n"
        "qc = pd.read_csv(m.PRECOORDINATE_QC_CSV, dtype={'pdb_id': str})\n"
        "pool = m.build_pool(qc); pool.to_csv(m.POOL_CSV, index=False)\n"
    )
    run_python(w, code, "18_pool")
    new = pd.read_csv(w / "data" / "manifests" / "stage3b_coordinate_pool.csv", dtype={"pdb_id": str})
    old = pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "stage3b_coordinate_pool.csv", dtype={"pdb_id": str})
    sel_new = set(new[new["coordinate_audit_selected"]]["pdb_id"])
    sel_old = set(old[old["coordinate_audit_selected"]]["pdb_id"])
    missing = sorted(p for p in sel_new if not mmcif_path(p).exists())
    res = {"pool_rows_v1": len(old), "pool_rows_new": len(new), "selected_pdbs_v1": len(sel_old),
           "selected_pdbs_new": len(sel_new), "added": sorted(sel_new - sel_old), "removed": sorted(sel_old - sel_new),
           "selected_but_not_cached": missing,
           "pool_identical_to_v1": new.equals(old)}
    if missing:
        raise SystemExit(f"CHECKPOINT A STOP: selected entries without local mmCIF: {missing}")
    return res


def checkpoint_b_selection(w: Path) -> dict:
    new = pd.read_csv(w / "data" / "manifests" / "final_structure_selection_audit.csv", dtype={"pdb_id": str})
    old = pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "final_structure_selection_audit.csv", dtype={"pdb_id": str})
    key = ["uniprot_id", "pdb_id", "polymer_entity_id"]
    a = new.set_index(key)
    b = old.set_index(key)
    sel_new = set(new[new["stage3b_coordinate_selected"]]["pdb_id"])
    sel_old = set(old[old["stage3b_coordinate_selected"]]["pdb_id"])
    rep_changed = []
    common = a.index.intersection(b.index)
    if "representative_instance_id" in a.columns:
        ra = a.loc[common, "representative_instance_id"].fillna("").astype(str)
        rb = b.loc[common, "representative_instance_id"].fillna("").astype(str)
        diff = ra != rb  # empty/NaN (no coordinates) on both sides is not a change
        rep_changed = [list(i) for i in common[diff.to_numpy()]]
    known = set(json.loads((w / "data" / "raw" / "api" / "rcsb" / "data" / "entry_container_identifiers_batch_9000.json")
                           .read_text())["ids_requested"])
    res = {"selected_pdbs_v1": len(sel_old), "selected_pdbs_new": len(sel_new),
           "added": sorted(sel_new - sel_old), "removed": sorted(sel_old - sel_new),
           "representative_instance_changes": rep_changed,
           "selected_pdbs_without_reconstructed_metadata": sorted(sel_new - known)}
    if res["selected_pdbs_without_reconstructed_metadata"]:
        raise SystemExit("CHECKPOINT B STOP: script 25 would need a network fetch for "
                         f"{res['selected_pdbs_without_reconstructed_metadata']}")
    return res


# ---------------------------------------------------------------------------


def run_python(w: Path, code: str, tag: str) -> None:
    env = {**os.environ, **NO_NETWORK_ENV, "PYTHONHASHSEED": "0"}
    log = w / "logs" / f"{tag}.log"
    with open(log, "w") as fh:
        r = subprocess.run([PY, "-c", code], cwd=w, env=env, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"step {tag} failed; see {log}")


def run_script(w: Path, script: str) -> None:
    env = {**os.environ, **NO_NETWORK_ENV, "PYTHONHASHSEED": "0"}
    log = w / "logs" / f"{script[:-3]}.log"
    with open(log, "w") as fh:
        r = subprocess.run([PY, f"scripts/{script}"], cwd=w, env=env, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"{script} failed (exit {r.returncode}); see {log}")


def export(variant: str) -> dict:
    """Copy every file the pipeline (re)wrote in the sandbox -- and nothing stale --
    to the versioned export location (data/rebuild_v2/..., reports/rebuild_v2/...)."""
    w = work_dir(variant)
    t0 = float((w / "logs" / "sandbox_built_at").read_text())
    exported = []
    for d, dst_root in [("data/manifests", EXPORT_DATA / "manifests"), ("data/processed", EXPORT_DATA / "processed"),
                        ("data/analysis", EXPORT_DATA / "analysis"), ("reports/tables", EXPORT_REPORTS / "tables")]:
        for f in sorted((w / d).rglob("*")):
            if f.is_file() and f.stat().st_mtime > t0:
                dst = dst_root / f.relative_to(w / d)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, dst)
                exported.append(str(dst.relative_to(PROJECT_ROOT)))
    (EXPORT_REPORTS / "foundation_export_manifest.json").write_text(json.dumps(
        {"variant": variant, "exported_files": exported,
         "checkpoints": json.loads((w / "logs" / "checkpoints.json").read_text())}, indent=2, default=str) + "\n")
    return {"n_exported": len(exported)}


def run(variant: str, from_step: str | None) -> dict:
    w = work_dir(variant) if from_step else build_sandbox(variant)
    if not from_step:
        import time
        (w / "logs" / "sandbox_built_at").write_text(repr(time.time()))
    started = from_step is None
    if variant == "v2":
        apply_sandbox_patches(w)
    cp_file = w / "logs" / "checkpoints.json"
    checkpoints = json.loads(cp_file.read_text()) if (from_step and cp_file.exists()) else {}
    for num, script in STEPS:
        if variant != "v2" and num in V2_ONLY_STEPS:
            continue
        if not started:
            started = num == from_step
            if not started:
                continue
        print(f"[{variant}] step {num}: {script}", flush=True)
        if script == "CHECKPOINT_A_POOL":
            checkpoints["A_pool"] = checkpoint_a_pool(w)
            cp_file.write_text(json.dumps(checkpoints, indent=2, default=str))
        elif script == "CHECKPOINT_B_SELECTION":
            checkpoints["B_selection"] = checkpoint_b_selection(w)
            cp_file.write_text(json.dumps(checkpoints, indent=2, default=str))
        else:
            run_script(w, script)
    (w / "logs" / "checkpoints.json").write_text(json.dumps(checkpoints, indent=2, default=str))
    return checkpoints


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--variant", choices=["v1_repro", "v2"], required=True)
    ap.add_argument("--from-step", default=None, help="resume an existing sandbox at this step number")
    ap.add_argument("--export", action="store_true", help="export rewritten files (v2 only)")
    ap.add_argument("--export-only", action="store_true")
    ap.add_argument("--run-single", default=None,
                    help="run one additional legacy script in an existing v2 sandbox (patches applied first)")
    args = ap.parse_args()
    if args.run_single:
        w = work_dir(args.variant)
        if args.variant == "v2":
            apply_sandbox_patches(w)
        run_script(w, args.run_single)
        print(f"ran {args.run_single} in {w}")
        return
    if not args.export_only:
        print(json.dumps(run(args.variant, args.from_step), indent=2, default=str))
    if (args.export or args.export_only) and args.variant == "v2":
        print(json.dumps(export(args.variant)))


if __name__ == "__main__":
    main()
