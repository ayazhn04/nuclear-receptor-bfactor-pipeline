"""Compute isolated-chain SASA for the corrected v2 representative instances.

Offline by default; --download restores missing RCSB files into a NEW cache.
Only files matching the frozen decompressed SHA256 are accepted. No upstream
table/cache is modified. Run with biotite, biopython, pandas and pyarrow installed.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/role_d/sasa_iqr"
CACHE = ROOT / "data/raw/external/role_d_sasa_mmcif"
SOURCE = ROOT / "data/rebuild_v2/processed/lbd_residue_bfactors_all.parquet"
# Tien et al. (2013), theoretical maximum ASA, Angstrom squared.
MAX_ASA = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    [129, 274, 195, 193, 167, 225, 223, 104, 224, 197, 201, 236, 224, 240, 159, 155, 172, 285, 263, 174],
))
SIDECHAIN = {
    "ALA": "CB", "ARG": "CB CG CD NE CZ NH1 NH2", "ASN": "CB CG OD1 ND2",
    "ASP": "CB CG OD1 OD2", "CYS": "CB SG", "GLN": "CB CG CD OE1 NE2",
    "GLU": "CB CG CD OE1 OE2", "GLY": "", "HIS": "CB CG ND1 CD2 CE1 NE2",
    "ILE": "CB CG1 CG2 CD1", "LEU": "CB CG CD1 CD2", "LYS": "CB CG CD CE NZ",
    "MET": "CB CG SD CE", "PHE": "CB CG CD1 CD2 CE1 CE2 CZ", "PRO": "CB CG CD",
    "SER": "CB OG", "THR": "CB OG1 CG2", "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "CB CG1 CG2",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_hashes() -> dict[str, str]:
    paths = [ROOT / "data/manifests/mmcif_download_manifest.csv",
             ROOT / "data/rebuild_v2/inputs/mmcif_download_additions_v2.csv"]
    result = {}
    for path in paths:
        if path.exists():
            data = pd.read_csv(path)
            if {"pdb_id", "sha256_decompressed"} <= set(data):
                result.update(data.set_index("pdb_id").sha256_decompressed.dropna().to_dict())
    return result


def restore(task):
    pdb, expected, download = task
    paths = [ROOT / f"data/raw/mmcif/{pdb}.cif.gz",
             ROOT / f"data/interim/raw_mmcif_v2_additions/{pdb}.cif.gz",
             CACHE / f"{pdb}.cif.gz"]
    path = next((p for p in paths if p.exists()), paths[-1])
    row = {"pdb_id": pdb, "url": f"https://files.rcsb.org/download/{pdb}.cif.gz",
           "expected_sha256_decompressed": expected, "status": "MISSING_CACHE"}
    try:
        if not path.exists() and download:
            import requests
            response = requests.get(row["url"], timeout=90)
            response.raise_for_status()
            # Validate gzip before retaining any download.
            gzip.decompress(response.content)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(response.content)
        if path.exists():
            digest = hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest()
            row.update(sha256_decompressed=digest, sha256_gzip=sha(path),
                       status="VERIFIED" if expected and digest == expected else "FROZEN_HASH_MISMATCH_OR_UNAVAILABLE",
                       path=str(path.relative_to(ROOT)))
    except Exception as exc:
        row.update(status="DOWNLOAD_OR_VALIDATION_FAILED", error=str(exc))
    return row


def select_atoms(atoms: pd.DataFrame, preferred: dict[int, str]) -> pd.DataFrame:
    """One coherent nonblank conformer per residue, shared blank atoms retained."""
    atoms = atoms.loc[(atoms.occupancy > 0) & ~atoms.element.isin(["H", "D"])].copy()
    pieces = []
    for seq, group in atoms.groupby("seq", sort=False):
        nonblank = group.loc[~group.alt.isin([".", "?", ""])]
        choice = preferred.get(int(seq), "")
        if choice not in set(nonblank.alt):
            ranking = nonblank.groupby("alt").occupancy.mean().sort_index().sort_values(ascending=False, kind="stable")
            choice = ranking.index[0] if len(ranking) else ""
        group = group.loc[group.alt.isin([".", "?", "", choice])]
        group = group.sort_values(["occupancy", "alt"], ascending=[False, True], kind="stable").drop_duplicates("atom")
        pieces.append(group)
    return pd.concat(pieces, ignore_index=True) if pieces else atoms.iloc[:0]


def calculate_chain(atoms: pd.DataFrame, mapping: pd.DataFrame, points: int) -> pd.DataFrame:
    import biotite.structure as bs
    from biotite.structure.info import vdw_radius_protor, vdw_radius_single

    preferred = dict(zip(mapping.entity_label_seq_id.dropna().astype(int),
                         mapping.loc[mapping.entity_label_seq_id.notna(), "selected_ca_altloc"].fillna("")))
    atoms = select_atoms(atoms, preferred)
    array = bs.AtomArray(len(atoms))
    array.coord = atoms[["x", "y", "z"]].to_numpy(dtype=float)
    array.res_id = atoms.seq.to_numpy(dtype=int)
    array.res_name = atoms.res.to_numpy(str)
    array.atom_name = atoms.atom.to_numpy(str)
    array.element = atoms.element.to_numpy(str)
    array.chain_id[:] = "A"
    radii = []
    fallback = 0
    for row in atoms.itertuples():
        try:
            radius = vdw_radius_protor(row.res, row.atom)
        except (KeyError, ValueError):
            radius = None
        if radius is None or not np.isfinite(radius):
            radius = vdw_radius_single(row.element)
            fallback += 1
        if radius is None or not np.isfinite(radius):
            raise ValueError(f"No radius for {row.element}/{row.res}/{row.atom}")
        radii.append(radius)
    area = bs.sasa(array, probe_radius=1.4, point_number=points, vdw_radii=np.array(radii), ignore_ions=False)
    atoms["atom_sasa"] = area
    groups = {int(k): v for k, v in atoms.groupby("seq")}
    rows = []
    for r in mapping.itertuples():
        row = {"uniprot_id": r.uniprot_id, "pdb_id": r.pdb_id, "instance_id": r.instance_id,
               "label_asym_id": r.label_asym_id, "model_num": r.model_num,
               "canonical_uniprot_position": r.canonical_uniprot_position,
               "entity_label_seq_id": r.entity_label_seq_id, "sasa_status": "NO_USABLE_CA",
               "chain_n_heavy_atoms": len(atoms), "chain_radius_fallback_atoms": fallback}
        group = groups.get(int(r.entity_label_seq_id)) if pd.notna(r.entity_label_seq_id) else None
        if r.ca_bfactor_usable and group is not None:
            names = set(group.atom)
            observed = group.res.iloc[0]
            row.update(observed_residue_name=observed, canonical_residue_name=r.canonical_residue_name)
            if observed not in MAX_ASA or group.res.nunique() != 1:
                row["sasa_status"] = "UNSUPPORTED_OR_MIXED_RESIDUE"
            elif not (set("N CA C O".split()) | set(SIDECHAIN[observed].split())) <= names:
                row["sasa_status"] = "INCOMPLETE_HEAVY_ATOMS"
            else:
                total = float(group.atom_sasa.sum())
                row.update(sasa_status="OK", sasa_a2=total, rsa=total / MAX_ASA[observed],
                           residue_matches_canonical=observed == r.canonical_residue_name)
        rows.append(row)
    return pd.DataFrame(rows)


def calculate_pdb(task):
    pdb, path, mapping, points = task
    from biotite.structure.io.pdbx import CIFFile
    try:
        with gzip.open(ROOT / path, "rt") as stream:
            data = CIFFile.read(stream).block["atom_site"]
        fields = dict(chain="label_asym_id", model="pdbx_PDB_model_num", seq="label_seq_id",
                      atom="label_atom_id", alt="label_alt_id", res="label_comp_id", element="type_symbol",
                      occupancy="occupancy", x="Cartn_x", y="Cartn_y", z="Cartn_z")
        atoms = pd.DataFrame({k: data[v].as_array(str) for k, v in fields.items()})
        atoms.seq = pd.to_numeric(atoms.seq, errors="coerce")
        atoms = atoms.loc[atoms.seq.notna()].copy()  # polymer only; no ligand/water
        atoms.occupancy = pd.to_numeric(atoms.occupancy)
        parts = []
        for (chain, model), group in mapping.groupby(["label_asym_id", "model_num"]):
            selected = atoms.loc[atoms.chain.eq(str(chain)) & atoms.model.eq(str(model))]
            parts.append(calculate_chain(selected, group, points))
        return pd.concat(parts, ignore_index=True), {"pdb_id": pdb, "status": "OK"}
    except Exception as exc:
        return pd.DataFrame(), {"pdb_id": pdb, "status": "CALCULATION_FAILED", "error": str(exc)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--points", type=int, default=1000)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_parquet(SOURCE)
    mapping = mapping.loc[mapping.is_representative_instance].copy()
    hashes = frozen_hashes()
    pdbs = sorted(mapping.loc[mapping.ca_bfactor_usable, "pdb_id"].unique())
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        downloads = list(pool.map(restore, [(pdb, hashes.get(pdb), args.download) for pdb in pdbs]))
    audit = pd.DataFrame(downloads)
    audit.to_csv(OUT / "coordinate_cache_audit.csv", index=False)
    verified = audit.loc[audit.status.eq("VERIFIED")]
    print(f"Verified frozen coordinate files: {len(verified)}/{len(audit)}", flush=True)
    tasks = [(r.pdb_id, r.path, mapping.loc[mapping.pdb_id.eq(r.pdb_id)], args.points) for r in verified.itertuples()]
    frames, statuses = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i, (frame, status) in enumerate(pool.map(calculate_pdb, tasks), 1):
            frames.append(frame)
            statuses.append(status)
            if i % 100 == 0:
                print(f"SASA completed: {i}/{len(tasks)}", flush=True)
    pd.DataFrame(statuses).to_csv(OUT / "calculation_audit.csv", index=False)
    if not frames or all(frame.empty for frame in frames):
        raise RuntimeError("No SASA calculated; inspect coordinate_cache_audit.csv / calculation_audit.csv")
    residues = pd.concat(frames, ignore_index=True)
    residues.to_parquet(OUT / "representative_instance_sasa.parquet", index=False)
    usable = residues.loc[residues.sasa_status.eq("OK")]
    keys = ["uniprot_id", "pdb_id", "canonical_uniprot_position"]
    entries = usable.groupby(keys, as_index=False).agg(
        sasa_a2=("sasa_a2", "median"), rsa=("rsa", "median"), n_instances=("instance_id", "nunique"))
    entries.to_parquet(OUT / "pdb_position_sasa.parquet", index=False)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Shrake-Rupley; isolated deposited polymer chain; all observed heavy atoms; no ligand/water/other chain/symmetry",
        "radii": "ProtOr; Single-element fallback for unclassified atoms (counts audited)",
        "probe_radius_angstrom": 1.4, "sphere_points": args.points,
        "rsa_reference": "Tien et al. 2013 theoretical maximum ASA; DOI 10.1371/journal.pone.0080635; observed residue type; not clipped at 1",
        "n_requested_pdb": len(pdbs), "n_verified_pdb": len(verified),
        "n_calculated_pdb": sum(s["status"] == "OK" for s in statuses),
        "instance_residue_status_counts": residues.sasa_status.value_counts().to_dict(),
        "input_sha256": {str(SOURCE.relative_to(ROOT)): sha(SOURCE)},
        "script_sha256": sha(Path(__file__)),
        "versions": {p: importlib.metadata.version(p) for p in ["biotite", "biopython", "numpy", "pandas", "pyarrow"]},
    }
    (OUT / "sasa_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
