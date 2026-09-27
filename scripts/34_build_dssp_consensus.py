"""
Structural-segment alignment, Part 3: per-receptor secondary-structure
consensus evidence.

No DSSP binary is installed in this environment (checked: mkdssp/dssp are
not on PATH) and none is installed here. Secondary-structure evidence
instead comes from the mmCIF `_struct_conf` (helices) and
`_struct_sheet_range` (beta strands) categories already deposited in every
cached mmCIF file -- the same underlying evidence (depositor/PDB-pipeline
secondary-structure assignment, typically DSSP- or PROMOTIF-derived at
deposition time), read directly with no new software. See
scripts/utils/mmcif_parse.py:get_secondary_structure_ranges.

Uses every representative coordinate instance (1840; independent of B-factor
normalizability -- secondary structure does not depend on B-factors), one
mmCIF parse per unique PDB entry.

Reads:
    data/processed/structures.csv
    data/analysis/normalized_lbd_bfactors_all.parquet (residue-level identity
        + entity_label_seq_id <-> canonical_uniprot_position mapping only;
        B-factor/Z columns are not used here)
    data/raw/mmcif/{PDB_ID}.cif.gz (already-cached; no download)

Writes:
    data/analysis/receptor_dssp_consensus.parquet
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.mmcif_parse import get_secondary_structure_ranges, parse_mmcif_gz  # noqa: E402

PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
RAW_MMCIF = PROJECT_ROOT / "data" / "raw" / "mmcif"

MASTER_CSV = PROCESSED / "structures.csv"
NORMALIZED_PARQUET = ANALYSIS / "normalized_lbd_bfactors_all.parquet"
OUT_PARQUET = ANALYSIS / "receptor_dssp_consensus.parquet"

STATES = ["HELIX", "SHEET", "COIL"]


def assign_state(label_seq_id: int, ranges: list[dict]) -> str:
    for r in ranges:
        if r["sse_type"] == "HELIX" and r["beg_label_seq_id"] <= label_seq_id <= r["end_label_seq_id"]:
            return "HELIX"
    for r in ranges:
        if r["sse_type"] == "SHEET" and r["beg_label_seq_id"] <= label_seq_id <= r["end_label_seq_id"]:
            return "SHEET"
    return "COIL"


def main() -> None:
    master = pd.read_csv(MASTER_CSV, dtype={"pdb_id": str})
    coord = master[master["has_stage3b_coordinate_data"]]
    assert len(coord) == 1840

    res = pd.read_parquet(
        NORMALIZED_PARQUET,
        columns=["instance_id", "uniprot_id", "nr_code", "pdb_id", "label_asym_id",
                 "canonical_uniprot_position", "canonical_residue_name", "entity_label_seq_id",
                 "is_construct_mapped", "ca_record_present"],
    )
    rep_ids = set(coord["representative_instance_id"])
    rep = res[res["instance_id"].isin(rep_ids)].copy()
    assert rep["instance_id"].nunique() == len(rep_ids) == 1840

    rows = []
    by_pdb = rep.groupby("pdb_id", sort=False)
    for n_done, (pdb_id, g) in enumerate(by_pdb, start=1):
        d = parse_mmcif_gz(RAW_MMCIF / f"{pdb_id}.cif.gz")
        ranges_by_asym: dict[str, list[dict]] = {}
        for r in get_secondary_structure_ranges(d):
            ranges_by_asym.setdefault(r["label_asym_id"], []).append(r)

        for instance_id, gi in g.groupby("instance_id", sort=False):
            asym = gi["label_asym_id"].iloc[0]
            ranges = ranges_by_asym.get(asym, [])
            observed = gi[gi["is_construct_mapped"] & gi["ca_record_present"]]
            for _, r in observed.iterrows():
                seq_id = int(r["entity_label_seq_id"])
                rows.append(
                    {
                        "uniprot_id": r["uniprot_id"], "canonical_uniprot_position": int(r["canonical_uniprot_position"]),
                        "instance_id": instance_id, "pdb_id": pdb_id, "state": assign_state(seq_id, ranges),
                    }
                )
        if n_done % 300 == 0:
            print(f"  parsed {n_done}/{by_pdb.ngroups} PDB entries")

    obs = pd.DataFrame(rows)
    print(f"Observed (residue x representative-instance) rows with secondary-structure evidence: {len(obs)}")

    names = rep.drop_duplicates(["uniprot_id", "canonical_uniprot_position"]).set_index(
        ["uniprot_id", "canonical_uniprot_position"])["canonical_residue_name"]
    nr_codes = rep.drop_duplicates("uniprot_id").set_index("uniprot_id")["nr_code"]

    key = ["uniprot_id", "canonical_uniprot_position"]
    counts = obs.groupby(key + ["state"]).size().unstack(fill_value=0).reindex(columns=STATES, fill_value=0)
    n_pdb = obs.groupby(key)["pdb_id"].nunique().rename("n_pdb_observed")
    n_inst = obs.groupby(key)["instance_id"].nunique().rename("n_dssp_observed")

    out = counts.join(n_pdb).join(n_inst).reset_index()
    total = out[STATES].sum(axis=1)
    out["consensus_dssp"] = out[STATES].idxmax(axis=1)
    out["consensus_fraction"] = out[STATES].max(axis=1) / total
    for s in STATES:
        out[f"{s.lower()}_fraction"] = out[s] / total
    out["canonical_residue_name"] = [names.get((u, p), "") for u, p in zip(out["uniprot_id"], out["canonical_uniprot_position"])]
    out["nr_code"] = out["uniprot_id"].map(nr_codes)

    out = out.rename(columns={"HELIX": "n_helix_votes", "SHEET": "n_sheet_votes", "COIL": "n_coil_votes"})
    cols = ["uniprot_id", "nr_code", "canonical_uniprot_position", "canonical_residue_name",
            "n_pdb_observed", "n_dssp_observed", "n_helix_votes", "n_sheet_votes", "n_coil_votes",
            "consensus_dssp", "consensus_fraction", "helix_fraction", "sheet_fraction", "coil_fraction"]
    out = out[cols].sort_values(["uniprot_id", "canonical_uniprot_position"]).reset_index(drop=True)

    assert not out.duplicated(key).any()
    # n_dssp_observed (instances) >= n_pdb_observed (unique PDBs) is expected: a small number of
    # receptors have >1 representative polymer-entity candidate in the same PDB entry (e.g. a
    # heterodimer construct deposited as two separate entities of the same receptor).
    assert (out["n_dssp_observed"] >= out["n_pdb_observed"]).all()
    assert np.isclose(out["helix_fraction"] + out["sheet_fraction"] + out["coil_fraction"], 1.0).all()

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT_PARQUET, index=False)
    print(f"Wrote {OUT_PARQUET.relative_to(PROJECT_ROOT)} ({len(out)} rows, {out['uniprot_id'].nunique()} receptors)")
    print(out["consensus_dssp"].value_counts().to_string())


if __name__ == "__main__":
    main()
