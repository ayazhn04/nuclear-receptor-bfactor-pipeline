"""
Structural-segment alignment, Parts 6-9: segment coverage QC, entry-level and
receptor-level segment flexibility features, and the primary/sensitivity
feature matrices.

Aggregation hierarchy (supervisor-approved, unchanged):
    raw B -> per-instance Z (already done, scripts/32) -> chain/instance ->
    PDB entry by MEDIAN (already done: data/analysis/pdb_position_bfactor_profiles.parquet)
    -> canonical structural segment (this script) -> receptor by MEDIAN (primary)
    / MEAN (sensitivity)

Segment definitions come only from named entries of
data/analysis/nr_lbd_structural_segment_map.csv (AMBIGUOUS/UNASSIGNED
positions are never pooled into a segment). Coverage is metadata and is kept
in its own columns/matrix, never mixed into a flexibility (Z) column, and
missing segment values are left as NaN, never 0 (0 is a meaningful z-score).

Reads:
    data/processed/structures.csv
    data/analysis/normalized_lbd_bfactors_all.parquet
    data/analysis/pdb_position_bfactor_profiles.parquet
    data/analysis/nr_lbd_structural_segment_map.csv
    data/manifests/nr_lbd_reference.csv

Writes:
    data/analysis/pdb_segment_coverage.parquet
    data/analysis/pdb_segment_flexibility_features.parquet
    data/analysis/receptor_segment_flexibility_profiles.parquet
    data/analysis/receptor_segment_feature_matrix_median.parquet
    data/analysis/receptor_segment_feature_matrix_mean.parquet
    data/analysis/pdb_segment_feature_matrix_median.parquet
    data/analysis/pdb_segment_feature_matrix_mean.parquet
    data/analysis/receptor_segment_coverage_matrix.parquet
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"

MASTER_CSV = PROCESSED / "structures.csv"
NORMALIZED_PARQUET = ANALYSIS / "normalized_lbd_bfactors_all.parquet"
PDB_POSITION_PARQUET = ANALYSIS / "pdb_position_bfactor_profiles.parquet"
SEGMENT_MAP_CSV = ANALYSIS / "nr_lbd_structural_segment_map.csv"
LBD_REFERENCE_CSV = MANIFESTS / "nr_lbd_reference.csv"

COVERAGE_PARQUET = ANALYSIS / "pdb_segment_coverage.parquet"
PDB_FEATURES_PARQUET = ANALYSIS / "pdb_segment_flexibility_features.parquet"
RECEPTOR_PROFILE_PARQUET = ANALYSIS / "receptor_segment_flexibility_profiles.parquet"
RECEPTOR_MATRIX_MEDIAN = ANALYSIS / "receptor_segment_feature_matrix_median.parquet"
RECEPTOR_MATRIX_MEAN = ANALYSIS / "receptor_segment_feature_matrix_mean.parquet"
PDB_MATRIX_MEDIAN = ANALYSIS / "pdb_segment_feature_matrix_median.parquet"
PDB_MATRIX_MEAN = ANALYSIS / "pdb_segment_feature_matrix_mean.parquet"
RECEPTOR_COVERAGE_MATRIX = ANALYSIS / "receptor_segment_coverage_matrix.parquet"


def load_segment_lookup(seg_map: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    named = seg_map[seg_map["structural_segment"] != ""].copy()
    seg_len = named.groupby(["uniprot_id", "structural_segment"]).size().rename("n_expected_residues")
    seg_class = named.groupby(["uniprot_id", "structural_segment"])["structural_class"].first()
    seg_order = named.groupby(["uniprot_id", "structural_segment"])["segment_order"].first()
    defs = pd.concat([seg_len, seg_class, seg_order], axis=1).reset_index()
    pos_to_seg = named.set_index(["uniprot_id", "canonical_uniprot_position"])["structural_segment"].to_dict()
    return defs, pos_to_seg


def build_coverage(master: pd.DataFrame, norm: pd.DataFrame, seg_map: pd.DataFrame, defs: pd.DataFrame) -> pd.DataFrame:
    reps = master.loc[master["has_stage3b_coordinate_data"], ["uniprot_id", "pdb_id", "representative_instance_id"]]
    n = norm[norm["instance_id"].isin(set(reps["representative_instance_id"])) & norm["is_construct_mapped"]].merge(
        seg_map[seg_map["structural_segment"] != ""][["uniprot_id", "canonical_uniprot_position", "structural_segment"]],
        on=["uniprot_id", "canonical_uniprot_position"], how="inner",
    )
    n = n.merge(reps.rename(columns={"representative_instance_id": "instance_id"}), on=["uniprot_id", "instance_id"], how="inner")

    rows = []
    for (uid, pdb, seg), g in n.groupby(["uniprot_id", "pdb_id", "structural_segment"]):
        exp = int(defs.loc[(defs["uniprot_id"] == uid) & (defs["structural_segment"] == seg), "n_expected_residues"].iloc[0])
        observed_positions = set(g.loc[g["ca_positive_occupancy"], "canonical_uniprot_position"])
        rows.append({"uniprot_id": uid, "pdb_id": pdb, "structural_segment": seg,
                     "n_expected_residues": exp, "n_observed_residues": len(observed_positions)})
    cov = pd.DataFrame(rows)
    cov["segment_coverage"] = cov["n_observed_residues"] / cov["n_expected_residues"]
    for thr in (0.50, 0.80, 0.90):
        cov[f"coverage_ge_0{int(thr * 100)}"] = cov["segment_coverage"] >= thr
    return cov


def build_pdb_features(master: pd.DataFrame, norm: pd.DataFrame, pdb_pos: pd.DataFrame,
                        seg_map: pd.DataFrame, defs: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    named = seg_map[seg_map["structural_segment"] != ""][["uniprot_id", "canonical_uniprot_position", "structural_segment"]]
    p = pdb_pos.merge(named, on=["uniprot_id", "canonical_uniprot_position"], how="inner")

    z = p.groupby(["uniprot_id", "pdb_id", "structural_segment"])["b_factor_z"].agg(
        mean_z="mean", median_z="median", sd_z=lambda s: s.std(ddof=0),
        q25_z=lambda s: s.quantile(0.25), q75_z=lambda s: s.quantile(0.75),
    ).reset_index()
    z["iqr_z"] = z["q75_z"] - z["q25_z"]

    state = p.groupby(["uniprot_id", "pdb_id", "structural_segment"])["ligand_state"].agg(
        lambda s: s.iloc[0] if s.nunique() == 1 else "AMBIGUOUS"
    ).rename("ligand_state").reset_index()

    reps = master.loc[master["has_stage3b_coordinate_data"], ["uniprot_id", "pdb_id", "representative_instance_id"]]
    n = norm[norm["instance_id"].isin(set(reps["representative_instance_id"])) & norm["ca_positive_occupancy"]].merge(
        named, on=["uniprot_id", "canonical_uniprot_position"], how="inner"
    ).merge(reps.rename(columns={"representative_instance_id": "instance_id"}), on=["uniprot_id", "instance_id"], how="inner")
    n_inst = n.groupby(["uniprot_id", "pdb_id", "structural_segment"])["instance_id"].nunique().rename("n_contributing_instances").reset_index()

    out = z.merge(state, on=["uniprot_id", "pdb_id", "structural_segment"], how="left")
    out = out.merge(n_inst, on=["uniprot_id", "pdb_id", "structural_segment"], how="left")
    out = out.merge(coverage[["uniprot_id", "pdb_id", "structural_segment", "n_expected_residues", "n_observed_residues", "segment_coverage"]],
                     on=["uniprot_id", "pdb_id", "structural_segment"], how="left")
    out = out.merge(defs[["uniprot_id", "structural_segment", "structural_class"]], on=["uniprot_id", "structural_segment"], how="left")
    meta = master.drop_duplicates("uniprot_id")[["uniprot_id", "nr_code", "common_name", "group"]]
    out = out.merge(meta, on="uniprot_id", how="left")
    cols = ["uniprot_id", "nr_code", "common_name", "group", "pdb_id", "structural_segment", "structural_class",
            "n_expected_residues", "n_observed_residues", "segment_coverage",
            "mean_z", "median_z", "sd_z", "q25_z", "q75_z", "iqr_z", "ligand_state", "n_contributing_instances"]
    out = out[cols].sort_values(["uniprot_id", "pdb_id", "structural_segment"]).reset_index(drop=True)
    assert not out.duplicated(["uniprot_id", "pdb_id", "structural_segment"]).any()
    return out


def _agg_state(df: pd.DataFrame, state: str, keys: list[str]) -> pd.DataFrame:
    d = df[df["ligand_state"] == state]
    g = d.groupby(keys)["median_z"]
    return pd.DataFrame({f"{state.lower()}_n_pdb": g.count(), f"{state.lower()}_median_of_entry_median_z": g.median(),
                          f"{state.lower()}_mean_of_entry_median_z": g.mean()})


def build_receptor_profile(pdb_feat: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    keys = ["uniprot_id", "structural_segment"]
    n_total = meta.set_index("uniprot_id")["n_pdbs_with_coordinate_data"]
    g_med = pdb_feat.groupby(keys)["median_z"]
    g_mean = pdb_feat.groupby(keys)["mean_z"]
    base = pd.DataFrame(
        {
            "n_pdb_entries_with_segment": g_med.count(),
            "median_of_entry_median_z": g_med.median(), "mean_of_entry_median_z": g_med.mean(),
            "sd_of_entry_median_z": g_med.std(ddof=0), "q25": g_med.quantile(0.25), "q75": g_med.quantile(0.75),
            "median_of_entry_mean_z": g_mean.median(), "mean_of_entry_mean_z": g_mean.mean(),
            "median_segment_coverage": pdb_feat.groupby(keys)["segment_coverage"].median(),
        }
    )
    base["iqr"] = base["q75"] - base["q25"]
    apo = _agg_state(pdb_feat, "APO", keys)
    holo = _agg_state(pdb_feat, "HOLO", keys)
    out = base.join(apo, how="left").join(holo, how="left").reset_index()
    out["n_pdb_entries_total"] = out["uniprot_id"].map(n_total)
    out["pdb_support_fraction"] = out["n_pdb_entries_with_segment"] / out["n_pdb_entries_total"]
    for c in ["apo_n_pdb", "holo_n_pdb"]:
        out[c] = out[c].fillna(0).astype(int)
    out = out.merge(meta[["uniprot_id", "nr_code", "common_name", "group"]], on="uniprot_id", how="left")
    seg_class = pdb_feat.drop_duplicates(["uniprot_id", "structural_segment"])[["uniprot_id", "structural_segment", "structural_class"]]
    out = out.merge(seg_class, on=["uniprot_id", "structural_segment"], how="left")
    cols = ["uniprot_id", "nr_code", "common_name", "group", "structural_segment", "structural_class",
            "n_pdb_entries_total", "n_pdb_entries_with_segment", "pdb_support_fraction", "median_segment_coverage",
            "median_of_entry_median_z", "mean_of_entry_median_z", "sd_of_entry_median_z", "q25", "q75", "iqr",
            "median_of_entry_mean_z", "mean_of_entry_mean_z",
            "apo_n_pdb", "apo_median_of_entry_median_z", "apo_mean_of_entry_median_z",
            "holo_n_pdb", "holo_median_of_entry_median_z", "holo_mean_of_entry_median_z"]
    return out[cols]


def main() -> None:
    master = pd.read_csv(MASTER_CSV, dtype={"pdb_id": str})
    norm = pd.read_parquet(
        NORMALIZED_PARQUET,
        columns=["instance_id", "uniprot_id", "canonical_uniprot_position", "is_construct_mapped", "ca_positive_occupancy"],
    )
    pdb_pos = pd.read_parquet(PDB_POSITION_PARQUET)
    seg_map = pd.read_csv(SEGMENT_MAP_CSV, keep_default_na=False)
    defs, _ = load_segment_lookup(seg_map)
    print(f"Declared (named) segments across all receptors: {len(defs)} (uniprot_id x structural_segment pairs)")

    print("Part 6: PDB x segment coverage...")
    coverage = build_coverage(master, norm, seg_map, defs)
    coverage.to_parquet(COVERAGE_PARQUET, index=False)
    print(f"  wrote {COVERAGE_PARQUET.relative_to(PROJECT_ROOT)} ({len(coverage)} rows)")

    print("Part 7: PDB-entry segment flexibility features...")
    pdb_feat = build_pdb_features(master, norm, pdb_pos, seg_map, defs, coverage)
    pdb_feat.to_parquet(PDB_FEATURES_PARQUET, index=False)
    print(f"  wrote {PDB_FEATURES_PARQUET.relative_to(PROJECT_ROOT)} ({len(pdb_feat)} rows)")

    print("Part 8: receptor-level segment profiles...")
    meta = master.drop_duplicates("uniprot_id")[["uniprot_id", "nr_code", "common_name", "group"]].copy()
    n_pdb_total = master.loc[master["has_stage3b_coordinate_data"]].groupby("uniprot_id")["pdb_id"].nunique()
    meta["n_pdbs_with_coordinate_data"] = meta["uniprot_id"].map(n_pdb_total).fillna(0).astype(int)
    rec_profile = build_receptor_profile(pdb_feat, meta)
    rec_profile.to_parquet(RECEPTOR_PROFILE_PARQUET, index=False)
    print(f"  wrote {RECEPTOR_PROFILE_PARQUET.relative_to(PROJECT_ROOT)} ({len(rec_profile)} rows)")

    print("Part 9: feature matrices (missing = NaN, never 0)...")
    rec_med = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="median_of_entry_median_z")
    rec_mean = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="mean_of_entry_mean_z")
    pdb_med = pdb_feat.pivot_table(index=["uniprot_id", "pdb_id"], columns="structural_segment", values="median_z")
    pdb_mean = pdb_feat.pivot_table(index=["uniprot_id", "pdb_id"], columns="structural_segment", values="mean_z")
    rec_cov = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="median_segment_coverage")

    rec_med.reset_index().to_parquet(RECEPTOR_MATRIX_MEDIAN, index=False)
    rec_mean.reset_index().to_parquet(RECEPTOR_MATRIX_MEAN, index=False)
    pdb_med.reset_index().to_parquet(PDB_MATRIX_MEDIAN, index=False)
    pdb_mean.reset_index().to_parquet(PDB_MATRIX_MEAN, index=False)
    rec_cov.reset_index().to_parquet(RECEPTOR_COVERAGE_MATRIX, index=False)
    print(f"  receptor x segment median matrix: {rec_med.shape}")
    print(f"  receptor x segment mean matrix:   {rec_mean.shape}")
    print(f"  pdb x segment median matrix:      {pdb_med.shape}")
    print(f"  pdb x segment mean matrix:        {pdb_mean.shape}")
    print(f"  receptor x segment coverage matrix (values are coverage, never a flexibility feature): {rec_cov.shape}")
    assert set(rec_med.columns) == set(rec_mean.columns) == set(rec_cov.columns)
    assert set(pdb_med.columns) == set(pdb_mean.columns)


if __name__ == "__main__":
    main()
