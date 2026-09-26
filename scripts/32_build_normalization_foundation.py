"""
Analysis foundation, Parts 3-8: within-instance Ca B-factor normalization and
descriptive receptor x position profiles, built from the FULL retained dataset.

QC thresholds (coverage, resolution, R-free, method, taxonomy, ligand state,
previous strict-subset membership) are metadata only: nothing is filtered out
anywhere in this script.

Inputs (tracked, no network, no raw mmCIF):
    data/processed/structures.csv                    (master, 2072 rows)
    data/processed/lbd_residue_bfactors_all.parquet  (771,542 residue rows)
    data/manifests/nr_lbd_reference.csv              (48 receptors)

Outputs:
    data/analysis/normalized_lbd_bfactors_all.parquet
    data/analysis/instance_bfactor_summary.csv
    data/analysis/pdb_position_bfactor_profiles.parquet
    data/analysis/receptor_position_flexibility_profiles.parquet
    data/analysis/receptor_flexibility_summary.csv
    data/processed/structures.csv   (ONLY adds representative_* descriptive columns;
                                     rows, order and all existing columns unchanged)

Normalization: z = (B - mean_instance_B) / population_SD_instance_B (ddof=0),
per coordinate instance, over rows whose selected raw Ca B-factor is usable and
finite. Raw B-factor summaries are descriptive QC only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"

MASTER_CSV = PROCESSED / "structures.csv"
RESIDUE_PARQUET = PROCESSED / "lbd_residue_bfactors_all.parquet"
LBD_REFERENCE_CSV = MANIFESTS / "nr_lbd_reference.csv"

NORMALIZED_PARQUET = ANALYSIS / "normalized_lbd_bfactors_all.parquet"
INSTANCE_SUMMARY_CSV = ANALYSIS / "instance_bfactor_summary.csv"
PDB_PROFILE_PARQUET = ANALYSIS / "pdb_position_bfactor_profiles.parquet"
RECEPTOR_PROFILE_PARQUET = ANALYSIS / "receptor_position_flexibility_profiles.parquet"
RECEPTOR_SUMMARY_CSV = ANALYSIS / "receptor_flexibility_summary.csv"

Z_HIGH = 1.5
SD_EPS = 1e-12
MIN_USABLE = 2

MASTER_ADDED_COLUMNS = [
    "representative_n_ca_usable",
    "representative_raw_b_mean_qc",
    "representative_raw_b_median_qc",
    "representative_raw_b_sd_qc",
    "representative_raw_b_iqr_qc",
    "representative_fraction_z_gt_1_5",
    "representative_fraction_z_lt_minus_1_5",
    "representative_normalization_status",
]

CANDIDATE_KEY = ["uniprot_id", "pdb_id", "polymer_entity_id"]


# --------------------------------------------------------------------------
# Part 3: normalized residue table
# --------------------------------------------------------------------------

def normalize_residues(res: pd.DataFrame) -> pd.DataFrame:
    b = res["selected_ca_b_iso_or_equiv"]
    usable = res["ca_bfactor_usable"].astype(bool) & b.notna() & np.isfinite(b.astype(float))
    # The two definitions of "usable" must agree exactly (Stage 3B.1 semantics).
    assert (usable == res["ca_bfactor_usable"].astype(bool)).all(), "ca_bfactor_usable disagrees with selected B"
    assert (usable == b.notna()).all(), "selected B present for a non-usable row or vice versa"

    out = res.copy()
    out["_b_usable"] = b.where(usable)
    g = out.groupby("instance_id", sort=False)["_b_usable"]
    stats = pd.DataFrame(
        {
            "instance_n_usable_b": g.count(),
            "instance_raw_b_mean": g.mean(),
            "instance_raw_b_median": g.median(),
            "instance_raw_b_sd": g.std(ddof=0),
            "instance_raw_b_q25": g.quantile(0.25),
            "instance_raw_b_q75": g.quantile(0.75),
        }
    )
    stats["instance_raw_b_iqr"] = stats["instance_raw_b_q75"] - stats["instance_raw_b_q25"]
    status = np.where(
        stats["instance_n_usable_b"] < MIN_USABLE,
        "NOT_NORMALIZABLE_FEWER_THAN_2_USABLE_B",
        np.where(~(stats["instance_raw_b_sd"] > SD_EPS), "NOT_NORMALIZABLE_ZERO_SD", "NORMALIZED"),
    )
    stats["normalization_status"] = status
    # Instances with <2 usable values: withhold dispersion statistics that are undefined/degenerate.
    for c in ["instance_raw_b_sd", "instance_raw_b_iqr"]:
        stats.loc[stats["instance_n_usable_b"] < MIN_USABLE, c] = np.nan

    out = out.join(stats, on="instance_id")
    ok = out["normalization_status"].eq("NORMALIZED") & out["_b_usable"].notna()
    out["b_factor_z"] = np.where(ok, (out["_b_usable"] - out["instance_raw_b_mean"]) / out["instance_raw_b_sd"], np.nan)
    out["row_normalization_status"] = np.where(
        out["b_factor_z"].notna(), "Z_COMPUTED",
        np.where(out["_b_usable"].isna(), "B_MISSING", "INSTANCE_NOT_NORMALIZABLE"),
    )
    out = out.drop(columns=["_b_usable"])

    # ---- assertions (Part 3) ----
    assert len(out) == len(res) == 771542, "row count changed"
    for c in res.columns:
        a, bcol = res[c], out[c]
        assert a.equals(bcol), f"existing column altered: {c}"
    assert out.loc[b.isna(), "b_factor_z"].isna().all(), "a missing raw B acquired a z value"
    fin = out["b_factor_z"].notna()
    assert np.isfinite(out.loc[fin, "b_factor_z"]).all() and fin.sum() > 0
    assert (b.notna() & np.isfinite(b.astype(float)))[fin].all(), "finite z without finite usable raw B"
    assert not np.isinf(out["b_factor_z"].astype(float)).any()
    zg = out.loc[fin].groupby("instance_id")["b_factor_z"]
    m, s = zg.mean(), zg.std(ddof=0)
    assert m.abs().max() < 1e-9, f"instance mean(z) not ~0: {m.abs().max()}"
    assert (s - 1).abs().max() < 1e-9, f"instance SD(z) not ~1: {(s - 1).abs().max()}"
    return out


# --------------------------------------------------------------------------
# Part 4: instance summary
# --------------------------------------------------------------------------

def build_instance_summary(norm: pd.DataFrame, master: pd.DataFrame) -> pd.DataFrame:
    ident = ["uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id", "entity_id",
             "instance_id", "label_asym_id", "auth_asym_id", "model_num", "is_representative_instance"]
    first = norm.groupby("instance_id", sort=False)[ident[:-1] + ["is_representative_instance"]].first().reset_index(drop=True)
    assert first["instance_id"].is_unique
    usable_flag = (norm["row_normalization_status"] != "B_MISSING").rename("_u")
    g = norm.assign(_u=usable_flag).groupby("instance_id", sort=False)
    counts = pd.DataFrame({"n_lbd_rows": g.size(), "n_ca_usable": g["_u"].sum().astype(int)})
    counts["n_ca_missing"] = counts["n_lbd_rows"] - counts["n_ca_usable"]
    counts["usable_fraction"] = counts["n_ca_usable"] / counts["n_lbd_rows"]

    inst = first.merge(counts, left_on="instance_id", right_index=True)
    raw = norm.groupby("instance_id", sort=False)[[
        "instance_raw_b_mean", "instance_raw_b_median", "instance_raw_b_sd",
        "instance_raw_b_q25", "instance_raw_b_q75", "instance_raw_b_iqr", "normalization_status",
    ]].first()
    raw.columns = ["raw_b_mean", "raw_b_median", "raw_b_sd", "raw_b_q25", "raw_b_q75", "raw_b_iqr", "normalization_status"]
    inst = inst.merge(raw, left_on="instance_id", right_index=True)

    zs = norm.dropna(subset=["b_factor_z"]).groupby("instance_id", sort=False)["b_factor_z"]
    zstat = pd.DataFrame(
        {
            "z_mean": zs.mean(), "z_sd": zs.std(ddof=0), "z_median": zs.median(),
            "z_q25": zs.quantile(0.25), "z_q75": zs.quantile(0.75),
            "fraction_z_gt_1_5": zs.apply(lambda s: (s > Z_HIGH).mean()),
            "fraction_z_lt_minus_1_5": zs.apply(lambda s: (s < -Z_HIGH).mean()),
        }
    )
    zstat["z_iqr"] = zstat["z_q75"] - zstat["z_q25"]
    inst = inst.merge(zstat, left_on="instance_id", right_index=True, how="left")

    meta_cols = CANDIDATE_KEY + [
        "group", "standardized_lbd_length", "resolution", "r_free", "positive_occupancy_lbd_coverage",
        "ligand_state", "ligand_annotation_status", "previous_primary_qc_member", "identity_class", "fusion_or_chimera",
    ]
    inst = inst.merge(master[meta_cols], on=CANDIDATE_KEY, how="left", validate="many_to_one")
    assert inst["group"].notna().all(), "instance without a master candidate row"
    assert (inst["standardized_lbd_length"].astype(int) == inst["n_lbd_rows"]).all()

    order = [
        "uniprot_id", "nr_code", "common_name", "group", "pdb_id", "polymer_entity_id", "entity_id",
        "instance_id", "label_asym_id", "auth_asym_id", "model_num", "is_representative_instance",
        "standardized_lbd_length", "n_lbd_rows", "n_ca_usable", "n_ca_missing", "usable_fraction",
        "raw_b_mean", "raw_b_median", "raw_b_sd", "raw_b_q25", "raw_b_q75", "raw_b_iqr",
        "z_mean", "z_sd", "z_median", "z_q25", "z_q75", "z_iqr",
        "fraction_z_gt_1_5", "fraction_z_lt_minus_1_5", "normalization_status",
        "resolution", "r_free", "positive_occupancy_lbd_coverage", "ligand_state", "ligand_annotation_status",
        "previous_primary_qc_member", "identity_class", "fusion_or_chimera",
    ]
    inst = inst[order].sort_values(["uniprot_id", "pdb_id", "polymer_entity_id", "instance_id", "model_num"]).reset_index(drop=True)
    assert len(inst) == norm["instance_id"].nunique()
    return inst


# --------------------------------------------------------------------------
# Part 5: master descriptive columns (membership unchanged)
# --------------------------------------------------------------------------

def add_master_columns(inst: pd.DataFrame) -> None:
    master = pd.read_csv(MASTER_CSV, dtype=str, keep_default_na=False)
    master = master.drop(columns=[c for c in MASTER_ADDED_COLUMNS if c in master.columns])
    n0, cols0 = len(master), list(master.columns)
    rep = inst.set_index("instance_id")
    mapping = {
        "representative_n_ca_usable": "n_ca_usable",
        "representative_raw_b_mean_qc": "raw_b_mean",
        "representative_raw_b_median_qc": "raw_b_median",
        "representative_raw_b_sd_qc": "raw_b_sd",
        "representative_raw_b_iqr_qc": "raw_b_iqr",
        "representative_fraction_z_gt_1_5": "fraction_z_gt_1_5",
        "representative_fraction_z_lt_minus_1_5": "fraction_z_lt_minus_1_5",
        "representative_normalization_status": "normalization_status",
    }
    ids = master["representative_instance_id"]
    new = {}
    for new_col, src in mapping.items():
        vals = ids.map(rep[src])  # blank id -> NaN
        if new_col == "representative_n_ca_usable":
            vals = vals.astype("Int64")
        new[new_col] = vals
    new_df = pd.DataFrame(new)
    coord = master["has_stage3b_coordinate_data"].eq("True")
    assert new_df.loc[coord, "representative_normalization_status"].notna().all()
    assert new_df.loc[~coord].isna().all().all(), "non-coordinate candidates must stay missing"
    out = pd.concat([master, new_df], axis=1)
    assert len(out) == n0 == 2072 and list(out.columns[: len(cols0)]) == cols0
    out.to_csv(MASTER_CSV, index=False)


# --------------------------------------------------------------------------
# Part 6: PDB-level observations (no pseudo-replication)
# --------------------------------------------------------------------------

def build_pdb_profiles(norm: pd.DataFrame, inst: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    reps = inst[inst["is_representative_instance"] & inst["normalization_status"].eq("NORMALIZED")]
    # PDB-level attributes across all normalizable representative candidates of receptor x PDB.
    rows = []
    for (uid, pdb), g in reps.groupby(["uniprot_id", "pdb_id"], sort=True):
        states = sorted(set(g["ligand_state"]))
        if len(states) == 1:
            state, conflict = states[0], False
        else:
            state, conflict = "AMBIGUOUS", True
        assert state in ("APO", "HOLO", "AMBIGUOUS"), state
        res_vals = g["resolution"].dropna().unique()
        assert len(res_vals) <= 1, f"{uid}/{pdb}: candidates disagree on resolution"
        rows.append(
            {
                "uniprot_id": uid, "pdb_id": pdb, "ligand_state": state, "ligand_state_conflict": conflict,
                "contributing_ligand_states_json": json.dumps(states),
                "n_representative_candidates_normalizable": len(g),
                "resolution": float(res_vals[0]) if len(res_vals) else np.nan,
                "any_previous_primary_qc_member": bool(g["previous_primary_qc_member"].any()),
                "max_positive_occupancy_lbd_coverage": float(g["positive_occupancy_lbd_coverage"].max()),
            }
        )
    pdb_attr = pd.DataFrame(rows)

    rep_ids = set(reps["instance_id"])
    obs = norm[norm["instance_id"].isin(rep_ids) & norm["b_factor_z"].notna()]
    grp = obs.groupby(["uniprot_id", "pdb_id", "canonical_uniprot_position"], sort=True)
    prof = grp.agg(
        b_factor_z=("b_factor_z", "median"),
        n_contributing_representative_candidates=("instance_id", "nunique"),
    ).reset_index()
    ident = norm.groupby("uniprot_id")[["nr_code", "common_name"]].first()
    grpmap = inst.groupby("uniprot_id")["group"].first()
    names = norm.drop_duplicates(["uniprot_id", "canonical_uniprot_position"]).set_index(
        ["uniprot_id", "canonical_uniprot_position"])["canonical_residue_name"]
    prof["nr_code"] = prof["uniprot_id"].map(ident["nr_code"])
    prof["common_name"] = prof["uniprot_id"].map(ident["common_name"])
    prof["group"] = prof["uniprot_id"].map(grpmap)
    prof["canonical_residue_name"] = [names[(u, p)] for u, p in zip(prof["uniprot_id"], prof["canonical_uniprot_position"])]
    prof = prof.merge(pdb_attr, on=["uniprot_id", "pdb_id"], how="left", validate="many_to_one")
    assert prof["ligand_state"].notna().all()
    assert not prof.duplicated(["uniprot_id", "pdb_id", "canonical_uniprot_position"]).any()
    cols = ["uniprot_id", "nr_code", "common_name", "group", "pdb_id", "canonical_uniprot_position",
            "canonical_residue_name", "b_factor_z", "n_contributing_representative_candidates",
            "ligand_state", "ligand_state_conflict", "contributing_ligand_states_json",
            "resolution", "any_previous_primary_qc_member", "max_positive_occupancy_lbd_coverage"]
    return prof[cols], pdb_attr


# --------------------------------------------------------------------------
# Part 7: receptor x position profile
# --------------------------------------------------------------------------

def _stats(df: pd.DataFrame, prefix: str, keys: list[str]) -> pd.DataFrame:
    g = df.groupby(keys)["b_factor_z"]
    return pd.DataFrame(
        {
            f"{prefix}n": g.count(),
            f"{prefix}mean_z": g.mean(),
            f"{prefix}median_z": g.median(),
            f"{prefix}sd_z": g.std(ddof=1),  # sample SD across PDB entries; NaN when n < 2
        }
    )


def build_receptor_profile(pdb_prof: pd.DataFrame, pdb_attr: pd.DataFrame, lbd_ref: pd.DataFrame,
                           norm: pd.DataFrame) -> pd.DataFrame:
    keys = ["uniprot_id", "canonical_uniprot_position"]
    n_total = pdb_attr.groupby("uniprot_id")["pdb_id"].nunique()
    g = pdb_prof.groupby(keys)["b_factor_z"]
    allst = pd.DataFrame(
        {
            "n_pdb_entries_observed_at_position": g.count(), "mean_z": g.mean(), "median_z": g.median(),
            "sd_z": g.std(ddof=1), "q25_z": g.quantile(0.25), "q75_z": g.quantile(0.75),
            "min_z": g.min(), "max_z": g.max(),
        }
    )
    allst["iqr_z"] = allst["q75_z"] - allst["q25_z"]
    apo = _stats(pdb_prof[pdb_prof["ligand_state"] == "APO"], "apo_", keys)
    holo = _stats(pdb_prof[pdb_prof["ligand_state"] == "HOLO"], "holo_", keys)

    names = norm.drop_duplicates(keys).set_index(keys)["canonical_residue_name"]
    rows = []
    for _, r in lbd_ref.iterrows():
        uid = r["uniprot_id"]
        if uid not in n_total.index:
            continue
        for pos in range(int(r["lbd_start"]), int(r["lbd_end"]) + 1):
            rows.append((uid, r["nr_code"], r["common_name"], r["group"], pos))
    grid = pd.DataFrame(rows, columns=["uniprot_id", "nr_code", "common_name", "group", "canonical_uniprot_position"])
    grid["canonical_residue_name"] = [names.get((u, p), "") for u, p in zip(grid["uniprot_id"], grid["canonical_uniprot_position"])]
    grid["n_pdb_entries_total_for_receptor"] = grid["uniprot_id"].map(n_total).astype(int)
    out = grid.join(allst, on=keys).join(apo, on=keys).join(holo, on=keys)
    for c in ["n_pdb_entries_observed_at_position", "apo_n", "holo_n"]:
        out[c] = out[c].fillna(0).astype(int)
    out["observation_fraction"] = out["n_pdb_entries_observed_at_position"] / out["n_pdb_entries_total_for_receptor"]
    both = (out["apo_n"] > 0) & (out["holo_n"] > 0)
    out["apo_holo_mean_delta"] = np.where(both, out["holo_mean_z"] - out["apo_mean_z"], np.nan)
    out["apo_holo_median_delta"] = np.where(both, out["holo_median_z"] - out["apo_median_z"], np.nan)
    cols = ["uniprot_id", "nr_code", "common_name", "group", "canonical_uniprot_position", "canonical_residue_name",
            "n_pdb_entries_total_for_receptor", "n_pdb_entries_observed_at_position", "observation_fraction",
            "mean_z", "median_z", "sd_z", "q25_z", "q75_z", "iqr_z", "min_z", "max_z",
            "apo_n", "apo_mean_z", "apo_median_z", "apo_sd_z",
            "holo_n", "holo_mean_z", "holo_median_z", "holo_sd_z",
            "apo_holo_mean_delta", "apo_holo_median_delta"]
    out = out[cols]
    assert not out.duplicated(keys).any()
    return out


# --------------------------------------------------------------------------
# Part 8: receptor-level summary (all 48 receptors)
# --------------------------------------------------------------------------

def build_receptor_summary(master: pd.DataFrame, inst: pd.DataFrame, pdb_attr: pd.DataFrame,
                           rec_prof: pd.DataFrame, lbd_ref: pd.DataFrame) -> pd.DataFrame:
    coord = master["has_stage3b_coordinate_data"].astype(str).eq("True")
    reps = inst[inst["is_representative_instance"] & inst["normalization_status"].eq("NORMALIZED")]
    rows = []
    for _, r in lbd_ref.iterrows():
        uid = r["uniprot_id"]
        m = master[master["uniprot_id"] == uid]
        mc = m[coord[m.index]]
        pa = pdb_attr[pdb_attr["uniprot_id"] == uid]
        rp = rec_prof[(rec_prof["uniprot_id"] == uid) & (rec_prof["n_pdb_entries_observed_at_position"] >= 1)]
        rr = reps[reps["uniprot_id"] == uid]
        n_norm = pa["pdb_id"].nunique()
        if len(m) == 0:
            status = "NO_STRUCTURES_IN_MASTER"
        elif len(mc) == 0:
            status = "NO_COORDINATE_DATA"
        elif n_norm == 0:
            status = "NO_NORMALIZABLE_INSTANCE"
        elif n_norm < 3:
            status = "PROFILE_AVAILABLE_LOW_N"
        else:
            status = "PROFILE_AVAILABLE"
        has = n_norm > 0
        med = rp["median_z"]
        rows.append(
            {
                "uniprot_id": uid, "nr_code": r["nr_code"], "common_name": r["common_name"], "group": r["group"],
                "n_master_candidates": len(m), "n_coordinate_candidates": len(mc),
                "n_unique_pdbs": m["pdb_id"].nunique(), "n_unique_coordinate_pdbs": mc["pdb_id"].nunique(),
                "n_normalizable_pdbs": n_norm,
                "median_lbd_observation_fraction": float(rr["usable_fraction"].median()) if has else np.nan,
                "n_positions_with_observation": int(len(rp)) if has else 0,
                "median_position_median_z": float(med.median()) if has and len(med) else np.nan,
                "iqr_of_position_median_z": float(med.quantile(0.75) - med.quantile(0.25)) if has and len(med) else np.nan,
                "fraction_positions_median_z_gt_1_0": float((med > 1.0).mean()) if has and len(med) else np.nan,
                "fraction_positions_median_z_lt_minus_1_0": float((med < -1.0).mean()) if has and len(med) else np.nan,
                "n_positions_ge3_pdbs": int((rp["n_pdb_entries_observed_at_position"] >= 3).sum()) if has else 0,
                "n_positions_ge5_pdbs": int((rp["n_pdb_entries_observed_at_position"] >= 5).sum()) if has else 0,
                "n_positions_ge10_pdbs": int((rp["n_pdb_entries_observed_at_position"] >= 10).sum()) if has else 0,
                "apo_pdb_count": int((pa["ligand_state"] == "APO").sum()),
                "holo_pdb_count": int((pa["ligand_state"] == "HOLO").sum()),
                "ambiguous_pdb_count": int((pa["ligand_state"] == "AMBIGUOUS").sum()),
                "profile_status": status,
            }
        )
    out = pd.DataFrame(rows)
    assert len(out) == 48 and out["uniprot_id"].is_unique
    return out


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(MASTER_CSV, dtype={"pdb_id": str})
    assert len(master) == 2072 and int(master["has_stage3b_coordinate_data"].sum()) == 1840
    res = pd.read_parquet(RESIDUE_PARQUET)
    assert len(res) == 771542
    lbd_ref = pd.read_csv(LBD_REFERENCE_CSV)

    print("Part 3: normalizing residues within instance...")
    norm = normalize_residues(res)
    norm.to_parquet(NORMALIZED_PARQUET, index=False)
    print(f"  wrote {NORMALIZED_PARQUET.relative_to(PROJECT_ROOT)} ({len(norm)} rows)")
    print(norm.groupby("instance_id")["normalization_status"].first().value_counts().to_string())

    print("Part 4: instance summary...")
    inst = build_instance_summary(norm, master)
    inst.to_csv(INSTANCE_SUMMARY_CSV, index=False)
    print(f"  wrote {INSTANCE_SUMMARY_CSV.relative_to(PROJECT_ROOT)} ({len(inst)} rows)")

    print("Part 5: descriptive representative columns on master (2072 rows retained)...")
    add_master_columns(inst)

    print("Part 6: PDB-level position observations...")
    pdb_prof, pdb_attr = build_pdb_profiles(norm, inst)
    pdb_prof.to_parquet(PDB_PROFILE_PARQUET, index=False)
    print(f"  wrote {PDB_PROFILE_PARQUET.relative_to(PROJECT_ROOT)} ({len(pdb_prof)} rows, {pdb_attr['pdb_id'].nunique()} PDBs)")

    print("Part 7: receptor x position profiles...")
    rec_prof = build_receptor_profile(pdb_prof, pdb_attr, lbd_ref, norm)
    rec_prof.to_parquet(RECEPTOR_PROFILE_PARQUET, index=False)
    print(f"  wrote {RECEPTOR_PROFILE_PARQUET.relative_to(PROJECT_ROOT)} ({len(rec_prof)} rows, {rec_prof['uniprot_id'].nunique()} receptors)")

    print("Part 8: receptor summary (48 receptors)...")
    master2 = pd.read_csv(MASTER_CSV, dtype={"pdb_id": str})
    summ = build_receptor_summary(master2, inst, pdb_attr, rec_prof, lbd_ref)
    summ.to_csv(RECEPTOR_SUMMARY_CSV, index=False)
    print(f"  wrote {RECEPTOR_SUMMARY_CSV.relative_to(PROJECT_ROOT)} ({len(summ)} rows)")
    print(summ["profile_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
