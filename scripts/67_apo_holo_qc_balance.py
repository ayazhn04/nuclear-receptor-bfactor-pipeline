"""
Part 11 of the APO vs HOLO flexibility stage: audit whether APO and HOLO
groups differ systematically in experimental quality (resolution, R-free,
LBD coverage) for every eligible receptor. No automatic regression
adjustment is performed at this stage -- receptors with substantial QC
imbalance are flagged, and a sensitivity re-run restricted to an
overlapping "reasonable quality" resolution band is reported where enough
structures remain on both sides. Nothing is deleted from master data.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

IMBALANCE_P_THRESHOLD = 0.05  # flag if resolution or coverage differs at this level
MIN_N_FOR_SENSITIVITY = 3


def main():
    eligibility = pd.read_csv(TABLES / "apo_holo_receptor_eligibility.csv")
    eligible_receptors = eligibility[eligibility["eligible_exploratory_3_3"]]["uniprot_id"].tolist()

    structures = pd.read_csv(PROCESSED / "structures.csv", dtype={"pdb_id": str}, keep_default_na=False)
    coord = structures[structures["coordinate_data_status"] == "COORDINATE_DATA_AVAILABLE"]
    coord_pdbs = coord.groupby(["uniprot_id", "nr_code", "pdb_id"]).agg(
        ligand_state=("ligand_state", "first"),
        resolution=("resolution", "first"),
        r_free=("r_free", "first"),
        positive_occupancy_lbd_coverage=("positive_occupancy_lbd_coverage", "first"),
    ).reset_index()
    coord_pdbs = coord_pdbs[coord_pdbs["uniprot_id"].isin(eligible_receptors)]
    coord_pdbs = coord_pdbs[coord_pdbs["ligand_state"].isin(["APO", "HOLO"])]

    rows = []
    for uid, g in coord_pdbs.groupby("uniprot_id"):
        nr_code = g["nr_code"].iloc[0]
        apo = g[g["ligand_state"] == "APO"]
        holo = g[g["ligand_state"] == "HOLO"]

        def summarize(metric):
            a = pd.to_numeric(apo[metric], errors="coerce").dropna()
            h = pd.to_numeric(holo[metric], errors="coerce").dropna()
            if len(a) < 2 or len(h) < 2:
                return dict(apo_median=np.nan, holo_median=np.nan, raw_p=np.nan, n_apo=len(a), n_holo=len(h))
            _, p = stats.mannwhitneyu(a, h, alternative="two-sided")
            return dict(apo_median=float(a.median()), holo_median=float(h.median()), raw_p=float(p),
                        n_apo=len(a), n_holo=len(h))

        res = summarize("resolution")
        rfree = summarize("r_free")
        cov = summarize("positive_occupancy_lbd_coverage")

        imbalanced = any(
            not np.isnan(d["raw_p"]) and d["raw_p"] < IMBALANCE_P_THRESHOLD
            for d in (res, rfree, cov)
        )

        rows.append({
            "uniprot_id": uid, "nr_code": nr_code,
            "apo_resolution_median": res["apo_median"], "holo_resolution_median": res["holo_median"],
            "resolution_mw_p": res["raw_p"],
            "apo_r_free_median": rfree["apo_median"], "holo_r_free_median": rfree["holo_median"],
            "r_free_mw_p": rfree["raw_p"],
            "apo_coverage_median": cov["apo_median"], "holo_coverage_median": cov["holo_median"],
            "coverage_mw_p": cov["raw_p"],
            "qc_imbalance_flagged": imbalanced,
        })

    qc = pd.DataFrame(rows).sort_values("uniprot_id").reset_index(drop=True)
    out = TABLES / "apo_holo_qc_balance.csv"
    qc.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(qc)} receptors)")
    print(qc.to_string(index=False))

    n_flagged = qc["qc_imbalance_flagged"].sum()
    print(f"\n{n_flagged}/{len(qc)} receptors flagged for QC imbalance (resolution, R-free, or "
          f"coverage differ APO vs HOLO at p<{IMBALANCE_P_THRESHOLD})")

    # ---- sensitivity: restrict to an overlapping "reasonable quality" resolution band ----
    print("\n--- Resolution-restricted sensitivity re-run ---")
    segment_effects = pd.read_parquet(ANALYSIS / "apo_holo_segment_effects.parquet")
    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")
    features = features[features["uniprot_id"].isin(eligible_receptors)]
    features = features[features["ligand_state"].isin(["APO", "HOLO"])]
    res_by_pdb = coord_pdbs.set_index(["uniprot_id", "pdb_id"])["resolution"]

    sens_rows = []
    for uid in qc.loc[qc["qc_imbalance_flagged"], "uniprot_id"]:
        g = coord_pdbs[coord_pdbs["uniprot_id"] == uid]
        apo_res = pd.to_numeric(g.loc[g["ligand_state"] == "APO", "resolution"], errors="coerce").dropna()
        holo_res = pd.to_numeric(g.loc[g["ligand_state"] == "HOLO", "resolution"], errors="coerce").dropna()
        if len(apo_res) == 0 or len(holo_res) == 0:
            continue
        lo = max(apo_res.min(), holo_res.min())
        hi = min(apo_res.max(), holo_res.max())
        if lo >= hi:
            sens_rows.append({"uniprot_id": uid, "overlap_band": None, "status": "NO_RESOLUTION_OVERLAP"})
            continue
        overlap_pdbs = set(g.loc[pd.to_numeric(g["resolution"], errors="coerce").between(lo, hi), "pdb_id"])
        sub = features[(features["uniprot_id"] == uid) & (features["pdb_id"].isin(overlap_pdbs))]
        for seg, sg in sub.groupby("structural_segment"):
            apo_vals = sg.loc[sg["ligand_state"] == "APO", "median_z"].dropna()
            holo_vals = sg.loc[sg["ligand_state"] == "HOLO", "median_z"].dropna()
            if len(apo_vals) < MIN_N_FOR_SENSITIVITY or len(holo_vals) < MIN_N_FOR_SENSITIVITY:
                continue
            _, p = stats.mannwhitneyu(holo_vals, apo_vals, alternative="two-sided")
            sens_rows.append({
                "uniprot_id": uid, "structural_segment": seg, "overlap_band": f"[{lo:.2f}, {hi:.2f}]",
                "n_apo_in_band": len(apo_vals), "n_holo_in_band": len(holo_vals),
                "delta_median_z_in_band": float(holo_vals.median() - apo_vals.median()),
                "raw_p_in_band": float(p), "status": "TESTED_IN_OVERLAP_BAND",
            })
    sens = pd.DataFrame(sens_rows)
    sens_out = TABLES / "apo_holo_qc_resolution_sensitivity.csv"
    sens.to_csv(sens_out, index=False)
    print(f"Wrote {sens_out.relative_to(PROJECT_ROOT)} ({len(sens)} rows)")
    if len(sens):
        print(sens.to_string(index=False))
    print("\nNo structures were deleted from master data for this sensitivity check -- it only "
          "restricts which entries feed the comparison, computed on a copy.")


if __name__ == "__main__":
    main()
