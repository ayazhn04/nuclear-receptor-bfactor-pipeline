"""
Part 9 of the APO vs HOLO flexibility stage: SECONDARY residue-level
analysis, restricted to receptors already eligible for the segment-level
analysis (Part 1/6, >=3 APO and >=3 HOLO PDB entries overall). Positions
are never pooled across receptors by UniProt integer -- every test and
every BH-FDR correction is performed strictly WITHIN one receptor.

Per-position eligibility mirrors the segment-level tiers (>=5/5 PRIMARY,
>=3/3 EXPLORATORY), again segment/position-specific: no imputation, no
zero-filling for missing residues.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

PRIMARY_MIN_N = 5
EXPLORATORY_MIN_N = 3


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    n_greater = sum((ai > b).sum() for ai in a)
    n_less = sum((ai < b).sum() for ai in a)
    return (n_greater - n_less) / (len(a) * len(b))


def main():
    eligibility = pd.read_csv(TABLES / "apo_holo_receptor_eligibility.csv")
    eligible_receptors = eligibility[eligibility["eligible_exploratory_3_3"]]["uniprot_id"].tolist()
    print(f"Residue-level analysis restricted to the {len(eligible_receptors)} segment-level-eligible "
          f"receptors: {eligible_receptors}")

    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")
    positions = positions[positions["uniprot_id"].isin(eligible_receptors)]
    positions = positions[positions["ligand_state"].isin(["APO", "HOLO"])]
    # a position observed inconsistently across contributing instances within one PDB
    # entry is flagged ligand_state_conflict -- exclude, never force it to a side
    if "ligand_state_conflict" in positions.columns:
        n_conflict = positions["ligand_state_conflict"].fillna(False).astype(bool).sum()
        positions = positions[~positions["ligand_state_conflict"].fillna(False).astype(bool)]
        print(f"Excluded {n_conflict} position observations with a ligand_state_conflict flag")

    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    seg_lookup = seg_map.set_index(["uniprot_id", "canonical_uniprot_position"])["structural_segment"]

    all_rows = []
    for uid, g_receptor in positions.groupby("uniprot_id"):
        nr_code = g_receptor["nr_code"].iloc[0]
        receptor_rows = []
        for pos, g in g_receptor.groupby("canonical_uniprot_position"):
            apo_vals = g.loc[g["ligand_state"] == "APO", "b_factor_z"].dropna().values
            holo_vals = g.loc[g["ligand_state"] == "HOLO", "b_factor_z"].dropna().values
            n_apo, n_holo = len(apo_vals), len(holo_vals)
            if n_apo >= PRIMARY_MIN_N and n_holo >= PRIMARY_MIN_N:
                tier = "PRIMARY"
            elif n_apo >= EXPLORATORY_MIN_N and n_holo >= EXPLORATORY_MIN_N:
                tier = "EXPLORATORY"
            else:
                continue

            apo_median, holo_median = float(np.median(apo_vals)), float(np.median(holo_vals))
            delta = holo_median - apo_median
            _, p = stats.mannwhitneyu(holo_vals, apo_vals, alternative="two-sided")
            cdelta = cliffs_delta(holo_vals, apo_vals)
            seg = seg_lookup.get((uid, pos), "")

            receptor_rows.append({
                "uniprot_id": uid, "nr_code": nr_code, "canonical_uniprot_position": int(pos),
                "apo_n": n_apo, "holo_n": n_holo,
                "apo_median_z": apo_median, "holo_median_z": holo_median, "delta_median_z": delta,
                "cliffs_delta": cdelta, "raw_p": float(p), "eligibility_tier": tier,
                "structural_segment": seg if seg else "",
            })

        if receptor_rows:
            r_df = pd.DataFrame(receptor_rows)
            r_df["q_within_receptor"] = benjamini_hochberg(r_df["raw_p"].values)
            all_rows.append(r_df)

    residue_effects = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame(columns=[
        "uniprot_id", "nr_code", "canonical_uniprot_position", "apo_n", "holo_n",
        "apo_median_z", "holo_median_z", "delta_median_z", "cliffs_delta", "raw_p",
        "q_within_receptor", "structural_segment",
    ])
    residue_effects["significant_within_receptor"] = residue_effects["q_within_receptor"] <= 0.05
    residue_effects = residue_effects.sort_values(["uniprot_id", "raw_p"]).reset_index(drop=True)

    out = ANALYSIS / "apo_holo_residue_effects.parquet"
    residue_effects.to_parquet(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(residue_effects)} position tests, "
          f"{residue_effects['uniprot_id'].nunique()} receptors)")

    n_nominal = (residue_effects["raw_p"] < 0.05).sum()
    n_sig = residue_effects["significant_within_receptor"].sum()
    print(f"{n_nominal} nominal (raw p<0.05); {n_sig} significant within-receptor FDR")
    if n_sig:
        print(residue_effects.loc[residue_effects["significant_within_receptor"],
              ["uniprot_id", "nr_code", "canonical_uniprot_position", "structural_segment",
               "delta_median_z", "cliffs_delta", "raw_p", "q_within_receptor"]].to_string(index=False))


if __name__ == "__main__":
    main()
