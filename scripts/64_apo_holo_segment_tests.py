"""
Parts 3-8 of the APO vs HOLO flexibility stage.

Analysis unit: PDB ENTRY (never chain, never residue for this primary
test). Outcome: entry-level segment median_z from
data/analysis/pdb_segment_flexibility_features.parquet -- itself already
built from the approved chain/instance -> per-instance Z -> PDB-entry-
median hierarchy. Eligibility is receptor AND segment specific (Part 6):
a receptor passing the overall 5/5 (or 3/3) tier can still have a
specific segment with too few observed APO/HOLO entries to test -- that
segment/receptor pair is simply not tested, never imputed, never treated
as 0.

Delta convention (Part 3): delta = HOLO - APO. Negative delta = segment
relatively more RIGID in holo structures; positive = relatively more
FLEXIBLE in holo structures. These are relative, within-instance
normalized B-factor comparisons -- never an absolute physical claim.

Multiple testing (Part 5): BH-FDR applied twice -- once within each
receptor across its own tested segments (q_within_receptor, the PRIMARY
significance criterion), and once globally across every receptor x
segment test in the whole study (q_global, reported alongside, never
substituted for the within-receptor criterion or for raw p).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

PRIMARY_MIN_N = 5
EXPLORATORY_MIN_N = 3
BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 0
FDR_ALPHA = 0.05


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
    """Cliff's delta comparing a (HOLO) vs b (APO), consistent with the
    delta = HOLO - APO direction: positive delta means HOLO values tend to
    be greater than APO values."""
    n_greater = sum((ai > b).sum() for ai in a)
    n_less = sum((ai < b).sum() for ai in a)
    return (n_greater - n_less) / (len(a) * len(b))


def bootstrap_delta_median_ci(apo: np.ndarray, holo: np.ndarray, seed=BOOTSTRAP_SEED, n_boot=BOOTSTRAP_N):
    rng = np.random.default_rng(seed)
    deltas = np.empty(n_boot)
    for i in range(n_boot):
        a = rng.choice(apo, size=len(apo), replace=True)
        h = rng.choice(holo, size=len(holo), replace=True)
        deltas[i] = np.median(h) - np.median(a)
    return float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))


def main():
    eligibility = pd.read_csv(TABLES / "apo_holo_receptor_eligibility.csv")
    candidate_receptors = eligibility[eligibility["eligible_exploratory_3_3"]]["uniprot_id"].tolist()
    print(f"{len(candidate_receptors)} receptors pass at least the EXPLORATORY (>=3/3) tier: "
          f"{candidate_receptors}")

    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")
    features = features[features["uniprot_id"].isin(candidate_receptors)]
    features = features[features["ligand_state"].isin(["APO", "HOLO"])]  # AMBIGUOUS excluded, never assigned either side

    rows = []
    for (uid, nr_code, seg), g in features.groupby(["uniprot_id", "nr_code", "structural_segment"]):
        apo_vals = g.loc[g["ligand_state"] == "APO", "median_z"].dropna().values
        holo_vals = g.loc[g["ligand_state"] == "HOLO", "median_z"].dropna().values
        n_apo, n_holo = len(apo_vals), len(holo_vals)

        if n_apo >= PRIMARY_MIN_N and n_holo >= PRIMARY_MIN_N:
            tier = "PRIMARY"
        elif n_apo >= EXPLORATORY_MIN_N and n_holo >= EXPLORATORY_MIN_N:
            tier = "EXPLORATORY"
        else:
            continue  # insufficient segment-specific coverage -- not tested, not imputed

        apo_median, holo_median = float(np.median(apo_vals)), float(np.median(holo_vals))
        apo_mean, holo_mean = float(np.mean(apo_vals)), float(np.mean(holo_vals))
        delta_median = holo_median - apo_median
        delta_mean = holo_mean - apo_mean

        U, p = stats.mannwhitneyu(holo_vals, apo_vals, alternative="two-sided")
        delta = cliffs_delta(holo_vals, apo_vals)
        ci_low, ci_high = bootstrap_delta_median_ci(apo_vals, holo_vals)

        rows.append({
            "uniprot_id": uid, "nr_code": nr_code, "structural_segment": seg,
            "apo_n": n_apo, "holo_n": n_holo,
            "apo_median_z": apo_median, "holo_median_z": holo_median, "delta_median_z": delta_median,
            "apo_mean_z": apo_mean, "holo_mean_z": holo_mean, "delta_mean_z": delta_mean,
            "bootstrap_ci_low": ci_low, "bootstrap_ci_high": ci_high,
            "cliffs_delta": delta, "mw_u": float(U), "raw_p": float(p),
            "eligibility_tier": tier,
        })

    effects = pd.DataFrame(rows)
    print(f"\n{len(effects)} receptor x segment tests run "
          f"({(effects['eligibility_tier'] == 'PRIMARY').sum()} PRIMARY tier, "
          f"{(effects['eligibility_tier'] == 'EXPLORATORY').sum()} EXPLORATORY tier)")

    # Part 5: BH-FDR within each receptor
    effects["q_within_receptor"] = np.nan
    for uid, idx in effects.groupby("uniprot_id").groups.items():
        pvals = effects.loc[idx, "raw_p"].values
        effects.loc[idx, "q_within_receptor"] = benjamini_hochberg(pvals)
    effects["significant_within_receptor"] = effects["q_within_receptor"] <= FDR_ALPHA

    # Part 5: SECONDARY global BH-FDR across every receptor x segment test
    effects["q_global"] = benjamini_hochberg(effects["raw_p"].values)
    effects["significant_global"] = effects["q_global"] <= FDR_ALPHA

    effects = effects.sort_values(["uniprot_id", "raw_p"]).reset_index(drop=True)

    parquet_out = ANALYSIS / "apo_holo_segment_effects.parquet"
    effects.to_parquet(parquet_out, index=False)
    print(f"\nWrote {parquet_out.relative_to(PROJECT_ROOT)}")

    csv_out = TABLES / "apo_holo_segment_tests.csv"
    effects.to_csv(csv_out, index=False)
    print(f"Wrote {csv_out.relative_to(PROJECT_ROOT)}")

    n_nominal = (effects["raw_p"] < 0.05).sum()
    n_sig_within = effects["significant_within_receptor"].sum()
    n_sig_global = effects["significant_global"].sum()
    print(f"\n{n_nominal} nominal (raw p<0.05); {n_sig_within} significant within-receptor FDR; "
          f"{n_sig_global} significant global FDR")
    if n_nominal > n_sig_within:
        print(f"NOTE: {n_nominal - n_sig_within} nominal result(s) do NOT survive within-receptor "
              "FDR correction and must not be described as significant.")

    if n_sig_within:
        print("\nWithin-receptor-FDR-significant results:")
        sig = effects[effects["significant_within_receptor"]]
        print(sig[["uniprot_id", "nr_code", "structural_segment", "apo_n", "holo_n",
                    "delta_median_z", "cliffs_delta", "raw_p", "q_within_receptor", "q_global"]].to_string(index=False))

    # ---- Part 8: receptor-level summary ----
    summary_rows = []
    for uid, g in effects.groupby("uniprot_id"):
        nr_code = g["nr_code"].iloc[0]
        n_tested = len(g)
        n_nominal_r = (g["raw_p"] < 0.05).sum()
        n_sig_within_r = g["significant_within_receptor"].sum()
        n_sig_global_r = g["significant_global"].sum()
        median_delta = float(g["delta_median_z"].median())
        frac_neg = float((g["delta_median_z"] < 0).mean())
        frac_pos = float((g["delta_median_z"] > 0).mean())
        strongest_rigid = g.loc[g["delta_median_z"].idxmin(), "structural_segment"] if n_tested else None
        strongest_flex = g.loc[g["delta_median_z"].idxmax(), "structural_segment"] if n_tested else None
        summary_rows.append({
            "uniprot_id": uid, "nr_code": nr_code,
            "n_tested_segments": n_tested, "n_nominal_p_lt_05": int(n_nominal_r),
            "n_significant_within_receptor_fdr": int(n_sig_within_r),
            "n_significant_global_fdr": int(n_sig_global_r),
            "median_delta_across_segments": median_delta,
            "fraction_segments_negative_delta": frac_neg,
            "fraction_segments_positive_delta": frac_pos,
            "strongest_rigidifying_segment": strongest_rigid,
            "strongest_flexibilizing_segment": strongest_flex,
        })
    receptor_summary = pd.DataFrame(summary_rows).sort_values("uniprot_id").reset_index(drop=True)
    receptor_out = TABLES / "apo_holo_receptor_summary.csv"
    receptor_summary.to_csv(receptor_out, index=False)
    print(f"\nWrote {receptor_out.relative_to(PROJECT_ROOT)} ({len(receptor_summary)} receptors)")
    print(receptor_summary.to_string(index=False))


if __name__ == "__main__":
    main()
