"""
Parts 3 and 5 of the functional-site specificity/validation stage:
structural-class-matched background, and the matched site-vs-background
paired statistics that use it.

Matched control hierarchy (fixed, tried in this exact order):
  1. same canonical structural segment, non-site residues of that receptor
  2. same structural class (HELIX/SHEET/LOOP/LOOP_UNCERTAIN), non-site
     residues of that receptor, any segment
  3. the original whole-LBD background (already computed, Part 8 of the
     prior stage) -- used only as a documented fallback, never silently
     conflated with tiers 1/2.
A matched comparison is only used when it has >= MIN_CONTROL_RESIDUES
control positions; otherwise the tier is skipped and the next one tried.

Only RESIDUE-level site definitions are matched here (Part 4: SEGMENT-
level proxies are never given a same-segment matched background --
comparing a segment against itself minus its own proxy "site" residues
would be circular, since the proxy IS the segment).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

MIN_CONTROL_RESIDUES = 5
MIN_PAIRS = 5


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


def rank_biserial_effect_size(diffs: np.ndarray) -> float:
    nz = diffs[diffs != 0]
    if len(nz) == 0:
        return 0.0
    ranks = stats.rankdata(np.abs(nz))
    w_pos = ranks[nz > 0].sum()
    w_neg = ranks[nz < 0].sum()
    return float((w_pos - w_neg) / (w_pos + w_neg))


def bootstrap_ci(diffs: np.ndarray, seed=0, n_boot=2000):
    if len(diffs) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.array([np.median(rng.choice(diffs, size=len(diffs), replace=True)) for _ in range(n_boot)])
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def build_matched_controls(definitions: pd.DataFrame, seg_map: pd.DataFrame) -> pd.DataFrame:
    """For every RESIDUE-level (uniprot_id, functional_site), return the matched
    control position set and which tier was used."""
    pos_seg = seg_map.set_index(["uniprot_id", "canonical_uniprot_position"])[["structural_segment", "structural_class"]]
    receptor_all_positions = seg_map[seg_map["structural_segment"] != ""].groupby("uniprot_id")[
        ["canonical_uniprot_position", "structural_segment", "structural_class"]
    ].apply(lambda g: g).reset_index(drop=True)

    residue_defs = definitions[definitions["definition_level"] == "RESIDUE"]
    rows = []
    for (uid, site), g in residue_defs.groupby(["uniprot_id", "functional_site"]):
        site_positions = set(g["canonical_uniprot_position"])
        all_rec_pos = seg_map[(seg_map["uniprot_id"] == uid) & (seg_map["structural_segment"] != "")]
        if len(all_rec_pos) == 0:
            continue

        # what segment(s)/class(es) do the site residues themselves sit in?
        site_annotated = all_rec_pos[all_rec_pos["canonical_uniprot_position"].isin(site_positions)]
        if len(site_annotated) == 0:
            continue  # site residues have no structural annotation at all -- cannot match

        site_segments = set(site_annotated["structural_segment"])
        site_classes = set(site_annotated["structural_class"])

        # tier 1: same segment(s), non-site residues
        tier1 = all_rec_pos[all_rec_pos["structural_segment"].isin(site_segments) &
                             ~all_rec_pos["canonical_uniprot_position"].isin(site_positions)]
        # tier 2: same class(es), non-site residues, any segment
        tier2 = all_rec_pos[all_rec_pos["structural_class"].isin(site_classes) &
                             ~all_rec_pos["canonical_uniprot_position"].isin(site_positions)]

        if len(tier1) >= MIN_CONTROL_RESIDUES:
            chosen_tier, control_positions = "1_same_segment", set(tier1["canonical_uniprot_position"])
        elif len(tier2) >= MIN_CONTROL_RESIDUES:
            chosen_tier, control_positions = "2_same_class", set(tier2["canonical_uniprot_position"])
        else:
            chosen_tier, control_positions = "3_whole_lbd_fallback", None  # handled downstream via existing background

        rows.append({
            "uniprot_id": uid, "functional_site": site,
            "site_segments": ",".join(sorted(site_segments)), "site_classes": ",".join(sorted(site_classes)),
            "matched_tier_used": chosen_tier,
            "n_control_residues": len(control_positions) if control_positions else None,
            "control_positions": control_positions,
        })
    return pd.DataFrame(rows)


def main():
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")

    matched_controls = build_matched_controls(definitions, seg_map)
    tier_counts = matched_controls["matched_tier_used"].value_counts()
    print("Matched-control tier usage (RESIDUE-level sites only):")
    print(tier_counts.to_string())

    control_out = TABLES / "functional_site_matched_control_tiers.csv"
    matched_controls.drop(columns=["control_positions"]).to_csv(control_out, index=False)
    print(f"\nWrote {control_out.relative_to(PROJECT_ROOT)}")

    # ---- PDB-level matched-background scores ----
    site_position_sets = {
        (uid, site): set(g["canonical_uniprot_position"])
        for (uid, site), g in definitions[definitions["definition_level"] == "RESIDUE"].groupby(["uniprot_id", "functional_site"])
    }
    control_lookup = {
        (r["uniprot_id"], r["functional_site"]): r["control_positions"]
        for _, r in matched_controls.iterrows() if r["control_positions"] is not None
    }
    # also need the original whole-LBD background for tier-3 fallback cases
    all_positions_by_receptor = {
        uid: set(g["canonical_uniprot_position"])
        for uid, g in seg_map[seg_map["structural_segment"] != ""].groupby("uniprot_id")
    }

    rows = []
    for (uid, pdb_id), g in positions.groupby(["uniprot_id", "pdb_id"]):
        nr_code = g["nr_code"].iloc[0]
        g = g.dropna(subset=["b_factor_z"])
        if len(g) == 0:
            continue
        observed = set(g["canonical_uniprot_position"])
        for (uid2, site), site_pos in site_position_sets.items():
            if uid2 != uid:
                continue
            observed_site = observed & site_pos
            if len(observed_site) == 0:
                continue

            key = (uid, site)
            tier_row = matched_controls[(matched_controls["uniprot_id"] == uid) & (matched_controls["functional_site"] == site)]
            if len(tier_row) == 0:
                continue
            tier = tier_row.iloc[0]["matched_tier_used"]
            if tier == "3_whole_lbd_fallback":
                control_pos = all_positions_by_receptor.get(uid, set()) - site_pos
            else:
                control_pos = control_lookup.get(key, set())

            observed_control = observed & control_pos
            if len(observed_control) < 3:
                continue

            site_vals = g.loc[g["canonical_uniprot_position"].isin(observed_site), "b_factor_z"]
            control_vals = g.loc[g["canonical_uniprot_position"].isin(observed_control), "b_factor_z"]

            rows.append({
                "uniprot_id": uid, "nr_code": nr_code, "pdb_id": pdb_id, "functional_site": site,
                "matched_tier_used": tier,
                "n_site_residues_observed": len(observed_site), "n_control_residues_observed": len(observed_control),
                "median_site_z": float(site_vals.median()), "median_matched_background_z": float(control_vals.median()),
                "paired_delta": float(site_vals.median() - control_vals.median()),
            })

    pdb_matched = pd.DataFrame(rows).sort_values(["uniprot_id", "functional_site", "pdb_id"]).reset_index(drop=True)
    out1 = ANALYSIS / "pdb_functional_site_matched_background.parquet"
    pdb_matched.to_parquet(out1, index=False)
    print(f"\nWrote {out1.relative_to(PROJECT_ROOT)} ({len(pdb_matched)} rows)")

    # ---- Wilcoxon signed-rank, paired by PDB entry, BH within + global ----
    test_rows = []
    for (uid, site), g in pdb_matched.groupby(["uniprot_id", "functional_site"]):
        if len(g) < MIN_PAIRS:
            continue
        nr_code = g["nr_code"].iloc[0]
        site_vals, bg_vals = g["median_site_z"].values, g["median_matched_background_z"].values
        diffs = site_vals - bg_vals
        if np.all(diffs == 0):
            continue
        stat, p = stats.wilcoxon(site_vals, bg_vals, alternative="two-sided", zero_method="wilcox")
        effect = rank_biserial_effect_size(diffs)
        ci_low, ci_high = bootstrap_ci(diffs)
        test_rows.append({
            "uniprot_id": uid, "nr_code": nr_code, "functional_site": site,
            "matched_tier_used": g["matched_tier_used"].iloc[0],
            "n_pairs": len(g), "median_site_z": float(np.median(site_vals)),
            "median_matched_background_z": float(np.median(bg_vals)),
            "median_paired_delta": float(np.median(diffs)),
            "bootstrap_ci_low": ci_low, "bootstrap_ci_high": ci_high,
            "wilcoxon_statistic": float(stat), "effect_size_rank_biserial": effect, "raw_p": float(p),
        })

    matched_tests = pd.DataFrame(test_rows)
    matched_tests["q_within_receptor"] = np.nan
    for uid, idx in matched_tests.groupby("uniprot_id").groups.items():
        matched_tests.loc[idx, "q_within_receptor"] = benjamini_hochberg(matched_tests.loc[idx, "raw_p"].values)
    matched_tests["significant_within_receptor"] = matched_tests["q_within_receptor"] <= 0.05
    matched_tests["q_global"] = benjamini_hochberg(matched_tests["raw_p"].values)
    matched_tests["significant_global"] = matched_tests["q_global"] <= 0.05
    matched_tests = matched_tests.sort_values(["uniprot_id", "raw_p"]).reset_index(drop=True)

    out2 = TABLES / "functional_site_matched_tests.csv"
    matched_tests.to_csv(out2, index=False)
    print(f"Wrote {out2.relative_to(PROJECT_ROOT)} ({len(matched_tests)} tests)")
    n_sig = matched_tests["significant_within_receptor"].sum()
    print(f"{n_sig}/{len(matched_tests)} significant within-receptor FDR (matched-background)")
    print(matched_tests["matched_tier_used"].value_counts().to_string())


if __name__ == "__main__":
    main()
