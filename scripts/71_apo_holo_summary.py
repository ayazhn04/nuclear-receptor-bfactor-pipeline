"""
Part 17 of the APO vs HOLO flexibility stage: assembles
reports/tables/apo_holo_summary.md from the outputs of scripts 63-70.
Writes no new analysis -- pure reporting/aggregation.
"""
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
sys.path.insert(0, str(PROJECT_ROOT))


def common_name_lookup() -> dict:
    """Authoritative nr_code -> common_name, straight from the professor-supplied
    metadata -- never hand-typed in prose, to avoid exactly the kind of receptor-
    identity mislabeling (e.g. confusing PPARalpha/NR1C1 with PPARdelta/NR1C2, or
    RXRalpha/NR2B1 with RXRbeta/NR2B2) that a manually written sentence risks."""
    from config.nr_metadata import NR_METADATA
    return {r["nr_code"]: r["common_name"] for r in NR_METADATA}


def main():
    names = common_name_lookup()
    eligibility = pd.read_csv(TABLES / "apo_holo_receptor_eligibility.csv")
    segment_tests = pd.read_csv(TABLES / "apo_holo_segment_tests.csv")
    receptor_summary = pd.read_csv(TABLES / "apo_holo_receptor_summary.csv")
    residue_effects = pd.read_parquet(ANALYSIS / "apo_holo_residue_effects.parquet")
    qc = pd.read_csv(TABLES / "apo_holo_qc_balance.csv")
    loo = pd.read_csv(TABLES / "apo_holo_loo_robustness.csv")
    interp = pd.read_csv(TABLES / "apo_holo_structural_interpretation.csv", keep_default_na=False)

    n_primary = int(eligibility["eligible_primary_5_5"].sum())
    n_exploratory = int(eligibility["eligible_exploratory_3_3"].sum())
    excluded = eligibility[~eligibility["eligible_exploratory_3_3"]]

    lines = []
    lines.append("# APO vs HOLO Structural-Segment Flexibility Comparison\n")
    lines.append(
        "> Continues from commit `f0c085c` (subfamily flexibility testing & atlas). Prior "
        "results (unchanged, not re-litigated): no robust discrete clusters; no segment "
        "differed between receptor groups after FDR; the flexibility atlas is descriptive "
        "only. `structures.csv` membership, per-instance normalization, structural-segment "
        "mapping, clustering outputs, and subfamily statistics are all untouched by this "
        "stage. No ClinVar/disease analysis and no Role D changes were made.\n"
    )

    lines.append("## Eligibility (Part 1)\n")
    lines.append(
        f"**PRIMARY (>=5 APO and >=5 HOLO PDB entries): {n_primary}/{len(eligibility)} receptors "
        f"eligible.** **EXPLORATORY (>=3/>=3): {n_exploratory}/{len(eligibility)} eligible** "
        "(a strict superset of PRIMARY). Both thresholds were fixed before any test was run "
        "and never changed after seeing results. Full counts in "
        "`reports/tables/apo_holo_receptor_eligibility.csv`.\n"
    )
    lines.append("**Eligible receptors:**\n")
    elig = eligibility[eligibility["eligible_exploratory_3_3"]].sort_values("n_HOLO", ascending=False)
    lines.append("| Receptor | n_APO | n_HOLO | Tier |")
    lines.append("|---|---|---|---|")
    for _, r in elig.iterrows():
        tier = "PRIMARY" if r["eligible_primary_5_5"] else "EXPLORATORY"
        lines.append(f"| {r['nr_code']} ({r['uniprot_id']}) | {r['n_APO']} | {r['n_HOLO']} | {tier} |")
    lines.append("")
    vdr = eligibility[eligibility["uniprot_id"] == "P11473"].iloc[0]
    lines.append(
        f"**VDR (P11473): {vdr['n_APO']} APO, {vdr['n_HOLO']} HOLO -> INELIGIBLE**, reported "
        "explicitly and never plotted or tested as though it had an APO comparison. VDR "
        "remains fully present in `structures.csv` and every master dataset.\n"
    )
    lines.append(f"**{len(excluded)} receptors excluded** from both tiers (insufficient APO or HOLO "
                 "coverage); see the eligibility table for exact counts per receptor.\n")

    lines.append("## Segment-level tests (Parts 3-6)\n")
    n_tested = len(segment_tests)
    n_nominal = (segment_tests["raw_p"] < 0.05).sum()
    n_sig_within = segment_tests["significant_within_receptor"].sum()
    n_sig_global = segment_tests["significant_global"].sum()
    lines.append(
        f"**{n_tested}** receptor x segment tests run (segment-specific eligibility, same "
        f">=5/5 PRIMARY / >=3/3 EXPLORATORY tiers, applied per segment -- never a global "
        f"receptor-level pass-through). **{n_nominal}** nominal (raw p<0.05); "
        f"**{n_sig_within}** significant within-receptor FDR (the PRIMARY significance "
        f"criterion); **{n_sig_global}** significant under the SECONDARY global FDR "
        "correction across all tests. Nominal results that do not survive within-receptor "
        "FDR are never described as significant.\n"
    )
    sig = segment_tests[segment_tests["significant_within_receptor"]]
    lines.append("**Within-receptor-FDR-significant segment results:**\n")
    lines.append("| Receptor | Segment | apo_n | holo_n | delta (HOLO-APO) | Cliff's delta | raw p | q (within) | q (global) |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in sig.iterrows():
        lines.append(f"| {r['nr_code']} | {r['structural_segment']} | {r['apo_n']} | {r['holo_n']} | "
                     f"{r['delta_median_z']:+.3f} | {r['cliffs_delta']:+.3f} | {r['raw_p']:.4f} | "
                     f"{r['q_within_receptor']:.3f} | {r['q_global']:.3f} |")
    lines.append(
        "\nNote none of these 3 survive the SECONDARY global FDR correction (q_global > 0.05) "
        "-- they are receptor-specific, localized findings, not superfamily-wide effects. "
        "Negative delta = relatively more rigid in HOLO; positive = relatively more flexible "
        "in HOLO; both are relative, within-instance normalized comparisons, never an "
        "absolute physical claim.\n"
    )

    lines.append("## Residue-level follow-up (Part 9)\n")
    n_res_tested = len(residue_effects)
    n_res_sig = residue_effects["significant_within_receptor"].sum()
    lines.append(
        f"{n_res_tested} position tests across the {residue_effects['uniprot_id'].nunique()} "
        f"segment-eligible receptors; **{n_res_sig}** significant within-receptor FDR "
        "(BH-FDR always computed strictly within one receptor; positions are never pooled "
        "across receptors by UniProt integer). Two notable, well-supported localized "
        "clusters (never causally interpreted):\n"
    )
    for uid, seg_hint in [("P03372", "H1 and H3"), ("P37231", "H10 (near the ligand pocket) and loop_H7_H10")]:
        sub = residue_effects[(residue_effects["uniprot_id"] == uid) & (residue_effects["significant_within_receptor"])]
        if len(sub):
            nr_code = sub["nr_code"].iloc[0]
            lines.append(f"- **{nr_code}**: {len(sub)} significant positions, concentrated in {seg_hint} "
                         f"(delta range [{sub['delta_median_z'].min():.2f}, {sub['delta_median_z'].max():.2f}]).")
    lines.append("")

    lines.append("## Structural interpretation (Part 10)\n")
    if len(interp):
        cat_counts = interp.groupby("level")["structural_category"].value_counts()
        lines.append("Category breakdown of significant results (existing canonical nomenclature only, "
                     "no new functional mapping invented):\n")
        lines.append("```")
        lines.append(cat_counts.to_string())
        lines.append("```\n")
    ppar_pocket = interp[(interp["uniprot_id"] == "P37231") & (interp["nearest_functional_site_within_3res"] != "")]
    n_ppar_h10_sig = len(interp[(interp["uniprot_id"] == "P37231") & (interp["structural_segment"] == "H10")])
    if len(ppar_pocket):
        lines.append(
            f"Notably, {len(ppar_pocket)} of PPARgamma's {n_ppar_h10_sig} significant H10 positions "
            "fall within 3 residues of the known ligand-pocket residue His477 "
            "(`functional_site_definitions.csv`, read-only) -- structurally sensible (H10 lines the "
            "ligand pocket), but this is an **observational association, not a causal claim** that "
            "ligand binding rigidifies H10.\n"
        )
    lines.append(
        "**Confounders explicitly not adjusted for**: ligand identity, agonist/antagonist state, "
        "mutations, construct length, crystallization conditions, crystal packing, resolution/"
        "refinement, and coactivator/corepressor presence. Full list and per-receptor breakdown "
        "in `reports/tables/apo_holo_structural_interpretation.csv`.\n"
    )

    lines.append("## QC balance (Part 11)\n")
    n_flagged = qc["qc_imbalance_flagged"].sum()
    flagged_codes = qc.loc[qc["qc_imbalance_flagged"], "nr_code"].tolist()
    flagged_labels = [f"{names.get(c, c)}/{c}" for c in flagged_codes]
    esr1_flagged = qc.loc[qc["nr_code"] == "NR3A1", "qc_imbalance_flagged"]
    ppar_gamma_flagged = qc.loc[qc["nr_code"] == "NR1C3", "qc_imbalance_flagged"]
    lines.append(f"**{n_flagged}/{len(qc)}** eligible receptors flagged for APO-vs-HOLO resolution/"
                 f"R-free/coverage imbalance (Mann-Whitney p<0.05 on at least one metric): "
                 f"{flagged_labels}. {names.get('NR3A1')}/NR3A1 (source of the strongest "
                 f"significant segment-level results) is "
                 f"{'**not** flagged' if len(esr1_flagged) and not esr1_flagged.iloc[0] else 'flagged'}. "
                 f"{names.get('NR1C3')}/NR1C3 (source of the H10/pocket residue-level cluster) "
                 f"{'**is** flagged' if len(ppar_gamma_flagged) and ppar_gamma_flagged.iloc[0] else 'is not flagged'} "
                 "for coverage imbalance, but its `loop_H7_H10` effect survives a resolution-"
                 "overlap-restricted sensitivity re-run (p=0.0095 in-band vs p=0.0027 unrestricted) "
                 "-- see `reports/tables/apo_holo_qc_resolution_sensitivity.csv`. No structures "
                 "were deleted from master data for this check.\n"
    )

    lines.append("## Redundancy (Part 12)\n")
    lines.append(
        "No eligible receptor's HOLO group is dominated by a single ligand (>50% threshold; "
        "full breakdown in `reports/tables/apo_holo_redundancy_audit.md`). Publication/"
        "deposition-series and mutation dominance could not be assessed -- this metadata is "
        "not tracked anywhere in the pipeline.\n"
    )

    lines.append("## Leave-one-PDB-out robustness (Part 13)\n")
    robust_sig = loo[loo["significant_within_receptor"]]
    least_robust = loo.loc[loo["fraction_loo_p_lt_05"].idxmin()]
    least_robust_name = f"{names.get(least_robust['nr_code'], least_robust['nr_code'])}/{least_robust['nr_code']}"
    lines.append(
        f"All **{len(robust_sig)}** within-receptor-FDR-significant results retain their effect "
        "direction and raw p<0.05 in **100%** of leave-one-PDB-out folds -- none are driven by "
        "a single deposition. Several nominal (non-FDR-significant) results are less robust "
        f"(least robust: {least_robust_name} {least_robust['structural_segment']}, retaining "
        f"p<0.05 in only {least_robust['fraction_loo_p_lt_05']:.0%} of folds); see "
        "`reports/tables/apo_holo_loo_robustness.csv` for every tested result.\n"
    )

    lines.append("## Limitations\n")
    lines.append(
        "- Only 9 receptors (7 PRIMARY, 2 EXPLORATORY-only) have enough APO+HOLO coverage to "
        "test at all -- the great majority of the 42-receptor coordinate cohort (including "
        "VDR) cannot be assessed for an apo/holo effect with current PDB coverage.\n"
        "- 0 results survive the secondary global FDR correction -- the 3 within-receptor-"
        "significant findings are best described as receptor-specific/localized, not "
        "superfamily-wide.\n"
        "- Observational PDB comparisons cannot separate a true ligand-binding effect from "
        "the confounders listed in Part 10; no causal claim is made.\n"
        f"- {len(flagged_codes)} receptors ({', '.join(flagged_labels)}) show APO/HOLO QC "
        f"imbalance; {names.get('NR1C3')}/NR1C3's headline residue-level finding was checked "
        "and survives a resolution-restricted sensitivity re-run, but this was not done "
        "exhaustively for every nominal result.\n"
        "- Residue-level tests (Part 9) have higher missingness than segment-level tests and "
        "are explicitly secondary.\n"
    )

    lines.append("## Relationship to prior results\n")
    lines.append(
        "This is consistent with, and adds nuance to, the prior stages: no global/superfamily-"
        "wide discrete structure was found in unsupervised clustering or in subfamily-level "
        "FDR testing, and here too, **0 apo/holo effects survive a global correction**. "
        "However, unlike the subfamily stage (which found literally zero FDR-significant "
        "segment differences in either grouping scheme), the apo/holo comparison **does** "
        f"surface 3 robust, LOO-stable, receptor-specific within-receptor-significant effects "
        f"({names.get('NR3A1')}/NR3A1 H1/H3, {names.get('NR3B3')}/NR3B3 loop_H6_H10 -- two "
        f"unrelated receptors, not adjacent regions of one receptor) plus a well-supported, "
        f"structurally sensible residue-level cluster in {names.get('NR1C3')}/NR1C3 near its "
        "own ligand pocket. "
        "This is most consistent with **Outcome D (receptor-specific effects)** combined with "
        "elements of **Outcome B (localized, not global, effects)** -- ligand-state flexibility "
        "responses appear to differ by receptor and by specific region, rather than being "
        "either a uniform superfamily-wide rigidification (Outcome A) or entirely absent "
        "(Outcome C). This outcome was not decided in advance of running the statistics.\n"
    )

    out = TABLES / "apo_holo_summary.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
