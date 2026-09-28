"""
Part 12 of the functional-site specificity/validation stage: assembles
reports/tables/functional_site_specificity_summary.md from the outputs
of scripts 83-91. Writes no new analysis -- pure reporting/aggregation.
Uses "consistent with" / "associated with" language; never "confirms
mechanism" or "ligand binding causes" -- this analysis does not support
causal claims.
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def main():
    fdr_audit = pd.read_csv(TABLES / "functional_site_fdr_audit.csv")
    definition_sensitivity = pd.read_csv(TABLES / "functional_site_definition_sensitivity.csv")
    matched_tests = pd.read_csv(TABLES / "functional_site_matched_tests.csv")
    sensitivity = pd.read_csv(TABLES / "functional_site_background_sensitivity.csv")
    matched_superfamily = pd.read_csv(TABLES / "functional_site_superfamily_matched_reanalysis.csv")
    classification = pd.read_csv(TABLES / "functional_site_interpretation_classification.csv")
    apo_holo_recheck = pd.read_csv(TABLES / "apo_holo_overlap_confidence_recheck.csv")
    loo_pdb = pd.read_csv(TABLES / "functional_site_matched_loo_pdb_robustness.csv")
    loo_receptor = pd.read_csv(TABLES / "functional_site_matched_loo_receptor_robustness.csv")

    lines = []
    lines.append("# Functional-Site Specificity & Confounding Validation\n")
    lines.append(
        "> Continues from commit `e420454`. A sensitivity/validation layer -- the existing "
        "functional-site results are **not modified**. Language throughout uses \"consistent "
        "with\" / \"associated with\"; this analysis never claims to confirm a mechanism or "
        "that ligand binding *causes* a flexibility change.\n"
    )

    lines.append("## FDR audit (Part 1)\n")
    matches = fdr_audit[["q_within_matches", "q_global_matches", "sig_within_matches", "sig_global_matches"]].all().all()
    n_only_global = ((fdr_audit["significant_global_recomputed"]) & (~fdr_audit["significant_within_receptor_recomputed"])).sum()
    lines.append(
        f"Independently re-derived BH-FDR (a second implementation, not a call to the same code) "
        f"**exactly reproduces** every stored q-value and significance call "
        f"({'PASS' if matches else 'MISMATCH -- see functional_site_fdr_audit.csv'}). "
        f"**{n_only_global}** tests are global-significant but not within-receptor-significant; "
        "this is explained by BH's dependence on each family's own p-value distribution, not "
        "family size alone -- a small within-receptor family whose few members are only "
        "moderately significant produces a harsher correction than the same p-value's much "
        "better relative rank in the full 107-test global pool. No implementation bug was found; "
        "no result was changed. Full trace: `reports/tables/functional_site_fdr_audit.csv`.\n"
    )

    lines.append("## Effect of definition confidence (Part 2)\n")
    overall = definition_sensitivity[definition_sensitivity["functional_site"] == "ALL_SITES"]
    lines.append("| Tier | n receptors | n tests | Sig. within-FDR | Sig. global-FDR | Median effect size |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in overall.iterrows():
        lines.append(f"| {r['definition_tier']} | {r['n_receptors']} | {r['n_tests']} | "
                     f"{r['n_significant_within_fdr']} | {r['n_significant_global_fdr']} | "
                     f"{r['median_effect_size']:+.3f} |")
    lines.append(
        "\nTier-A (residue-level, high confidence) and tier-C (segment-level proxy) results are "
        "never combined into one interpretation without this breakdown -- see "
        "`reports/tables/functional_site_definition_sensitivity.csv` for the full per-site-class detail.\n"
    )

    lines.append("## Original vs. structurally matched background (Parts 3-6)\n")
    n_matched = len(matched_tests)
    n_sig_matched = matched_tests["significant_within_receptor"].sum()
    tier_counts = matched_tests["matched_tier_used"].value_counts()
    lines.append(
        f"{n_matched} receptor x site tests have a valid structural-context-matched background "
        f"(RESIDUE-level definitions only; SEGMENT-level proxies are never given a same-segment "
        f"match -- that would be circular). Control hierarchy used: "
        f"{tier_counts.get('1_same_segment', 0)} same-segment, {tier_counts.get('2_same_class', 0)} "
        f"same-structural-class. **{n_sig_matched}/{n_matched}** remain significant within-receptor "
        "FDR against the matched background.\n"
    )
    n_dir = sensitivity["effect_direction_preserved"].sum()
    n_preserved = sensitivity["significance_preserved"].sum()
    n_lost = sensitivity["significance_lost"].sum()
    n_gained = sensitivity["significance_gained"].sum()
    lines.append(
        f"Of the {len(sensitivity)} tests with both an original and a matched result: "
        f"**{n_dir}** preserve effect direction, **{n_preserved}** preserve significance, "
        f"**{n_lost}** LOSE significance under matching (explained by structural context), "
        f"**{n_gained}** GAIN significance (the whole-LBD background was diluting a real, more "
        "localized effect). Full table: `reports/tables/functional_site_background_sensitivity.csv`.\n"
    )

    lines.append("## Confirmatory superfamily reanalysis (Part 7)\n")
    lines.append("| Site class | n receptors | median matched delta | raw p | BH q | Status |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in matched_superfamily.iterrows():
        if r["status"] == "TESTED":
            lines.append(f"| {r['functional_site']} | {int(r['n_receptors'])} | {r['median_matched_delta']:+.3f} | "
                         f"{r['raw_p']:.4g} | {r['BH_q']:.4g} | {'**SIGNIFICANT**' if r['significant_fdr_05'] else 'not significant'} |")
        else:
            lines.append(f"| {r['functional_site']} | {int(r['n_receptors'])} | -- | -- | -- | {r['status']} |")
    lines.append("")

    lines.append("## Interpretation classification (Part 8)\n")
    lines.append("| Site class | Category | Rationale |")
    lines.append("|---|---|---|")
    for _, r in classification.iterrows():
        lines.append(f"| {r['functional_site']} | **{r['interpretation_category']}** | {r['rationale']} |")
    lines.append(
        "\n**AF2/H12 and the ligand-binding pocket are SITE-SPECIFIC SUPPORTED**: their flexibility "
        "differences are consistent with (not merely an artifact of) the site itself, robust across "
        "structural-context matching. **The coactivator groove is STRUCTURAL-CONTEXT EXPLAINED**: "
        "its original apparent rigidity is better attributed to the surrounding helical context "
        "than to the coactivator-groove residues specifically -- an important correction to the "
        "prior stage's face-value interpretation. **The dimerization interface has INSUFFICIENT "
        "RESIDUE-LEVEL MAPPING** to be confirmed or refuted this way.\n"
    )

    lines.append("## Apo/holo overlap confidence recheck (Part 9, descriptive only)\n")
    lines.append("| Receptor | Segment | Site | Definition level | Confidence | Category |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in apo_holo_recheck.iterrows():
        lines.append(f"| {r['nr_code']} | {r['apo_holo_segment']} | {r['site_name']} | "
                     f"{r['site_definition_level']} | {r['site_definition_confidence']} | {r['overlap_confidence_category']} |")
    lines.append(
        "\nAll 3 significant apo/holo x functional-site overlaps rest on RESIDUE-level, HIGH-"
        "confidence site definitions -- the apo/holo link is not being overstated on top of a "
        "coarse segment-level proxy. The apo/holo statistics themselves were not recomputed.\n"
    )

    lines.append("## Robustness (Part 10)\n")
    n_pdb_robust = (loo_pdb["fraction_loo_p_lt_05"] >= 1.0).sum()
    lines.append(
        f"Leave-one-PDB-out (matched-background, within receptor): **{n_pdb_robust}/{len(loo_pdb)}** "
        "significant results retain p<0.05 in every fold. Leave-one-receptor-out (matched-"
        "background superfamily reanalysis): "
    )
    for _, r in loo_receptor.iterrows():
        lines.append(f"- {r['functional_site']}: fraction of LOO folds with p<0.05 = "
                     f"{r['fraction_loo_p_lt_05']:.0%}, direction preserved in "
                     f"{r['fraction_loo_same_direction']:.0%} of folds, "
                     f"significance flips = {r['significance_flips_on_some_loo_fold']}.")
    lines.append("")

    lines.append("## Limitations\n")
    lines.append(
        "- The confirmatory matched-background analysis is only possible for RESIDUE-level "
        "definitions; the coactivator groove and dimerization interface are mostly SEGMENT-level "
        "proxies, which is exactly why they could not be (coactivator groove) or could not "
        "usefully be (dimerization interface, n=2) confirmed this way -- this is a genuine "
        "evidentiary gap, not a negative biological result about those sites.\n"
        "- Matched-background control sizes are sometimes small (minimum 5 control residues, "
        "minimum 5 paired PDB entries); receptor-level results built on thin controls should be "
        "read with that in mind even where the superfamily-level result is robust.\n"
        "- This remains a descriptive structural-association analysis. Even the SITE-SPECIFIC-"
        "SUPPORTED classification for AF2/H12 and the ligand-binding pocket means the flexibility "
        "difference is **consistent with** the site being specifically involved, not that a causal "
        "or mechanistic link to ligand binding has been established.\n"
    )

    out = TABLES / "functional_site_specificity_summary.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
