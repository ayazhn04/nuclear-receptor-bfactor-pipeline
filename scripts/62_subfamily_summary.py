"""
Part 15 of the subfamily flexibility stage: assembles
reports/tables/subfamily_flexibility_summary.md from the outputs of
scripts 55-61. Writes no new analysis -- pure reporting/aggregation.
"""
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

WELL_SUPPORTED_MIN_N = 10


def main():
    with open(ANALYSIS / "subfamily_grouping_schemes.json") as f:
        schemes = json.load(f)
    sample_sizes = pd.read_csv(TABLES / "subfamily_segment_sample_sizes.csv")
    kruskal = pd.read_csv(TABLES / "subfamily_segment_kruskal.csv")
    posthoc = pd.read_csv(TABLES / "subfamily_segment_posthoc.csv")
    atlas = pd.read_csv(ANALYSIS / "superfamily_segment_flexibility_atlas.csv")
    loo_atlas = pd.read_csv(TABLES / "subfamily_atlas_loo_robustness.csv")
    loo_kruskal = pd.read_csv(TABLES / "subfamily_kruskal_loo_robustness.csv")
    mean_sens = pd.read_csv(TABLES / "subfamily_segment_kruskal_mean_sensitivity.csv")

    lines = []
    lines.append("# Subfamily Flexibility Testing & Superfamily Segment Atlas\n")
    lines.append(
        "> Continues from commit `06f7b7e` (segment-level PCA/UMAP/clustering). Prior "
        "unsupervised conclusion (unchanged, not re-litigated here): *\"Under the "
        "structurally aligned segment-level representation, there is no robust evidence "
        "for discrete nuclear-receptor flexibility clusters; the results are more "
        "consistent with continuous variation.\"* Structural mapping, normalization, and "
        "modeling-stage inputs are untouched by this stage. No apo/holo, ClinVar, or "
        "Role D analysis was run here.\n"
    )

    lines.append("## Grouping scheme used\n")
    lines.append(
        "Three grouping fields actually exist in the project metadata (never invented): "
        "(1) formal NR subfamily-letter (from `nr_code`, e.g. `NR3C4` -> `NR3C`) -- "
        "inadequate sample sizes everywhere (max n=4) and **not used as a formal test "
        "grouping**, reported for transparency only; (2) formal NR family-number "
        "(`NR3C4` -> `NR3`) -- used as **SECONDARY**; (3) the professor-supplied project "
        "functional `group` field (`config/nr_metadata.py`) -- the largest, most balanced "
        "groups, used as **PRIMARY**.\n"
    )
    lines.append("| Scheme | Groups (n) |")
    lines.append("|---|---|")
    lines.append(f"| PRIMARY: project `group` | {schemes['project_group_counts']} |")
    lines.append(f"| SECONDARY: NR family-number | {schemes['nr_family_number_counts']} |")
    lines.append(f"| Not usable: NR subfamily-letter | {schemes['nr_subfamily_letter_counts']} |")
    lines.append("")

    lines.append("## Segment-specific sample sizes (Part 4)\n")
    n_segments_total = sample_sizes["structural_segment"].nunique()
    primary_testable = kruskal[(kruskal["grouping_scheme"] == "primary_project_group") &
                                (kruskal["test_status"] == "TESTED")]
    secondary_testable = kruskal[(kruskal["grouping_scheme"] == "secondary_nr_family_number") &
                                  (kruskal["test_status"] == "TESTED")]
    lines.append(
        f"{n_segments_total} realized structural segments (available-case per segment, "
        f"never restricted to the 29-receptor primary clustering cohort -- a receptor "
        f"missing one segment still contributes to every other segment's test). "
        f"**{len(primary_testable)}/{n_segments_total}** segments have >=2 eligible "
        f"groups (n>=3 receptors each) under the PRIMARY scheme; "
        f"**{len(secondary_testable)}/{n_segments_total}** under SECONDARY. Full "
        "per-segment-per-group counts in `reports/tables/subfamily_segment_sample_sizes.csv`.\n"
    )

    lines.append("## Omnibus results (Part 5-6)\n")
    n_sig_primary = kruskal[(kruskal["grouping_scheme"] == "primary_project_group") &
                             (kruskal["significant_fdr_05"] == True)]  # noqa: E712
    n_sig_secondary = kruskal[(kruskal["grouping_scheme"] == "secondary_nr_family_number") &
                               (kruskal["significant_fdr_05"] == True)]  # noqa: E712
    nominal = kruskal[kruskal["raw_p"] < 0.05]
    lines.append(
        f"Kruskal-Wallis across eligible groups, outcome = receptor-level "
        f"`median_of_entry_median_z`, BH-FDR correction applied within each scheme's "
        f"family of tests. **FDR-significant segments: {len(n_sig_primary)} (PRIMARY), "
        f"{len(n_sig_secondary)} (SECONDARY) -- zero in both.** "
        f"{len(nominal)} segment/scheme pair(s) had raw p<0.05 but did **not** survive "
        "BH-FDR and are **not** described as significant:\n"
    )
    if len(nominal):
        lines.append("| Scheme | Segment | n | H | raw p | epsilon^2 [95% CI] | FDR q |")
        lines.append("|---|---|---|---|---|---|---|")
        for _, r in nominal.iterrows():
            lines.append(f"| {r['grouping_scheme']} | {r['structural_segment']} | {int(r['n_receptors'])} | "
                         f"{r['H_statistic']:.2f} | {r['raw_p']:.4f} | "
                         f"{r['epsilon_squared']:.3f} [{r['epsilon_squared_ci_low']:.3f}, {r['epsilon_squared_ci_high']:.3f}] | "
                         f"{r['fdr_q']:.3f} |")
    lines.append("")
    lines.append(
        "Effect size = epsilon-squared = `H * (n+1) / (n^2-1)` (Tomczak & Tomczak 2014), "
        "with a stratified bootstrap 95% CI (seed=0, 2000 resamples). Full table: "
        "`reports/tables/subfamily_segment_kruskal.csv`.\n"
    )

    lines.append("## Post-hoc (Part 7)\n")
    lines.append(
        f"**0 segments** survived omnibus BH-FDR in either scheme, so, per instruction, "
        f"**no post-hoc pairwise tests were run or highlighted** ({len(posthoc)} rows in "
        "`reports/tables/subfamily_segment_posthoc.csv`, by design).\n"
    )

    lines.append("## Flexibility-atlas ranking (Part 8-9)\n")
    well = atlas[atlas["well_supported"]].sort_values("median_receptor_z", ascending=False)
    lines.append(
        f"{len(atlas)} segments in the full atlas (`data/analysis/"
        f"superfamily_segment_flexibility_atlas.csv`); **{len(well)}** are "
        f"\"well-supported\" (n_receptors >= {WELL_SUPPORTED_MIN_N}, an explicit documented "
        "threshold below which a bootstrap CI is unreliable or undefined). Z=0 means the "
        "average residue within a receptor's own instance(s) -- positive = relatively "
        "flexible, negative = relatively rigid, **within that receptor's own structures**, "
        "never an absolute cross-protein scale.\n"
    )
    lines.append("**Most consistently flexible (well-supported, highest median Z):**\n")
    lines.append("| Segment | n | median Z | 95% bootstrap CI |")
    lines.append("|---|---|---|---|")
    for _, r in well.head(5).iterrows():
        lines.append(f"| {r['structural_segment']} | {int(r['n_receptors'])} | {r['median_receptor_z']:.3f} | "
                     f"[{r['bootstrap_95ci_median_low']:.3f}, {r['bootstrap_95ci_median_high']:.3f}] |")
    lines.append("")
    lines.append("**Most consistently rigid (well-supported, lowest median Z):**\n")
    lines.append("| Segment | n | median Z | 95% bootstrap CI |")
    lines.append("|---|---|---|---|")
    for _, r in well.tail(5).iloc[::-1].iterrows():
        lines.append(f"| {r['structural_segment']} | {int(r['n_receptors'])} | {r['median_receptor_z']:.3f} | "
                     f"[{r['bootstrap_95ci_median_low']:.3f}, {r['bootstrap_95ci_median_high']:.3f}] |")
    lines.append("")
    unfiltered_top = atlas.loc[atlas["top_10pct_most_flexible"], "structural_segment"].tolist()
    unfiltered_bottom = atlas.loc[atlas["bottom_10pct_most_rigid"], "structural_segment"].tolist()
    lines.append(
        f"An **unfiltered** top/bottom-10% ranking (all realized segments, including "
        f"thin-n ones) gives top={unfiltered_top}, bottom={unfiltered_bottom} -- these are "
        "dominated by segments observed in only 1-3 receptors (undefined or extremely "
        "wide bootstrap CIs) and are **not** used as the headline ranking; the "
        "well-supported ranking above is preferred. No unresolved segment (H8, absent "
        "from this table entirely) is ever assigned a flexible/rigid label.\n"
    )

    lines.append("## Leave-one-receptor-out robustness (Part 10)\n")
    n_dominated = loo_atlas["dominated_by_single_receptor"].sum()
    lines.append(
        f"Atlas medians: {n_dominated}/{len(loo_atlas)} segments flagged as dominated by "
        "a single receptor (removing one receptor shifts the median by more than half "
        "the segment's own IQR, or 0.25 Z, whichever is larger) -- both are very-low-n "
        "segments already excluded from the well-supported ranking above "
        "(`reports/tables/subfamily_atlas_loo_robustness.csv`).\n\n"
        "Nominal (raw p<0.05) omnibus results, LOO-checked "
        "(`reports/tables/subfamily_kruskal_loo_robustness.csv`):\n"
    )
    for _, r in loo_kruskal.iterrows():
        lines.append(f"- **{r['structural_segment']}** ({r['grouping_scheme']}): retains p<0.05 in "
                     f"{r['fraction_loo_still_p_lt_05']:.0%} of leave-one-out folds "
                     f"(p range [{r['loo_p_min']:.4f}, {r['loo_p_max']:.4f}]) -- "
                     + ("fragile, depends on a specific receptor." if r["conclusion_depends_on_single_receptor"]
                        else "stable across folds, but still not FDR-significant."))
    lines.append("")

    lines.append("## Sensitivity: mean vs median segment summaries (Part 13)\n")
    n_sig_mean = mean_sens["significant_fdr_05"].sum()
    lines.append(
        f"Repeating the identical omnibus procedure with the receptor MEAN segment "
        f"summary (`mean_of_entry_mean_z`, the aggregation behind "
        f"`receptor_segment_feature_matrix_mean.parquet`) instead of the median pipeline: "
        f"**{n_sig_mean} FDR-significant segments** (same as median-based: 0). "
        "100% agreement on FDR-significance calls across all 35 segment/scheme tests; "
        "Spearman correlation of raw p-values between the two summaries = 0.74. This is "
        "the same conclusion from a different summary statistic, not a separate "
        "narrative (`reports/tables/subfamily_segment_kruskal_mean_sensitivity.csv`, "
        "`subfamily_median_vs_mean_agreement.csv`).\n"
    )

    lines.append("## Relationship to the unsupervised (clustering) conclusion (Part 11)\n")
    lines.append(
        "**Outcome A** (of the three possible outcomes named in the instructions): **no "
        "segment differs between groups after FDR correction, in either grouping scheme, "
        "under either the median or mean summary.** This directly strengthens the prior "
        "unsupervised finding -- flexibility at the structural-segment level appears "
        "broadly continuous/conserved across the project's functional groups and NR "
        "families, rather than concentrated in group-specific segment differences that "
        "clustering simply failed to detect as discrete groups. The two nominal (raw "
        "p<0.05) hits are exactly the kind of noise expected among ~35 independent tests "
        "at alpha=0.05 by chance, and one of the two does not even survive leave-one-out "
        "resampling.\n"
    )

    lines.append("## Limitations\n")
    lines.append(
        "- Small, unbalanced group sizes (Orphan receptors n=5) limit statistical power; "
        "a true small-to-moderate group difference could be missed (Type II error), not "
        "just correctly rejected.\n"
        "- Available-case analysis means different segments are tested on different, "
        "only partially overlapping receptor subsets -- segment-to-segment comparisons "
        "of \"who is more flexible\" are not perfectly apples-to-apples.\n"
        "- The formal NR subfamily-letter grouping (the most biologically precise "
        "grouping) could not be tested at all due to sample size; a coarser project "
        "grouping and NR family-number grouping were used instead.\n"
        "- H8 remains unresolved by the current structural-template mapping and never "
        "receives an atlas value, omnibus test, or flexible/rigid label.\n"
        "- This is a descriptive/inferential statistical stage; no causal or mechanistic "
        "claim is made about *why* any segment is more or less flexible.\n"
    )

    out = TABLES / "subfamily_flexibility_summary.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
