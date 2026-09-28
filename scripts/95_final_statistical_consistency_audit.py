"""
Part 4 of the freeze/consolidate stage: cross-check every major claim in
the project's final Markdown summaries against a value INDEPENDENTLY
RECOMPUTED from its machine-readable source table (never a hand-typed
expected value) -- significant counts, q<0.05 counts, receptor/sample
counts, PCA variance, cluster metrics, atlas rankings, apo/holo names and
deltas, functional-site classifications.

Each check: (claim description, recomputed value, exact substring the
recomputed value should produce in the target .md file, PASS/FAIL).
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

results = []


def check(claim, md_file, expected_substring, source_file, method):
    text = (PROJECT_ROOT / md_file).read_text()
    passed = expected_substring in text
    results.append({
        "claim": claim, "md_file": md_file, "recomputed_expected_substring": expected_substring,
        "source_file": source_file, "recomputation_method": method,
        "status": "PASS" if passed else "FAIL",
    })
    return passed


def main():
    # ---- PCA (modeling_summary.md) ----
    pca_var = pd.read_csv(TABLES / "pca_variance.csv")
    pc1 = pca_var.iloc[0]["explained_variance_ratio"]
    check("PC1 explained variance", "reports/tables/modeling_summary.md", f"{pc1:.1%}".replace("%", ""),
          "reports/tables/pca_variance.csv", "pca_var.iloc[0]['explained_variance_ratio']")

    # ---- Clustering (modeling_summary.md) ----
    dbscan = pd.read_csv(TABLES / "clustering_dbscan_grid.csv")
    max_clusters = int(dbscan[dbscan["input"] == "pca_primary"]["n_clusters_found"].max())
    check("DBSCAN max clusters found (PCA-primary)", "reports/tables/modeling_summary.md", f"**{max_clusters}**",
          "reports/tables/clustering_dbscan_grid.csv", "max(n_clusters_found) where input==pca_primary")

    clustering_val = pd.read_csv(TABLES / "clustering_validation.csv")
    best = clustering_val[clustering_val["input"] == "pca_primary"].dropna(subset=["silhouette"]).sort_values(
        "silhouette", ascending=False).iloc[0]
    check("Best silhouette value (PCA-primary)", "reports/tables/modeling_summary.md",
          f"silhouette={best['silhouette']:.3f}", "reports/tables/clustering_validation.csv",
          "max(silhouette) where input==pca_primary")

    # ---- Subfamily testing (subfamily_flexibility_summary.md) ----
    kruskal = pd.read_csv(TABLES / "subfamily_segment_kruskal.csv")
    n_sig_primary = int((kruskal[kruskal["grouping_scheme"] == "primary_project_group"]["significant_fdr_05"] == True).sum())  # noqa: E712
    check("Subfamily PRIMARY-grouping FDR-significant segment count",
          "reports/tables/subfamily_flexibility_summary.md", f"significant segments: {n_sig_primary} (PRIMARY)",
          "reports/tables/subfamily_segment_kruskal.csv", "sum(significant_fdr_05) where scheme==primary_project_group")

    atlas = pd.read_csv(ANALYSIS / "superfamily_segment_flexibility_atlas.csv")
    well = atlas[atlas["well_supported"]].sort_values("median_receptor_z", ascending=False)
    top_segment = well.iloc[0]["structural_segment"]
    check("Atlas most-flexible well-supported segment (rank 1)",
          "reports/tables/subfamily_flexibility_summary.md", top_segment,
          "data/analysis/superfamily_segment_flexibility_atlas.csv",
          "sort_values(median_receptor_z, desc) among well_supported, iloc[0]")

    # ---- apo/holo (apo_holo_summary.md) ----
    elig = pd.read_csv(TABLES / "apo_holo_receptor_eligibility.csv")
    n_primary = int(elig["eligible_primary_5_5"].sum())
    n_exploratory = int(elig["eligible_exploratory_3_3"].sum())
    check("apo/holo PRIMARY-eligible receptor count", "reports/tables/apo_holo_summary.md",
          f"PRIMARY (>=5 APO and >=5 HOLO PDB entries): {n_primary}/{len(elig)} receptors eligible",
          "reports/tables/apo_holo_receptor_eligibility.csv", "sum(eligible_primary_5_5)")
    check("apo/holo EXPLORATORY-eligible receptor count", "reports/tables/apo_holo_summary.md",
          f"EXPLORATORY (>=3/>=3): {n_exploratory}/{len(elig)} eligible",
          "reports/tables/apo_holo_receptor_eligibility.csv", "sum(eligible_exploratory_3_3)")

    seg_tests = pd.read_csv(TABLES / "apo_holo_segment_tests.csv")
    n_within = int(seg_tests["significant_within_receptor"].sum())
    n_global = int(seg_tests["significant_global"].sum())
    check("apo/holo within-receptor FDR significant count", "reports/tables/apo_holo_summary.md",
          f"**{n_within}** significant within-receptor FDR",
          "reports/tables/apo_holo_segment_tests.csv", "sum(significant_within_receptor)")
    check("apo/holo global FDR significant count", "reports/tables/apo_holo_summary.md",
          f"**{n_global}** significant under the SECONDARY global FDR",
          "reports/tables/apo_holo_segment_tests.csv", "sum(significant_global)")

    sig = seg_tests[seg_tests["significant_within_receptor"]]
    for _, r in sig.iterrows():
        check(f"apo/holo delta for {r['nr_code']}/{r['structural_segment']}", "reports/tables/apo_holo_summary.md",
              f"{r['delta_median_z']:+.3f}", "reports/tables/apo_holo_segment_tests.csv",
              f"row where nr_code=={r['nr_code']} & structural_segment=={r['structural_segment']}, delta_median_z")

    # ---- functional sites (current_functional_site_summary.md) ----
    within = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv")
    n_within_fs = int(within["significant_within_receptor"].sum())
    check("Functional-site within-receptor FDR significant count",
          "reports/tables/current_functional_site_summary.md", f"**{n_within_fs}**",
          "reports/tables/current_functional_site_within_receptor_tests.csv", "sum(significant_within_receptor)")

    superfamily = pd.read_csv(TABLES / "current_functional_site_superfamily_tests.csv")
    for _, r in superfamily.iterrows():
        check(f"Superfamily site-class delta for {r['functional_site']}",
              "reports/tables/current_functional_site_summary.md", f"{r['median_delta_vs_background']:+.3f}",
              "reports/tables/current_functional_site_superfamily_tests.csv",
              f"row where functional_site=={r['functional_site']}, median_delta_vs_background")

    # ---- specificity validation (functional_site_specificity_summary.md) ----
    classif = pd.read_csv(TABLES / "functional_site_interpretation_classification.csv")
    for _, r in classif.iterrows():
        check(f"Interpretation classification for {r['functional_site']}",
              "reports/tables/functional_site_specificity_summary.md", r["interpretation_category"],
              "reports/tables/functional_site_interpretation_classification.csv", "interpretation_category")

    audit_df = pd.DataFrame(results)
    n_pass = (audit_df["status"] == "PASS").sum()
    n_fail = (audit_df["status"] == "FAIL").sum()

    lines = ["# Final Statistical Consistency Audit\n",
             f"> {len(audit_df)} claims cross-checked: every value below was INDEPENDENTLY RECOMPUTED from its "
             "machine-readable source table (never hand-typed) and confirmed to appear correctly in the "
             f"corresponding Markdown summary. **{n_pass}/{len(audit_df)} PASS, {n_fail} FAIL.**\n"]
    lines.append("| Claim | Recomputed value | Target file | Status |")
    lines.append("|---|---|---|---|")
    for _, r in audit_df.iterrows():
        status = "✅ PASS" if r["status"] == "PASS" else "❌ FAIL"
        lines.append(f"| {r['claim']} | `{r['recomputed_expected_substring']}` | {r['md_file']} | {status} |")

    if n_fail:
        lines.append("\n## Failures require investigation\n")
        for _, r in audit_df[audit_df["status"] == "FAIL"].iterrows():
            lines.append(f"- **{r['claim']}**: expected substring `{r['recomputed_expected_substring']}` "
                         f"(from {r['source_file']}, via {r['recomputation_method']}) not found in {r['md_file']}")

    out_md = PROJECT_ROOT / "reports" / "tables" / "final_statistical_consistency_audit.md"
    out_md.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_md.relative_to(PROJECT_ROOT)}")
    print(f"\n{n_pass}/{len(audit_df)} PASS, {n_fail} FAIL")
    if n_fail:
        print(audit_df[audit_df["status"] == "FAIL"].to_string(index=False))


if __name__ == "__main__":
    main()
