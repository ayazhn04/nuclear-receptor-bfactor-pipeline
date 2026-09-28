"""
Part 2 of the functional-site specificity/validation stage: stratify the
existing within-receptor test results by site-definition quality, so
tier-A residue-level and tier-C segment-level proxy definitions are never
silently combined into one interpretation.

Tier A: RESIDUE-level definitions with HIGH confidence.
Tier B: RESIDUE-level definitions with MODERATE (or lower) confidence, or
        SEGMENT-level definitions with HIGH confidence ("partial").
Tier C: SEGMENT-level definitions with MODERATE/LOW confidence (the
        weakest evidence: an approximate structural-region proxy, not a
        residue-specific functional annotation).
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def classify_tier(definition_level: str, confidence: str) -> str:
    if definition_level == "RESIDUE" and confidence == "HIGH":
        return "A_residue_high_confidence"
    if definition_level == "RESIDUE" or (definition_level == "SEGMENT" and confidence == "HIGH"):
        return "B_partial_moderate"
    return "C_segment_level_proxy"


def main():
    within = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv")

    # within_receptor_tests.csv doesn't carry definition_level/confidence directly --
    # join from the receptor-level flexibility table, which does
    receptor_site = pd.read_parquet(ANALYSIS / "receptor_functional_site_flexibility.parquet")
    within = within.merge(
        receptor_site[["uniprot_id", "functional_site", "definition_level", "site_definition_confidence"]],
        on=["uniprot_id", "functional_site"], how="left",
    )
    within["definition_tier"] = within.apply(
        lambda r: classify_tier(r["definition_level"], r["site_definition_confidence"]), axis=1
    )

    rows = []
    for tier, g in within.groupby("definition_tier"):
        rows.append({
            "definition_tier": tier,
            "n_receptors": g["uniprot_id"].nunique(),
            "n_tests": len(g),
            "n_significant_within_fdr": int(g["significant_within_receptor"].sum()),
            "n_significant_global_fdr": int(g["significant_global"].sum()),
            "median_effect_size": float(g["effect_size_rank_biserial"].median()),
            "median_raw_p": float(g["raw_p"].median()),
        })
    stratified = pd.DataFrame(rows).sort_values("definition_tier")

    # also stratify per site class within each tier for the summary's finer detail
    rows2 = []
    for (tier, site), g in within.groupby(["definition_tier", "functional_site"]):
        rows2.append({
            "definition_tier": tier, "functional_site": site,
            "n_receptors": g["uniprot_id"].nunique(), "n_tests": len(g),
            "n_significant_within_fdr": int(g["significant_within_receptor"].sum()),
            "n_significant_global_fdr": int(g["significant_global"].sum()),
            "median_effect_size": float(g["effect_size_rank_biserial"].median()),
        })
    stratified_by_site = pd.DataFrame(rows2).sort_values(["functional_site", "definition_tier"])

    out = TABLES / "functional_site_definition_sensitivity.csv"
    combined = pd.concat([
        stratified.assign(functional_site="ALL_SITES"),
        stratified_by_site,
    ], ignore_index=True)
    combined = combined[["definition_tier", "functional_site", "n_receptors", "n_tests",
                          "n_significant_within_fdr", "n_significant_global_fdr",
                          "median_effect_size", "median_raw_p"]]
    combined.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(combined)} rows)")
    print("\nOverall by tier:")
    print(stratified.to_string(index=False))
    print("\nBy tier x site class:")
    print(stratified_by_site.to_string(index=False))

    # also write per-row tier assignment for downstream scripts to reuse
    tier_lookup_out = ANALYSIS / "functional_site_definition_tiers.parquet"
    within[["uniprot_id", "nr_code", "functional_site", "definition_level",
            "site_definition_confidence", "definition_tier"]].to_parquet(tier_lookup_out, index=False)
    print(f"\nWrote {tier_lookup_out.relative_to(PROJECT_ROOT)} (reusable tier lookup)")


if __name__ == "__main__":
    main()
