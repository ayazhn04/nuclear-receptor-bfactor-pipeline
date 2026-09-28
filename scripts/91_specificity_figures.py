"""
Part 11 of the functional-site specificity/validation stage: reproducible
figures under reports/figures/functional_sites_validation/.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "functional_sites_validation"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def fig_original_vs_matched_scatter(sensitivity: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(6, 6))
    colors = np.where(sensitivity["significance_preserved"], "#2ca02c",
              np.where(sensitivity["significance_lost"], "#d62728",
              np.where(sensitivity["significance_gained"], "#1f77b4", "#888888")))
    ax.scatter(sensitivity["original_delta"], sensitivity["matched_delta"], c=colors, s=40, alpha=0.8)
    lim = max(sensitivity["original_delta"].abs().max(), sensitivity["matched_delta"].abs().max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], color="gray", linestyle="--", linewidth=0.7)
    ax.axhline(0, color="lightgray", linewidth=0.5); ax.axvline(0, color="lightgray", linewidth=0.5)
    ax.set_xlabel("original delta (site - whole-LBD background)")
    ax.set_ylabel("matched delta (site - structural-context-matched background)")
    ax.set_title("Original vs. structurally matched-background effect\n"
                 "green=sig both, red=sig lost, blue=sig gained, gray=neither")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "original_vs_matched_effect_scatter.png", dpi=150)
    plt.close(fig)


def fig_effect_by_confidence(within: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(7, 5))
    tiers = ["A_residue_high_confidence", "B_partial_moderate", "C_segment_level_proxy"]
    data = [within.loc[within["definition_tier"] == t, "effect_size_rank_biserial"].dropna().values for t in tiers]
    present = [(t, d) for t, d in zip(tiers, data) if len(d) > 0]
    ax.boxplot([d for _, d in present], tick_labels=[t for t, _ in present])
    ax.axhline(0, color="lightgray", linewidth=0.7)
    ax.set_ylabel("effect size (rank-biserial), original whole-LBD-background test")
    ax.set_title("Effect size by site-definition confidence tier")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "effect_by_definition_confidence.png", dpi=150)
    plt.close(fig)


def fig_matched_superfamily_forest(matched_superfamily: pd.DataFrame):
    df = matched_superfamily[matched_superfamily["status"] == "TESTED"].sort_values("median_matched_delta")
    fig, ax = plt.subplots(figsize=(10, 3.5))
    y = np.arange(len(df))
    colors = ["#2ca02c" if s else "#888888" for s in df["significant_fdr_05"]]
    ax.errorbar(df["median_matched_delta"], y,
                xerr=[df["median_matched_delta"] - df["bootstrap_ci_low"],
                      df["bootstrap_ci_high"] - df["median_matched_delta"]],
                fmt="none", ecolor="#cccccc", capsize=3)
    ax.scatter(df["median_matched_delta"], y, color=colors, s=60, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(df["functional_site"])
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_xlabel("median matched-background delta + bootstrap 95% CI")
    ax.set_title("Superfamily CONFIRMATORY reanalysis (structural-context-matched background)\n"
                 "HIGH/MODERATE-confidence RESIDUE-level definitions only; green = FDR-significant")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "superfamily_matched_forest_plot.png", dpi=150)
    plt.close(fig)


def fig_classification_counts(classification: pd.DataFrame):
    counts = classification["interpretation_category"].value_counts()
    fig, ax = plt.subplots(figsize=(7, 4))
    colors_map = {
        "SITE_SPECIFIC_SUPPORTED": "#2ca02c", "STRUCTURAL_CONTEXT_EXPLAINED": "#d62728",
        "INSUFFICIENT_RESIDUE_LEVEL_MAPPING": "#888888", "MIXED": "#ff7f0e",
    }
    ax.bar(counts.index, counts.values, color=[colors_map.get(c, "#333333") for c in counts.index])
    ax.set_ylabel("n site classes")
    ax.set_title("Specificity validation outcome by site class")
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "interpretation_category_counts.png", dpi=150)
    plt.close(fig)


def fig_paired_examples(pdb_matched: pd.DataFrame):
    counts = pdb_matched.groupby(["uniprot_id", "functional_site"]).size()
    well = counts[counts >= 15].sort_values(ascending=False)
    if len(well) == 0:
        print("No receptor x site pair with >=15 matched-background PDB entries -- skipping paired example figure")
        return
    top = well.index[0]
    uid, site = top
    g = pdb_matched[(pdb_matched["uniprot_id"] == uid) & (pdb_matched["functional_site"] == site)]
    nr_code = g["nr_code"].iloc[0]

    fig, ax = plt.subplots(figsize=(5, 5))
    for _, r in g.iterrows():
        ax.plot([0, 1], [r["median_matched_background_z"], r["median_site_z"]], color="gray", alpha=0.4, linewidth=0.8)
    ax.scatter(np.zeros(len(g)), g["median_matched_background_z"], color="#888888", zorder=3, label="matched background")
    ax.scatter(np.ones(len(g)), g["median_site_z"], color="#d62728", zorder=3, label=site)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["matched\nbackground", site])
    ax.set_ylabel("median Z per PDB entry")
    ax.set_title(f"{nr_code} ({uid}): {site} vs. structurally matched background\n"
                 f"{len(g)} paired PDB entries, tier={g['matched_tier_used'].iloc[0]}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"paired_matched_example_{nr_code}_{site}.png", dpi=150)
    plt.close(fig)


def main():
    sensitivity = pd.read_csv(TABLES / "functional_site_background_sensitivity.csv")
    within = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv").merge(
        pd.read_parquet(ANALYSIS / "functional_site_definition_tiers.parquet")[
            ["uniprot_id", "functional_site", "definition_tier"]],
        on=["uniprot_id", "functional_site"], how="left",
    )
    matched_superfamily = pd.read_csv(TABLES / "functional_site_superfamily_matched_reanalysis.csv")
    classification = pd.read_csv(TABLES / "functional_site_interpretation_classification.csv")
    pdb_matched = pd.read_parquet(ANALYSIS / "pdb_functional_site_matched_background.parquet")

    fig_original_vs_matched_scatter(sensitivity)
    fig_effect_by_confidence(within)
    fig_matched_superfamily_forest(matched_superfamily)
    fig_classification_counts(classification)
    fig_paired_examples(pdb_matched)

    print(f"Wrote figures to {FIG_DIR.relative_to(PROJECT_ROOT)}/:")
    for f in sorted(FIG_DIR.glob("*.png")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
