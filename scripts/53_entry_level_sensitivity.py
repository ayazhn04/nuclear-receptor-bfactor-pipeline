"""
Part 13 of the segment-level PCA/clustering modeling stage.

The primary scientific unit throughout this stage is the RECEPTOR (scripts
48-52). This script is a SECONDARY sensitivity look at the 1789-entry
PDB-level matrix, used only to examine within-receptor conformational
heterogeneity -- never to claim subfamily separation, and never treating
PDB entries as independent receptor-family observations (a receptor with
hundreds of entries, e.g. PPARgamma/ESR1, would otherwise dominate any
superfamily-level pattern purely by structure count).

Produces: an entry-count-per-receptor table (the overrepresentation this
guards against), and a descriptive PCA embedding of PDB entries on the
same 8 primary segment features, colored by receptor. No clustering
battery is run here and no biological conclusion is drawn from it.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "modeling"


def main():
    with open(ANALYSIS / "modeling_cohort_matrix_summary.json") as f:
        cohort_summary = json.load(f)
    primary_features = cohort_summary["primary_features"]

    entry = pd.read_parquet(ANALYSIS / "pdb_segment_feature_matrix_median.parquet")
    cohort = pd.read_csv(ANALYSIS / "modeling_receptor_cohort.csv")

    # --- overrepresentation table (all 42 coordinate receptors, all PDB entries) ---
    counts = entry.groupby("uniprot_id").size().rename("n_pdb_entries")
    counts = counts.reindex(cohort["uniprot_id"]).fillna(0).astype(int)
    overrep = cohort.set_index("uniprot_id")[["nr_code", "group", "included_primary"]].join(counts)
    overrep = overrep.sort_values("n_pdb_entries", ascending=False)
    overrep_out = TABLES / "entry_level_overrepresentation.csv"
    overrep.reset_index().to_csv(overrep_out, index=False)
    print(f"Wrote {overrep_out.relative_to(PROJECT_ROOT)}")
    print(f"PDB-entry count per receptor: min={counts.min()}, median={counts.median():.0f}, "
          f"max={counts.max()} ({overrep.index[0]}/{overrep.iloc[0]['nr_code']})")
    print(f"Top 5 most-represented receptors account for "
          f"{overrep['n_pdb_entries'].head(5).sum() / overrep['n_pdb_entries'].sum():.1%} of all {int(counts.sum())} entries "
          f"-- this is exactly why the receptor-level median matrix, not this entry-level matrix, "
          f"is the primary clustering input.")

    # --- descriptive entry-level embedding (secondary only) ---
    entry_primary = entry[["uniprot_id", "pdb_id"] + primary_features].copy()
    entry_primary = entry_primary.dropna(subset=primary_features, how="all")
    # require every primary feature observed for this specific entry (a stricter,
    # honest complete-case view at entry level -- no imputation invented here)
    complete_entries = entry_primary.dropna(subset=primary_features, how="any")
    print(f"\n{len(complete_entries)}/{len(entry)} PDB entries have all {len(primary_features)} "
          f"primary features observed (complete-case, no entry-level imputation performed)")

    if len(complete_entries) >= 10:
        X = complete_entries[primary_features].values
        X = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
        pca = PCA(n_components=2, random_state=0)
        scores = pca.fit_transform(X)
        plot_df = complete_entries[["uniprot_id", "pdb_id"]].copy()
        plot_df["PC1"], plot_df["PC2"] = scores[:, 0], scores[:, 1]
        plot_df = plot_df.merge(cohort[["uniprot_id", "nr_code", "group"]], on="uniprot_id", how="left")

        fig, ax = plt.subplots(figsize=(8, 6))
        receptors = plot_df["uniprot_id"].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(receptors)))
        for uid, c in zip(receptors, colors):
            sub = plot_df[plot_df["uniprot_id"] == uid]
            ax.scatter(sub["PC1"], sub["PC2"], color=c, s=15, alpha=0.7,
                       label=f"{sub['nr_code'].iloc[0]} (n={len(sub)})")
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} var)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} var)")
        ax.set_title(f"Entry-level PCA (n={len(complete_entries)} PDB entries), colored by receptor\n"
                     "SECONDARY/descriptive only -- not used to claim subfamily separation "
                     "(receptor overrepresentation is visible by construction)")
        ax.legend(fontsize=6, ncol=2, loc="center left", bbox_to_anchor=(1.0, 0.5))
        fig.tight_layout()
        fig.savefig(FIG_DIR / "entry_level_pca_by_receptor.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {(FIG_DIR / 'entry_level_pca_by_receptor.png').relative_to(PROJECT_ROOT)}")
    else:
        print("Fewer than 10 complete-case PDB entries -- skipping the descriptive entry-level "
              "PCA figure (would not be a meaningful embedding).")

    print("\nNo clustering was run at the entry level, and no subfamily-separation claim is made "
          "from this secondary view -- receptor is the unit of analysis for all formal conclusions.")


if __name__ == "__main__":
    main()
