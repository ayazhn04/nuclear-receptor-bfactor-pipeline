"""
Part 6 of the segment-level PCA/clustering modeling stage.

Runs PCA on the scaled PRIMARY receptor matrix
(receptor_segment_matrix_primary_scaled.parquet, 29 receptors x 8 segment
features). n is only ~29 receptors and the feature count is small (8), so
n_components is never hard-coded to 50 -- it is min(n_samples, n_features)
minus 1 in practice via sklearn's default full SVD, i.e. up to 8 components
here (bounded by the smaller of n_samples-1 and n_features). Group labels
are NEVER used for fitting -- PCA is fit only on the 8 selected flexibility
features -- and are joined back in afterward purely for visualization.
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
FIG_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 0  # PCA itself is deterministic (full SVD); kept for any solver fallback


def load_scaled_matrix():
    scaled = pd.read_parquet(ANALYSIS / "receptor_segment_matrix_primary_scaled.parquet")
    scaled = scaled.set_index("uniprot_id")
    return scaled


def load_cohort_metadata():
    cohort = pd.read_csv(ANALYSIS / "modeling_receptor_cohort.csv")
    return cohort.set_index("uniprot_id")


def run_pca(scaled: pd.DataFrame):
    n_samples, n_features = scaled.shape
    max_components = min(n_samples, n_features)
    pca = PCA(n_components=max_components, svd_solver="full", random_state=RANDOM_SEED)
    scores = pca.fit_transform(scaled.values)
    return pca, scores, max_components


def main():
    scaled = load_scaled_matrix()
    cohort = load_cohort_metadata()
    feature_names = list(scaled.columns)
    n_samples, n_features = scaled.shape
    print(f"PCA input: {n_samples} receptors (primary cohort) x {n_features} scaled segment features")
    print(f"Features (fit only on these -- group/subfamily labels NOT used for fitting): {feature_names}")

    pca, scores, max_components = run_pca(scaled)
    explained = pca.explained_variance_ratio_
    cumulative = np.cumsum(explained)
    eigenvalues = pca.explained_variance_

    variance_table = pd.DataFrame({
        "pc": [f"PC{i+1}" for i in range(max_components)],
        "eigenvalue": eigenvalues,
        "explained_variance_ratio": explained,
        "cumulative_explained_variance": cumulative,
    })
    variance_out = TABLES / "pca_variance.csv"
    variance_table.to_csv(variance_out, index=False)
    print(f"\nWrote {variance_out.relative_to(PROJECT_ROOT)}")
    print(variance_table.to_string(index=False))

    n_for = {}
    for target in (0.70, 0.80, 0.90, 0.95):
        idx = np.searchsorted(cumulative, target) + 1
        idx = min(idx, max_components)
        n_for[target] = int(idx)
    print("\nComponents needed for variance-explained targets:")
    for target, n in n_for.items():
        print(f"  {target:.0%}: {n} PC(s) (cumulative = {cumulative[n-1]:.1%})")

    # receptor scores
    scores_df = pd.DataFrame(scores, index=scaled.index,
                              columns=[f"PC{i+1}" for i in range(max_components)])
    scores_df = scores_df.join(cohort[["nr_code", "group"]])
    scores_out = ANALYSIS / "pca_receptor_scores.csv"
    scores_df.reset_index().rename(columns={"index": "uniprot_id"}).to_csv(scores_out, index=False)
    print(f"\nWrote {scores_out.relative_to(PROJECT_ROOT)}")

    # loadings
    loadings = pd.DataFrame(
        pca.components_.T, index=feature_names,
        columns=[f"PC{i+1}" for i in range(max_components)],
    )
    loadings_out = ANALYSIS / "pca_loadings.csv"
    loadings.reset_index().rename(columns={"index": "structural_segment"}).to_csv(loadings_out, index=False)
    print(f"Wrote {loadings_out.relative_to(PROJECT_ROOT)}")

    strongest = loadings["PC1"].abs().sort_values(ascending=False)
    print("\nStrongest |PC1| loadings:")
    print(strongest.to_string())

    # ---- figures ----
    # scree plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(range(1, max_components + 1), explained, "o-", color="#1f77b4")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Explained variance ratio")
    ax.set_title("PCA scree plot (primary cohort, n=%d, %d features)" % (n_samples, n_features))
    ax.set_xticks(range(1, max_components + 1))
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pca_scree.png", dpi=150)
    plt.close(fig)

    # cumulative variance
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(range(1, max_components + 1), cumulative, "o-", color="#2ca02c")
    for target in (0.70, 0.80, 0.90, 0.95):
        ax.axhline(target, color="gray", linestyle="--", linewidth=0.7)
    ax.set_xlabel("Number of PCs")
    ax.set_ylabel("Cumulative explained variance")
    ax.set_title("PCA cumulative explained variance")
    ax.set_xticks(range(1, max_components + 1))
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pca_cumulative_variance.png", dpi=150)
    plt.close(fig)

    # PC1 vs PC2 by receptor
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(scores_df["PC1"], scores_df["PC2"], color="#444444")
    for uid, row in scores_df.iterrows():
        ax.annotate(row["nr_code"], (row["PC1"], row["PC2"]), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({explained[0]:.1%} var)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1%} var)")
    ax.set_title("PC1 vs PC2, labeled by receptor (nr_code)")
    ax.axhline(0, color="lightgray", linewidth=0.5)
    ax.axvline(0, color="lightgray", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pca_pc1_pc2_by_receptor.png", dpi=150)
    plt.close(fig)

    # PC1 vs PC2 by group (metadata overlay only, NOT used in fitting)
    fig, ax = plt.subplots(figsize=(7, 6))
    groups = scores_df["group"].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))
    for g, c in zip(groups, colors):
        sub = scores_df[scores_df["group"] == g]
        ax.scatter(sub["PC1"], sub["PC2"], color=c, label=g, s=50)
    ax.set_xlabel(f"PC1 ({explained[0]:.1%} var)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1%} var)")
    ax.set_title("PC1 vs PC2, colored by project NR group (metadata overlay, not fit input)")
    ax.axhline(0, color="lightgray", linewidth=0.5)
    ax.axvline(0, color="lightgray", linewidth=0.5)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pca_pc1_pc2_by_group.png", dpi=150)
    plt.close(fig)

    # loadings heatmap
    fig, ax = plt.subplots(figsize=(8, 5))
    n_show = min(5, max_components)
    im = ax.imshow(loadings.iloc[:, :n_show].values, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(n_show))
    ax.set_xticklabels([f"PC{i+1}" for i in range(n_show)])
    ax.set_yticks(range(len(feature_names)))
    ax.set_yticklabels(feature_names, fontsize=8)
    ax.set_title("PCA loadings heatmap (first %d PCs)" % n_show)
    fig.colorbar(im, ax=ax, label="loading")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pca_loadings_heatmap.png", dpi=150)
    plt.close(fig)

    print(f"\nWrote 5 figures to {FIG_DIR.relative_to(PROJECT_ROOT)}/")
    print("\nNote: visual PC1/PC2 overlap or separation is descriptive only, "
          "not formal clustering evidence -- that comes in the clustering stage.")

    summary = {
        "n_samples": n_samples,
        "n_features": n_features,
        "max_components": max_components,
        "n_components_for_variance": {f"{k:.0%}": v for k, v in n_for.items()},
        "pc1_explained_variance": float(explained[0]),
        "pc2_explained_variance": float(explained[1]) if max_components > 1 else None,
    }
    summary_out = ANALYSIS / "pca_summary.json"
    summary_out.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {summary_out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
