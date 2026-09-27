"""
Part 7 of the segment-level PCA/clustering modeling stage.

UMAP is used here strictly for VISUALIZATION, never as proof of discrete
clusters. Run on (a) the scaled primary matrix directly and (b), as a
sensitivity check, on the retained PCA representation (components covering
>=80% cumulative variance, from Part 6). Because n is small (~29
receptors), a parameter grid (n_neighbors x min_dist x seed) is swept and
a stability table reports whether the global visual pattern changes
substantially across settings -- operationalized as the Spearman
correlation between each setting's pairwise-embedding-distance matrix and
a fixed reference setting's. A UMAP setting is never chosen because it
"looks more separated" -- the default representative figure uses the
most literature-typical middle setting (n_neighbors=15, min_dist=0.1),
not a cherry-picked one.

Cluster-colored UMAP figures are produced later (Part 9 script) only
after clustering exists.
"""
import json
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from umap import UMAP

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "modeling"
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_NEIGHBORS_GRID = [5, 10, 15]
MIN_DIST_GRID = [0.1, 0.3, 0.5]
SEEDS = [0, 1, 2]
REFERENCE_SETTING = (15, 0.1, 0)  # (n_neighbors, min_dist, seed) used as the stability reference
DEFAULT_SETTING = (15, 0.1, 0)    # used for the representative saved figures


def load_inputs():
    scaled = pd.read_parquet(ANALYSIS / "receptor_segment_matrix_primary_scaled.parquet").set_index("uniprot_id")
    cohort = pd.read_csv(ANALYSIS / "modeling_receptor_cohort.csv").set_index("uniprot_id")
    pca_scores = pd.read_csv(ANALYSIS / "pca_receptor_scores.csv").set_index("uniprot_id")
    with open(ANALYSIS / "pca_summary.json") as f:
        pca_summary = json.load(f)
    n_pcs_80 = pca_summary["n_components_for_variance"]["80%"]
    pca_repr = pca_scores[[f"PC{i+1}" for i in range(n_pcs_80)]]
    return scaled, cohort, pca_repr, n_pcs_80


def run_umap_grid(X: np.ndarray, ids, label: str) -> tuple[pd.DataFrame, dict]:
    n_samples = X.shape[0]
    embeddings = {}
    rows = []
    for n_neighbors, min_dist, seed in product(N_NEIGHBORS_GRID, MIN_DIST_GRID, SEEDS):
        nn = min(n_neighbors, n_samples - 1)
        reducer = UMAP(n_neighbors=nn, min_dist=min_dist, n_components=2,
                        random_state=seed, metric="euclidean")
        emb = reducer.fit_transform(X)
        embeddings[(n_neighbors, min_dist, seed)] = emb
        rows.append({
            "input": label, "n_neighbors_requested": n_neighbors,
            "n_neighbors_used": nn, "min_dist": min_dist, "seed": seed,
        })

    ref_key = REFERENCE_SETTING
    if ref_key not in embeddings:
        ref_key = next(iter(embeddings))
    ref_dist = pdist(embeddings[ref_key])

    for row in rows:
        key = (row["n_neighbors_requested"], row["min_dist"], row["seed"])
        d = pdist(embeddings[key])
        if np.std(d) == 0 or np.std(ref_dist) == 0:
            rho = np.nan
        else:
            rho, _ = spearmanr(d, ref_dist)
        row["spearman_corr_to_reference_setting"] = rho

    # within-setting stability across seeds (same n_neighbors/min_dist, different seed)
    for n_neighbors, min_dist in product(N_NEIGHBORS_GRID, MIN_DIST_GRID):
        seed_embs = [embeddings[(n_neighbors, min_dist, s)] for s in SEEDS if (n_neighbors, min_dist, s) in embeddings]
        if len(seed_embs) < 2:
            continue
        dists = [pdist(e) for e in seed_embs]
        pairwise_rhos = []
        for i in range(len(dists)):
            for j in range(i + 1, len(dists)):
                if np.std(dists[i]) == 0 or np.std(dists[j]) == 0:
                    continue
                r, _ = spearmanr(dists[i], dists[j])
                pairwise_rhos.append(r)
        mean_rho = float(np.mean(pairwise_rhos)) if pairwise_rhos else np.nan
        for row in rows:
            if row["n_neighbors_requested"] == n_neighbors and row["min_dist"] == min_dist:
                row["mean_spearman_corr_across_seeds"] = mean_rho

    table = pd.DataFrame(rows)
    return table, embeddings


def plot_default(embeddings: dict, cohort: pd.DataFrame, ids, label: str, suffix: str):
    key = DEFAULT_SETTING if DEFAULT_SETTING in embeddings else next(iter(embeddings))
    emb = embeddings[key]
    df = pd.DataFrame(emb, columns=["UMAP1", "UMAP2"], index=ids).join(cohort[["nr_code", "group"]])

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(df["UMAP1"], df["UMAP2"], color="#444444")
    for uid, row in df.iterrows():
        ax.annotate(row["nr_code"], (row["UMAP1"], row["UMAP2"]), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_title(f"UMAP ({label}), n_neighbors={key[0]}, min_dist={key[1]}, seed={key[2]} -- by receptor")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"umap_{suffix}_by_receptor.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    groups = df["group"].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))
    for g, c in zip(groups, colors):
        sub = df[df["group"] == g]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], color=c, label=g, s=50)
    ax.set_title(f"UMAP ({label}), n_neighbors={key[0]}, min_dist={key[1]}, seed={key[2]} -- by group (metadata overlay)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"umap_{suffix}_by_group.png", dpi=150)
    plt.close(fig)

    return df, key


def main():
    scaled, cohort, pca_repr, n_pcs_80 = load_inputs()
    ids = scaled.index.tolist()
    print(f"UMAP grid: n_neighbors={N_NEIGHBORS_GRID}, min_dist={MIN_DIST_GRID}, seeds={SEEDS} "
          f"({len(N_NEIGHBORS_GRID) * len(MIN_DIST_GRID) * len(SEEDS)} runs per input)")

    print("\n--- Primary run: scaled 8-feature matrix ---")
    table_feat, emb_feat = run_umap_grid(scaled.values, ids, "scaled_feature_matrix")

    print(f"--- Sensitivity run: PCA representation ({n_pcs_80} PCs, >=80% variance) ---")
    table_pca, emb_pca = run_umap_grid(pca_repr.values, ids, "pca_representation")

    stability = pd.concat([table_feat, table_pca], ignore_index=True)
    stability_out = TABLES / "umap_stability.csv"
    stability.to_csv(stability_out, index=False)
    print(f"\nWrote {stability_out.relative_to(PROJECT_ROOT)}")

    agg = stability.groupby("input")["spearman_corr_to_reference_setting"].agg(["mean", "min", "max"])
    print("\nStability of global pattern vs the reference setting (n_neighbors=15, min_dist=0.1, seed=0):")
    print(agg.to_string())
    seed_agg = stability.drop_duplicates(["input", "n_neighbors_requested", "min_dist"]).groupby("input")[
        "mean_spearman_corr_across_seeds"].agg(["mean", "min"])
    print("\nStability across random seeds at fixed (n_neighbors, min_dist):")
    print(seed_agg.to_string())

    df_feat, key_feat = plot_default(emb_feat, cohort, ids, "scaled feature matrix", "feature_matrix")
    df_pca, key_pca = plot_default(emb_pca, cohort, ids, "PCA representation", "pca_representation")
    print(f"\nWrote representative figures for default setting {DEFAULT_SETTING} to {FIG_DIR.relative_to(PROJECT_ROOT)}/")

    # persist default embeddings for the clustering-colored figure later
    df_feat.reset_index().rename(columns={"index": "uniprot_id"})[
        ["uniprot_id", "UMAP1", "UMAP2"]
    ].to_csv(ANALYSIS / "umap_default_embedding_feature_matrix.csv", index=False)
    df_pca.reset_index().rename(columns={"index": "uniprot_id"})[
        ["uniprot_id", "UMAP1", "UMAP2"]
    ].to_csv(ANALYSIS / "umap_default_embedding_pca_representation.csv", index=False)

    overall_mean_stability = float(stability["spearman_corr_to_reference_setting"].mean())
    verdict = (
        "STABLE: global pattern correlates strongly (mean rho >= 0.7) across parameter settings"
        if overall_mean_stability >= 0.7 else
        "UNSTABLE: global pattern correlates weakly across settings; UMAP layout should not be "
        "over-interpreted as a robust global structure"
    )
    print(f"\nOverall UMAP stability verdict: {verdict} (mean rho = {overall_mean_stability:.3f})")

    summary = {
        "n_neighbors_grid": N_NEIGHBORS_GRID, "min_dist_grid": MIN_DIST_GRID, "seeds": SEEDS,
        "reference_setting": REFERENCE_SETTING, "default_setting": DEFAULT_SETTING,
        "overall_mean_spearman_to_reference": overall_mean_stability,
        "verdict": verdict,
    }
    (ANALYSIS / "umap_stability_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {(ANALYSIS / 'umap_stability_summary.json').relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
