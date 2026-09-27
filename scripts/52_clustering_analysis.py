"""
Parts 8-11 of the segment-level PCA/clustering modeling stage.

Part 8: PRIMARY clustering input = PCA scores retaining components for
>=80% cumulative variance (from Part 6, 5 PCs). A sensitivity run repeats
the identical battery directly on the scaled 8-feature matrix. Clustering
is never run on 2D UMAP coordinates -- UMAP is visualization-only.

Part 9: full planned battery (K-means, Agglomerative/Ward, GMM, Spectral,
DBSCAN), k = 2..min(10, n_receptors-1) for the four k-based methods, a
documented eps/min_samples grid for DBSCAN (never forced to find clusters).

Part 10: Silhouette / Davies-Bouldin / Calinski-Harabasz for every
applicable solution; Gap Statistic for K-means; AIC/BIC for GMM; pairwise
Adjusted Rand Index / Rand Index across methods at matched k; stability
checks across imputation method, feature-support threshold, median-vs-mean
matrix, and PCA-space vs scaled-feature-space.

Part 11: no method is picked because it "looks nicest" -- every metric for
every method/k/input is reported in reports/tables/clustering_validation.csv,
and the continuum-vs-discrete-clusters verdict in the modeling summary is
derived directly from these numbers, not asserted first.

Random state is fixed everywhere a method uses one, so outputs (including
DBSCAN, which has no randomness at all) are deterministic across reruns.
"""
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering, DBSCAN
from sklearn.mixture import GaussianMixture
from sklearn.metrics import (
    silhouette_score, davies_bouldin_score, calinski_harabasz_score,
    adjusted_rand_score, rand_score,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "modeling"

RANDOM_STATE = 0
GAP_N_REFS = 20  # bootstrap reference datasets for the gap statistic


# ---------------------------------------------------------------- inputs ---
def load_primary_input():
    pca_scores = pd.read_csv(ANALYSIS / "pca_receptor_scores.csv").set_index("uniprot_id")
    with open(ANALYSIS / "pca_summary.json") as f:
        pca_summary = json.load(f)
    n_pcs_80 = pca_summary["n_components_for_variance"]["80%"]
    X = pca_scores[[f"PC{i+1}" for i in range(n_pcs_80)]]
    return X


def load_sensitivity_inputs():
    """Returns a dict of {label: DataFrame(index=uniprot_id)} for every
    sensitivity comparison requested in Part 10."""
    inputs = {}
    inputs["scaled_feature_matrix"] = pd.read_parquet(
        ANALYSIS / "receptor_segment_matrix_primary_scaled.parquet"
    ).set_index("uniprot_id")
    inputs["knn_imputed_scaled"] = pd.read_parquet(
        ANALYSIS / "receptor_segment_matrix_sensitivity_knn_scaled.parquet"
    ).set_index("uniprot_id")
    inputs["meanmatrix_scaled"] = pd.read_parquet(
        ANALYSIS / "receptor_segment_matrix_sensitivity_meanmatrix_scaled.parquet"
    ).set_index("uniprot_id")
    complete_case_scaled_path = ANALYSIS / "receptor_segment_matrix_sensitivity_complete_case_scaled.parquet"
    if complete_case_scaled_path.exists():
        inputs["complete_case_scaled"] = pd.read_parquet(complete_case_scaled_path).set_index("uniprot_id")
    return inputs


# --------------------------------------------------------- gap statistic ---
def gap_statistic(X: np.ndarray, k_range, n_refs=GAP_N_REFS, seed=RANDOM_STATE):
    rng = np.random.default_rng(seed)
    mins, maxs = X.min(axis=0), X.max(axis=0)

    def dispersion(data, k):
        km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(data)
        labels = km.labels_
        # standard Wk = sum over clusters of within-cluster pairwise distance / (2*n_c)
        W = 0.0
        for c in np.unique(labels):
            members = data[labels == c]
            n_c = len(members)
            if n_c < 2:
                continue
            D = squareform(pdist(members))
            W += D.sum() / (2 * n_c)
        return W

    results = []
    for k in k_range:
        Wk = dispersion(X, k)
        ref_logs = []
        for _ in range(n_refs):
            ref = rng.uniform(mins, maxs, size=X.shape)
            ref_logs.append(np.log(max(dispersion(ref, k), 1e-12)))
        ref_logs = np.array(ref_logs)
        gap = ref_logs.mean() - np.log(max(Wk, 1e-12))
        sk = ref_logs.std(ddof=0) * np.sqrt(1 + 1.0 / n_refs)
        results.append({"k": k, "gap": gap, "sk": sk, "log_Wk": np.log(max(Wk, 1e-12)),
                         "ref_log_Wk_mean": ref_logs.mean()})
    df = pd.DataFrame(results)
    # Tibshirani rule: smallest k such that gap(k) >= gap(k+1) - sk(k+1)
    chosen_k = None
    for i in range(len(df) - 1):
        if df.loc[i, "gap"] >= df.loc[i + 1, "gap"] - df.loc[i + 1, "sk"]:
            chosen_k = int(df.loc[i, "k"])
            break
    return df, chosen_k


# ------------------------------------------------------- clustering runs ---
def run_kbased_battery(X: pd.DataFrame, k_range, input_label: str):
    ids = X.index.tolist()
    Xv = X.values
    n = Xv.shape[0]
    rows = []
    labels_store = {}

    for k in k_range:
        if k >= n:
            continue

        km = KMeans(n_clusters=k, n_init=10, random_state=RANDOM_STATE).fit(Xv)
        labels_store[("kmeans", k)] = pd.Series(km.labels_, index=ids)

        agg = AgglomerativeClustering(n_clusters=k, linkage="ward").fit(Xv)
        labels_store[("agglomerative_ward", k)] = pd.Series(agg.labels_, index=ids)

        gmm = GaussianMixture(n_components=k, covariance_type="diag", n_init=5,
                               random_state=RANDOM_STATE, reg_covar=1e-4).fit(Xv)
        gmm_labels = gmm.predict(Xv)
        labels_store[("gmm", k)] = pd.Series(gmm_labels, index=ids)

        n_neighbors_spec = max(2, min(10, n - 1))
        try:
            spec = SpectralClustering(n_clusters=k, affinity="nearest_neighbors",
                                       n_neighbors=n_neighbors_spec, assign_labels="kmeans",
                                       random_state=RANDOM_STATE).fit(Xv)
            labels_store[("spectral", k)] = pd.Series(spec.labels_, index=ids)
        except Exception as e:
            print(f"  spectral clustering failed at k={k} ({input_label}): {e}")

        for method_key in ("kmeans", "agglomerative_ward", "gmm", "spectral"):
            if (method_key, k) not in labels_store:
                continue
            labels = labels_store[(method_key, k)].values
            n_unique = len(np.unique(labels))
            row = {"input": input_label, "method": method_key, "k": k, "n_clusters_found": n_unique}
            if n_unique >= 2:
                row["silhouette"] = silhouette_score(Xv, labels)
                row["davies_bouldin"] = davies_bouldin_score(Xv, labels)
                row["calinski_harabasz"] = calinski_harabasz_score(Xv, labels)
            else:
                row["silhouette"] = np.nan
                row["davies_bouldin"] = np.nan
                row["calinski_harabasz"] = np.nan
            if method_key == "gmm":
                row["aic"] = gmm.aic(Xv)
                row["bic"] = gmm.bic(Xv)
            rows.append(row)

    validation = pd.DataFrame(rows)

    gap_df, gap_chosen_k = gap_statistic(Xv, [k for k in k_range if k < n])
    gap_df["input"] = input_label
    gap_df["method"] = "kmeans_gap_statistic"

    return validation, gap_df, gap_chosen_k, labels_store


def run_dbscan_grid(X: pd.DataFrame, input_label: str):
    ids = X.index.tolist()
    Xv = X.values
    D = squareform(pdist(Xv))
    np.fill_diagonal(D, np.nan)
    nn_dist = np.nanmin(D, axis=1)
    eps_grid = sorted(set(np.round(np.percentile(pdist(Xv), [10, 20, 30, 40, 50, 60, 70]), 3)))
    min_samples_grid = [2, 3, 4, 5]

    rows = []
    labels_store = {}
    for eps, min_samples in itertools.product(eps_grid, min_samples_grid):
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(Xv)
        labels = db.labels_
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        noise_fraction = float((labels == -1).mean())
        row = {
            "input": input_label, "eps": eps, "min_samples": min_samples,
            "n_clusters_found": n_clusters, "noise_fraction": noise_fraction,
        }
        if n_clusters >= 2:
            non_noise = labels != -1
            if non_noise.sum() >= 2 and len(set(labels[non_noise])) >= 2:
                row["silhouette"] = silhouette_score(Xv[non_noise], labels[non_noise])
                row["davies_bouldin"] = davies_bouldin_score(Xv[non_noise], labels[non_noise])
                row["calinski_harabasz"] = calinski_harabasz_score(Xv[non_noise], labels[non_noise])
        rows.append(row)
        labels_store[(eps, min_samples)] = pd.Series(labels, index=ids)

    return pd.DataFrame(rows), labels_store


# ------------------------------------------------------- method agreement --
def compute_method_agreement(labels_store: dict, k_range, input_label: str):
    methods = ["kmeans", "agglomerative_ward", "gmm", "spectral"]
    rows = []
    for k in k_range:
        available = [m for m in methods if (m, k) in labels_store]
        for m1, m2 in itertools.combinations(available, 2):
            l1 = labels_store[(m1, k)]
            l2 = labels_store[(m2, k)]
            common = l1.index.intersection(l2.index)
            ari = adjusted_rand_score(l1.loc[common], l2.loc[common])
            ri = rand_score(l1.loc[common], l2.loc[common])
            rows.append({"input": input_label, "k": k, "method_a": m1, "method_b": m2,
                         "adjusted_rand_index": ari, "rand_index": ri})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ main ---
def main():
    primary_X = load_primary_input()
    n = len(primary_X)
    k_range = list(range(2, min(10, n - 1) + 1))
    print(f"Primary clustering input: PCA scores (>=80% variance), {primary_X.shape[0]} receptors x "
          f"{primary_X.shape[1]} PCs. k range = {k_range}")

    all_validation = []
    all_gap = []
    all_agreement = []
    all_dbscan = []
    gap_chosen = {}
    labels_by_input = {}

    print("\n--- PRIMARY: PCA-space clustering ---")
    val_pca, gap_pca, gap_k_pca, labels_pca = run_kbased_battery(primary_X, k_range, "pca_primary")
    all_validation.append(val_pca)
    all_gap.append(gap_pca)
    gap_chosen["pca_primary"] = gap_k_pca
    labels_by_input["pca_primary"] = labels_pca
    agree_pca = compute_method_agreement(labels_pca, k_range, "pca_primary")
    all_agreement.append(agree_pca)
    dbscan_pca, dbscan_labels_pca = run_dbscan_grid(primary_X, "pca_primary")
    all_dbscan.append(dbscan_pca)

    sensitivity_inputs = load_sensitivity_inputs()
    for label, X in sensitivity_inputs.items():
        n_s = len(X)
        k_range_s = list(range(2, min(10, n_s - 1) + 1))
        print(f"\n--- SENSITIVITY: {label} ({n_s} receptors x {X.shape[1]} features), k range = {k_range_s} ---")
        val, gap, gap_k, labels = run_kbased_battery(X, k_range_s, label)
        all_validation.append(val)
        all_gap.append(gap)
        gap_chosen[label] = gap_k
        labels_by_input[label] = labels
        all_agreement.append(compute_method_agreement(labels, k_range_s, label))
        dbscan_df, dbscan_labels = run_dbscan_grid(X, label)
        all_dbscan.append(dbscan_df)

    validation = pd.concat(all_validation, ignore_index=True)
    gap_all = pd.concat(all_gap, ignore_index=True)
    agreement = pd.concat(all_agreement, ignore_index=True)
    dbscan_all = pd.concat(all_dbscan, ignore_index=True)

    validation_out = TABLES / "clustering_validation.csv"
    validation.to_csv(validation_out, index=False)
    gap_out = TABLES / "clustering_gap_statistic.csv"
    gap_all.to_csv(gap_out, index=False)
    agreement_out = TABLES / "clustering_method_agreement.csv"
    agreement.to_csv(agreement_out, index=False)
    dbscan_out = TABLES / "clustering_dbscan_grid.csv"
    dbscan_all.to_csv(dbscan_out, index=False)

    print(f"\nWrote {validation_out.relative_to(PROJECT_ROOT)} ({len(validation)} rows)")
    print(f"Wrote {gap_out.relative_to(PROJECT_ROOT)} ({len(gap_all)} rows)")
    print(f"Wrote {agreement_out.relative_to(PROJECT_ROOT)} ({len(agreement)} rows)")
    print(f"Wrote {dbscan_out.relative_to(PROJECT_ROOT)} ({len(dbscan_all)} rows)")

    print("\nGap-statistic chosen k per input (Tibshirani rule):")
    for label, k in gap_chosen.items():
        print(f"  {label}: k = {k}")

    # ---- sensitivity: cross-input ARI at matched k for the PRIMARY method (kmeans) ----
    sens_rows = []
    primary_ids = set(primary_X.index)
    for label, labels_dict in labels_by_input.items():
        if label == "pca_primary":
            continue
        for k in k_range:
            if ("kmeans", k) not in labels_dict or ("kmeans", k) not in labels_by_input["pca_primary"]:
                continue
            l_sens = labels_dict[("kmeans", k)]
            l_prim = labels_by_input["pca_primary"][("kmeans", k)]
            common = sorted(set(l_sens.index) & set(l_prim.index) & primary_ids)
            if len(common) < max(4, k + 1):
                continue
            ari = adjusted_rand_score(l_prim.loc[common], l_sens.loc[common])
            sens_rows.append({
                "comparison": f"pca_primary_vs_{label}", "method": "kmeans", "k": k,
                "n_common_receptors": len(common), "adjusted_rand_index": ari,
            })
    sensitivity_df = pd.DataFrame(sens_rows)
    sensitivity_out = TABLES / "clustering_sensitivity.csv"
    sensitivity_df.to_csv(sensitivity_out, index=False)
    print(f"\nWrote {sensitivity_out.relative_to(PROJECT_ROOT)} ({len(sensitivity_df)} rows)")
    if len(sensitivity_df):
        mean_ari_by_comparison = sensitivity_df.groupby("comparison")["adjusted_rand_index"].mean()
        print(mean_ari_by_comparison.to_string())

    # ---- descriptive continuum-vs-cluster signal from the numbers actually computed ----
    primary_silhouettes = validation[(validation["input"] == "pca_primary")].dropna(subset=["silhouette"])
    max_sil = primary_silhouettes["silhouette"].max() if len(primary_silhouettes) else np.nan
    mean_ari_across_methods = agreement[agreement["input"] == "pca_primary"]["adjusted_rand_index"].mean()
    print(f"\nDescriptive signal (from computed numbers, not asserted): "
          f"max silhouette (PCA-primary, any method/k) = {max_sil:.3f}; "
          f"mean cross-method ARI (PCA-primary, all k) = {mean_ari_across_methods:.3f}")

    # ---- figures: k-selection curves for the primary input ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for method in ["kmeans", "agglomerative_ward", "gmm", "spectral"]:
        sub = val_pca[val_pca["method"] == method]
        axes[0].plot(sub["k"], sub["silhouette"], "o-", label=method)
        axes[1].plot(sub["k"], sub["davies_bouldin"], "o-", label=method)
        axes[2].plot(sub["k"], sub["calinski_harabasz"], "o-", label=method)
    axes[0].set_title("Silhouette (higher better)"); axes[0].set_xlabel("k")
    axes[1].set_title("Davies-Bouldin (lower better)"); axes[1].set_xlabel("k")
    axes[2].set_title("Calinski-Harabasz (higher better)"); axes[2].set_xlabel("k")
    axes[0].legend(fontsize=7)
    fig.suptitle("Clustering validation metrics vs k -- PCA-primary input")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "clustering_validation_metrics.png", dpi=150)
    plt.close(fig)

    # ---- cluster-colored UMAP figure (Part 7 deferred requirement) ----
    best_row = primary_silhouettes.sort_values("silhouette", ascending=False).iloc[0] if len(primary_silhouettes) else None
    if best_row is not None:
        best_method, best_k = best_row["method"], int(best_row["k"])
        best_labels = labels_by_input["pca_primary"][(best_method, best_k)]
        umap_path = ANALYSIS / "umap_default_embedding_pca_representation.csv"
        if umap_path.exists():
            umap_df = pd.read_csv(umap_path).set_index("uniprot_id")
            umap_df = umap_df.join(best_labels.rename("cluster"))
            fig, ax = plt.subplots(figsize=(7, 6))
            for c in sorted(umap_df["cluster"].dropna().unique()):
                sub = umap_df[umap_df["cluster"] == c]
                ax.scatter(sub["UMAP1"], sub["UMAP2"], label=f"cluster {int(c)}", s=50)
            ax.set_title(f"UMAP (PCA representation), colored by best-silhouette solution\n"
                         f"({best_method}, k={best_k}, silhouette={best_row['silhouette']:.3f}) "
                         f"-- descriptive only")
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(FIG_DIR / "umap_pca_representation_by_cluster.png", dpi=150)
            plt.close(fig)
        print(f"\nBest single silhouette across the whole PCA-primary battery: "
              f"{best_method}, k={best_k}, silhouette={best_row['silhouette']:.3f} "
              f"(reported for completeness -- NOT selected as 'the' answer; see Part 11)")

    summary = {
        "k_range": k_range,
        "gap_chosen_k": gap_chosen,
        "max_silhouette_pca_primary": float(max_sil) if not np.isnan(max_sil) else None,
        "mean_cross_method_ari_pca_primary": float(mean_ari_across_methods) if not np.isnan(mean_ari_across_methods) else None,
    }
    (ANALYSIS / "clustering_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nWrote {(ANALYSIS / 'clustering_summary.json').relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
