"""
Part 12 of the subfamily flexibility stage: reproducible figures under
reports/figures/subfamily_stats/. All cross-group inferential plots use
RECEPTOR points (never the 1789 PDB entries). The PCA plot reuses the
already-fit modeling-stage PCA scores (data/analysis/pca_receptor_scores.csv)
purely for annotation/interpretation -- it is NOT refit here.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "subfamily_stats"
FIG_DIR.mkdir(parents=True, exist_ok=True)

WELL_SUPPORTED_MIN_N = 10


def fig_atlas_ci(atlas: pd.DataFrame):
    well = atlas[atlas["well_supported"]].sort_values("median_receptor_z")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(well))))
    y = np.arange(len(well))
    ax.errorbar(
        well["median_receptor_z"], y,
        xerr=[well["median_receptor_z"] - well["bootstrap_95ci_median_low"],
              well["bootstrap_95ci_median_high"] - well["median_receptor_z"]],
        fmt="o", color="#1f77b4", ecolor="#888888", capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels([f"{s} (n={n})" for s, n in zip(well["structural_segment"], well["n_receptors"])], fontsize=8)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Median receptor Z (within-receptor relative flexibility)")
    ax.set_title(f"Segment flexibility atlas: median Z + 95% bootstrap CI\n"
                 f"(well-supported segments only, n_receptors >= {WELL_SUPPORTED_MIN_N})")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "segment_atlas_median_ci.png", dpi=150)
    plt.close(fig)


def fig_receptor_segment_heatmap(profiles: pd.DataFrame, atlas: pd.DataFrame):
    well_segments = atlas.loc[atlas["well_supported"], "structural_segment"].tolist()
    pivot = profiles.pivot(index="uniprot_id", columns="structural_segment", values="median_of_entry_median_z")
    pivot = pivot[well_segments]
    # order receptors by nr_code for readability
    meta = pd.read_csv(ANALYSIS / "subfamily_receptor_metadata.csv").set_index("uniprot_id")
    pivot = pivot.reindex(meta.sort_values("nr_code").index).dropna(how="all")
    labels = meta.loc[pivot.index, "nr_code"]

    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(well_segments)), max(6, 0.22 * len(pivot))))
    masked = np.ma.masked_invalid(pivot.values)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="#e0e0e0")  # missing shown explicitly, not blended into the color scale
    vmax = np.nanmax(np.abs(pivot.values))
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(well_segments)))
    ax.set_xticklabels(well_segments, rotation=90, fontsize=7)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_title("Receptor x segment median Z (well-supported segments)\ngray = missing (never imputed for this view)")
    fig.colorbar(im, ax=ax, label="median receptor Z", shrink=0.6)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "receptor_segment_heatmap.png", dpi=150)
    plt.close(fig)


def fig_violin_by_group(profiles: pd.DataFrame, meta: pd.DataFrame, kruskal: pd.DataFrame):
    tested = kruskal[(kruskal["grouping_scheme"] == "primary_project_group") &
                      (kruskal["test_status"] == "TESTED")].sort_values("raw_p")
    segs = tested["structural_segment"].tolist()
    n = len(segs)
    if n == 0:
        return
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.2 * nrows), squeeze=False)
    profiles_core = profiles[["uniprot_id", "structural_segment", "median_of_entry_median_z"]]
    merged = profiles_core.merge(meta[["uniprot_id", "group"]], on="uniprot_id", how="left")
    groups = sorted(merged["group"].dropna().unique())
    for i, seg in enumerate(segs):
        ax = axes[i // ncols][i % ncols]
        seg_df = merged[merged["structural_segment"] == seg]
        data = [seg_df.loc[seg_df["group"] == g, "median_of_entry_median_z"].dropna().values for g in groups]
        parts = ax.violinplot([d for d in data if len(d) > 0], showmedians=True)
        present_groups = [g for g, d in zip(groups, data) if len(d) > 0]
        for xi, d in enumerate(data, start=1):
            if len(d) == 0:
                continue
            jitter = np.random.default_rng(0).normal(0, 0.04, size=len(d))
            ax.scatter(np.full(len(d), xi) + jitter, d, color="black", s=10, alpha=0.6, zorder=3)
        ax.set_xticks(range(1, len(present_groups) + 1))
        ax.set_xticklabels([g[:12] for g in present_groups], rotation=30, fontsize=6, ha="right")
        row = kruskal[(kruskal["grouping_scheme"] == "primary_project_group") & (kruskal["structural_segment"] == seg)].iloc[0]
        ax.set_title(f"{seg}\nq={row['fdr_q']:.2f}", fontsize=8)
        ax.axhline(0, color="lightgray", linewidth=0.5)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle("Per-segment receptor Z by project group (PRIMARY grouping) -- "
                 "individual receptor points shown; none FDR-significant", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "segment_by_group_violin.png", dpi=150)
    plt.close(fig)


def fig_effect_sizes(kruskal: pd.DataFrame):
    primary = kruskal[(kruskal["grouping_scheme"] == "primary_project_group") &
                       (kruskal["test_status"] == "TESTED")].sort_values("epsilon_squared", ascending=True)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(primary))))
    y = np.arange(len(primary))
    ax.errorbar(
        primary["epsilon_squared"], y,
        xerr=[primary["epsilon_squared"] - primary["epsilon_squared_ci_low"],
              primary["epsilon_squared_ci_high"] - primary["epsilon_squared"]],
        fmt="s", color="#d62728", ecolor="#999999", capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(primary["structural_segment"], fontsize=8)
    colors = ["green" if sig else "black" for sig in primary["significant_fdr_05"]]
    for yi, c in zip(y, colors):
        ax.get_yticklabels()[yi].set_color(c)
    ax.set_xlabel("epsilon-squared (Kruskal-Wallis rank-based effect size)")
    ax.set_title("Omnibus effect size per segment, PRIMARY grouping (project functional group)\n"
                 "(none reach FDR significance; green would mark significant)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "effect_size_by_segment.png", dpi=150)
    plt.close(fig)


def fig_pca_annotated(kruskal: pd.DataFrame):
    pca_path = ANALYSIS / "pca_receptor_scores.csv"
    if not pca_path.exists():
        print("modeling-stage PCA scores not found -- skipping annotated PCA figure")
        return
    scores = pd.read_csv(pca_path)  # already-fit scores, NOT refit here
    fig, ax = plt.subplots(figsize=(7, 6))
    groups = scores["group"].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))
    for g, c in zip(groups, colors):
        sub = scores[scores["group"] == g]
        ax.scatter(sub["PC1"], sub["PC2"], color=c, label=g, s=50)
    for _, row in scores.iterrows():
        ax.annotate(row["nr_code"], (row["PC1"], row["PC2"]), fontsize=6, xytext=(3, 3), textcoords="offset points")
    nominal = kruskal[(kruskal["raw_p"] < 0.05)]
    note = ("No segment reached FDR significance in either grouping scheme.\n"
            f"Nominal (uncorrected) hits: {', '.join(nominal['structural_segment'].unique())} "
            "(not robust under BH-FDR or, for one, leave-one-out)")
    ax.set_title("PCA (modeling stage, reused not refit) annotated with subfamily-test outcome")
    ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    ax.legend(fontsize=8)
    ax.axhline(0, color="lightgray", linewidth=0.5)
    ax.axvline(0, color="lightgray", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pca_annotated_with_subfamily_results.png", dpi=150)
    plt.close(fig)


def main():
    profiles = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    meta = pd.read_csv(ANALYSIS / "subfamily_receptor_metadata.csv")
    atlas = pd.read_csv(ANALYSIS / "superfamily_segment_flexibility_atlas.csv")
    kruskal = pd.read_csv(TABLES / "subfamily_segment_kruskal.csv")

    fig_atlas_ci(atlas)
    print(f"Wrote {(FIG_DIR / 'segment_atlas_median_ci.png').relative_to(PROJECT_ROOT)}")

    fig_receptor_segment_heatmap(profiles, atlas)
    print(f"Wrote {(FIG_DIR / 'receptor_segment_heatmap.png').relative_to(PROJECT_ROOT)}")

    fig_violin_by_group(profiles, meta, kruskal)
    print(f"Wrote {(FIG_DIR / 'segment_by_group_violin.png').relative_to(PROJECT_ROOT)}")

    fig_effect_sizes(kruskal)
    print(f"Wrote {(FIG_DIR / 'effect_size_by_segment.png').relative_to(PROJECT_ROOT)}")

    fig_pca_annotated(kruskal)
    print(f"Wrote {(FIG_DIR / 'pca_annotated_with_subfamily_results.png').relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
