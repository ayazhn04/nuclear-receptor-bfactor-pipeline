"""
Structural-segment alignment, Part 12: analysis/QC figures (not final
publication figures). Reproducible: reads only already-generated
data/analysis/ tables; no data is computed here that isn't already saved.

Writes to reports/figures/alignment_qc/:
    01_segment_coverage_heatmap.png
    02_segment_pdb_support_by_receptor.png
    03_segment_missingness_summary.png
    04_example_lbd_schematic_<receptor>.png  (VDR, PPARGAMMA, ESR1, GR, NR4A1)
    05_segment_median_z_distributions.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "alignment_qc"

CANONICAL_VOCAB = (
    ["H1", "loop_H1_H2", "H2", "loop_H2_H3", "H3", "loop_H3_H4", "H4", "loop_H4_H5", "H5",
     "loop_H5_beta_sheet", "beta_sheet", "loop_beta_sheet_H6", "H6", "loop_H6_H7", "H7",
     "loop_H7_H8", "H8", "loop_H8_H9", "H9", "loop_H9_H10", "H10", "loop_H10_H11", "H11",
     "loop_H11_H12", "H12"]
)

ANCHORS = {"P11473": "VDR", "P37231": "PPARGAMMA", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}


def ordered_cols(cols) -> list[str]:
    rank = {s: i for i, s in enumerate(CANONICAL_VOCAB)}
    return sorted(cols, key=lambda c: rank.get(c, 99))


def fig1_coverage_heatmap(cov: pd.DataFrame, meta: dict) -> None:
    cov = cov.set_index("uniprot_id")
    cols = ordered_cols(cov.columns)
    cov = cov[cols]
    labels = [meta.get(u, u) for u in cov.index]
    fig, ax = plt.subplots(figsize=(max(8, len(cols) * 0.5), max(6, len(cov) * 0.28)))
    im = ax.imshow(cov.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=6)
    ax.set_title("Median per-segment residue coverage by receptor (QC; metadata, not a flexibility feature)")
    fig.colorbar(im, ax=ax, label="median segment_coverage")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "01_segment_coverage_heatmap.png", dpi=150)
    plt.close(fig)


def fig2_pdb_support(rec_prof: pd.DataFrame, meta: dict) -> None:
    p = rec_prof.pivot(index="uniprot_id", columns="structural_segment", values="n_pdb_entries_with_segment")
    cols = ordered_cols(p.columns)
    p = p[cols].fillna(0)
    labels = [meta.get(u, u) for u in p.index]
    fig, ax = plt.subplots(figsize=(max(8, len(cols) * 0.5), max(6, len(p) * 0.28)))
    im = ax.imshow(p.values, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=6)
    ax.set_title("Number of PDB entries supporting each segment, by receptor")
    fig.colorbar(im, ax=ax, label="n_pdb_entries_with_segment")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "02_segment_pdb_support_by_receptor.png", dpi=150)
    plt.close(fig)


def fig3_missingness_summary(matrix: pd.DataFrame, lbd_ref: pd.DataFrame) -> None:
    idx = matrix.set_index("uniprot_id")
    total = len(lbd_ref)
    cols = ordered_cols(idx.columns)
    represented = idx[cols].notna().sum()
    fig, ax = plt.subplots(figsize=(max(8, len(cols) * 0.4), 4))
    ax.bar(range(len(cols)), represented.values, color="#4c72b0")
    ax.axhline(total, color="grey", linestyle="--", linewidth=1, label=f"all {total} receptors")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=7)
    ax.set_ylabel("receptors represented")
    ax.set_title("Segment-feature-matrix missingness summary (out of 48 project receptors)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_segment_missingness_summary.png", dpi=150)
    plt.close(fig)


def fig4_schematics(seg_map: pd.DataFrame, lbd_ref: pd.DataFrame) -> None:
    class_color = {"H1": "#d62728", "H12": "#9467bd", "beta_sheet": "#2ca02c"}

    def color_for(seg: str, cls: str) -> str:
        if seg in class_color:
            return class_color[seg]
        if cls == "HELIX":
            return "#1f77b4"
        if cls == "LOOP":
            return "#c7c7c7"
        return "#eeeeee"

    for uid, name in ANCHORS.items():
        g = seg_map[seg_map["uniprot_id"] == uid].sort_values("canonical_uniprot_position")
        if len(g) == 0:
            continue
        row = lbd_ref[lbd_ref["uniprot_id"] == uid].iloc[0]
        fig, ax = plt.subplots(figsize=(12, 1.8))
        for _, r in g.iterrows():
            ax.add_patch(plt.Rectangle((r["canonical_uniprot_position"], 0), 1, 1,
                                        color=color_for(r["structural_segment"], r["structural_class"]), linewidth=0))
        for seg in g[g["structural_segment"] != ""]["structural_segment"].unique():
            gg = g[g["structural_segment"] == seg]
            mid = (gg["canonical_uniprot_position"].min() + gg["canonical_uniprot_position"].max()) / 2
            ax.text(mid, 1.15, seg, ha="center", va="bottom", fontsize=6, rotation=90 if len(seg) > 3 else 0)
        ax.set_xlim(row["lbd_start"] - 2, row["lbd_end"] + 2)
        ax.set_ylim(0, 2.2)
        ax.set_yticks([])
        ax.set_xlabel("canonical UniProt position")
        ax.set_title(f"{name} ({uid}) standardized LBD segment map (red=H1, blue=other helix, green=beta_sheet, purple=H12, grey=loop)", fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"04_example_lbd_schematic_{name}.png", dpi=150)
        plt.close(fig)


def fig5_z_distributions(rec_prof: pd.DataFrame) -> None:
    segs = ["H1", "H12", "beta_sheet"]
    data = [rec_prof.loc[rec_prof["structural_segment"] == s, "median_of_entry_median_z"].dropna() for s in segs]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.boxplot(data, tick_labels=segs, showmeans=True)
    ax.axhline(0, color="grey", linestyle=":", linewidth=1)
    ax.set_ylabel("receptor median-of-PDB-median normalized Z")
    ax.set_title("Distribution of segment median Z across receptors (well-supported segments only)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "05_segment_median_z_distributions.png", dpi=150)
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    cov = pd.read_parquet(ANALYSIS / "receptor_segment_coverage_matrix.parquet")
    matrix = pd.read_parquet(ANALYSIS / "receptor_segment_feature_matrix_median.parquet")
    rec_prof = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    lbd_ref = pd.read_csv(MANIFESTS / "nr_lbd_reference.csv")
    meta = dict(zip(lbd_ref["uniprot_id"], lbd_ref["common_name"]))

    fig1_coverage_heatmap(cov, meta)
    fig2_pdb_support(rec_prof, meta)
    fig3_missingness_summary(matrix, lbd_ref)
    fig4_schematics(seg_map, lbd_ref)
    fig5_z_distributions(rec_prof)
    print("Wrote figures to", FIG_DIR.relative_to(PROJECT_ROOT))
    for p in sorted(FIG_DIR.glob("*.png")):
        print(" ", p.name)


if __name__ == "__main__":
    main()
