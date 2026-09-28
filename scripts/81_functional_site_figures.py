"""
Part 16 of the current functional-site flexibility stage: reproducible
figures under reports/figures/functional_sites_current/. Every figure
title/caption is explicit that these are CURRENT full-dataset results,
distinct from the historical Role D analysis (reports/ROLE_D_PROGRESS_
REPORT.md, reports/figures/functional_site_qc/ if present) -- never
implied to be the same numbers.
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
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "functional_sites_current"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SUBTITLE = "CURRENT full-dataset analysis (this stage) -- distinct from the historical Role D results"


def fig_receptor_site_heatmap(receptor_site: pd.DataFrame):
    pivot = receptor_site.pivot(index="nr_code", columns="functional_site", values="median_site_z")
    fig, ax = plt.subplots(figsize=(6, max(6, 0.22 * len(pivot))))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="#e0e0e0")
    vmax = np.nanmax(np.abs(pivot.values))
    im = ax.imshow(np.ma.masked_invalid(pivot.values), aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(pivot.shape[1])); ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0])); ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_title(f"Receptor x functional-site median Z\n{SUBTITLE}\ngray = site UNRESOLVED for that receptor (never zero-filled)", fontsize=9)
    fig.colorbar(im, ax=ax, label="median site Z", shrink=0.6)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "receptor_site_median_z_heatmap.png", dpi=150)
    plt.close(fig)


def fig_site_background_delta_heatmap(receptor_site: pd.DataFrame):
    pivot = receptor_site.pivot(index="nr_code", columns="functional_site", values="median_delta_site_vs_background")
    fig, ax = plt.subplots(figsize=(6, max(6, 0.22 * len(pivot))))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="#e0e0e0")
    vmax = np.nanmax(np.abs(pivot.values))
    im = ax.imshow(np.ma.masked_invalid(pivot.values), aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(pivot.shape[1])); ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0])); ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_title(f"Receptor x site delta (site - background median Z)\n{SUBTITLE}", fontsize=9)
    fig.colorbar(im, ax=ax, label="delta (site - background)", shrink=0.6)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "receptor_site_background_delta_heatmap.png", dpi=150)
    plt.close(fig)


def fig_site_class_forest(superfamily: pd.DataFrame):
    df = superfamily.sort_values("median_delta_vs_background")
    fig, ax = plt.subplots(figsize=(9, 3.5))
    y = np.arange(len(df))
    colors = ["#2ca02c" if s else "#888888" for s in df["significant_fdr_05"]]
    ax.errorbar(df["median_delta_vs_background"], y,
                xerr=[df["median_delta_vs_background"] - df["bootstrap_ci_low"],
                      df["bootstrap_ci_high"] - df["median_delta_vs_background"]],
                fmt="none", ecolor="#cccccc", capsize=3)
    ax.scatter(df["median_delta_vs_background"], y, color=colors, s=60, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(df["functional_site"])
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_xlabel("median receptor-level delta (site - background) + bootstrap 95% CI")
    ax.set_title(f"Superfamily site-class comparison\n{SUBTITLE}\n(green = FDR-significant vs zero)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "site_class_forest_plot.png", dpi=150)
    plt.close(fig)


def fig_site_distribution(receptor_site: pd.DataFrame, site: str, fname: str):
    g = receptor_site[receptor_site["functional_site"] == site]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(g["median_site_z"], bins=15, color="#1f77b4", edgecolor="white")
    ax.axvline(g["median_site_z"].median(), color="red", linestyle="--", label="median")
    ax.set_xlabel("receptor median site Z"); ax.set_ylabel("n receptors")
    ax.set_title(f"{site} flexibility distribution across receptors (n={len(g)})\n{SUBTITLE}", fontsize=9)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / fname, dpi=150)
    plt.close(fig)


def fig_site_vs_background_paired(pdb_site: pd.DataFrame):
    counts = pdb_site.groupby(["uniprot_id", "functional_site"]).size()
    well = counts[counts >= 15].reset_index()
    well.columns = ["uniprot_id", "functional_site", "n"]
    if len(well) == 0:
        print("No receptor x site pair with >=15 PDB entries -- skipping paired plot")
        return
    top = well.sort_values("n", ascending=False).iloc[0]
    uid, site = top["uniprot_id"], top["functional_site"]
    g = pdb_site[(pdb_site["uniprot_id"] == uid) & (pdb_site["functional_site"] == site)]
    nr_code = g["nr_code"].iloc[0]

    fig, ax = plt.subplots(figsize=(5, 5))
    for _, r in g.iterrows():
        ax.plot([0, 1], [r["comparison_non_site_median_z"], r["median_site_z"]], color="gray", alpha=0.4, linewidth=0.8)
    ax.scatter(np.zeros(len(g)), g["comparison_non_site_median_z"], color="#888888", zorder=3, label="background")
    ax.scatter(np.ones(len(g)), g["median_site_z"], color="#d62728", zorder=3, label=site)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["background", site])
    ax.set_ylabel("median Z per PDB entry")
    ax.set_title(f"{nr_code} ({uid}): {site} vs background, {len(g)} paired PDB entries\n{SUBTITLE}", fontsize=9)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"paired_site_vs_background_{nr_code}_{site}.png", dpi=150)
    plt.close(fig)


def fig_apo_holo_site_schematic(overlap: pd.DataFrame):
    sig = overlap[overlap["apo_holo_significant_within_receptor"]]
    if len(sig) == 0:
        print("No significant apo/holo results -- skipping overlap schematic")
        return
    fig, ax = plt.subplots(figsize=(8, max(3, 0.4 * len(sig))))
    y = np.arange(len(sig))
    colors = ["#2ca02c" if o else "#888888" for o in sig["functional_site_overlap"]]
    ax.barh(y, sig["apo_holo_delta_median_z"], color=colors)
    labels = [f"{n}/{s}" + (f" -> {site}" if site else " (no site overlap)")
              for n, s, site in zip(sig["nr_code"], sig["apo_holo_segment"], sig["site_name"])]
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_xlabel("apo/holo delta median Z (HOLO - APO)")
    ax.set_title(f"Apo/holo significant segments x functional-site overlap (descriptive only)\n{SUBTITLE}", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "apo_holo_site_overlap_schematic.png", dpi=150)
    plt.close(fig)


def fig_anchor_schematics(definitions: pd.DataFrame):
    anchors = {"P11473": "VDR", "P37231": "PPARgamma", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    colors = {"AF2_H12": "#9467bd", "LIGAND_BINDING_POCKET": "#1f77b4",
              "COACTIVATOR_GROOVE": "#ff7f0e", "DIMERIZATION_INTERFACE": "#2ca02c"}
    for uid, name in anchors.items():
        rec_def = definitions[definitions["uniprot_id"] == uid]
        if len(rec_def) == 0:
            continue
        rec_seg = seg_map[(seg_map["uniprot_id"] == uid) & (seg_map["structural_segment"] != "")]
        if len(rec_seg) == 0:
            continue
        lo, hi = rec_seg["canonical_uniprot_position"].min(), rec_seg["canonical_uniprot_position"].max()

        fig, ax = plt.subplots(figsize=(12, 2.5))
        ax.set_xlim(lo - 5, hi + 5); ax.set_ylim(0, 1)
        for i, (site, c) in enumerate(colors.items()):
            positions = sorted(rec_def.loc[rec_def["functional_site"] == site, "canonical_uniprot_position"])
            y = 0.15 + i * 0.2
            for p in positions:
                ax.plot([p, p], [y, y + 0.15], color=c, linewidth=1)
            ax.text(hi + 8, y + 0.07, site, fontsize=7, va="center", color=c)
        ax.set_yticks([])
        ax.set_xlabel("canonical UniProt position")
        ax.set_title(f"{name} ({uid}): current functional-site map\n{SUBTITLE}", fontsize=9)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"anchor_schematic_{name}.png", dpi=150)
        plt.close(fig)


def main():
    receptor_site = pd.read_parquet(ANALYSIS / "receptor_functional_site_flexibility.parquet")
    pdb_site = pd.read_parquet(ANALYSIS / "pdb_functional_site_flexibility.parquet")
    superfamily = pd.read_csv(TABLES / "current_functional_site_superfamily_tests.csv")
    overlap = pd.read_csv(TABLES / "apo_holo_functional_site_overlap.csv", keep_default_na=False)
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")

    fig_receptor_site_heatmap(receptor_site)
    fig_site_background_delta_heatmap(receptor_site)
    fig_site_class_forest(superfamily)
    fig_site_distribution(receptor_site, "AF2_H12", "af2_h12_distribution.png")
    fig_site_distribution(receptor_site, "LIGAND_BINDING_POCKET", "lbp_distribution.png")
    fig_site_vs_background_paired(pdb_site)
    fig_apo_holo_site_schematic(overlap)
    fig_anchor_schematics(definitions)

    print(f"Wrote figures to {FIG_DIR.relative_to(PROJECT_ROOT)}/:")
    for f in sorted(FIG_DIR.glob("*.png")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
