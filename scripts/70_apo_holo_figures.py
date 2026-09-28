"""
Part 14 of the APO vs HOLO flexibility stage: reproducible figures under
reports/figures/apo_holo/. VDR is never plotted as though it has an
APO comparison (it is excluded from every APO/HOLO figure -- it has 0 APO
entries).
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
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "apo_holo"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def fig_counts_by_receptor(eligibility: pd.DataFrame):
    df = eligibility.sort_values("n_total_coordinate_pdbs", ascending=True)
    fig, ax = plt.subplots(figsize=(8, max(6, 0.28 * len(df))))
    y = np.arange(len(df))
    ax.barh(y, df["n_APO"], color="#1f77b4", label="APO")
    ax.barh(y, df["n_HOLO"], left=df["n_APO"], color="#d62728", label="HOLO")
    ax.barh(y, df["n_AMBIGUOUS"], left=df["n_APO"] + df["n_HOLO"], color="#999999", label="AMBIGUOUS")
    colors = ["green" if p else "black" for p in df["eligible_primary_5_5"]]
    ax.set_yticks(y)
    ax.set_yticklabels(df["nr_code"], fontsize=7)
    for tick, c in zip(ax.get_yticklabels(), colors):
        tick.set_color(c)
    ax.set_xlabel("PDB entries")
    ax.set_title("APO/HOLO/AMBIGUOUS PDB counts by receptor\n(green label = PRIMARY-eligible >=5/5; VDR = 0 APO, 48 HOLO)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "apo_holo_counts_by_receptor.png", dpi=150)
    plt.close(fig)


def fig_delta_heatmap(effects: pd.DataFrame):
    pivot = effects.pivot_table(index="nr_code", columns="structural_segment", values="delta_median_z")
    fig, ax = plt.subplots(figsize=(max(8, 0.5 * pivot.shape[1]), max(4, 0.4 * pivot.shape[0])))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="#e0e0e0")
    vmax = np.nanmax(np.abs(pivot.values))
    masked = np.ma.masked_invalid(pivot.values)
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title("Receptor x segment delta (HOLO - APO median Z)\ngray = not tested (insufficient segment-specific coverage, never imputed)")
    fig.colorbar(im, ax=ax, label="delta median Z (HOLO - APO)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "receptor_segment_delta_heatmap.png", dpi=150)
    plt.close(fig)


def fig_forest(effects: pd.DataFrame):
    df = effects.sort_values("delta_median_z")
    fig, ax = plt.subplots(figsize=(8, max(6, 0.2 * len(df))))
    y = np.arange(len(df))
    colors = ["#2ca02c" if s else "#888888" for s in df["significant_within_receptor"]]
    ax.errorbar(df["delta_median_z"], y,
                xerr=[df["delta_median_z"] - df["bootstrap_ci_low"], df["bootstrap_ci_high"] - df["delta_median_z"]],
                fmt="none", ecolor="#cccccc", capsize=2, zorder=1)
    ax.scatter(df["delta_median_z"], y, color=colors, s=25, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{n} / {s}" for n, s in zip(df["nr_code"], df["structural_segment"])], fontsize=5)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_xlabel("delta median Z (HOLO - APO) + bootstrap 95% CI")
    ax.set_title("All tested receptor x segment effects\n(green = within-receptor FDR significant)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "effect_forest_plot.png", dpi=150)
    plt.close(fig)


def fig_per_receptor_segments(features: pd.DataFrame, eligibility: pd.DataFrame):
    eligible = eligibility[eligibility["eligible_exploratory_3_3"]]["uniprot_id"].tolist()
    for uid in eligible:
        sub = features[(features["uniprot_id"] == uid) & (features["ligand_state"].isin(["APO", "HOLO"]))]
        segs = sorted(sub["structural_segment"].unique())
        nr_code = sub["nr_code"].iloc[0]
        n = len(segs)
        if n == 0:
            continue
        ncols = 4
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.8 * nrows), squeeze=False)
        for i, seg in enumerate(segs):
            ax = axes[i // ncols][i % ncols]
            seg_df = sub[sub["structural_segment"] == seg]
            apo = seg_df.loc[seg_df["ligand_state"] == "APO", "median_z"].dropna().values
            holo = seg_df.loc[seg_df["ligand_state"] == "HOLO", "median_z"].dropna().values
            data = [d for d in (apo, holo) if len(d) > 0]
            labels = [lbl for lbl, d in zip(["APO", "HOLO"], (apo, holo)) if len(d) > 0]
            if data:
                ax.violinplot(data, showmedians=True)
                for xi, d in enumerate(data, start=1):
                    jitter = np.random.default_rng(0).normal(0, 0.04, size=len(d))
                    ax.scatter(np.full(len(d), xi) + jitter, d, color="black", s=8, alpha=0.6)
                ax.set_xticks(range(1, len(labels) + 1))
                ax.set_xticklabels([f"{l} (n={len(d)})" for l, d in zip(labels, data)], fontsize=6)
            ax.set_title(seg, fontsize=8)
            ax.axhline(0, color="lightgray", linewidth=0.5)
        for j in range(n, nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.suptitle(f"{nr_code} ({uid}): APO vs HOLO per segment", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"per_receptor_{nr_code}.png", dpi=150)
        plt.close(fig)


def fig_residue_profile(residue_effects: pd.DataFrame):
    for uid in residue_effects["uniprot_id"].unique():
        sub = residue_effects[residue_effects["uniprot_id"] == uid].sort_values("canonical_uniprot_position")
        nr_code = sub["nr_code"].iloc[0]
        fig, ax = plt.subplots(figsize=(10, 3.5))
        colors = ["#2ca02c" if s else "#999999" for s in sub["significant_within_receptor"]]
        ax.scatter(sub["canonical_uniprot_position"], sub["delta_median_z"], c=colors, s=12)
        ax.axhline(0, color="black", linewidth=0.6)
        ax.set_xlabel("canonical UniProt position")
        ax.set_ylabel("delta median Z (HOLO - APO)")
        ax.set_title(f"{nr_code} ({uid}): residue-level delta profile\n(green = within-receptor FDR significant)")
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"residue_profile_{nr_code}.png", dpi=150)
        plt.close(fig)


def build_pymol_table(residue_effects: pd.DataFrame, segment_effects: pd.DataFrame):
    sig_res = residue_effects[residue_effects["significant_within_receptor"]].copy()
    sig_res["color"] = np.where(sig_res["delta_median_z"] < 0, "blue", "red")
    sig_res["pymol_select_command"] = sig_res.apply(
        lambda r: f"select sig_{r['nr_code']}_{r['canonical_uniprot_position']}, "
                  f"(chain A and resi {r['canonical_uniprot_position']})  # UniProt numbering; "
                  f"map to auth_seq_id per receptor/PDB before use", axis=1)
    out_cols = ["uniprot_id", "nr_code", "canonical_uniprot_position", "structural_segment",
                "delta_median_z", "q_within_receptor", "color", "pymol_select_command"]
    out = TABLES / "apo_holo_pymol_ready_table.csv"
    sig_res[out_cols].to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(sig_res)} rows) -- "
          "canonical UniProt numbering; must be re-mapped to each structure's own "
          "auth_seq_id before loading into PyMOL (not done here -- no new residue "
          "mapping is invented in this figures script).")


def main():
    eligibility = pd.read_csv(TABLES / "apo_holo_receptor_eligibility.csv")
    segment_effects = pd.read_parquet(ANALYSIS / "apo_holo_segment_effects.parquet")
    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")
    residue_effects = pd.read_parquet(ANALYSIS / "apo_holo_residue_effects.parquet")

    assert (eligibility.loc[eligibility["uniprot_id"] == "P11473", "n_APO"] == 0).all(), \
        "VDR must have 0 APO entries -- sanity check before generating figures that exclude it"

    fig_counts_by_receptor(eligibility)
    print(f"Wrote {(FIG_DIR / 'apo_holo_counts_by_receptor.png').relative_to(PROJECT_ROOT)}")

    fig_delta_heatmap(segment_effects)
    print(f"Wrote {(FIG_DIR / 'receptor_segment_delta_heatmap.png').relative_to(PROJECT_ROOT)}")

    fig_forest(segment_effects)
    print(f"Wrote {(FIG_DIR / 'effect_forest_plot.png').relative_to(PROJECT_ROOT)}")

    fig_per_receptor_segments(features, eligibility)
    print(f"Wrote per-receptor figures to {FIG_DIR.relative_to(PROJECT_ROOT)}/per_receptor_*.png")

    fig_residue_profile(residue_effects)
    print(f"Wrote residue-profile figures to {FIG_DIR.relative_to(PROJECT_ROOT)}/residue_profile_*.png")

    build_pymol_table(residue_effects, segment_effects)


if __name__ == "__main__":
    main()
