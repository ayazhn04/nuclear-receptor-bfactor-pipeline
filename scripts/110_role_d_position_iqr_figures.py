"""Figure for the Role D v2 contextual-heterogeneity analysis."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "reports/role_d/contextual_iqr/role_d_functional_position_iqr.csv"
OUTPUT = ROOT / "reports/role_d/contextual_iqr/role_d_functional_position_iqr_scatter.png"


def main() -> None:
    data = pd.read_csv(INPUT)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    names = {"AF2_H12": "AF-2/H12", "LIGAND_BINDING_POCKET": "Ligand pocket"}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharex=True, sharey=True)
    for ax, (site, label) in zip(axes, names.items()):
        d = data[data.functional_site.eq(site)]
        sparse = d[~d.eligible_primary.astype(bool)]
        eligible = d[d.eligible_primary.astype(bool)]
        ax.scatter(sparse.median_z, sparse.iqr_z, s=10, color="#bdbdbd", alpha=.55, label="<5 PDB entries")
        ax.scatter(eligible.median_z, eligible.iqr_z, s=13, color="#2166ac", alpha=.7, label="≥5 PDB entries")
        ax.axvline(0, color="black", linewidth=.7)
        ax.set_title(label)
        ax.set_xlabel("Median PDB-entry Z")
    axes[0].set_ylabel("Between-PDB position IQR of Z")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Functional residues: average relative mobility versus contextual heterogeneity")
    fig.text(.5, .01, "Each point is one canonical residue; IQR is descriptive crystallographic-disorder heterogeneity, not solution dynamics.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .06, 1, .93))
    fig.savefig(OUTPUT, dpi=300)
    plt.close(fig)
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
