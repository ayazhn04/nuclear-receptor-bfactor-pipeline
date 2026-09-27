"""
Structural-segment alignment v3, Part 11: final validation table for the 5
anchor receptors against the v3 (fully curated, cross-template-rescued)
segment map, plus regenerated schematics.

Reads:
    data/analysis/nr_lbd_structural_segment_map.csv
    data/analysis/canonical_lbd_template_segments_v2.csv
    data/analysis/receptor_dssp_consensus.parquet
    reports/tables/functional_site_definitions.csv
    data/manifests/nr_lbd_reference.csv

Writes:
    reports/tables/structural_segment_validation_v3.csv
    reports/figures/alignment_qc/v3_schematic_<receptor>.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"
TABLES = PROJECT_ROOT / "reports" / "tables"
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "alignment_qc"

SEGMENT_MAP_CSV = ANALYSIS / "nr_lbd_structural_segment_map.csv"
TEMPLATE_V2_CSV = ANALYSIS / "canonical_lbd_template_segments_v2.csv"
CONSENSUS_PARQUET = ANALYSIS / "receptor_dssp_consensus.parquet"
FUNCTIONAL_SITE_CSV = TABLES / "functional_site_definitions.csv"

VALIDATION_V3_CSV = TABLES / "structural_segment_validation_v3.csv"

ANCHORS = {"P11473": "VDR", "P37231": "PPARgamma", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}
CANONICAL_ALL = ["H1", "H2", "H3", "H4", "H5", "beta_sheet", "H6", "H7", "H8", "H9", "H10", "H11", "H12"]


def build_validation(seg_map: pd.DataFrame, consensus: pd.DataFrame, fsite: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for uid, name in ANCHORS.items():
        g = seg_map[(seg_map["uniprot_id"] == uid) & (seg_map["structural_segment"] != "")]
        for seg, gg in g.groupby("structural_segment"):
            method = gg["assignment_method"].iloc[0]
            c = consensus[(consensus["uniprot_id"] == uid) & consensus["canonical_uniprot_position"].between(gg["canonical_uniprot_position"].min(), gg["canonical_uniprot_position"].max())]
            ss_agreement = f"{(c['consensus_dssp'] == 'HELIX').mean():.0%} consensus HELIX" if len(c) and gg['structural_class'].iloc[0] == "HELIX" else (
                f"{(c['consensus_dssp'] == 'SHEET').mean():.0%} consensus SHEET" if len(c) else "no consensus evidence")
            transfer_agreement = "own curated template (no transfer)" if method == "own_canonical_template" else gg["evidence_notes"].iloc[0][:160]
            fs = fsite[(fsite["uniprot_id"] == uid) & fsite["uniprot_position"].between(gg["canonical_uniprot_position"].min(), gg["canonical_uniprot_position"].max())]
            func_consistency = "; ".join(f"{r['site_name']}:{r['uniprot_position']}({r['expected_residue_name']})" for _, r in fs.iterrows()) or "no repo functional-site residue in this range"
            rows.append({
                "uniprot_id": uid, "nr_code": name, "structural_segment": seg,
                "canonical_uniprot_start": int(gg["canonical_uniprot_position"].min()),
                "canonical_uniprot_end": int(gg["canonical_uniprot_position"].max()),
                "n_residues": len(gg),
                "supporting_reference_structure": gg["reference_pdb"].iloc[0],
                "publication_evidence": gg["evidence_notes"].iloc[0][:200],
                "secondary_structure_agreement": ss_agreement,
                "structural_transfer_agreement": transfer_agreement,
                "confidence": gg["assignment_confidence"].iloc[0],
            })
    return pd.DataFrame(rows).sort_values(["uniprot_id", "canonical_uniprot_start"])


def build_schematic(uid: str, name: str, seg_map: pd.DataFrame, fsite: pd.DataFrame, lbd_ref: pd.DataFrame) -> None:
    class_color = {"H1": "#d62728", "H12": "#9467bd", "beta_sheet": "#2ca02c", "H3": "#ff7f0e", "H5": "#8c564b",
                   "H7": "#e377c2", "H11": "#17becf", "H6": "#bcbd22"}

    def color_for(seg: str, cls: str) -> str:
        if seg in class_color:
            return class_color[seg]
        if cls == "HELIX":
            return "#1f77b4"
        if cls == "LOOP":
            return "#c7c7c7"
        return "#eeeeee"

    g = seg_map[seg_map["uniprot_id"] == uid].sort_values("canonical_uniprot_position")
    row = lbd_ref[lbd_ref["uniprot_id"] == uid].iloc[0]
    fig, ax = plt.subplots(figsize=(13, 2.6))
    for _, r in g.iterrows():
        ax.add_patch(plt.Rectangle((r["canonical_uniprot_position"], 1), 1, 1,
                                    color=color_for(r["structural_segment"], r["structural_class"]), linewidth=0))
    for seg in g[g["structural_segment"] != ""]["structural_segment"].unique():
        gg = g[g["structural_segment"] == seg]
        mid = (gg["canonical_uniprot_position"].min() + gg["canonical_uniprot_position"].max()) / 2
        ax.text(mid, 2.15, seg, ha="center", va="bottom", fontsize=6, rotation=90)
    fs = fsite[fsite["uniprot_id"] == uid]
    for _, r in fs.iterrows():
        ax.plot([r["uniprot_position"]], [0.5], marker="v", color="black", markersize=4)
    ax.set_xlim(row["lbd_start"] - 2, row["lbd_end"] + 2)
    ax.set_ylim(0, 3.2)
    ax.set_yticks([])
    ax.set_xlabel("canonical UniProt position")
    n_named = g["structural_segment"].ne("").sum()
    ax.set_title(f"{name} ({uid}) v3 curated segment map -- {g[g.structural_segment!=''].structural_segment.nunique()} segments, {n_named} positions named; triangles = repo functional-site residues", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"v3_schematic_{name}.png", dpi=150)
    plt.close(fig)


def main() -> None:
    seg_map = pd.read_csv(SEGMENT_MAP_CSV, keep_default_na=False)
    consensus = pd.read_parquet(CONSENSUS_PARQUET)
    fsite = pd.read_csv(FUNCTIONAL_SITE_CSV)
    lbd_ref = pd.read_csv(MANIFESTS / "nr_lbd_reference.csv")

    validation = build_validation(seg_map, consensus, fsite)
    TABLES.mkdir(parents=True, exist_ok=True)
    validation.to_csv(VALIDATION_V3_CSV, index=False)
    print(f"Wrote {VALIDATION_V3_CSV.relative_to(PROJECT_ROOT)} ({len(validation)} rows)")

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for uid, name in ANCHORS.items():
        build_schematic(uid, name, seg_map, fsite, lbd_ref)
        n_seg = validation[validation["uniprot_id"] == uid]["structural_segment"].nunique()
        print(f"  {name}: {n_seg} curated segments -> v3_schematic_{name}.png")


if __name__ == "__main__":
    main()
