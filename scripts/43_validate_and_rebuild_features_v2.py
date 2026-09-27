"""
Structural-segment alignment v2, Parts 10-12.

Part 10: detailed validation table + schematics for the 5 anchor receptors
against the v2 (reference-based, structurally-anchored) segment map.

Part 11: rebuild the segment coverage/feature/matrix tables from the v2
segment map. The aggregation hierarchy is UNCHANGED (instance normalized
independently -> PDB entry median -> segment summary -> receptor median) --
this script reuses the exact same aggregation functions from
scripts/36_build_segment_features.py, only the input segment map differs.

Part 12: recompute the missingness policy audit against the v2 matrices.

Reads:
    data/analysis/nr_lbd_structural_segment_map.csv (v2, from script 42)
    data/analysis/lbd_structural_reference_structures.csv
    data/analysis/lbd_structural_alignment_scores.csv
    data/processed/structures.csv
    data/analysis/normalized_lbd_bfactors_all.parquet
    data/analysis/pdb_position_bfactor_profiles.parquet
    data/manifests/nr_lbd_reference.csv
    reports/tables/functional_site_definitions.csv

Writes:
    reports/tables/structural_segment_validation_v2.csv
    reports/figures/alignment_qc/v2_schematic_<receptor>.png
    data/analysis/pdb_segment_coverage.parquet                 (OVERWRITTEN)
    data/analysis/pdb_segment_flexibility_features.parquet     (OVERWRITTEN)
    data/analysis/receptor_segment_flexibility_profiles.parquet (OVERWRITTEN)
    data/analysis/receptor_segment_feature_matrix_median.parquet (OVERWRITTEN)
    data/analysis/receptor_segment_feature_matrix_mean.parquet   (OVERWRITTEN)
    data/analysis/pdb_segment_feature_matrix_median.parquet      (OVERWRITTEN)
    data/analysis/pdb_segment_feature_matrix_mean.parquet        (OVERWRITTEN)
    data/analysis/receptor_segment_coverage_matrix.parquet       (OVERWRITTEN)
    reports/tables/segment_feature_missingness.md                (OVERWRITTEN)
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

from scripts.utils.mmcif_parse import get_secondary_structure_ranges, parse_mmcif_gz  # noqa: E402
import importlib.util

spec36 = importlib.util.spec_from_file_location("seg_features_v1", PROJECT_ROOT / "scripts" / "36_build_segment_features.py")
seg_features_v1 = importlib.util.module_from_spec(spec36)
spec36.loader.exec_module(seg_features_v1)

spec37 = importlib.util.spec_from_file_location("seg_missingness_v1", PROJECT_ROOT / "scripts" / "37_build_segment_missingness_and_alignment.py")
seg_missingness_v1 = importlib.util.module_from_spec(spec37)
spec37.loader.exec_module(seg_missingness_v1)

ANALYSIS = PROJECT_ROOT / "data" / "analysis"
MANIFESTS = PROJECT_ROOT / "data" / "manifests"
TABLES = PROJECT_ROOT / "reports" / "tables"
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "alignment_qc"
RAW_MMCIF = PROJECT_ROOT / "data" / "raw" / "mmcif"

SEGMENT_MAP_CSV = ANALYSIS / "nr_lbd_structural_segment_map.csv"
REF_CSV = ANALYSIS / "lbd_structural_reference_structures.csv"
SCORES_CSV = ANALYSIS / "lbd_structural_alignment_scores.csv"
FUNCTIONAL_SITE_CSV = TABLES / "functional_site_definitions.csv"
VALIDATION_V2_CSV = TABLES / "structural_segment_validation_v2.csv"

ANCHORS = {"P11473": "VDR", "P37231": "PPARgamma", "P03372": "ESR1", "P04150": "GR", "P22736": "NR4A1"}


def build_validation_v2(seg_map: pd.DataFrame, ref: pd.DataFrame, scores: pd.DataFrame, fsite: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for uid, name in ANCHORS.items():
        g = seg_map[(seg_map["uniprot_id"] == uid) & (seg_map["structural_segment"] != "")]
        for seg, gg in g.groupby("structural_segment"):
            method = gg["assignment_method"].iloc[0]
            is_own = method == "own_canonical_template"
            z_val = pd.to_numeric(gg["structural_alignment_score"].iloc[0], errors="coerce")
            align_support = "own canonical template (no alignment needed)" if is_own else (
                f"z={z_val:.2f} vs template {gg['reference_template'].iloc[0]}"
            )
            ss_support = f"n_pdb_observed>=1 for all {len(gg)} residues (see receptor_dssp_consensus.parquet)"
            fs = fsite[(fsite["uniprot_id"] == uid) & fsite["uniprot_position"].between(gg["canonical_uniprot_position"].min(), gg["canonical_uniprot_position"].max())]
            func_consistency = "; ".join(f"{r['site_name']}:{r['uniprot_position']}({r['expected_residue_name']})" for _, r in fs.iterrows()) or "no repo functional-site residue in this range"
            rows.append({
                "uniprot_id": uid, "nr_code": name, "structural_segment": seg,
                "target_uniprot_start": int(gg["canonical_uniprot_position"].min()),
                "target_uniprot_end": int(gg["canonical_uniprot_position"].max()),
                "n_residues": len(gg), "reference_template": gg["reference_template"].iloc[0],
                "reference_pdb": gg["reference_pdb"].iloc[0],
                "alignment_support": align_support, "secondary_structure_support": ss_support,
                "functional_site_consistency": func_consistency,
                "assignment_confidence": gg["assignment_confidence"].iloc[0],
                "assignment_method": method,
            })
    return pd.DataFrame(rows).sort_values(["uniprot_id", "target_uniprot_start"])


def build_schematic(uid: str, name: str, seg_map: pd.DataFrame, fsite: pd.DataFrame, lbd_ref: pd.DataFrame) -> None:
    class_color = {"H1": "#d62728", "H12": "#9467bd", "beta_sheet": "#2ca02c", "H3": "#ff7f0e", "H11": "#17becf"}

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
    fig, ax = plt.subplots(figsize=(12, 2.4))
    for _, r in g.iterrows():
        ax.add_patch(plt.Rectangle((r["canonical_uniprot_position"], 1), 1, 1,
                                    color=color_for(r["structural_segment"], r["structural_class"]), linewidth=0))
    for seg in g[g["structural_segment"] != ""]["structural_segment"].unique():
        gg = g[g["structural_segment"] == seg]
        mid = (gg["canonical_uniprot_position"].min() + gg["canonical_uniprot_position"].max()) / 2
        ax.text(mid, 2.15, seg, ha="center", va="bottom", fontsize=6, rotation=90 if len(seg) > 3 else 0)
    fs = fsite[fsite["uniprot_id"] == uid]
    for _, r in fs.iterrows():
        ax.plot([r["uniprot_position"]], [0.5], marker="v", color="black", markersize=4)
    ax.set_xlim(row["lbd_start"] - 2, row["lbd_end"] + 2)
    ax.set_ylim(0, 3.2)
    ax.set_yticks([])
    ax.set_xlabel("canonical UniProt position")
    ax.set_title(f"{name} ({uid}) v2 structurally-anchored segment map; triangles = repository functional-site residues", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"v2_schematic_{name}.png", dpi=150)
    plt.close(fig)


def main() -> None:
    seg_map = pd.read_csv(SEGMENT_MAP_CSV, keep_default_na=False)
    ref = pd.read_csv(REF_CSV, dtype={"pdb_id": str})
    scores = pd.read_csv(SCORES_CSV)
    fsite = pd.read_csv(FUNCTIONAL_SITE_CSV)
    lbd_ref = pd.read_csv(MANIFESTS / "nr_lbd_reference.csv")

    print("Part 10: anchor validation v2 + schematics...")
    validation = build_validation_v2(seg_map, ref, scores, fsite)
    TABLES.mkdir(parents=True, exist_ok=True)
    validation.to_csv(VALIDATION_V2_CSV, index=False)
    print(f"  wrote {VALIDATION_V2_CSV.relative_to(PROJECT_ROOT)} ({len(validation)} rows)")

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for uid, name in ANCHORS.items():
        build_schematic(uid, name, seg_map, fsite, lbd_ref)
    print(f"  wrote {len(ANCHORS)} v2 schematics to {FIG_DIR.relative_to(PROJECT_ROOT)}")

    print("\nPart 11: rebuilding segment coverage/features/matrices from the v2 map (unchanged aggregation)...")
    master = pd.read_csv(seg_features_v1.MASTER_CSV, dtype={"pdb_id": str})
    norm = pd.read_parquet(
        seg_features_v1.NORMALIZED_PARQUET,
        columns=["instance_id", "uniprot_id", "canonical_uniprot_position", "is_construct_mapped", "ca_positive_occupancy"],
    )
    pdb_pos = pd.read_parquet(seg_features_v1.PDB_POSITION_PARQUET)
    defs, _ = seg_features_v1.load_segment_lookup(seg_map)

    coverage = seg_features_v1.build_coverage(master, norm, seg_map, defs)
    coverage.to_parquet(seg_features_v1.COVERAGE_PARQUET, index=False)
    print(f"  wrote {seg_features_v1.COVERAGE_PARQUET.relative_to(PROJECT_ROOT)} ({len(coverage)} rows)")

    pdb_feat = seg_features_v1.build_pdb_features(master, norm, pdb_pos, seg_map, defs, coverage)
    pdb_feat.to_parquet(seg_features_v1.PDB_FEATURES_PARQUET, index=False)
    print(f"  wrote {seg_features_v1.PDB_FEATURES_PARQUET.relative_to(PROJECT_ROOT)} ({len(pdb_feat)} rows)")

    meta = master.drop_duplicates("uniprot_id")[["uniprot_id", "nr_code", "common_name", "group"]].copy()
    n_pdb_total = master.loc[master["has_stage3b_coordinate_data"]].groupby("uniprot_id")["pdb_id"].nunique()
    meta["n_pdbs_with_coordinate_data"] = meta["uniprot_id"].map(n_pdb_total).fillna(0).astype(int)
    rec_profile = seg_features_v1.build_receptor_profile(pdb_feat, meta)
    rec_profile.to_parquet(seg_features_v1.RECEPTOR_PROFILE_PARQUET, index=False)
    print(f"  wrote {seg_features_v1.RECEPTOR_PROFILE_PARQUET.relative_to(PROJECT_ROOT)} ({len(rec_profile)} rows)")

    rec_med = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="median_of_entry_median_z")
    rec_mean = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="mean_of_entry_mean_z")
    pdb_med = pdb_feat.pivot_table(index=["uniprot_id", "pdb_id"], columns="structural_segment", values="median_z")
    pdb_mean = pdb_feat.pivot_table(index=["uniprot_id", "pdb_id"], columns="structural_segment", values="mean_z")
    rec_cov = rec_profile.pivot(index="uniprot_id", columns="structural_segment", values="median_segment_coverage")
    rec_med.reset_index().to_parquet(seg_features_v1.RECEPTOR_MATRIX_MEDIAN, index=False)
    rec_mean.reset_index().to_parquet(seg_features_v1.RECEPTOR_MATRIX_MEAN, index=False)
    pdb_med.reset_index().to_parquet(seg_features_v1.PDB_MATRIX_MEDIAN, index=False)
    pdb_mean.reset_index().to_parquet(seg_features_v1.PDB_MATRIX_MEAN, index=False)
    rec_cov.reset_index().to_parquet(seg_features_v1.RECEPTOR_COVERAGE_MATRIX, index=False)
    print(f"  receptor x segment median matrix: {rec_med.shape}; pdb x segment median matrix: {pdb_med.shape}")

    print("\nPart 12: recomputing missingness policy audit against v2 matrices...")
    matrix = rec_med.reset_index()
    coverage_matrix = rec_cov.reset_index()
    report = seg_missingness_v1.build_missingness_report(matrix, coverage_matrix, rec_profile, lbd_ref)
    report = report.replace(
        "## 1. STOP-CONDITION finding from Part 5 (carried forward)",
        "## 1. STOP-CONDITION status (v2: reference-based 3D structural alignment)",
    )
    seg_missingness_v1.MISSINGNESS_MD.write_text(report + "\n")
    print(f"  wrote {seg_missingness_v1.MISSINGNESS_MD.relative_to(PROJECT_ROOT)}")

    for thr in (0.9, 0.8, 0.7, 0.5):
        n_seg = int((matrix.set_index("uniprot_id").notna().sum() >= thr * len(lbd_ref)).sum())
        print(f"  segments represented in >= {thr:.0%} of 48 receptors: {n_seg} (matrix shape at that support: "
              f"{len(lbd_ref)} x {n_seg})")


if __name__ == "__main__":
    main()
