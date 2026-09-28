"""
Part 10 of the APO vs HOLO flexibility stage: structural interpretation of
statistically supported results, using ONLY existing structural/functional
annotations already in this repository (canonical H1-H12/beta-sheet/loop
nomenclature established in the v3 structural-segment mapping, and the
Role-D `functional_site_definitions.csv`, read-only, never modified). No
new functional mapping is invented here.

Category assignment is a direct, literal reading of the canonical segment
name already used throughout this project (H12 = AF-2 activation helix,
H3/H4/H5 = coactivator/pocket-adjacent helices, H9/H10/H11 =
dimerization-interface helices, beta_sheet = the small beta-hairpin,
loop_* = major loops) -- not a new inference.

This script explicitly enumerates possible confounders (ligand identity,
agonist/antagonist state, mutations, construct length, crystallization
conditions, crystal packing, resolution/refinement, coactivator/
corepressor presence) and does not claim ligand binding CAUSES any
observed difference -- these are observational PDB comparisons only.
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

CONFOUNDERS = [
    "ligand identity (different HOLO structures may carry chemically distinct ligands)",
    "agonist vs antagonist/inverse-agonist pharmacology",
    "point mutations or engineered stabilizing mutations in the construct",
    "construct length / boundary differences between depositions",
    "crystallization conditions (buffer, pH, additives)",
    "crystal packing contacts (can artificially rigidify or perturb surface loops)",
    "resolution and refinement protocol differences",
    "presence or absence of a bound coactivator/corepressor peptide",
]

CATEGORY_MAP = {
    "H12": "H12 / AF-2 activation helix",
    "H3": "H3/H4/H5 coactivator/pocket region", "H4": "H3/H4/H5 coactivator/pocket region",
    "H5": "H3/H4/H5 coactivator/pocket region",
    "H9": "H9/H10/H11 dimerization region", "H10": "H9/H10/H11 dimerization region",
    "H11": "H9/H10/H11 dimerization region",
    "beta_sheet": "beta-sheet",
}


def categorize(segment: str) -> str:
    if segment in CATEGORY_MAP:
        return CATEGORY_MAP[segment]
    if segment.startswith("loop_"):
        return "major loop"
    if segment.startswith("H"):
        return "other helix (H1/H2/H6/H7/H8, not in the named functional categories above)"
    return "other"


def nearest_functional_site(uid: str, pos: int, sites: pd.DataFrame, window: int = 3):
    sub = sites[sites["uniprot_id"] == uid].copy()
    if len(sub) == 0:
        return None
    sub["distance"] = (sub["uniprot_position"] - pos).abs()
    nearest = sub.sort_values("distance").iloc[0]
    if nearest["distance"] <= window:
        return f"{nearest['site_name']} (residue {int(nearest['uniprot_position'])}, {int(nearest['distance'])} away)"
    return None


def main():
    sites = pd.read_csv(PROJECT_ROOT / "reports" / "tables" / "functional_site_definitions.csv")

    segment_effects = pd.read_parquet(ANALYSIS / "apo_holo_segment_effects.parquet")
    residue_effects = pd.read_parquet(ANALYSIS / "apo_holo_residue_effects.parquet")

    seg_sig = segment_effects[segment_effects["significant_within_receptor"]].copy()
    seg_sig["structural_category"] = seg_sig["structural_segment"].apply(categorize)
    seg_sig["nearest_functional_site_within_3res"] = ""  # segment-level has no single position

    res_sig = residue_effects[residue_effects["significant_within_receptor"]].copy()
    res_sig["structural_category"] = res_sig["structural_segment"].apply(
        lambda s: categorize(s) if s else "position not confidently mapped to a named segment"
    )
    res_sig["nearest_functional_site_within_3res"] = res_sig.apply(
        lambda r: nearest_functional_site(r["uniprot_id"], r["canonical_uniprot_position"], sites) or "", axis=1
    )

    print("=" * 70)
    print("SEGMENT-LEVEL significant results, structural category:")
    print("=" * 70)
    if len(seg_sig):
        print(seg_sig[["uniprot_id", "nr_code", "structural_segment", "structural_category",
                        "delta_median_z", "q_within_receptor"]].to_string(index=False))
    else:
        print("(none)")

    print("\n" + "=" * 70)
    print("RESIDUE-LEVEL significant results, structural category + nearest known functional site:")
    print("=" * 70)
    if len(res_sig):
        print(res_sig[["uniprot_id", "nr_code", "canonical_uniprot_position", "structural_segment",
                        "structural_category", "delta_median_z", "q_within_receptor",
                        "nearest_functional_site_within_3res"]].to_string(index=False))
    else:
        print("(none)")

    cat_counts_seg = seg_sig["structural_category"].value_counts()
    cat_counts_res = res_sig["structural_category"].value_counts()
    print("\nSegment-level category counts:")
    print(cat_counts_seg.to_string())
    print("\nResidue-level category counts:")
    print(cat_counts_res.to_string())

    out = TABLES / "apo_holo_structural_interpretation.csv"
    combined = pd.concat([
        seg_sig[["uniprot_id", "nr_code", "structural_segment", "structural_category",
                 "delta_median_z", "q_within_receptor", "nearest_functional_site_within_3res"]].assign(level="segment"),
        res_sig[["uniprot_id", "nr_code", "structural_segment", "structural_category",
                 "delta_median_z", "q_within_receptor", "nearest_functional_site_within_3res"]].assign(level="residue"),
    ], ignore_index=True)
    combined.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(combined)} rows)")

    print("\n" + "=" * 70)
    print("CONFOUNDERS (observational data -- causal claims are NOT made):")
    print("=" * 70)
    for c in CONFOUNDERS:
        print(f"  - {c}")
    print("\nNone of these confounders are adjusted for at this stage; every result above is a "
          "descriptive, observational PDB-level comparison, not a causal claim about ligand binding.")


if __name__ == "__main__":
    main()
