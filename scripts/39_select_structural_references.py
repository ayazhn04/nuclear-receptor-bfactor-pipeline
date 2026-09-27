"""
Structural-segment alignment v2, Part 2: select ONE primary reference LBD
structure per receptor with coordinate data, for establishing the
receptor-specific canonical segment map only. This is NOT an analysis
filter -- it does not touch structures.csv or any flexibility-feature
table; master-inventory membership is untouched.

Selection priority (in order, each a tie-break on the previous):
    1. human / human-derived identity (identity_class == HUMAN_OR_HUMAN_DERIVED)
    2. maximal standardized-LBD positive-occupancy coordinate coverage
    3. maximal B-factor-usable coverage (structural completeness proxy)
    4. best (lowest) resolution, then best (lowest) R-free, where available
    5. not a fusion/chimera construct (preferred, not required, if the only
       candidates are fusions)
    6. one clear, already-deterministically-selected LBD chain (reuses this
       receptor x PDB x entity candidate's own representative_instance_id)

Reads:
    data/processed/structures.csv

Writes:
    data/analysis/lbd_structural_reference_structures.csv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"

MASTER_CSV = PROCESSED / "structures.csv"
OUT_CSV = ANALYSIS / "lbd_structural_reference_structures.csv"


def select_reference(g: pd.DataFrame) -> pd.Series:
    human = g[g["identity_class"] == "HUMAN_OR_HUMAN_DERIVED"]
    pool = human if len(human) else g
    non_fusion = pool[~pool["fusion_or_chimera"]]
    if len(non_fusion):
        pool = non_fusion
    ranked = pool.sort_values(
        by=["positive_occupancy_lbd_coverage", "bfactor_usable_lbd_coverage", "resolution", "r_free"],
        ascending=[False, False, True, True],
        na_position="last",
    )
    return ranked.iloc[0]


def main() -> None:
    master = pd.read_csv(MASTER_CSV, dtype={"pdb_id": str})
    coord = master[master["has_stage3b_coordinate_data"]].copy()

    rows = []
    for uid, g in coord.groupby("uniprot_id"):
        r = select_reference(g)
        n_human = int((g["identity_class"] == "HUMAN_OR_HUMAN_DERIVED").sum())
        n_nonfusion = int((~g["fusion_or_chimera"]).sum())
        reasons = [
            f"identity: {'human/human-derived' if r['identity_class'] == 'HUMAN_OR_HUMAN_DERIVED' else r['identity_class'] + ' (no human candidate available)'} "
            f"({n_human}/{len(g)} candidates human)",
            f"positive_occupancy_lbd_coverage={r['positive_occupancy_lbd_coverage']:.4f} (best available among the filtered pool)",
            f"bfactor_usable_lbd_coverage={r['bfactor_usable_lbd_coverage']:.4f}",
            f"resolution={r['resolution']}, r_free={r['r_free']}",
            f"fusion_or_chimera={r['fusion_or_chimera']} ({n_nonfusion}/{len(g)} candidates non-fusion)",
            f"chain: representative_instance_id={r['representative_instance_id']} (label_asym_id={r['label_asym_id']}, auth_asym_id={r['auth_asym_id']})",
        ]
        rows.append(
            {
                "uniprot_id": uid, "nr_code": r["nr_code"], "common_name": r["common_name"],
                "pdb_id": r["pdb_id"], "polymer_entity_id": r["polymer_entity_id"],
                "instance_id": r["representative_instance_id"],
                "label_asym_id": r["label_asym_id"], "auth_asym_id": r["auth_asym_id"],
                "standardized_lbd_start": r["standardized_lbd_start"], "standardized_lbd_end": r["standardized_lbd_end"],
                "positive_occupancy_lbd_coverage": r["positive_occupancy_lbd_coverage"],
                "bfactor_usable_lbd_coverage": r["bfactor_usable_lbd_coverage"],
                "resolution": r["resolution"], "r_free": r["r_free"],
                "identity_class": r["identity_class"], "fusion_or_chimera": bool(r["fusion_or_chimera"]),
                "n_candidates_considered": len(g), "n_human_candidates": n_human, "n_non_fusion_candidates": n_nonfusion,
                "selection_reason": "; ".join(reasons),
            }
        )

    out = pd.DataFrame(rows).sort_values("uniprot_id").reset_index(drop=True)
    assert out["uniprot_id"].is_unique
    assert len(out) == coord["uniprot_id"].nunique()

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV.relative_to(PROJECT_ROOT)} ({len(out)} rows)")
    print(f"  non-human reference chosen for: {sorted(out.loc[out['identity_class'] != 'HUMAN_OR_HUMAN_DERIVED', 'nr_code'].tolist())}")
    print(f"  fusion/chimera reference chosen for: {sorted(out.loc[out['fusion_or_chimera'], 'nr_code'].tolist())}")


if __name__ == "__main__":
    main()
