"""
Parts 7-8 of the current functional-site flexibility stage: PDB-level site
flexibility. For every receptor x PDB entry x available functional site,
summarize the already-approved per-instance-Z -> PDB-entry-median position
profile (data/analysis/pdb_position_bfactor_profiles.parquet) over the
site's defined positions, and compare against the non-site "background"
(every other standardized-LBD position observed in that same entry).

A site score is never computed from zero observed site residues (Part 8:
"Do not use a site score if there are no observed site residues"). An
explicit site-coverage criterion (fraction of the site's defined positions
actually observed in that PDB entry) is required and reported alongside
every row -- this matters especially for SEGMENT-level sites, where the
full nominal span is rarely all observed in every structure.
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"

MIN_SITE_COVERAGE = 0.0  # every row reports coverage explicitly; rows with 0 observed
                          # site residues are dropped entirely (never scored), not just
                          # flagged -- see build_site_flexibility()


def main():
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")
    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")

    site_position_sets = {
        (uid, site): set(g["canonical_uniprot_position"])
        for (uid, site), g in definitions.groupby(["uniprot_id", "functional_site"])
    }
    n_expected = {k: len(v) for k, v in site_position_sets.items()}
    def_level = definitions.groupby(["uniprot_id", "functional_site"])["definition_level"].first().to_dict()
    def_conf = definitions.groupby(["uniprot_id", "functional_site"])["assignment_confidence"] \
        .agg(lambda s: s.mode().iloc[0]).to_dict()

    rows = []
    for (uid, pdb_id), g in positions.groupby(["uniprot_id", "pdb_id"]):
        nr_code = g["nr_code"].iloc[0]
        g = g.dropna(subset=["b_factor_z"])
        if len(g) == 0:
            continue
        all_positions_observed = set(g["canonical_uniprot_position"])
        for (uid2, site), site_positions in site_position_sets.items():
            if uid2 != uid:
                continue
            observed_site = all_positions_observed & site_positions
            if len(observed_site) == 0:
                continue  # never score a site with zero observed residues in this entry
            site_vals = g.loc[g["canonical_uniprot_position"].isin(observed_site), "b_factor_z"]
            background_vals = g.loc[~g["canonical_uniprot_position"].isin(site_positions), "b_factor_z"]
            if len(background_vals) == 0:
                continue  # need a background to compare against too

            rows.append({
                "uniprot_id": uid, "nr_code": nr_code, "pdb_id": pdb_id, "functional_site": site,
                "definition_level": def_level.get((uid, site), ""),
                "site_definition_confidence": def_conf.get((uid, site), ""),
                "n_site_residues_expected": n_expected[(uid, site)],
                "n_site_residues_observed": len(observed_site),
                "site_coverage": len(observed_site) / n_expected[(uid, site)],
                "median_site_z": float(site_vals.median()),
                "mean_site_z": float(site_vals.mean()),
                "comparison_non_site_median_z": float(background_vals.median()),
                "n_background_residues_observed": len(background_vals),
                "delta_site_vs_background": float(site_vals.median() - background_vals.median()),
            })

    pdb_site = pd.DataFrame(rows).sort_values(["uniprot_id", "functional_site", "pdb_id"]).reset_index(drop=True)
    out = ANALYSIS / "pdb_functional_site_flexibility.parquet"
    pdb_site.to_parquet(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(pdb_site)} rows, "
          f"{pdb_site['uniprot_id'].nunique()} receptors, {pdb_site['functional_site'].nunique()} site classes)")

    print("\nRows per site class:")
    print(pdb_site["functional_site"].value_counts().to_string())
    print("\nMedian site_coverage by site class:")
    print(pdb_site.groupby("functional_site")["site_coverage"].median().to_string())
    print(f"\nAssertion: no row has n_site_residues_observed == 0 -> "
          f"{(pdb_site['n_site_residues_observed'] == 0).sum() == 0}")


if __name__ == "__main__":
    main()
