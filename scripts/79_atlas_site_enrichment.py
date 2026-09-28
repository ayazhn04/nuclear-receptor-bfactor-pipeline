"""
Part 14 of the current functional-site flexibility stage: descriptively,
and then formally, ask whether the most flexible atlas segments
(data/analysis/superfamily_segment_flexibility_atlas.csv, from the
subfamily/atlas stage, commit f0c085c) are enriched for functional-site
overlap, and whether rigid core helices are depleted.

The universe and the threshold are defined BEFORE testing (Part 14's own
instruction: "do not select thresholds after seeing the result"): the
universe is every well-supported atlas segment (n_receptors >= 10, the
same explicit threshold already established and used in the atlas stage,
never re-chosen here), and "flexible"/"rigid" halves are defined by a
median split of that universe's median_receptor_z (not an arbitrary
post-hoc cutoff) -- a segment counts as "functional-site-associated" if
ANY of its residues, in ANY receptor, was ever assigned to a current
functional site (segment- or residue-level).
"""
from pathlib import Path

import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def main():
    atlas = pd.read_csv(ANALYSIS / "superfamily_segment_flexibility_atlas.csv")
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")

    # universe: well-supported segments only, exactly the threshold already established
    # and used throughout the atlas stage -- never re-picked here
    universe = atlas[atlas["well_supported"]].copy()
    print(f"Universe: {len(universe)} well-supported atlas segments (n_receptors >= 10, "
          "threshold fixed in the prior atlas stage, not re-chosen here)")

    site_associated_segments = set(definitions.loc[definitions["structural_segment"] != "", "structural_segment"].unique())
    # AF2_H12 and LBP rows are RESIDUE-level with structural_segment == "" in the
    # definitions table for LBP (pocket residues aren't tagged to one named segment) --
    # cross-reference those positions against the segment map directly instead.
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    pos_to_seg = seg_map.set_index(["uniprot_id", "canonical_uniprot_position"])["structural_segment"]
    lbp_af2 = definitions[definitions["functional_site"].isin(["AF2_H12", "LIGAND_BINDING_POCKET"])]
    for _, r in lbp_af2.iterrows():
        seg = pos_to_seg.get((r["uniprot_id"], r["canonical_uniprot_position"]), "")
        if seg:
            site_associated_segments.add(seg)

    universe["functional_site_associated"] = universe["structural_segment"].isin(site_associated_segments)
    print(f"{universe['functional_site_associated'].sum()}/{len(universe)} universe segments are "
          f"functional-site-associated: {sorted(universe.loc[universe['functional_site_associated'], 'structural_segment'])}")

    median_z = universe["median_receptor_z"].median()
    universe["flexible_half"] = universe["median_receptor_z"] >= median_z
    print(f"\nMedian split at median_receptor_z = {median_z:.3f} (fixed by the data's own median, "
          "not chosen after inspecting the enrichment result)")

    table = pd.crosstab(universe["flexible_half"], universe["functional_site_associated"])
    print("\n2x2 contingency table (rows=flexible half, cols=site-associated):")
    print(table.to_string())

    if table.shape == (2, 2):
        odds_ratio, p = stats.fisher_exact(table.values)
        print(f"\nFisher's exact test: odds ratio = {odds_ratio:.3f}, p = {p:.4f}")
    else:
        odds_ratio, p = None, None
        print("\nContingency table is degenerate (one cell empty) -- Fisher's exact test not computed.")

    out = TABLES / "atlas_functional_site_enrichment.csv"
    universe.to_csv(out, index=False)
    print(f"\nWrote {out.relative_to(PROJECT_ROOT)} ({len(universe)} rows)")

    summary = pd.DataFrame([{
        "universe_size": len(universe),
        "n_site_associated": int(universe["functional_site_associated"].sum()),
        "median_split_threshold": median_z,
        "odds_ratio": odds_ratio, "fisher_p": p,
    }])
    summary_out = TABLES / "atlas_functional_site_enrichment_test.csv"
    summary.to_csv(summary_out, index=False)
    print(f"Wrote {summary_out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
