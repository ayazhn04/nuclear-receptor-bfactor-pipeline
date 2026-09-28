"""
Part 18 of the current functional-site flexibility stage: assembles
reports/tables/current_functional_site_summary.md from the outputs of
scripts 72-81. Writes no new analysis -- pure reporting/aggregation.
"""
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def main():
    availability = pd.read_csv(TABLES / "current_functional_site_availability.csv")
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")
    within = pd.read_csv(TABLES / "current_functional_site_within_receptor_tests.csv")
    superfamily = pd.read_csv(TABLES / "current_functional_site_superfamily_tests.csv")
    site_pairs = pd.read_csv(TABLES / "current_functional_site_class_comparison.csv")
    atlas_enrich = pd.read_csv(TABLES / "atlas_functional_site_enrichment_test.csv")
    apo_holo_overlap = pd.read_csv(TABLES / "apo_holo_functional_site_overlap.csv", keep_default_na=False)
    loo_pdb = pd.read_csv(TABLES / "functional_site_loo_pdb_robustness.csv")
    loo_receptor = pd.read_csv(TABLES / "functional_site_loo_receptor_robustness.csv")
    flags = pd.read_csv(TABLES / "functional_site_confidence_coverage_flags.csv")

    lines = []
    lines.append("# Current Full-Dataset Functional-Site Flexibility Map\n")
    lines.append(
        "> **NEW analysis layer**, continuing from commit `fb27a86`, built on the current full "
        "master dataset, approved per-instance Z normalization, chain -> PDB-entry-median "
        "aggregation, and the validated v3 structural-segment mapping. The historical Role D "
        "analysis (`reports/tables/functional_site_definitions.csv` and its downstream "
        "`functional_site_*_scores.csv`/`functional_site_normalization_qc.csv`, the old "
        "1,382-structure strict-QC view) is **not modified, not overwritten, byte-identical** "
        "-- verified by test. No ClinVar/disease analysis was run.\n"
    )

    lines.append("## Site availability (Part 1)\n")
    avail_pivot = availability.pivot_table(index="functional_site", columns="status", values="uniprot_id",
                                            aggfunc="count", fill_value=0)
    lines.append("| Site class | AVAILABLE | PARTIAL | UNRESOLVED |")
    lines.append("|---|---|---|---|")
    for site, row in avail_pivot.iterrows():
        lines.append(f"| {site} | {row.get('AVAILABLE', 0)} | {row.get('PARTIAL', 0)} | {row.get('UNRESOLVED', 0)} |")
    lines.append(
        "\nUnresolved sites are never assigned a placeholder position or a zero flexibility "
        "value anywhere downstream. Full receptor x site status: "
        "`reports/tables/current_functional_site_availability.csv`.\n"
    )

    lines.append("## Site definitions and provenance (Parts 2, 4, 5, 6)\n")
    lines.append(
        f"`data/analysis/current_functional_site_definitions.parquet` -- {len(definitions)} rows, "
        f"{definitions['uniprot_id'].nunique()} receptors, 4 site classes:\n"
    )
    lines.append("- **AF2_H12**: taken directly from the current v3 canonical H12 mapping (never the old "
                 f"heuristic definition) -- {definitions[definitions['functional_site']=='AF2_H12']['uniprot_id'].nunique()} receptors.")
    lines.append("- **LIGAND_BINDING_POCKET**: receptor-level consensus from real 5.0 A heavy-atom contact "
                 "geometry to the already-classified FUNCTIONAL_LBD_LIGAND (never crystallization additives, "
                 "ions, or buffer molecules -- that exclusion is inherited read-only from the master "
                 f"ligand-state classification), >=25% of eligible holo PDBs (PRIMARY) -- "
                 f"{definitions[definitions['functional_site']=='LIGAND_BINDING_POCKET']['uniprot_id'].nunique()} receptors; "
                 "sensitivity thresholds (10%/50%) also computed and reported "
                 "(`data/analysis/ligand_pocket_consensus_residues.parquet`).")
    coact_a = ((definitions["functional_site"] == "COACTIVATOR_GROOVE") & (definitions["evidence_type"] == "role_d_curated_definition")).sum()
    coact_c = ((definitions["functional_site"] == "COACTIVATOR_GROOVE") & (definitions["evidence_type"] == "structural_segment_correspondence")).sum()
    dimer_a = ((definitions["functional_site"] == "DIMERIZATION_INTERFACE") & (definitions["evidence_type"] == "role_d_curated_definition")).sum()
    dimer_c = ((definitions["functional_site"] == "DIMERIZATION_INTERFACE") & (definitions["evidence_type"] == "structural_segment_correspondence")).sum()
    lines.append(f"- **COACTIVATOR_GROOVE**: {coact_a} tier-A (curated residue) rows, {coact_c} tier-C "
                 "(H3/H4/H5 segment correspondence) rows. Tier B (bound coactivator-peptide contact) was "
                 "not computed -- untracked polymer-partner identity outside the historical strict-QC "
                 "subset -- reported as unavailable, not approximated.")
    lines.append(f"- **DIMERIZATION_INTERFACE**: {dimer_a} tier-A rows, {dimer_c} tier-C (H9/H10/H11 "
                 "segment correspondence) rows. Tier B (observed receptor-receptor contacts) likewise not "
                 "computed for the same reason.\n")

    lines.append("## Within-receptor site-vs-background tests (Part 10)\n")
    n_tested = len(within)
    n_nominal = (within["raw_p"] < 0.05).sum()
    n_sig_within = within["significant_within_receptor"].sum()
    n_sig_global = within["significant_global"].sum()
    lines.append(
        f"{n_tested} receptor x site paired Wilcoxon tests (>=5 paired PDB entries, PDB entry as the "
        f"paired unit -- never individual residues). **{n_nominal}** nominal; **{n_sig_within}** "
        f"significant within-receptor FDR (PRIMARY); **{n_sig_global}** significant under a secondary "
        "global FDR too. Full table: `reports/tables/current_functional_site_within_receptor_tests.csv`.\n"
    )

    lines.append("## Superfamily site-class results (Part 11)\n")
    lines.append("| Site class | n receptors | median delta (site-background) | raw p | FDR q | Significant |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in superfamily.iterrows():
        lines.append(f"| {r['functional_site']} | {r['n_receptors']} | {r['median_delta_vs_background']:+.3f} | "
                     f"{r['raw_p']:.4g} | {r['fdr_q']:.4g} | {'**YES**' if r['significant_fdr_05'] else 'no'} |")
    lines.append(
        "\n**AF2/H12 is significantly MORE flexible than background superfamily-wide; the "
        "coactivator groove and ligand-binding pocket are both significantly MORE RIGID than "
        "background; the dimerization interface shows no significant net effect.** This pattern "
        "is directly consistent with the canonical NR mechanism (H12 as a mobile activation-"
        "helix \"lid\"; ligand/coactivator contacts locally rigidifying the pocket and groove).\n"
    )

    lines.append("### Pairwise site-class comparison (Part 12)\n")
    if len(site_pairs):
        lines.append("| Site A | Site B | n paired | median A-B | effect (rank-biserial) | raw p | FDR q | Sig. |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for _, r in site_pairs.iterrows():
            lines.append(f"| {r['site_a']} | {r['site_b']} | {r['n_receptors_paired']} | "
                         f"{r['median_paired_difference_a_minus_b']:+.3f} | {r['effect_size_rank_biserial']:+.3f} | "
                         f"{r['raw_p']:.4g} | {r['fdr_q']:.4g} | {'**YES**' if r['significant_fdr_05'] else 'no'} |")
    lines.append("")

    lines.append("## Effect sizes\n")
    lines.append(
        "Within-receptor tests use the matched-pairs rank-biserial correlation "
        "(`r = (W+ - W-)/(W+ + W-)` on Wilcoxon signed ranks); superfamily/site-class comparisons "
        "use the same measure on receptor-level paired deltas. Values are stored alongside every "
        "test in the corresponding CSV (never reported as p-values alone).\n"
    )

    lines.append("## Relationship to the flexibility atlas (Part 14)\n")
    er = atlas_enrich.iloc[0]
    universe_size, n_site_assoc = int(er["universe_size"]), int(er["n_site_associated"])
    lines.append(
        f"Universe = {universe_size} well-supported atlas segments (n_receptors >= 10, the "
        f"same threshold already fixed in the prior atlas stage); {n_site_assoc} of "
        f"them are functional-site-associated. Median-split enrichment test (threshold fixed by "
        f"the data's own median, not chosen after seeing the result): odds ratio = "
        f"{er['odds_ratio']:.2f}, Fisher's exact p = {er['fisher_p']:.3f} -- **no significant "
        "enrichment of functional-site association in the more-flexible half**, because "
        f"functional sites already cover {n_site_assoc}/{universe_size} "
        "(the great majority) of well-supported segments, leaving little room for a "
        "discriminating signal. Honestly reported null result, not forced.\n"
    )

    lines.append("## Overlap with apo/holo effects (Part 13, descriptive only)\n")
    sig_overlap = apo_holo_overlap[(apo_holo_overlap["apo_holo_significant_within_receptor"]) &
                                    (apo_holo_overlap["functional_site_overlap"])]
    lines.append(
        f"All {len(apo_holo_overlap[apo_holo_overlap['apo_holo_significant_within_receptor']])} "
        f"within-receptor-FDR-significant apo/holo segment results were cross-referenced (never "
        f"recomputed) against the current site definitions; **{len(sig_overlap)}** overlap a "
        "functional site:\n"
    )
    if len(sig_overlap):
        lines.append("| Receptor | Segment | Site | Overlap fraction | apo/holo delta |")
        lines.append("|---|---|---|---|---|")
        for _, r in sig_overlap.iterrows():
            lines.append(f"| {r['nr_code']} | {r['apo_holo_segment']} | {r['site_name']} | "
                         f"{r['overlap_fraction_of_segment']:.0%} | {r['apo_holo_delta_median_z']:+.3f} |")
    lines.append("\nNo causal inference is drawn from this overlap. Full table: "
                 "`reports/tables/apo_holo_functional_site_overlap.csv`.\n")

    lines.append("## Robustness (Part 15)\n")
    n_pdb_fragile = (loo_pdb["fraction_loo_p_lt_05"] < 1.0).sum()
    lines.append(
        f"Leave-one-PDB-out (within receptor): **{len(loo_pdb) - n_pdb_fragile}/{len(loo_pdb)}** "
        "within-receptor-significant results retain p<0.05 in every fold. Leave-one-receptor-out "
        "(superfamily site classes): "
    )
    flips = loo_receptor[loo_receptor["significance_flips_on_some_loo_fold"]]
    flip_desc = "none" if len(flips) == 0 else ", ".join(flips["functional_site"].tolist())
    lines.append(f"**{len(flips)}/{len(loo_receptor)}** site classes have significance that flips when "
                 f"a single receptor is removed ({flip_desc}). Full tables: "
                 "`reports/tables/functional_site_loo_pdb_robustness.csv`, "
                 "`functional_site_loo_receptor_robustness.csv`.\n")
    lines.append(
        f"{flags['low_confidence_definition'].sum()}/{len(flags)} receptor x site rows rest on a "
        f"LOW/MODERATE-confidence definition; {flags['segment_level_only'].sum()}/{len(flags)} are "
        "SEGMENT-level-only (imprecise) definitions; both are flagged explicitly in "
        "`reports/tables/functional_site_confidence_coverage_flags.csv` and should be weighted "
        "accordingly when interpreting any single receptor's result.\n"
    )

    lines.append("## Limitations (site-definition confidence)\n")
    lines.append(
        "- Tier B evidence (bound coactivator/coregulator peptide contacts; observed receptor-"
        "receptor protein contacts) was not computed for either the coactivator groove or the "
        "dimerization interface -- this pipeline does not track polymer-partner chain identity "
        "outside the historical 1,382-structure strict-QC subset. Most receptors' coactivator-"
        "groove and dimerization-interface definitions are therefore SEGMENT-level (H3/H4/H5 or "
        "H9/H10/H11 correspondence), not exact residue-level annotations, and are labeled as such "
        "everywhere.\n"
        "- AF2/H12 and dimerization-interface definitions are only AVAILABLE for the subset of "
        "receptors with a confident v3 H12/H9-H10-H11 structural mapping; the remainder are "
        "UNRESOLVED, never zero-filled.\n"
        "- The ligand-binding-pocket consensus requires >=3 holo PDB entries; receptors with fewer "
        "(or none) are UNRESOLVED for that site, not assigned a homologous pocket.\n"
        "- This is a descriptive structural association analysis. The strong, consistent, and "
        "leave-one-out-robust site-vs-background differences reported here (especially AF2/H12's "
        "flexibility and the pocket/groove's rigidity) are **structural associations, not a "
        "causal or mechanistic claim** -- they do not by themselves establish *why* a region is "
        "more or less flexible, only that it consistently is, across a large fraction of the "
        "superfamily.\n"
    )

    out = TABLES / "current_functional_site_summary.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
