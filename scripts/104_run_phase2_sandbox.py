"""
Phase 2 of the v2 rebuild: re-run the legacy analyses (scripts 48-95) on the
corrected v2 foundation, inside the existing v2 sandbox (data/interim/rebuild_v2/,
built by scripts/101 + scripts/103).

The working-tree scripts are never edited. Every approved adaptation (G34,
A01-A28, S1-S3) is applied as a logged textual patch to the SANDBOX copy of the
legacy script; each patch must match its target exactly once. The report
generators 54/62/71/82/92 are not run: scripts/105 writes the same Markdown files
from the v2 result tables only (A24), and the unchanged script 95 then checks them.

Safety
    * before the run, every v1 copy in the sandbox's analysis/report folders is
      removed (only the v2 foundation and the static inputs listed in STATIC_INPUTS
      remain), so a script can never silently read a v1 result;
    * every file each script opens is logged (Python audit hook) and checked: a read
      is allowed only if the file is part of the v2 foundation, was written earlier
      in this run, or is a listed static input. Anything else stops the run;
    * network access is blocked.

Variants
    v2                patched run on the v2 foundation (the Phase 2 result)
    v1_repro --legacy unpatched legacy 48-95 on the v1 reproduction sandbox, compared
                      with the frozen tracked outputs (proves the sandbox reproduces the
                      released analyses before any adaptation is applied)

Usage:
    python scripts/104_run_phase2_sandbox.py --variant v2 [--export]
    python scripts/104_run_phase2_sandbox.py --variant v1_repro --legacy --compare
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
EXPORT_DATA = PROJECT_ROOT / "data" / "rebuild_v2"
EXPORT_REPORTS = PROJECT_ROOT / "reports" / "rebuild_v2"
NO_NETWORK_ENV = {"HTTP_PROXY": "http://127.0.0.1:9", "HTTPS_PROXY": "http://127.0.0.1:9",
                  "ALL_PROXY": "http://127.0.0.1:9", "http_proxy": "http://127.0.0.1:9",
                  "https_proxy": "http://127.0.0.1:9", "NO_PROXY": "", "no_proxy": "", "PYTHONHASHSEED": "0"}

# folders whose v1 copies are removed before the run (they hold analysis results)
CLEAN_DIRS = ["data/analysis", "reports/tables", "reports/figures"]
# v1 files that are legitimate static inputs of Phase 2 (not analysis results)
STATIC_INPUTS = {
    "reports/tables/functional_site_definitions.csv",   # Role D literature table (UniProt numbering; no TRalpha rows)
    "data/manifests/rcsb_discovery_anomalies.csv",       # discovery stage; discovery candidates unchanged (2,072)
    "data/manifests/nr_metadata_validated.csv",          # the 48-receptor project definition (unchanged by design)
}
STATIC_PREFIXES = ("config/", "scripts/", "data/raw/")
SCAN_PREFIXES = ("data/manifests/", "data/processed/", "README.md", "TEAM_HANDOFF.md", "reports/")
# legacy fidelity run only: frozen v1 inputs of Phase 2 produced by stages 34-47 (not rebuilt by 101)
LEGACY_V1_INPUTS = {f"data/analysis/{f}" for f in [
    "nr_lbd_structural_segment_map.csv", "receptor_dssp_consensus.parquet",
    "pdb_segment_coverage.parquet", "pdb_segment_flexibility_features.parquet",
    "receptor_segment_flexibility_profiles.parquet", "receptor_segment_feature_matrix_median.parquet",
    "receptor_segment_feature_matrix_mean.parquet", "pdb_segment_feature_matrix_median.parquet",
    "pdb_segment_feature_matrix_mean.parquet", "receptor_segment_coverage_matrix.parquet"]}
NEW_UTILS = ["canonical_helix.py", "apo_holo_scope.py", "matching_context.py"]

SUMMARY_SCRIPT = "105_phase2_v2_summaries.py"
REPORT_SCRIPTS_REPLACED = {"54", "62", "71", "82", "92"}  # A24: replaced by SUMMARY_SCRIPT
S3_SENSITIVITY_ENV = "S3_SENSITIVITY_EXCLUDE_REVIEW_REQUIRED"

STEPS = [
    "48_feature_support_audit.py", "49_build_modeling_cohort_and_matrix.py", "50_pca_analysis.py",
    "51_umap_analysis.py", "52_clustering_analysis.py", "53_entry_level_sensitivity.py", "54_modeling_summary.py",
    "55_subfamily_group_definition.py", "56_subfamily_kruskal_omnibus.py", "57_subfamily_posthoc.py",
    "58_flexibility_atlas.py", "59_loo_robustness.py", "60_mean_sensitivity_omnibus.py", "61_subfamily_figures.py",
    "62_subfamily_summary.py",
    "63_apo_holo_eligibility.py", "64_apo_holo_segment_tests.py", "65_apo_holo_residue_tests.py",
    "66_apo_holo_structural_interpretation.py", "67_apo_holo_qc_balance.py", "68_apo_holo_redundancy_audit.py",
    "69_apo_holo_loo_robustness.py", "70_apo_holo_figures.py", "71_apo_holo_summary.py",
    "72_ligand_pocket_geometry.py", "73_current_functional_site_definitions.py", "74_pdb_site_flexibility.py",
    "75_receptor_site_summary.py", "76_within_receptor_site_tests.py", "77_superfamily_site_tests.py",
    "78_apo_holo_site_overlap.py", "79_atlas_site_enrichment.py", "80_site_robustness.py",
    "81_functional_site_figures.py", "82_functional_site_summary.py", "83_functional_site_fdr_audit.py",
    "84_definition_quality_stratification.py", "85_matched_background_tests.py",
    "86_background_sensitivity_comparison.py", "87_superfamily_matched_reanalysis.py",
    "88_apo_holo_overlap_confidence_recheck.py", "89_matched_robustness.py", "90_interpretation_classification.py",
    "91_specificity_figures.py", "92_specificity_summary.py",
    "93_final_receptor_identity_audit.py", "94_final_project_numbers.py", "95_final_statistical_consistency_audit.py",
]
# v2 run: report generators replaced; S3 sensitivity re-runs 87/89/90 with the flag set
S3_SENSITIVITY_STEPS = ["87_superfamily_matched_reanalysis.py", "89_matched_robustness.py",
                        "90_interpretation_classification.py"]

ROOT_LINE = "PROJECT_ROOT = Path(__file__).resolve().parents[1]"


def _imports(*mods: str) -> str:
    return (ROOT_LINE + "\nimport sys as _sys  # sandbox patch (Phase 2): v2 helpers\n"
            "_sys.path.insert(0, str(PROJECT_ROOT))\n"
            + "".join(f"from scripts.utils import {m}  # noqa: E402\n" for m in mods))


def _s3_suffix_lines() -> str:
    return ('import os as _os  # sandbox patch (S3 sensitivity output naming)\n'
            f'_S3 = _os.environ.get("{S3_SENSITIVITY_ENV}") == "1"\n'
            '_SFX = "_sensitivity_excl_review_required" if _S3 else ""\n')


ELIGIBLE_S3 = ('(matched["site_definition_confidence"].isin({conf}) | '
               '(matched["definition_source"] == "PROFESSOR_GROUND_TRUTH"))')

LEGACY_CONFIDENT_LIST = '["HIGH", "MODERATE"]'

# (script, old, new, adaptation id, reason)
PATCHES: list[tuple[str, str, str, str, str]] = [
    # ---- A02 (48): professor vocabulary and professor-status reporting columns
    ("48_feature_support_audit.py", ROOT_LINE, _imports("canonical_helix"), "A02", "import professor vocabulary"),
    ("48_feature_support_audit.py",
     'CANONICAL_CORE_SEGMENTS = [f"H{i}" for i in range(1, 13)] + ["beta_sheet"]',
     "CANONICAL_CORE_SEGMENTS = list(canonical_helix.ELEMENT_ORDER)  # sandbox patch (A02): professor elements",
     "A02", "audited vocabulary = professor element order (no beta_sheet, which is not a professor element)"),
    ("48_feature_support_audit.py",
     '        conf_high_frac = float(conf_counts.get("HIGH", 0) / n_conf_rows) if n_conf_rows else np.nan\n'
     '        conf_moderate_frac = float(conf_counts.get("MODERATE", 0) / n_conf_rows) if n_conf_rows else np.nan\n'
     '        conf_low_frac = float(conf_counts.get("LOW", 0) / n_conf_rows) if n_conf_rows else np.nan\n',
     '        # sandbox patch (A02): professor assignment status (verbatim) instead of v3 confidence\n'
     '        st_pub = float(conf_counts.get("published_boundaries", 0) / n_conf_rows) if n_conf_rows else np.nan\n'
     '        st_hom = float(conf_counts.get("structural_homology", 0) / n_conf_rows) if n_conf_rows else np.nan\n'
     '        st_rev = float(conf_counts.get("review_required", 0) / n_conf_rows) if n_conf_rows else np.nan\n',
     "A02", "report professor status fractions"),
    ("48_feature_support_audit.py",
     '            "confidence_high_fraction": conf_high_frac,\n'
     '            "confidence_moderate_fraction": conf_moderate_frac,\n'
     '            "confidence_low_fraction": conf_low_frac,\n'
     '            "n_positions_own_canonical_template": n_own_template,\n'
     '            "n_positions_structural_alignment_transfer": n_alignment_transfer,\n'
     '            "n_positions_complementary_structure_rescue": n_complementary_rescue,\n',
     '            "professor_status_published_boundaries_fraction": st_pub,\n'
     '            "professor_status_structural_homology_fraction": st_hom,\n'
     '            "professor_status_review_required_fraction": st_rev,\n'
     '            "n_positions_professor_canonical_table": int(method_counts.get("PROFESSOR_CANONICAL_TABLE_FINAL", 0)),\n',
     "A02", "v3 method/confidence columns replaced by professor provenance columns"),
    # ---- A06 (58)
    ("58_flexibility_atlas.py",
     '    assert "H8" not in profiles["structural_segment"].unique(), (\n'
     '        "H8 was never confidently resolved -- it must never receive an atlas value"\n'
     '    )\n',
     "    # sandbox patch (A06): v3-only invariant removed; the professor map defines H8\n",
     "A06", "remove v3 invariant"),
    # ---- S1 / A07 (63-70): TRACE_LBD entries out of apo/holo scope, both sides
    ("63_apo_holo_eligibility.py", ROOT_LINE, _imports("apo_holo_scope"), "A07", "import S1 helper"),
    ("63_apo_holo_eligibility.py",
     '    ).reset_index()\n\n    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n',
     '    ).reset_index()\n    coord_pdbs = apo_holo_scope.apply(coord_pdbs)  # sandbox patch (S1)\n\n'
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n'
     '    features = apo_holo_scope.apply(features)  # sandbox patch (S1)\n',
     "A07", "S1 applied to master entry states and to the element features"),
    ("63_apo_holo_eligibility.py",
     '            "n_APO": n_apo, "n_HOLO": n_holo, "n_AMBIGUOUS": n_ambiguous, "n_NOT_ASSESSABLE": n_not_assessable,\n',
     '            "n_APO": n_apo, "n_HOLO": n_holo, "n_AMBIGUOUS": n_ambiguous, "n_NOT_ASSESSABLE": n_not_assessable,\n'
     '            "n_excluded_TRACE_LBD": int(counts.get(apo_holo_scope.EXCLUDED_STATE, 0)),  # sandbox patch (S1)\n',
     "A07", "report the excluded entries"),
    ("64_apo_holo_segment_tests.py", ROOT_LINE, _imports("apo_holo_scope"), "A07", "import S1 helper"),
    ("64_apo_holo_segment_tests.py",
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n',
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n'
     '    features = apo_holo_scope.apply(features)  # sandbox patch (S1)\n', "A07", "S1"),
    ("65_apo_holo_residue_tests.py", ROOT_LINE, _imports("apo_holo_scope", "canonical_helix"), "A07/A08", "imports"),
    ("65_apo_holo_residue_tests.py",
     '    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")\n',
     '    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")\n'
     '    positions = apo_holo_scope.apply(positions)  # sandbox patch (S1)\n', "A07", "S1"),
    ("65_apo_holo_residue_tests.py",
     '    seg_lookup = seg_map.set_index(["uniprot_id", "canonical_uniprot_position"])["structural_segment"]\n',
     '    # sandbox patch (A08): membership lookup; a dual professor position gets its joined DISPLAY label\n'
     '    _members = canonical_helix.position_memberships(\n'
     '        seg_map[["uniprot_id", "canonical_uniprot_position", "structural_segment"]].itertuples(index=False))\n'
     '    seg_lookup = {k: canonical_helix.display_label(v) for k, v in _members.items()}\n',
     "A08", "one test row per position; joined label is display only"),
    ("66_apo_holo_structural_interpretation.py",
     '    "H11": "H9/H10/H11 dimerization region",\n',
     '    "H11": "H9/H10/H11 dimerization region",\n'
     '    "H10/H11": "H9/H10/H11 dimerization region",  # sandbox patch (A10): professor combined label\n',
     "A10", "professor H10/H11 belongs to the dimerization category"),
    ("66_apo_holo_structural_interpretation.py",
     "def categorize(segment: str) -> str:\n    if segment in CATEGORY_MAP:\n",
     "def categorize(segment: str) -> str:\n"
     '    if "|" in segment:  # sandbox patch (A10): dual professor position -> category of each element\n'
     '        return " | ".join(dict.fromkeys(categorize(s) for s in segment.split("|")))\n'
     "    if segment in CATEGORY_MAP:\n",
     "A10", "joined display labels are split, never treated as an element"),
    ("66_apo_holo_structural_interpretation.py",
     '        lambda s: categorize(s) if s else "position not confidently mapped to a named segment"\n',
     '        lambda s: categorize(s) if s else "between professor helices (no professor element; derived position)"  # sandbox patch (A10)\n',
     "A10", "v2 meaning of an unlabelled position"),
    ("67_apo_holo_qc_balance.py", ROOT_LINE, _imports("apo_holo_scope"), "A07", "import S1 helper"),
    ("67_apo_holo_qc_balance.py",
     '    coord_pdbs = coord_pdbs[coord_pdbs["uniprot_id"].isin(eligible_receptors)]\n',
     '    coord_pdbs = apo_holo_scope.apply(coord_pdbs)  # sandbox patch (S1)\n'
     '    coord_pdbs = coord_pdbs[coord_pdbs["uniprot_id"].isin(eligible_receptors)]\n', "A07", "S1"),
    ("67_apo_holo_qc_balance.py",
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n',
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n'
     '    features = apo_holo_scope.apply(features)  # sandbox patch (S1)\n', "A07", "S1"),
    ("68_apo_holo_redundancy_audit.py", ROOT_LINE, _imports("apo_holo_scope"), "A07", "import S1 helper"),
    ("68_apo_holo_redundancy_audit.py",
     '    coord_pdbs = coord.groupby(["uniprot_id", "pdb_id"]).agg(ligand_state=("ligand_state", "first")).reset_index()\n',
     '    coord_pdbs = coord.groupby(["uniprot_id", "pdb_id"]).agg(ligand_state=("ligand_state", "first")).reset_index()\n'
     '    coord_pdbs = apo_holo_scope.apply(coord_pdbs)  # sandbox patch (S1)\n', "A07", "S1"),
    ("69_apo_holo_loo_robustness.py", ROOT_LINE, _imports("apo_holo_scope"), "A07", "import S1 helper"),
    ("69_apo_holo_loo_robustness.py",
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n',
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n'
     '    features = apo_holo_scope.apply(features)  # sandbox patch (S1)\n', "A07", "S1"),
    ("70_apo_holo_figures.py", ROOT_LINE, _imports("apo_holo_scope"), "A07", "import S1 helper"),
    ("70_apo_holo_figures.py",
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n',
     '    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")\n'
     '    features = apo_holo_scope.apply(features)  # sandbox patch (S1)\n', "A07", "S1"),
    # ---- A11/A12/A13/A26 (73)
    ("73_current_functional_site_definitions.py",
     'DIMERIZATION_SEGMENTS = ["H9", "H10", "H11"]',
     'DIMERIZATION_SEGMENTS = ["H9", "H10", "H11", "H10/H11"]  # sandbox patch (A11): professor combined H10/H11',
     "A11", "dimer proxy covers the professor combined element"),
    ("73_current_functional_site_definitions.py",
     '"assignment_confidence": r["assignment_confidence"], "evidence_type": "current_structural_segment_map",',
     '"assignment_confidence": r["assignment_confidence"], "evidence_type": "current_structural_segment_map",\n'
     '            "definition_source": "PROFESSOR_GROUND_TRUTH",  # sandbox patch (S3/A13)',
     "A13", "provenance flag; professor status kept verbatim in assignment_confidence"),
    ("73_current_functional_site_definitions.py",
     '"AF-2/H12 taken directly from the current v3 canonical H12 mapping, never the old heuristic definition."',
     '"AF-2/H12 = professor canonical H12 (final professor table), never the old heuristic definition."',
     "A12", "provenance text"),
    ("73_current_functional_site_definitions.py",
     '"assignment_confidence": confidence, "evidence_type": "geometric_ligand_contact_consensus",',
     '"assignment_confidence": confidence, "evidence_type": "geometric_ligand_contact_consensus",\n'
     '            "definition_source": "GEOMETRIC_LIGAND_CONTACT",  # sandbox patch (S3/A13)',
     "A13", "provenance flag"),
    ("73_current_functional_site_definitions.py",
     '            "evidence_type": "role_d_curated_definition",\n',
     '            "evidence_type": "role_d_curated_definition",\n'
     '            "definition_source": "ROLE_D_CURATED",  # sandbox patch (S3/A13)\n',
     "A13", "provenance flag"),
    ("73_current_functional_site_definitions.py",
     '"assignment_confidence": r["assignment_confidence"], "evidence_type": "structural_segment_correspondence",',
     '"assignment_confidence": r["assignment_confidence"], "evidence_type": "structural_segment_correspondence",\n'
     '            "definition_source": "PROFESSOR_GROUND_TRUTH",  # sandbox patch (S3/A13)',
     "A13", "provenance flag"),
    ("73_current_functional_site_definitions.py",
     '                "n_positions": len(sub),\n',
     '                "n_positions": int(sub["canonical_uniprot_position"].nunique()),  # sandbox patch (A26): unique positions\n'
     '                "n_definition_rows": len(sub),\n',
     "A26", "a dual-label position counts once per site"),
    ("73_current_functional_site_definitions.py",
     'positions (from current v3 mapping)")', 'positions (professor H12)")', "A12", "print text"),
    # ---- A13 pass-through (74, 75)
    ("74_pdb_site_flexibility.py",
     '        .agg(lambda s: s.mode().iloc[0]).to_dict()\n',
     '        .agg(lambda s: s.mode().iloc[0]).to_dict()\n'
     '    # sandbox patch (S3/A13): provenance pass-through; professor statuses stay verbatim\n'
     '    def_src = definitions.groupby(["uniprot_id", "functional_site"])["definition_source"] \\\n'
     '        .agg(lambda s: "|".join(sorted(set(s)))).to_dict()\n'
     '    def_rr = {k: bool(((g["definition_source"] == "PROFESSOR_GROUND_TRUTH") &\n'
     '                       (g["assignment_confidence"] == "review_required")).any())\n'
     '              for k, g in definitions.groupby(["uniprot_id", "functional_site"])}\n',
     "A13", "pass definition_source and the professor review flag through"),
    ("74_pdb_site_flexibility.py",
     '                "site_definition_confidence": def_conf.get((uid, site), ""),\n',
     '                "site_definition_confidence": def_conf.get((uid, site), ""),\n'
     '                "definition_source": def_src.get((uid, site), ""),\n'
     '                "professor_review_required": def_rr.get((uid, site), False),\n',
     "A13", "pass-through"),
    ("75_receptor_site_summary.py",
     '            "definition_level": g["definition_level"].mode().iloc[0],\n',
     '            "definition_level": g["definition_level"].mode().iloc[0],\n'
     '            "definition_source": g["definition_source"].iloc[0],  # sandbox patch (S3/A13)\n'
     '            "professor_review_required": bool(g["professor_review_required"].any()),\n',
     "A13", "pass-through"),
    # ---- A09 (79)
    ("79_atlas_site_enrichment.py", ROOT_LINE, _imports("canonical_helix"), "A09", "import"),
    ("79_atlas_site_enrichment.py",
     '    pos_to_seg = seg_map.set_index(["uniprot_id", "canonical_uniprot_position"])["structural_segment"]\n',
     '    # sandbox patch (A09): membership lookup; each element of a dual position is added\n'
     '    pos_to_seg = canonical_helix.position_memberships(\n'
     '        seg_map[["uniprot_id", "canonical_uniprot_position", "structural_segment"]].itertuples(index=False))\n',
     "A09", "membership lookup"),
    ("79_atlas_site_enrichment.py",
     '        seg = pos_to_seg.get((r["uniprot_id"], r["canonical_uniprot_position"]), "")\n'
     '        if seg:\n'
     '            site_associated_segments.add(seg)\n',
     '        for seg in pos_to_seg.get((r["uniprot_id"], int(r["canonical_uniprot_position"])), ()):\n'
     '            site_associated_segments.add(seg)\n',
     "A09", "never the joined string"),
    # ---- A17 (80)
    ("80_site_robustness.py",
     '    flags["low_confidence_definition"] = flags["site_definition_confidence"].isin(["LOW", "MODERATE"])\n',
     '    # sandbox patch (S3/A17): professor ground truth is not a low-confidence definition\n'
     '    flags["low_confidence_definition"] = flags["site_definition_confidence"].isin(["LOW", "MODERATE"]) & \\\n'
     '        (flags["definition_source"] != "PROFESSOR_GROUND_TRUTH")\n',
     "A17", "explicit"),
    ("80_site_robustness.py",
     '                  "low_coverage", "segment_level_only"]]\n',
     '                  "low_coverage", "segment_level_only", "definition_source", "professor_review_required"]]\n',
     "A17", "separate professor review flag"),
    # ---- A14 (84)
    ("84_definition_quality_stratification.py",
     '        receptor_site[["uniprot_id", "functional_site", "definition_level", "site_definition_confidence"]],\n',
     '        receptor_site[["uniprot_id", "functional_site", "definition_level", "site_definition_confidence",\n'
     '                       "definition_source", "professor_review_required"]],\n',
     "A14", "provenance columns"),
    ("84_definition_quality_stratification.py",
     '        lambda r: classify_tier(r["definition_level"], r["site_definition_confidence"]), axis=1\n',
     '        # sandbox patch (S3/A14): professor ground truth ranks with HIGH; status kept verbatim\n'
     '        lambda r: classify_tier(r["definition_level"], "HIGH" if r["definition_source"] == "PROFESSOR_GROUND_TRUTH"\n'
     '                                else r["site_definition_confidence"]), axis=1\n',
     "A14", "tier rule"),
    ("84_definition_quality_stratification.py",
     '            "site_definition_confidence", "definition_tier"]].to_parquet(tier_lookup_out, index=False)\n',
     '            "site_definition_confidence", "definition_source", "professor_review_required",\n'
     '            "definition_tier"]].to_parquet(tier_lookup_out, index=False)\n',
     "A14", "keep provenance in the lookup"),
    # ---- S2 / A18 / A27 (85)
    ("85_matched_background_tests.py", ROOT_LINE, _imports("matching_context"), "A18", "import"),
    ("85_matched_background_tests.py",
     '    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)\n',
     '    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)\n'
     '    # sandbox patch (S2/A18): tier-1 unit = professor element or derived inter-element stretch\n'
     '    # (matching device only); tier-2 class = independent script-34 local class (HELIX/SHEET/COIL)\n'
     '    seg_map = matching_context.build_matching_frame(\n'
     '        seg_map, pd.read_parquet(ANALYSIS / "receptor_dssp_consensus.parquet"))\n',
     "A18", "approved S2 hierarchy"),
    ("85_matched_background_tests.py",
     '        site_classes = set(site_annotated["structural_class"])\n',
     '        site_classes = set(site_annotated["structural_class"]) - {""}  # sandbox patch (S2): NO_OBSERVATION never matches\n',
     "A18", "NO_OBSERVATION is not a class"),
    ("85_matched_background_tests.py",
     "        if len(tier1) >= MIN_CONTROL_RESIDUES:\n",
     '        if tier1["canonical_uniprot_position"].nunique() >= MIN_CONTROL_RESIDUES:  # sandbox patch (A27)\n',
     "A27", "unique positions"),
    ("85_matched_background_tests.py",
     "        elif len(tier2) >= MIN_CONTROL_RESIDUES:\n",
     '        elif tier2["canonical_uniprot_position"].nunique() >= MIN_CONTROL_RESIDUES:  # sandbox patch (A27)\n',
     "A27", "unique positions"),
    ("85_matched_background_tests.py",
     '            "site_segments": ",".join(sorted(site_segments)), "site_classes": ",".join(sorted(site_classes)),\n',
     '            # sandbox patch (S2): matching strata (elements or derived stretches) and local classes\n'
     '            "site_matching_strata": ",".join(sorted(site_segments)),\n'
     '            "site_local_ss_classes": ",".join(sorted(site_classes)),\n',
     "A18", "column names say what they hold"),
    ("85_matched_background_tests.py",
     "  1. same canonical structural segment, non-site residues of that receptor\n"
     "  2. same structural class (HELIX/SHEET/LOOP/LOOP_UNCERTAIN), non-site\n"
     "     residues of that receptor, any segment\n",
     "  1. same professor element (dual positions belong to both); for residues outside\n"
     "     professor helices the same derived inter-element stretch (matching device only)\n"
     "  2. same independent local class (script-34 consensus HELIX/SHEET/COIL), non-site\n"
     "     residues of that receptor, any stratum  [sandbox patch S2]\n",
     "A18", "docstring states the v2 hierarchy"),
    # ---- A15 (87, 89) + S3 sensitivity naming; A19 (90)
    ("87_superfamily_matched_reanalysis.py", "from pathlib import Path\n", "from pathlib import Path\n" + _s3_suffix_lines(),
     "A15", "sensitivity naming"),
    ("87_superfamily_matched_reanalysis.py",
     '        receptor_site[["uniprot_id", "functional_site", "site_definition_confidence"]],\n',
     '        receptor_site[["uniprot_id", "functional_site", "site_definition_confidence",\n'
     '                       "definition_source", "professor_review_required"]],\n', "A15", "provenance columns"),
    ("87_superfamily_matched_reanalysis.py",
     '    confident = matched[matched["site_definition_confidence"].isin(CONFIDENT_TIERS)]\n',
     '    # sandbox patch (S3/A15): professor ground truth eligible regardless of the legacy filter\n'
     f'    confident = matched[{ELIGIBLE_S3.format(conf="CONFIDENT_TIERS")}]\n'
     '    if _S3:  # labelled SENSITIVITY only: drop definitions resting on review_required professor labels\n'
     '        confident = confident[~confident["professor_review_required"].astype(bool)]\n',
     "A15", "S3 eligibility"),
    ("87_superfamily_matched_reanalysis.py",
     '    out = TABLES / "functional_site_superfamily_matched_reanalysis.csv"\n',
     '    out = TABLES / f"functional_site_superfamily_matched_reanalysis{_SFX}.csv"\n', "A15", "sensitivity naming"),
    ("89_matched_robustness.py", "from pathlib import Path\n", "from pathlib import Path\n" + _s3_suffix_lines(),
     "A15", "sensitivity naming"),
    ("89_matched_robustness.py",
     '    superfamily = pd.read_csv(TABLES / "functional_site_superfamily_matched_reanalysis.csv")\n',
     '    superfamily = pd.read_csv(TABLES / f"functional_site_superfamily_matched_reanalysis{_SFX}.csv")\n',
     "A15", "sensitivity naming"),
    ("89_matched_robustness.py",
     '        receptor_site[["uniprot_id", "functional_site", "site_definition_confidence"]],\n',
     '        receptor_site[["uniprot_id", "functional_site", "site_definition_confidence",\n'
     '                       "definition_source", "professor_review_required"]],\n', "A15", "provenance columns"),
    ("89_matched_robustness.py",
     '    confident = matched[matched["site_definition_confidence"].isin(["HIGH", "MODERATE"])]\n',
     '    # sandbox patch (S3/A15)\n'
     '    confident = matched[' + ELIGIBLE_S3.format(conf=LEGACY_CONFIDENT_LIST) + ']\n'
     '    if _S3:\n'
     '        confident = confident[~confident["professor_review_required"].astype(bool)]\n',
     "A15", "S3 eligibility"),
    ("89_matched_robustness.py",
     '    out1 = TABLES / "functional_site_matched_loo_pdb_robustness.csv"\n',
     '    out1 = TABLES / f"functional_site_matched_loo_pdb_robustness{_SFX}.csv"\n', "A15", "sensitivity naming"),
    ("89_matched_robustness.py",
     '    out2 = TABLES / "functional_site_matched_loo_receptor_robustness.csv"\n',
     '    out2 = TABLES / f"functional_site_matched_loo_receptor_robustness{_SFX}.csv"\n', "A15", "sensitivity naming"),
    ("90_interpretation_classification.py", "from pathlib import Path\n", "from pathlib import Path\n" + _s3_suffix_lines(),
     "A19", "sensitivity naming"),
    ("90_interpretation_classification.py",
     '    superfamily_matched = pd.read_csv(TABLES / "functional_site_superfamily_matched_reanalysis.csv")\n',
     '    superfamily_matched = pd.read_csv(TABLES / f"functional_site_superfamily_matched_reanalysis{_SFX}.csv")\n',
     "A19", "sensitivity naming"),
    ("90_interpretation_classification.py",
     '"Too few HIGH/MODERATE-confidence RESIDUE-level definitions to test the matched-background superfamily reanalysis at all."',
     '"Too few eligible (HIGH/MODERATE or professor ground truth) RESIDUE-level definitions to test the '
     'matched-background superfamily reanalysis at all."', "A19", "wording"),
    ("90_interpretation_classification.py",
     '    out = TABLES / "functional_site_interpretation_classification.csv"\n',
     '    out = TABLES / f"functional_site_interpretation_classification{_SFX}.csv"\n', "A19", "sensitivity naming"),
    # ---- A16 (88)
    ("88_apo_holo_overlap_confidence_recheck.py",
     '                category = "RESIDUE_LEVEL_HIGH_CONFIDENCE"\n',
     '                category = "RESIDUE_LEVEL_HIGH_CONFIDENCE"\n'
     '            elif level == "RESIDUE" and conf_row.iloc[0]["definition_source"] == "PROFESSOR_GROUND_TRUTH":\n'
     '                category = "RESIDUE_LEVEL_PROFESSOR_GROUND_TRUTH"  # sandbox patch (S3/A16); status not relabelled\n',
     "A16", "label"),
    # ---- A23/A26 (94)
    ("94_final_project_numbers.py",
     '    add("dataset", "project_receptors", structures["uniprot_id"].nunique(), "config/nr_metadata.py",\n'
     '        "len(NR_METADATA)", "48 professor-supplied nuclear receptors")\n',
     '    # sandbox patch (metric label): the legacy "project_receptors" value is the number of receptors\n'
     '    # with >= 1 discovery candidate, not the project definition; same count, honest name\n'
     '    add("dataset", "receptors_with_any_discovery_candidate", structures["uniprot_id"].nunique(),\n'
     '        "data/processed/structures.csv", "uniprot_id.nunique()",\n'
     '        "receptors with >= 1 RCSB discovery candidate; NOT the project receptor list (see project_receptors_defined)")\n'
     '    add("dataset", "project_receptors_defined",\n'
     '        len(pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "nr_metadata_validated.csv")),\n'
     '        "data/manifests/nr_metadata_validated.csv", "len(df)", "the professor-supplied project receptor list")\n',
     "LABEL", "misleading legacy metric name corrected; underlying count unchanged"),
    ("94_final_project_numbers.py",
     '    add("structural_mapping", "unresolved_element_H8_receptor_count", 0,\n'
     '        "data/analysis/nr_lbd_structural_segment_map.csv", "structural_segment == \'H8\' count",\n'
     '        "H8 never confidently resolved by the current mapping; never fabricated or zero-filled")\n',
     '    # sandbox patch (A23): computed, not a v1 literal\n'
     '    add("structural_mapping", "receptors_with_professor_H8",\n'
     '        int(named.loc[named["structural_segment"] == "H8", "uniprot_id"].nunique()),\n'
     '        "data/analysis/nr_lbd_structural_segment_map.csv", "uniprot_id.nunique() where structural_segment == \'H8\'")\n',
     "A23", "H8 count computed"),
    ("94_final_project_numbers.py",
     '        f"method={best[\'method\']}, k={int(best[\'k\'])}; this solution is a degenerate 1-vs-28 split, not evidence of real clusters")\n',
     '        f"method={best[\'method\']}, k={int(best[\'k\'])}")  # sandbox patch (A23): v1 interpretation removed\n',
     "A23", "v1 note removed"),
    ("94_final_project_numbers.py",
     '    add("pca_clustering", "supported_conclusion", "no_robust_discrete_clusters_continuum_consistent",\n'
     '        "reports/tables/modeling_summary.md", "Part 11 formal conclusion",\n'
     '        "no method picked because it looked nicest; degenerate best-silhouette solution; DBSCAN finds <=1 cluster")\n',
     "    # sandbox patch (A23): hard-coded v1 conclusion removed (numbers only; see modeling_summary.md)\n",
     "A23", "v1 verdict removed"),
    ("94_final_project_numbers.py",
     '        "reports/tables/subfamily_segment_kruskal.csv", "significant_fdr_05.sum()",\n'
     '        "0 segments differ between groups after BH-FDR correction, in either grouping scheme")\n',
     '        "reports/tables/subfamily_segment_kruskal.csv", "significant_fdr_05.sum()")  # sandbox patch (A23)\n',
     "A23", "v1 note removed"),
    ("94_final_project_numbers.py",
     '        "VDR is INELIGIBLE for apo/holo comparison (0 APO, 48 HOLO); never removed from structures.csv")\n',
     '        f"VDR apo/holo eligible: {bool(vdr[\'eligible_exploratory_3_3\'])} ({int(vdr[\'n_APO\'])} APO, '
     '{int(vdr[\'n_HOLO\'])} HOLO); never removed from structures.csv")  # sandbox patch (A23)\n',
     "A23", "note derived"),
    ("94_final_project_numbers.py",
     '    add("functional_sites", "definition_rows", len(defs), "data/analysis/current_functional_site_definitions.parquet", "len(df)")\n',
     '    add("functional_sites", "definition_rows", len(defs), "data/analysis/current_functional_site_definitions.parquet", "len(df)")\n'
     '    add("functional_sites", "definition_unique_positions",  # sandbox patch (A26)\n'
     '        int(defs[["uniprot_id", "functional_site", "canonical_uniprot_position"]].drop_duplicates().shape[0]),\n'
     '        "data/analysis/current_functional_site_definitions.parquet", "unique (uniprot_id, site, position)")\n',
     "A26", "unique positions next to rows"),
]


# Labelled sensitivity (not the primary policy): functional-site chain without TRACE_LBD entries.
# Applied on top of PATCHES in a separate copy of the v2 sandbox.
TRACE_SENS_STEPS = [s for s in STEPS if "72" <= s[:2] <= "91" and s[:2] not in REPORT_SCRIPTS_REPLACED]
TRACE_PATCHES: list[tuple[str, str, str, str, str]] = [
    ("72_ligand_pocket_geometry.py",
     '    functional = inventory[inventory["nonpolymer_relevance_class"] == "FUNCTIONAL_LBD_LIGAND"].copy()\n',
     '    functional = inventory[inventory["nonpolymer_relevance_class"] == "FUNCTIONAL_LBD_LIGAND"].copy()\n'
     '    from scripts.utils import apo_holo_scope  # sandbox patch (TRACE_LBD site sensitivity)\n'
     '    functional = apo_holo_scope.drop(functional)\n',
     "TRACE_SENS", "no pocket contacts from TRACE_LBD-only entries"),
    ("74_pdb_site_flexibility.py",
     '    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")\n',
     '    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")\n'
     '    import sys as _sys; _sys.path.insert(0, str(PROJECT_ROOT))  # sandbox patch (TRACE_LBD site sensitivity)\n'
     '    from scripts.utils import apo_holo_scope\n'
     '    positions = apo_holo_scope.drop(positions)\n',
     "TRACE_SENS", "no site rows from TRACE_LBD-only entries"),
    ("85_matched_background_tests.py",
     '    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")\n',
     '    positions = pd.read_parquet(ANALYSIS / "pdb_position_bfactor_profiles.parquet")\n'
     '    from scripts.utils import apo_holo_scope  # sandbox patch (TRACE_LBD site sensitivity)\n'
     '    positions = apo_holo_scope.drop(positions)\n',
     "TRACE_SENS", "no matched rows from TRACE_LBD-only entries"),
]


def run_trace_sensitivity() -> dict:
    """Copy the finished v2 sandbox, remove the functional-site chain outputs, apply TRACE_PATCHES
    and re-run steps 72-91 with the same stale-read guard."""
    src, w = work_dir("v2"), work_dir("v2_trace_sens")
    if w.exists():
        shutil.rmtree(w)
    shutil.copytree(src, w, symlinks=True)
    primary = json.loads((src / "logs" / "phase2_run.json").read_text())
    chain = {s["step"] for s in primary["steps"] if s["script"] in TRACE_SENS_STEPS}
    chain_out = {f for s in primary["steps"] if s["step"] in chain for f in s["written"]}
    for f in chain_out:
        (w / f).unlink(missing_ok=True)
    log = []
    for script, old, new, aid, reason in TRACE_PATCHES:
        f = w / "scripts" / script
        text = f.read_text()
        if text.count(old) != 1:
            raise SystemExit(f"trace patch target not found exactly once in {script}")
        f.write_text(text.replace(old, new))
        log.append({"script": script, "adaptation": aid, "reason": reason, "status": "APPLIED"})
    (w / "logs" / "sandbox_patches_trace_sensitivity.json").write_text(json.dumps(log, indent=2))
    foundation = set(json.loads((w / "logs" / "phase2_foundation_snapshot.json").read_text()))
    produced = set(primary["produced"]) - chain_out  # earlier Phase 2 outputs (apo/holo, atlas, ...)
    (w / "logs" / "phase2_io").mkdir(exist_ok=True)
    steps = []
    for script in TRACE_SENS_STEPS:
        tag = f"{script[:2]}_trace_sens"
        print(f"[v2 TRACE_LBD sensitivity] {tag}", flush=True)
        t0 = time.time()
        io = run_step(w, script, tag)
        check_reads(w, io, produced, foundation, tag)
        written = {k for k, m in files_under(w, ["data", "reports"]).items() if m >= t0 - 1}
        produced |= written
        steps.append({"step": tag, "written": sorted(written)})
    out = EXPORT_REPORTS / "tables" / "trace_lbd_site_sensitivity"
    out.mkdir(parents=True, exist_ok=True)
    exported = []
    for s in steps:
        for rel in s["written"]:
            if rel.endswith((".csv", ".parquet")) and rel.startswith(("data/analysis/", "reports/tables/")):
                shutil.copy2(w / rel, out / Path(rel).name)
                exported.append(Path(rel).name)
    (w / "logs" / "phase2_trace_sens_run.json").write_text(json.dumps({"steps": steps}, indent=1))
    return {"n_steps": len(steps), "n_exported_tables": len(set(exported))}


def work_dir(variant: str) -> Path:
    return PROJECT_ROOT / "data" / "interim" / f"rebuild_{variant}"


def apply_patches(w: Path) -> list[dict]:
    log = []
    for script, old, new, aid, reason in PATCHES:
        f = w / "scripts" / script
        text = f.read_text()
        if new in text:
            status = "ALREADY_APPLIED"
        elif text.count(old) == 1:
            f.write_text(text.replace(old, new))
            status = "APPLIED"
        else:
            raise SystemExit(f"patch {aid} target not found exactly once in {script}: {old[:120]!r}")
        log.append({"script": script, "adaptation": aid, "reason": reason, "replaced": old, "with": new, "status": status})
    (w / "logs" / "sandbox_patches_phase2.json").write_text(json.dumps(log, indent=2))
    return log


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def files_under(w: Path, rels) -> dict[str, float]:
    out = {}
    for rel in rels:
        for f in (w / rel).rglob("*"):
            if f.is_file():
                out[str(f.relative_to(w))] = f.stat().st_mtime
    return out


def prepare(w: Path, legacy: bool) -> dict:
    """Snapshot the foundation once, then remove every non-foundation, non-static file
    from the analysis/report folders (v1 copies and outputs of any earlier Phase 2 run)."""
    built = float((w / "logs" / "sandbox_built_at").read_text())
    snap_path = w / "logs" / "phase2_foundation_snapshot.json"
    if not snap_path.exists():
        allf = files_under(w, ["data", "reports"])
        snap = sorted(k for k, m in allf.items() if m > built and not k.startswith("data/raw/"))
        snap_path.write_text(json.dumps(snap, indent=1))
    foundation = set(json.loads(snap_path.read_text()))
    removed = []
    for rel, _ in files_under(w, CLEAN_DIRS).items():
        if rel in foundation or rel in STATIC_INPUTS or (legacy and rel in LEGACY_V1_INPUTS):
            continue
        (w / rel).unlink()
        removed.append(rel)
    for d in ["scripts/utils"]:
        (w / d).mkdir(parents=True, exist_ok=True)
    if not legacy:
        for u in NEW_UTILS:
            shutil.copy2(PROJECT_ROOT / "scripts" / "utils" / u, w / "scripts" / "utils" / u)
        if (PROJECT_ROOT / "scripts" / SUMMARY_SCRIPT).exists():
            shutil.copy2(PROJECT_ROOT / "scripts" / SUMMARY_SCRIPT, w / "scripts" / SUMMARY_SCRIPT)
    io_dir = w / "logs" / "phase2_io"
    if io_dir.exists():
        shutil.rmtree(io_dir)
    io_dir.mkdir(parents=True)
    return {"foundation_files": len(foundation), "removed_stale_or_previous": len(removed)}


HOOK = r'''
import os, sys, json, runpy
_log = open(os.environ["PHASE2_IO_LOG"], "a")
_root = os.getcwd()
_W = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND | os.O_TRUNC
def _hook(event, args):
    if event != "open" or not args:
        return
    p = args[0]
    if isinstance(p, int):
        return
    try:
        p = os.fsdecode(p)
    except TypeError:
        return
    ap = os.path.abspath(p)
    if not ap.startswith(_root + os.sep):
        return
    mode = args[1] if len(args) > 1 else None
    flags = args[2] if len(args) > 2 else 0
    write = (isinstance(mode, str) and any(c in mode for c in "wax+")) or (isinstance(flags, int) and bool(flags & _W))
    if os.path.isdir(ap):
        return
    _log.write(json.dumps({"path": os.path.relpath(ap, _root), "write": write}) + "\n")
    _log.flush()
sys.addaudithook(_hook)
script = sys.argv[1]
sys.argv = [script]
runpy.run_path(script, run_name="__main__")
'''


def run_step(w: Path, script: str, tag: str, extra_env: dict | None = None) -> list[dict]:
    io_log = w / "logs" / "phase2_io" / f"{tag}.jsonl"
    env = {**os.environ, **NO_NETWORK_ENV, "PHASE2_IO_LOG": str(io_log), "MPLBACKEND": "Agg", **(extra_env or {})}
    log = w / "logs" / f"phase2_{tag}.log"
    with open(log, "w") as fh:
        r = subprocess.run([PY, "-c", HOOK, f"scripts/{script}"], cwd=w, env=env, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"{script} failed (exit {r.returncode}); see {log}")
    return [json.loads(l) for l in io_log.read_text().splitlines()] if io_log.exists() else []


def check_reads(w: Path, io: list[dict], produced: set[str], foundation: set[str], tag: str,
                legacy: bool = False) -> list[str]:
    bad = []
    allowed = STATIC_INPUTS | (LEGACY_V1_INPUTS if legacy else set())
    # script 93 is an identity AUDIT that scans every manifest/processed file present; scanning
    # is its purpose, not consuming a result, so those reads are logged but not stopped
    scan = SCAN_PREFIXES if tag.startswith("93") else ()
    for rel in sorted({e["path"] for e in io if not e["write"]}):
        if rel in produced or rel in foundation or rel in allowed or rel.startswith(STATIC_PREFIXES + scan):
            continue
        if rel.startswith("logs/") or "__pycache__" in rel or rel.endswith(".pyc"):
            continue
        bad.append(rel)
    if bad:
        raise SystemExit(f"STALE READ STOP in {tag}: {bad}")
    return bad


def run(variant: str, legacy: bool, until: str | None = None) -> dict:
    w = work_dir(variant)
    if not legacy:
        patches = apply_patches(w)
    prep = prepare(w, legacy)
    foundation = set(json.loads((w / "logs" / "phase2_foundation_snapshot.json").read_text()))
    t_start = time.time()
    produced: set[str] = set()
    steps = []
    plan = []
    for s in STEPS:
        if not legacy and s[:2] in REPORT_SCRIPTS_REPLACED:
            if s.startswith("92"):
                for sx in S3_SENSITIVITY_STEPS:  # labelled S3 sensitivity (never replaces the primary)
                    plan.append((sx, f"{sx[:2]}_s3_sensitivity", {S3_SENSITIVITY_ENV: "1"}))
                plan.append((SUMMARY_SCRIPT, "105_summaries", None))  # all five summaries, after the last table
            continue
        plan.append((s, s[:-3], None))
    for script, tag, env in plan:
        if until and tag.startswith(until):
            print(f"stopping before {tag} (--until)")
            break
        print(f"[{variant}{' legacy' if legacy else ''}] {tag}", flush=True)
        t0 = time.time()
        io = run_step(w, script, tag, env)
        check_reads(w, io, produced, foundation, tag, legacy)
        written = {k for k, m in files_under(w, ["data", "reports"]).items() if m >= t0 - 1}
        written |= {e["path"] for e in io if e["write"]}
        written = {p for p in written if (w / p).exists() and not p.startswith("logs/")}
        produced |= written
        steps.append({"step": tag, "script": script, "seconds": round(time.time() - t0, 1),
                      "n_reads": len({e['path'] for e in io if not e['write']}), "written": sorted(written)})
    summary = {"variant": variant, "legacy_unpatched": legacy, "prepare": prep,
               "n_patches": 0 if legacy else len(patches), "started": t_start, "steps": steps,
               "produced": sorted(produced)}
    (w / "logs" / "phase2_run.json").write_text(json.dumps(summary, indent=1))
    return summary


def compare_with_frozen(w: Path) -> dict:
    """legacy run on v1 inputs vs frozen tracked outputs (same relative paths)."""
    run_info = json.loads((w / "logs" / "phase2_run.json").read_text())
    tracked = set(subprocess.run(["git", "ls-files"], cwd=PROJECT_ROOT, capture_output=True, text=True,
                                 check=True).stdout.split())
    rows = []
    for rel in run_info["produced"]:
        if rel not in tracked:
            rows.append({"file": rel, "status": "NOT_TRACKED_IN_RELEASE"})
            continue
        same = sha256(w / rel) == sha256(PROJECT_ROOT / rel)
        rows.append({"file": rel, "status": "IDENTICAL" if same else "DIFFERENT"})
    out = {"n_files": len(rows), "counts": {s: sum(r["status"] == s for r in rows) for s in
                                             ["IDENTICAL", "DIFFERENT", "NOT_TRACKED_IN_RELEASE"]},
           "files": rows}
    (w / "logs" / "phase2_v1_compare.json").write_text(json.dumps(out, indent=1))
    return out


def export(variant: str) -> dict:
    w = work_dir(variant)
    run_info = json.loads((w / "logs" / "phase2_run.json").read_text())
    dst_for = {"data/analysis/": EXPORT_DATA / "analysis", "reports/tables/": EXPORT_REPORTS / "tables",
               "reports/figures/": EXPORT_REPORTS / "figures", "data/manifests/": EXPORT_DATA / "manifests",
               "data/processed/": EXPORT_DATA / "processed"}
    exported = []
    for rel in run_info["produced"]:
        for pre, root in dst_for.items():
            if rel.startswith(pre):
                dst = root / rel[len(pre):]
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(w / rel, dst)
                exported.append(str(dst.relative_to(PROJECT_ROOT)))
                break
        else:
            raise SystemExit(f"no export location for {rel}")
    # A28: the independent local-structure consensus used by script 85
    src = w / "data" / "analysis" / "receptor_dssp_consensus.parquet"
    shutil.copy2(src, EXPORT_DATA / "analysis" / src.name)
    exported.append(f"data/rebuild_v2/analysis/{src.name}")
    shutil.copy2(w / "logs" / "sandbox_patches_phase2.json", EXPORT_REPORTS / "sandbox_patches_phase2.json")
    (EXPORT_REPORTS / "phase2_export_manifest.json").write_text(json.dumps(
        {"variant": variant, "exported_files": sorted(exported),
         "steps": [{k: v for k, v in s.items() if k != "written"} for s in run_info["steps"]]}, indent=2) + "\n")
    return {"n_exported": len(exported)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--variant", choices=["v2", "v1_repro"], required=True)
    ap.add_argument("--legacy", action="store_true", help="unpatched legacy scripts (fidelity check)")
    ap.add_argument("--compare", action="store_true", help="compare a legacy v1 run with the frozen outputs")
    ap.add_argument("--export", action="store_true")
    ap.add_argument("--export-only", action="store_true")
    ap.add_argument("--trace-sensitivity", action="store_true",
                    help="labelled sensitivity: functional-site chain 72-91 without TRACE_LBD entries")
    ap.add_argument("--until", default=None, help="stop before the step whose tag starts with this")
    args = ap.parse_args()
    if args.trace_sensitivity:
        print(json.dumps(run_trace_sensitivity()))
        return
    if args.variant == "v2" and args.legacy:
        raise SystemExit("the v2 run is always patched")
    if not args.export_only:
        s = run(args.variant, args.legacy, args.until)
        print(json.dumps({k: v for k, v in s.items() if k not in ("steps", "produced")}, indent=2))
    if args.compare:
        c = compare_with_frozen(work_dir(args.variant))
        print(json.dumps(c["counts"]))
    if (args.export or args.export_only) and args.variant == "v2":
        print(json.dumps(export(args.variant)))


if __name__ == "__main__":
    main()
