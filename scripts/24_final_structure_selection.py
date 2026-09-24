"""
Final Data Release, Step 1-4: freeze the final structure-selection policy,
select exactly one representative instance per (project receptor x PDB
entry x receptor polymer entity), build the complete 2072-row selection
audit lineage, and compute primary/sensitivity set membership flags.

Uses only already-frozen/approved manifests (Stage 2 candidate inventory,
Stage 3A precoordinate QC, Stage 3B corrected instance QC + observations
parquet, Stage 2.1 taxonomy and multi-UniProt fusion audits). Performs NO
new network access and NO mmCIF re-parsing.

Writes:
    data/manifests/final_representative_instance_selection.csv
    data/manifests/final_structure_selection_audit.csv
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.missing_run_analysis import analyze_missing_runs  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
INTERIM_QC_DIR = PROJECT_ROOT / "data" / "interim" / "qc"

POOL_CSV = MANIFESTS_DIR / "stage3b_coordinate_pool.csv"
PRECOORDINATE_QC_CSV = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"
INSTANCE_QC_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_instance_qc.csv"
OBSERVATIONS_PARQUET = INTERIM_QC_DIR / "stage3b_lbd_ca_observations.parquet"
LBD_REFERENCE_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"
MULTI_UNIPROT_AUDIT_CSV = MANIFESTS_DIR / "rcsb_multi_uniprot_entity_audit.csv"
TAXONOMY_AUDIT_CSV = MANIFESTS_DIR / "rcsb_source_taxonomy_audit.csv"

REP_SELECTION_CSV = MANIFESTS_DIR / "final_representative_instance_selection.csv"
SELECTION_AUDIT_CSV = MANIFESTS_DIR / "final_structure_selection_audit.csv"

SELECTION_POLICY_VERSION = "final-structure-selection-v1"

PRIMARY_COMPLETENESS_THRESHOLD = 0.90

# --- Section E: explicit, individually-documented identity audit ---------
# Only the 9 "unusual taxonomy" candidates (source_includes_homo_sapiens_9606
# == False) need an override; every other candidate is HUMAN_OR_HUMAN_DERIVED
# by construction (its deposited source organism list literally includes
# Homo sapiens). Keyed by (pdb_id, polymer_entity_id); each entry documents
# the specific evidence used, per the Stage 2.1 taxonomy audit
# (data/manifests/rcsb_source_taxonomy_audit.csv).
IDENTITY_AUDIT_OVERRIDES: dict[tuple[str, str], tuple[str, str]] = {
    ("3ZQT", "3ZQT_1"): (
        "NON_HUMAN_CONSTRUCT",
        "Deposited source organism is Escherichia coli with no synthetic-construct or "
        "ancestral-reconstruction annotation; Stage 2.1 taxonomy audit (EXPLAINED) confirmed "
        "this is genuinely non-human, not an expression-host misattribution.",
    ),
    ("5CJ6", "5CJ6_2"): (
        "HUMAN_OR_HUMAN_DERIVED",
        "Source organism is 'synthetic construct' (taxid 32630); UniProt mapping is confirmed "
        "and unique to the human receptor accession (P10275) — a synthetically produced "
        "peptide identical to the human sequence, not a different-organism source.",
    ),
    ("3CBM", "3CBM_2"): (
        "HUMAN_OR_HUMAN_DERIVED",
        "Same synthetic-construct rationale as 5CJ6_2: UniProt mapping confirmed and unique "
        "to the human receptor accession (P03372).",
    ),
    ("8ARW", "8ARW_2"): (
        "AMBIGUOUS_IDENTITY",
        "Stage 2.1 taxonomy audit flagged REVIEW_NEEDED: source organism Escherichia coli "
        "does not match a recognized expression-host artifact pattern and was not further "
        "resolved; genuine biological identity remains unconfirmed.",
    ),
    ("5CBX", "5CBX_1"): (
        "ANCESTRAL_RECONSTRUCTION",
        "Entity description 'AncGR DNA Binding Domain' plus unclassified source taxid is the "
        "established naming/annotation pattern for a computationally reconstructed ancestral "
        "glucocorticoid-receptor DBD, not an extant-organism source.",
    ),
    ("5CBY", "5CBY_1"): (
        "ANCESTRAL_RECONSTRUCTION",
        "Entity description 'AncGR2 DNA Binding Domain' — same ancestral-reconstruction "
        "pattern as 5CBX_1.",
    ),
    ("5CBZ", "5CBZ_1"): (
        "ANCESTRAL_RECONSTRUCTION",
        "Entity description 'AncMR DNA Binding Domain' — same ancestral-reconstruction "
        "pattern as 5CBX_1/5CBY_1 (ancestral mineralocorticoid-receptor DBD).",
    ),
    ("5CC0", "5CC0_1"): (
        "ANCESTRAL_RECONSTRUCTION",
        "Entity description 'AncSR2 DNA Binding Domain' indicates a computationally "
        "reconstructed ancestral steroid-receptor DBD despite the 'synthetic construct' "
        "source-organism label; the ancestral-reconstruction naming takes precedence.",
    ),
    ("1VJB", "1VJB_1"): (
        "NON_HUMAN_ORTHOLOG",
        "Source organism Mus musculus (taxid 10090), confirmed by Stage 2.1 taxonomy audit as "
        "genuinely non-human — this is the mouse ortholog of ERRgamma, not a human-derived "
        "construct.",
    ),
}

EXCLUSION_REASON_VOCAB = {
    "NON_XRAY", "NO_LBD_OVERLAP", "LBD_COVERAGE_LT_0_90",
    "RESOLUTION_LT_1_8", "RESOLUTION_GT_3_5", "RESOLUTION_MISSING",
    "RFREE_GT_0_30", "RFREE_MISSING",
    "NON_HUMAN_ORTHOLOG", "ANCESTRAL_RECONSTRUCTION", "AMBIGUOUS_IDENTITY", "NON_HUMAN_CONSTRUCT",
}


def classify_identity(pdb_id: str, polymer_entity_id: str, source_includes_human: bool) -> tuple[str, str]:
    if source_includes_human:
        return "HUMAN_OR_HUMAN_DERIVED", "Deposited source organism list includes Homo sapiens (9606)."
    override = IDENTITY_AUDIT_OVERRIDES.get((pdb_id, polymer_entity_id))
    if override is None:
        raise ValueError(
            f"Unaudited unusual-taxonomy candidate {pdb_id}/{polymer_entity_id}: "
            "source_includes_homo_sapiens_9606 is False but no explicit identity audit entry "
            "exists. Every unusual-taxonomy candidate must be individually audited (Section E) "
            "before this script can run — add an entry to IDENTITY_AUDIT_OVERRIDES."
        )
    return override


def compute_instance_run_stats(obs: pd.DataFrame, lbd_lengths: dict[str, int]) -> pd.DataFrame:
    rows = []
    for instance_id, group in obs.groupby("instance_id", sort=False):
        uniprot_id = group["uniprot_id"].iloc[0]
        lbd_length = lbd_lengths[uniprot_id]
        stats = analyze_missing_runs(group, lbd_length)
        stats["instance_id"] = instance_id
        rows.append(stats)
    return pd.DataFrame(rows)


TIE_BREAK_STEPS = [
    ("highest positive-occupancy full standardized-LBD Ca coverage", "observed_full_lbd_ca_coverage", False),
    ("highest B-factor-usable full standardized-LBD coverage", "bfactor_usable_full_lbd_coverage", False),
    ("fewest internal mapped-but-unmodeled standardized-LBD residues", "internal_missing_residue_count", True),
    ("fewest zero-occupancy-only standardized-LBD residues", "mapped_with_zero_occupancy_only_ca_residues", True),
    ("fewest total mapped-but-unmodeled standardized-LBD residues", "total_mapped_unmodeled_residues", True),
    ("lexicographically smallest label_asym_id", "label_asym_id", True),
    ("lexicographically smallest auth_asym_id", "auth_asym_id", True),
]


def select_representative(candidate_instances: pd.DataFrame) -> tuple[pd.Series, str]:
    """Apply the predeclared deterministic tie-break cascade. Returns the
    selected row plus a human-readable description of the step that broke
    the tie (or 'NO_TIE_SINGLE_INSTANCE' / 'NO_TIE_FIRST_CRITERION')."""
    remaining = candidate_instances.copy()
    tie_break_step = "NO_TIE_SINGLE_INSTANCE" if len(remaining) == 1 else None

    for step_idx, (description, column, ascending) in enumerate(TIE_BREAK_STEPS):
        if len(remaining) == 1:
            break
        best = remaining[column].min() if ascending else remaining[column].max()
        narrowed = remaining[remaining[column] == best]
        if len(narrowed) < len(remaining) and tie_break_step is None:
            tie_break_step = description
        remaining = narrowed

    if tie_break_step is None:
        tie_break_step = "NO_TIE_FIRST_CRITERION" if len(candidate_instances) > 1 else "NO_TIE_SINGLE_INSTANCE"

    selected = remaining.iloc[0]
    return selected, tie_break_step


def build_representative_selection(instance_qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (uniprot_id, pdb_id, polymer_entity_id), group in instance_qc.groupby(
        ["uniprot_id", "pdb_id", "polymer_entity_id"], sort=False
    ):
        selected, tie_break_step = select_representative(group)
        rows.append(
            {
                "uniprot_id": uniprot_id, "nr_code": selected["nr_code"], "common_name": selected["common_name"],
                "pdb_id": pdb_id, "polymer_entity_id": polymer_entity_id, "entity_id": selected["entity_id"],
                "candidate_instance_count": len(group),
                "all_candidate_instance_ids_json": json.dumps(sorted(group["instance_id"].tolist())),
                "all_candidate_qc_json": json.dumps(
                    group[
                        [
                            "instance_id", "label_asym_id", "auth_asym_id",
                            "observed_full_lbd_ca_coverage", "bfactor_usable_full_lbd_coverage",
                            "internal_missing_residue_count", "mapped_with_zero_occupancy_only_ca_residues",
                            "total_mapped_unmodeled_residues",
                        ]
                    ].to_dict(orient="records")
                ),
                "selected_instance_id": selected["instance_id"],
                "selected_label_asym_id": selected["label_asym_id"],
                "selected_auth_asym_id": selected["auth_asym_id"],
                "selected_model_num": selected["model_num"],
                "selected_observed_full_lbd_ca_coverage": selected["observed_full_lbd_ca_coverage"],
                "selected_bfactor_usable_full_lbd_coverage": selected["bfactor_usable_full_lbd_coverage"],
                "tie_break_step": tie_break_step,
                "selection_reason": (
                    f"{len(group)} candidate instance(s) for this receptor x PDB entry x "
                    f"polymer entity; selected {selected['instance_id']} via: {tie_break_step}."
                ),
            }
        )
    return pd.DataFrame(rows)


def build_selection_audit(
    pool: pd.DataFrame,
    precoordinate_qc: pd.DataFrame,
    instance_qc: pd.DataFrame,
    rep_selection: pd.DataFrame,
    multi_uniprot_ids: set[str],
    inventory_human: pd.DataFrame,
) -> pd.DataFrame:
    pq_indexed = precoordinate_qc.set_index(["uniprot_id", "polymer_entity_id"])
    rep_indexed = rep_selection.set_index(["uniprot_id", "pdb_id", "polymer_entity_id"])
    human_map = inventory_human.set_index(["pdb_id", "polymer_entity_id"])["source_includes_homo_sapiens_9606"].to_dict()

    rows = []
    for _, prow in pool.iterrows():
        key_pq = (prow["uniprot_id"], prow["polymer_entity_id"])
        pq = pq_indexed.loc[key_pq] if key_pq in pq_indexed.index else None

        key_rep = (prow["uniprot_id"], prow["pdb_id"], prow["polymer_entity_id"])
        has_rep = key_rep in rep_indexed.index
        rep = rep_indexed.loc[key_rep] if has_rep else None

        source_includes_human = bool(human_map.get((prow["pdb_id"], prow["polymer_entity_id"]), True))
        identity_class, identity_reason = classify_identity(
            prow["pdb_id"], prow["polymer_entity_id"], source_includes_human
        )
        fusion_or_chimera = prow["polymer_entity_id"] in multi_uniprot_ids

        experimental_methods = pq["experimental_methods"] if pq is not None else ""
        resolution = pq["resolution_max"] if pq is not None else None
        r_free = pq["r_free_max"] if pq is not None else None
        resolution_status = pq["resolution_project_range_status"] if pq is not None else "MISSING"
        rfree_status = pq["r_free_status"] if pq is not None else "MISSING"
        lbd_start = pq["lbd_start"] if pq is not None else None
        lbd_end = pq["lbd_end"] if pq is not None else None
        lbd_length = pq["lbd_length"] if pq is not None else None

        if has_rep:
            positive_occ_cov = float(rep["selected_observed_full_lbd_ca_coverage"])
            bfactor_cov = float(rep["selected_bfactor_usable_full_lbd_coverage"])
            representative_instance_id = rep["selected_instance_id"]
            label_asym_id = rep["selected_label_asym_id"]
            auth_asym_id = rep["selected_auth_asym_id"]
        else:
            positive_occ_cov = None
            bfactor_cov = None
            representative_instance_id = ""
            label_asym_id = ""
            auth_asym_id = ""

        completeness_ge_090 = (positive_occ_cov is not None) and (positive_occ_cov >= 0.90)
        completeness_ge_095 = (positive_occ_cov is not None) and (positive_occ_cov >= 0.95)
        completeness_ge_099 = (positive_occ_cov is not None) and (positive_occ_cov >= 0.99)
        completeness_ge_080 = (positive_occ_cov is not None) and (positive_occ_cov >= 0.80)
        completeness_eq_100 = (positive_occ_cov is not None) and (positive_occ_cov >= 1.00)

        primary_method_pass = bool(prow["method_is_xray"])
        primary_completeness_pass = bool(completeness_ge_090)
        primary_resolution_pass = resolution_status == "PASS_PROJECT_RANGE"
        primary_rfree_pass = rfree_status == "PASS_ALL_LE_0_30"
        primary_identity_pass = identity_class == "HUMAN_OR_HUMAN_DERIVED"

        reasons: list[str] = []
        if not has_rep:
            if not prow["method_is_xray"]:
                reasons.append("NON_XRAY")
            if not prow["mapping_overlaps_lbd"]:
                reasons.append("NO_LBD_OVERLAP")
            primary_selection_status = "NO_COORDINATE_CANDIDATE"
        else:
            if not primary_method_pass:
                reasons.append("NON_XRAY")
            if not primary_completeness_pass:
                reasons.append("LBD_COVERAGE_LT_0_90")
            if resolution_status == "ULTRAHIGH_LT_1_8":
                reasons.append("RESOLUTION_LT_1_8")
            elif resolution_status == "TOO_LOW_RESOLUTION_GT_3_5":
                reasons.append("RESOLUTION_GT_3_5")
            elif resolution_status == "MISSING":
                reasons.append("RESOLUTION_MISSING")
            if rfree_status == "FAIL_ALL_GT_0_30":
                reasons.append("RFREE_GT_0_30")
            elif rfree_status == "MISSING":
                reasons.append("RFREE_MISSING")
            if identity_class in ("NON_HUMAN_ORTHOLOG", "ANCESTRAL_RECONSTRUCTION", "AMBIGUOUS_IDENTITY", "NON_HUMAN_CONSTRUCT"):
                reasons.append(identity_class)

            primary_selection_status = "SELECTED_PRIMARY" if not reasons else "EXCLUDED_PRIMARY"

        assert set(reasons) <= EXCLUSION_REASON_VOCAB, f"unexpected reason code(s): {reasons}"

        rows.append(
            {
                "uniprot_id": prow["uniprot_id"], "nr_code": prow["nr_code"], "common_name": prow["common_name"],
                "group": prow["group"],
                "pdb_id": prow["pdb_id"], "polymer_entity_id": prow["polymer_entity_id"], "entity_id": prow["entity_id"],
                "stage3b_coordinate_selected": bool(prow["coordinate_audit_selected"]),
                "representative_instance_id": representative_instance_id,
                "label_asym_id": label_asym_id, "auth_asym_id": auth_asym_id,
                "experimental_methods": experimental_methods, "resolution": resolution, "r_free": r_free,
                "standardized_lbd_start": lbd_start, "standardized_lbd_end": lbd_end, "standardized_lbd_length": lbd_length,
                "positive_occupancy_lbd_coverage": positive_occ_cov,
                "bfactor_usable_lbd_coverage": bfactor_cov,
                "completeness_ge_090": completeness_ge_090, "completeness_ge_095": completeness_ge_095,
                "completeness_ge_099": completeness_ge_099, "completeness_ge_080": completeness_ge_080,
                "completeness_eq_100": completeness_eq_100,
                "identity_class": identity_class,
                "fusion_or_chimera": fusion_or_chimera,
                "primary_method_pass": primary_method_pass,
                "primary_completeness_pass": primary_completeness_pass,
                "primary_resolution_pass": primary_resolution_pass,
                "primary_rfree_pass": primary_rfree_pass,
                "primary_identity_pass": primary_identity_pass,
                "primary_selection_status": primary_selection_status,
                "primary_exclusion_reasons_json": json.dumps(reasons),
                "in_primary_set": primary_selection_status == "SELECTED_PRIMARY",
                "in_completeness_95_sensitivity": bool(has_rep and completeness_ge_095 and primary_method_pass and primary_resolution_pass and primary_rfree_pass and primary_identity_pass),
                "in_completeness_99_sensitivity": bool(has_rep and completeness_ge_099 and primary_method_pass and primary_resolution_pass and primary_rfree_pass and primary_identity_pass),
                "in_completeness_80_sensitivity": bool(has_rep and completeness_ge_080 and primary_method_pass and primary_resolution_pass and primary_rfree_pass and primary_identity_pass),
                "in_ultrahigh_resolution_sensitivity": bool(has_rep and resolution_status == "ULTRAHIGH_LT_1_8" and primary_completeness_pass and primary_method_pass and primary_rfree_pass and primary_identity_pass),
                "in_no_fusion_sensitivity": bool(primary_selection_status == "SELECTED_PRIMARY" and not fusion_or_chimera),
                "selection_policy_version": SELECTION_POLICY_VERSION,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    print("Loading frozen manifests (no recomputation of Stage 2/3A/3B results)...")
    pool = pd.read_csv(POOL_CSV, dtype={"pdb_id": str})
    precoordinate_qc = pd.read_csv(PRECOORDINATE_QC_CSV, dtype={"pdb_id": str})
    instance_qc = pd.read_csv(INSTANCE_QC_CSV, dtype={"pdb_id": str})
    obs = pd.read_parquet(OBSERVATIONS_PARQUET)
    lbd_reference = pd.read_csv(LBD_REFERENCE_CSV)
    multi_uniprot_audit = pd.read_csv(MULTI_UNIPROT_AUDIT_CSV, dtype={"pdb_id": str})
    inventory_human = pd.read_csv(
        PROJECT_ROOT / "data" / "manifests" / "rcsb_candidate_inventory.csv", dtype={"pdb_id": str}
    )[["pdb_id", "polymer_entity_id", "source_includes_homo_sapiens_9606"]].drop_duplicates(
        subset=["pdb_id", "polymer_entity_id"]
    )

    lbd_lengths = dict(zip(lbd_reference["uniprot_id"], lbd_reference["lbd_length"]))

    print("Computing internal-vs-terminal missing-run statistics for all 3159 measured instances...")
    run_stats = compute_instance_run_stats(obs, lbd_lengths)
    instance_qc = instance_qc.merge(run_stats, on="instance_id", how="left")
    instance_qc["total_mapped_unmodeled_residues"] = (
        instance_qc["mapped_lbd_residues"] - instance_qc["positive_occupancy_lbd_ca_residues"]
    )
    assert instance_qc["internal_missing_residue_count"].notna().all()

    print("Selecting representative instance per (receptor x PDB entry x polymer entity)...")
    rep_selection = build_representative_selection(instance_qc)
    rep_selection.to_csv(REP_SELECTION_CSV, index=False)
    print(f"Wrote {REP_SELECTION_CSV.relative_to(PROJECT_ROOT)} ({len(rep_selection)} rows)")
    print(f"  {(rep_selection['candidate_instance_count'] > 1).sum()} candidates required tie-break selection among >1 instance")
    print(f"  tie-break step distribution:\n{rep_selection['tie_break_step'].value_counts().to_string()}")

    multi_uniprot_ids = set(multi_uniprot_audit["polymer_entity_id"])

    print("\nBuilding full 2072-row final structure selection audit...")
    audit = build_selection_audit(pool, precoordinate_qc, instance_qc, rep_selection, multi_uniprot_ids, inventory_human)
    assert len(audit) == 2072
    audit.to_csv(SELECTION_AUDIT_CSV, index=False)
    print(f"Wrote {SELECTION_AUDIT_CSV.relative_to(PROJECT_ROOT)} ({len(audit)} rows)")
    print(audit["primary_selection_status"].value_counts().to_string())
    print()
    print("Primary set: ", audit["in_primary_set"].sum(), "structures,",
          audit.loc[audit["in_primary_set"], "pdb_id"].nunique(), "unique PDB entries,",
          audit.loc[audit["in_primary_set"], "uniprot_id"].nunique(), "receptors")


if __name__ == "__main__":
    main()
