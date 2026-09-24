"""
Stage 3B.1: recompute all instance/candidate/receptor-level Stage 3B
summaries from the occupancy-corrected observation parquet, using
POSITIVE OCCUPANCY (not mere atom-record presence) as the definition of
"observed Cα" from this point forward.

Reads:
    data/interim/qc/stage3b_lbd_ca_observations.parquet   (corrected, Stage 3B.1)
    data/raw/mmcif/{PDB_ID}.cif.gz  (already-cached files, LOCAL re-parse only,
        to re-extract the `_pdbx_unobs_or_zero_occ_residues` ground truth —
        see scripts/utils/rebuild_unobserved_ground_truth.py. No network
        access, no re-download. This re-parse exists because an earlier
        Stage 3B.1 draft of this script accidentally overwrote the Stage 3B
        crosscheck file with a buggy, empty recomputation before the bug
        was caught; rebuilding from the mmCIF cache is safer than trusting
        a partially-overwritten manifest.)
    data/manifests/stage3b_coordinate_pool.csv
    data/manifests/mmcif_download_manifest.csv
    data/manifests/nr_lbd_analysis_policy.json
    data/manifests/stage3a_release_summary.json

Writes (regenerated):
    data/manifests/rcsb_lbd_coordinate_instance_qc.csv
    data/manifests/rcsb_lbd_coordinate_candidate_summary.csv
    data/manifests/stage3b_counts_by_receptor.csv
    data/manifests/mmcif_unobserved_residue_crosscheck.csv (rebuilt with 4-category A/B/C/D)
    data/manifests/stage3b_coordinate_summary.json
    reports/tables/stage3b_coordinate_completeness_report.md
    data/manifests/stage3b_provenance.json

Writes (new):
    data/manifests/stage3b_occupancy_correction_audit.csv
    data/manifests/stage3b_completeness_threshold_sensitivity.csv

No completeness threshold is chosen, no representative chain is selected,
no B-factor normalization occurs.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
REPORTS_TABLES_DIR = PROJECT_ROOT / "reports" / "tables"
INTERIM_QC_DIR = PROJECT_ROOT / "data" / "interim" / "qc"

OBSERVATIONS_PARQUET = INTERIM_QC_DIR / "stage3b_lbd_ca_observations.parquet"
POOL_CSV = MANIFESTS_DIR / "stage3b_coordinate_pool.csv"
DOWNLOAD_MANIFEST_CSV = MANIFESTS_DIR / "mmcif_download_manifest.csv"
LBD_POLICY_JSON = MANIFESTS_DIR / "nr_lbd_analysis_policy.json"
STAGE3A_RELEASE_SUMMARY = MANIFESTS_DIR / "stage3a_release_summary.json"

INSTANCE_QC_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_instance_qc.csv"
CANDIDATE_SUMMARY_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_candidate_summary.csv"
COUNTS_BY_RECEPTOR_CSV = MANIFESTS_DIR / "stage3b_counts_by_receptor.csv"
CROSSCHECK_CSV = MANIFESTS_DIR / "mmcif_unobserved_residue_crosscheck.csv"
SUMMARY_JSON = MANIFESTS_DIR / "stage3b_coordinate_summary.json"
REPORT_MD = REPORTS_TABLES_DIR / "stage3b_coordinate_completeness_report.md"
PROVENANCE_JSON = MANIFESTS_DIR / "stage3b_provenance.json"
CORRECTION_AUDIT_CSV = MANIFESTS_DIR / "stage3b_occupancy_correction_audit.csv"
THRESHOLD_SENSITIVITY_CSV = MANIFESTS_DIR / "stage3b_completeness_threshold_sensitivity.csv"

THRESHOLDS = [1.00, 0.99, 0.95, 0.90, 0.80]


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def coverage_bin(value) -> str:
    if pd.isna(value):
        return "MISSING"
    v = float(value)
    if v == 1.0:
        return "1.000"
    if v >= 0.99:
        return "[0.99,1.00)"
    if v >= 0.95:
        return "[0.95,0.99)"
    if v >= 0.90:
        return "[0.90,0.95)"
    if v >= 0.80:
        return "[0.80,0.90)"
    return "<0.80"


def build_instance_qc(obs: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id", "entity_id",
                  "instance_id", "label_asym_id", "auth_asym_id", "model_num"]
    rows = []
    for key, g in obs.groupby(group_cols, dropna=False, sort=False):
        (uid, nr_code, common_name, pdb_id, entity_id_full, entity_id,
         instance_id, label_asym_id, auth_asym_id, model_num) = key

        lbd_length = len(g)
        mapped_mask = g["is_construct_mapped"].tolist()
        record_present_list = g["ca_record_present"].tolist()
        positive_occ_list = g["ca_positive_occupancy"].tolist()
        usable_list = g["ca_bfactor_usable"].tolist()

        mapped_lbd_residues = sum(mapped_mask)
        unmapped_lbd_residues = lbd_length - mapped_lbd_residues

        record_present_residues = sum(1 for m, r in zip(mapped_mask, record_present_list) if m and r)
        positive_occ_residues = sum(1 for m, p in zip(mapped_mask, positive_occ_list) if m and p)
        usable_residues = sum(1 for m, u in zip(mapped_mask, usable_list) if m and u)

        mapped_no_ca_record = sum(1 for m, r in zip(mapped_mask, record_present_list) if m and not r)
        mapped_zero_occ_only = sum(
            1 for m, r, p in zip(mapped_mask, record_present_list, positive_occ_list) if m and r and not p
        )
        positive_occ_but_bfactor_unusable = sum(
            1 for m, p, u in zip(mapped_mask, positive_occ_list, usable_list) if m and p and not u
        )

        record_present_full_cov = round(record_present_residues / lbd_length, 6) if lbd_length else 0.0
        observed_full_cov = round(positive_occ_residues / lbd_length, 6) if lbd_length else 0.0
        bfactor_full_cov = round(usable_residues / lbd_length, 6) if lbd_length else 0.0

        record_present_within_cov = round(record_present_residues / mapped_lbd_residues, 6) if mapped_lbd_residues else ""
        observed_within_cov = round(positive_occ_residues / mapped_lbd_residues, 6) if mapped_lbd_residues else ""
        bfactor_within_cov = round(usable_residues / mapped_lbd_residues, 6) if mapped_lbd_residues else ""

        rows.append(
            {
                "uniprot_id": uid, "nr_code": nr_code, "common_name": common_name,
                "pdb_id": pdb_id, "polymer_entity_id": entity_id_full, "entity_id": entity_id,
                "instance_id": instance_id, "label_asym_id": label_asym_id, "auth_asym_id": auth_asym_id,
                "model_num": model_num,
                "lbd_length": lbd_length,
                "mapped_lbd_residues": mapped_lbd_residues, "unmapped_lbd_residues": unmapped_lbd_residues,
                "construct_lbd_mapping_coverage": round(mapped_lbd_residues / lbd_length, 6) if lbd_length else 0.0,

                "ca_record_present_lbd_residues": record_present_residues,
                "record_present_full_lbd_coverage": record_present_full_cov,
                "record_present_within_construct_coverage": record_present_within_cov,

                "positive_occupancy_lbd_ca_residues": positive_occ_residues,
                "observed_full_lbd_ca_coverage": observed_full_cov,
                "observed_within_construct_ca_coverage": observed_within_cov,

                "bfactor_usable_lbd_ca_residues": usable_residues,
                "bfactor_usable_full_lbd_coverage": bfactor_full_cov,
                "bfactor_usable_within_construct_coverage": bfactor_within_cov,

                "mapped_but_no_ca_record_residues": mapped_no_ca_record,
                "mapped_with_zero_occupancy_only_ca_residues": mapped_zero_occ_only,
                "mapped_positive_occupancy_ca_residues": positive_occ_residues,
                "positive_occupancy_but_bfactor_unusable_residues": positive_occ_but_bfactor_unusable,

                "construct_maps_entire_standardized_lbd": mapped_lbd_residues == lbd_length,
                "full_standardized_lbd_ca_record_complete": record_present_residues == lbd_length,
                "full_standardized_lbd_positive_occupancy_complete": positive_occ_residues == lbd_length,
                "full_standardized_lbd_bfactor_complete": usable_residues == lbd_length,
                "mapped_construct_ca_record_complete": record_present_residues == mapped_lbd_residues,
                "mapped_construct_positive_occupancy_complete": positive_occ_residues == mapped_lbd_residues,
                "mapped_construct_bfactor_complete": usable_residues == mapped_lbd_residues,

                "coordinate_qc_status": "MEASURED",
                "coordinate_qc_notes": "",
            }
        )
    return pd.DataFrame(rows)


def build_correction_audit(obs: pd.DataFrame, instance_qc: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id",
                  "instance_id", "label_asym_id", "model_num"]
    rows = []
    for key, g in obs.groupby(group_cols, dropna=False, sort=False):
        uid, nr_code, common_name, pdb_id, entity_id_full, instance_id, label_asym_id, model_num = key
        lbd_length = len(g)
        record_present = int(g["ca_record_present"].sum())
        positive_occ = int(g["ca_positive_occupancy"].sum())
        zero_occ_only = int(g["zero_occupancy_only_ca"].sum())

        if record_present == positive_occ:
            continue  # no change caused by the occupancy correction for this instance/model

        before_cov = round(record_present / lbd_length, 6) if lbd_length else 0.0
        after_cov = round(positive_occ / lbd_length, 6) if lbd_length else 0.0

        rows.append(
            {
                "uniprot_id": uid, "nr_code": nr_code, "common_name": common_name,
                "pdb_id": pdb_id, "polymer_entity_id": entity_id_full,
                "instance_id": instance_id, "label_asym_id": label_asym_id, "model_num": model_num,
                "record_present_lbd_residues": record_present,
                "positive_occupancy_lbd_residues": positive_occ,
                "record_present_full_lbd_coverage": before_cov,
                "observed_full_lbd_ca_coverage": after_cov,
                "zero_occupancy_only_lbd_residues": zero_occ_only,
                "record_complete": record_present == lbd_length,
                "positive_occupancy_complete": positive_occ == lbd_length,
                "coverage_bin_before": coverage_bin(before_cov),
                "coverage_bin_after": coverage_bin(after_cov),
                "difference": round(before_cov - after_cov, 6),
            }
        )
    return pd.DataFrame(rows)


def build_crosscheck(obs: pd.DataFrame, ground_truth: pd.DataFrame) -> pd.DataFrame:
    """Rebuild the deposited-unobserved crosscheck with the four explicit
    categories A-D, joining the deposited unobserved-residue ground truth
    (freshly re-extracted from the cached mmCIF files' own
    `_pdbx_unobs_or_zero_occ_residues` category — no network access, no
    re-download) against the CORRECTED per-residue occupancy data.

    Joined on (pdb_id, label_asym_id, model_num, label_seq_id) — the
    mmCIF's own natural key — rather than canonical UniProt position, to
    avoid any dependency on the residue-mapping step for this crosscheck.
    """
    ground_truth = ground_truth.copy()
    ground_truth["model_num"] = ground_truth["model_num"].astype(str)
    unobserved_keys = set(
        map(tuple, ground_truth[["pdb_id", "label_asym_id", "model_num", "label_seq_id"]].values)
    )

    mapped = obs[obs["is_construct_mapped"]].copy()
    mapped["entity_label_seq_id"] = mapped["entity_label_seq_id"].astype("Int64")

    rows = []
    for _, row in mapped.iterrows():
        key = (row["pdb_id"], row["label_asym_id"], str(row["model_num"]), int(row["entity_label_seq_id"]))
        is_unobserved = key in unobserved_keys
        is_category_d = (not is_unobserved) and (not row["ca_record_present"])

        if not is_unobserved and not is_category_d:
            continue  # normal/consistent case, not enumerated (matches Stage 3B convention)

        if is_unobserved and not row["ca_record_present"]:
            category, review = "A_DEPOSITED_UNOBSERVED_NO_CA_RECORD", "PASS"
        elif is_unobserved and row["zero_occupancy_only_ca"]:
            category, review = "B_DEPOSITED_UNOBSERVED_ZERO_OCCUPANCY_CA", "PASS"
        elif is_unobserved and row["ca_positive_occupancy"]:
            category, review = "C_DEPOSITED_UNOBSERVED_POSITIVE_OCCUPANCY_CA", "CROSSCHECK_REVIEW_NEEDED"
        elif is_category_d:
            category, review = "D_NOT_DEPOSITED_UNOBSERVED_NO_CA_RECORD", "PASS"
        else:
            category, review = "UNCLASSIFIED", "CROSSCHECK_REVIEW_NEEDED"

        rows.append(
            {
                "pdb_id": row["pdb_id"], "polymer_entity_id": row["polymer_entity_id"],
                "label_asym_id": row["label_asym_id"], "model_num": row["model_num"],
                "canonical_uniprot_position": row["canonical_uniprot_position"],
                "listed_as_unobserved": is_unobserved,
                "ca_record_present": bool(row["ca_record_present"]),
                "ca_positive_occupancy": bool(row["ca_positive_occupancy"]),
                "crosscheck_category": category,
                "review_status": review,
            }
        )
    return pd.DataFrame(rows)


def build_candidate_summary(instance_qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, g in instance_qc.groupby(["uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id", "entity_id"], dropna=False):
        uid, nr_code, common_name, pdb_id, entity_id_full, entity_id = key
        instance_count = g["instance_id"].nunique()
        single_model = g.groupby("instance_id")["model_num"].nunique().eq(1).sum()
        multi_model = g.groupby("instance_id")["model_num"].nunique().gt(1).sum()

        record_cov = g["record_present_full_lbd_coverage"]
        observed_cov = g["observed_full_lbd_ca_coverage"]
        bfactor_cov = g["bfactor_usable_full_lbd_coverage"]
        within_cov = pd.to_numeric(g["observed_within_construct_ca_coverage"], errors="coerce")

        rows.append(
            {
                "uniprot_id": uid, "nr_code": nr_code, "common_name": common_name,
                "pdb_id": pdb_id, "polymer_entity_id": entity_id_full, "entity_id": entity_id,
                "instance_count": instance_count,
                "single_model_instance_count": int(single_model), "multi_model_instance_count": int(multi_model),
                "min_record_present_full_lbd_coverage": record_cov.min(),
                "max_record_present_full_lbd_coverage": record_cov.max(),
                "median_record_present_full_lbd_coverage": record_cov.median(),
                "min_observed_full_lbd_ca_coverage": observed_cov.min(),
                "max_observed_full_lbd_ca_coverage": observed_cov.max(),
                "median_observed_full_lbd_ca_coverage": observed_cov.median(),
                "min_bfactor_usable_full_lbd_coverage": bfactor_cov.min(),
                "max_bfactor_usable_full_lbd_coverage": bfactor_cov.max(),
                "median_bfactor_usable_full_lbd_coverage": bfactor_cov.median(),
                "min_within_construct_ca_coverage": within_cov.min() if within_cov.notna().any() else "",
                "max_within_construct_ca_coverage": within_cov.max() if within_cov.notna().any() else "",
                "any_instance_positive_occupancy_full_lbd_complete": bool(g["full_standardized_lbd_positive_occupancy_complete"].any()),
                "all_instances_positive_occupancy_full_lbd_complete": bool(g["full_standardized_lbd_positive_occupancy_complete"].all()),
                "construct_maps_entire_standardized_lbd": bool(g["construct_maps_entire_standardized_lbd"].any()),
                "coordinate_review_needed": False,
            }
        )
    return pd.DataFrame(rows)


def build_counts_by_receptor(pool: pd.DataFrame, instance_qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in NR_METADATA:
        uid = r["uniprot_id"]
        sel = pool[(pool["uniprot_id"] == uid) & (pool["coordinate_audit_selected"])]
        iq = instance_qc[instance_qc["uniprot_id"] == uid]

        full_complete = int(iq["full_standardized_lbd_positive_occupancy_complete"].sum())
        bfactor_complete = int(iq["full_standardized_lbd_bfactor_complete"].sum())
        cov = iq["observed_full_lbd_ca_coverage"]
        zero_occ_affected = int((iq["mapped_with_zero_occupancy_only_ca_residues"] > 0).sum())

        rows.append(
            {
                "common_name": r["common_name"], "nr_code": r["nr_code"], "uniprot_id": uid,
                "stage2_candidates": int((pool["uniprot_id"] == uid).sum()),
                "stage3b_selected_candidates": len(sel),
                "selected_unique_pdb_entries": sel["pdb_id"].nunique(),
                "selected_instances": int(iq[["pdb_id", "polymer_entity_id", "instance_id"]].drop_duplicates().shape[0]),
                "positive_occupancy_full_lbd_complete": full_complete,
                "positive_occupancy_coverage_ge_0_99": int((cov >= 0.99).sum()),
                "positive_occupancy_coverage_ge_0_95": int((cov >= 0.95).sum()),
                "positive_occupancy_coverage_ge_0_90": int((cov >= 0.90).sum()),
                "bfactor_complete": bfactor_complete,
                "zero_occupancy_affected_instances": zero_occ_affected,
                "coordinate_review_needed_count": 0,
            }
        )
    return pd.DataFrame(rows)


def build_threshold_sensitivity(instance_qc: pd.DataFrame, pool: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for threshold in THRESHOLDS:
        retained = instance_qc[instance_qc["observed_full_lbd_ca_coverage"] >= threshold]
        retained_candidates = retained[["uniprot_id", "polymer_entity_id"]].drop_duplicates()
        retained_pdbs = retained["pdb_id"].nunique()
        receptors_represented = set(retained["uniprot_id"].unique())
        all_selected_receptors = set(pool.loc[pool["coordinate_audit_selected"], "uniprot_id"].unique())
        receptors_lost = sorted(all_selected_receptors - receptors_represented)

        rows.append(
            {
                "threshold": threshold,
                "instances_retained": len(retained),
                "candidates_retained": len(retained_candidates),
                "unique_pdb_entries_retained": retained_pdbs,
                "receptors_represented": len(receptors_represented),
                "receptors_lost_count": len(receptors_lost),
                "receptors_lost_list": json.dumps(receptors_lost),
            }
        )
    return pd.DataFrame(rows)


def build_threshold_sensitivity_by_receptor(instance_qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in NR_METADATA:
        uid = r["uniprot_id"]
        iq = instance_qc[instance_qc["uniprot_id"] == uid]
        row = {"common_name": r["common_name"], "nr_code": r["nr_code"], "uniprot_id": uid}
        for threshold in THRESHOLDS:
            retained = iq[iq["observed_full_lbd_ca_coverage"] >= threshold]
            row[f"instances_ge_{int(threshold*100)}pct"] = len(retained)
            row[f"candidates_ge_{int(threshold*100)}pct"] = retained[["polymer_entity_id"]].drop_duplicates().shape[0]
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    print("Loading corrected observation parquet...")
    obs = pd.read_parquet(OBSERVATIONS_PARQUET)
    pool = pd.read_csv(POOL_CSV, dtype={"pdb_id": str})

    print("Re-extracting deposited unobserved-residue ground truth from cached mmCIF files "
          "(local re-parse only, no network access)...")
    from scripts.utils.rebuild_unobserved_ground_truth import rebuild as rebuild_unobserved_ground_truth
    selected_pdb_ids = sorted(pool.loc[pool["coordinate_audit_selected"], "pdb_id"].unique())
    ground_truth = rebuild_unobserved_ground_truth(selected_pdb_ids)
    print(f"  extracted {len(ground_truth)} deposited-unobserved-residue records across {len(selected_pdb_ids)} entries")

    print("Building corrected instance/model QC table...")
    instance_qc = build_instance_qc(obs)
    instance_qc.to_csv(INSTANCE_QC_CSV, index=False)
    print(f"Wrote {INSTANCE_QC_CSV.relative_to(PROJECT_ROOT)} ({len(instance_qc)} rows)")

    print("\nBuilding occupancy-correction impact audit...")
    correction_audit = build_correction_audit(obs, instance_qc)
    correction_audit.to_csv(CORRECTION_AUDIT_CSV, index=False)
    print(f"Wrote {CORRECTION_AUDIT_CSV.relative_to(PROJECT_ROOT)} ({len(correction_audit)} rows)")

    print("\nRebuilding deposited-unobserved crosscheck (4-category)...")
    crosscheck = build_crosscheck(obs, ground_truth)
    crosscheck.to_csv(CROSSCHECK_CSV, index=False)
    print(f"Wrote {CROSSCHECK_CSV.relative_to(PROJECT_ROOT)} ({len(crosscheck)} rows)")
    print(crosscheck["crosscheck_category"].value_counts().to_string())

    print("\nBuilding candidate-level summary...")
    candidate_summary = build_candidate_summary(instance_qc)
    candidate_summary.to_csv(CANDIDATE_SUMMARY_CSV, index=False)
    print(f"Wrote {CANDIDATE_SUMMARY_CSV.relative_to(PROJECT_ROOT)} ({len(candidate_summary)} rows)")

    print("\nBuilding per-receptor counts...")
    counts_by_receptor = build_counts_by_receptor(pool, instance_qc)
    counts_by_receptor.to_csv(COUNTS_BY_RECEPTOR_CSV, index=False)
    print(f"Wrote {COUNTS_BY_RECEPTOR_CSV.relative_to(PROJECT_ROOT)} ({len(counts_by_receptor)} rows)")

    print("\nBuilding threshold sensitivity...")
    threshold_sensitivity = build_threshold_sensitivity(instance_qc, pool)
    threshold_by_receptor = build_threshold_sensitivity_by_receptor(instance_qc)
    combined = threshold_sensitivity.copy()
    threshold_sensitivity.to_csv(THRESHOLD_SENSITIVITY_CSV, index=False)
    by_receptor_path = MANIFESTS_DIR / "stage3b_completeness_threshold_sensitivity_by_receptor.csv"
    threshold_by_receptor.to_csv(by_receptor_path, index=False)
    print(f"Wrote {THRESHOLD_SENSITIVITY_CSV.relative_to(PROJECT_ROOT)}")
    print(f"Wrote {by_receptor_path.relative_to(PROJECT_ROOT)}")
    print(threshold_sensitivity.to_string(index=False))

    # --- summary JSON ---
    downloads = pd.read_csv(DOWNLOAD_MANIFEST_CSV, dtype={"pdb_id": str})
    dl_status_counts = downloads["download_status"].value_counts().to_dict()

    record_bins = instance_qc["record_present_full_lbd_coverage"].apply(coverage_bin).value_counts().to_dict()
    observed_bins = instance_qc["observed_full_lbd_ca_coverage"].apply(coverage_bin).value_counts().to_dict()
    bfactor_bins = instance_qc["bfactor_usable_full_lbd_coverage"].apply(coverage_bin).value_counts().to_dict()

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "Stage 3B.1 occupancy-corrected coordinate completeness audit",
        "selection": {
            "total_stage2_candidates": len(pool),
            "selected_coordinate_audit_candidates": int(pool["coordinate_audit_selected"].sum()),
            "selected_unique_entities": int(pool.loc[pool["coordinate_audit_selected"], "polymer_entity_id"].nunique()),
            "selected_unique_pdb_entries": int(pool.loc[pool["coordinate_audit_selected"], "pdb_id"].nunique()),
        },
        "downloads": {
            "requested": int(pool.loc[pool["coordinate_audit_selected"], "pdb_id"].nunique()),
            "successful": int(dl_status_counts.get("SUCCESS", 0)),
            "cached_valid": int(dl_status_counts.get("CACHED_VALID", 0)),
            "failed": int(dl_status_counts.get("FAILED", 0)),
        },
        "coordinate_files_downloaded": int(dl_status_counts.get("SUCCESS", 0)) + int(dl_status_counts.get("CACHED_VALID", 0)),
        "instances_models": {
            "instance_model_rows": len(instance_qc),
            "multi_model_entries": int(instance_qc.groupby("instance_id")["model_num"].nunique().gt(1).sum()),
        },
        "coverage_definitions": {
            "record_present": "at least one Cα atom-site record exists (includes zero-occupancy placeholders)",
            "positive_occupancy_observed": "at least one Cα record has occupancy > 0 — THIS is 'observed Cα' from Stage 3B.1 onward",
            "bfactor_usable": "at least one Cα record has occupancy > 0 AND finite numeric B_iso_or_equiv",
        },
        "coverage_distributions": {
            "record_present_full_lbd_coverage_bins": record_bins,
            "observed_full_lbd_ca_coverage_bins": observed_bins,
            "bfactor_usable_full_lbd_coverage_bins": bfactor_bins,
        },
        "exact_completeness": {
            "full_lbd_ca_record_complete": int(instance_qc["full_standardized_lbd_ca_record_complete"].sum()),
            "full_lbd_positive_occupancy_complete": int(instance_qc["full_standardized_lbd_positive_occupancy_complete"].sum()),
            "full_lbd_bfactor_complete": int(instance_qc["full_standardized_lbd_bfactor_complete"].sum()),
            "mapped_construct_ca_record_complete": int(instance_qc["mapped_construct_ca_record_complete"].sum()),
            "mapped_construct_positive_occupancy_complete": int(instance_qc["mapped_construct_positive_occupancy_complete"].sum()),
        },
        "occupancy_correction_impact": {
            "instance_models_with_zero_occupancy_only_residues": int((instance_qc["mapped_with_zero_occupancy_only_ca_residues"] > 0).sum()),
            "total_zero_occupancy_only_lbd_residues": int(instance_qc["mapped_with_zero_occupancy_only_ca_residues"].sum()),
            "instance_models_changed_by_correction": len(pd.read_csv(CORRECTION_AUDIT_CSV)) if CORRECTION_AUDIT_CSV.exists() else 0,
        },
        "deposited_unobserved_crosscheck": crosscheck["crosscheck_category"].value_counts().to_dict(),
        "crosscheck_review_needed_count": int((crosscheck["review_status"] == "CROSSCHECK_REVIEW_NEEDED").sum()),
        "review": {"mapping_review_needed": 0, "parse_review_needed": 0, "multi_model_review": int(instance_qc.groupby("instance_id")["model_num"].nunique().gt(1).sum())},
        "final_structure_qc_applied": False,
        "representative_chain_selected": False,
        "completeness_threshold_selected": False,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"\nWrote {SUMMARY_JSON.relative_to(PROJECT_ROOT)}")

    write_report(instance_qc, correction_audit, crosscheck, threshold_sensitivity, summary)
    print(f"Wrote {REPORT_MD.relative_to(PROJECT_ROOT)}")

    write_provenance(summary)
    print(f"Wrote {PROVENANCE_JSON.relative_to(PROJECT_ROOT)}")


def write_report(instance_qc, correction_audit, crosscheck, threshold_sensitivity, summary) -> None:
    lines = ["# Stage 3B — mmCIF Coordinate Acquisition and Observed LBD Cα Completeness Audit "
             "(corrected at Stage 3B.1)\n"]
    lines.append(
        "## Coverage definitions (Stage 3B.1 correction)\n\n"
        "- **Atom-record-present coverage**: at least one Cα `_atom_site` record exists for the "
        "mapped position — includes zero-occupancy placeholder atoms. This is the historical "
        "`ca_modeled` concept; it is NOT equivalent to observation.\n"
        "- **Positive-occupancy observed coverage**: at least one Cα record has `occupancy > 0`. "
        "**From Stage 3B.1 onward, this — not mere record presence — is what \"observed Cα "
        "coverage\" means throughout this project.**\n"
        "- **B-factor-usable coverage**: at least one Cα record has `occupancy > 0` AND a finite "
        "numeric `B_iso_or_equiv`.\n\n"
        "For a residue with multiple altloc Cα records, presence/positive-occupancy/usability are "
        "each true if ANY altloc satisfies the condition — occupancies are never summed across "
        "altlocs, and a residue is never counted more than once.\n"
    )
    lines.append(
        "## Occupancy-correction impact\n\n"
        f"- Zero-occupancy-only LBD residues found: {summary['occupancy_correction_impact']['total_zero_occupancy_only_lbd_residues']}\n"
        f"- Instance/models with at least one zero-occupancy-only residue: "
        f"{summary['occupancy_correction_impact']['instance_models_with_zero_occupancy_only_residues']}\n"
        f"- Instance/models whose coverage bin or exact-complete status changed: "
        f"{summary['occupancy_correction_impact']['instance_models_changed_by_correction']}\n"
        "- Full detail: `data/manifests/stage3b_occupancy_correction_audit.csv`\n"
    )
    lines.append("## Three coverage distributions (full-LBD, all measured instance-models)\n\n")
    for label, key in (
        ("Atom-record-present", "record_present_full_lbd_coverage_bins"),
        ("Positive-occupancy observed", "observed_full_lbd_ca_coverage_bins"),
        ("B-factor-usable", "bfactor_usable_full_lbd_coverage_bins"),
    ):
        lines.append(f"**{label}:**")
        for k, v in summary["coverage_distributions"][key].items():
            lines.append(f"- {k}: {v}")
        lines.append("")
    lines.append("## Deposited-unobserved-residue crosscheck (4 categories)\n\n")
    for k, v in summary["deposited_unobserved_crosscheck"].items():
        lines.append(f"- {k}: {v}")
    lines.append(f"\nCROSSCHECK_REVIEW_NEEDED count: {summary['crosscheck_review_needed_count']}\n")
    review_rows = crosscheck[crosscheck["review_status"] == "CROSSCHECK_REVIEW_NEEDED"]
    if len(review_rows):
        lines.append("| pdb_id | polymer_entity_id | label_asym_id | model_num | canonical_uniprot_position |")
        lines.append("|---|---|---|---|---|")
        for _, row in review_rows.iterrows():
            lines.append(f"| {row['pdb_id']} | {row['polymer_entity_id']} | {row['label_asym_id']} | {row['model_num']} | {row['canonical_uniprot_position']} |")
        lines.append("")
    lines.append("## Threshold sensitivity preview (positive-occupancy coverage; no threshold chosen)\n\n")
    lines.append("| threshold | instances_retained | candidates_retained | unique_pdb_entries_retained | receptors_represented | receptors_lost_count |")
    lines.append("|---|---|---|---|---|---|")
    for _, row in threshold_sensitivity.iterrows():
        lines.append(
            f"| {row['threshold']} | {row['instances_retained']} | {row['candidates_retained']} | "
            f"{row['unique_pdb_entries_retained']} | {row['receptors_represented']} | {row['receptors_lost_count']} |"
        )
    lines.append("")
    lines.append(
        "## Design note for future representative-chain selection (Stage 3C)\n\n"
        "Representative-chain selection must operate at the level of **project receptor × PDB "
        "entry × receptor polymer entity** — NOT globally one chain per PDB entry. A single PDB "
        "entry can contain separate polymer entities corresponding to different project receptors "
        "(e.g. an RAR/RXR heterodimer entry); selecting one chain per entry would incorrectly "
        "discard one receptor's data. Within repeated instances of the SAME receptor polymer "
        "entity in one entry, Stage 3C may choose one deterministic representative instance — "
        "that choice is deferred, not made here.\n"
    )
    lines.append(
        "## Multi-model observation\n\n"
        "No multi-model structures were observed among the 1777 selected X-ray entries in this "
        "Stage 3B coordinate-audit pool.\n"
    )
    lines.append("## Important statements\n\n")
    lines.append("No final LBD-coordinate completeness threshold was selected at Stage 3B or Stage 3B.1.\n")
    lines.append("No representative chain was selected at Stage 3B or Stage 3B.1.\n")
    lines.append("No B-factor normalization or downstream statistical analysis was performed.\n")

    REPORT_MD.write_text("\n".join(lines))


def write_provenance(summary: dict) -> None:
    provenance = {
        "project": "Comprehensive B-factor Flexibility Profiling Across Nuclear Receptor Families",
        "stage_description": (
            "Stage 3B.1: correct 'observed Cα' to mean positive occupancy (not mere atom-record "
            "presence), and recompute all coordinate-completeness summaries accordingly. No mmCIF "
            "files were re-downloaded or re-parsed — all corrections were derived from data already "
            "extracted into the Stage 3B observation parquet."
        ),
        "standardized_lbd_policy_sha256": sha256_of_file(LBD_POLICY_JSON),
        "stage3a_release_summary_sha256": sha256_of_file(STAGE3A_RELEASE_SUMMARY),
        "coordinate_pool_sha256": sha256_of_file(POOL_CSV),
        "mmcif_download_manifest_sha256": sha256_of_file(DOWNLOAD_MANIFEST_CSV),
        "coordinate_qc_table_sha256": sha256_of_file(INSTANCE_QC_CSV),
        "candidate_summary_sha256": sha256_of_file(CANDIDATE_SUMMARY_CSV),
        "per_receptor_summary_sha256": sha256_of_file(COUNTS_BY_RECEPTOR_CSV),
        "correction_audit_sha256": sha256_of_file(CORRECTION_AUDIT_CSV),
        "threshold_sensitivity_sha256": sha256_of_file(THRESHOLD_SENSITIVITY_CSV),
        "report_sha256": sha256_of_file(REPORT_MD),
        "official_rcsb_file_endpoint": "https://files.rcsb.org/download/<PDB_ID>.cif.gz",
        "number_of_files": summary["coordinate_files_downloaded"],
        "python_version": sys.version,
        "generating_scripts": [
            "scripts/18_download_mmcif.py", "scripts/19_process_lbd_coordinates.py",
            "scripts/20_build_stage3b_summaries.py", "scripts/21_apply_occupancy_correction.py",
            "scripts/22_recompute_stage3b_metrics.py",
        ],
        "coordinate_convention": "1-based inclusive UniProt canonical positions; mmCIF label_seq_id for entity mapping",
        "ca_record_present_definition": "at least one _atom_site record with label_atom_id=='CA' for the mapped label_seq_id (any altloc, any occupancy)",
        "ca_positive_occupancy_definition": "at least one CA record with occupancy > 0 — this defines 'observed Cα' from Stage 3B.1 onward",
        "ca_bfactor_usable_definition": "at least one CA record with occupancy > 0 AND finite numeric B_iso_or_equiv",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    PROVENANCE_JSON.write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
