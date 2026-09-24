"""
Stage 3B, step 2: parse every successfully downloaded selected mmCIF file
and build the per-residue observed-Cα observation table plus all
instance/candidate/receptor-level coordinate-completeness summaries.

Reads:
    data/manifests/stage3b_coordinate_pool.csv
    data/manifests/mmcif_download_manifest.csv
    data/manifests/nr_lbd_reference.csv
    data/manifests/rcsb_sequence_coordinate_alignments.csv
    data/manifests/rcsb_polymer_instances.csv
    data/manifests/rcsb_multi_uniprot_entity_audit.csv
    data/manifests/rcsb_source_taxonomy_audit.csv
    data/manifests/rcsb_precoordinate_qc.csv
    data/raw/mmcif/{PDB_ID}.cif.gz

Writes:
    data/interim/qc/stage3b_lbd_ca_observations.parquet   (gitignored; large)
    data/manifests/rcsb_lbd_coordinate_instance_qc.csv
    data/manifests/rcsb_lbd_coordinate_candidate_summary.csv
    data/manifests/mmcif_unobserved_residue_crosscheck.csv

Canonical mapping uses ONLY label_seq_id (never auth_seq_id) — see
scripts/utils/residue_mapping.py. No representative chain is chosen, no
completeness threshold is applied, no B-factor normalization occurs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.mmcif_parse import AtomSiteIndex, get_entry_metadata, get_unobserved_residues, parse_mmcif_gz  # noqa: E402
from scripts.utils.residue_mapping import build_entity_to_uniprot_map, invert_mapping  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
INTERIM_QC_DIR = PROJECT_ROOT / "data" / "interim" / "qc"
RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"

POOL_CSV = MANIFESTS_DIR / "stage3b_coordinate_pool.csv"
DOWNLOAD_MANIFEST_CSV = MANIFESTS_DIR / "mmcif_download_manifest.csv"
LBD_REFERENCE_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"
ALIGNMENTS_CSV = MANIFESTS_DIR / "rcsb_sequence_coordinate_alignments.csv"
INSTANCES_CSV = MANIFESTS_DIR / "rcsb_polymer_instances.csv"
MULTI_UNIPROT_CSV = MANIFESTS_DIR / "rcsb_multi_uniprot_entity_audit.csv"
TAXONOMY_AUDIT_CSV = MANIFESTS_DIR / "rcsb_source_taxonomy_audit.csv"
PRECOORDINATE_QC_CSV = MANIFESTS_DIR / "rcsb_precoordinate_qc.csv"

OBSERVATIONS_PARQUET = INTERIM_QC_DIR / "stage3b_lbd_ca_observations.parquet"
INSTANCE_QC_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_instance_qc.csv"
CANDIDATE_SUMMARY_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_candidate_summary.csv"
UNOBSERVED_CROSSCHECK_CSV = MANIFESTS_DIR / "mmcif_unobserved_residue_crosscheck.csv"

BATCH_SIZE_PDBS = 100

OBS_SCHEMA = pa.schema(
    [
        ("uniprot_id", pa.string()), ("nr_code", pa.string()), ("common_name", pa.string()),
        ("pdb_id", pa.string()), ("polymer_entity_id", pa.string()), ("entity_id", pa.string()),
        ("instance_id", pa.string()), ("label_asym_id", pa.string()), ("auth_asym_id", pa.string()),
        ("model_num", pa.string()),
        ("canonical_uniprot_position", pa.int64()), ("lbd_relative_position", pa.int64()),
        ("entity_label_seq_id", pa.float64()),  # nullable int
        ("is_construct_mapped", pa.bool_()),
        ("ca_modeled", pa.bool_()), ("ca_bfactor_usable", pa.bool_()),
        ("ca_atom_record_count", pa.int64()),
        ("label_comp_ids_json", pa.string()), ("auth_comp_ids_json", pa.string()),
        ("group_pdb_values_json", pa.string()),
        ("auth_seq_ids_json", pa.string()), ("insertion_codes_json", pa.string()),
        ("altloc_ids_json", pa.string()), ("occupancies_json", pa.string()), ("b_iso_values_json", pa.string()),
        ("mapping_status", pa.string()), ("observation_status", pa.string()),
    ]
)


def load_reference_data():
    pool = pd.read_csv(POOL_CSV, dtype={"pdb_id": str, "entity_id": str})
    selected = pool[pool["coordinate_audit_selected"]].copy()

    downloads = pd.read_csv(DOWNLOAD_MANIFEST_CSV, dtype={"pdb_id": str})
    lbd_ref = pd.read_csv(LBD_REFERENCE_CSV).set_index("uniprot_id")
    alignments = pd.read_csv(ALIGNMENTS_CSV, dtype={"pdb_id": str, "entity_id": str})
    alignments = alignments[alignments["is_stage2_candidate"]].set_index(["uniprot_id", "polymer_entity_id"])
    instances = pd.read_csv(INSTANCES_CSV, dtype={"pdb_id": str})
    instances = instances[instances["instance_metadata_status"] == "PASS"]
    multi_uniprot = set(pd.read_csv(MULTI_UNIPROT_CSV)["polymer_entity_id"]) if MULTI_UNIPROT_CSV.exists() else set()
    taxonomy = pd.read_csv(TAXONOMY_AUDIT_CSV, dtype={"pdb_id": str}).set_index("polymer_entity_id")["taxonomy_audit_status"].to_dict()

    return selected, downloads, lbd_ref, alignments, instances, multi_uniprot, taxonomy


def json_field(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def process_one_pdb(
    pdb_id: str,
    candidates: pd.DataFrame,
    lbd_ref: pd.DataFrame,
    alignments: pd.DataFrame,
    instances: pd.DataFrame,
) -> tuple[list[dict], list[dict], list[dict]]:
    """Returns (observation_rows, instance_qc_rows, unobserved_crosscheck_rows) for all candidates in this PDB."""
    mmcif_path = RAW_MMCIF_DIR / f"{pdb_id.upper()}.cif.gz"
    obs_rows: list[dict] = []
    instance_qc_rows: list[dict] = []
    crosscheck_rows: list[dict] = []

    try:
        d = parse_mmcif_gz(mmcif_path)
    except Exception as exc:  # noqa: BLE001 - any parse failure must be captured, not crash the batch
        for _, cand in candidates.iterrows():
            instance_qc_rows.append(
                _blank_instance_qc_row(cand, lbd_ref, coordinate_qc_status="PARSE_REVIEW_NEEDED", notes=f"mmCIF parse failed: {exc}")
            )
        return obs_rows, instance_qc_rows, crosscheck_rows

    entry_meta = get_entry_metadata(d)
    idx = AtomSiteIndex(d)
    unobs = get_unobserved_residues(d)
    unobs_set = {(r["label_asym_id"], r["model_num"], r["label_seq_id"]) for r in unobs}

    for _, cand in candidates.iterrows():
        uid, entity_id_full = cand["uniprot_id"], cand["polymer_entity_id"]

        if uid not in lbd_ref.index:
            instance_qc_rows.append(_blank_instance_qc_row(cand, lbd_ref, "MAPPING_REVIEW_NEEDED", "No LBD reference row."))
            continue
        ref_row = lbd_ref.loc[uid]
        lbd_start, lbd_end = int(ref_row["lbd_start"]), int(ref_row["lbd_end"])
        lbd_length = lbd_end - lbd_start + 1

        align_key = (uid, entity_id_full)
        if align_key not in alignments.index:
            instance_qc_rows.append(_blank_instance_qc_row(cand, lbd_ref, "MAPPING_REVIEW_NEEDED", "No Sequence Coordinates alignment found."))
            continue
        align_row = alignments.loc[align_key]
        regions = json.loads(align_row["aligned_regions_json"])
        map_result = build_entity_to_uniprot_map(regions)
        if map_result["status"] != "OK":
            instance_qc_rows.append(_blank_instance_qc_row(cand, lbd_ref, "MAPPING_REVIEW_NEEDED", map_result["notes"]))
            continue
        entity_to_uniprot = map_result["mapping"]
        inv, ambiguous = invert_mapping(entity_to_uniprot)
        if ambiguous:
            instance_qc_rows.append(
                _blank_instance_qc_row(cand, lbd_ref, "MAPPING_REVIEW_NEEDED", f"Ambiguous UniProt positions: {sorted(ambiguous)}")
            )
            continue

        entity_instances = instances[instances["polymer_entity_id"] == entity_id_full]
        if len(entity_instances) == 0:
            instance_qc_rows.append(_blank_instance_qc_row(cand, lbd_ref, "MAPPING_REVIEW_NEEDED", "No polymer instances found for this entity."))
            continue

        for _, inst in entity_instances.iterrows():
            label_asym_id, auth_asym_id, instance_id = inst["label_asym_id"], inst["auth_asym_id"], inst["instance_id"]
            models_for_asym = sorted({m for (a, m) in idx.index.keys() if a == label_asym_id})

            if not models_for_asym:
                instance_qc_rows.append(
                    _blank_instance_qc_row(cand, lbd_ref, "PARSE_REVIEW_NEEDED",
                                            f"No CA atoms found for label_asym_id={label_asym_id!r} in this entry.",
                                            instance_id=instance_id, label_asym_id=label_asym_id, auth_asym_id=auth_asym_id)
                )
                continue

            qc_status = "MULTI_MODEL_REVIEW" if len(models_for_asym) > 1 else "MEASURED"

            for model_num in models_for_asym:
                anisotrop_ca_count = 0
                for rel_pos, canonical_pos in enumerate(range(lbd_start, lbd_end + 1), start=1):
                    label_seq_id = inv.get(canonical_pos)
                    is_mapped = label_seq_id is not None
                    ca_rows = idx.get_residue_ca_rows(label_asym_id, model_num, label_seq_id) if is_mapped else []
                    ca_modeled = len(ca_rows) > 0
                    ca_usable = any(
                        r["B_iso_or_equiv"] is not None and r["occupancy"] is not None and r["occupancy"] > 0
                        for r in ca_rows
                    )
                    if entry_meta["anisotrop_category_present"] and any(
                        r["atom_id"] in entry_meta["anisotrop_atom_ids"] for r in ca_rows
                    ):
                        anisotrop_ca_count += 1

                    obs_rows.append(
                        {
                            "uniprot_id": uid, "nr_code": cand["nr_code"], "common_name": cand["common_name"],
                            "pdb_id": pdb_id, "polymer_entity_id": entity_id_full, "entity_id": cand["entity_id"],
                            "instance_id": instance_id, "label_asym_id": label_asym_id, "auth_asym_id": auth_asym_id,
                            "model_num": model_num,
                            "canonical_uniprot_position": canonical_pos, "lbd_relative_position": rel_pos,
                            "entity_label_seq_id": float(label_seq_id) if is_mapped else None,
                            "is_construct_mapped": is_mapped,
                            "ca_modeled": ca_modeled, "ca_bfactor_usable": ca_usable,
                            "ca_atom_record_count": len(ca_rows),
                            "label_comp_ids_json": json_field([r["label_comp_id"] for r in ca_rows]),
                            "auth_comp_ids_json": json_field([r["auth_comp_id"] for r in ca_rows]),
                            "group_pdb_values_json": json_field([r["group_PDB"] for r in ca_rows]),
                            "auth_seq_ids_json": json_field([r["auth_seq_id"] for r in ca_rows]),
                            "insertion_codes_json": json_field([r["pdbx_PDB_ins_code"] for r in ca_rows]),
                            "altloc_ids_json": json_field([r["label_alt_id"] for r in ca_rows]),
                            "occupancies_json": json_field([r["occupancy"] for r in ca_rows]),
                            "b_iso_values_json": json_field([r["B_iso_or_equiv"] for r in ca_rows]),
                            "mapping_status": "OK",
                            "observation_status": "OK",
                        }
                    )

                    if is_mapped and (label_asym_id, model_num, label_seq_id) in unobs_set:
                        crosscheck_rows.append(
                            {
                                "pdb_id": pdb_id, "polymer_entity_id": entity_id_full, "label_asym_id": label_asym_id,
                                "model_num": model_num, "label_seq_id": label_seq_id,
                                "canonical_uniprot_position": canonical_pos,
                                "listed_as_unobserved": True, "has_ca_atom_record": ca_modeled,
                                "consistent": not ca_modeled,
                            }
                        )
                    elif is_mapped and ca_modeled:
                        pass  # not in unobserved table and modeled: expected, no crosscheck row needed
                    elif is_mapped and not ca_modeled:
                        # mapped, not modeled, but ALSO absent from the deposited unobserved-residue table.
                        crosscheck_rows.append(
                            {
                                "pdb_id": pdb_id, "polymer_entity_id": entity_id_full, "label_asym_id": label_asym_id,
                                "model_num": model_num, "label_seq_id": label_seq_id,
                                "canonical_uniprot_position": canonical_pos,
                                "listed_as_unobserved": False, "has_ca_atom_record": False,
                                "consistent": False,
                            }
                        )

                mapped_mask = [r["is_construct_mapped"] for r in obs_rows[-lbd_length:]]
                ca_modeled_list = [r["ca_modeled"] for r in obs_rows[-lbd_length:]]
                ca_usable_list = [r["ca_bfactor_usable"] for r in obs_rows[-lbd_length:]]
                canonical_positions = [r["canonical_uniprot_position"] for r in obs_rows[-lbd_length:]]

                qc_row = _compute_instance_qc(
                    cand, ref_row, lbd_start, lbd_end, lbd_length,
                    instance_id, label_asym_id, auth_asym_id, model_num,
                    mapped_mask, ca_modeled_list, ca_usable_list, canonical_positions,
                    obs_rows[-lbd_length:], qc_status, pdb_id, entry_meta, anisotrop_ca_count,
                )
                instance_qc_rows.append(qc_row)

    return obs_rows, instance_qc_rows, crosscheck_rows


def _blank_instance_qc_row(cand, lbd_ref, coordinate_qc_status, notes, instance_id="", label_asym_id="", auth_asym_id=""):
    uid = cand["uniprot_id"]
    lbd_start = lbd_end = lbd_length = ""
    if uid in lbd_ref.index:
        ref_row = lbd_ref.loc[uid]
        lbd_start, lbd_end, lbd_length = ref_row["lbd_start"], ref_row["lbd_end"], ref_row["lbd_length"]
    return {
        "uniprot_id": uid, "nr_code": cand["nr_code"], "common_name": cand["common_name"], "group": cand["group"],
        "pdb_id": cand["pdb_id"], "polymer_entity_id": cand["polymer_entity_id"], "entity_id": cand["entity_id"],
        "instance_id": instance_id, "label_asym_id": label_asym_id, "auth_asym_id": auth_asym_id, "model_num": "",
        "lbd_start": lbd_start, "lbd_end": lbd_end, "lbd_length": lbd_length,
        "mapped_lbd_residues": "", "construct_lbd_mapping_coverage": "",
        "unmapped_lbd_residues": "",
        "modeled_lbd_ca_residues": "", "observed_full_lbd_ca_coverage": "", "observed_within_construct_ca_coverage": "",
        "mapped_but_unmodeled_lbd_residues": "",
        "bfactor_usable_lbd_ca_residues": "", "bfactor_usable_full_lbd_coverage": "", "bfactor_usable_within_construct_coverage": "",
        "modeled_but_bfactor_unusable_residues": "",
        "n_terminal_unmapped_count": "", "c_terminal_unmapped_count": "",
        "n_terminal_mapped_but_unmodeled_count": "", "c_terminal_mapped_but_unmodeled_count": "",
        "internal_unmodeled_residue_count": "", "internal_missing_run_count": "", "longest_internal_missing_run": "",
        "alternate_ca_residue_count": "", "zero_occupancy_ca_residue_count": "",
        "anisotrop_category_present": "", "lbd_ca_with_anisotrop_count": "",
        "mmcif_sha256": "",
        "coordinate_qc_status": coordinate_qc_status, "coordinate_qc_notes": notes,
    }


def _compute_instance_qc(
    cand, ref_row, lbd_start, lbd_end, lbd_length,
    instance_id, label_asym_id, auth_asym_id, model_num,
    mapped_mask, ca_modeled_list, ca_usable_list, canonical_positions,
    residue_rows, qc_status, pdb_id, entry_meta, anisotrop_ca_count,
):
    mapped_lbd_residues = sum(mapped_mask)
    unmapped_lbd_residues = lbd_length - mapped_lbd_residues
    modeled_lbd_ca_residues = sum(ca_modeled_list)
    bfactor_usable_lbd_ca_residues = sum(ca_usable_list)
    mapped_but_unmodeled = sum(1 for m, c in zip(mapped_mask, ca_modeled_list) if m and not c)
    modeled_but_unusable = sum(1 for c, u in zip(ca_modeled_list, ca_usable_list) if c and not u)

    construct_coverage = round(mapped_lbd_residues / lbd_length, 6) if lbd_length else 0.0
    full_ca_coverage = round(modeled_lbd_ca_residues / lbd_length, 6) if lbd_length else 0.0
    within_construct_ca_coverage = round(modeled_lbd_ca_residues / mapped_lbd_residues, 6) if mapped_lbd_residues else ""
    bfactor_full_coverage = round(bfactor_usable_lbd_ca_residues / lbd_length, 6) if lbd_length else 0.0
    bfactor_within_construct = round(bfactor_usable_lbd_ca_residues / mapped_lbd_residues, 6) if mapped_lbd_residues else ""

    # Terminal vs internal missingness topology.
    n = len(mapped_mask)
    n_term_unmapped = 0
    for m in mapped_mask:
        if not m:
            n_term_unmapped += 1
        else:
            break
    c_term_unmapped = 0
    for m in reversed(mapped_mask):
        if not m:
            c_term_unmapped += 1
        else:
            break

    # Mapped-but-unmodeled at termini: scan from the first mapped residue onward.
    mapped_unmodeled_flags = [m and not c for m, c in zip(mapped_mask, ca_modeled_list)]
    first_mapped_idx = next((i for i, m in enumerate(mapped_mask) if m), None)
    last_mapped_idx = next((i for i in range(n - 1, -1, -1) if mapped_mask[i]), None)

    n_term_mapped_unmodeled = 0
    if first_mapped_idx is not None:
        i = first_mapped_idx
        while i < n and mapped_unmodeled_flags[i]:
            n_term_mapped_unmodeled += 1
            i += 1
    c_term_mapped_unmodeled = 0
    if last_mapped_idx is not None:
        i = last_mapped_idx
        while i >= 0 and mapped_unmodeled_flags[i]:
            c_term_mapped_unmodeled += 1
            i -= 1

    # Internal unmodeled = mapped-but-unmodeled positions strictly between
    # the first and last "modeled" position (i.e. excluding the terminal
    # runs already counted above).
    internal_unmodeled_positions = []
    if first_mapped_idx is not None and last_mapped_idx is not None:
        interior_start = first_mapped_idx + n_term_mapped_unmodeled
        interior_end = last_mapped_idx - c_term_mapped_unmodeled
        for i in range(interior_start, interior_end + 1):
            if mapped_unmodeled_flags[i]:
                internal_unmodeled_positions.append(i)

    internal_unmodeled_count = len(internal_unmodeled_positions)
    # Count contiguous runs.
    runs = []
    if internal_unmodeled_positions:
        run_start = internal_unmodeled_positions[0]
        prev = run_start
        for pos in internal_unmodeled_positions[1:]:
            if pos == prev + 1:
                prev = pos
            else:
                runs.append(prev - run_start + 1)
                run_start = pos
                prev = pos
        runs.append(prev - run_start + 1)
    internal_missing_run_count = len(runs)
    longest_internal_missing_run = max(runs) if runs else 0

    # Altloc / occupancy diagnostics from residue rows.
    alternate_ca_count = 0
    zero_occ_count = 0
    for row in residue_rows:
        altlocs = json.loads(row["altloc_ids_json"])
        occs = json.loads(row["occupancies_json"])
        if len(altlocs) > 1 and any(a not in (".", "?", "") for a in altlocs):
            alternate_ca_count += 1
        if row["ca_modeled"] and occs and all((o is not None and o == 0) for o in occs):
            zero_occ_count += 1

    # Anisotropic cross-reference: computed by the caller (which still has
    # atom-level `atom_id` access) via direct atom_id membership in
    # `_atom_site_anisotrop.id` — never fabricated when the category is
    # absent (left as an explicit 0 only when the category truly is absent,
    # not as a stand-in for "not computed").
    anisotrop_present = entry_meta["anisotrop_category_present"]
    lbd_ca_with_anisotrop = anisotrop_ca_count if anisotrop_present else ""

    mmcif_path = RAW_MMCIF_DIR / f"{pdb_id.upper()}.cif.gz"
    mmcif_sha256 = hashlib.sha256(mmcif_path.read_bytes()).hexdigest() if mmcif_path.exists() else ""

    return {
        "uniprot_id": cand["uniprot_id"], "nr_code": cand["nr_code"], "common_name": cand["common_name"], "group": cand["group"],
        "pdb_id": pdb_id, "polymer_entity_id": cand["polymer_entity_id"], "entity_id": cand["entity_id"],
        "instance_id": instance_id, "label_asym_id": label_asym_id, "auth_asym_id": auth_asym_id, "model_num": model_num,
        "lbd_start": lbd_start, "lbd_end": lbd_end, "lbd_length": lbd_length,
        "mapped_lbd_residues": mapped_lbd_residues, "construct_lbd_mapping_coverage": construct_coverage,
        "unmapped_lbd_residues": unmapped_lbd_residues,
        "modeled_lbd_ca_residues": modeled_lbd_ca_residues,
        "observed_full_lbd_ca_coverage": full_ca_coverage, "observed_within_construct_ca_coverage": within_construct_ca_coverage,
        "mapped_but_unmodeled_lbd_residues": mapped_but_unmodeled,
        "bfactor_usable_lbd_ca_residues": bfactor_usable_lbd_ca_residues,
        "bfactor_usable_full_lbd_coverage": bfactor_full_coverage, "bfactor_usable_within_construct_coverage": bfactor_within_construct,
        "modeled_but_bfactor_unusable_residues": modeled_but_unusable,
        "n_terminal_unmapped_count": n_term_unmapped, "c_terminal_unmapped_count": c_term_unmapped,
        "n_terminal_mapped_but_unmodeled_count": n_term_mapped_unmodeled, "c_terminal_mapped_but_unmodeled_count": c_term_mapped_unmodeled,
        "internal_unmodeled_residue_count": internal_unmodeled_count,
        "internal_missing_run_count": internal_missing_run_count, "longest_internal_missing_run": longest_internal_missing_run,
        "alternate_ca_residue_count": alternate_ca_count, "zero_occupancy_ca_residue_count": zero_occ_count,
        "anisotrop_category_present": anisotrop_present, "lbd_ca_with_anisotrop_count": lbd_ca_with_anisotrop,
        "mmcif_sha256": mmcif_sha256,
        "coordinate_qc_status": qc_status, "coordinate_qc_notes": "",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", type=str, default=None, help="Comma-separated PDB IDs to restrict to (pilot mode).")
    args = parser.parse_args()

    selected, downloads, lbd_ref, alignments, instances, multi_uniprot, taxonomy = load_reference_data()
    usable_pdb_ids = set(downloads.loc[downloads["download_status"].isin(["SUCCESS", "CACHED_VALID"]), "pdb_id"])
    failed_pdb_ids = set(downloads.loc[downloads["download_status"] == "FAILED", "pdb_id"])

    if args.only:
        only_set = set(args.only.split(","))
        selected = selected[selected["pdb_id"].isin(only_set)]

    grouped = selected.groupby("pdb_id")
    pdb_ids_to_process = sorted(grouped.groups.keys())
    print(f"Processing {len(pdb_ids_to_process)} unique selected PDB entries "
          f"({len([p for p in pdb_ids_to_process if p in usable_pdb_ids])} with usable downloads, "
          f"{len([p for p in pdb_ids_to_process if p in failed_pdb_ids])} with failed downloads).")

    INTERIM_QC_DIR.mkdir(parents=True, exist_ok=True)
    writer = None
    all_instance_qc_rows: list[dict] = []
    all_crosscheck_rows: list[dict] = []

    for batch_start in range(0, len(pdb_ids_to_process), BATCH_SIZE_PDBS):
        batch_pdb_ids = pdb_ids_to_process[batch_start: batch_start + BATCH_SIZE_PDBS]
        batch_obs_rows: list[dict] = []

        for pdb_id in batch_pdb_ids:
            candidates = grouped.get_group(pdb_id)
            if pdb_id in failed_pdb_ids:
                for _, cand in candidates.iterrows():
                    all_instance_qc_rows.append(
                        _blank_instance_qc_row(cand, lbd_ref, "PARSE_REVIEW_NEEDED", "mmCIF download failed for this entry.")
                    )
                continue
            if pdb_id not in usable_pdb_ids:
                for _, cand in candidates.iterrows():
                    all_instance_qc_rows.append(
                        _blank_instance_qc_row(cand, lbd_ref, "PARSE_REVIEW_NEEDED", "mmCIF not downloaded in this run.")
                    )
                continue

            obs_rows, instance_qc_rows, crosscheck_rows = process_one_pdb(pdb_id, candidates, lbd_ref, alignments, instances)
            batch_obs_rows.extend(obs_rows)
            all_instance_qc_rows.extend(instance_qc_rows)
            all_crosscheck_rows.extend(crosscheck_rows)

        if batch_obs_rows:
            table = pa.Table.from_pylist(batch_obs_rows, schema=OBS_SCHEMA)
            if writer is None:
                writer = pq.ParquetWriter(OBSERVATIONS_PARQUET, OBS_SCHEMA)
            writer.write_table(table)

        print(f"  batch {batch_start // BATCH_SIZE_PDBS + 1}: processed {len(batch_pdb_ids)} PDBs, "
              f"{len(batch_obs_rows)} observation rows "
              f"({batch_start + len(batch_pdb_ids)}/{len(pdb_ids_to_process)} total)")

    if writer is not None:
        writer.close()
        print(f"\nWrote {OBSERVATIONS_PARQUET.relative_to(PROJECT_ROOT)}")
    else:
        print("\nNo observation rows generated (empty pilot set or all failed).")

    instance_qc_df = pd.DataFrame(all_instance_qc_rows)
    if args.only and INSTANCE_QC_CSV.exists():
        existing = pd.read_csv(INSTANCE_QC_CSV, dtype={"pdb_id": str})
        existing = existing[~existing["pdb_id"].isin(set(args.only.split(",")))]
        instance_qc_df = pd.concat([existing, instance_qc_df], ignore_index=True)
    instance_qc_df.to_csv(INSTANCE_QC_CSV, index=False)
    print(f"Wrote {INSTANCE_QC_CSV.relative_to(PROJECT_ROOT)} ({len(instance_qc_df)} rows)")

    crosscheck_df = pd.DataFrame(all_crosscheck_rows)
    if args.only and UNOBSERVED_CROSSCHECK_CSV.exists():
        existing = pd.read_csv(UNOBSERVED_CROSSCHECK_CSV, dtype={"pdb_id": str})
        existing = existing[~existing["pdb_id"].isin(set(args.only.split(",")))]
        crosscheck_df = pd.concat([existing, crosscheck_df], ignore_index=True)
    crosscheck_df.to_csv(UNOBSERVED_CROSSCHECK_CSV, index=False)
    print(f"Wrote {UNOBSERVED_CROSSCHECK_CSV.relative_to(PROJECT_ROOT)} ({len(crosscheck_df)} rows)")

    if len(crosscheck_df):
        print(crosscheck_df.groupby(["listed_as_unobserved", "has_ca_atom_record"]).size().to_string())


if __name__ == "__main__":
    main()
