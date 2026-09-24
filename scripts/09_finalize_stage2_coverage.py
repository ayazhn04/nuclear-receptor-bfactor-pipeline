"""
Stage 2.2, step 2: robustify the Stage 2.1 coverage computation (interval
UNION instead of naive length sum), compare against the official RCSB
Sequence Coordinates values, and finalize the candidate inventory with both
official and computed coverage fields side by side.

Reads:
    data/manifests/rcsb_polymer_entity_uniprot_mappings.csv   (Stage 2.1)
    data/manifests/rcsb_sequence_coordinate_alignments.csv    (Stage 2.2 step 1)
    data/manifests/rcsb_polymer_entities.csv
    data/manifests/rcsb_candidate_inventory.csv
    data/manifests/rcsb_source_taxonomy_audit.csv
    data/manifests/nr_metadata_uniprot_audit.csv               (Stage 1: reference lengths)

Writes (updated in place):
    data/manifests/rcsb_polymer_entity_uniprot_mappings.csv   (adds union-based computed_* + overlap flags)
    data/manifests/rcsb_candidate_inventory.csv                (adds rcsb_* official + computed_* fields)
    data/manifests/rcsb_source_taxonomy_audit.csv               (adds qc_suitability = NOT_EVALUATED)
    data/manifests/rcsb_discovery_summary.json
    data/manifests/rcsb_stage2_provenance.json
    reports/tables/rcsb_discovery_report.md

Writes (new):
    data/manifests/rcsb_coverage_validation.csv
    data/manifests/rcsb_short_mapping_audit.csv

Tolerance (documented, not arbitrary): both RCSB's official coverage and
our computed coverage are exact rational numbers (aligned residues / total
residues), each independently rounded to a fixed number of decimal places.
The only expected source of difference between two CORRECT computations is
that rounding, which is bounded well under 1e-3 for any sequence longer
than a few residues. We therefore use:
    MATCH                  : abs difference <= 1e-4
    MINOR_NUMERIC_DIFFERENCE: 1e-4 < abs difference <= 1e-2  (rounding-scale)
    REVIEW_NEEDED           : abs difference > 1e-2, or either value missing
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

from scripts.utils.rcsb_extract import compute_union_coverage  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
REPORTS_TABLES_DIR = PROJECT_ROOT / "reports" / "tables"

MAPPINGS_CSV = MANIFESTS_DIR / "rcsb_polymer_entity_uniprot_mappings.csv"
ALIGNMENTS_CSV = MANIFESTS_DIR / "rcsb_sequence_coordinate_alignments.csv"
POLYMER_ENTITIES_CSV = MANIFESTS_DIR / "rcsb_polymer_entities.csv"
CANDIDATE_INVENTORY_CSV = MANIFESTS_DIR / "rcsb_candidate_inventory.csv"
TAXONOMY_AUDIT_CSV = MANIFESTS_DIR / "rcsb_source_taxonomy_audit.csv"
UNIPROT_AUDIT_CSV = MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv"
DISCOVERY_BY_RECEPTOR_CSV = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
SEQ_COORD_MANIFEST_CSV = MANIFESTS_DIR / "rcsb_sequence_coordinates_manifest.csv"
SEARCH_MANIFEST_CSV = MANIFESTS_DIR / "rcsb_search_manifest.csv"
SUMMARY_JSON = MANIFESTS_DIR / "rcsb_discovery_summary.json"
PROVENANCE_JSON = MANIFESTS_DIR / "rcsb_stage2_provenance.json"
REPORT_MD = REPORTS_TABLES_DIR / "rcsb_discovery_report.md"

COVERAGE_VALIDATION_CSV = MANIFESTS_DIR / "rcsb_coverage_validation.csv"
SHORT_MAPPING_AUDIT_CSV = MANIFESTS_DIR / "rcsb_short_mapping_audit.csv"

MATCH_TOLERANCE = 1e-4
MINOR_DIFFERENCE_TOLERANCE = 1e-2
SHORT_COVERAGE_THRESHOLD = 0.05
SHORT_LENGTH_THRESHOLD = 30


def json_field(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def load_reference_lengths() -> dict[str, int]:
    df = pd.read_csv(UNIPROT_AUDIT_CSV)
    return {
        row["resolved_primary_accession"]: int(row["sequence_length"])
        for _, row in df.iterrows()
        if pd.notna(row.get("sequence_length")) and row.get("resolved_primary_accession")
    }


def robustify_mapping_coverage(mapping_df: pd.DataFrame, entity_lengths: dict[str, int], reference_lengths: dict[str, int]) -> pd.DataFrame:
    computed_cols = {
        "aligned_region_count": [],
        "reference_intervals_overlap": [],
        "entity_intervals_overlap": [],
        "computed_reference_sequence_coverage": [],
        "computed_entity_sequence_coverage": [],
    }
    for _, row in mapping_df.iterrows():
        regions = json.loads(row["alignment_regions_json"])
        entity_length = entity_lengths.get(row["polymer_entity_id"])
        ref_length = reference_lengths.get(row["mapped_uniprot_accession"])
        result = compute_union_coverage(regions, entity_length, ref_length)
        for k in computed_cols:
            computed_cols[k].append(result[k])

    out = mapping_df.copy()
    for k, v in computed_cols.items():
        out[k] = v
    # Drop the old (Stage 2.1, naive-sum) coverage columns — superseded.
    out = out.drop(columns=["reference_sequence_coverage", "entity_sequence_coverage"], errors="ignore")
    return out


def build_coverage_validation(
    mapping_df: pd.DataFrame,
    alignments_df: pd.DataFrame,
    discovery_df: pd.DataFrame,
) -> pd.DataFrame:
    # One row per Stage 2 candidate (discovery_df), joined to its own
    # (query accession) computed mapping and official alignment.
    query_mappings = mapping_df[mapping_df["is_query_receptor_accession"]]
    query_mappings = query_mappings.drop_duplicates(subset=["polymer_entity_id", "mapped_uniprot_accession"])

    merged = discovery_df[["uniprot_id", "polymer_entity_id"]].drop_duplicates().merge(
        query_mappings[[
            "polymer_entity_id", "mapped_uniprot_accession",
            "computed_reference_sequence_coverage", "computed_entity_sequence_coverage",
            "reference_intervals_overlap", "entity_intervals_overlap",
        ]],
        left_on=["polymer_entity_id", "uniprot_id"], right_on=["polymer_entity_id", "mapped_uniprot_accession"],
        how="left",
    )
    merged = merged.merge(
        alignments_df[["uniprot_id", "polymer_entity_id", "rcsb_reference_sequence_coverage", "rcsb_entity_sequence_coverage"]],
        on=["uniprot_id", "polymer_entity_id"], how="left",
    )

    rows = []
    for _, row in merged.iterrows():
        rcsb_ref = row["rcsb_reference_sequence_coverage"]
        computed_ref = row["computed_reference_sequence_coverage"]
        rcsb_ent = row["rcsb_entity_sequence_coverage"]
        computed_ent = row["computed_entity_sequence_coverage"]

        ref_diff = abs(rcsb_ref - computed_ref) if pd.notna(rcsb_ref) and pd.notna(computed_ref) else None
        ent_diff = abs(rcsb_ent - computed_ent) if pd.notna(rcsb_ent) and pd.notna(computed_ent) else None

        if ref_diff is None or ent_diff is None:
            status, reason = "REVIEW_NEEDED", (
                f"Missing value(s) for comparison: rcsb_reference={rcsb_ref}, computed_reference={computed_ref}, "
                f"rcsb_entity={rcsb_ent}, computed_entity={computed_ent}."
            )
        elif ref_diff <= MATCH_TOLERANCE and ent_diff <= MATCH_TOLERANCE:
            status, reason = "MATCH", f"Both differences <= {MATCH_TOLERANCE}."
        elif ref_diff <= MINOR_DIFFERENCE_TOLERANCE and ent_diff <= MINOR_DIFFERENCE_TOLERANCE:
            status, reason = "MINOR_NUMERIC_DIFFERENCE", (
                f"reference_abs_difference={ref_diff:.6f}, entity_abs_difference={ent_diff:.6f} "
                f"— within rounding-scale tolerance ({MINOR_DIFFERENCE_TOLERANCE})."
            )
        else:
            status, reason = "REVIEW_NEEDED", (
                f"reference_abs_difference={ref_diff:.6f}, entity_abs_difference={ent_diff:.6f} "
                f"exceeds rounding-scale tolerance ({MINOR_DIFFERENCE_TOLERANCE})."
            )

        rows.append(
            {
                "uniprot_id": row["uniprot_id"],
                "polymer_entity_id": row["polymer_entity_id"],
                "rcsb_reference_sequence_coverage": rcsb_ref,
                "computed_reference_sequence_coverage": computed_ref,
                "reference_abs_difference": ref_diff,
                "rcsb_entity_sequence_coverage": rcsb_ent,
                "computed_entity_sequence_coverage": computed_ent,
                "entity_abs_difference": ent_diff,
                "reference_intervals_overlap": row["reference_intervals_overlap"],
                "entity_intervals_overlap": row["entity_intervals_overlap"],
                "validation_status": status,
                "validation_reason": reason,
            }
        )
    return pd.DataFrame(rows)


def update_candidate_inventory(
    inventory_df: pd.DataFrame,
    mapping_df: pd.DataFrame,
    alignments_df: pd.DataFrame,
) -> pd.DataFrame:
    query_mappings = mapping_df[mapping_df["is_query_receptor_accession"]].drop_duplicates(
        subset=["polymer_entity_id", "mapped_uniprot_accession"]
    )
    computed = query_mappings.rename(columns={"mapped_uniprot_accession": "uniprot_id"})[
        ["polymer_entity_id", "uniprot_id", "computed_reference_sequence_coverage", "computed_entity_sequence_coverage"]
    ]
    official = alignments_df[["polymer_entity_id", "uniprot_id", "rcsb_reference_sequence_coverage", "rcsb_entity_sequence_coverage"]]

    # Drop both the pre-Stage-2.2 column names AND this function's own
    # output columns, so reruns (e.g. Stage 2.3 rerunning Stage 2.2's
    # finalize script) don't collide with columns it wrote last time.
    out = inventory_df.drop(
        columns=[
            "query_uniprot_reference_sequence_coverage", "query_uniprot_entity_sequence_coverage",
            "computed_reference_sequence_coverage", "computed_entity_sequence_coverage",
            "rcsb_reference_sequence_coverage", "rcsb_entity_sequence_coverage",
        ],
        errors="ignore",
    )
    out = out.merge(computed, on=["polymer_entity_id", "uniprot_id"], how="left")
    out = out.merge(official, on=["polymer_entity_id", "uniprot_id"], how="left")
    return out


def build_short_mapping_audit(
    discovery_df: pd.DataFrame,
    mapping_df: pd.DataFrame,
    entities_df: pd.DataFrame,
) -> pd.DataFrame:
    query_mappings = mapping_df[mapping_df["is_query_receptor_accession"]].drop_duplicates(
        subset=["polymer_entity_id", "mapped_uniprot_accession"]
    )
    entities_by_id = entities_df.set_index("polymer_entity_id")

    rows = []
    for _, row in query_mappings.iterrows():
        regions = json.loads(row["alignment_regions_json"])
        entity_union_length = sum(r["entity_end_seq_id"] - r["entity_beg_seq_id"] + 1 for r in regions) if regions else 0
        # Re-derive union (not naive sum) length for reporting mapped_length.
        merged_entity_intervals = []
        for beg, end in sorted((r["entity_beg_seq_id"], r["entity_end_seq_id"]) for r in regions):
            if merged_entity_intervals and beg <= merged_entity_intervals[-1][1] + 1:
                merged_entity_intervals[-1] = (merged_entity_intervals[-1][0], max(merged_entity_intervals[-1][1], end))
            else:
                merged_entity_intervals.append((beg, end))
        mapped_length = sum(e - b + 1 for b, e in merged_entity_intervals)

        ref_coverage = row["computed_reference_sequence_coverage"]
        is_short = (pd.notna(ref_coverage) and ref_coverage < SHORT_COVERAGE_THRESHOLD) or (mapped_length < SHORT_LENGTH_THRESHOLD)
        if not is_short:
            continue

        entity_row = entities_by_id.loc[row["polymer_entity_id"]] if row["polymer_entity_id"] in entities_by_id.index else None
        disc_row = discovery_df[
            (discovery_df["polymer_entity_id"] == row["polymer_entity_id"]) & (discovery_df["uniprot_id"] == row["query_project_uniprot"])
        ]
        if len(disc_row) == 0:
            continue
        disc_row = disc_row.iloc[0]

        rows.append(
            {
                "uniprot_id": disc_row["uniprot_id"],
                "nr_code": disc_row["nr_code"],
                "common_name": disc_row["common_name"],
                "pdb_id": disc_row["pdb_id"],
                "polymer_entity_id": row["polymer_entity_id"],
                "entity_sequence_length": entity_row["entity_sequence_length"] if entity_row is not None else "",
                "mapped_length": mapped_length,
                "reference_coverage": ref_coverage,
                "entity_coverage": row["computed_entity_sequence_coverage"],
                "entity_description": entity_row["entity_description"] if entity_row is not None else "",
            }
        )

    columns = [
        "uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id",
        "entity_sequence_length", "mapped_length", "reference_coverage", "entity_coverage", "entity_description",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows)[columns].sort_values("mapped_length")


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def update_report_with_stage_2_2_section(stage_2_2_summary: dict) -> None:
    """Append a Stage 2.2 section to the existing report rather than
    rewriting the Stage 2 / Stage 2.1 historical sections above it."""
    existing_text = REPORT_MD.read_text()
    marker = "\n## Stage 2.2 — Official Sequence Coordinates Validation\n"
    if marker in existing_text:
        existing_text = existing_text.split(marker)[0].rstrip() + "\n"

    lines = [existing_text.rstrip(), "", marker.strip(), ""]
    lines.append(
        "Three official, distinct RCSB services were used across Stage 2 / 2.1 / 2.2 / 2.3, each for a "
        "different purpose:\n"
        "- **Search API** (`search.rcsb.org/rcsbsearch/v2/query`) = candidate discovery (Stage 2)\n"
        "- **Data API** (`data.rcsb.org/graphql`) = entity/entry metadata enrichment (Stage 2, 2.1)\n"
        "- **Sequence Coordinates API** (`sequence-coordinates.rcsb.org/graphql`) = independent "
        "sequence-alignment/coverage validation (Stage 2.2)\n"
    )
    lines.append(
        f"All 48 project UniProt accessions were queried once each (UNIPROT -> PDB_ENTITY), rather "
        f"than once per polymer entity. {stage_2_2_summary['total_experimental_target_alignments']} "
        f"experimental-PDB-entity target alignments were returned (computed structure models such as "
        f"AlphaFold were identified by target ID and excluded); "
        f"{stage_2_2_summary['stage2_candidates_matched']} of the 2072 Stage 2 candidates were matched, "
        f"{stage_2_2_summary['stage2_candidates_missing']} were missing, and "
        f"{stage_2_2_summary['targets_outside_stage2_scope_found']} target(s) fell outside the Stage 2 "
        f"Search API candidate scope. Stage 2.3 inspected these against authoritative RCSB Data API "
        f"metadata (not titles) and found they are legitimately out of scope: their "
        f"`structure_determination_methodology` is `integrative`, not `experimental` — Stage 2's Search "
        f"query explicitly restricts to `experimental`, while the Sequence Coordinates API applies no "
        f"such filter. They were NOT recently deposited (2018/2019) and their presence is not evidence "
        f"of any indexing-timing discrepancy between RCSB services. See "
        f"`rcsb_sequence_coordinates_extra_targets.csv` for the full per-target classification; they "
        f"were not added to the Stage 2 candidate inventory.\n"
    )
    lines.append(
        f"Coverage agreement (official RCSB `query_coverage`/`target_coverage` vs. this project's "
        f"independently computed union-interval coverage): "
        f"MATCH={stage_2_2_summary['coverage_validation_match']}, "
        f"MINOR_NUMERIC_DIFFERENCE={stage_2_2_summary['coverage_validation_minor_numeric_difference']}, "
        f"REVIEW_NEEDED={stage_2_2_summary['coverage_validation_review_needed']} "
        f"(max reference diff={stage_2_2_summary['max_reference_abs_difference']:.2e}, "
        f"max entity diff={stage_2_2_summary['max_entity_abs_difference']:.2e}). "
        f"Reference-interval overlaps found: {stage_2_2_summary['mappings_with_overlapping_reference_intervals']}; "
        f"entity-interval overlaps found: {stage_2_2_summary['mappings_with_overlapping_entity_intervals']} "
        f"(out of 2118 normalized mapping rows checked) — where zero, the interval-union computation "
        f"reduces exactly to the prior simple sum, i.e. the Stage 2.1 naive-sum coverage was already "
        f"numerically correct for every such row; the union logic is retained as the robust general "
        f"method going forward.\n"
    )
    lines.append(
        f"{stage_2_2_summary['short_mappings_length_below_30aa']} candidates were flagged as short "
        f"mappings (reference coverage < 0.05 and/or aligned length < 30 residues) in "
        f"`rcsb_short_mapping_audit.csv` — informational only, for Stage 3 to identify peptide/"
        f"fragment constructs that are not full LBD structures.\n"
    )
    lines.append(
        "`rcsb_candidate_inventory.csv` now carries BOTH official "
        "(`rcsb_reference_sequence_coverage`, `rcsb_entity_sequence_coverage`) and independently "
        "computed (`computed_reference_sequence_coverage`, `computed_entity_sequence_coverage`) "
        "coverage fields side by side; the official value is null (never guessed) where the Sequence "
        "Coordinates API did not return a match.\n"
    )
    lines.append(
        "The Stage 2.1 taxonomy audit's `EXPLAINED` status describes only that the unusual source "
        "annotation has been characterized, NOT that the structure is suitable for the final human "
        "NR B-factor dataset — all 9 taxonomy-audit rows now carry `qc_suitability = NOT_EVALUATED`; "
        "no inclusion/exclusion decision has been made.\n"
    )
    lines.append("No final structural QC filters were applied at Stage 2, 2.1, 2.2, or 2.3.\n")
    lines.append("No mmCIF coordinate files were downloaded at Stage 2, 2.1, 2.2, or 2.3.\n")

    REPORT_MD.write_text("\n".join(lines))


def update_provenance(stage_2_2_summary: dict) -> None:
    provenance = json.loads(PROVENANCE_JSON.read_text())
    provenance["description"] = (
        "Stage 2 (discovery) + Stage 2.1 (coverage-field correction, source-taxonomy and "
        "multi-UniProt-mapping audit) + Stage 2.2 (independent validation of sequence coverage "
        "against the official RCSB Sequence Coordinates API, and an interval-union robustness fix). "
        "No final structural QC filtering, no coordinate downloads."
    )
    provenance["rcsb_sequence_coordinates_api_endpoint"] = "https://sequence-coordinates.rcsb.org/graphql"
    provenance["stage_2_2_sequence_coordinates_validation"] = stage_2_2_summary
    provenance["generated_manifests"] = {
        "rcsb_search_manifest.csv": sha256_of_file(SEARCH_MANIFEST_CSV),
        "rcsb_discovery_by_receptor.csv": sha256_of_file(DISCOVERY_BY_RECEPTOR_CSV),
        "rcsb_polymer_entities.csv": sha256_of_file(POLYMER_ENTITIES_CSV),
        "rcsb_candidate_inventory.csv": sha256_of_file(CANDIDATE_INVENTORY_CSV),
        "rcsb_polymer_entity_uniprot_mappings.csv": sha256_of_file(MAPPINGS_CSV),
        "rcsb_source_taxonomy_audit.csv": sha256_of_file(TAXONOMY_AUDIT_CSV),
        "rcsb_sequence_coordinates_manifest.csv": sha256_of_file(SEQ_COORD_MANIFEST_CSV),
        "rcsb_sequence_coordinate_alignments.csv": sha256_of_file(ALIGNMENTS_CSV),
        "rcsb_coverage_validation.csv": sha256_of_file(COVERAGE_VALIDATION_CSV),
        "rcsb_short_mapping_audit.csv": sha256_of_file(SHORT_MAPPING_AUDIT_CSV),
    }
    provenance["generating_scripts"] = [
        "scripts/04_discover_rcsb_candidates.py",
        "scripts/05_fetch_rcsb_metadata.py",
        "scripts/06_build_rcsb_inventory.py",
        "scripts/07_correct_rcsb_coverage_and_audit.py",
        "scripts/08_fetch_sequence_coordinates.py",
        "scripts/09_finalize_stage2_coverage.py",
    ]
    provenance["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    PROVENANCE_JSON.write_text(json.dumps(provenance, indent=2) + "\n")


def main() -> None:
    mapping_df = pd.read_csv(MAPPINGS_CSV)
    alignments_df = pd.read_csv(ALIGNMENTS_CSV, dtype={"pdb_id": str})
    entities_df = pd.read_csv(POLYMER_ENTITIES_CSV, dtype={"pdb_id": str})
    inventory_df = pd.read_csv(CANDIDATE_INVENTORY_CSV, dtype={"pdb_id": str})
    taxonomy_audit_df = pd.read_csv(TAXONOMY_AUDIT_CSV)
    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR_CSV, dtype={"pdb_id": str, "entity_id": str})

    reference_lengths = load_reference_lengths()
    entity_lengths = dict(zip(entities_df["polymer_entity_id"], entities_df["entity_sequence_length"]))

    print("Robustifying mapping coverage with interval UNION logic...")
    mapping_df = robustify_mapping_coverage(mapping_df, entity_lengths, reference_lengths)
    n_ref_overlap = int(mapping_df["reference_intervals_overlap"].sum())
    n_ent_overlap = int(mapping_df["entity_intervals_overlap"].sum())
    if n_ref_overlap == 0 and n_ent_overlap == 0:
        print("  No overlapping aligned regions found in ANY of the "
              f"{len(mapping_df)} mapping rows — union coverage equals the "
              "prior naive-sum coverage for every row.")
    else:
        print(f"  Overlapping intervals found: reference={n_ref_overlap}, entity={n_ent_overlap} "
              f"— union coverage differs from naive sum for these rows.")
    mapping_df.to_csv(MAPPINGS_CSV, index=False)
    print(f"Updated {MAPPINGS_CSV.relative_to(PROJECT_ROOT)}")

    print("\nBuilding coverage validation (official vs computed)...")
    validation_df = build_coverage_validation(mapping_df, alignments_df, discovery_df)
    validation_df.to_csv(COVERAGE_VALIDATION_CSV, index=False)
    print(f"Wrote {COVERAGE_VALIDATION_CSV.relative_to(PROJECT_ROOT)} ({len(validation_df)} rows)")
    print(validation_df["validation_status"].value_counts().to_string())

    print("\nUpdating candidate inventory with official + computed coverage fields...")
    inventory_df = update_candidate_inventory(inventory_df, mapping_df, alignments_df)
    inventory_df.to_csv(CANDIDATE_INVENTORY_CSV, index=False)
    print(f"Updated {CANDIDATE_INVENTORY_CSV.relative_to(PROJECT_ROOT)} ({len(inventory_df)} rows)")

    print("\nAdding qc_suitability = NOT_EVALUATED to taxonomy audit...")
    taxonomy_audit_df["qc_suitability"] = "NOT_EVALUATED"
    taxonomy_audit_df.to_csv(TAXONOMY_AUDIT_CSV, index=False)
    print(f"Updated {TAXONOMY_AUDIT_CSV.relative_to(PROJECT_ROOT)} ({len(taxonomy_audit_df)} rows)")

    print("\nBuilding short-mapping audit...")
    short_mapping_df = build_short_mapping_audit(discovery_df, mapping_df, entities_df)
    short_mapping_df.to_csv(SHORT_MAPPING_AUDIT_CSV, index=False)
    print(f"Wrote {SHORT_MAPPING_AUDIT_CSV.relative_to(PROJECT_ROOT)} ({len(short_mapping_df)} rows)")

    # --- update summary json (add Stage 2.2 section) ---
    summary = json.loads(SUMMARY_JSON.read_text())
    seq_coord_manifest_df = pd.read_csv(SEQ_COORD_MANIFEST_CSV)
    status_counts = validation_df["validation_status"].value_counts().to_dict()
    summary["stage_2_2_sequence_coordinates_validation"] = {
        "successful_receptor_queries": int(seq_coord_manifest_df["success"].sum()),
        "failed_receptor_queries": int((~seq_coord_manifest_df["success"]).sum()),
        "total_experimental_target_alignments": int(len(alignments_df)),
        "stage2_candidates_matched": int(seq_coord_manifest_df["stage2_targets_matched"].fillna(0).sum()),
        "stage2_candidates_missing": int(
            seq_coord_manifest_df["missing_stage2_targets"].apply(lambda v: len(json.loads(v)) if isinstance(v, str) and v else 0).sum()
        ),
        "targets_outside_stage2_scope_found": int(
            seq_coord_manifest_df["targets_outside_stage2_scope"].apply(lambda v: len(json.loads(v)) if isinstance(v, str) and v else 0).sum()
        ),
        "coverage_validation_match": int(status_counts.get("MATCH", 0)),
        "coverage_validation_minor_numeric_difference": int(status_counts.get("MINOR_NUMERIC_DIFFERENCE", 0)),
        "coverage_validation_review_needed": int(status_counts.get("REVIEW_NEEDED", 0)),
        "max_reference_abs_difference": float(validation_df["reference_abs_difference"].max(skipna=True)) if validation_df["reference_abs_difference"].notna().any() else None,
        "max_entity_abs_difference": float(validation_df["entity_abs_difference"].max(skipna=True)) if validation_df["entity_abs_difference"].notna().any() else None,
        "mappings_with_overlapping_reference_intervals": n_ref_overlap,
        "mappings_with_overlapping_entity_intervals": n_ent_overlap,
        "short_mappings_reference_coverage_below_0_05": int((short_mapping_df["reference_coverage"] < SHORT_COVERAGE_THRESHOLD).sum()) if len(short_mapping_df) else 0,
        "short_mappings_length_below_30aa": int((short_mapping_df["mapped_length"] < SHORT_LENGTH_THRESHOLD).sum()) if len(short_mapping_df) else 0,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"\nUpdated {SUMMARY_JSON.relative_to(PROJECT_ROOT)}")

    print("\n--- Stage 2.2 validation summary ---")
    for k, v in summary["stage_2_2_sequence_coordinates_validation"].items():
        print(f"{k}: {v}")

    update_report_with_stage_2_2_section(summary["stage_2_2_sequence_coordinates_validation"])
    print(f"\nUpdated {REPORT_MD.relative_to(PROJECT_ROOT)}")

    update_provenance(summary["stage_2_2_sequence_coordinates_validation"])
    print(f"Updated {PROVENANCE_JSON.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
