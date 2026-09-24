"""
Stage 3A, step 2: map every Stage 2 candidate (receptor x polymer entity)
against its receptor's authoritative LBD (and, diagnostically, DBD)
reference interval, using the official RCSB Sequence Coordinates alignment
data already validated in Stage 2.2.

Reads:
    data/manifests/nr_lbd_reference.csv
    data/manifests/rcsb_sequence_coordinate_alignments.csv   (query_begin/end are
        already in the QUERYING receptor's own UniProt coordinates — no
        cross-fusion-partner ambiguity, since each row was produced by a
        UNIPROT -> PDB_ENTITY query for that specific accession)
    data/manifests/rcsb_discovery_by_receptor.csv
    data/manifests/rcsb_short_mapping_audit.csv

Writes:
    data/manifests/rcsb_candidate_lbd_mapping.csv    (2072 rows)
    data/manifests/rcsb_short_mapping_lbd_audit.csv

IMPORTANT: `lbd_mapping_coverage` is a SEQUENCE MAPPING / CONSTRUCT
coverage metric — the fraction of the canonical LBD interval that this
polymer entity's construct/alignment spans. It says nothing about how many
of those residues have resolved atomic coordinates in the deposited
structure (that requires the mmCIF file, not fetched at Stage 3A).

Stage 3A.2 policy (Option A): `lbd_start`/`lbd_end` here are the
STANDARDIZED project-wide interval (UniProt/PROSITE NR-LBD profile), now
labeled explicitly via `lbd_definition`/`lbd_coordinate_system` on every
row. This is the sole interval used for `lbd_mapping_category` and
downstream prefiltering — the previous cohort's legacy VDR-only window is
tracked separately (data/manifests/nr_lbd_boundary_audit.csv,
rcsb_lbd_boundary_sensitivity.csv) and never substituted here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.intervals import intersect_with_reference, merge_intervals, total_length  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
LBD_REFERENCE_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"
ALIGNMENTS_CSV = MANIFESTS_DIR / "rcsb_sequence_coordinate_alignments.csv"
DISCOVERY_CSV = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
SHORT_MAPPING_CSV = MANIFESTS_DIR / "rcsb_short_mapping_audit.csv"

OUTPUT_CSV = MANIFESTS_DIR / "rcsb_candidate_lbd_mapping.csv"
SHORT_MAPPING_LBD_CSV = MANIFESTS_DIR / "rcsb_short_mapping_lbd_audit.csv"

LBD_DEFINITION = "UNIPROT_PROSITE_NR_LBD"
LBD_COORDINATE_SYSTEM = "UNIPROT_CANONICAL_1_BASED_INCLUSIVE"

OUTPUT_COLUMNS = [
    "uniprot_id", "nr_code", "common_name", "group",
    "pdb_id", "polymer_entity_id", "entity_id",
    "lbd_start", "lbd_end", "lbd_length", "lbd_definition", "lbd_coordinate_system",
    "mapped_reference_residues_total", "mapped_lbd_residues", "lbd_mapping_coverage",
    "mapped_non_lbd_residues",
    "mapping_overlaps_lbd", "mapping_fully_spans_lbd_reference_interval",
    "mapped_dbd_residues", "dbd_mapping_coverage", "mapping_overlaps_dbd",
    "lbd_mapping_category",
    "alignment_regions_json",
    "mapping_status", "mapping_notes",
]


def categorize(coverage: float) -> str:
    if coverage == 0:
        return "NO_LBD_OVERLAP"
    if coverage < 0.25:
        return "TRACE_LBD"
    if coverage < 0.80:
        return "PARTIAL_LBD"
    if coverage < 0.95:
        return "SUBSTANTIAL_LBD"
    return "NEAR_COMPLETE_LBD"


def main() -> None:
    lbd_ref = pd.read_csv(LBD_REFERENCE_CSV).set_index("uniprot_id")
    alignments = pd.read_csv(ALIGNMENTS_CSV, dtype={"pdb_id": str, "entity_id": str})
    alignments = alignments[alignments["is_stage2_candidate"]]
    discovery = pd.read_csv(DISCOVERY_CSV, dtype={"pdb_id": str, "entity_id": str})

    alignments_by_key = alignments.set_index(["uniprot_id", "polymer_entity_id"])

    rows = []
    for _, disc_row in discovery.iterrows():
        uniprot_id = disc_row["uniprot_id"]
        entity_id_full = disc_row["polymer_entity_id"]
        key = (uniprot_id, entity_id_full)

        ref_row = lbd_ref.loc[uniprot_id] if uniprot_id in lbd_ref.index else None
        align_row = alignments_by_key.loc[key] if key in alignments_by_key.index else None

        if ref_row is None or align_row is None:
            rows.append(
                {
                    "uniprot_id": uniprot_id, "nr_code": disc_row["nr_code"],
                    "common_name": disc_row["common_name"], "group": disc_row["group"],
                    "pdb_id": disc_row["pdb_id"], "polymer_entity_id": entity_id_full, "entity_id": disc_row["entity_id"],
                    "lbd_start": "", "lbd_end": "", "lbd_length": "",
                    "lbd_definition": LBD_DEFINITION, "lbd_coordinate_system": LBD_COORDINATE_SYSTEM,
                    "mapped_reference_residues_total": "", "mapped_lbd_residues": "", "lbd_mapping_coverage": "",
                    "mapped_non_lbd_residues": "",
                    "mapping_overlaps_lbd": "", "mapping_fully_spans_lbd_reference_interval": "",
                    "mapped_dbd_residues": "", "dbd_mapping_coverage": "", "mapping_overlaps_dbd": "",
                    "lbd_mapping_category": "REVIEW_NEEDED",
                    "alignment_regions_json": "[]",
                    "mapping_status": "REVIEW_NEEDED",
                    "mapping_notes": (
                        "Missing LBD reference row" if ref_row is None else
                        "No matching official Sequence Coordinates alignment found for this "
                        "(uniprot_id, polymer_entity_id) pair."
                    ),
                }
            )
            continue

        regions = json.loads(align_row["aligned_regions_json"])
        query_intervals = [(r["query_begin"], r["query_end"]) for r in regions]
        merged_query = merge_intervals(query_intervals)
        mapped_reference_total = total_length(merged_query)

        lbd_start, lbd_end, lbd_length = ref_row["lbd_start"], ref_row["lbd_end"], ref_row["lbd_length"]
        lbd_overlap_intervals = intersect_with_reference(merged_query, int(lbd_start), int(lbd_end))
        mapped_lbd_residues = total_length(lbd_overlap_intervals)
        lbd_coverage = round(mapped_lbd_residues / lbd_length, 6) if lbd_length else 0.0
        mapped_non_lbd = mapped_reference_total - mapped_lbd_residues

        overlaps_lbd = mapped_lbd_residues > 0
        # "Fully spans" means the LBD reference interval is entirely
        # contained within the union of aligned query intervals (not that
        # every residue in it has resolved coordinates).
        fully_spans_lbd = any(
            beg <= lbd_start and end >= lbd_end for beg, end in merged_query
        )

        dbd_status = ref_row["dbd_status"]
        if dbd_status == "PASS":
            dbd_start, dbd_end = int(ref_row["dbd_start"]), int(ref_row["dbd_end"])
            dbd_overlap_intervals = intersect_with_reference(merged_query, dbd_start, dbd_end)
            mapped_dbd_residues = total_length(dbd_overlap_intervals)
            dbd_length = dbd_end - dbd_start + 1
            dbd_coverage = round(mapped_dbd_residues / dbd_length, 6) if dbd_length else 0.0
            overlaps_dbd = mapped_dbd_residues > 0
        else:
            mapped_dbd_residues = ""
            dbd_coverage = ""
            overlaps_dbd = ""

        category = categorize(lbd_coverage)

        rows.append(
            {
                "uniprot_id": uniprot_id, "nr_code": disc_row["nr_code"],
                "common_name": disc_row["common_name"], "group": disc_row["group"],
                "pdb_id": disc_row["pdb_id"], "polymer_entity_id": entity_id_full, "entity_id": disc_row["entity_id"],
                "lbd_start": lbd_start, "lbd_end": lbd_end, "lbd_length": lbd_length,
                "lbd_definition": LBD_DEFINITION, "lbd_coordinate_system": LBD_COORDINATE_SYSTEM,
                "mapped_reference_residues_total": mapped_reference_total,
                "mapped_lbd_residues": mapped_lbd_residues,
                "lbd_mapping_coverage": lbd_coverage,
                "mapped_non_lbd_residues": mapped_non_lbd,
                "mapping_overlaps_lbd": overlaps_lbd,
                "mapping_fully_spans_lbd_reference_interval": fully_spans_lbd,
                "mapped_dbd_residues": mapped_dbd_residues,
                "dbd_mapping_coverage": dbd_coverage,
                "mapping_overlaps_dbd": overlaps_dbd,
                "lbd_mapping_category": category,
                "alignment_regions_json": json.dumps(regions, separators=(",", ":")),
                "mapping_status": "PASS",
                "mapping_notes": "",
            }
        )

    df = pd.DataFrame(rows)[OUTPUT_COLUMNS]
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {OUTPUT_CSV.relative_to(PROJECT_ROOT)} ({len(df)} rows)")
    print(df["lbd_mapping_category"].value_counts().to_string())

    # --- short-mapping LBD audit ---
    short_df = pd.read_csv(SHORT_MAPPING_CSV, dtype={"pdb_id": str})
    merged = short_df.merge(
        df[[
            "uniprot_id", "polymer_entity_id", "lbd_mapping_coverage", "mapping_overlaps_lbd",
            "dbd_mapping_coverage", "mapping_overlaps_dbd", "lbd_mapping_category",
        ]],
        on=["uniprot_id", "polymer_entity_id"], how="left",
    )
    short_lbd_columns = [
        "uniprot_id", "nr_code", "common_name", "pdb_id", "polymer_entity_id",
        "mapped_length", "reference_coverage",
        "mapping_overlaps_lbd", "lbd_mapping_coverage",
        "mapping_overlaps_dbd", "dbd_mapping_coverage",
        "lbd_mapping_category", "entity_description",
    ]
    merged[short_lbd_columns].to_csv(SHORT_MAPPING_LBD_CSV, index=False)
    print(f"Wrote {SHORT_MAPPING_LBD_CSV.relative_to(PROJECT_ROOT)} ({len(merged)} rows)")
    n_short_overlap_lbd = int(merged["mapping_overlaps_lbd"].sum())
    print(f"Short mappings overlapping LBD: {n_short_overlap_lbd} / {len(merged)}")


if __name__ == "__main__":
    main()
