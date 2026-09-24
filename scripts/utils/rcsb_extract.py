"""Field extraction, flag computation, and mapping-verification logic for
raw RCSB Data API GraphQL responses. Defensive: every extractor returns
None/empty when a field is absent, never inventing a value.
"""

from __future__ import annotations

from typing import Any

from scripts.utils.intervals import merge_intervals as _merge_intervals
from scripts.utils.intervals import overlap_fraction as _interval_overlap_fraction

EXPERIMENTAL_METHOD_XRAY_TOKENS = {"X-RAY DIFFRACTION"}
EXPERIMENTAL_METHOD_CRYOEM_TOKENS = {"ELECTRON MICROSCOPY"}
EXPERIMENTAL_METHOD_NMR_TOKENS = {
    "SOLUTION NMR", "SOLID-STATE NMR",
}
HOMO_SAPIENS_TAXON_ID = 9606


def parse_polymer_entity_identifier(identifier: str) -> tuple[str, str]:
    """'1DB1_1' -> ('1DB1', '1'). Raises ValueError if the identifier does
    not have the expected <PDBID>_<ENTITY_ID> form."""
    if "_" not in identifier:
        raise ValueError(f"Polymer entity identifier missing '_': {identifier!r}")
    pdb_id, entity_id = identifier.rsplit("_", 1)
    if not pdb_id or not entity_id:
        raise ValueError(f"Could not parse polymer entity identifier: {identifier!r}")
    return pdb_id, entity_id


def extract_polymer_entity_fields(raw_entity: dict[str, Any]) -> dict[str, Any]:
    """Extract Stage-2-required fields from one polymer_entities[] element
    of a raw RCSB Data API GraphQL response. `raw_entity` may be None if
    the Data API returned null for a requested ID (e.g. suspended/obsolete
    entity) — callers must check for that before calling this function.

    Does NOT compute coverage — see extract_uniprot_mappings() (Stage 2.1)
    for the per-mapping alignment/coverage table. `all_mapped_uniprot_accessions`
    here is the flat identity list only.
    """
    container = raw_entity.get("rcsb_polymer_entity_container_identifiers") or {}
    entity_poly = raw_entity.get("entity_poly") or {}
    polymer_entity = raw_entity.get("rcsb_polymer_entity") or {}
    source_organisms = raw_entity.get("rcsb_entity_source_organism") or []

    ref_seq_ids = container.get("reference_sequence_identifiers") or []
    mapped_accessions = sorted(
        {
            r["database_accession"]
            for r in ref_seq_ids
            if r.get("database_name") == "UniProt" and r.get("database_accession")
        }
    )

    source_organism_names = sorted({o["scientific_name"] for o in source_organisms if o.get("scientific_name")})
    source_taxonomy_ids = sorted({o["ncbi_taxonomy_id"] for o in source_organisms if o.get("ncbi_taxonomy_id") is not None})

    return {
        "entry_id": container.get("entry_id"),
        "entity_id": container.get("entity_id"),
        "entity_description": polymer_entity.get("pdbx_description"),
        "polymer_type": entity_poly.get("type"),
        "entity_sequence_length": entity_poly.get("rcsb_sample_sequence_length"),
        "label_asym_ids": container.get("asym_ids") or [],
        "auth_asym_ids": container.get("auth_asym_ids") or [],
        "all_mapped_uniprot_accessions": mapped_accessions,
        "source_organism_names": source_organism_names,
        "source_taxonomy_ids": source_taxonomy_ids,
    }


def extract_uniprot_mappings(raw_entity: dict[str, Any], entity_sequence_length: int | None) -> list[dict[str, Any]]:
    """Stage 2.1 correction: one row per (polymer_entity, UniProt reference
    mapping), from `rcsb_polymer_entity_align` — the only field in the RCSB
    Data API schema that carries per-mapping alignment detail (confirmed via
    GraphQL introspection: `RcsbPolymerEntityAlign` has no precomputed
    coverage-fraction field, only `aligned_regions` with
    entity_beg_seq_id/ref_beg_seq_id/length).

    For each region we store the raw begin+length exactly as returned, and
    additionally derive entity_end_seq_id = entity_beg_seq_id + length - 1
    (and the equivalent for ref_end_seq_id) — unambiguous from begin+length,
    not a guess.

    entity_sequence_coverage = (sum of aligned region lengths for this
    mapping) / entity_sequence_length — always computable, since
    entity_sequence_length comes from this same entity's own record.

    reference_sequence_coverage is NOT computed here (this function has no
    access to the reference accession's total sequence length); the caller
    fills it in when the reference accession's length is known (i.e. for
    the project's own 48 validated receptors, from the Stage 1 UniProt
    audit) and leaves it None otherwise — never guessed.
    """
    align_records = raw_entity.get("rcsb_polymer_entity_align") or []
    mappings = []
    for a in align_records:
        raw_regions = a.get("aligned_regions") or []
        regions = []
        total_aligned_length = 0
        for r in raw_regions:
            length = r.get("length")
            entity_beg = r.get("entity_beg_seq_id")
            ref_beg = r.get("ref_beg_seq_id")
            region = {
                "entity_beg_seq_id": entity_beg,
                "ref_beg_seq_id": ref_beg,
                "length": length,
                "entity_end_seq_id": (entity_beg + length - 1) if (entity_beg is not None and length is not None) else None,
                "ref_end_seq_id": (ref_beg + length - 1) if (ref_beg is not None and length is not None) else None,
            }
            regions.append(region)
            if length is not None:
                total_aligned_length += length

        entity_coverage = (
            round(total_aligned_length / entity_sequence_length, 6)
            if entity_sequence_length else None
        )

        mappings.append(
            {
                "mapped_database_name": a.get("reference_database_name"),
                "mapped_uniprot_accession": a.get("reference_database_accession"),
                "aligned_regions": regions,
                "total_aligned_length": total_aligned_length if raw_regions else None,
                "entity_sequence_coverage": entity_coverage,
            }
        )
    return mappings


def compute_mapping_status(query_uniprot_id: str, extracted_entity: dict[str, Any]) -> tuple[bool, str, str]:
    """Returns (mapping_contains_query_uniprot, mapping_status, mapping_notes)."""
    mapped = extracted_entity["all_mapped_uniprot_accessions"]
    contains_query = query_uniprot_id in mapped

    if contains_query and len(mapped) == 1:
        return True, "CONFIRMED", ""
    if contains_query and len(mapped) > 1:
        return (
            True,
            "CONFIRMED_MULTI_MAPPED",
            f"Entity maps to multiple UniProt accessions: {mapped}",
        )
    if not contains_query and mapped:
        return (
            False,
            "REVIEW_NEEDED",
            f"Entity was discovered via query accession {query_uniprot_id!r} but its Data API "
            f"metadata maps only to {mapped} — expected mapping not confirmed.",
        )
    return (
        False,
        "REVIEW_NEEDED",
        f"Entity was discovered via query accession {query_uniprot_id!r} but Data API "
        f"metadata reports no UniProt reference-sequence mapping at all.",
    )


def extract_entry_fields(raw_entry: dict[str, Any]) -> dict[str, Any]:
    """Extract Stage-2-required fields from one entries[] element of a raw
    RCSB Data API GraphQL response. `raw_entry` may be None if the Data API
    returned null for a requested ID."""
    struct = raw_entry.get("struct") or {}
    entry_info = raw_entry.get("rcsb_entry_info") or {}
    exptl = raw_entry.get("exptl") or []
    refine = raw_entry.get("refine") or []
    accession_info = raw_entry.get("rcsb_accession_info") or {}
    symmetry = raw_entry.get("symmetry") or {}
    citation = raw_entry.get("rcsb_primary_citation") or {}

    experimental_methods = sorted({e["method"] for e in exptl if e.get("method")})
    r_work_values = [r["ls_R_factor_R_work"] for r in refine if r.get("ls_R_factor_R_work") is not None]
    r_free_values = [r["ls_R_factor_R_free"] for r in refine if r.get("ls_R_factor_R_free") is not None]

    return {
        "structure_title": struct.get("title"),
        "structure_determination_methodology": entry_info.get("structure_determination_methodology"),
        "experimental_methods": experimental_methods,
        "resolution_combined": entry_info.get("resolution_combined") or [],
        "r_work": r_work_values,
        "r_free": r_free_values,
        "deposit_date": accession_info.get("deposit_date"),
        "initial_release_date": accession_info.get("initial_release_date"),
        "revision_date": accession_info.get("revision_date"),
        "space_group": symmetry.get("space_group_name_H_M"),
        "polymer_entity_count": entry_info.get("polymer_entity_count"),
        "protein_entity_count": entry_info.get("polymer_entity_count_protein"),
        "nonpolymer_entity_count": entry_info.get("nonpolymer_entity_count"),
        "deposited_polymer_instance_count": entry_info.get("deposited_polymer_entity_instance_count"),
        "primary_citation_title": citation.get("title"),
        "primary_citation_doi": citation.get("pdbx_database_id_DOI"),
        "primary_citation_year": citation.get("year"),
    }


def compute_union_coverage(
    aligned_regions: list[dict[str, Any]],
    entity_sequence_length: int | None,
    reference_sequence_length: int | None,
) -> dict[str, Any]:
    """Stage 2.2 robustness fix: compute coverage from the UNION of aligned
    region intervals on each coordinate system, not a naive sum of region
    lengths (which double-counts if regions overlap).

    Returns aligned_region_count, whether entity/reference intervals
    overlapped, and the union-based coverage fractions (None where the
    corresponding total length is unknown).
    """
    entity_intervals = [
        (r["entity_beg_seq_id"], r["entity_end_seq_id"])
        for r in aligned_regions
        if r.get("entity_beg_seq_id") is not None and r.get("entity_end_seq_id") is not None
    ]
    ref_intervals = [
        (r["ref_beg_seq_id"], r["ref_end_seq_id"])
        for r in aligned_regions
        if r.get("ref_beg_seq_id") is not None and r.get("ref_end_seq_id") is not None
    ]

    def raw_span_sum(intervals: list[tuple[int, int]]) -> int:
        return sum(end - beg + 1 for beg, end in intervals)

    entity_merged = _merge_intervals(entity_intervals)
    ref_merged = _merge_intervals(ref_intervals)

    entity_union_length = raw_span_sum(entity_merged)
    ref_union_length = raw_span_sum(ref_merged)

    entity_overlap = raw_span_sum(entity_intervals) != entity_union_length
    ref_overlap = raw_span_sum(ref_intervals) != ref_union_length

    return {
        "aligned_region_count": len(aligned_regions),
        "entity_intervals_overlap": entity_overlap,
        "reference_intervals_overlap": ref_overlap,
        "computed_entity_sequence_coverage": (
            round(entity_union_length / entity_sequence_length, 6) if entity_sequence_length else None
        ),
        "computed_reference_sequence_coverage": (
            round(ref_union_length / reference_sequence_length, 6) if reference_sequence_length else None
        ),
    }



def classify_multi_uniprot_entity(
    project_regions: list[dict[str, Any]],
    other_mappings: list[dict[str, Any]],
) -> tuple[str, str]:
    """Conservative classification of a multi-UniProt-mapped polymer entity,
    based ONLY on entity-sequence-coordinate overlap between the project
    receptor's aligned region(s) and each other mapping's aligned region(s)
    — never inferred from entity_description text alone.

    FUSION_OR_CHIMERA: exactly one other mapping, and its entity-sequence
      span is essentially disjoint (<=5% overlap) from the project
      receptor's span — different parts of the same chain map to different
      proteins, consistent with an engineered fusion/chimera construct.
    MULTI_MAPPING_SAME_BIOLOGICAL_PROTEIN: exactly one other mapping, and
      its entity-sequence span overlaps the project receptor's span by
      >=80% — the two accessions are being attached to essentially the same
      stretch of residues (e.g. redundant/alternate reference entries).
    AMBIGUOUS_REVIEW_NEEDED: anything else (more than one other mapping,
      partial overlap, or missing region data).
    """
    if len(other_mappings) != 1:
        return (
            "AMBIGUOUS_REVIEW_NEEDED",
            f"{len(other_mappings)} additional UniProt mapping(s) besides the project receptor "
            f"— classification rule only handles exactly one additional mapping.",
        )

    project_intervals = _merge_intervals(
        [(r["entity_beg_seq_id"], r["entity_end_seq_id"]) for r in project_regions
         if r.get("entity_beg_seq_id") is not None and r.get("entity_end_seq_id") is not None]
    )
    other = other_mappings[0]
    other_intervals = _merge_intervals(
        [(r["entity_beg_seq_id"], r["entity_end_seq_id"]) for r in other["aligned_regions"]
         if r.get("entity_beg_seq_id") is not None and r.get("entity_end_seq_id") is not None]
    )

    if not project_intervals or not other_intervals:
        return (
            "AMBIGUOUS_REVIEW_NEEDED",
            "Missing or unparseable entity-sequence alignment coordinates for one of the mappings.",
        )

    overlap_fraction = _interval_overlap_fraction(project_intervals, other_intervals)

    if overlap_fraction <= 0.05:
        return (
            "FUSION_OR_CHIMERA",
            f"Project receptor's aligned entity-sequence range(s) {project_intervals} and "
            f"{other['mapped_uniprot_accession']}'s range(s) {other_intervals} overlap by only "
            f"{overlap_fraction:.1%} — disjoint regions of the same chain, consistent with an "
            f"engineered fusion/chimera construct.",
        )
    if overlap_fraction >= 0.80:
        return (
            "MULTI_MAPPING_SAME_BIOLOGICAL_PROTEIN",
            f"Project receptor's range(s) {project_intervals} and "
            f"{other['mapped_uniprot_accession']}'s range(s) {other_intervals} overlap by "
            f"{overlap_fraction:.1%} of the shorter span — both accessions map to essentially "
            f"the same residues.",
        )
    return (
        "AMBIGUOUS_REVIEW_NEEDED",
        f"Partial overlap ({overlap_fraction:.1%}) between project receptor range(s) "
        f"{project_intervals} and {other['mapped_uniprot_accession']}'s range(s) {other_intervals} "
        f"— does not clearly fit fusion or same-protein pattern.",
    )


def compute_method_flags(experimental_methods: list[str]) -> dict[str, bool]:
    methods_upper = {m.upper() for m in experimental_methods}
    is_xray = bool(methods_upper & EXPERIMENTAL_METHOD_XRAY_TOKENS)
    is_cryo_em = bool(methods_upper & EXPERIMENTAL_METHOD_CRYOEM_TOKENS)
    is_nmr = bool(methods_upper & EXPERIMENTAL_METHOD_NMR_TOKENS)
    known_tokens = EXPERIMENTAL_METHOD_XRAY_TOKENS | EXPERIMENTAL_METHOD_CRYOEM_TOKENS | EXPERIMENTAL_METHOD_NMR_TOKENS
    is_other_experimental = bool(methods_upper - known_tokens)
    return {
        "is_xray": is_xray,
        "is_cryo_em": is_cryo_em,
        "is_nmr": is_nmr,
        "is_other_experimental": is_other_experimental,
    }
