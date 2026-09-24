"""Exact polymer-entity `label_seq_id` <-> canonical UniProt residue mapper,
built from the already-validated official RCSB Sequence Coordinates
alignment (Stage 2.2), for the QUERY receptor accession only.

Critical design rule (Stage 3B Section 7): canonical residue mapping is
built from `label_seq_id` (entity sequence position), never `auth_seq_id`
(author/depositor numbering, which can be discontinuous, offset, or reused
across chains). `auth_seq_id` is carried through observation tables only as
display/provenance metadata.

The alignment's `target_begin`/`target_end` (from Sequence Coordinates,
`to: PDB_ENTITY`) are entity sequence positions using the SAME 1-based
numbering as `_atom_site.label_seq_id` — RCSB's own polymer entity sequence
numbering, confirmed by construction (the alignment was queried against the
entity's own `pdbx_seq_one_letter_code_can`).
"""

from __future__ import annotations

from typing import Any


def build_entity_to_uniprot_map(aligned_regions: list[dict[str, Any]]) -> dict[str, Any]:
    """Build an unambiguous label_seq_id -> canonical UniProt position map.

    Each region in `aligned_regions` (from
    rcsb_sequence_coordinate_alignments.csv `aligned_regions_json`) has
    query_begin/query_end (UniProt) and target_begin/target_end (entity
    label_seq_id), 1-based inclusive, from the SAME alignment record.

    Returns:
        {
            "mapping": {label_seq_id (int): uniprot_position (int), ...},
            "status": "OK" | "REVIEW_NEEDED",
            "notes": str,
        }
    Never interpolates across gaps or extrapolates outside an aligned
    region. If the same label_seq_id would map to more than one UniProt
    position (ambiguous), or a region's reference/entity aligned lengths
    disagree, the whole map is marked REVIEW_NEEDED and returned empty
    (callers must not silently use a partial/ambiguous map).
    """
    mapping: dict[int, int] = {}
    notes: list[str] = []

    for region in aligned_regions:
        q_beg, q_end = region.get("query_begin"), region.get("query_end")
        t_beg, t_end = region.get("target_begin"), region.get("target_end")

        if q_beg is None or q_end is None or t_beg is None or t_end is None:
            return {"mapping": {}, "status": "REVIEW_NEEDED", "notes": f"Incomplete region: {region}"}
        if q_end < q_beg or t_end < t_beg:
            return {"mapping": {}, "status": "REVIEW_NEEDED", "notes": f"Invalid interval order in region: {region}"}

        ref_length = q_end - q_beg + 1
        entity_length = t_end - t_beg + 1
        if ref_length != entity_length:
            return {
                "mapping": {}, "status": "REVIEW_NEEDED",
                "notes": f"Reference-aligned length ({ref_length}) != entity-aligned length ({entity_length}) "
                         f"for region {region} — cannot map 1:1 without guessing indel placement.",
            }

        for offset in range(ref_length):
            label_seq_id = t_beg + offset
            uniprot_pos = q_beg + offset
            if label_seq_id in mapping and mapping[label_seq_id] != uniprot_pos:
                return {
                    "mapping": {}, "status": "REVIEW_NEEDED",
                    "notes": f"Ambiguous mapping: label_seq_id {label_seq_id} maps to both "
                             f"{mapping[label_seq_id]} and {uniprot_pos} across overlapping regions.",
                }
            mapping[label_seq_id] = uniprot_pos

    return {"mapping": mapping, "status": "OK", "notes": "; ".join(notes)}


def invert_mapping(entity_to_uniprot: dict[int, int]) -> tuple[dict[int, int], set[int]]:
    """uniprot_position -> label_seq_id, plus the set of any UniProt
    positions claimed by more than one entity label_seq_id (should not
    happen for a well-formed 1:1 alignment, but is detected and reported
    explicitly rather than silently resolved by overwriting)."""
    inverted: dict[int, int] = {}
    ambiguous: set[int] = set()
    for label_seq_id, uniprot_pos in sorted(entity_to_uniprot.items()):
        if uniprot_pos in inverted:
            ambiguous.add(uniprot_pos)
            continue
        inverted[uniprot_pos] = label_seq_id
    return inverted, ambiguous
