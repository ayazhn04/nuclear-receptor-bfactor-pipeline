"""Defensive field extraction and conservative identity checks over a raw
UniProtKB JSON entry (rest.uniprot.org/uniprotkb/{accession}.json schema).

Every extractor returns None/"" when a field is absent — nothing here
invents a value that was not present in the API response.
"""

from __future__ import annotations

import re
import unicodedata
from typing import Any

EXPECTED_TAXON_ID = 9606
EXPECTED_ORGANISM = "Homo sapiens"

# Greek-letter <-> Latin-letter equivalents commonly used interchangeably in
# receptor nomenclature (e.g. "ERα" vs "ERA", "PPARγ" vs "PPARG").
_GREEK_TO_LATIN = str.maketrans(
    {
        "α": "A", "β": "B", "γ": "G", "δ": "D", "ε": "E",
    }
)


def _get(d: Any, *path, default=None):
    cur = d
    for key in path:
        if cur is None:
            return default
        if isinstance(key, int):
            if not isinstance(cur, list) or len(cur) <= key:
                return default
            cur = cur[key]
        else:
            if not isinstance(cur, dict) or key not in cur:
                return default
            cur = cur[key]
    return cur if cur is not None else default


def extract_fields(raw: dict[str, Any]) -> dict[str, Any]:
    """Extract the Stage-1-required fields from a raw UniProtKB JSON entry."""
    entry_type = _get(raw, "entryType", default="")
    reviewed_status = (
        "reviewed" if "reviewed" in entry_type.lower() and "unreviewed" not in entry_type.lower()
        else "unreviewed" if "unreviewed" in entry_type.lower()
        else "unknown" if entry_type else ""
    )

    recommended_name = _get(raw, "proteinDescription", "recommendedName", "fullName", "value")
    if recommended_name is None:
        # Unreviewed (TrEMBL) entries often lack recommendedName and instead
        # carry submissionNames.
        recommended_name = _get(raw, "proteinDescription", "submissionNames", 0, "fullName", "value")

    genes = _get(raw, "genes", default=[])
    primary_gene_name = _get(genes, 0, "geneName", "value") if isinstance(genes, list) else None
    gene_synonyms = []
    if isinstance(genes, list) and len(genes) > 0:
        synonyms = genes[0].get("synonyms", []) if isinstance(genes[0], dict) else []
        gene_synonyms = [s.get("value") for s in synonyms if isinstance(s, dict) and s.get("value")]

    sequence_length = _get(raw, "sequence", "length")
    sequence_version = _get(raw, "entryAudit", "sequenceVersion")
    entry_version = _get(raw, "entryAudit", "entryVersion")
    first_public_date = _get(raw, "entryAudit", "firstPublicDate")
    last_sequence_update_date = _get(raw, "entryAudit", "lastSequenceUpdateDate")
    last_annotation_update_date = _get(raw, "entryAudit", "lastAnnotationUpdateDate")

    return {
        "resolved_primary_accession": _get(raw, "primaryAccession"),
        "uniprot_entry_name": _get(raw, "uniProtkbId"),
        "entry_type": entry_type,
        "reviewed_status": reviewed_status,
        "organism_scientific_name": _get(raw, "organism", "scientificName"),
        "organism_common_name": _get(raw, "organism", "commonName"),
        "organism_taxon_id": _get(raw, "organism", "taxonId"),
        "recommended_protein_name": recommended_name,
        "primary_gene_name": primary_gene_name,
        "gene_synonyms": ";".join(gene_synonyms) if gene_synonyms else "",
        "sequence_length": sequence_length,
        "sequence_version": sequence_version,
        "entry_version": entry_version,
        "first_public_date": first_public_date,
        "last_sequence_update_date": last_sequence_update_date,
        "last_annotation_update_date": last_annotation_update_date,
    }


def _normalize_symbol(text: str) -> str:
    """Normalize a receptor symbol/name for loose overlap comparison only.

    Uppercases, strips whitespace/hyphens, and maps Greek letters to their
    common Latin-letter equivalents (α->A etc.) so e.g. "RXRα" and "RXRA"
    compare equal. This is a coarse overlap heuristic, NOT a claim of
    biological equivalence — see check_identity_heuristics().
    """
    text = text.translate(_GREEK_TO_LATIN)
    text = unicodedata.normalize("NFKD", text)
    text = re.sub(r"[^A-Za-z0-9]", "", text)
    return text.upper()


def naming_overlap_found(professor_common_name: str, professor_nr_code: str, extracted: dict[str, Any]) -> bool:
    """Conservative, low-false-positive overlap check between the professor's
    common_name/nr_code and UniProt's gene name / gene synonyms / protein name.

    This intentionally does NOT use exact-string equality against common_name
    to declare a mismatch (per Stage 1 instructions). It only asks whether
    *any* normalized candidate from the professor record appears as a
    normalized substring/equal match against *any* UniProt-supplied gene
    symbol, synonym, or the recommended protein name. Absence of overlap is
    reported as REVIEW_NEEDED (never FAIL) precisely because this heuristic
    can miss legitimate matches (e.g. purely descriptive UniProt names).
    """
    professor_candidates = {
        _normalize_symbol(professor_common_name),
        _normalize_symbol(professor_nr_code),
    }
    professor_candidates = {c for c in professor_candidates if c}

    uniprot_candidates: set[str] = set()
    if extracted.get("primary_gene_name"):
        uniprot_candidates.add(_normalize_symbol(extracted["primary_gene_name"]))
    for syn in (extracted.get("gene_synonyms") or "").split(";"):
        if syn:
            uniprot_candidates.add(_normalize_symbol(syn))

    protein_name_normalized = _normalize_symbol(extracted.get("recommended_protein_name") or "")

    for pc in professor_candidates:
        if not pc:
            continue
        for uc in uniprot_candidates:
            if pc == uc:
                return True
        # Fall back to substring containment against the protein full name
        # only for reasonably specific (non-trivial) candidates, to avoid
        #1-2 character false positives.
        if len(pc) >= 3 and pc in protein_name_normalized:
            return True

    return False


def check_identity(professor_record: dict[str, Any], extracted: dict[str, Any]) -> dict[str, Any]:
    requested_accession = professor_record["uniprot_id"]
    resolved = extracted.get("resolved_primary_accession")

    organism_name = extracted.get("organism_scientific_name") or ""
    taxon_id = extracted.get("organism_taxon_id")
    sequence_length = extracted.get("sequence_length")

    checks = {
        "accession_resolved": resolved is not None,
        "primary_accession_matches_professor": resolved == requested_accession,
        "organism_is_homo_sapiens": organism_name == EXPECTED_ORGANISM,
        "organism_taxon_is_9606": taxon_id == EXPECTED_TAXON_ID,
        "entry_is_reviewed": extracted.get("reviewed_status") == "reviewed",
        "gene_name_available": bool(extracted.get("primary_gene_name")),
        "protein_name_available": bool(extracted.get("recommended_protein_name")),
        "sequence_available": sequence_length is not None,
        "sequence_length_positive": isinstance(sequence_length, int) and sequence_length > 0,
    }
    checks["naming_overlap_found"] = naming_overlap_found(
        professor_record["common_name"], professor_record["nr_code"], extracted
    )
    return checks


def decide_audit_status(
    fetch_success: bool,
    http_status: int | None,
    checks: dict[str, Any] | None,
) -> tuple[str, str]:
    """Return (audit_status, audit_reason). Conservative: only FAIL on clear,
    objective problems; REVIEW_NEEDED for anything requiring a human look."""
    if not fetch_success or checks is None:
        return (
            "FAIL",
            f"UniProt request did not succeed (http_status={http_status}); "
            f"accession could not be resolved after retries.",
        )

    if not checks["accession_resolved"]:
        return "FAIL", "Response received but no primaryAccession present in the record."

    if not checks["organism_taxon_is_9606"] or not checks["organism_is_homo_sapiens"]:
        return (
            "FAIL",
            "Resolved UniProt entry organism does not match Homo sapiens / taxon 9606 — "
            "clear species mismatch.",
        )

    if not checks["sequence_available"] or not checks["sequence_length_positive"]:
        return "FAIL", "Resolved entry has no usable (positive-length) sequence."

    reasons = []
    if not checks["primary_accession_matches_professor"]:
        reasons.append(
            "Requested accession resolved to a DIFFERENT current primary accession "
            "(possible merge/demerge/redirect) — professor value preserved, not overwritten."
        )
    if not checks["entry_is_reviewed"]:
        reasons.append("Entry is not a reviewed (Swiss-Prot) record — lower curation confidence.")
    if not checks["gene_name_available"]:
        reasons.append("No primary gene name present in the UniProt record.")
    if not checks["protein_name_available"]:
        reasons.append("No recommended/submitted protein name present in the UniProt record.")
    if not checks["naming_overlap_found"]:
        reasons.append(
            "Heuristic overlap check found no match between professor common_name/nr_code "
            "and UniProt gene name, gene synonyms, or protein full name — may still be "
            "correct (heuristic is conservative and can miss legitimate matches); needs "
            "manual confirmation."
        )

    if reasons:
        return "REVIEW_NEEDED", " | ".join(reasons)

    return "PASS", ""
