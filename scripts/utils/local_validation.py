"""Shared local-validation helpers for the professor-supplied NR metadata.

No network access. Pure, deterministic checks over the in-memory record list.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Any

REQUIRED_FIELDS = ["uniprot_id", "nr_code", "common_name", "search_terms", "group"]
NON_BLANK_STRING_FIELDS = ["uniprot_id", "nr_code", "common_name", "group"]

# Standard UniProtKB accession format (6-character form used by all Swiss-Prot
# entries referenced here). See https://www.uniprot.org/help/accession_numbers
UNIPROT_ACCESSION_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$"
    r"|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$"
    r"|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]$"
)

# NR nomenclature format: "NR" + subfamily digit (0-9) + group letter(s) + gene number
NR_CODE_RE = re.compile(r"^NR[0-9][A-Z]{1,2}[0-9]{1,2}$")


def validate_record(record: dict[str, Any]) -> dict[str, Any]:
    """Run all local checks on a single record. Returns a flat dict of booleans/values."""
    result: dict[str, Any] = {}

    missing_fields = [f for f in REQUIRED_FIELDS if f not in record]
    result["has_all_required_fields"] = len(missing_fields) == 0
    result["missing_fields"] = ";".join(missing_fields)

    for field in NON_BLANK_STRING_FIELDS:
        value = record.get(field)
        is_non_blank = isinstance(value, str) and len(value.strip()) > 0
        result[f"{field}_non_blank"] = is_non_blank

    search_terms = record.get("search_terms")
    result["search_terms_is_nonempty_list_of_str"] = (
        isinstance(search_terms, list)
        and len(search_terms) > 0
        and all(isinstance(t, str) and len(t.strip()) > 0 for t in search_terms)
    )

    uniprot_id = record.get("uniprot_id", "")
    result["uniprot_id_format_plausible"] = bool(
        isinstance(uniprot_id, str) and UNIPROT_ACCESSION_RE.match(uniprot_id)
    )

    nr_code = record.get("nr_code", "")
    result["nr_code_format_plausible"] = bool(
        isinstance(nr_code, str) and NR_CODE_RE.match(nr_code)
    )

    result["record_locally_valid"] = all(
        [
            result["has_all_required_fields"],
            result["uniprot_id_non_blank"],
            result["nr_code_non_blank"],
            result["common_name_non_blank"],
            result["group_non_blank"],
            result["search_terms_is_nonempty_list_of_str"],
            result["uniprot_id_format_plausible"],
            result["nr_code_format_plausible"],
        ]
    )

    return result


def find_duplicates(values: list[str]) -> list[str]:
    """Return values that appear more than once, in first-seen order, no repeats."""
    counts = Counter(values)
    seen = []
    for v in values:
        if counts[v] > 1 and v not in seen:
            seen.append(v)
    return seen


def find_exact_duplicate_records(records: list[dict[str, Any]]) -> list[int]:
    """Return record_indexes of records that are exact duplicates of an earlier record."""
    seen_canonical: dict[tuple, int] = {}
    duplicate_indexes: list[int] = []
    for idx, record in enumerate(records):
        canonical = (
            record.get("uniprot_id"),
            record.get("nr_code"),
            record.get("common_name"),
            tuple(sorted(record.get("search_terms", []))),
            record.get("group"),
        )
        if canonical in seen_canonical:
            duplicate_indexes.append(idx)
        else:
            seen_canonical[canonical] = idx
    return duplicate_indexes


def build_local_validation_summary(records: list[dict[str, Any]]) -> dict[str, Any]:
    uniprot_ids = [r.get("uniprot_id") for r in records]
    nr_codes = [r.get("nr_code") for r in records]
    common_names = [r.get("common_name") for r in records]

    per_record = [validate_record(r) for r in records]

    return {
        "total_records": len(records),
        "expected_record_count": 48,
        "record_count_matches_expected": len(records) == 48,
        "unique_uniprot_id_count": len(set(uniprot_ids)),
        "unique_nr_code_count": len(set(nr_codes)),
        "duplicate_uniprot_ids": find_duplicates(uniprot_ids),
        "duplicate_nr_codes": find_duplicates(nr_codes),
        "duplicate_common_names": find_duplicates(common_names),
        "exact_duplicate_record_indexes": find_exact_duplicate_records(records),
        "records_missing_required_fields": [
            i for i, r in enumerate(per_record) if not r["has_all_required_fields"]
        ],
        "records_with_blank_required_string_field": [
            i
            for i, r in enumerate(per_record)
            if not (
                r["uniprot_id_non_blank"]
                and r["nr_code_non_blank"]
                and r["common_name_non_blank"]
                and r["group_non_blank"]
            )
        ],
        "records_with_invalid_search_terms": [
            i for i, r in enumerate(per_record) if not r["search_terms_is_nonempty_list_of_str"]
        ],
        "records_with_implausible_uniprot_format": [
            i for i, r in enumerate(per_record) if not r["uniprot_id_format_plausible"]
        ],
        "records_with_implausible_nr_code_format": [
            i for i, r in enumerate(per_record) if not r["nr_code_format_plausible"]
        ],
        "group_counts": dict(Counter(r.get("group") for r in records)),
        "all_records_locally_valid": all(r["record_locally_valid"] for r in per_record),
    }
