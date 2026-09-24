"""Minimal, reproducible batching client for the official RCSB Data API
GraphQL endpoint (https://data.rcsb.org/graphql).

Batches IDs conservatively (default 150/request), retries transient errors,
and reduces batch size on repeated failure rather than silently dropping
IDs. Every batch attempt (success or failure) is reported back as a
structured result.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

import requests

RCSB_DATA_API_URL = "https://data.rcsb.org/graphql"
USER_AGENT = (
    "BIOL363-NR-Bfactor-DataSupply/0.1 "
    "(Nazarbayev University coursework; reproducible-research script; "
    "contact: a.uzbekbay@gmail.com)"
)
REQUEST_TIMEOUT_SECONDS = 60
MAX_ATTEMPTS = 3
BACKOFF_BASE_SECONDS = 3.0
MIN_SECONDS_BETWEEN_REQUESTS = 0.5
DEFAULT_BATCH_SIZE = 150
MIN_BATCH_SIZE = 25

RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}

POLYMER_ENTITY_QUERY = """
query ($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    rcsb_polymer_entity_container_identifiers {
      entry_id
      entity_id
      asym_ids
      auth_asym_ids
      reference_sequence_identifiers {
        database_accession
        database_name
      }
    }
    entity_poly {
      rcsb_sample_sequence_length
      type
    }
    rcsb_polymer_entity {
      pdbx_description
    }
    rcsb_entity_source_organism {
      scientific_name
      ncbi_taxonomy_id
    }
    rcsb_polymer_entity_align {
      reference_database_name
      reference_database_accession
      aligned_regions {
        entity_beg_seq_id
        ref_beg_seq_id
        length
      }
    }
  }
}
"""

POLYMER_ENTITY_INSTANCES_QUERY = """
query ($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    polymer_entity_instances {
      rcsb_id
      rcsb_polymer_entity_instance_container_identifiers {
        entry_id
        entity_id
        asym_id
        auth_asym_id
      }
    }
  }
}
"""

HOST_ORGANISM_SUPPLEMENT_QUERY = """
query ($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    rcsb_entity_host_organism {
      scientific_name
      ncbi_taxonomy_id
    }
  }
}
"""

ENTRY_CONTAINER_IDENTIFIERS_QUERY = """
query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_entry_container_identifiers {
      polymer_entity_ids
      non_polymer_entity_ids
      water_entity_ids
    }
  }
}
"""

NONPOLYMER_ENTITY_QUERY = """
query ($ids: [String!]!) {
  nonpolymer_entities(entity_ids: $ids) {
    rcsb_id
    rcsb_nonpolymer_entity_container_identifiers {
      entry_id
      entity_id
      asym_ids
      auth_asym_ids
      nonpolymer_comp_id
    }
    pdbx_entity_nonpoly {
      comp_id
      name
      entity_id
    }
    rcsb_nonpolymer_entity_keywords {
      text
    }
    rcsb_nonpolymer_entity_annotation {
      annotation_id
      type
      name
      description
      provenance_source
    }
    nonpolymer_comp {
      chem_comp {
        id
        name
        type
        formula
        formula_weight
        three_letter_code
      }
      rcsb_chem_comp_info {
        atom_count
        atom_count_heavy
      }
    }
    nonpolymer_entity_instances {
      rcsb_id
      rcsb_nonpolymer_entity_instance_container_identifiers {
        asym_id
        auth_asym_id
        auth_seq_id
        comp_id
        entry_id
      }
    }
  }
}
"""

ENTRY_QUERY = """
query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct {
      title
    }
    rcsb_entry_info {
      experimental_method
      resolution_combined
      structure_determination_methodology
      polymer_entity_count
      polymer_entity_count_protein
      nonpolymer_entity_count
      deposited_polymer_entity_instance_count
    }
    exptl {
      method
    }
    refine {
      ls_R_factor_R_work
      ls_R_factor_R_free
    }
    rcsb_accession_info {
      initial_release_date
      deposit_date
      revision_date
    }
    symmetry {
      space_group_name_H_M
    }
    rcsb_primary_citation {
      title
      pdbx_database_id_DOI
      year
    }
  }
}
"""


@dataclass
class RCSBBatchResult:
    batch_type: str
    batch_number: int
    ids_requested: list[str]
    request_url: str
    http_status: int | None
    retrieved_at_utc: str
    success: bool
    error_message: str
    raw_json: dict[str, Any] | None = None
    attempts_made: int = 0


class RateLimiter:
    def __init__(self, min_interval_seconds: float) -> None:
        self.min_interval_seconds = min_interval_seconds
        self._last_request_time: float | None = None

    def wait(self) -> None:
        if self._last_request_time is not None:
            elapsed = time.monotonic() - self._last_request_time
            remaining = self.min_interval_seconds - elapsed
            if remaining > 0:
                time.sleep(remaining)
        self._last_request_time = time.monotonic()


_rate_limiter = RateLimiter(MIN_SECONDS_BETWEEN_REQUESTS)


def _post_graphql(query: str, ids: list[str]) -> tuple[int | None, dict | None, str, int]:
    """Single GraphQL POST with bounded retry/backoff. Returns
    (http_status, parsed_json_or_None, error_message, attempts_made)."""
    headers = {"Content-Type": "application/json", "User-Agent": USER_AGENT}
    payload = {"query": query, "variables": {"ids": ids}}

    last_status: int | None = None
    last_error = ""

    for attempt in range(1, MAX_ATTEMPTS + 1):
        _rate_limiter.wait()
        try:
            response = requests.post(
                RCSB_DATA_API_URL, json=payload, headers=headers, timeout=REQUEST_TIMEOUT_SECONDS
            )
        except requests.exceptions.Timeout:
            last_error = f"Timeout after {REQUEST_TIMEOUT_SECONDS}s (attempt {attempt}/{MAX_ATTEMPTS})"
            last_status = None
        except requests.exceptions.RequestException as exc:
            last_error = f"Network error on attempt {attempt}/{MAX_ATTEMPTS}: {exc}"
            last_status = None
        else:
            last_status = response.status_code
            if response.status_code == 200:
                try:
                    parsed = response.json()
                except ValueError as exc:
                    return response.status_code, None, f"HTTP 200 but invalid JSON: {exc}", attempt
                if "errors" in parsed and parsed.get("errors"):
                    return (
                        response.status_code,
                        parsed,
                        f"GraphQL response contained errors: {parsed['errors']}",
                        attempt,
                    )
                return response.status_code, parsed, "", attempt

            if response.status_code not in RETRYABLE_STATUS_CODES:
                return response.status_code, None, f"HTTP {response.status_code}: {response.text[:500]}", attempt

            last_error = f"HTTP {response.status_code} (retryable) on attempt {attempt}/{MAX_ATTEMPTS}"

        if attempt < MAX_ATTEMPTS:
            time.sleep(BACKOFF_BASE_SECONDS * (2 ** (attempt - 1)))

    return last_status, None, f"Exhausted {MAX_ATTEMPTS} attempts; last error: {last_error}", MAX_ATTEMPTS


def fetch_batches(
    batch_type: str,
    all_ids: list[str],
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> list[RCSBBatchResult]:
    """Fetch metadata for `all_ids` (polymer entity IDs or entry IDs) in
    conservative batches. On a batch failure, retries once at half the
    batch size (down to MIN_BATCH_SIZE) before giving up and recording the
    failure — never silently drops IDs from the audit trail.
    """
    if batch_type == "polymer_entity":
        query = POLYMER_ENTITY_QUERY
    elif batch_type == "entry":
        query = ENTRY_QUERY
    elif batch_type == "host_organism_supplement":
        query = HOST_ORGANISM_SUPPLEMENT_QUERY
    elif batch_type == "polymer_entity_instances":
        query = POLYMER_ENTITY_INSTANCES_QUERY
    elif batch_type == "entry_container_identifiers":
        query = ENTRY_CONTAINER_IDENTIFIERS_QUERY
    elif batch_type == "nonpolymer_entity":
        query = NONPOLYMER_ENTITY_QUERY
    else:
        raise ValueError(f"Unknown batch_type: {batch_type}")

    results: list[RCSBBatchResult] = []
    batch_number = 0
    idx = 0
    current_batch_size = batch_size

    while idx < len(all_ids):
        batch_number += 1
        chunk = all_ids[idx : idx + current_batch_size]
        retrieved_at_utc = datetime.now(timezone.utc).isoformat()
        status, parsed, error, attempts = _post_graphql(query, chunk)

        if parsed is not None and not error:
            results.append(
                RCSBBatchResult(
                    batch_type=batch_type,
                    batch_number=batch_number,
                    ids_requested=chunk,
                    request_url=RCSB_DATA_API_URL,
                    http_status=status,
                    retrieved_at_utc=retrieved_at_utc,
                    success=True,
                    error_message="",
                    raw_json=parsed,
                    attempts_made=attempts,
                )
            )
            idx += current_batch_size
            continue

        # Batch failed (timeout, non-retryable error, or GraphQL errors).
        # If we still have room to shrink, retry this same chunk smaller
        # rather than silently skipping it.
        if current_batch_size > MIN_BATCH_SIZE:
            current_batch_size = max(MIN_BATCH_SIZE, current_batch_size // 2)
            results.append(
                RCSBBatchResult(
                    batch_type=batch_type,
                    batch_number=batch_number,
                    ids_requested=chunk,
                    request_url=RCSB_DATA_API_URL,
                    http_status=status,
                    retrieved_at_utc=retrieved_at_utc,
                    success=False,
                    error_message=f"{error} — reducing batch size to {current_batch_size} and retrying this chunk",
                    raw_json=parsed,
                    attempts_made=attempts,
                )
            )
            continue  # retry same idx with smaller current_batch_size

        # Already at minimum batch size and still failing: record the
        # permanent failure for this chunk and move on (do not spin forever).
        results.append(
            RCSBBatchResult(
                batch_type=batch_type,
                batch_number=batch_number,
                ids_requested=chunk,
                request_url=RCSB_DATA_API_URL,
                http_status=status,
                retrieved_at_utc=retrieved_at_utc,
                success=False,
                error_message=error,
                raw_json=parsed,
                attempts_made=attempts,
            )
        )
        idx += current_batch_size

    return results
