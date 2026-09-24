"""Client for the official RCSB Sequence Coordinates GraphQL service
(https://sequence-coordinates.rcsb.org/graphql).

Schema confirmed via live introspection on 2026-09-23 (see
reports/tables/rcsb_discovery_report.md "Methods" for the recorded
evidence):

    alignments(from: SequenceReference!, queryId: String!, to: SequenceReference!, range: [...]): SequenceAlignments
    SequenceReference enum includes UNIPROT and PDB_ENTITY (confirmed).
    SequenceAlignments { query_sequence, target_alignments: [TargetAlignments] }
    TargetAlignments { target_id, orientation, coverage: Coverage, aligned_regions: [AlignedRegions] }
    Coverage { query_coverage, query_length, target_coverage, target_length }
    AlignedRegions { query_begin, query_end, target_begin, target_end }

Confirmed semantics (from = UNIPROT, to = PDB_ENTITY), via a manually
understandable example (VDR/P11473, target 1KB2_3):
    query_length = 427  == the full UniProt (P11473) sequence length
    target_length = 110 == the PDB polymer entity's own sequence length
    query_coverage = 0.2576 == 110/427 (aligned length / UniProt length)
    target_coverage = 1.0   == 110/110 (aligned length / entity length)
  => query_coverage is therefore the REFERENCE (UniProt) sequence coverage,
     and target_coverage is the ENTITY (PDB polymer entity) sequence
     coverage — confirmed by the numbers, not assumed from field names alone.

`target_id` values include BOTH experimental PDB entities (e.g. "1KB2_3")
and computed structure models (e.g. "AF_AFP11473F1_1" for AlphaFold) when
querying to: PDB_ENTITY — confirmed empirically. Callers must filter to
4-character PDB ID targets to match Stage 2's experimental-only scope.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

import requests

SEQUENCE_COORDINATES_URL = "https://sequence-coordinates.rcsb.org/graphql"
USER_AGENT = (
    "BIOL363-NR-Bfactor-DataSupply/0.1 "
    "(Nazarbayev University coursework; reproducible-research script; "
    "contact: a.uzbekbay@gmail.com)"
)
REQUEST_TIMEOUT_SECONDS = 30
MAX_ATTEMPTS = 3
BACKOFF_BASE_SECONDS = 2.0
MIN_SECONDS_BETWEEN_REQUESTS = 0.5

RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}

ALIGNMENTS_QUERY = """
query ($id: String!) {
  alignments(from: UNIPROT, to: PDB_ENTITY, queryId: $id) {
    query_sequence
    target_alignments {
      target_id
      orientation
      coverage {
        query_coverage
        query_length
        target_coverage
        target_length
      }
      aligned_regions {
        query_begin
        query_end
        target_begin
        target_end
      }
    }
  }
}
"""


@dataclass
class SequenceCoordinatesResult:
    requested_accession: str
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


def fetch_alignments_for_uniprot(accession: str) -> SequenceCoordinatesResult:
    """Fetch UNIPROT -> PDB_ENTITY alignments for one UniProt accession.
    Never raises for ordinary HTTP/network failures — captured in the
    returned result."""
    headers = {"Content-Type": "application/json", "User-Agent": USER_AGENT}
    payload = {"query": ALIGNMENTS_QUERY, "variables": {"id": accession}}

    last_status: int | None = None
    last_error = ""

    for attempt in range(1, MAX_ATTEMPTS + 1):
        _rate_limiter.wait()
        retrieved_at_utc = datetime.now(timezone.utc).isoformat()
        try:
            response = requests.post(
                SEQUENCE_COORDINATES_URL, json=payload, headers=headers, timeout=REQUEST_TIMEOUT_SECONDS
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
                    return SequenceCoordinatesResult(
                        requested_accession=accession, request_url=SEQUENCE_COORDINATES_URL,
                        http_status=response.status_code, retrieved_at_utc=retrieved_at_utc,
                        success=False, error_message=f"HTTP 200 but invalid JSON: {exc}", attempts_made=attempt,
                    )
                if "errors" in parsed and parsed.get("errors"):
                    return SequenceCoordinatesResult(
                        requested_accession=accession, request_url=SEQUENCE_COORDINATES_URL,
                        http_status=response.status_code, retrieved_at_utc=retrieved_at_utc,
                        success=False, error_message=f"GraphQL errors: {parsed['errors']}",
                        raw_json=parsed, attempts_made=attempt,
                    )
                if parsed.get("data", {}).get("alignments") is None:
                    # Syntactically valid response but no alignments object at
                    # all for this accession — a genuine zero-target result,
                    # not a failure. Represented with an empty target list.
                    parsed = {"data": {"alignments": {"query_sequence": None, "target_alignments": []}}}
                return SequenceCoordinatesResult(
                    requested_accession=accession, request_url=SEQUENCE_COORDINATES_URL,
                    http_status=response.status_code, retrieved_at_utc=retrieved_at_utc,
                    success=True, error_message="", raw_json=parsed, attempts_made=attempt,
                )

            if response.status_code not in RETRYABLE_STATUS_CODES:
                return SequenceCoordinatesResult(
                    requested_accession=accession, request_url=SEQUENCE_COORDINATES_URL,
                    http_status=response.status_code, retrieved_at_utc=retrieved_at_utc,
                    success=False, error_message=f"HTTP {response.status_code} (non-retryable): {response.text[:500]}",
                    attempts_made=attempt,
                )
            last_error = f"HTTP {response.status_code} (retryable) on attempt {attempt}/{MAX_ATTEMPTS}"

        if attempt < MAX_ATTEMPTS:
            time.sleep(BACKOFF_BASE_SECONDS * (2 ** (attempt - 1)))

    return SequenceCoordinatesResult(
        requested_accession=accession, request_url=SEQUENCE_COORDINATES_URL,
        http_status=last_status, retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
        success=False, error_message=f"Exhausted {MAX_ATTEMPTS} attempts; last error: {last_error}",
        attempts_made=MAX_ATTEMPTS,
    )
