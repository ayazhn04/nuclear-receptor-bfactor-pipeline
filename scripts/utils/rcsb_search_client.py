"""Minimal, reproducible client for the official RCSB Search API
(https://search.rcsb.org/rcsbsearch/v2/query).

Discovers current-holdings, experimental polymer entities whose
reference-sequence mapping includes a given UniProt accession.

Design mirrors scripts/utils/uniprot_client.py: HTTPS only, finite timeout,
bounded retry/backoff for transient errors, polite pacing, and every
attempt (success or failure) is reported back as a structured result —
nothing is silently dropped.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

import requests

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
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


def build_uniprot_polymer_entity_query(accession: str) -> dict[str, Any]:
    """Build the RCSB Search API request payload for discovering current
    EXPERIMENTAL polymer entities whose reference-sequence mapping includes
    `accession` (a validated UniProt accession).

    Criteria used (confirmed against the live schema during the Stage 2
    pilot, see reports/tables/rcsb_discovery_report.md "Methods"):
      - rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers
          .database_accession == accession
      - ...database_name == "UniProt"  (defends against a same-string
        accession collision in a different reference database)
      - rcsb_entry_info.structure_determination_methodology == "experimental"
        (excludes computed structure models, e.g. AlphaFold, per Stage 2
        scope rules — NOT an experimental-method restriction; X-ray, cryo-EM,
        NMR, and other experimental methods are all retained here)
    """
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": (
                            "rcsb_polymer_entity_container_identifiers."
                            "reference_sequence_identifiers.database_accession"
                        ),
                        "operator": "exact_match",
                        "value": accession,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": (
                            "rcsb_polymer_entity_container_identifiers."
                            "reference_sequence_identifiers.database_name"
                        ),
                        "operator": "exact_match",
                        "value": "UniProt",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.structure_determination_methodology",
                        "operator": "exact_match",
                        "value": "experimental",
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {"return_all_hits": True},
    }


@dataclass
class RCSBSearchResult:
    requested_accession: str
    request_url: str
    request_payload: dict[str, Any]
    http_status: int | None
    retrieved_at_utc: str
    success: bool
    error_message: str
    raw_json: dict[str, Any] | None = None
    identifiers: list[str] = field(default_factory=list)
    total_count: int | None = None
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


def search_polymer_entities_by_uniprot(accession: str) -> RCSBSearchResult:
    """Query the RCSB Search API for current experimental polymer entities
    mapped to `accession`. A RCSB Search API "no hits" response is HTTP 200
    with an empty result_set (or, for some server versions, HTTP 204) — both
    are treated as a SUCCESSFUL zero-result query, distinct from a failed
    request. Never raises for ordinary HTTP/network failures.
    """
    payload = build_uniprot_polymer_entity_query(accession)
    headers = {"Content-Type": "application/json", "User-Agent": USER_AGENT}

    last_status: int | None = None
    last_error = ""

    for attempt in range(1, MAX_ATTEMPTS + 1):
        _rate_limiter.wait()
        retrieved_at_utc = datetime.now(timezone.utc).isoformat()
        try:
            response = requests.post(
                RCSB_SEARCH_URL, json=payload, headers=headers, timeout=REQUEST_TIMEOUT_SECONDS
            )
        except requests.exceptions.Timeout:
            last_error = f"Timeout after {REQUEST_TIMEOUT_SECONDS}s (attempt {attempt}/{MAX_ATTEMPTS})"
            last_status = None
        except requests.exceptions.RequestException as exc:
            last_error = f"Network error on attempt {attempt}/{MAX_ATTEMPTS}: {exc}"
            last_status = None
        else:
            last_status = response.status_code

            if response.status_code == 204:
                # RCSB returns 204 No Content for a syntactically valid
                # query with zero matches (server-version dependent).
                return RCSBSearchResult(
                    requested_accession=accession,
                    request_url=RCSB_SEARCH_URL,
                    request_payload=payload,
                    http_status=response.status_code,
                    retrieved_at_utc=retrieved_at_utc,
                    success=True,
                    error_message="",
                    raw_json={"query_id": None, "result_type": "polymer_entity", "total_count": 0, "result_set": []},
                    identifiers=[],
                    total_count=0,
                    attempts_made=attempt,
                )

            if response.status_code == 200:
                content_type = response.headers.get("Content-Type", "")
                if "json" not in content_type.lower():
                    return RCSBSearchResult(
                        requested_accession=accession,
                        request_url=RCSB_SEARCH_URL,
                        request_payload=payload,
                        http_status=response.status_code,
                        retrieved_at_utc=retrieved_at_utc,
                        success=False,
                        error_message=f"HTTP 200 but unexpected Content-Type '{content_type}'",
                        attempts_made=attempt,
                    )
                try:
                    parsed = response.json()
                except ValueError as exc:
                    return RCSBSearchResult(
                        requested_accession=accession,
                        request_url=RCSB_SEARCH_URL,
                        request_payload=payload,
                        http_status=response.status_code,
                        retrieved_at_utc=retrieved_at_utc,
                        success=False,
                        error_message=f"HTTP 200 but response body was not valid JSON: {exc}",
                        attempts_made=attempt,
                    )
                identifiers = [r["identifier"] for r in parsed.get("result_set", [])]
                return RCSBSearchResult(
                    requested_accession=accession,
                    request_url=RCSB_SEARCH_URL,
                    request_payload=payload,
                    http_status=response.status_code,
                    retrieved_at_utc=retrieved_at_utc,
                    success=True,
                    error_message="",
                    raw_json=parsed,
                    identifiers=identifiers,
                    total_count=parsed.get("total_count", len(identifiers)),
                    attempts_made=attempt,
                )

            if response.status_code not in RETRYABLE_STATUS_CODES:
                body_snippet = response.text[:500]
                return RCSBSearchResult(
                    requested_accession=accession,
                    request_url=RCSB_SEARCH_URL,
                    request_payload=payload,
                    http_status=response.status_code,
                    retrieved_at_utc=retrieved_at_utc,
                    success=False,
                    error_message=f"HTTP {response.status_code} (non-retryable): {body_snippet}",
                    attempts_made=attempt,
                )

            last_error = f"HTTP {response.status_code} (retryable) on attempt {attempt}/{MAX_ATTEMPTS}"

        if attempt < MAX_ATTEMPTS:
            time.sleep(BACKOFF_BASE_SECONDS * (2 ** (attempt - 1)))

    return RCSBSearchResult(
        requested_accession=accession,
        request_url=RCSB_SEARCH_URL,
        request_payload=payload,
        http_status=last_status,
        retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
        success=False,
        error_message=f"Exhausted {MAX_ATTEMPTS} attempts; last error: {last_error}",
        attempts_made=MAX_ATTEMPTS,
    )
