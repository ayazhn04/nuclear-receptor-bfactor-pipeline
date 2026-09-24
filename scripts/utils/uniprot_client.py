"""Minimal, reproducible client for the official UniProt REST API.

Design goals (per project data-safety requirements):
- HTTPS only, official rest.uniprot.org endpoint.
- Finite timeout, bounded retries with backoff for transient errors only.
- Polite rate limiting between requests.
- Never silently drop a failed request — every attempt is reported back
  to the caller as a structured result, success or failure.
- Never fabricate a response, checksum, or field value.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

import requests

UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"
USER_AGENT = (
    "BIOL363-NR-Bfactor-DataSupply/0.1 "
    "(Nazarbayev University coursework; reproducible-research script; "
    "contact: a.uzbekbay@gmail.com)"
)
REQUEST_TIMEOUT_SECONDS = 20
MAX_ATTEMPTS = 3
BACKOFF_BASE_SECONDS = 2.0
MIN_SECONDS_BETWEEN_REQUESTS = 1.0

# Status codes worth retrying: transient server-side / rate-limit conditions.
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}


@dataclass
class UniProtFetchResult:
    requested_accession: str
    request_url: str
    http_status: int | None
    retrieved_at_utc: str
    success: bool
    error_message: str
    raw_json: dict[str, Any] | None = None
    redirect_history: list[str] = field(default_factory=list)
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


def fetch_uniprot_entry(accession: str) -> UniProtFetchResult:
    """Fetch a single UniProtKB entry as JSON, with retry/backoff.

    Never raises for ordinary HTTP/network failures — those are captured in
    the returned UniProtFetchResult. Only truly unexpected local errors
    (e.g. programming bugs) would propagate.
    """
    url = f"{UNIPROT_API_BASE}/{accession}.json"
    headers = {"Accept": "application/json", "User-Agent": USER_AGENT}

    last_status: int | None = None
    last_error = ""
    redirect_history: list[str] = []

    for attempt in range(1, MAX_ATTEMPTS + 1):
        _rate_limiter.wait()
        retrieved_at_utc = datetime.now(timezone.utc).isoformat()
        try:
            response = requests.get(
                url, headers=headers, timeout=REQUEST_TIMEOUT_SECONDS, allow_redirects=True
            )
        except requests.exceptions.Timeout:
            last_error = f"Timeout after {REQUEST_TIMEOUT_SECONDS}s (attempt {attempt}/{MAX_ATTEMPTS})"
            last_status = None
        except requests.exceptions.RequestException as exc:
            last_error = f"Network error on attempt {attempt}/{MAX_ATTEMPTS}: {exc}"
            last_status = None
        else:
            last_status = response.status_code
            redirect_history = [r.url for r in response.history]

            if response.status_code == 200:
                content_type = response.headers.get("Content-Type", "")
                if "json" not in content_type.lower():
                    return UniProtFetchResult(
                        requested_accession=accession,
                        request_url=url,
                        http_status=response.status_code,
                        retrieved_at_utc=retrieved_at_utc,
                        success=False,
                        error_message=(
                            f"HTTP 200 but unexpected Content-Type '{content_type}' "
                            f"(not JSON) — refusing to parse"
                        ),
                        redirect_history=redirect_history,
                        attempts_made=attempt,
                    )
                try:
                    parsed = response.json()
                except json.JSONDecodeError as exc:
                    return UniProtFetchResult(
                        requested_accession=accession,
                        request_url=url,
                        http_status=response.status_code,
                        retrieved_at_utc=retrieved_at_utc,
                        success=False,
                        error_message=f"HTTP 200 but response body was not valid JSON: {exc}",
                        redirect_history=redirect_history,
                        attempts_made=attempt,
                    )
                return UniProtFetchResult(
                    requested_accession=accession,
                    request_url=url,
                    http_status=response.status_code,
                    retrieved_at_utc=retrieved_at_utc,
                    success=True,
                    error_message="",
                    raw_json=parsed,
                    redirect_history=redirect_history,
                    attempts_made=attempt,
                )

            if response.status_code == 404:
                return UniProtFetchResult(
                    requested_accession=accession,
                    request_url=url,
                    http_status=response.status_code,
                    retrieved_at_utc=retrieved_at_utc,
                    success=False,
                    error_message=(
                        "HTTP 404 Not Found — accession not resolvable via current "
                        "UniProt REST API (may be obsolete, deleted, or never valid)"
                    ),
                    redirect_history=redirect_history,
                    attempts_made=attempt,
                )

            if response.status_code not in RETRYABLE_STATUS_CODES:
                # Non-retryable, non-200, non-404 status (e.g. 400 Bad Request,
                # 300 Multiple Choices for demerged accessions).
                body_snippet = response.text[:500]
                return UniProtFetchResult(
                    requested_accession=accession,
                    request_url=url,
                    http_status=response.status_code,
                    retrieved_at_utc=retrieved_at_utc,
                    success=False,
                    error_message=(
                        f"HTTP {response.status_code} (non-retryable): {body_snippet}"
                    ),
                    redirect_history=redirect_history,
                    attempts_made=attempt,
                )

            last_error = f"HTTP {response.status_code} (retryable) on attempt {attempt}/{MAX_ATTEMPTS}"

        if attempt < MAX_ATTEMPTS:
            backoff = BACKOFF_BASE_SECONDS * (2 ** (attempt - 1))
            time.sleep(backoff)

    return UniProtFetchResult(
        requested_accession=accession,
        request_url=url,
        http_status=last_status,
        retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
        success=False,
        error_message=f"Exhausted {MAX_ATTEMPTS} attempts; last error: {last_error}",
        redirect_history=redirect_history,
        attempts_made=MAX_ATTEMPTS,
    )
