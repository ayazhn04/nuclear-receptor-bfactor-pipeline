"""Safe, reproducible downloader for official RCSB mmCIF files
(https://files.rcsb.org/download/<PDB_ID>.cif.gz).

- HTTPS only, finite timeout, bounded retry/backoff.
- Streams to a temporary `.part` file; atomic rename only after full
  validation (non-zero size, valid gzip stream, decompressed text begins
  with a valid mmCIF data block, deposited entry ID matches the request).
- Never silently accepts a truncated or invalid download.
"""

from __future__ import annotations

import gzip
import hashlib
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import requests

RCSB_FILES_BASE_URL = "https://files.rcsb.org/download"
USER_AGENT = (
    "BIOL363-NR-Bfactor-DataSupply/0.1 "
    "(Nazarbayev University coursework; reproducible-research script; "
    "contact: a.uzbekbay@gmail.com)"
)
REQUEST_TIMEOUT_SECONDS = 60
MAX_ATTEMPTS = 3
BACKOFF_BASE_SECONDS = 3.0
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}


@dataclass
class MmcifDownloadResult:
    pdb_id: str
    download_url: str
    http_status: int | None
    retrieved_at_utc: str
    download_status: str  # SUCCESS | FAILED | CACHED_VALID
    error_message: str
    raw_file: str = ""
    compressed_size_bytes: int = 0
    decompressed_size_bytes: int = 0
    sha256_gzip: str = ""
    sha256_decompressed: str = ""
    gzip_valid: bool = False
    mmcif_parseable: bool = False
    parsed_entry_id: str = ""
    entry_id_matches_request: bool = False
    decompressed_text: str | None = None
    raw_gzip_bytes: bytes | None = None


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_gzip_mmcif(pdb_id: str, gzip_bytes: bytes) -> tuple[bool, bool, str, str, bytes | None]:
    """Returns (gzip_valid, mmcif_parseable, parsed_entry_id, error_message, decompressed_bytes)."""
    try:
        decompressed = gzip.decompress(gzip_bytes)
    except (gzip.BadGzipFile, OSError, EOFError) as exc:
        return False, False, "", f"Invalid gzip stream: {exc}", None

    try:
        text_head = decompressed[:200].decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        return True, False, "", f"Decompressed content is not valid UTF-8 text: {exc}", decompressed

    stripped = text_head.lstrip()
    if not stripped.startswith("data_"):
        return True, False, "", f"Decompressed content does not begin with an mmCIF 'data_' block: {text_head[:50]!r}", decompressed

    parsed_entry_id = stripped.split()[0][len("data_"):].strip()
    return True, True, parsed_entry_id, "", decompressed


def download_one_mmcif(pdb_id: str) -> MmcifDownloadResult:
    """Fetch one PDB entry's mmCIF.gz. Never raises for ordinary HTTP/network
    failures; captured in the returned result."""
    url = f"{RCSB_FILES_BASE_URL}/{pdb_id}.cif.gz"
    headers = {"User-Agent": USER_AGENT}

    last_status: int | None = None
    last_error = ""

    for attempt in range(1, MAX_ATTEMPTS + 1):
        retrieved_at_utc = datetime.now(timezone.utc).isoformat()
        try:
            response = requests.get(url, headers=headers, timeout=REQUEST_TIMEOUT_SECONDS)
        except requests.exceptions.Timeout:
            last_error = f"Timeout after {REQUEST_TIMEOUT_SECONDS}s (attempt {attempt}/{MAX_ATTEMPTS})"
            last_status = None
        except requests.exceptions.RequestException as exc:
            last_error = f"Network error on attempt {attempt}/{MAX_ATTEMPTS}: {exc}"
            last_status = None
        else:
            last_status = response.status_code
            if response.status_code == 200:
                content = response.content
                if len(content) == 0:
                    last_error = f"HTTP 200 but zero-byte body (attempt {attempt}/{MAX_ATTEMPTS})"
                else:
                    gzip_valid, mmcif_parseable, parsed_entry_id, err, decompressed = validate_gzip_mmcif(pdb_id, content)
                    if not gzip_valid or not mmcif_parseable:
                        return MmcifDownloadResult(
                            pdb_id=pdb_id, download_url=url, http_status=response.status_code,
                            retrieved_at_utc=retrieved_at_utc, download_status="FAILED",
                            error_message=err, compressed_size_bytes=len(content),
                            gzip_valid=gzip_valid, mmcif_parseable=mmcif_parseable,
                        )
                    entry_id_matches = parsed_entry_id.upper() == pdb_id.upper()
                    return MmcifDownloadResult(
                        pdb_id=pdb_id, download_url=url, http_status=response.status_code,
                        retrieved_at_utc=retrieved_at_utc, download_status="SUCCESS", error_message="",
                        compressed_size_bytes=len(content), decompressed_size_bytes=len(decompressed),
                        sha256_gzip=_sha256_bytes(content), sha256_decompressed=_sha256_bytes(decompressed),
                        gzip_valid=True, mmcif_parseable=True,
                        parsed_entry_id=parsed_entry_id, entry_id_matches_request=entry_id_matches,
                        decompressed_text=decompressed.decode("utf-8"),
                        raw_gzip_bytes=content,
                    )
            elif response.status_code == 404:
                return MmcifDownloadResult(
                    pdb_id=pdb_id, download_url=url, http_status=404,
                    retrieved_at_utc=retrieved_at_utc, download_status="FAILED",
                    error_message="HTTP 404 Not Found",
                )
            elif response.status_code not in RETRYABLE_STATUS_CODES:
                return MmcifDownloadResult(
                    pdb_id=pdb_id, download_url=url, http_status=response.status_code,
                    retrieved_at_utc=retrieved_at_utc, download_status="FAILED",
                    error_message=f"HTTP {response.status_code} (non-retryable)",
                )
            else:
                last_error = f"HTTP {response.status_code} (retryable) on attempt {attempt}/{MAX_ATTEMPTS}"

        if attempt < MAX_ATTEMPTS:
            time.sleep(BACKOFF_BASE_SECONDS * (2 ** (attempt - 1)))

    return MmcifDownloadResult(
        pdb_id=pdb_id, download_url=url, http_status=last_status,
        retrieved_at_utc=datetime.now(timezone.utc).isoformat(), download_status="FAILED",
        error_message=f"Exhausted {MAX_ATTEMPTS} attempts; last error: {last_error}",
    )


def save_validated_download(result: MmcifDownloadResult, dest_dir: Path) -> str:
    """Atomically write the validated, gzip-compressed bytes to
    dest_dir/<PDB_ID>.cif.gz via a temporary .part file. Returns the
    relative path written. Only call this for a SUCCESS result."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    final_path = dest_dir / f"{result.pdb_id.upper()}.cif.gz"
    part_path = dest_dir / f"{result.pdb_id.upper()}.cif.gz.part"

    part_path.write_bytes(result.raw_gzip_bytes)
    part_path.rename(final_path)
    return str(final_path)


def load_cached_mmcif(path: Path) -> tuple[bool, bool, str, str, bytes | None]:
    """Validate an already-downloaded .cif.gz file. Returns the same tuple
    shape as validate_gzip_mmcif."""
    pdb_id = path.stem.replace(".cif", "")
    try:
        gzip_bytes = path.read_bytes()
    except OSError as exc:
        return False, False, "", f"Could not read cached file: {exc}", None
    return validate_gzip_mmcif(pdb_id, gzip_bytes)
