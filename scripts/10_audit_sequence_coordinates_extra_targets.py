"""
Stage 2.3: correct the Stage 2.2 interpretation of Sequence Coordinates
targets that fall outside the Stage 2 Search API candidate set.

Stage 2.2 labeled these "unexpected targets" and speculated they reflected
an indexing-timing lag between RCSB microservices. That speculation was
NOT verified against actual RCSB metadata and has been shown to be wrong:
programmatic inspection (this script) confirms both underlying PDB entries
(8ZZR, 8ZZZ) are INTEGRATIVE structures (rcsb_entry_info.structure_
determination_methodology == "integrative", not "experimental"), deposited
in 2018/2019 — not recently. They were correctly excluded from Stage 2
because Stage 2's Search API query explicitly filters on
structure_determination_methodology == "experimental"; the Sequence
Coordinates API applies no such filter, so it surfaces them anyway.

Reads:
    data/manifests/rcsb_sequence_coordinate_alignments.csv
    data/manifests/rcsb_discovery_by_receptor.csv
    data/raw/api/rcsb/data/entries_batch_*.json   (reused where the PDB ID is already cached)

Writes:
    data/raw/api/rcsb/data/extra_targets_entries.json   (gitignored; minimal targeted fetch)
    data/manifests/rcsb_sequence_coordinates_extra_targets.csv
    data/manifests/rcsb_sequence_coordinates_manifest.csv   (terminology correction, in place)

This script performs NO coordinate downloads and makes at most one small,
targeted Data API request (only for PDB IDs not already covered by the
Stage 2 raw batches).
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.rcsb_data_client import RCSB_DATA_API_URL, fetch_batches  # noqa: E402
from scripts.utils.rcsb_extract import extract_entry_fields  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "data"

ALIGNMENTS_CSV = MANIFESTS_DIR / "rcsb_sequence_coordinate_alignments.csv"
DISCOVERY_BY_RECEPTOR_CSV = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
SEQ_COORD_MANIFEST_CSV = MANIFESTS_DIR / "rcsb_sequence_coordinates_manifest.csv"
EXTRA_TARGETS_CSV = MANIFESTS_DIR / "rcsb_sequence_coordinates_extra_targets.csv"
EXTRA_TARGETS_RAW_PATH = RAW_DATA_DIR / "extra_targets_entries.json"

EXPERIMENTAL_METHODOLOGY = "experimental"


def load_cached_entry_records() -> dict[str, dict]:
    records: dict[str, dict] = {}
    for path in sorted(RAW_DATA_DIR.glob("entries_batch_*.json")):
        wrapper = json.loads(path.read_text())
        if not wrapper.get("success"):
            continue
        entries = ((wrapper.get("response") or {}).get("data") or {}).get("entries") or []
        for e in entries:
            if e and e.get("rcsb_id"):
                records[e["rcsb_id"]] = e
    return records


def fetch_extra_entry_metadata(pdb_ids: list[str], refresh: bool) -> dict[str, dict]:
    """Fetch entry metadata for `pdb_ids` not already covered by cached
    Stage 2 batches, via one small targeted Data API request. Cached in a
    single wrapper file; reused unless --refresh."""
    if not refresh and EXTRA_TARGETS_RAW_PATH.exists():
        wrapper = json.loads(EXTRA_TARGETS_RAW_PATH.read_text())
        cached_ids = set(wrapper.get("ids_requested", []))
        if cached_ids >= set(pdb_ids):
            entries = ((wrapper.get("response") or {}).get("data") or {}).get("entries") or []
            print(f"  extra_targets_entries: {len(pdb_ids)} IDs already covered by cache — reusing.")
            return {e["rcsb_id"]: e for e in entries if e and e.get("rcsb_id")}

    print(f"  extra_targets_entries: fetching entry metadata for {len(pdb_ids)} PDB ID(s) via minimal targeted request...")
    results = fetch_batches("entry", sorted(pdb_ids), batch_size=max(len(pdb_ids), 1))
    RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)

    records: dict[str, dict] = {}
    all_ids_requested: list[str] = []
    combined_entries: list[dict] = []
    success = True
    for r in results:
        all_ids_requested.extend(r.ids_requested)
        success = success and r.success
        if r.success and r.raw_json:
            entries = ((r.raw_json.get("data") or {}).get("entries")) or []
            combined_entries.extend(entries)
            for e in entries:
                if e and e.get("rcsb_id"):
                    records[e["rcsb_id"]] = e

    wrapper = {
        "batch_type": "extra_targets_entries",
        "ids_requested": all_ids_requested,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "success": success,
        "response": {"data": {"entries": combined_entries}},
    }
    EXTRA_TARGETS_RAW_PATH.write_text(json.dumps(wrapper, indent=2, ensure_ascii=False))
    return records


def classify_target(entry_fields: dict | None) -> tuple[str, str]:
    """Returns (scope_exclusion_reason, review_status). Classification is
    based ONLY on authoritative RCSB metadata fields, never on the entry
    title."""
    if entry_fields is None:
        return "REVIEW_NEEDED", "No Data API entry record could be retrieved for this PDB ID."

    methodology = (entry_fields.get("structure_determination_methodology") or "").lower()
    methods = entry_fields.get("experimental_methods") or []

    if methodology == EXPERIMENTAL_METHODOLOGY:
        # This would be unexpected — an "experimental" entry outside the
        # Stage 2 candidate set would need real investigation, not a
        # confident classification.
        return "OTHER_SEARCH_SCOPE_MISMATCH", (
            "REVIEW_NEEDED: structure_determination_methodology is 'experimental' yet this target "
            "was absent from the Stage 2 Search API result — does not fit the integrative/computed "
            "exclusion pattern; requires manual investigation."
        )
    if methodology == "integrative":
        return "INTEGRATIVE_STRUCTURE", "PASS"
    if methodology in ("computational", "computed"):
        return "COMPUTED_STRUCTURE_MODEL", "PASS"
    if not methods and not methodology:
        return "NON_EXPERIMENTAL_MODEL", "REVIEW_NEEDED: no methodology or experimental method reported at all."
    return "OTHER_SEARCH_SCOPE_MISMATCH", (
        f"REVIEW_NEEDED: structure_determination_methodology={methodology!r} does not match a "
        f"recognized exclusion category (integrative/computational)."
    )


def main() -> None:
    alignments_df = pd.read_csv(ALIGNMENTS_CSV, dtype={"pdb_id": str})
    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR_CSV, dtype={"pdb_id": str, "entity_id": str})
    stage2_entity_ids = set(discovery_df["polymer_entity_id"])

    extra_targets = alignments_df[~alignments_df["polymer_entity_id"].isin(stage2_entity_ids)].copy()
    print(f"Found {len(extra_targets)} Sequence Coordinates target(s) outside the Stage 2 candidate inventory.")

    if len(extra_targets) == 0:
        pd.DataFrame(
            columns=[
                "uniprot_id", "nr_code", "common_name", "polymer_entity_id", "pdb_id", "entity_id",
                "structure_title", "deposit_date", "release_date",
                "structure_determination_methodology", "experimental_methods",
                "stage2_scope_status", "scope_exclusion_reason", "review_status", "review_reason",
            ]
        ).to_csv(EXTRA_TARGETS_CSV, index=False)
        print(f"Wrote {EXTRA_TARGETS_CSV.relative_to(PROJECT_ROOT)} (0 rows)")
        return

    unique_pdb_ids = sorted(extra_targets["pdb_id"].unique())
    cached = load_cached_entry_records()
    still_needed = [pid for pid in unique_pdb_ids if pid not in cached]
    fetched = fetch_extra_entry_metadata(still_needed, refresh=False) if still_needed else {}
    entry_records = {**cached, **fetched}

    rows = []
    for _, row in extra_targets.iterrows():
        raw = entry_records.get(row["pdb_id"])
        extracted = extract_entry_fields(raw) if raw is not None else None

        scope_reason, review = classify_target(extracted)

        rows.append(
            {
                "uniprot_id": row["uniprot_id"],
                "nr_code": row["nr_code"],
                "common_name": row["common_name"],
                "polymer_entity_id": row["polymer_entity_id"],
                "pdb_id": row["pdb_id"],
                "entity_id": row["entity_id"],
                "structure_title": extracted["structure_title"] if extracted else "",
                "deposit_date": extracted["deposit_date"] if extracted else "",
                "release_date": extracted["initial_release_date"] if extracted else "",
                "structure_determination_methodology": extracted["structure_determination_methodology"] if extracted else "",
                "experimental_methods": json.dumps(extracted["experimental_methods"]) if extracted else "[]",
                "stage2_scope_status": "OUTSIDE_STAGE2_SEARCH_SCOPE",
                "scope_exclusion_reason": scope_reason,
                "review_status": "PASS" if review == "PASS" else "REVIEW_NEEDED",
                "review_reason": "" if review == "PASS" else review,
            }
        )

    extra_targets_df = pd.DataFrame(rows)
    extra_targets_df.to_csv(EXTRA_TARGETS_CSV, index=False)
    print(f"Wrote {EXTRA_TARGETS_CSV.relative_to(PROJECT_ROOT)} ({len(extra_targets_df)} rows)")
    for _, r in extra_targets_df.iterrows():
        print(f"  {r['common_name']:10s} {r['polymer_entity_id']:10s} methodology={r['structure_determination_methodology']!r} "
              f"deposit={r['deposit_date']} release={r['release_date']} -> {r['scope_exclusion_reason']}")

    # --- correct manifest terminology (rename unexpected_targets column) ---
    manifest_df = pd.read_csv(SEQ_COORD_MANIFEST_CSV)
    if "unexpected_targets" in manifest_df.columns:
        manifest_df = manifest_df.rename(columns={"unexpected_targets": "targets_outside_stage2_scope"})
        manifest_df.to_csv(SEQ_COORD_MANIFEST_CSV, index=False)
        print(f"\nRenamed 'unexpected_targets' -> 'targets_outside_stage2_scope' in "
              f"{SEQ_COORD_MANIFEST_CSV.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
