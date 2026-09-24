"""
Stage 3A, step 3: fetch authoritative polymer-instance (chain) metadata for
every unique Stage 2 candidate polymer entity, via the official RCSB Data
API, and resolve label_asym_id <-> auth_asym_id unambiguously (never by
array-position pairing of the Stage 2 entity-level lists).

Reads:
    data/manifests/rcsb_discovery_by_receptor.csv   (2072 unique polymer entity IDs)
    data/manifests/rcsb_polymer_entities.csv          (entity-level label/auth lists, for the anomaly audit)

Writes:
    data/raw/api/rcsb/instances/instances_batch_XXXX.json   (gitignored)
    data/manifests/rcsb_instance_api_manifest.csv
    data/manifests/rcsb_polymer_instances.csv                 (one row per instance)
    data/manifests/rcsb_instance_mapping_anomalies.csv

Idempotent by default: cached raw batch files are reused unless --refresh.
No coordinate files are downloaded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.rcsb_data_client import RCSB_DATA_API_URL  # noqa: E402
from scripts.utils.rcsb_instance_client import fetch_instances_for_entities  # noqa: E402

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
RAW_INSTANCES_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "rcsb" / "instances"

DISCOVERY_CSV = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
POLYMER_ENTITIES_CSV = MANIFESTS_DIR / "rcsb_polymer_entities.csv"

API_MANIFEST_CSV = MANIFESTS_DIR / "rcsb_instance_api_manifest.csv"
INSTANCES_CSV = MANIFESTS_DIR / "rcsb_polymer_instances.csv"
ANOMALIES_CSV = MANIFESTS_DIR / "rcsb_instance_mapping_anomalies.csv"

DEFAULT_BATCH_SIZE = 150


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_cached_entity_ids() -> set[str]:
    covered: set[str] = set()
    for path in sorted(RAW_INSTANCES_DIR.glob("instances_batch_*.json")):
        try:
            wrapper = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        covered.update(wrapper.get("ids_requested", []))
    return covered


def load_cached_entity_records() -> dict[str, dict]:
    records: dict[str, dict] = {}
    for path in sorted(RAW_INSTANCES_DIR.glob("instances_batch_*.json")):
        wrapper = json.loads(path.read_text())
        if not wrapper.get("success"):
            continue
        entities = ((wrapper.get("response") or {}).get("data") or {}).get("polymer_entities") or []
        for e in entities:
            if e and e.get("rcsb_id"):
                records[e["rcsb_id"]] = e
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Force re-fetching every batch from the live API.")
    args = parser.parse_args()

    discovery_df = pd.read_csv(DISCOVERY_CSV, dtype={"pdb_id": str, "entity_id": str})
    unique_entity_ids = sorted(discovery_df["polymer_entity_id"].unique())
    print(f"Unique Stage 2 polymer entities to fetch instances for: {len(unique_entity_ids)}")

    RAW_INSTANCES_DIR.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    if not args.refresh:
        cached_ids = load_cached_entity_ids()
        missing = [i for i in unique_entity_ids if i not in cached_ids]
    else:
        missing = unique_entity_ids
        for path in RAW_INSTANCES_DIR.glob("instances_batch_*.json"):
            path.unlink()

    if missing:
        print(f"Fetching {len(missing)} of {len(unique_entity_ids)} entities not yet cached...")
        results = fetch_instances_for_entities(missing, batch_size=DEFAULT_BATCH_SIZE)
        for r in results:
            raw_path = RAW_INSTANCES_DIR / f"instances_batch_{r.batch_number:04d}.json"
            wrapper = {
                "batch_type": "polymer_entity_instances",
                "batch_number": r.batch_number,
                "ids_requested": r.ids_requested,
                "retrieved_at_utc": r.retrieved_at_utc,
                "http_status": r.http_status,
                "success": r.success,
                "error_message": r.error_message,
                "response": r.raw_json,
            }
            raw_path.write_bytes(json.dumps(wrapper, indent=2, ensure_ascii=False).encode("utf-8"))
            status_label = "OK" if r.success else "FAILED"
            print(f"  batch {r.batch_number:3d} ({len(r.ids_requested):3d} ids) -> {status_label}"
                  + (f"  [{r.error_message}]" if not r.success else ""))
    else:
        print("All entities already covered by cached batches — reusing.")

    for path in sorted(RAW_INSTANCES_DIR.glob("instances_batch_*.json")):
        wrapper = json.loads(path.read_text())
        raw_bytes = path.read_bytes()
        manifest_rows.append(
            {
                "batch_number": wrapper["batch_number"],
                "number_requested": len(wrapper["ids_requested"]),
                "request_url": RCSB_DATA_API_URL,
                "http_status": wrapper["http_status"] if wrapper["http_status"] is not None else "",
                "retrieved_at_utc": wrapper["retrieved_at_utc"],
                "raw_response_file": str(path.relative_to(PROJECT_ROOT)),
                "raw_response_sha256": sha256_bytes(raw_bytes),
                "success": wrapper["success"],
                "error_message": wrapper["error_message"],
            }
        )
    manifest_df = pd.DataFrame(manifest_rows).sort_values("batch_number")
    manifest_df.to_csv(API_MANIFEST_CSV, index=False)
    print(f"\nWrote {API_MANIFEST_CSV.relative_to(PROJECT_ROOT)} ({len(manifest_df)} batch rows)")

    entity_records = load_cached_entity_records()
    requested_but_missing = [e for e in unique_entity_ids if e not in entity_records]

    instance_rows = []
    for entity_id in unique_entity_ids:
        raw = entity_records.get(entity_id)
        if raw is None:
            instance_rows.append(
                {
                    "pdb_id": entity_id.split("_", 1)[0] if "_" in entity_id else "",
                    "polymer_entity_id": entity_id,
                    "entity_id": entity_id.split("_", 1)[1] if "_" in entity_id else "",
                    "instance_id": "",
                    "label_asym_id": "", "auth_asym_id": "",
                    "instance_metadata_status": "REVIEW_NEEDED",
                    "instance_metadata_notes": "No Data API record returned for this entity in any batch.",
                }
            )
            continue

        instances = raw.get("polymer_entity_instances") or []
        if not instances:
            instance_rows.append(
                {
                    "pdb_id": entity_id.split("_", 1)[0],
                    "polymer_entity_id": entity_id,
                    "entity_id": entity_id.split("_", 1)[1],
                    "instance_id": "",
                    "label_asym_id": "", "auth_asym_id": "",
                    "instance_metadata_status": "REVIEW_NEEDED",
                    "instance_metadata_notes": "Entity record returned but zero polymer_entity_instances listed.",
                }
            )
            continue

        for inst in instances:
            ids = inst.get("rcsb_polymer_entity_instance_container_identifiers") or {}
            instance_rows.append(
                {
                    "pdb_id": ids.get("entry_id", ""),
                    "polymer_entity_id": entity_id,
                    "entity_id": ids.get("entity_id", ""),
                    "instance_id": inst.get("rcsb_id", ""),
                    "label_asym_id": ids.get("asym_id", ""),
                    "auth_asym_id": ids.get("auth_asym_id", ""),
                    "instance_metadata_status": "PASS",
                    "instance_metadata_notes": "",
                }
            )

    instances_df = pd.DataFrame(instance_rows)
    instances_df.to_csv(INSTANCES_CSV, index=False)
    print(f"Wrote {INSTANCES_CSV.relative_to(PROJECT_ROOT)} ({len(instances_df)} instance rows)")

    # --- instance mapping anomaly audit (Section 15) ---
    entities_df = pd.read_csv(POLYMER_ENTITIES_CSV)
    entities_by_id = entities_df.set_index("polymer_entity_id")

    anomaly_rows = []
    counts_per_entity = instances_df[instances_df["instance_metadata_status"] == "PASS"].groupby("polymer_entity_id").size()

    for entity_id, group in instances_df.groupby("polymer_entity_id"):
        pass_group = group[group["instance_metadata_status"] == "PASS"]

        dup_label = pass_group[pass_group.duplicated("label_asym_id", keep=False) & (pass_group["label_asym_id"] != "")]
        for _, row in dup_label.drop_duplicates("label_asym_id").iterrows():
            anomaly_rows.append(
                {"anomaly_type": "DUPLICATE_LABEL_ASYM_ID", "polymer_entity_id": entity_id,
                 "pdb_id": row["pdb_id"], "details": f"label_asym_id {row['label_asym_id']!r} appears more than once."}
            )

        dup_auth = pass_group[pass_group.duplicated("auth_asym_id", keep=False) & (pass_group["auth_asym_id"] != "")]
        for _, row in dup_auth.drop_duplicates("auth_asym_id").iterrows():
            anomaly_rows.append(
                {"anomaly_type": "DUPLICATE_AUTH_ASYM_ID", "polymer_entity_id": entity_id,
                 "pdb_id": row["pdb_id"], "details": f"auth_asym_id {row['auth_asym_id']!r} appears more than once."}
            )

        missing_auth = pass_group[pass_group["auth_asym_id"] == ""]
        for _, row in missing_auth.iterrows():
            anomaly_rows.append(
                {"anomaly_type": "MISSING_AUTH_ASYM_ID", "polymer_entity_id": entity_id,
                 "pdb_id": row["pdb_id"], "details": f"instance {row['instance_id']} has no auth_asym_id."}
            )
        missing_label = pass_group[pass_group["label_asym_id"] == ""]
        for _, row in missing_label.iterrows():
            anomaly_rows.append(
                {"anomaly_type": "MISSING_LABEL_ASYM_ID", "polymer_entity_id": entity_id,
                 "pdb_id": row["pdb_id"], "details": f"instance {row['instance_id']} has no label_asym_id."}
            )

        # Cross-check against Stage 2 entity-level lists.
        if entity_id in entities_by_id.index:
            entity_row = entities_by_id.loc[entity_id]
            stage2_label = set(json.loads(entity_row["label_asym_ids"])) if entity_row["label_asym_ids"] else set()
            stage2_auth = set(json.loads(entity_row["auth_asym_ids"])) if entity_row["auth_asym_ids"] else set()
            instance_label = set(pass_group["label_asym_id"]) - {""}
            instance_auth = set(pass_group["auth_asym_id"]) - {""}
            if stage2_label != instance_label:
                anomaly_rows.append(
                    {"anomaly_type": "LABEL_ASYM_SET_MISMATCH_VS_STAGE2_ENTITY", "polymer_entity_id": entity_id,
                     "pdb_id": entity_row["pdb_id"],
                     "details": f"Stage2 entity label_asym_ids={sorted(stage2_label)} vs instance-level={sorted(instance_label)}"}
                )
            if stage2_auth != instance_auth:
                anomaly_rows.append(
                    {"anomaly_type": "AUTH_ASYM_SET_MISMATCH_VS_STAGE2_ENTITY", "polymer_entity_id": entity_id,
                     "pdb_id": entity_row["pdb_id"],
                     "details": f"Stage2 entity auth_asym_ids={sorted(stage2_auth)} vs instance-level={sorted(instance_auth)}"}
                )

    for entity_id in requested_but_missing:
        anomaly_rows.append(
            {"anomaly_type": "NO_DATA_API_RECORD", "polymer_entity_id": entity_id,
             "pdb_id": entity_id.split("_", 1)[0] if "_" in entity_id else "",
             "details": "Entity was requested but never returned by the Data API in any batch."}
        )

    anomalies_df = pd.DataFrame(anomaly_rows, columns=["anomaly_type", "polymer_entity_id", "pdb_id", "details"])
    anomalies_df.to_csv(ANOMALIES_CSV, index=False)
    print(f"Wrote {ANOMALIES_CSV.relative_to(PROJECT_ROOT)} ({len(anomalies_df)} anomaly rows)")

    print(f"\n--- Instance fetch summary ---")
    print(f"unique entities processed: {len(unique_entity_ids)}")
    print(f"total instance rows: {len(instances_df)}")
    print(f"entities with exactly 1 instance: {int((counts_per_entity == 1).sum())}")
    print(f"entities with >1 instance: {int((counts_per_entity > 1).sum())}")
    print(f"max instances for one entity: {int(counts_per_entity.max()) if len(counts_per_entity) else 0}")
    print(f"anomalies: {anomalies_df['anomaly_type'].value_counts().to_dict() if len(anomalies_df) else {}}")


if __name__ == "__main__":
    main()
