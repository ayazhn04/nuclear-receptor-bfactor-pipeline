"""Thin wrapper around the shared RCSB Data API batching engine
(scripts/utils/rcsb_data_client.py) for polymer-entity-instance metadata.

Schema confirmed via live GraphQL introspection (2026-09-23):
    CorePolymerEntity.polymer_entity_instances: [CorePolymerEntityInstance]
    CorePolymerEntityInstance.rcsb_polymer_entity_instance_container_identifiers
        { entry_id, entity_id, asym_id, auth_asym_id }
    Instance rcsb_id format: "<ENTRY_ID>.<ASYM_ID>" (label asym), e.g. "1KB2.D".

Pilot-verified against a known multi-chain entity (1KB2_3, a VDR DBD dimer):
returns 2 instances with an UNAMBIGUOUS, authoritative label<->auth asym_id
pairing (asym_id "C" -> auth_asym_id "A", asym_id "D" -> auth_asym_id "B") —
obtained from true instance-level records, never by pairing the Stage 2
entity-level `label_asym_ids`/`auth_asym_ids` arrays by position.

This query is nested under `polymer_entities(entity_ids: [...])`, so the
existing Stage 2 polymer-entity IDs can be reused directly — no separate
instance-ID discovery step is needed.
"""

from __future__ import annotations

from scripts.utils.rcsb_data_client import RCSBBatchResult, fetch_batches


def fetch_instances_for_entities(entity_ids: list[str], batch_size: int = 150) -> list[RCSBBatchResult]:
    return fetch_batches("polymer_entity_instances", entity_ids, batch_size=batch_size)
