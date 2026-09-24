"""Extraction of authoritative nuclear-receptor LBD/DBD domain intervals
from a raw UniProtKB JSON entry (the same schema already used in Stage 1).

Coordinate convention: 1-based, inclusive (UniProt's own convention —
`location.start.value` and `location.end.value` are the first and last
residue positions of the feature, both included).

LBD identification rule (Stage 3A Section 2): only a `features[]` entry
with `type == "Domain"` and `description == "NR LBD"` is accepted
automatically. Verified empirically (2026-09-23) against all 48 cached
Stage 1 UniProt responses: every one of the 48 project receptors has
EXACTLY one such feature, evidenced by PROSITE-ProRule (e.g. PRU01189 for
VDR) — no fallback source was required for any receptor.

DBD identification rule (Stage 3A Section 5, diagnostic only): a
`features[]` entry with `type == "DNA binding"` and
`description == "Nuclear receptor"` (PROSITE-ProRule, e.g. PRU00407).
Verified empirically: present for 46/48 receptors. The 2 without one
(DAX-1 / P51843, SHP / Q15466) are real, well-documented atypical orphan
nuclear receptors that lack a classical zinc-finger DNA-binding domain —
this is genuine biology, not missing data, and no interval is invented for
them.
"""

from __future__ import annotations

from typing import Any


def _find_domain_features(raw: dict[str, Any], feature_type: str, description: str) -> list[dict[str, Any]]:
    matches = []
    for f in raw.get("features", []):
        if f.get("type") == feature_type and f.get("description") == description:
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            if start is None or end is None:
                continue
            evidences = f.get("evidences", [])
            evidence_str = "; ".join(
                f"{e.get('source', '?')}:{e.get('id', '?')}" for e in evidences
            ) if evidences else ""
            matches.append({"start": start, "end": end, "evidence": evidence_str})
    return matches


def extract_lbd_feature(raw: dict[str, Any]) -> dict[str, Any] | None:
    """Returns {'start', 'end', 'evidence'} for the single accepted 'NR LBD'
    UniProt Domain feature, or None if zero or more than one was found
    (ambiguous cases must not be silently resolved by the caller)."""
    matches = _find_domain_features(raw, "Domain", "NR LBD")
    if len(matches) != 1:
        return None
    return matches[0]


def extract_dbd_feature(raw: dict[str, Any]) -> dict[str, Any] | None:
    """Returns {'start', 'end', 'evidence'} for the single 'DNA binding' /
    'Nuclear receptor' UniProt feature, or None if absent or ambiguous.
    Absence is a valid, real outcome for atypical orphan receptors (e.g.
    DAX-1, SHP) — callers must not infer a DBD interval in that case."""
    matches = _find_domain_features(raw, "DNA binding", "Nuclear receptor")
    if len(matches) != 1:
        return None
    return matches[0]


def extract_entry_audit(raw: dict[str, Any]) -> dict[str, Any]:
    audit = raw.get("entryAudit", {})
    return {
        "entry_version": audit.get("entryVersion"),
        "sequence_version": audit.get("sequenceVersion"),
    }


def extract_sequence_length(raw: dict[str, Any]) -> int | None:
    return raw.get("sequence", {}).get("length")
