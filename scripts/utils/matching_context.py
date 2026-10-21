"""
Matched-background context for script 85 (approved decision S2, Phase 2 of the
v2 rebuild).

Canonical identity comes from the professor map; local structure comes from the
independent script-34 consensus (HELIX / SHEET / COIL). Neither is derived from
the other.

matching_stratum (tier 1)
    the professor element at a position inside a professor helix (a dual position
    keeps one row per element, so it belongs to both); otherwise the derived
    stretch between the flanking professor elements, "INTER:<before>-<after>"
    ("Nterm"/"Cterm" at the span ends). The stretch is a MATCHING DEVICE ONLY,
    never a canonical or biological label, and is not written to any element table.
local_ss_class (tier 2)
    script-34 consensus class; "" where no representative structure observed the
    position (NO_OBSERVATION), so such positions never match on class.
Receptors without a professor map get no rows (they cannot be matched).
"""

from __future__ import annotations

import pandas as pd

INTER_PREFIX = "INTER:"
NO_OBSERVATION = "NO_OBSERVATION"
LOCAL_CLASSES = {"HELIX", "SHEET", "COIL"}
KEY = ["uniprot_id", "canonical_uniprot_position"]


def inter_element_strata(seg_map: pd.DataFrame) -> pd.DataFrame:
    """One row per inter-helix position of a mapped receptor: KEY + matching_stratum."""
    rows = []
    mapped = seg_map[seg_map["structural_class"] != "NO_PROFESSOR_MAPPING"]
    for uid, g in mapped.groupby("uniprot_id"):
        labels = g.groupby("canonical_uniprot_position")["structural_segment"].agg(lambda s: [x for x in s if x])
        items = list(labels.items())
        for i, (pos, lab) in enumerate(items):
            if lab:
                continue
            before = next((l[-1] for _, l in reversed(items[:i]) if l), "Nterm")
            after = next((l[0] for _, l in items[i + 1:] if l), "Cterm")
            rows.append({"uniprot_id": uid, "canonical_uniprot_position": int(pos),
                         "matching_stratum": f"{INTER_PREFIX}{before}-{after}"})
    return pd.DataFrame(rows, columns=KEY + ["matching_stratum"])


def build_matching_frame(seg_map: pd.DataFrame, ss_consensus: pd.DataFrame) -> pd.DataFrame:
    """Rows usable by script 85, in its own column names: structural_segment holds the
    matching stratum and structural_class the local class. Dual positions keep 2 rows."""
    helix = seg_map[seg_map["structural_segment"] != ""][KEY + ["structural_segment"]].rename(
        columns={"structural_segment": "matching_stratum"})
    frame = pd.concat([helix, inter_element_strata(seg_map)], ignore_index=True)
    cls = ss_consensus[KEY + ["consensus_dssp"]]
    unknown = set(cls["consensus_dssp"]) - LOCAL_CLASSES
    if unknown:
        raise ValueError(f"unexpected local classes in script-34 consensus: {unknown}")
    frame = frame.merge(cls, on=KEY, how="left", validate="many_to_one")
    frame["local_ss_class"] = frame["consensus_dssp"].fillna("")
    frame = frame.drop(columns="consensus_dssp")
    frame["canonical_uniprot_position"] = frame["canonical_uniprot_position"].astype(int)
    return frame.rename(columns={"matching_stratum": "structural_segment", "local_ss_class": "structural_class"}) \
        .sort_values(KEY + ["structural_segment"]).reset_index(drop=True)
