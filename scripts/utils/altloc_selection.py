"""
Final Data Release, Step 12: deterministic raw Ca altloc handoff rule.

For teammate convenience, exactly ONE recommended raw Ca record is chosen
per residue. This is NOT normalized B-factor analysis — B-factor magnitude
never participates in the choice, only occupancy and altloc identity.

Rule, applied only among Ca records with occupancy > 0 AND finite
B_iso_or_equiv:
    1. highest occupancy
    2. if tied: blank/no-altloc ('.' or '') preferred
    3. then altloc 'A'
    4. then lexicographically smallest altloc identifier

If no record satisfies occupancy > 0 and finite B_iso_or_equiv, all
selected fields are None.
"""

from __future__ import annotations

import math


def _altloc_priority(altloc: str) -> tuple[int, str]:
    if altloc in (".", "", "?"):
        return (0, "")
    if altloc == "A":
        return (1, "")
    return (2, altloc)


def select_representative_ca_record(
    occupancies: list, b_isos: list, altlocs: list
) -> dict:
    """occupancies/b_isos/altlocs are parallel lists (one entry per Ca atom
    record at this residue, across all altlocs) as stored in the Stage 3B
    observations parquet's occupancies_json/b_iso_values_json/altloc_ids_json
    columns."""
    candidates = []
    for occ, b_iso, altloc in zip(occupancies, b_isos, altlocs):
        if occ is None or b_iso is None:
            continue
        if not (occ > 0):
            continue
        if not math.isfinite(b_iso):
            continue
        candidates.append((occ, altloc, b_iso))

    if not candidates:
        return {"altloc": None, "occupancy": None, "b_iso_or_equiv": None}

    candidates.sort(key=lambda c: (-c[0], _altloc_priority(c[1])))
    best_occ, best_altloc, best_b_iso = candidates[0]
    return {"altloc": best_altloc, "occupancy": best_occ, "b_iso_or_equiv": best_b_iso}
