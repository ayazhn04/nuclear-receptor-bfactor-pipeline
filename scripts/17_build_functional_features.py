"""
Stage 3A.2: capture authoritative functional-feature intervals (UniProt
`Motif` features, e.g. AF-2/9aaTAD transactivation motifs, LXXLL
coactivator-binding motifs, nuclear localization/export signals) separately
from the standardized LBD analysis interval, so functionally important
residues are never silently lost just because they sit outside the
standardized profile boundary — and never silently used to expand it
either.

Also preserves the VDR terminal Helix feature (416-422), commonly referred
to as "H12" in nuclear-receptor literature, though UniProt itself does not
apply that label to this feature (evidenced only as an unnamed
PDB-derived secondary-structure Helix, source PDB:3B0T) — recorded exactly
as UniProt provides it, not renamed to assert an identity UniProt doesn't
state.

Only features already present in the cached Stage 1 raw UniProt responses
are used. No literature curation, no invented coordinates: where a
receptor has no UniProt Motif feature, it simply has no row here.

Reads:
    config/nr_metadata.py
    data/raw/api/uniprot/{uniprot_id}.json
    data/manifests/nr_lbd_reference.csv   (for within/outside comparison only — not modified)

Writes:
    data/manifests/nr_functional_features.csv
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.nr_metadata import NR_METADATA  # noqa: E402

RAW_UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "api" / "uniprot"
MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
LBD_REFERENCE_CSV = MANIFESTS_DIR / "nr_lbd_reference.csv"
OUTPUT_CSV = MANIFESTS_DIR / "nr_functional_features.csv"

COLUMNS = [
    "uniprot_id", "nr_code", "common_name",
    "feature_name", "feature_type",
    "feature_start", "feature_end",
    "source", "source_identifier",
    "within_standardized_lbd", "extends_outside_standardized_lbd",
    "notes",
]

# VDR-specific terminal-helix exception (Section 6): UniProt records this
# secondary-structure Helix but assigns it no biological name. We preserve
# it explicitly, cautiously worded, without asserting a name UniProt does
# not itself state.
VDR_TERMINAL_HELIX_UNIPROT_ID = "P11473"


def evidence_str(feature: dict) -> str:
    evidences = feature.get("evidences", [])
    return "; ".join(f"{e.get('source', '?')}:{e.get('id', '?')}" for e in evidences) if evidences else ""


def main() -> None:
    lbd_ref = pd.read_csv(LBD_REFERENCE_CSV).set_index("uniprot_id")

    rows = []
    for record in NR_METADATA:
        uid = record["uniprot_id"]
        raw_path = RAW_UNIPROT_DIR / f"{uid}.json"
        if not raw_path.exists():
            continue
        raw = json.loads(raw_path.read_text())

        lbd_row = lbd_ref.loc[uid] if uid in lbd_ref.index else None
        lbd_start = int(lbd_row["lbd_start"]) if lbd_row is not None and pd.notna(lbd_row["lbd_start"]) else None
        lbd_end = int(lbd_row["lbd_end"]) if lbd_row is not None and pd.notna(lbd_row["lbd_end"]) else None

        for f in raw.get("features", []):
            if f.get("type") != "Motif":
                continue
            start = f.get("location", {}).get("start", {}).get("value")
            end = f.get("location", {}).get("end", {}).get("value")
            if start is None or end is None:
                continue

            if lbd_start is not None and lbd_end is not None:
                within = start >= lbd_start and end <= lbd_end
                outside = end > lbd_end or start < lbd_start
            else:
                within, outside = "", "REVIEW_NEEDED"

            evidences = f.get("evidences", [])
            source = evidences[0].get("source", "UniProt") if evidences else "UniProt"

            rows.append(
                {
                    "uniprot_id": uid, "nr_code": record["nr_code"], "common_name": record["common_name"],
                    "feature_name": f.get("description", ""), "feature_type": "Motif",
                    "feature_start": start, "feature_end": end,
                    "source": f"UniProt ({source})", "source_identifier": evidence_str(f),
                    "within_standardized_lbd": within, "extends_outside_standardized_lbd": outside,
                    "notes": "",
                }
            )

        # VDR terminal helix exception.
        if uid == VDR_TERMINAL_HELIX_UNIPROT_ID:
            for f in raw.get("features", []):
                if f.get("type") == "Helix" and f.get("location", {}).get("start", {}).get("value") == 416 \
                        and f.get("location", {}).get("end", {}).get("value") == 422:
                    within = lbd_start is not None and 416 >= lbd_start and 422 <= lbd_end
                    rows.append(
                        {
                            "uniprot_id": uid, "nr_code": record["nr_code"], "common_name": record["common_name"],
                            "feature_name": "Terminal helix (commonly referred to as H12 in NR literature; "
                                            "UniProt does not itself apply this label)",
                            "feature_type": "Helix",
                            "feature_start": 416, "feature_end": 422,
                            "source": "UniProt (PDB)", "source_identifier": evidence_str(f),
                            "within_standardized_lbd": within, "extends_outside_standardized_lbd": not within,
                            "notes": (
                                "Preserved per Stage 3A.2 Section 6 (functional-feature separation) so this "
                                "region is not lost from downstream functional analysis even though the "
                                "standardized LBD interval (127-423) already fully contains it."
                            ),
                        }
                    )

    df = pd.DataFrame(rows, columns=COLUMNS)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {OUTPUT_CSV.relative_to(PROJECT_ROOT)} ({len(df)} rows)")
    print(f"Receptors with at least one functional feature: {df['uniprot_id'].nunique()} / 48")
    print(f"Features extending outside the standardized LBD interval: {int((df['extends_outside_standardized_lbd'] == True).sum())}")


if __name__ == "__main__":
    main()
