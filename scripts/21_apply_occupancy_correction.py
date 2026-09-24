"""
Stage 3B.1: correct the semantics of "observed Cα" from mere atom-record
presence to POSITIVE OCCUPANCY, without re-parsing any mmCIF file.

Background: Stage 3B's `ca_modeled` field only checked whether a Cα
`_atom_site` record existed, conflating genuinely observed atoms with
zero-occupancy placeholder atoms (270 such cases found in the Stage 3B
cross-check). This script derives the three explicitly separate Cα states
from data ALREADY extracted into the existing observation parquet
(`occupancies_json`, `b_iso_values_json`, `altloc_ids_json` per residue) —
no mmCIF file needs to be re-downloaded or re-parsed.

Three states (all altloc-aware; a residue is never double-counted):
    A. ca_record_present      = at least one Cα atom-site record exists
                                 (this IS the old `ca_modeled`; kept as an
                                 explicit alias, documented as "atom record
                                 present, including zero-occupancy
                                 placeholders" — no longer used alone to
                                 mean "observed")
    B. ca_positive_occupancy  = at least one Cα record has occupancy > 0
                                 (THIS is now what "observed Cα" means)
    C. ca_bfactor_usable      = at least one Cα record has occupancy > 0
                                 AND a finite numeric B_iso_or_equiv
                                 (unchanged definition from Stage 3B —
                                 it already required positive occupancy)

Reads/writes in place:
    data/interim/qc/stage3b_lbd_ca_observations.parquet
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

INTERIM_QC_DIR = PROJECT_ROOT / "data" / "interim" / "qc"
OBSERVATIONS_PARQUET = INTERIM_QC_DIR / "stage3b_lbd_ca_observations.parquet"

BATCH_ROWS = 200_000


def compute_occupancy_fields(occupancies_json: str, b_iso_values_json: str) -> dict:
    occupancies = json.loads(occupancies_json) if occupancies_json else []
    b_isos = json.loads(b_iso_values_json) if b_iso_values_json else []

    ca_record_present = len(occupancies) > 0
    positive_occ_flags = [o is not None and o > 0 for o in occupancies]
    ca_positive_occupancy = any(positive_occ_flags)
    positive_occupancy_ca_altloc_count = sum(positive_occ_flags)

    ca_bfactor_usable = any(
        (o is not None and o > 0) and (b is not None and math.isfinite(b))
        for o, b in zip(occupancies, b_isos)
    )
    zero_occupancy_only_ca = ca_record_present and not ca_positive_occupancy

    return {
        "ca_record_present": ca_record_present,
        "ca_positive_occupancy": ca_positive_occupancy,
        "ca_bfactor_usable": ca_bfactor_usable,
        "zero_occupancy_only_ca": zero_occupancy_only_ca,
        "positive_occupancy_ca_altloc_count": positive_occupancy_ca_altloc_count,
    }


def main() -> None:
    print(f"Reading {OBSERVATIONS_PARQUET.relative_to(PROJECT_ROOT)} in batches...")
    parquet_file = pq.ParquetFile(OBSERVATIONS_PARQUET)
    total_rows = parquet_file.metadata.num_rows
    print(f"Total rows: {total_rows}")

    out_schema = None
    writer = None
    tmp_path = OBSERVATIONS_PARQUET.with_suffix(".parquet.tmp")

    rows_done = 0
    mismatch_count = 0

    for batch in parquet_file.iter_batches(batch_size=BATCH_ROWS):
        df = batch.to_pandas()

        computed = df.apply(
            lambda r: compute_occupancy_fields(r["occupancies_json"], r["b_iso_values_json"]), axis=1
        )
        computed_df = pd.DataFrame(list(computed))

        # Sanity: ca_record_present must equal the existing ca_modeled for
        # every row (they are defined identically) — verify rather than
        # assume, since a silent mismatch here would indicate the original
        # extraction was inconsistent with the raw JSON now being re-derived
        # from.
        mismatches = (computed_df["ca_record_present"] != df["ca_modeled"]).sum()
        mismatch_count += int(mismatches)

        df["ca_record_present"] = computed_df["ca_record_present"]
        df["ca_positive_occupancy"] = computed_df["ca_positive_occupancy"]
        df["ca_bfactor_usable"] = computed_df["ca_bfactor_usable"]
        df["zero_occupancy_only_ca"] = computed_df["zero_occupancy_only_ca"]
        df["positive_occupancy_ca_altloc_count"] = computed_df["positive_occupancy_ca_altloc_count"]

        table = pa.Table.from_pandas(df, preserve_index=False)
        if writer is None:
            out_schema = table.schema
            writer = pq.ParquetWriter(tmp_path, out_schema)
        writer.write_table(table)

        rows_done += len(df)
        print(f"  processed {rows_done}/{total_rows} rows")

    if writer is not None:
        writer.close()

    if mismatch_count > 0:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"STOPPING: {mismatch_count} rows had ca_record_present != legacy ca_modeled — "
            f"the two are defined identically and must always agree. Not overwriting the "
            f"observation table with inconsistent data."
        )

    tmp_path.replace(OBSERVATIONS_PARQUET)
    print(f"\nWrote corrected {OBSERVATIONS_PARQUET.relative_to(PROJECT_ROOT)} "
          f"({rows_done} rows, 0 ca_record_present/ca_modeled mismatches)")


if __name__ == "__main__":
    main()
