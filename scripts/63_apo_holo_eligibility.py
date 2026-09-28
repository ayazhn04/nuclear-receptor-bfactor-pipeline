"""
Part 1 of the APO vs HOLO flexibility stage (continuing from commit
f0c085c on analysis/normalization-v1).

Uses the current ligand-state annotations (data/processed/structures.csv)
and the already-approved PDB-level aggregation
(data/analysis/pdb_segment_flexibility_features.parquet, itself built from
the approved chain/instance -> per-instance Z -> PDB-entry-median
hierarchy) to audit, per receptor, how many coordinate PDB entries are
APO, HOLO, AMBIGUOUS, or NOT_ASSESSABLE for a Z-score apo/holo comparison.

NOT_ASSESSABLE here means: the PDB entry has Stage-3B coordinate data and
a resolved ligand_state in structures.csv, but it could not be Z-normalized
(fewer than 2 usable B-factor values, or zero population SD) and therefore
never appears in pdb_segment_flexibility_features.parquet at all -- it is
not usable for this specific comparison, independent of its ligand state.

Two eligibility tiers, fixed BEFORE any test is run (never changed after
seeing significance):
  PRIMARY:     >=5 APO and >=5 HOLO PDB entries
  EXPLORATORY: >=3 APO and >=3 HOLO PDB entries
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

PRIMARY_MIN_APO = 5
PRIMARY_MIN_HOLO = 5
EXPLORATORY_MIN_APO = 3
EXPLORATORY_MIN_HOLO = 3


def main():
    structures = pd.read_csv(PROCESSED / "structures.csv", dtype={"pdb_id": str}, keep_default_na=False)
    coord = structures[structures["coordinate_data_status"] == "COORDINATE_DATA_AVAILABLE"]
    # one row per receptor x PDB entry (a receptor's candidates can include >1 polymer
    # entity per PDB in rare fusion cases; the PDB entry is the analysis unit here)
    coord_pdbs = coord.groupby(["uniprot_id", "nr_code", "common_name", "group", "pdb_id"]).agg(
        ligand_state=("ligand_state", "first")
    ).reset_index()

    features = pd.read_parquet(ANALYSIS / "pdb_segment_flexibility_features.parquet")
    normalizable_pdbs = features[["uniprot_id", "pdb_id", "ligand_state"]].drop_duplicates(subset=["uniprot_id", "pdb_id"])
    normalizable_pdbs = normalizable_pdbs.rename(columns={"ligand_state": "normalized_ligand_state"})

    merged = coord_pdbs.merge(normalizable_pdbs, on=["uniprot_id", "pdb_id"], how="left")
    # sanity: where both are present, ligand_state must agree (same source truth)
    both = merged.dropna(subset=["normalized_ligand_state"])
    mismatches = both[both["ligand_state"] != both["normalized_ligand_state"]]
    if len(mismatches):
        print(f"WARNING: {len(mismatches)} PDB entries have disagreeing ligand_state between "
              "structures.csv and pdb_segment_flexibility_features.parquet -- treating as AMBIGUOUS")
    merged["effective_state"] = merged["normalized_ligand_state"]
    merged.loc[mismatches.index, "effective_state"] = "AMBIGUOUS"
    merged.loc[merged["normalized_ligand_state"].isna(), "effective_state"] = "NOT_ASSESSABLE"

    rows = []
    for (uid, nr_code, common_name, group), g in merged.groupby(["uniprot_id", "nr_code", "common_name", "group"]):
        counts = g["effective_state"].value_counts()
        n_apo = int(counts.get("APO", 0))
        n_holo = int(counts.get("HOLO", 0))
        n_ambiguous = int(counts.get("AMBIGUOUS", 0))
        n_not_assessable = int(counts.get("NOT_ASSESSABLE", 0))
        primary_eligible = n_apo >= PRIMARY_MIN_APO and n_holo >= PRIMARY_MIN_HOLO
        exploratory_eligible = n_apo >= EXPLORATORY_MIN_APO and n_holo >= EXPLORATORY_MIN_HOLO
        rows.append({
            "uniprot_id": uid, "nr_code": nr_code, "common_name": common_name, "group": group,
            "n_total_coordinate_pdbs": len(g),
            "n_APO": n_apo, "n_HOLO": n_holo, "n_AMBIGUOUS": n_ambiguous, "n_NOT_ASSESSABLE": n_not_assessable,
            "eligible_primary_5_5": primary_eligible,
            "eligible_exploratory_3_3": exploratory_eligible,
        })

    eligibility = pd.DataFrame(rows).sort_values("uniprot_id").reset_index(drop=True)
    out = TABLES / "apo_holo_receptor_eligibility.csv"
    eligibility.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(eligibility)} receptors)")

    n_primary = eligibility["eligible_primary_5_5"].sum()
    n_exploratory = eligibility["eligible_exploratory_3_3"].sum()
    print(f"\nPRIMARY tier (>= {PRIMARY_MIN_APO} APO and >= {PRIMARY_MIN_HOLO} HOLO): "
          f"{n_primary}/{len(eligibility)} receptors eligible")
    print(f"EXPLORATORY tier (>= {EXPLORATORY_MIN_APO} APO and >= {EXPLORATORY_MIN_HOLO} HOLO): "
          f"{n_exploratory}/{len(eligibility)} receptors eligible")
    if n_primary < 5:
        print(f"NOTE: the PRIMARY tier leaves only {n_primary} receptor(s) -- reported explicitly. "
              "The EXPLORATORY (>=3/3) tier is used as a genuinely exploratory supplement, "
              "not a silent replacement for the pre-declared PRIMARY threshold.")

    vdr = eligibility[eligibility["uniprot_id"] == "P11473"]
    if len(vdr):
        r = vdr.iloc[0]
        print(f"\nVDR (P11473): {r['n_APO']} APO, {r['n_HOLO']} HOLO -> "
              f"{'ELIGIBLE' if r['eligible_primary_5_5'] else 'INELIGIBLE'} for apo/holo comparison "
              "(reported explicitly; not removed from structures.csv or any master dataset).")

    print("\nPRIMARY-eligible receptors:")
    print(eligibility.loc[eligibility["eligible_primary_5_5"], ["uniprot_id", "nr_code", "n_APO", "n_HOLO"]].to_string(index=False))
    print("\nEXPLORATORY-only receptors (eligible at 3/3 but not 5/5):")
    excl_only = eligibility[eligibility["eligible_exploratory_3_3"] & ~eligibility["eligible_primary_5_5"]]
    print(excl_only[["uniprot_id", "nr_code", "n_APO", "n_HOLO"]].to_string(index=False))


if __name__ == "__main__":
    main()
