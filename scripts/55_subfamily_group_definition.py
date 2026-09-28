"""
Part 3-4 of the subfamily flexibility / segment-atlas stage (continuing
from commit 06f7b7e on analysis/normalization-v1).

Part 3: inspect and report the EXACT existing grouping fields in the
project metadata -- never invent or relabel groups. Three grouping
schemes actually exist in this project's data:
  1. `nr_code` subfamily-letter level (e.g. NR3C4 -> "NR3C") -- the
     finest formal NR nomenclature grouping.
  2. `nr_code` family-number level (e.g. NR3C4 -> "NR3") -- a coarser
     formal NR nomenclature grouping.
  3. `group` -- the professor-supplied project-level functional grouping
     (config/nr_metadata.py: "Steroid hormone receptors" / "Nutritional
     sensors" / "Orphan receptors").
The formal subfamily-letter grouping (1) has inadequate sample sizes for
every group (max n=4, most n<=3) among the 42 coordinate receptors, so it
cannot serve as a usable PRIMARY or SECONDARY test grouping and is
reported for transparency only. Between the two grouping schemes that DO
have workable sample sizes, the project-supplied functional `group` (3)
has the largest, most balanced groups (25/12/5) and is used as PRIMARY;
the coarser formal NR family-number grouping (2) is used as SECONDARY
(with its own small groups -- NR4, NR5, NR0 -- flagged, not merged away).

Part 4: for every realized structural segment and every group in both
schemes, record eligibility (>=3 receptors observed for that segment in
that group) -- entirely per-segment, never a global receptor exclusion.
"""
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

MIN_N_PER_GROUP = 3
MIN_ELIGIBLE_GROUPS = 2


def load_receptor_metadata() -> pd.DataFrame:
    rfs = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")
    coord = rfs[rfs["n_coordinate_candidates"] > 0][
        ["uniprot_id", "nr_code", "common_name", "group", "n_coordinate_candidates"]
    ].copy()
    coord["nr_subfamily_letter"] = coord["nr_code"].str.extract(r"^(NR\d+[A-Z])")
    coord["nr_family_number"] = coord["nr_code"].str.extract(r"^(NR\d+)")
    return coord.reset_index(drop=True)


def report_grouping_schemes(meta: pd.DataFrame):
    print("=" * 70)
    print("Part 3: exact existing grouping fields (42 coordinate receptors)")
    print("=" * 70)
    print("\n1. Formal NR subfamily-letter (from nr_code, e.g. NR3C4 -> NR3C):")
    counts1 = meta["nr_subfamily_letter"].value_counts().sort_index()
    print(counts1.to_string())
    print(f"   -> max group n = {counts1.max()}, groups with n>=3: {(counts1 >= 3).sum()}/{len(counts1)}")
    print("   VERDICT: inadequate sample sizes almost everywhere -- reported for "
          "transparency only, NOT used as a formal test grouping.")

    print("\n2. Formal NR family-number (from nr_code, e.g. NR3C4 -> NR3):")
    counts2 = meta["nr_family_number"].value_counts().sort_index()
    print(counts2.to_string())
    print(f"   -> groups with n>=3: {(counts2 >= 3).sum()}/{len(counts2)} "
          f"(small groups NR4/NR5/NR0 flagged, not merged)")
    print("   VERDICT: usable as SECONDARY grouping.")

    print("\n3. Project-supplied functional `group` (config/nr_metadata.py):")
    counts3 = meta["group"].value_counts()
    print(counts3.to_string())
    print(f"   -> groups with n>=3: {(counts3 >= 3).sum()}/{len(counts3)}")
    print("   VERDICT: largest, most balanced sample sizes -- used as PRIMARY grouping.")

    schemes_out = ANALYSIS / "subfamily_grouping_schemes.json"
    schemes_out.write_text(json.dumps({
        "nr_subfamily_letter_counts": counts1.to_dict(),
        "nr_family_number_counts": counts2.to_dict(),
        "project_group_counts": counts3.to_dict(),
        "primary_scheme": "group",
        "secondary_scheme": "nr_family_number",
        "subfamily_letter_scheme_usable": False,
        "min_n_per_group": MIN_N_PER_GROUP,
        "min_eligible_groups": MIN_ELIGIBLE_GROUPS,
    }, indent=2) + "\n")
    print(f"\nWrote {schemes_out.relative_to(PROJECT_ROOT)}")


def build_eligibility_table(meta: pd.DataFrame, profiles: pd.DataFrame, scheme_col: str, scheme_name: str) -> pd.DataFrame:
    """Per-segment, per-group sample sizes and eligibility -- available-case only,
    no imputation, no global receptor exclusion. A receptor missing segment X but
    observed for segment Y contributes to Y's counts regardless."""
    all_receptors = meta[["uniprot_id", scheme_col]].rename(columns={scheme_col: "group_value"})
    all_segments = sorted(profiles["structural_segment"].unique())

    rows = []
    for seg in all_segments:
        observed_ids = set(profiles.loc[profiles["structural_segment"] == seg, "uniprot_id"])
        for group_value, group_members in all_receptors.groupby("group_value"):
            n_total = len(group_members)
            n_observed = group_members["uniprot_id"].isin(observed_ids).sum()
            n_missing = n_total - n_observed
            rows.append({
                "grouping_scheme": scheme_name,
                "structural_segment": seg,
                "group": group_value,
                "n_total_group_receptors": n_total,
                "n_observed": int(n_observed),
                "n_missing": int(n_missing),
                "eligible_for_test": bool(n_observed >= MIN_N_PER_GROUP),
            })
    return pd.DataFrame(rows)


def main():
    meta = load_receptor_metadata()
    report_grouping_schemes(meta)

    profiles = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")

    primary = build_eligibility_table(meta, profiles, "group", "primary_project_group")
    secondary = build_eligibility_table(meta, profiles, "nr_family_number", "secondary_nr_family_number")
    combined = pd.concat([primary, secondary], ignore_index=True)

    out = TABLES / "subfamily_segment_sample_sizes.csv"
    combined.to_csv(out, index=False)
    print(f"\nPart 4: wrote {out.relative_to(PROJECT_ROOT)} ({len(combined)} rows)")

    for scheme_name, table in [("primary_project_group", primary), ("secondary_nr_family_number", secondary)]:
        n_eligible_segments = table[table["eligible_for_test"]].groupby("structural_segment")["group"].nunique()
        n_segments_testable = (n_eligible_segments >= MIN_ELIGIBLE_GROUPS).sum()
        print(f"\n{scheme_name}: {n_segments_testable}/{len(sorted(profiles['structural_segment'].unique()))} "
              f"segments have >= {MIN_ELIGIBLE_GROUPS} eligible groups (n>={MIN_N_PER_GROUP} each) -- testable")

    meta_out = ANALYSIS / "subfamily_receptor_metadata.csv"
    meta.to_csv(meta_out, index=False)
    print(f"\nWrote {meta_out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
