"""
Part 9 of the functional-site specificity/validation stage: for the
existing significant apo/holo x functional-site overlaps (already
computed, never rerun -- data/analysis/apo_holo_segment_effects.parquet
and reports/tables/apo_holo_functional_site_overlap.csv, both read-only
here), report whether the overlapping site annotation is residue-level
high confidence, a segment-level proxy, or ambiguous -- so the apo/holo
functional-site link is never overstated as a precise residue-level
mechanistic association when it is really only a coarse segment-level
correspondence.
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def main():
    overlap = pd.read_csv(TABLES / "apo_holo_functional_site_overlap.csv", keep_default_na=False)
    sig_overlap = overlap[(overlap["apo_holo_significant_within_receptor"]) & (overlap["functional_site_overlap"])]

    receptor_site = pd.read_parquet(ANALYSIS / "receptor_functional_site_flexibility.parquet")

    rows = []
    for _, r in sig_overlap.iterrows():
        conf_row = receptor_site[(receptor_site["uniprot_id"] == r["uniprot_id"]) &
                                  (receptor_site["functional_site"] == r["site_name"])]
        if len(conf_row) == 0:
            level, conf, category = "", "", "AMBIGUOUS_NO_DEFINITION_FOUND"
        else:
            level = conf_row.iloc[0]["definition_level"]
            conf = conf_row.iloc[0]["site_definition_confidence"]
            if level == "RESIDUE" and conf == "HIGH":
                category = "RESIDUE_LEVEL_HIGH_CONFIDENCE"
            elif level == "RESIDUE":
                category = "RESIDUE_LEVEL_MODERATE_OR_LOWER"
            elif level == "SEGMENT":
                category = "SEGMENT_LEVEL_PROXY"
            else:
                category = "AMBIGUOUS"

        rows.append({
            "uniprot_id": r["uniprot_id"], "nr_code": r["nr_code"], "apo_holo_segment": r["apo_holo_segment"],
            "site_name": r["site_name"], "apo_holo_delta_median_z": r["apo_holo_delta_median_z"],
            "apo_holo_q_within_receptor": r["apo_holo_q_within_receptor"],
            "site_definition_level": level, "site_definition_confidence": conf,
            "overlap_confidence_category": category,
        })

    recheck = pd.DataFrame(rows)
    out = TABLES / "apo_holo_overlap_confidence_recheck.csv"
    recheck.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(recheck)} rows)")
    print(recheck.to_string(index=False))
    print(
        "\nNo apo/holo statistic was recomputed here -- this is purely a confidence-labeling pass "
        "over the already-stored overlap table, to prevent overstating a segment-level-proxy "
        "overlap as a precise residue-level mechanistic link."
    )


if __name__ == "__main__":
    main()
