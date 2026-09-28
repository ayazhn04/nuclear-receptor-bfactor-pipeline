"""
Part 13 of the current functional-site flexibility stage: a DESCRIPTIVE
overlap between the already-completed apo/holo segment-level results
(data/analysis/apo_holo_segment_effects.parquet, commit 4212a5f) and the
current functional-site definitions. The apo/holo tests are NOT rerun
here -- this script only reads their stored output and asks whether a
tested apo/holo segment corresponds (by canonical position overlap) to a
functional site. No causal inference is drawn from the overlap.
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"


def main():
    apo_holo = pd.read_parquet(ANALYSIS / "apo_holo_segment_effects.parquet")
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)
    definitions = pd.read_parquet(ANALYSIS / "current_functional_site_definitions.parquet")

    # positions belonging to each tested apo/holo (receptor, segment)
    seg_positions = seg_map[seg_map["structural_segment"] != ""].groupby(
        ["uniprot_id", "structural_segment"]
    )["canonical_uniprot_position"].apply(set).to_dict()

    site_positions = definitions.groupby(["uniprot_id", "functional_site"])["canonical_uniprot_position"].apply(set).to_dict()

    rows = []
    for _, r in apo_holo.iterrows():
        uid, seg = r["uniprot_id"], r["structural_segment"]
        seg_pos = seg_positions.get((uid, seg), set())
        overlaps = []
        for (uid2, site), pos_set in site_positions.items():
            if uid2 != uid:
                continue
            overlap = seg_pos & pos_set
            if overlap:
                overlaps.append((site, len(overlap), len(seg_pos)))
        if overlaps:
            for site, n_overlap, n_seg in overlaps:
                rows.append({
                    "uniprot_id": uid, "nr_code": r["nr_code"], "apo_holo_segment": seg,
                    "apo_holo_delta_median_z": r["delta_median_z"], "apo_holo_raw_p": r["raw_p"],
                    "apo_holo_q_within_receptor": r["q_within_receptor"],
                    "apo_holo_significant_within_receptor": r["significant_within_receptor"],
                    "functional_site_overlap": True, "site_name": site,
                    "n_overlapping_positions": n_overlap, "n_segment_positions": n_seg,
                    "overlap_fraction_of_segment": n_overlap / n_seg if n_seg else None,
                })
        else:
            rows.append({
                "uniprot_id": uid, "nr_code": r["nr_code"], "apo_holo_segment": seg,
                "apo_holo_delta_median_z": r["delta_median_z"], "apo_holo_raw_p": r["raw_p"],
                "apo_holo_q_within_receptor": r["q_within_receptor"],
                "apo_holo_significant_within_receptor": r["significant_within_receptor"],
                "functional_site_overlap": False, "site_name": "",
                "n_overlapping_positions": 0, "n_segment_positions": len(seg_pos),
                "overlap_fraction_of_segment": 0.0,
            })

    overlap_df = pd.DataFrame(rows)
    out = TABLES / "apo_holo_functional_site_overlap.csv"
    overlap_df.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(PROJECT_ROOT)} ({len(overlap_df)} rows)")

    sig_with_overlap = overlap_df[overlap_df["apo_holo_significant_within_receptor"] & overlap_df["functional_site_overlap"]]
    print(f"\n{len(sig_with_overlap)} within-receptor-FDR-significant apo/holo segment results overlap "
          "a current functional-site definition:")
    if len(sig_with_overlap):
        print(sig_with_overlap[["nr_code", "apo_holo_segment", "site_name", "overlap_fraction_of_segment",
                                 "apo_holo_delta_median_z", "apo_holo_q_within_receptor"]].to_string(index=False))
    print("\nNo causal inference is drawn from this overlap -- it is a descriptive cross-reference only, "
          "and the apo/holo statistics themselves were not recomputed.")


if __name__ == "__main__":
    main()
