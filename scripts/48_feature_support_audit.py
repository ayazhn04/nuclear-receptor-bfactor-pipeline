"""
Part 1-2 of the segment-level PCA/clustering modeling stage (continuing from
commit 7d280bc on analysis/normalization-v1).

Part 1: freeze the modeling inputs. Record SHA256 hashes of the five
upstream files this whole modeling stage is built on; they are never
modified by any script in this stage (48-54). Any downstream test can
recompute and compare against data/analysis/modeling_frozen_input_hashes.json
to prove that.

Part 2: for every segment column realized in the v3 feature matrix, compute
an objective support audit (no hard-coded segment names) and derive
candidate feature sets at >=70/60/50/40% receptor support.

Primary biological unit = RECEPTOR (the 42 coordinate-bearing receptors).
The PDB-entry-level matrix is explicitly NOT used here -- it is reserved
for the Part 13 secondary/sensitivity analysis so that receptors with
hundreds of PDB entries (e.g. ESR1/GR) cannot dominate the superfamily
picture.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

# Part 1 -- frozen inputs (never modified by scripts 48-54)
FROZEN_FILES = [
    ANALYSIS / "receptor_segment_feature_matrix_median.parquet",
    ANALYSIS / "receptor_segment_feature_matrix_mean.parquet",
    ANALYSIS / "receptor_segment_coverage_matrix.parquet",
    ANALYSIS / "receptor_segment_flexibility_profiles.parquet",
    ANALYSIS / "nr_lbd_structural_segment_map.csv",
]
FROZEN_HASHES_PATH = ANALYSIS / "modeling_frozen_input_hashes.json"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def freeze_inputs() -> dict:
    hashes = {str(f.relative_to(PROJECT_ROOT)): _sha256(f) for f in FROZEN_FILES}
    FROZEN_HASHES_PATH.write_text(json.dumps(hashes, indent=2, sort_keys=True) + "\n")
    print(f"Part 1: froze {len(hashes)} modeling input file hashes -> "
          f"{FROZEN_HASHES_PATH.relative_to(PROJECT_ROOT)}")
    for k, v in hashes.items():
        print(f"  {k}: {v[:16]}...")
    return hashes


def get_coordinate_receptors() -> pd.DataFrame:
    rfs = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")
    coord = rfs[rfs["n_coordinate_candidates"] > 0][
        ["uniprot_id", "nr_code", "common_name", "group", "n_coordinate_candidates"]
    ].reset_index(drop=True)
    return coord


def build_feature_support_audit(coord_receptors: pd.DataFrame) -> pd.DataFrame:
    median = pd.read_parquet(ANALYSIS / "receptor_segment_feature_matrix_median.parquet")
    profiles = pd.read_parquet(ANALYSIS / "receptor_segment_flexibility_profiles.parquet")
    seg_map = pd.read_csv(ANALYSIS / "nr_lbd_structural_segment_map.csv", keep_default_na=False)

    n_coord = len(coord_receptors)
    coord_ids = coord_receptors["uniprot_id"].tolist()

    segment_cols = [c for c in median.columns if c != "uniprot_id"]
    # reindex the median matrix onto the full 42-receptor coordinate cohort so a
    # receptor entirely absent from the matrix (Q15466/SHP -- no usable CE
    # alignment to any anchor template) is counted as unsupported for every
    # segment, not silently dropped from the denominator.
    median_full = median.set_index("uniprot_id").reindex(coord_ids)

    named_seg_map = seg_map[seg_map["structural_segment"] != ""]

    # The core canonical H1-H12/beta_sheet vocabulary this whole project stage has been
    # using (v3 template curation, Part 9 acceptance table). A segment absent from the
    # feature-matrix columns entirely (currently only H8) is NOT silently invisible here
    # -- it is still audited explicitly, at 0 support, so "not fabricated" is verifiable
    # rather than merely "never mentioned."
    CANONICAL_CORE_SEGMENTS = [f"H{i}" for i in range(1, 13)] + ["beta_sheet"]
    all_audited_segments = list(segment_cols)
    for seg in CANONICAL_CORE_SEGMENTS:
        if seg not in all_audited_segments:
            all_audited_segments.append(seg)

    rows = []
    for seg in all_audited_segments:
        if seg in median_full.columns:
            vals = median_full[seg]
            n_supported = int(vals.notna().sum())
        else:
            # segment never realized as a feature-matrix column at all (currently only
            # H8) -- by construction that means zero receptors have a value for it;
            # audited explicitly rather than silently omitted.
            n_supported = 0
        fraction_supported = n_supported / n_coord

        prof_seg = profiles[profiles["structural_segment"] == seg]
        median_n_pdb = float(prof_seg["n_pdb_entries_with_segment"].median()) if len(prof_seg) else np.nan
        median_coverage = float(prof_seg["median_segment_coverage"].median()) if len(prof_seg) else np.nan

        smap_seg = named_seg_map[named_seg_map["structural_segment"] == seg]
        conf_counts = smap_seg["assignment_confidence"].value_counts()
        n_conf_rows = len(smap_seg)
        conf_high_frac = float(conf_counts.get("HIGH", 0) / n_conf_rows) if n_conf_rows else np.nan
        conf_moderate_frac = float(conf_counts.get("MODERATE", 0) / n_conf_rows) if n_conf_rows else np.nan
        conf_low_frac = float(conf_counts.get("LOW", 0) / n_conf_rows) if n_conf_rows else np.nan

        method_counts = smap_seg["assignment_method"].value_counts()
        n_own_template = int(method_counts.get("own_canonical_template", 0))
        n_alignment_transfer = int(method_counts.get("structural_alignment_transfer", 0))
        n_complementary_rescue = int(method_counts.get("complementary_structure_rescue", 0))

        if len(smap_seg):
            structural_class = smap_seg["structural_class"].mode().iloc[0]
        elif seg.startswith("H") and seg[1:].isdigit():
            # true by the H<n> naming convention itself, not fabricated biological
            # evidence -- every canonical H-numbered segment is a helix by definition,
            # independent of whether any receptor's structure was mapped to it
            structural_class = "HELIX"
        elif seg == "beta_sheet":
            structural_class = "SHEET"
        else:
            structural_class = "UNKNOWN"

        rows.append({
            "structural_segment": seg,
            "structural_class": structural_class,
            "n_coordinate_receptors_supported": n_supported,
            "n_coordinate_receptors_total": n_coord,
            "fraction_supported": fraction_supported,
            "median_n_pdb_entries_contributing": median_n_pdb,
            "median_segment_residue_coverage": median_coverage,
            "confidence_high_fraction": conf_high_frac,
            "confidence_moderate_fraction": conf_moderate_frac,
            "confidence_low_fraction": conf_low_frac,
            "n_positions_own_canonical_template": n_own_template,
            "n_positions_structural_alignment_transfer": n_alignment_transfer,
            "n_positions_complementary_structure_rescue": n_complementary_rescue,
        })

    audit = pd.DataFrame(rows).sort_values(
        ["fraction_supported", "structural_segment"], ascending=[False, True]
    ).reset_index(drop=True)

    out = TABLES / "modeling_feature_support.csv"
    audit.to_csv(out, index=False)
    print(f"Part 2: wrote {out.relative_to(PROJECT_ROOT)} ({len(audit)} segment columns audited)")
    return audit


def build_candidate_matrices(coord_receptors: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame:
    median = pd.read_parquet(ANALYSIS / "receptor_segment_feature_matrix_median.parquet")
    coord_ids = coord_receptors["uniprot_id"].tolist()
    n_coord = len(coord_ids)
    median_full = median.set_index("uniprot_id").reindex(coord_ids)

    thresholds = {"A": 0.70, "B": 0.60, "C": 0.50, "D": 0.40}
    rows = []
    for label, thresh in thresholds.items():
        features = audit.loc[audit["fraction_supported"] >= thresh, "structural_segment"].tolist()
        features = sorted(features)  # deterministic, never relies on audit row order
        sub = median_full[features] if features else median_full[[]]

        n_features = len(features)
        if n_features == 0:
            n_complete = 0
            total_missingness = np.nan
            per_receptor_missingness = pd.Series(dtype=float)
        else:
            observed = sub.notna()
            n_complete = int((observed.sum(axis=1) == n_features).sum())
            total_missingness = float(1.0 - observed.values.sum() / (n_coord * n_features))
            per_receptor_missingness = 1.0 - observed.sum(axis=1) / n_features

        def n_receptors_with_at_least(frac):
            if n_features == 0:
                return 0
            observed_frac = 1.0 - per_receptor_missingness
            return int((observed_frac >= frac - 1e-12).sum())

        rows.append({
            "candidate_set": label,
            "support_threshold": thresh,
            "n_features": n_features,
            "feature_names_json": json.dumps(features),
            "n_receptors_total": n_coord,
            "n_receptors_with_complete_data": n_complete,
            "total_matrix_missingness_fraction": total_missingness,
            "mean_per_receptor_missingness_fraction": (
                float(per_receptor_missingness.mean()) if n_features else np.nan
            ),
            "median_per_receptor_missingness_fraction": (
                float(per_receptor_missingness.median()) if n_features else np.nan
            ),
            "n_receptors_ge50pct_features_observed": n_receptors_with_at_least(0.50),
            "n_receptors_ge70pct_features_observed": n_receptors_with_at_least(0.70),
            "n_receptors_ge80pct_features_observed": n_receptors_with_at_least(0.80),
            "n_receptors_eq100pct_features_observed": n_receptors_with_at_least(1.00),
        })

    candidates = pd.DataFrame(rows)
    out = TABLES / "modeling_matrix_candidates.csv"
    candidates.to_csv(out, index=False)
    print(f"Part 2: wrote {out.relative_to(PROJECT_ROOT)}")
    print(candidates[["candidate_set", "support_threshold", "n_features",
                       "n_receptors_with_complete_data",
                       "total_matrix_missingness_fraction"]].to_string(index=False))
    return candidates


def main():
    freeze_inputs()
    coord_receptors = get_coordinate_receptors()
    print(f"\n{len(coord_receptors)} coordinate-bearing receptors (primary biological unit).")
    audit = build_feature_support_audit(coord_receptors)
    build_candidate_matrices(coord_receptors, audit)


if __name__ == "__main__":
    main()
