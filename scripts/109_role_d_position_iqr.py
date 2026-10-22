"""Role D v2: descriptive between-PDB variability at functional residues.

The observation unit is the established PDB-entry median of instance-normalized
Ca B-factor Z.  Position IQR therefore describes heterogeneous relative
crystallographic disorder across deposited structures; it is not solution
dynamics, an intrinsic-flexibility score, or a hypothesis test.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "data/rebuild_v2/analysis"
TABLES = ROOT / "reports/role_d/contextual_iqr"
PROFILE = ANALYSIS / "receptor_position_professor_profile_v2.csv"
CURRENT_SITES = ANALYSIS / "current_functional_site_definitions.parquet"
MIN_PDB = 5


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    profile = pd.read_csv(PROFILE, keep_default_na=False)
    numeric = ["iqr_z", "q25_z", "q75_z", "median_z", "sd_z", "n_pdb_entries_observed_at_position"]
    for col in numeric:
        profile[col] = pd.to_numeric(profile[col], errors="coerce")
    profile["review_required_flag"] = profile["review_required_flag"].astype(str).str.strip().str.casefold().eq("true")
    profile["eligible_primary"] = profile.n_pdb_entries_observed_at_position.ge(MIN_PDB) & profile.iqr_z.notna()
    baseline = profile[profile.eligible_primary].groupby("uniprot_id").iqr_z.median()
    profile["receptor_median_position_iqr_z"] = profile.uniprot_id.map(baseline)
    profile["excess_iqr_z"] = np.where(profile.eligible_primary, profile.iqr_z - profile.receptor_median_position_iqr_z, np.nan)
    profile["iqr_eligibility"] = np.where(
        profile.n_pdb_entries_observed_at_position.eq(0), "NO_PDB_OBSERVATION",
        np.where(profile.eligible_primary, "ELIGIBLE_PRIMARY", f"INSUFFICIENT_PDB_SUPPORT_LT_{MIN_PDB}"),
    )
    profile["context_stratum_display"] = np.where(
        profile.canonical_elements.ne(""), profile.canonical_elements,
        np.where(profile.position_status.eq("INTER_HELIX_UNASSIGNED_DERIVED"), "INTER_HELIX_DERIVED", "NO_PROFESSOR_MAPPING"),
    )

    # Ligand-pocket residues retain the current geometry-derived v2 definition.
    defs = pd.read_parquet(CURRENT_SITES)
    pocket = defs[(defs.functional_site == "LIGAND_BINDING_POCKET") & (defs.definition_level == "RESIDUE")].copy()
    # AF2/H12 is rebuilt directly from corrected professor membership, avoiding
    # stale phase-2 H12 endpoints outside the corrected analysis span.
    h12 = profile[profile.canonical_elements.str.split("|").map(lambda labels: "H12" in labels)][
        ["uniprot_id", "canonical_uniprot_position"]
    ].copy()
    h12["functional_site"] = "AF2_H12"
    h12["definition_level"] = "RESIDUE"
    h12["evidence_type"] = "professor_canonical_h12_v2"
    h12["source_or_method"] = "corrected professor canonical H12 membership"
    selected = pd.concat([pocket, h12], ignore_index=True, sort=False).drop(columns="nr_code", errors="ignore")
    selected = selected.rename(columns={"assignment_status": "site_assignment_status"})
    selected.canonical_uniprot_position = selected.canonical_uniprot_position.astype(int)
    profile = profile.rename(columns={"assignment_status": "professor_assignment_status"})
    positions = selected.merge(profile, on=["uniprot_id", "canonical_uniprot_position"], validate="many_to_one")
    assert len(positions) == len(selected), "functional residue absent from corrected v2 profile"
    assert not positions.duplicated(["uniprot_id", "functional_site", "canonical_uniprot_position"]).any()
    finite = positions.iqr_z.notna()
    assert np.allclose(positions.loc[finite, "iqr_z"], positions.loc[finite, "q75_z"] - positions.loc[finite, "q25_z"])
    positions["functional_position_label"] = (
        positions.nr_code + " " + positions.canonical_residue_name + positions.canonical_uniprot_position.astype(str)
    )
    positions = positions.sort_values(["functional_site", "uniprot_id", "canonical_uniprot_position"]).reset_index(drop=True)

    keys = ["uniprot_id", "nr_code", "common_name", "functional_site"]
    positions["primary_iqr_z"] = positions["iqr_z"].where(positions["eligible_primary"])
    summary = positions.groupby(keys, sort=True).agg(
        n_defined_positions=("canonical_uniprot_position", "size"),
        n_positions_observed=("has_observation", "sum"),
        n_positions_eligible_primary=("eligible_primary", "sum"),
        median_position_iqr_z=("primary_iqr_z", "median"),
        median_position_iqr_z_all_observed=("iqr_z", "median"),
        median_excess_iqr_z=("excess_iqr_z", "median"),
    ).reset_index()
    summary["interpretation"] = "DESCRIPTIVE_ONLY_BETWEEN_PDB_POSITION_VARIABILITY"
    return positions, summary


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    positions, summary = build()
    positions.to_csv(TABLES / "role_d_functional_position_iqr.csv", index=False, lineterminator="\n")
    summary.to_csv(TABLES / "role_d_functional_site_iqr_summary.csv", index=False, lineterminator="\n")
    print(f"Wrote {len(positions)} functional residues and {len(summary)} receptor-site summaries.")


if __name__ == "__main__":
    main()
