"""Normalize LBD C-alpha B-factors and summarize predefined functional sites.

The script uses only primary selected structures and canonical UniProt
positions from the final handoff. It creates descriptive site-level results;
it does not infer ligand pharmacology or perform apo-versus-holo tests.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "reports" / "tables"

STRUCTURES_CSV = PROCESSED_DIR / "structures.csv"
RESIDUE_MAP_PARQUET = PROCESSED_DIR / "final_lbd_residue_map.parquet"
LBD_REFERENCE_CSV = PROJECT_ROOT / "data" / "manifests" / "nr_lbd_reference.csv"
SITE_DEFINITIONS_CSV = TABLES_DIR / "functional_site_definitions.csv"

STRUCTURE_KEYS = ["uniprot_id", "pdb_id", "polymer_entity_id", "instance_id"]


def require_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def main() -> None:
    structures = pd.read_csv(STRUCTURES_CSV)
    residue_map = pd.read_parquet(RESIDUE_MAP_PARQUET)
    sites = pd.read_csv(SITE_DEFINITIONS_CSV)
    lbd_reference = pd.read_csv(LBD_REFERENCE_CSV)

    require_columns(structures, STRUCTURE_KEYS, "structures.csv")
    require_columns(
        residue_map,
        STRUCTURE_KEYS + ["canonical_uniprot_position", "ca_bfactor_usable", "selected_ca_b_iso_or_equiv"],
        "final_lbd_residue_map.parquet",
    )
    require_columns(
        sites,
        ["uniprot_id", "site_id", "uniprot_position", "expected_residue_name"],
        "functional-site definitions",
    )
    require_columns(lbd_reference, ["uniprot_id", "lbd_start", "lbd_end"], "LBD reference")

    if sites.duplicated(["uniprot_id", "site_id", "uniprot_position"]).any():
        raise ValueError("Functional-site definitions contain duplicate residue assignments.")

    site_bounds = sites.merge(lbd_reference[["uniprot_id", "lbd_start", "lbd_end"]], on="uniprot_id", how="left")
    outside_lbd = site_bounds.loc[
        (site_bounds["uniprot_position"] < site_bounds["lbd_start"])
        | (site_bounds["uniprot_position"] > site_bounds["lbd_end"])
        | site_bounds["lbd_start"].isna()
    ]
    if not outside_lbd.empty:
        raise ValueError("Functional-site positions outside the standardized LBD:\n" + outside_lbd.to_string(index=False))

    canonical_residues = residue_map[
        ["uniprot_id", "canonical_uniprot_position", "canonical_residue_name"]
    ].drop_duplicates()
    conflicting_residues = canonical_residues.duplicated(
        ["uniprot_id", "canonical_uniprot_position"], keep=False
    )
    if conflicting_residues.any():
        raise ValueError("Residue map contains conflicting canonical residue identities.")
    site_identities = sites.merge(
        canonical_residues,
        left_on=["uniprot_id", "uniprot_position"],
        right_on=["uniprot_id", "canonical_uniprot_position"],
        how="left",
    )
    identity_mismatches = site_identities.loc[
        site_identities["canonical_residue_name"].notna()
        & site_identities["canonical_residue_name"].ne(site_identities["expected_residue_name"])
    ]
    if not identity_mismatches.empty:
        raise ValueError(
            "Functional-site residue identities do not match the canonical sequence:\n"
            + identity_mismatches[
                ["uniprot_id", "site_id", "uniprot_position", "expected_residue_name", "canonical_residue_name"]
            ].to_string(index=False)
        )

    usable = residue_map.loc[
        residue_map["ca_bfactor_usable"].eq(True)
        & residue_map["selected_ca_b_iso_or_equiv"].notna()
    ].copy()
    usable["selected_ca_b_iso_or_equiv"] = pd.to_numeric(usable["selected_ca_b_iso_or_equiv"], errors="raise")

    # Normalize only within one selected receptor instance; raw B scales differ between entries.
    instance_group = usable.groupby(STRUCTURE_KEYS)["selected_ca_b_iso_or_equiv"]
    usable["structure_mean_b"] = instance_group.transform("mean")
    usable["structure_sd_b"] = instance_group.transform(lambda values: values.std(ddof=0))
    usable["b_factor_z"] = np.where(
        usable["structure_sd_b"] > 0,
        (usable["selected_ca_b_iso_or_equiv"] - usable["structure_mean_b"]) / usable["structure_sd_b"],
        np.nan,
    )
    usable = usable.dropna(subset=["b_factor_z"])

    normalization_qc = (
        usable.groupby(STRUCTURE_KEYS, as_index=False)
        .agg(
            n_usable_positions=("b_factor_z", "size"),
            raw_b_mean=("structure_mean_b", "first"),
            raw_b_sd=("structure_sd_b", "first"),
            z_mean=("b_factor_z", "mean"),
            z_sd=("b_factor_z", lambda values: values.std(ddof=0)),
        )
    )

    # Multiple selected entities from one PDB entry are reduced before a PDB becomes an observation.
    entry_position = (
        usable.groupby(["uniprot_id", "pdb_id", "canonical_uniprot_position"], as_index=False)
        .agg(b_factor_z=("b_factor_z", "median"), n_selected_instances=("instance_id", "nunique"))
    )
    receptor_entries = structures[["uniprot_id", "pdb_id"]].drop_duplicates()
    site_entry_grid = sites.merge(receptor_entries, on="uniprot_id", how="inner")
    site_entry = site_entry_grid.merge(
        entry_position,
        left_on=["uniprot_id", "pdb_id", "uniprot_position"],
        right_on=["uniprot_id", "pdb_id", "canonical_uniprot_position"],
        how="left",
    )
    site_entry_scores = (
        site_entry.groupby(["uniprot_id", "nr_code", "common_name", "site_id", "site_name", "pdb_id"], as_index=False)
        .agg(
            site_mean_z=("b_factor_z", "mean"),
            n_site_residues=("uniprot_position", "size"),
            n_available_residues=("b_factor_z", "count"),
        )
    )
    site_entry_scores["site_coverage"] = (
        site_entry_scores["n_available_residues"] / site_entry_scores["n_site_residues"]
    )
    # A site score must represent the same predefined residue set in every PDB entry.
    site_entry_scores.loc[site_entry_scores["site_coverage"].lt(1), "site_mean_z"] = np.nan

    receptor_site_scores = (
        site_entry_scores.groupby(
            ["uniprot_id", "nr_code", "common_name", "site_id", "site_name"], as_index=False
        )
        .agg(
            receptor_median_site_z=("site_mean_z", "median"),
            n_pdb_entries=("pdb_id", "nunique"),
            n_entries_with_site_data=("site_mean_z", "count"),
            median_site_coverage=("site_coverage", "median"),
        )
    )
    site_catalog = sites[["uniprot_id", "nr_code", "common_name", "site_id", "site_name"]].drop_duplicates()
    receptor_site_scores = site_catalog.merge(
        receptor_site_scores,
        on=["uniprot_id", "nr_code", "common_name", "site_id", "site_name"],
        how="left",
    )
    receptor_site_scores[["n_pdb_entries", "n_entries_with_site_data"]] = receptor_site_scores[
        ["n_pdb_entries", "n_entries_with_site_data"]
    ].fillna(0).astype(int)

    normalization_qc.to_csv(TABLES_DIR / "functional_site_normalization_qc.csv", index=False)
    site_entry_scores.to_csv(TABLES_DIR / "functional_site_entry_scores.csv", index=False)
    receptor_site_scores.to_csv(TABLES_DIR / "functional_site_receptor_scores.csv", index=False)

    print(f"Normalized {normalization_qc.shape[0]} selected receptor instances.")
    print(f"Wrote {len(site_entry_scores)} PDB-entry site scores.")
    print(f"Wrote {len(receptor_site_scores)} receptor-level site summaries.")


if __name__ == "__main__":
    main()
