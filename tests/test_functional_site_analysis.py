from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
TABLES_DIR = PROJECT_ROOT / "reports" / "tables"


def test_functional_site_residues_match_canonical_sequence():
    sites = pd.read_csv(TABLES_DIR / "functional_site_definitions.csv")
    residue_map = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "final_lbd_residue_map_primary_qc_subset.parquet")
    canonical = residue_map[
        ["uniprot_id", "canonical_uniprot_position", "canonical_residue_name"]
    ].drop_duplicates()
    checked = sites.merge(
        canonical,
        left_on=["uniprot_id", "uniprot_position"],
        right_on=["uniprot_id", "canonical_uniprot_position"],
        how="inner",
    )

    assert len(checked) == len(sites[sites["uniprot_id"] != "P11473"])
    assert checked["expected_residue_name"].eq(checked["canonical_residue_name"]).all()


def test_pparg_sites_use_canonical_pparg2_numbering():
    sites = pd.read_csv(TABLES_DIR / "functional_site_definitions.csv")
    pparg = sites[sites["uniprot_id"] == "P37231"]

    assert set(pparg.loc[pparg["site_id"] == "ligand_pocket", "uniprot_position"]) == {
        317,
        351,
        477,
        501,
    }
    assert set(pparg.loc[pparg["site_id"] == "af2_coactivator_groove", "uniprot_position"]) == {
        329,
        499,
    }


def test_only_complete_sites_receive_scores():
    entries = pd.read_csv(TABLES_DIR / "functional_site_entry_scores.csv")

    assert entries.loc[entries["site_mean_z"].notna(), "site_coverage"].eq(1).all()
    assert entries.loc[entries["site_coverage"].lt(1), "site_mean_z"].isna().all()
