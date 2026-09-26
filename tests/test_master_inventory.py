"""
Master-inventory revision tests: QC criteria are metadata/analysis flags,
never deletion criteria. Offline; uses only already-generated files.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

MASTER_CSV = PROCESSED_DIR / "structures.csv"
SUBSET_CSV = PROCESSED_DIR / "structures_primary_qc_subset.csv"
SUBSET_PARQUET = PROCESSED_DIR / "final_lbd_residue_map_primary_qc_subset.parquet"
RESIDUE_ALL = PROCESSED_DIR / "lbd_residue_bfactors_all.parquet"
OBS_PARQUET = PROJECT_ROOT / "data" / "interim" / "qc" / "stage3b_lbd_ca_observations.parquet"
POOL_CSV = MANIFESTS_DIR / "stage3b_coordinate_pool.csv"
INSTANCE_QC_CSV = MANIFESTS_DIR / "rcsb_lbd_coordinate_instance_qc.csv"

# SHA256 of the strict subset files exactly as committed in a42e875.
SUBSET_CSV_SHA256 = "5222cfd462cb558c421f2a6780952612e8b7d4bedd6a7930de287afeb90ed58b"
SUBSET_PARQUET_SHA256 = "ce118a4d049d5b4123daf460bffaca965493cb1154071a25bfe15afccadfacbc"


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@pytest.fixture(scope="module")
def master():
    return pd.read_csv(MASTER_CSV, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def residues():
    return pd.read_parquet(RESIDUE_ALL)


def test_master_has_exactly_2072_rows(master):
    assert len(master) == 2072


def test_every_stage2_candidate_appears_exactly_once(master):
    pool = pd.read_csv(POOL_CSV, dtype={"pdb_id": str})
    key = ["uniprot_id", "pdb_id", "polymer_entity_id"]
    assert not master.duplicated(subset=key).any()
    assert set(map(tuple, master[key].values)) == set(map(tuple, pool[key].values))


def test_required_master_columns_present(master):
    required = """uniprot_id nr_code common_name group pdb_id polymer_entity_id entity_id experimental_method
    has_stage3b_coordinate_data method_is_xray mapping_overlaps_lbd standardized_lbd_start standardized_lbd_end
    standardized_lbd_length construct_lbd_mapping_coverage positive_occupancy_lbd_coverage bfactor_usable_lbd_coverage
    coverage_ge_080 coverage_ge_090 coverage_ge_095 coverage_ge_099 coverage_eq_100 resolution
    resolution_project_range_status r_free r_free_status identity_class fusion_or_chimera taxonomy_audit_status
    previous_primary_qc_member previous_primary_exclusion_reasons_json representative_instance_id label_asym_id
    auth_asym_id ligand_state""".split()
    assert set(required) <= set(master.columns)


def test_vdr_represented_in_master_and_residues(master, residues):
    v = master[master["uniprot_id"] == "P11473"]
    assert len(v) == 52 and v["pdb_id"].nunique() == 52
    assert int(v["has_stage3b_coordinate_data"].sum()) == 48
    vr = residues[residues["uniprot_id"] == "P11473"]
    assert len(vr) == 48 * 297


def test_no_structure_removed_solely_for_coverage(master):
    """Every row with coordinate data is present regardless of coverage;
    rows failing 0.90 are retained with the flag False, not dropped."""
    with_coord = master[master["has_stage3b_coordinate_data"]]
    assert len(with_coord) == 1840
    assert (~with_coord["coverage_ge_090"].astype(bool)).sum() > 0
    low = with_coord[~with_coord["coverage_ge_090"].astype(bool)]
    assert low["positive_occupancy_lbd_coverage"].notna().all()


def test_rows_without_coordinate_data_have_explicit_missing_and_status(master):
    no = master[~master["has_stage3b_coordinate_data"]]
    assert len(no) == 232
    assert no["positive_occupancy_lbd_coverage"].isna().all()
    assert no["coverage_ge_090"].isna().all()
    assert (no["coordinate_data_status"] == "NO_STAGE3B_COORDINATE_DATA").all()


def test_coverage_flags_consistent_with_coverage(master):
    w = master[master["has_stage3b_coordinate_data"]]
    c = w["positive_occupancy_lbd_coverage"].astype(float)
    assert (w["coverage_ge_090"].astype(bool) == (c >= 0.90)).all()
    assert (w["coverage_ge_080"].astype(bool) == (c >= 0.80)).all()
    assert (w["coverage_eq_100"].astype(bool) == (c >= 1.0)).all()


def test_nested_view_counts(master):
    assert len(master) == 2072
    assert int(master["has_stage3b_coordinate_data"].sum()) == 1840
    assert int(master["previous_primary_qc_member"].sum()) == 1382


def test_old_subset_preserved_exactly():
    assert _sha(SUBSET_CSV) == SUBSET_CSV_SHA256
    assert _sha(SUBSET_PARQUET) == SUBSET_PARQUET_SHA256
    old = subprocess.run(
        ["git", "show", "a42e875:data/processed/structures.csv"], cwd=PROJECT_ROOT, capture_output=True
    ).stdout
    assert hashlib.sha256(old).hexdigest() == SUBSET_CSV_SHA256


def test_previous_membership_matches_preserved_subset(master):
    sub = pd.read_csv(SUBSET_CSV, dtype={"pdb_id": str})
    key = ["uniprot_id", "pdb_id", "polymer_entity_id"]
    m = master[master["previous_primary_qc_member"]]
    assert set(map(tuple, m[key].values)) == set(map(tuple, sub[key].values))
    assert len(sub) == 1382


def test_all_stage3b_residue_observations_represented(residues):
    if not OBS_PARQUET.exists():
        pytest.skip("development-only Stage 3B interim parquet is not tracked in the clean release repository")
    obs = pd.read_parquet(OBS_PARQUET)
    assert len(residues) == len(obs) == 771542
    assert residues["instance_id"].nunique() == pd.read_csv(INSTANCE_QC_CSV)["instance_id"].nunique() == 3159
    k = ["instance_id", "lbd_relative_position"]
    assert not residues.duplicated(subset=k).any()
    assert set(map(tuple, residues[k].values)) == set(map(tuple, obs[k].values))


def test_residues_not_filtered_by_coverage(residues, master):
    low_ids = set(master.loc[master["coverage_ge_090"].astype("boolean") == False, "representative_instance_id"])  # noqa: E712
    assert len(low_ids) > 0
    assert low_ids <= set(residues["instance_id"])


def test_positive_occupancy_semantics_unchanged(residues):
    if not OBS_PARQUET.exists():
        pytest.skip("development-only Stage 3B interim parquet is not tracked in the clean release repository")
    obs = pd.read_parquet(OBS_PARQUET)
    for c in ["ca_record_present", "ca_positive_occupancy", "ca_bfactor_usable", "zero_occupancy_only_ca"]:
        assert (residues[c].values == obs[c].values).all()
    # observed => record present; zero-occupancy-only is never observed
    assert not (residues["ca_positive_occupancy"] & ~residues["ca_record_present"]).any()
    assert not (residues["zero_occupancy_only_ca"] & residues["ca_positive_occupancy"]).any()
    assert not (residues["ca_bfactor_usable"] & ~residues["ca_positive_occupancy"]).any()


def test_missing_lbd_residues_explicit_and_selected_altloc_valid(residues):
    missing = residues[~residues["ca_positive_occupancy"]]
    assert len(missing) > 0
    assert missing["selected_ca_occupancy"].isna().all()
    usable = residues[residues["selected_ca_occupancy"].notna()]
    assert (usable["selected_ca_occupancy"] > 0).all()
    assert usable["ca_bfactor_usable"].all()


def test_no_normalization_columns_in_residue_table(residues):
    bad = {"z_score", "normalized_b", "b_factor_zscore", "cluster_label"}
    assert bad.isdisjoint({c.lower() for c in residues.columns})


def test_frozen_stage_summaries_unchanged():
    r = subprocess.run(
        ["git", "diff", "--quiet", "a42e875", "--",
         "data/manifests/stage2_release_summary.json", "data/manifests/stage3a_release_summary.json",
         "data/manifests/stage3b_release_summary.json", "data/manifests/nr_lbd_reference.csv",
         "data/manifests/rcsb_lbd_coordinate_instance_qc.csv", "config/nr_metadata.py"],
        cwd=PROJECT_ROOT,
    )
    assert r.returncode == 0
    # The clean GitHub release intentionally omits development-history tags.
    # Frozen-artifact immutability is anchored to the clean baseline commit.
    assert subprocess.run(
        ["git", "cat-file", "-e", "a42e875^{commit}"],
        cwd=PROJECT_ROOT,
        capture_output=True,
    ).returncode == 0


# ---------------------------------------------------------------------------
# Ligand annotation extended to all coordinate candidates (QC = metadata only)
# ---------------------------------------------------------------------------

LIGAND_CSV = MANIFESTS_DIR / "final_ligand_annotations.csv"
RESIDUE_ALL_SHA256 = "28cc0f7feff59d27cf6b724becad23391068c6ed74181db3a6b9923750f81a27"


def test_all_1840_coordinate_candidates_have_explicit_ligand_assessment(master):
    lig = pd.read_csv(LIGAND_CSV, dtype={"pdb_id": str})
    key = ["uniprot_id", "pdb_id", "polymer_entity_id"]
    assert len(lig) == 1840 and not lig.duplicated(subset=key).any()
    coord = master[master["has_stage3b_coordinate_data"]]
    assert set(map(tuple, coord[key].values)) == set(map(tuple, lig[key].values))
    assert lig["ligand_state"].isin(["APO", "HOLO", "AMBIGUOUS"]).all()
    assert coord["ligand_state"].isin(["APO", "HOLO", "AMBIGUOUS"]).all()
    assert coord["ligand_annotation_status"].isin(["RESOLVED", "REVIEW_NEEDED"]).all()


def test_non_coordinate_candidates_retained_not_labelled_apo(master):
    no = master[~master["has_stage3b_coordinate_data"]]
    assert len(no) == 232
    assert (no["ligand_state"] == "NOT_ASSESSABLE_NO_COORDINATE_DATA").all()
    assert (no["ligand_annotation_status"] == "NOT_ASSESSABLE_NO_COORDINATE_DATA").all()


def test_ligand_state_never_removes_candidates(master):
    assert len(master) == 2072
    assert master.groupby("ligand_state").size().sum() == 2072


def test_low_coverage_and_nonstrict_candidates_are_annotated(master):
    coord = master[master["has_stage3b_coordinate_data"]]
    low = coord[~coord["previous_primary_qc_member"]]
    assert len(low) == 458
    assert low["ligand_state"].isin(["APO", "HOLO", "AMBIGUOUS"]).all()


def test_strict_subset_ligand_assessments_unchanged_in_master(master):
    sub = pd.read_csv(SUBSET_CSV, dtype={"pdb_id": str})
    m = master.merge(sub[["uniprot_id", "pdb_id", "polymer_entity_id", "ligand_state"]],
                     on=["uniprot_id", "pdb_id", "polymer_entity_id"], suffixes=("", "_subset"))
    assert len(m) == 1382
    assert (m["ligand_state"] == m["ligand_state_subset"]).all()


def test_vdr_ligand_assessed_regardless_of_coverage(master):
    v = master[(master["uniprot_id"] == "P11473") & master["has_stage3b_coordinate_data"]]
    assert len(v) == 48
    assert (~v["previous_primary_qc_member"]).all()
    assert v["ligand_state"].isin(["APO", "HOLO", "AMBIGUOUS"]).all()


def test_holo_has_functional_ligand_and_apo_none_all_candidates():
    lig = pd.read_csv(LIGAND_CSV, dtype={"pdb_id": str})
    n = lig["functional_ligand_component_ids_json"].apply(lambda j: len(json.loads(j)))
    assert (n[lig["ligand_state"] == "HOLO"] >= 1).all()
    assert (n[lig["ligand_state"] == "APO"] == 0).all()
    assert (lig.loc[lig["ligand_state"] == "AMBIGUOUS", "ligand_annotation_status"] == "REVIEW_NEEDED").all()


def test_residue_parquet_unchanged_by_ligand_extension():
    assert _sha(RESIDUE_ALL) == RESIDUE_ALL_SHA256
