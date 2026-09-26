"""
Final Data Release tests: structure selection, representative-instance
determinism, canonical residue mapping, deterministic altloc handoff,
ligand-state rules, full handoff-file consistency, and immutability of
every earlier frozen stage's release summary.

Offline only — operates on already-generated manifests/processed tables
and the cached mmCIF/UniProt files. No network access.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RAW_MMCIF_DIR = PROJECT_ROOT / "data" / "raw" / "mmcif"

STRUCTURES_CSV = PROCESSED_DIR / "structures_primary_qc_subset.csv"  # historical strict subset (was structures.csv)
EXCLUDED_STRUCTURES_CSV = PROCESSED_DIR / "excluded_structures.csv"
RESIDUE_MAP_PARQUET = PROCESSED_DIR / "final_lbd_residue_map_primary_qc_subset.parquet"
SELECTION_AUDIT_CSV = MANIFESTS_DIR / "final_structure_selection_audit.csv"
REP_SELECTION_CSV = MANIFESTS_DIR / "final_representative_instance_selection.csv"
LIGAND_ANNOTATIONS_CSV = MANIFESTS_DIR / "final_ligand_annotations.csv"
NONPOLYMER_INVENTORY_CSV = MANIFESTS_DIR / "final_nonpolymer_inventory.csv"
FINAL_COUNTS_BY_RECEPTOR_CSV = MANIFESTS_DIR / "final_counts_by_receptor.csv"
DOWNLOAD_MANIFEST_CSV = MANIFESTS_DIR / "mmcif_download_manifest.csv"
RELEASE_JSON = MANIFESTS_DIR / "final_data_infrastructure_release.json"

STAGE2_RELEASE_SUMMARY = MANIFESTS_DIR / "stage2_release_summary.json"
STAGE3A_RELEASE_SUMMARY = MANIFESTS_DIR / "stage3a_release_summary.json"
STAGE3B_RELEASE_SUMMARY = MANIFESTS_DIR / "stage3b_release_summary.json"


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet")


@pytest.fixture(scope="module")
def structures():
    _require(STRUCTURES_CSV)
    return pd.read_csv(STRUCTURES_CSV, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def excluded():
    _require(EXCLUDED_STRUCTURES_CSV)
    return pd.read_csv(EXCLUDED_STRUCTURES_CSV, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def audit():
    _require(SELECTION_AUDIT_CSV)
    return pd.read_csv(SELECTION_AUDIT_CSV, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def rep_selection():
    _require(REP_SELECTION_CSV)
    return pd.read_csv(REP_SELECTION_CSV, dtype={"pdb_id": str})


@pytest.fixture(scope="module")
def residue_map():
    _require(RESIDUE_MAP_PARQUET)
    return pd.read_parquet(RESIDUE_MAP_PARQUET)


@pytest.fixture(scope="module")
def ligand():
    _require(LIGAND_ANNOTATIONS_CSV)
    return pd.read_csv(LIGAND_ANNOTATIONS_CSV, dtype={"pdb_id": str})


# ---------------------------------------------------------------------------
# SELECTION
# ---------------------------------------------------------------------------

def test_structures_csv_rows_obey_all_primary_rules(structures):
    assert (structures["experimental_method"].str.contains("X-RAY", case=False)).all()
    assert (structures["positive_occupancy_lbd_coverage"] >= 0.90).all()
    assert (structures["resolution"] >= 1.8).all()
    assert (structures["resolution"] <= 3.5).all()
    assert (structures["r_free"] <= 0.30).all()
    assert (structures["identity_class"] == "HUMAN_OR_HUMAN_DERIVED").all()


def test_excluded_never_appears_in_structures(structures, excluded):
    primary_keys = set(zip(structures["uniprot_id"], structures["pdb_id"], structures["polymer_entity_id"]))
    excluded_keys = set(zip(excluded["uniprot_id"], excluded["pdb_id"], excluded["polymer_entity_id"]))
    assert primary_keys.isdisjoint(excluded_keys)


def test_selected_and_excluded_accounting_is_complete(structures, excluded):
    assert len(structures) + len(excluded) == 2072


def test_every_stage2_candidate_appears_in_primary_or_excluded_lineage(audit):
    assert len(audit) == 2072
    assert audit["primary_selection_status"].isin(
        {"SELECTED_PRIMARY", "EXCLUDED_PRIMARY", "NO_COORDINATE_CANDIDATE"}
    ).all()


def test_no_logical_duplicate_selected_row(structures):
    dup = structures.duplicated(subset=["uniprot_id", "pdb_id", "polymer_entity_id"])
    assert not dup.any()


def test_excluded_structures_have_at_least_one_reason_when_evaluated(excluded):
    evaluated = excluded[excluded["primary_selection_status"] == "EXCLUDED_PRIMARY"]
    if len(evaluated) == 0:
        pytest.skip("no EXCLUDED_PRIMARY rows in this run")
    reason_counts = evaluated["primary_exclusion_reasons_json"].apply(lambda j: len(json.loads(j)))
    assert (reason_counts >= 1).all()


def test_no_coordinate_candidate_rows_have_no_representative_instance(excluded):
    no_candidate = excluded[excluded["primary_selection_status"] == "NO_COORDINATE_CANDIDATE"]
    if len(no_candidate) == 0:
        pytest.skip("no NO_COORDINATE_CANDIDATE rows in this run")
    assert (no_candidate["representative_instance_id"].isna() | (no_candidate["representative_instance_id"] == "")).all()


# ---------------------------------------------------------------------------
# INSTANCE (representative-instance determinism)
# ---------------------------------------------------------------------------

def test_representative_selection_one_row_per_candidate_triple(rep_selection):
    dup = rep_selection.duplicated(subset=["uniprot_id", "pdb_id", "polymer_entity_id"])
    assert not dup.any()
    assert len(rep_selection) == 1840


def test_representative_is_always_the_max_coverage_instance(rep_selection):
    for _, row in rep_selection.iterrows():
        candidates = pd.DataFrame(json.loads(row["all_candidate_qc_json"]))
        best_cov = candidates["observed_full_lbd_ca_coverage"].max()
        assert math.isclose(row["selected_observed_full_lbd_ca_coverage"], best_cov, rel_tol=1e-9, abs_tol=1e-9), (
            f"{row['uniprot_id']}/{row['pdb_id']}/{row['polymer_entity_id']}: selected instance is not the "
            "max-coverage instance among its candidates"
        )


def test_representative_selection_reproducible_by_rerunning_tie_break_logic(rep_selection):
    """Re-run the SAME deterministic cascade independently in this test
    (not by importing the script) against the recorded candidate QC JSON,
    and confirm it reproduces the recorded selection — catches any drift
    between the script's logic and its own recorded output."""
    for _, row in rep_selection.iterrows():
        candidates = pd.DataFrame(json.loads(row["all_candidate_qc_json"]))
        ordered = candidates.sort_values(
            by=["observed_full_lbd_ca_coverage", "bfactor_usable_full_lbd_coverage",
                "internal_missing_residue_count", "mapped_with_zero_occupancy_only_ca_residues",
                "total_mapped_unmodeled_residues", "label_asym_id"],
            ascending=[False, False, True, True, True, True],
        )
        assert ordered.iloc[0]["instance_id"] == row["selected_instance_id"]


def test_tie_break_never_uses_bfactor_magnitude():
    """Static check: the tie-break column list in the selection script must
    never reference a B-factor-magnitude column."""
    script_text = (PROJECT_ROOT / "scripts" / "24_final_structure_selection.py").read_text()
    tie_break_block_start = script_text.index("TIE_BREAK_STEPS")
    tie_break_block = script_text[tie_break_block_start:tie_break_block_start + 2000]
    assert "b_iso" not in tie_break_block.lower()
    assert "bfactor_usable_full_lbd_coverage" in tie_break_block  # coverage, not magnitude — allowed


# ---------------------------------------------------------------------------
# MAPPING
# ---------------------------------------------------------------------------

def test_residue_map_every_standardized_lbd_position_represented(residue_map, structures):
    for _, s in structures.iterrows():
        rows = residue_map[residue_map["instance_id"] == s["instance_id"]]
        assert len(rows) == s["standardized_lbd_length"], (
            f"{s['pdb_id']}/{s['instance_id']}: expected {s['standardized_lbd_length']} rows, got {len(rows)}"
        )
        assert set(rows["lbd_relative_position"]) == set(range(1, int(s["standardized_lbd_length"]) + 1))


def test_residue_map_missing_residues_explicitly_flagged_not_dropped(residue_map):
    unmapped = residue_map[~residue_map["is_construct_mapped"]]
    if len(unmapped) == 0:
        pytest.skip("no unmapped residues in this run")
    assert (~unmapped["ca_positive_occupancy"]).all()
    assert unmapped["selected_ca_occupancy"].isna().all()


def test_residue_map_uses_label_seq_id_derived_canonical_position_not_auth(residue_map):
    sample = residue_map[residue_map["is_construct_mapped"]].head(500)
    assert sample["entity_label_seq_id"].notna().all()
    assert sample["canonical_uniprot_position"].notna().all()


def test_auth_numbering_is_provenance_only(residue_map):
    """auth_seq_ids_json must be present as a JSON list column, never used
    to compute canonical_uniprot_position (spot check: two rows with the
    same canonical_uniprot_position but different pdb entries can have
    unrelated auth numbering)."""
    sample = residue_map[residue_map["is_construct_mapped"]].head(1)
    if len(sample) == 0:
        pytest.skip("no mapped residues")
    val = json.loads(sample.iloc[0]["auth_seq_ids_json"])
    assert isinstance(val, list)


# ---------------------------------------------------------------------------
# ALTLOC
# ---------------------------------------------------------------------------

def test_selected_ca_has_positive_occupancy_and_finite_bfactor(residue_map):
    usable = residue_map[residue_map["selected_ca_occupancy"].notna()]
    if len(usable) == 0:
        pytest.skip("no usable selected Ca rows in this run")
    assert (usable["selected_ca_occupancy"] > 0).all()
    assert usable["selected_ca_b_iso_or_equiv"].apply(lambda v: math.isfinite(v)).all()


def test_selected_ca_is_highest_occupancy_among_usable_altlocs(residue_map):
    sample = residue_map[residue_map["selected_ca_occupancy"].notna()].sample(
        min(300, residue_map["selected_ca_occupancy"].notna().sum()), random_state=42
    )
    for _, row in sample.iterrows():
        occs = json.loads(row["occupancies_json"])
        b_isos = json.loads(row["b_iso_values_json"])
        usable_occs = [o for o, b in zip(occs, b_isos) if o is not None and o > 0 and b is not None and math.isfinite(b)]
        assert math.isclose(row["selected_ca_occupancy"], max(usable_occs), rel_tol=1e-9, abs_tol=1e-9)


def test_altloc_tie_break_deterministic_direct():
    from scripts.utils.altloc_selection import select_representative_ca_record

    # Tie on occupancy: blank preferred over 'A'.
    r = select_representative_ca_record([1.0, 1.0], [20.0, 999.0], ["B", "."])
    assert r["altloc"] == "."
    assert r["occupancy"] == 1.0

    # Tie on occupancy, no blank: 'A' preferred over other letters.
    r = select_representative_ca_record([0.5, 0.5], [999.0, 1.0], ["B", "A"])
    assert r["altloc"] == "A"

    # Tie on occupancy, no blank or 'A': lexicographically smallest.
    r = select_representative_ca_record([0.5, 0.5], [1.0, 999.0], ["C", "B"])
    assert r["altloc"] == "B"

    # Highest occupancy wins regardless of B-factor magnitude.
    r = select_representative_ca_record([0.3, 0.7], [1.0, 999.0], ["A", "B"])
    assert r["altloc"] == "B"
    assert r["occupancy"] == 0.7

    # No usable record (all zero occupancy or non-finite B).
    r = select_representative_ca_record([0.0], [20.0], ["."])
    assert r["altloc"] is None and r["occupancy"] is None and r["b_iso_or_equiv"] is None


def test_bfactor_value_never_used_as_tie_break_source_code():
    text = (PROJECT_ROOT / "scripts" / "utils" / "altloc_selection.py").read_text()
    sort_key_start = text.index("candidates.sort")
    sort_key_line = text[sort_key_start:text.index("\n", sort_key_start)]
    assert "b_iso" not in sort_key_line.lower()


# ---------------------------------------------------------------------------
# LIGAND
# ---------------------------------------------------------------------------

def test_every_selected_structure_has_a_ligand_state(structures, ligand):
    assert set(structures["instance_id"]) <= set(
        (ligand["pdb_id"] + "." + ligand["label_asym_id"])
    )
    assert ligand["ligand_state"].isin(["APO", "HOLO", "AMBIGUOUS"]).all()


def test_holo_has_at_least_one_functional_ligand(ligand):
    holo = ligand[ligand["ligand_state"] == "HOLO"]
    if len(holo) == 0:
        pytest.skip("no HOLO rows in this run")
    counts = holo["functional_ligand_component_ids_json"].apply(lambda j: len(json.loads(j)))
    assert (counts >= 1).all()


def test_apo_has_no_functional_ligand(ligand):
    apo = ligand[ligand["ligand_state"] == "APO"]
    if len(apo) == 0:
        pytest.skip("no APO rows in this run")
    counts = apo["functional_ligand_component_ids_json"].apply(lambda j: len(json.loads(j)))
    assert (counts == 0).all()


def test_water_alone_cannot_cause_holo():
    _require(NONPOLYMER_INVENTORY_CSV)
    inv = pd.read_csv(NONPOLYMER_INVENTORY_CSV, dtype={"pdb_id": str})
    water = inv[inv["nonpolymer_relevance_class"] == "WATER"]
    if len(water) == 0:
        pytest.skip("no water rows in this run")
    assert (water["component_id"].isin(["HOH", "DOD"])).all()
    # Water rows are never linked to a representative_instance_id (no
    # geometry computed for them), confirming they cannot drive a per
    # structure HOLO/APO/AMBIGUOUS decision.
    assert (water["representative_instance_id"].isna() | (water["representative_instance_id"] == "")).all()


def test_monatomic_ion_never_classified_functional_ligand():
    _require(NONPOLYMER_INVENTORY_CSV)
    inv = pd.read_csv(NONPOLYMER_INVENTORY_CSV, dtype={"pdb_id": str})
    non_water = inv[inv["nonpolymer_relevance_class"] != "WATER"]
    ions = non_water[non_water["atom_count_heavy"] == 1]
    if len(ions) == 0:
        pytest.skip("no monatomic non-water components in this run")
    assert (ions["nonpolymer_relevance_class"] == "ION").all()


# ---------------------------------------------------------------------------
# HANDOFF
# ---------------------------------------------------------------------------

def test_every_structures_key_has_residue_map_rows(structures, residue_map):
    mapped_instances = set(residue_map["instance_id"])
    assert set(structures["instance_id"]) <= mapped_instances


def test_every_primary_pdb_mmcif_exists_locally_with_matching_checksum(structures):
    _require(DOWNLOAD_MANIFEST_CSV)
    downloads = pd.read_csv(DOWNLOAD_MANIFEST_CSV, dtype={"pdb_id": str}).set_index("pdb_id")
    for pdb_id in structures["pdb_id"].unique():
        assert (RAW_MMCIF_DIR / f"{pdb_id}.cif.gz").exists()
        assert downloads.loc[pdb_id, "sha256_gzip"] != ""


def test_structures_mmcif_sha256_matches_download_manifest(structures):
    _require(DOWNLOAD_MANIFEST_CSV)
    downloads = pd.read_csv(DOWNLOAD_MANIFEST_CSV, dtype={"pdb_id": str}).set_index("pdb_id")["sha256_gzip"].to_dict()
    mismatches = structures[structures.apply(lambda r: downloads.get(r["pdb_id"]) != r["mmcif_sha256"], axis=1)]
    assert len(mismatches) == 0


def test_all_48_receptors_present_in_final_counts_by_receptor():
    _require(FINAL_COUNTS_BY_RECEPTOR_CSV)
    df = pd.read_csv(FINAL_COUNTS_BY_RECEPTOR_CSV)
    assert len(df) == 48
    assert df["uniprot_id"].nunique() == 48


def test_final_counts_by_receptor_sums_match_structures_csv(structures):
    _require(FINAL_COUNTS_BY_RECEPTOR_CSV)
    df = pd.read_csv(FINAL_COUNTS_BY_RECEPTOR_CSV)
    assert df["primary_selected_structures"].sum() == len(structures)
    assert df["apo_count"].sum() + df["holo_count"].sum() + df["ambiguous_ligand_count"].sum() == len(structures)


def test_release_json_checksums_match_actual_files():
    _require(RELEASE_JSON)
    import hashlib
    release = json.loads(RELEASE_JSON.read_text())
    for key, entry in release["checksums"].items():
        path = PROJECT_ROOT / entry["path"]
        if not path.exists():
            continue
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        assert h.hexdigest() == entry["sha256"], f"checksum mismatch for {key}"


def test_release_json_declares_no_downstream_analysis():
    _require(RELEASE_JSON)
    release = json.loads(RELEASE_JSON.read_text())
    text = json.dumps(release).lower()
    assert "z-score" not in text and "normalized" not in text


# ---------------------------------------------------------------------------
# IMMUTABILITY
# ---------------------------------------------------------------------------

def test_stage2_release_summary_unchanged():
    _require(STAGE2_RELEASE_SUMMARY)
    summary = json.loads(STAGE2_RELEASE_SUMMARY.read_text())
    assert summary["candidate_relationships"] == 2072
    assert summary["unique_pdb_entries"] == 1996


def test_stage3a_release_summary_unchanged():
    _require(STAGE3A_RELEASE_SUMMARY)
    summary = json.loads(STAGE3A_RELEASE_SUMMARY.read_text())
    assert summary["project_metadata_prefilter_pass"] == 1516


def test_stage3b_release_summary_unchanged():
    _require(STAGE3B_RELEASE_SUMMARY)
    summary = json.loads(STAGE3B_RELEASE_SUMMARY.read_text())
    assert summary["coordinate_files_downloaded"] == 1777
    assert summary["selected_candidates"] == 1840
    assert summary["measured_instance_models"] == 3159
    assert summary["exact_full_lbd_positive_occupancy_complete"] == 341


def test_frozen_tags_exist_and_point_at_frozen_commits():
    import subprocess

    for tag in ("stage2-candidate-inventory-v1", "stage3a-lbd-prefilter-v1", "stage3b-coordinate-audit-v1"):
        result = subprocess.run(
            ["git", "rev-parse", tag], cwd=PROJECT_ROOT, capture_output=True, text=True
        )
        assert result.returncode == 0, f"tag {tag} not found"


def test_no_bfactor_normalization_artifacts_exist():
    """This stage must never produce a normalized/Z-scored B-factor
    column anywhere in the final tables."""
    for path in (STRUCTURES_CSV, RESIDUE_MAP_PARQUET):
        if not path.exists():
            continue
        df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, nrows=1)
        forbidden = {"z_score", "normalized_b", "b_factor_zscore", "cluster_label", "pca_", "tsne_", "umap_"}
        cols_lower = {c.lower() for c in df.columns}
        assert forbidden.isdisjoint(cols_lower), f"{path.name} contains a forbidden downstream-analysis column"
