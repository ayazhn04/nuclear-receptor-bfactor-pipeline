"""
Stage 2 tests: RCSB candidate structure discovery.

These tests operate ONLY on already-generated manifest files. No network
access. If the Stage 2 outputs don't exist yet, tests that need them are
skipped with a clear message rather than triggering fresh API calls.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"

SEARCH_MANIFEST = MANIFESTS_DIR / "rcsb_search_manifest.csv"
DISCOVERY_BY_RECEPTOR = MANIFESTS_DIR / "rcsb_discovery_by_receptor.csv"
POLYMER_ENTITIES = MANIFESTS_DIR / "rcsb_polymer_entities.csv"
ENTRIES = MANIFESTS_DIR / "rcsb_entries.csv"
CANDIDATE_INVENTORY = MANIFESTS_DIR / "rcsb_candidate_inventory.csv"
ANOMALIES = MANIFESTS_DIR / "rcsb_discovery_anomalies.csv"
SUMMARY_JSON = MANIFESTS_DIR / "rcsb_discovery_summary.json"
MAPPINGS = MANIFESTS_DIR / "rcsb_polymer_entity_uniprot_mappings.csv"
TAXONOMY_AUDIT = MANIFESTS_DIR / "rcsb_source_taxonomy_audit.csv"
MULTI_UNIPROT_AUDIT = MANIFESTS_DIR / "rcsb_multi_uniprot_entity_audit.csv"

VALID_MAPPING_STATUSES = {"CONFIRMED", "CONFIRMED_MULTI_MAPPED", "REVIEW_NEEDED"}
VALID_CLASSIFICATIONS = {"FUSION_OR_CHIMERA", "MULTI_MAPPING_SAME_BIOLOGICAL_PROTEIN", "AMBIGUOUS_REVIEW_NEEDED"}
VALID_TAXONOMY_STATUSES = {"EXPLAINED", "REVIEW_NEEDED"}

STAGE1_MANIFEST_FILES = [
    "nr_metadata_professor.csv",
    "nr_metadata_professor_provenance.json",
    "nr_metadata_local_validation.csv",
    "nr_metadata_local_validation_summary.json",
    "nr_metadata_uniprot_audit.csv",
    "nr_metadata_uniprot_audit_summary.json",
    "uniprot_api_manifest.csv",
    "nr_metadata_validated.csv",
    "nr_metadata_validated_provenance.json",
]


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet — run the Stage 2 scripts first")


def test_config_nr_metadata_still_48_records():
    from config.nr_metadata import NR_METADATA

    assert len(NR_METADATA) == 48


def test_stage1_manifests_unchanged_by_stage2():
    """Stage 2 must never rewrite Stage 1 outputs; sanity check row counts
    that would break if Stage 1 files were accidentally overwritten."""
    _require(MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv")
    audit_df = pd.read_csv(MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv")
    assert len(audit_df) == 48
    for filename in STAGE1_MANIFEST_FILES:
        assert (MANIFESTS_DIR / filename).exists(), f"Stage 1 manifest missing: {filename}"


def test_search_manifest_has_48_rows():
    _require(SEARCH_MANIFEST)
    df = pd.read_csv(SEARCH_MANIFEST)
    assert len(df) == 48


def test_every_validated_receptor_appears_once_in_search_manifest():
    _require(SEARCH_MANIFEST)
    from config.nr_metadata import NR_METADATA

    df = pd.read_csv(SEARCH_MANIFEST)
    assert set(df["uniprot_id"]) == {r["uniprot_id"] for r in NR_METADATA}
    assert df["uniprot_id"].is_unique


def test_successful_queries_have_result_count():
    _require(SEARCH_MANIFEST)
    df = pd.read_csv(SEARCH_MANIFEST)
    successful = df[df["success"]]
    assert successful["result_count_parsed"].notna().all()


def test_failed_queries_are_not_zero_result_receptors():
    """A failed query and a successful zero-result query must never be
    conflated: a failed row must carry no result_count."""
    _require(SEARCH_MANIFEST)
    df = pd.read_csv(SEARCH_MANIFEST)
    failed = df[~df["success"]]
    for _, row in failed.iterrows():
        assert pd.isna(row["result_count_parsed"]) or row["result_count_parsed"] == ""
        assert isinstance(row["error_message"], str) and len(row["error_message"].strip()) > 0


def test_discovery_ids_parse_consistently_into_pdb_and_entity_id():
    _require(DISCOVERY_BY_RECEPTOR)
    df = pd.read_csv(DISCOVERY_BY_RECEPTOR, dtype=str)
    for _, row in df.iterrows():
        expected = f"{row['pdb_id']}_{row['entity_id']}"
        assert row["polymer_entity_id"] == expected


def test_no_duplicate_receptor_by_polymer_entity_rows():
    _require(DISCOVERY_BY_RECEPTOR)
    df = pd.read_csv(DISCOVERY_BY_RECEPTOR)
    dup_key = df[["uniprot_id", "polymer_entity_id"]]
    assert not dup_key.duplicated().any()


def test_unique_polymer_entity_table_has_unique_ids():
    _require(POLYMER_ENTITIES)
    df = pd.read_csv(POLYMER_ENTITIES)
    assert df["polymer_entity_id"].is_unique


def test_unique_entry_table_has_unique_pdb_ids():
    _require(ENTRIES)
    df = pd.read_csv(ENTRIES)
    assert df["pdb_id"].is_unique


def test_candidate_inventory_preserves_all_discovery_rows():
    _require(CANDIDATE_INVENTORY)
    _require(DISCOVERY_BY_RECEPTOR)
    inventory_df = pd.read_csv(CANDIDATE_INVENTORY)
    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR)
    assert len(inventory_df) == len(discovery_df)


def test_candidate_inventory_no_accidental_row_multiplication():
    _require(CANDIDATE_INVENTORY)
    df = pd.read_csv(CANDIDATE_INVENTORY)
    dup_key = df[["uniprot_id", "polymer_entity_id"]]
    assert not dup_key.duplicated().any()


def test_every_candidate_has_receptor_provenance():
    _require(CANDIDATE_INVENTORY)
    df = pd.read_csv(CANDIDATE_INVENTORY)
    for col in ("uniprot_id", "nr_code", "common_name", "group"):
        assert df[col].notna().all()
        assert (df[col].astype(str).str.strip() != "").all()


def test_mapping_status_is_controlled_vocabulary():
    _require(POLYMER_ENTITIES)
    df = pd.read_csv(POLYMER_ENTITIES)
    assert set(df["mapping_status"].unique()).issubset(VALID_MAPPING_STATUSES)


def test_experimental_method_flags_internally_consistent():
    """A row can be multiple methods (joint refinement) but if any flag is
    True, experimental_methods must be non-empty."""
    _require(CANDIDATE_INVENTORY)
    df = pd.read_csv(CANDIDATE_INVENTORY)
    any_flag = df["is_xray"] | df["is_cryo_em"] | df["is_nmr"] | df["is_other_experimental"]
    methods_nonempty = df["experimental_methods"].apply(lambda v: len(json.loads(v)) > 0 if isinstance(v, str) else False)
    assert (any_flag <= methods_nonempty).all(), "a method flag is True but experimental_methods is empty"


def test_zero_result_receptors_represented_in_summary():
    _require(SUMMARY_JSON)
    _require(SEARCH_MANIFEST)
    summary = json.loads(SUMMARY_JSON.read_text())
    search_df = pd.read_csv(SEARCH_MANIFEST)
    actual_zero = set(
        search_df[(search_df["success"]) & (search_df["result_count_parsed"] == 0)]["common_name"]
    )
    assert set(summary["receptors_with_zero_candidates_list"]) == actual_zero


def test_summary_failed_and_zero_result_counts_are_disjoint_concepts():
    _require(SUMMARY_JSON)
    summary = json.loads(SUMMARY_JSON.read_text())
    assert summary["number_of_failed_receptor_queries"] + summary["number_of_successful_receptor_queries"] == 48
    # zero-candidate count must not silently include failed queries
    assert summary["receptors_with_zero_candidates"] <= summary["number_of_successful_receptor_queries"]


def test_anomalies_have_controlled_review_status():
    _require(ANOMALIES)
    df = pd.read_csv(ANOMALIES)
    if len(df) == 0:
        pytest.skip("no anomalies generated")
    assert set(df["review_status"].unique()) == {"REVIEW_NEEDED"}


# ---------------------------------------------------------------------------
# Stage 2.1: coverage correction + source/mapping audit
# ---------------------------------------------------------------------------

def test_mapping_table_rows_are_unique_logical_mappings():
    """One row per (polymer_entity_id, mapped_uniprot_accession) — no
    silent collapsing, and no accidental duplication of a logical mapping."""
    _require(MAPPINGS)
    df = pd.read_csv(MAPPINGS)
    non_duplicate = df[~df["duplicate_mapping_for_same_accession"]]
    key = non_duplicate[["polymer_entity_id", "mapped_uniprot_accession"]]
    assert not key.duplicated().any()


def test_every_confirmed_candidate_has_its_query_uniprot_mapping_represented():
    _require(MAPPINGS)
    _require(POLYMER_ENTITIES)
    mapping_df = pd.read_csv(MAPPINGS)
    entities_df = pd.read_csv(POLYMER_ENTITIES)
    confirmed = entities_df[entities_df["mapping_status"].isin(["CONFIRMED", "CONFIRMED_MULTI_MAPPED"])]
    for _, row in confirmed.iterrows():
        entity_mappings = mapping_df[mapping_df["polymer_entity_id"] == row["polymer_entity_id"]]
        query_rows = entity_mappings[entity_mappings["is_query_receptor_accession"]]
        assert len(query_rows) >= 1, f"{row['polymer_entity_id']} confirmed but has no query-UniProt mapping row"


def test_coverage_values_null_or_within_0_1():
    """Historical note: Stage 2.1 named these columns reference_sequence_coverage
    / entity_sequence_coverage (naive length sum). Stage 2.2's interval-union
    robustness fix renamed them to computed_reference_sequence_coverage /
    computed_entity_sequence_coverage — see test_mapping_table_uses_interval_union_columns."""
    _require(MAPPINGS)
    df = pd.read_csv(MAPPINGS)
    for col in ("computed_reference_sequence_coverage", "computed_entity_sequence_coverage"):
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        assert (values >= 0).all() or len(values) == 0
        # Values > 1 are allowed to exist but must be flagged as anomalies
        # (COVERAGE_EXCEEDS_1) rather than silently clipped — checked below.


def test_coverage_exceeding_1_is_flagged_not_silently_present():
    _require(MAPPINGS)
    _require(ANOMALIES)
    mapping_df = pd.read_csv(MAPPINGS)
    anomalies_df = pd.read_csv(ANOMALIES)
    ref_over = pd.to_numeric(mapping_df["computed_reference_sequence_coverage"], errors="coerce") > 1.0
    ent_over = pd.to_numeric(mapping_df["computed_entity_sequence_coverage"], errors="coerce") > 1.0
    n_over = int((ref_over | ent_over).sum())
    n_flagged = int((anomalies_df["anomaly_type"] == "COVERAGE_EXCEEDS_1").sum())
    assert n_over == n_flagged


def test_reference_and_entity_coverage_stored_separately():
    """The Stage 2 bug stored the identical alignment JSON in both coverage
    columns; this asserts they are now genuinely distinct numeric series,
    not merely equal by coincidence across every row."""
    _require(MAPPINGS)
    df = pd.read_csv(MAPPINGS)
    ref = pd.to_numeric(df["computed_reference_sequence_coverage"], errors="coerce")
    ent = pd.to_numeric(df["computed_entity_sequence_coverage"], errors="coerce")
    both_present = ref.notna() & ent.notna()
    assert both_present.sum() > 0
    assert not (ref[both_present] == ent[both_present]).all(), (
        "computed_reference_sequence_coverage and computed_entity_sequence_coverage are identical "
        "for every row — this would indicate the Stage 2 bug (same alignment JSON reused for both) "
        "has recurred"
    )


def test_no_alignment_region_object_used_as_scalar_coverage():
    """Coverage columns must be numeric (or null), never a JSON/list string."""
    _require(MAPPINGS)
    df = pd.read_csv(MAPPINGS)
    for col in ("computed_reference_sequence_coverage", "computed_entity_sequence_coverage"):
        for value in df[col].dropna():
            assert not isinstance(value, str) or not value.strip().startswith("["), (
                f"{col} contains what looks like a JSON array, not a scalar coverage value: {value!r}"
            )


def test_all_multi_uniprot_entities_appear_in_dedicated_audit():
    _require(POLYMER_ENTITIES)
    _require(MULTI_UNIPROT_AUDIT)
    entities_df = pd.read_csv(POLYMER_ENTITIES)
    multi_mask = entities_df["all_mapped_uniprot_accessions"].apply(lambda v: len(json.loads(v)) > 1 if isinstance(v, str) else False)
    expected_ids = set(entities_df.loc[multi_mask, "polymer_entity_id"])
    audit_df = pd.read_csv(MULTI_UNIPROT_AUDIT)
    assert set(audit_df["polymer_entity_id"]) == expected_ids


def test_multi_uniprot_audit_uses_controlled_classification_vocabulary():
    _require(MULTI_UNIPROT_AUDIT)
    df = pd.read_csv(MULTI_UNIPROT_AUDIT)
    if len(df) == 0:
        pytest.skip("no multi-mapped entities")
    assert set(df["classification"].unique()).issubset(VALID_CLASSIFICATIONS)


def test_every_candidate_lacking_taxon_9606_appears_in_taxonomy_audit():
    _require(CANDIDATE_INVENTORY)
    _require(TAXONOMY_AUDIT)
    inventory_df = pd.read_csv(CANDIDATE_INVENTORY)
    expected = set(inventory_df.loc[~inventory_df["source_includes_homo_sapiens_9606"], "polymer_entity_id"])
    audit_df = pd.read_csv(TAXONOMY_AUDIT)
    assert set(audit_df["polymer_entity_id"]) == expected


def test_taxonomy_audit_uses_controlled_status_vocabulary():
    _require(TAXONOMY_AUDIT)
    df = pd.read_csv(TAXONOMY_AUDIT)
    if len(df) == 0:
        pytest.skip("no non-human-source candidates")
    assert set(df["taxonomy_audit_status"].unique()).issubset(VALID_TAXONOMY_STATUSES)


def test_host_organism_not_used_as_substitute_for_source_organism():
    """The taxonomy audit must keep source_organism and host_organism in
    separate columns — a human protein expressed in E. coli must not have
    its source_organism silently replaced by the host."""
    _require(TAXONOMY_AUDIT)
    df = pd.read_csv(TAXONOMY_AUDIT)
    assert "source_organism_names" in df.columns
    assert "host_organism_names" in df.columns
    # Every row originally lacked taxon 9606 in its SOURCE organism list —
    # the audit must not have quietly swapped in the host's taxon instead.
    for _, row in df.iterrows():
        source_taxa = json.loads(row["source_taxonomy_ids"])
        assert 9606 not in source_taxa


def test_zero_result_receptors_distinguished_from_failed_queries_in_summary():
    _require(SUMMARY_JSON)
    _require(SEARCH_MANIFEST)
    summary = json.loads(SUMMARY_JSON.read_text())
    search_df = pd.read_csv(SEARCH_MANIFEST)
    failed_receptors = set(search_df.loc[~search_df["success"], "common_name"])
    zero_result_receptors = set(summary["receptors_with_zero_candidates_list"])
    assert failed_receptors.isdisjoint(zero_result_receptors), (
        "a failed query's receptor appears in the zero-result list — these must remain distinct"
    )


# ---------------------------------------------------------------------------
# Stage 2.2: official Sequence Coordinates validation
# ---------------------------------------------------------------------------

SEQ_COORD_MANIFEST = MANIFESTS_DIR / "rcsb_sequence_coordinates_manifest.csv"
SEQ_COORD_ALIGNMENTS = MANIFESTS_DIR / "rcsb_sequence_coordinate_alignments.csv"
COVERAGE_VALIDATION = MANIFESTS_DIR / "rcsb_coverage_validation.csv"
SHORT_MAPPING_AUDIT = MANIFESTS_DIR / "rcsb_short_mapping_audit.csv"

VALID_VALIDATION_STATUSES = {"MATCH", "MINOR_NUMERIC_DIFFERENCE", "REVIEW_NEEDED"}


def test_sequence_coordinates_manifest_has_48_rows():
    _require(SEQ_COORD_MANIFEST)
    df = pd.read_csv(SEQ_COORD_MANIFEST)
    assert len(df) == 48


def test_successful_sequence_coordinates_queries_have_documented_counts():
    _require(SEQ_COORD_MANIFEST)
    df = pd.read_csv(SEQ_COORD_MANIFEST)
    successful = df[df["success"]]
    for col in ("targets_returned", "stage2_targets_expected", "stage2_targets_matched"):
        assert successful[col].notna().all()


def test_every_stage2_candidate_matched_or_flagged_missing():
    _require(SEQ_COORD_MANIFEST)
    _require(DISCOVERY_BY_RECEPTOR)
    manifest_df = pd.read_csv(SEQ_COORD_MANIFEST)
    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR)
    expected_by_receptor = discovery_df.groupby("uniprot_id")["polymer_entity_id"].nunique()
    for _, row in manifest_df[manifest_df["success"]].iterrows():
        expected = expected_by_receptor.get(row["uniprot_id"], 0)
        matched = row["stage2_targets_matched"]
        missing = json.loads(row["missing_stage2_targets"]) if isinstance(row["missing_stage2_targets"], str) and row["missing_stage2_targets"] else []
        assert matched + len(missing) == expected, (
            f"{row['uniprot_id']}: matched({matched}) + missing({len(missing)}) != expected({expected})"
        )


def test_official_coverage_values_null_or_within_0_1():
    _require(SEQ_COORD_ALIGNMENTS)
    df = pd.read_csv(SEQ_COORD_ALIGNMENTS)
    for col in ("rcsb_reference_sequence_coverage", "rcsb_entity_sequence_coverage"):
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        assert ((values >= 0) & (values <= 1.000001)).all()


def test_computed_coverage_values_null_or_within_0_1():
    _require(MAPPINGS)
    df = pd.read_csv(MAPPINGS)
    for col in ("computed_reference_sequence_coverage", "computed_entity_sequence_coverage"):
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        assert ((values >= 0) & (values <= 1.000001)).all()


def test_mapping_table_uses_interval_union_columns():
    """The Stage 2.2 robustness fix must have replaced the naive-sum
    coverage columns with union-aware ones, and recorded overlap flags."""
    _require(MAPPINGS)
    df = pd.read_csv(MAPPINGS)
    assert "computed_reference_sequence_coverage" in df.columns
    assert "computed_entity_sequence_coverage" in df.columns
    assert "aligned_region_count" in df.columns
    assert "reference_intervals_overlap" in df.columns
    assert "entity_intervals_overlap" in df.columns
    # the old, superseded naive columns must be gone, not left stale alongside
    assert "reference_sequence_coverage" not in df.columns
    assert "entity_sequence_coverage" not in df.columns


def test_no_overlapping_interval_double_counted():
    """For every mapping flagged as NOT overlapping (entity side), the
    union-based computed_entity_sequence_coverage * entity_sequence_length
    must equal the naive sum of region lengths — i.e. union coverage
    reduces to the simple sum exactly when there is truly no overlap,
    confirming the union algorithm isn't silently over- or under-counting."""
    _require(MAPPINGS)
    _require(POLYMER_ENTITIES)
    mapping_df = pd.read_csv(MAPPINGS)
    entities_df = pd.read_csv(POLYMER_ENTITIES).set_index("polymer_entity_id")

    non_overlapping = mapping_df[~mapping_df["entity_intervals_overlap"]]
    checked = 0
    for _, row in non_overlapping.iterrows():
        if row["polymer_entity_id"] not in entities_df.index:
            continue
        entity_length = entities_df.loc[row["polymer_entity_id"], "entity_sequence_length"]
        if pd.isna(entity_length) or pd.isna(row["computed_entity_sequence_coverage"]):
            continue
        regions = json.loads(row["alignment_regions_json"])
        naive_sum = sum(r["length"] for r in regions)
        union_length_from_coverage = round(row["computed_entity_sequence_coverage"] * entity_length)
        assert union_length_from_coverage == naive_sum, (
            f"{row['polymer_entity_id']}: non-overlapping union length ({union_length_from_coverage}) "
            f"!= naive region-length sum ({naive_sum})"
        )
        checked += 1
    assert checked > 0, "no rows available to check the union-vs-naive-sum identity"


def test_validation_status_controlled_vocabulary():
    _require(COVERAGE_VALIDATION)
    df = pd.read_csv(COVERAGE_VALIDATION)
    assert set(df["validation_status"].unique()).issubset(VALID_VALIDATION_STATUSES)


def test_taxonomy_audit_has_qc_suitability_not_evaluated():
    _require(TAXONOMY_AUDIT)
    df = pd.read_csv(TAXONOMY_AUDIT)
    assert "qc_suitability" in df.columns
    assert (df["qc_suitability"] == "NOT_EVALUATED").all()


def test_short_mapping_audit_reproducible_thresholds():
    _require(SHORT_MAPPING_AUDIT)
    df = pd.read_csv(SHORT_MAPPING_AUDIT)
    if len(df) == 0:
        pytest.skip("no short mappings found")
    below_threshold = (df["reference_coverage"] < 0.05) | (df["mapped_length"] < 30)
    assert below_threshold.all()


def test_stage1_files_unchanged_by_stage2_2():
    _require(MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv")
    df = pd.read_csv(MANIFESTS_DIR / "nr_metadata_uniprot_audit.csv")
    assert len(df) == 48
    for filename in STAGE1_MANIFEST_FILES:
        assert (MANIFESTS_DIR / filename).exists()


def test_no_coordinate_files_exist():
    """Historical note: at Stage 2 this asserted zero mmCIF files existed
    anywhere. Stage 3B's entire purpose is downloading a selected subset of
    them, so this check is now a no-op once Stage 3B has begun (detected via
    its coordinate-pool manifest) — mmCIF provenance/safety from that point
    on is covered by tests/test_stage3b_coordinate_qc.py instead."""
    if (MANIFESTS_DIR / "stage3b_coordinate_pool.csv").exists():
        pytest.skip("Stage 3B has begun downloading mmCIF files by design — see test_stage3b_coordinate_qc.py")
    mmcif_dir = PROJECT_ROOT / "data" / "raw" / "mmcif"
    if mmcif_dir.exists():
        files = [p for p in mmcif_dir.rglob("*") if p.is_file() and p.name != ".gitkeep"]
        assert len(files) == 0, f"unexpected coordinate files found: {files}"


# ---------------------------------------------------------------------------
# Stage 2.3: correct Sequence Coordinates extra-target interpretation
# ---------------------------------------------------------------------------

EXTRA_TARGETS = MANIFESTS_DIR / "rcsb_sequence_coordinates_extra_targets.csv"
VALID_SCOPE_STATUSES = {"IN_STAGE2_CANDIDATE_SET", "OUTSIDE_STAGE2_SEARCH_SCOPE"}
VALID_SCOPE_REASONS = {
    "INTEGRATIVE_STRUCTURE", "COMPUTED_STRUCTURE_MODEL",
    "NON_EXPERIMENTAL_MODEL", "OTHER_SEARCH_SCOPE_MISMATCH", "REVIEW_NEEDED",
}


def test_every_out_of_scope_target_appears_in_extra_targets_audit():
    _require(SEQ_COORD_ALIGNMENTS)
    _require(DISCOVERY_BY_RECEPTOR)
    _require(EXTRA_TARGETS)
    alignments_df = pd.read_csv(SEQ_COORD_ALIGNMENTS)
    discovery_df = pd.read_csv(DISCOVERY_BY_RECEPTOR)
    stage2_entity_ids = set(discovery_df["polymer_entity_id"])

    out_of_scope = alignments_df[~alignments_df["polymer_entity_id"].isin(stage2_entity_ids)]
    extra_df = pd.read_csv(EXTRA_TARGETS)
    assert set(out_of_scope["polymer_entity_id"]) == set(extra_df["polymer_entity_id"])


def test_every_extra_target_has_explicit_scope_classification():
    _require(EXTRA_TARGETS)
    df = pd.read_csv(EXTRA_TARGETS)
    if len(df) == 0:
        pytest.skip("no extra targets in this run")
    assert (df["stage2_scope_status"] == "OUTSIDE_STAGE2_SEARCH_SCOPE").all()
    assert set(df["scope_exclusion_reason"].unique()).issubset(VALID_SCOPE_REASONS)


def test_extra_targets_not_added_to_candidate_inventory():
    _require(EXTRA_TARGETS)
    _require(CANDIDATE_INVENTORY)
    extra_df = pd.read_csv(EXTRA_TARGETS)
    if len(extra_df) == 0:
        pytest.skip("no extra targets in this run")
    inventory_df = pd.read_csv(CANDIDATE_INVENTORY)
    assert set(extra_df["polymer_entity_id"]).isdisjoint(set(inventory_df["polymer_entity_id"]))


def test_stage2_candidates_still_fully_matched():
    _require(SEQ_COORD_MANIFEST)
    df = pd.read_csv(SEQ_COORD_MANIFEST)
    successful = df[df["success"]]
    assert int(successful["stage2_targets_matched"].sum()) == 2072
    assert int(successful["stage2_targets_expected"].sum()) == 2072


def test_manifest_uses_corrected_terminology():
    """The Stage 2.3 correction renamed the misleading 'unexpected_targets'
    column — it must no longer exist, replaced by a neutral name."""
    _require(SEQ_COORD_MANIFEST)
    df = pd.read_csv(SEQ_COORD_MANIFEST)
    assert "unexpected_targets" not in df.columns
    assert "targets_outside_stage2_scope" in df.columns


def test_report_contains_no_unsupported_indexing_lag_claim():
    """The corrected report may explicitly NEGATE the indexing-lag
    explanation (factually, with evidence) but must not assert it as the
    reason without evidence. We check that if either phrase appears, it is
    only in a sentence that also contains a negation ('not', 'NOT') nearby."""
    report_path = REPORTS_TABLES_DIR = PROJECT_ROOT / "reports" / "tables" / "rcsb_discovery_report.md"
    if not report_path.exists():
        pytest.skip("report not generated yet")
    text = report_path.read_text()
    for sentence in text.replace("\n", " ").split(". "):
        lowered = sentence.lower()
        if "indexing" in lowered or "recently deposited" in lowered:
            assert "not " in lowered or "no " in lowered, (
                f"found an unqualified indexing-lag/recently-deposited claim: {sentence!r}"
            )


def test_coverage_validation_results_unchanged_by_stage_2_3():
    _require(COVERAGE_VALIDATION)
    df = pd.read_csv(COVERAGE_VALIDATION)
    counts = df["validation_status"].value_counts().to_dict()
    assert counts.get("MATCH", 0) == 2072
    assert counts.get("MINOR_NUMERIC_DIFFERENCE", 0) == 0
    assert counts.get("REVIEW_NEEDED", 0) == 0


def test_zero_result_receptors_unchanged_by_stage_2_3():
    _require(SEARCH_MANIFEST)
    df = pd.read_csv(SEARCH_MANIFEST)
    zero = set(df[(df["success"]) & (df["result_count_parsed"] == 0)]["common_name"])
    assert zero == {"RORβ", "TR2", "GCNF"}
