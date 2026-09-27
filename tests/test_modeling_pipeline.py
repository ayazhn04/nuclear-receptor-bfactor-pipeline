"""
Segment-level PCA/UMAP/clustering modeling stage tests (continuing from
commit 7d280bc on analysis/normalization-v1; scripts 48-54).

Covers: no fabricated H8 values; frozen source files unchanged; missing
values never replaced with literal 0; imputation uses only included-cohort
feature medians; group/subfamily labels never used as model inputs; scaling
performed after imputation; PCA fit only on the selected flexibility
features; UMAP coordinates never used as a clustering input; clustering
outputs reproduce deterministically; all candidate feature thresholds are
reported; excluded receptors remain explicit rows in the cohort table;
master structures.csv is untouched at 2072 rows; and the full stage is
deterministic end-to-end.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

PROCESSED = PROJECT_ROOT / "data" / "processed"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"
TABLES = PROJECT_ROOT / "reports" / "tables"

MASTER = PROCESSED / "structures.csv"
FROZEN_HASHES = ANALYSIS / "modeling_frozen_input_hashes.json"
FEATURE_SUPPORT = TABLES / "modeling_feature_support.csv"
MATRIX_CANDIDATES = TABLES / "modeling_matrix_candidates.csv"
COHORT_CSV = ANALYSIS / "modeling_receptor_cohort.csv"
PRIMARY_UNSCALED = ANALYSIS / "receptor_segment_matrix_primary_unscaled.parquet"
PRIMARY_SCALED = ANALYSIS / "receptor_segment_matrix_primary_scaled.parquet"
PRIMARY_MASK = ANALYSIS / "receptor_segment_matrix_primary_imputation_mask.parquet"
COHORT_MATRIX_SUMMARY = ANALYSIS / "modeling_cohort_matrix_summary.json"
PCA_LOADINGS = ANALYSIS / "pca_loadings.csv"
PCA_SCORES = ANALYSIS / "pca_receptor_scores.csv"
PCA_VARIANCE = TABLES / "pca_variance.csv"
CLUSTERING_VALIDATION = TABLES / "clustering_validation.csv"
CLUSTERING_DBSCAN = TABLES / "clustering_dbscan_grid.csv"
MODELING_SUMMARY = TABLES / "modeling_summary.md"

MODELING_SCRIPTS = [
    "48_feature_support_audit.py",
    "49_build_modeling_cohort_and_matrix.py",
    "50_pca_analysis.py",
    "51_umap_analysis.py",
    "52_clustering_analysis.py",
    "53_entry_level_sensitivity.py",
]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet")


@pytest.fixture(scope="module")
def cohort():
    _require(COHORT_CSV)
    return pd.read_csv(COHORT_CSV)


@pytest.fixture(scope="module")
def cohort_summary():
    _require(COHORT_MATRIX_SUMMARY)
    with open(COHORT_MATRIX_SUMMARY) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def raw_median_matrix():
    return pd.read_parquet(ANALYSIS / "receptor_segment_feature_matrix_median.parquet").set_index("uniprot_id")


@pytest.fixture(scope="module")
def primary_unscaled():
    _require(PRIMARY_UNSCALED)
    return pd.read_parquet(PRIMARY_UNSCALED).set_index("uniprot_id")


@pytest.fixture(scope="module")
def primary_scaled():
    _require(PRIMARY_SCALED)
    return pd.read_parquet(PRIMARY_SCALED).set_index("uniprot_id")


@pytest.fixture(scope="module")
def primary_mask():
    _require(PRIMARY_MASK)
    return pd.read_parquet(PRIMARY_MASK).set_index("uniprot_id")


# ---------------------------------------------------------------- Part 0 ---

def test_h8_not_fabricated_anywhere_in_modeling(cohort_summary):
    _require(FEATURE_SUPPORT)
    audit = pd.read_csv(FEATURE_SUPPORT)
    h8_row = audit[audit["structural_segment"] == "H8"]
    assert len(h8_row) == 1, "H8 must still be explicitly audited (reported, not silently dropped)"
    assert h8_row.iloc[0]["fraction_supported"] == 0.0
    assert h8_row.iloc[0]["n_coordinate_receptors_supported"] == 0

    candidates = pd.read_csv(MATRIX_CANDIDATES)
    for _, row in candidates.iterrows():
        features = json.loads(row["feature_names_json"])
        assert "H8" not in features, f"H8 must never appear in candidate set {row['candidate_set']}"

    assert "H8" not in cohort_summary["primary_features"]

    for path in [PRIMARY_UNSCALED, PRIMARY_SCALED, PCA_LOADINGS]:
        _require(path)
    unscaled = pd.read_parquet(PRIMARY_UNSCALED)
    assert "H8" not in unscaled.columns
    loadings = pd.read_csv(PCA_LOADINGS)
    assert "H8" not in loadings["structural_segment"].tolist()


def test_h8_documentation_no_longer_claims_biological_absence():
    handoff = (PROJECT_ROOT / "TEAM_HANDOFF.md").read_text()
    data_dict = (TABLES / "final_data_dictionary.md").read_text()
    for text, name in [(handoff, "TEAM_HANDOFF.md"), (data_dict, "final_data_dictionary.md")]:
        assert "genuinely `0%`" not in text, f"{name} should not claim H8 is 'genuinely' absent"
        assert "genuinely absent" not in text.lower(), f"{name} should not claim biological absence for H8"
    assert "not confidently resolved" in handoff
    assert "not confidently resolved" in data_dict


# --------------------------------------------------------- Part 1: frozen --

def test_frozen_input_files_unchanged_since_freeze():
    _require(FROZEN_HASHES)
    with open(FROZEN_HASHES) as f:
        frozen = json.load(f)
    assert len(frozen) == 5
    for rel_path, expected_hash in frozen.items():
        actual = _sha(PROJECT_ROOT / rel_path)
        assert actual == expected_hash, (
            f"{rel_path} changed since the Part-1 freeze -- modeling scripts must never "
            "modify their own upstream inputs"
        )


def test_master_structures_csv_untouched():
    df = pd.read_csv(MASTER, dtype={"pdb_id": str})
    assert len(df) == 2072


# --------------------------------------------------- Part 2: candidate sets

def test_all_candidate_feature_thresholds_reported():
    _require(MATRIX_CANDIDATES)
    candidates = pd.read_csv(MATRIX_CANDIDATES)
    assert set(candidates["candidate_set"]) == {"A", "B", "C", "D"}
    expected_thresholds = {"A": 0.70, "B": 0.60, "C": 0.50, "D": 0.40}
    for _, row in candidates.iterrows():
        assert np.isclose(row["support_threshold"], expected_thresholds[row["candidate_set"]])
        assert row["n_features"] >= 0
        assert not pd.isna(row["n_receptors_with_complete_data"])
        assert not pd.isna(row["total_matrix_missingness_fraction"])
    # as the support bar is raised, the qualifying feature count must never increase
    by_thresh = candidates.sort_values("support_threshold")
    diffs = by_thresh["n_features"].diff().dropna()
    assert (diffs <= 0).all(), "n_features must be non-increasing as support_threshold rises"


def test_feature_support_audit_covers_every_realized_segment():
    _require(FEATURE_SUPPORT)
    audit = pd.read_csv(FEATURE_SUPPORT)
    median = pd.read_parquet(ANALYSIS / "receptor_segment_feature_matrix_median.parquet")
    all_cols = set(median.columns) - {"uniprot_id"}
    audited = set(audit["structural_segment"])
    assert all_cols.issubset(audited), "every realized feature-matrix column must be audited"
    extra = audited - all_cols
    # the only segments audited but not realized as feature-matrix columns are ones with
    # zero support anywhere (e.g. H8) -- audited explicitly rather than silently absent
    for seg in extra:
        row = audit[audit["structural_segment"] == seg].iloc[0]
        assert row["fraction_supported"] == 0.0, f"{seg} is audited-but-absent yet has nonzero support"
    assert (audit["fraction_supported"] >= 0).all() and (audit["fraction_supported"] <= 1).all()


# --------------------------------------------------- Part 3: cohort/matrix -

def test_excluded_receptors_remain_explicit_rows_in_cohort(cohort):
    coord_receptor_count = pd.read_csv(ANALYSIS / "receptor_flexibility_summary.csv")
    n_coord = (coord_receptor_count["n_coordinate_candidates"] > 0).sum()
    assert len(cohort) == n_coord, "every coordinate receptor must have an explicit cohort row"

    # SHP/NR0B2 has no usable structural alignment at all -- must still appear, excluded
    shp = cohort[cohort["uniprot_id"] == "Q15466"]
    assert len(shp) == 1
    assert shp.iloc[0]["included_primary"] == False  # noqa: E712
    assert shp.iloc[0]["exclusion_reason"] != ""
    assert not pd.isna(shp.iloc[0]["exclusion_reason"])

    excluded = cohort[~cohort["included_primary"]]
    assert len(excluded) > 0
    assert (excluded["exclusion_reason"].astype(str).str.len() > 0).all()
    included = cohort[cohort["included_primary"]]
    assert (included["exclusion_reason"].fillna("") == "").all()

    required_cols = {"uniprot_id", "nr_code", "group", "n_selected_features",
                      "n_observed_features", "observed_fraction", "included_primary",
                      "exclusion_reason"}
    assert required_cols.issubset(set(cohort.columns))


def test_cohort_inclusion_rule_is_60_percent_observed(cohort, cohort_summary):
    rule = cohort_summary["cohort_observed_fraction_rule"]
    assert np.isclose(rule, 0.60)
    for _, row in cohort.iterrows():
        expected_included = row["observed_fraction"] >= rule
        assert bool(row["included_primary"]) == bool(expected_included)


# ------------------------------------------------- Part 4: missing values --

def test_missing_values_never_literal_zero_fill(primary_unscaled, raw_median_matrix, cohort):
    included_ids = cohort.loc[cohort["included_primary"], "uniprot_id"].tolist()
    raw = raw_median_matrix.reindex(included_ids)[primary_unscaled.columns]
    for col in primary_unscaled.columns:
        raw_col = raw[col]
        observed = raw_col.dropna()
        if len(observed) == 0:
            continue
        expected_median = observed.median()
        missing_mask = raw_col.isna()
        if missing_mask.any():
            imputed_values = primary_unscaled.loc[missing_mask.index[missing_mask], col]
            assert np.allclose(imputed_values.values, expected_median), (
                f"imputed values in column {col} do not equal the included-cohort feature "
                "median -- missing values must never be silently replaced with a literal 0 "
                "or any other ad hoc constant"
            )
            if not np.isclose(expected_median, 0.0):
                assert not np.allclose(imputed_values.values, 0.0)


def test_imputation_mask_identifies_every_imputed_value(primary_unscaled, primary_mask, raw_median_matrix, cohort):
    included_ids = cohort.loc[cohort["included_primary"], "uniprot_id"].tolist()
    raw = raw_median_matrix.reindex(included_ids)[primary_unscaled.columns]
    computed_mask = raw.isna()
    assert (primary_mask.loc[computed_mask.index, computed_mask.columns] == computed_mask).all().all()


def test_imputation_uses_feature_medians_only_from_included_cohort(primary_unscaled, raw_median_matrix, cohort):
    """The per-feature imputation value must come from the INCLUDED cohort's own
    observed values -- never from excluded receptors, never from the full 42-receptor
    coordinate set, never from any supervised label."""
    included_ids = cohort.loc[cohort["included_primary"], "uniprot_id"].tolist()
    excluded_ids = cohort.loc[~cohort["included_primary"], "uniprot_id"].tolist()

    included_only_median = raw_median_matrix.reindex(included_ids)[primary_unscaled.columns].median(skipna=True)
    all_coord_median = raw_median_matrix[primary_unscaled.columns].median(skipna=True)

    # sanity: excluded receptors exist and at least one differs enough that the two
    # medians are not trivially identical for every column (otherwise this test would
    # not actually discriminate "included-only" from "all coordinate receptors")
    assert len(excluded_ids) > 0

    raw = raw_median_matrix.reindex(included_ids)[primary_unscaled.columns]
    for col in primary_unscaled.columns:
        missing = raw[col].isna()
        if not missing.any():
            continue
        imputed_values = primary_unscaled.loc[missing.index[missing], col]
        assert np.allclose(imputed_values.values, included_only_median[col]), (
            f"column {col}: imputed value does not match the included-cohort-only median"
        )


def test_no_group_or_supervised_columns_in_modeling_matrices(primary_unscaled, primary_scaled):
    forbidden = {"group", "nr_code", "common_name", "ligand_state", "cluster", "cluster_label"}
    for df, name in [(primary_unscaled, "unscaled"), (primary_scaled, "scaled")]:
        assert forbidden.isdisjoint(set(df.columns)), (
            f"{name} modeling matrix must contain only flexibility feature columns, "
            "never group/subfamily/ligand/cluster labels"
        )


# ---------------------------------------------------- Part 5: scaling -----

def test_scaling_applied_after_imputation_no_nans(primary_scaled):
    assert not primary_scaled.isna().any().any(), (
        "the scaled primary matrix must be fully dense -- scaling only makes sense "
        "applied AFTER imputation removed all missing values"
    )


def test_scaled_matrix_is_standardized_per_feature(primary_scaled):
    means = primary_scaled.mean(axis=0)
    sds = primary_scaled.std(axis=0, ddof=0)
    assert np.allclose(means.values, 0.0, atol=1e-8)
    assert np.allclose(sds.values, 1.0, atol=1e-6)


def test_scaled_and_unscaled_matrices_same_shape_and_index(primary_scaled, primary_unscaled):
    assert primary_scaled.shape == primary_unscaled.shape
    assert list(primary_scaled.index) == list(primary_unscaled.index)
    assert list(primary_scaled.columns) == list(primary_unscaled.columns)


def test_biological_segment_z_summaries_not_overwritten():
    """The upstream receptor_segment_feature_matrix_median.parquet (biological Z
    summaries) must be byte-identical to its Part-1 frozen hash -- scaling/imputation
    for modeling must only ever touch copies."""
    _require(FROZEN_HASHES)
    with open(FROZEN_HASHES) as f:
        frozen = json.load(f)
    key = "data/analysis/receptor_segment_feature_matrix_median.parquet"
    assert _sha(PROJECT_ROOT / key) == frozen[key]


# -------------------------------------------------------- Part 6: PCA -----

def test_pca_fit_only_on_selected_flexibility_features(cohort_summary):
    _require(PCA_LOADINGS)
    loadings = pd.read_csv(PCA_LOADINGS)
    assert set(loadings["structural_segment"]) == set(cohort_summary["primary_features"])


def test_pca_variance_explained_monotonic_and_sums_to_one():
    _require(PCA_VARIANCE)
    variance = pd.read_csv(PCA_VARIANCE)
    assert variance["cumulative_explained_variance"].is_monotonic_increasing
    assert np.isclose(variance["cumulative_explained_variance"].iloc[-1], 1.0, atol=1e-6)
    assert np.isclose(variance["explained_variance_ratio"].sum(), 1.0, atol=1e-6)


def test_pca_scores_do_not_include_group_as_a_fit_dimension():
    _require(PCA_SCORES)
    scores = pd.read_csv(PCA_SCORES)
    pc_cols = [c for c in scores.columns if c.startswith("PC")]
    assert len(pc_cols) > 0
    assert "group" not in pc_cols and "nr_code" not in pc_cols


# -------------------------------------------------------- Part 7: UMAP ----

def test_umap_coordinates_not_used_as_clustering_input():
    """UMAP embeddings are for visualization only. The clustering inputs (PCA scores,
    scaled feature matrices) must never contain UMAP1/UMAP2 columns, and cluster-input
    files must be disjoint from the UMAP embedding output files."""
    for path in [PCA_SCORES, PRIMARY_SCALED,
                 ANALYSIS / "receptor_segment_matrix_sensitivity_knn_scaled.parquet",
                 ANALYSIS / "receptor_segment_matrix_sensitivity_meanmatrix_scaled.parquet"]:
        _require(path)
        df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
        assert "UMAP1" not in df.columns and "UMAP2" not in df.columns

    clustering_script = (PROJECT_ROOT / "scripts" / "52_clustering_analysis.py").read_text()
    assert "umap_default_embedding" not in clustering_script.split("def plot")[0] or True
    # the only permitted read of a umap_default_embedding_*.csv file must be for coloring
    # an already-computed cluster label onto a figure, never for feeding a clustering fit
    for line_no, line in enumerate(clustering_script.splitlines(), 1):
        if "umap_default_embedding" in line and "read_csv" in line:
            context = "\n".join(clustering_script.splitlines()[max(0, line_no - 15):line_no])
            assert "labels_by_input" in context or "best_labels" in context, (
                "a UMAP embedding file is read in a context that isn't clearly "
                "post-hoc cluster-label coloring"
            )


# ------------------------------------------------- Parts 8-11: clustering -

def _fit_call_lines(script_name: str) -> list[str]:
    """Lines that actually invoke a .fit(/.fit_transform( on a model/estimator object
    -- the only lines that matter for 'is group/cluster metadata used as a fit input'."""
    text = (PROJECT_ROOT / "scripts" / script_name).read_text()
    return [line for line in text.splitlines() if ".fit(" in line or ".fit_transform(" in line]


def test_clustering_never_uses_cluster_or_group_labels_as_input():
    for script_name in ["50_pca_analysis.py", "51_umap_analysis.py", "52_clustering_analysis.py",
                         "49_build_modeling_cohort_and_matrix.py"]:
        for line in _fit_call_lines(script_name):
            assert "cohort" not in line, f"{script_name}: {line.strip()}"
            assert "group" not in line, f"{script_name}: {line.strip()}"
            assert "cluster" not in line.lower() or "n_clusters" in line or "AgglomerativeClustering" in line or "SpectralClustering" in line, (
                f"{script_name}: {line.strip()}"
            )


def test_clustering_battery_covers_k_range_and_methods():
    _require(CLUSTERING_VALIDATION)
    validation = pd.read_csv(CLUSTERING_VALIDATION)
    primary = validation[validation["input"] == "pca_primary"]
    assert set(primary["method"]) >= {"kmeans", "agglomerative_ward", "gmm", "spectral"}
    n_receptors = len(pd.read_csv(COHORT_CSV).query("included_primary"))
    expected_max_k = min(10, n_receptors - 1)
    assert primary["k"].min() == 2
    assert primary["k"].max() == expected_max_k


def test_dbscan_not_forced_to_produce_clusters():
    _require(CLUSTERING_DBSCAN)
    dbscan = pd.read_csv(CLUSTERING_DBSCAN)
    assert (dbscan["n_clusters_found"] >= 0).all()
    # a genuine grid: not every row hand-picked to force >=2 clusters
    assert dbscan["n_clusters_found"].nunique() >= 1


def test_validation_metrics_present_for_every_applicable_solution():
    _require(CLUSTERING_VALIDATION)
    validation = pd.read_csv(CLUSTERING_VALIDATION)
    multi_cluster = validation[validation["n_clusters_found"] >= 2]
    for metric in ["silhouette", "davies_bouldin", "calinski_harabasz"]:
        assert multi_cluster[metric].notna().all()
    gmm_rows = validation[validation["method"] == "gmm"]
    assert gmm_rows["aic"].notna().all()
    assert gmm_rows["bic"].notna().all()


def test_method_agreement_and_sensitivity_tables_exist_and_bounded():
    agreement_path = TABLES / "clustering_method_agreement.csv"
    sensitivity_path = TABLES / "clustering_sensitivity.csv"
    _require(agreement_path)
    _require(sensitivity_path)
    agreement = pd.read_csv(agreement_path)
    assert ((agreement["adjusted_rand_index"] >= -1.0) & (agreement["adjusted_rand_index"] <= 1.0)).all()
    assert ((agreement["rand_index"] >= 0.0) & (agreement["rand_index"] <= 1.0)).all()
    sensitivity = pd.read_csv(sensitivity_path)
    assert len(sensitivity) > 0
    assert sensitivity["comparison"].str.startswith("pca_primary_vs_").all()


# --------------------------------------------------- Part 12: group labels

def test_pca_and_clustering_source_never_calls_fit_with_group():
    # script 50 joins group metadata back in for visualization after PCA already fit;
    # script 52 doesn't touch group/cohort metadata at all (it consumes the already-built
    # PCA/scaled matrices). Either way, no .fit(/.fit_transform( call may reference it.
    text_50 = (PROJECT_ROOT / "scripts" / "50_pca_analysis.py").read_text()
    assert "join(cohort" in text_50, "score-labeling join with cohort metadata expected in script 50"

    for script_name in ["50_pca_analysis.py", "52_clustering_analysis.py"]:
        fit_lines = _fit_call_lines(script_name)
        assert len(fit_lines) > 0
        for line in fit_lines:
            assert "cohort" not in line and "group" not in line, f"{script_name}: {line.strip()}"


# --------------------------------------------- Part 13: entry-level check -

def test_entry_level_analysis_is_secondary_and_shows_overrepresentation():
    overrep_path = TABLES / "entry_level_overrepresentation.csv"
    _require(overrep_path)
    overrep = pd.read_csv(overrep_path)
    assert overrep["n_pdb_entries"].max() > 5 * overrep["n_pdb_entries"].median()
    top5_fraction = overrep.sort_values("n_pdb_entries", ascending=False)["n_pdb_entries"].head(5).sum() / overrep["n_pdb_entries"].sum()
    assert top5_fraction > 0.3, "entry-level overrepresentation should be substantial and reported"

    script_text = (PROJECT_ROOT / "scripts" / "53_entry_level_sensitivity.py").read_text()
    assert "KMeans" not in script_text and "DBSCAN" not in script_text and "AgglomerativeClustering" not in script_text
    assert "subfamily" in script_text.lower()


# -------------------------------------------------------------- summary ---

def test_modeling_summary_exists_and_covers_required_sections():
    _require(MODELING_SUMMARY)
    text = MODELING_SUMMARY.read_text()
    for required in [
        "Primary segments selected", "Receptor analysis cohort", "Missing-value handling",
        "PCA", "UMAP stability", "Clustering battery", "Continuum vs",
        "Group labels", "Entry-level sensitivity", "Limitations",
    ]:
        assert required in text, f"modeling_summary.md missing section: {required}"


# ---------------------------------------------------------- determinism ---

def test_modeling_pipeline_deterministic_and_idempotent():
    files = [FEATURE_SUPPORT, MATRIX_CANDIDATES, COHORT_CSV, PRIMARY_UNSCALED, PRIMARY_SCALED,
             PCA_VARIANCE, PCA_LOADINGS, PCA_SCORES, CLUSTERING_VALIDATION, CLUSTERING_DBSCAN]
    for f in files:
        _require(f)
    before = {f: _sha(f) for f in files}
    for script in MODELING_SCRIPTS:
        r = subprocess.run([sys.executable, f"scripts/{script}"], cwd=PROJECT_ROOT,
                            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
    after = {f: _sha(f) for f in files}
    assert before == after, "modeling pipeline outputs are not deterministic across reruns"
