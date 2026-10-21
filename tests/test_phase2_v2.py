"""
Phase 2 (corrected v2 rerun) checks: the S1 scope helper, the S2 matching context,
the logged sandbox patches, and the exported v2 results.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils import apo_holo_scope, matching_context  # noqa: E402

V2 = PROJECT_ROOT / "data" / "rebuild_v2"
V2A = V2 / "analysis"
V2R = PROJECT_ROOT / "reports" / "rebuild_v2"
V2T = V2R / "tables"
SEGMAP = V2A / "nr_lbd_structural_segment_map.csv"
SS = V2A / "receptor_dssp_consensus.parquet"

needs_v2 = pytest.mark.skipif(not (V2R / "phase2_export_manifest.json").exists(), reason="Phase 2 not exported")


# ---- S1 -------------------------------------------------------------------

@pytest.mark.skipif(not (V2 / "processed" / "structures.csv").exists(), reason="v2 foundation not built")
def test_trace_scope_uses_only_all_trace_entries(tmp_path):
    root = V2
    s = pd.read_csv(root / "processed" / "structures.csv", dtype={"pdb_id": str}, keep_default_na=False, low_memory=False)
    coord = s[s["coordinate_data_status"] == "COORDINATE_DATA_AVAILABLE"]
    cats = coord.groupby(["uniprot_id", "pdb_id"])["lbd_mapping_category"].agg(set)
    expected = {k for k, v in cats.items() if v == {"TRACE_LBD"}}
    assert ("P04150", "6YMO") in expected and ("P62508", "6Y58") in expected
    # RORgamma entries with an LBD entity plus a peptide entity are NOT excluded
    assert ("P51449", "4WQP") not in expected
    df = pd.DataFrame({"uniprot_id": ["P04150", "P51449"], "pdb_id": ["6YMO", "4WQP"], "ligand_state": ["APO", "AMBIGUOUS"]})
    # the helper reads <project>/data/processed/structures.csv: give it the v2 table in a temp layout
    (tmp_path / "data" / "processed").mkdir(parents=True)
    (tmp_path / "data" / "processed" / "structures.csv").write_bytes((root / "processed" / "structures.csv").read_bytes())
    out = apo_holo_scope.apply(df, tmp_path)
    assert list(out["ligand_state"]) == [apo_holo_scope.EXCLUDED_STATE, "AMBIGUOUS"]
    assert list(df["ligand_state"]) == ["APO", "AMBIGUOUS"]  # input untouched


@needs_v2
def test_s1_applied_symmetrically_in_eligibility():
    e = pd.read_csv(V2T / "apo_holo_receptor_eligibility.csv")
    gr = e[e["uniprot_id"] == "P04150"].iloc[0]
    assert gr["n_APO"] == 0 and not gr["eligible_exploratory_3_3"]
    assert int(e["n_excluded_TRACE_LBD"].sum()) > 0
    feats = pd.read_parquet(V2A / "pdb_segment_flexibility_features.parquet")
    tests = pd.read_csv(V2T / "apo_holo_segment_tests.csv")
    # every tested receptor is eligible; no excluded entry can be counted on either side
    assert set(tests["uniprot_id"]) <= set(e.loc[e["eligible_exploratory_3_3"], "uniprot_id"])
    assert len(feats) > 0


# ---- S2 -------------------------------------------------------------------

@pytest.mark.skipif(not (SEGMAP.exists() and SS.exists()), reason="v2 inputs missing")
def test_matching_frame():
    seg = pd.read_csv(SEGMAP, keep_default_na=False)
    ss = pd.read_parquet(SS)
    f = matching_context.build_matching_frame(seg, ss)
    mapped = seg[seg["structural_class"] != "NO_PROFESSOR_MAPPING"]
    # every mapped span position appears; dual positions twice; unmapped receptors never
    assert len(f) == len(mapped)
    assert set(f["uniprot_id"]) == set(mapped["uniprot_id"])
    assert f.duplicated(["uniprot_id", "canonical_uniprot_position"]).sum() == 2
    inter = f[f["structural_segment"].str.startswith(matching_context.INTER_PREFIX)]
    assert len(inter) == int((mapped["structural_class"] == "INTER_HELIX_DERIVED").sum())
    assert inter["structural_segment"].str.match(r"^INTER:[^-]+-[^-]+$").all()
    # local class is independent and NO_OBSERVATION is never a class
    assert set(f["structural_class"]) <= {"HELIX", "SHEET", "COIL", ""}
    obs = set(zip(ss["uniprot_id"], ss["canonical_uniprot_position"]))
    blank = f[f["structural_class"] == ""]
    assert not any(k in obs for k in zip(blank["uniprot_id"], blank["canonical_uniprot_position"]))
    # RORgamma gap 467-469 is an inter-element stretch inside H10/H11, never part of it
    ror = f[(f["uniprot_id"] == "P51449") & f["canonical_uniprot_position"].between(467, 469)]
    assert set(ror["structural_segment"]) == {"INTER:H10/H11-H10/H11"}


# ---- patches and outputs ----------------------------------------------------

@needs_v2
def test_patch_log_complete():
    log = json.loads((V2R / "sandbox_patches_phase2.json").read_text())
    assert log and all(p["status"] in ("APPLIED", "ALREADY_APPLIED") for p in log)
    ids = {p["adaptation"] for p in log}
    for a in ["A02", "A06", "A07", "A08", "A09", "A10", "A11", "A13", "A14", "A15", "A16", "A17", "A18", "A23", "A26", "A27"]:
        assert any(a in i for i in ids), a


@needs_v2
def test_phase2_validation_all_pass():
    v = json.loads((V2R / "phase2_validation_v2.json").read_text())
    assert v["n_fail"] == 0, [c for c in v["checks"] if c["status"] == "FAIL"]


@needs_v2
def test_imputed_cells_documented_and_not_zero():
    imp = pd.read_csv(V2T / "modeling_imputed_cells.csv")
    mask = pd.read_parquet(V2A / "receptor_segment_matrix_primary_imputation_mask.parquet").set_index("uniprot_id")
    assert len(imp) == int(mask.values.sum())
    assert set(imp["value_status"]) == {"IMPUTED_PLACEHOLDER_NOT_AN_OBSERVATION"}
    med = pd.read_parquet(V2A / "receptor_segment_feature_matrix_median.parquet").set_index("uniprot_id")
    for _, r in imp.iterrows():
        assert pd.isna(med.loc[r["uniprot_id"], r["feature"]])


@needs_v2
def test_s3_primary_and_sensitivity_are_separate_files():
    prim = V2T / "functional_site_superfamily_matched_reanalysis.csv"
    sens = V2T / "functional_site_superfamily_matched_reanalysis_sensitivity_excl_review_required.csv"
    assert prim.exists() and sens.exists()
    defs = pd.read_parquet(V2A / "current_functional_site_definitions.parquet")
    prof = defs[defs["definition_source"] == "PROFESSOR_GROUND_TRUTH"]
    # professor statuses are kept verbatim, never relabelled HIGH
    assert set(prof["assignment_confidence"]) <= {"published_boundaries", "structural_homology", "review_required"}


@needs_v2
def test_no_v3_columns_in_feature_support():
    sup = pd.read_csv(V2T / "modeling_feature_support.csv")
    assert not any(c.startswith("confidence_") for c in sup.columns)
    assert "beta_sheet" not in set(sup["structural_segment"])


@needs_v2
def test_discovery_metric_not_mislabelled_as_project_receptors():
    fp = pd.read_csv(V2T / "final_project_numbers.csv").set_index("metric")["value"]
    assert "project_receptors" not in fp.index
    assert int(fp["project_receptors_defined"]) == 48
    assert int(fp["receptors_with_any_discovery_candidate"]) == pd.read_csv(
        V2 / "processed" / "structures.csv", low_memory=False)["uniprot_id"].nunique()


@pytest.mark.skipif(not (V2T / "trace_lbd_site_sensitivity_comparison.csv").exists(), reason="sensitivity not run")
def test_trace_lbd_site_sensitivity_changes_no_headline():
    c = pd.read_csv(V2T / "trace_lbd_site_sensitivity_comparison.csv")
    head = c[c["quantity"].str.startswith(("classification", "superfamily", "matched superfamily",
                                           "within-receptor FDR", "matched-background FDR"))]
    assert len(head) and int(head["changed"].sum()) == 0
