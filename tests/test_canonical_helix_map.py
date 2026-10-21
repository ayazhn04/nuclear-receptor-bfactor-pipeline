"""
Professor canonical helix map tests -- FINAL policy (professor decision
2026-10-20): the professor table is ground truth and is transferred exactly.

Covers source preservation, interval expansion, identity/TRalpha (P10827-2 as
final numbering reference), use of every professor label incl. review_required,
both labels at overlapping professor intervals, derived inter-helix positions
(never loop names, never v3), receptors without mapping, analysis-span
containment, the audit-history tables (v3 comparison, endpoint concordance,
pre-build LBD window audit) and determinism / frozen-file protection.
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils import canonical_helix as ch  # noqa: E402

SOURCE_DIR = PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table"
SOURCE_CSV = SOURCE_DIR / "NR_canonical_helix_boundaries.csv"
SOURCE_DOCX = SOURCE_DIR / "Assignment of conventional LBD helix identities across hNRs.docx"
ANALYSIS = PROJECT_ROOT / "data" / "analysis"          # frozen legacy layer (read-only)
STAGE_DATA = PROJECT_ROOT / "data" / "canonical_remap"
TABLES = PROJECT_ROOT / "reports" / "canonical_remap"
SPAN = PROJECT_ROOT / "data" / "rebuild_v2" / "inputs" / "analysis_span_v2.csv"
INTERVALS = STAGE_DATA / "nr_lbd_canonical_helix_intervals_professor_2026.csv"
MAP = STAGE_DATA / "nr_lbd_canonical_helix_map_professor_2026.csv"
SUMMARY = STAGE_DATA / "canonical_helix_map_summary.json"
IDENTITY = TABLES / "canonical_helix_identity_isoform_audit.csv"
COVERAGE = TABLES / "canonical_helix_receptor_coverage_audit.csv"
CONFLICTS = TABLES / "canonical_helix_conflicts.csv"
V3_ELEMENTS = TABLES / "canonical_helix_v3_comparison_elements.csv"
V3_RECEPTOR = TABLES / "canonical_helix_v3_comparison_by_receptor.csv"
ENDPOINTS = TABLES / "canonical_helix_endpoint_concordance.csv"
ENDPOINTS_BY_RECEPTOR = TABLES / "canonical_helix_endpoint_concordance_by_receptor.csv"
LBD_WINDOW = TABLES / "lbd_window_interval_audit_v1.csv"
LBD_WINDOW_RECEPTORS = TABLES / "lbd_window_receptor_summary_v1.csv"
BUILD_OUTPUTS = [INTERVALS, MAP, SUMMARY, IDENTITY, COVERAGE, CONFLICTS, V3_ELEMENTS, V3_RECEPTOR]
ENDPOINT_OUTPUTS = [ENDPOINTS, ENDPOINTS_BY_RECEPTOR]
HAVE_RAW = (PROJECT_ROOT / "data" / "raw" / "mmcif").is_dir() and (PROJECT_ROOT / "data" / "raw" / "api" / "uniprot").is_dir()

SOURCE_SHA256 = "ae1fb1005e8e4a10d7cf568ac3657e5bf692c62f5eac7572117b558b92e4268c"
DOCX_SHA256 = "6d12b046a40a856c6ba81cc8004ff567162bfadf5543619f859e3154fec95758"
V3_MAP_SHA256 = "b1b18be7d6b2f6450da6c494cb19ef168a0d0f92ec050295ffa4f3c92d8faa2a"
FROZEN_SHA256 = {
    "data/processed/structures.csv": "b9aacc22ea0fcf803b576e883d5b4d0e4b3c86d2b9239a6e3eee385500142246",
    "data/processed/lbd_residue_bfactors_all.parquet": "28cc0f7feff59d27cf6b724becad23391068c6ed74181db3a6b9923750f81a27",
    "data/analysis/normalized_lbd_bfactors_all.parquet": "620fe63a74a10953a1c6b2b2d216d2d0e5290f870d1a686037aa78ea66f9b904",
    "data/analysis/pdb_position_bfactor_profiles.parquet": "6dfa602e05a0fdf41d3ee117f10dbc0ded87aaf0f1e16f76edcd858d8010beee",
    "data/analysis/receptor_position_flexibility_profiles.parquet": "05e3ba6bd2ab76493376e197b401264fddc817bd6aaa5fbde4c32f832aa2b222",
}
NOT_IN_SOURCE = {"NR1F2": "RORβ", "NR2C1": "TR2", "NR2F1": "COUP-TFI", "NR6A1": "GCNF", "NR0B1": "DAX-1",
                 "NR0B2": "SHP"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), PROJECT_ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def source():
    return pd.read_csv(SOURCE_CSV, encoding=ch.SOURCE_ENCODING, dtype=str, keep_default_na=False)


@pytest.fixture(scope="module")
def intervals():
    return pd.read_csv(INTERVALS, keep_default_na=False)


@pytest.fixture(scope="module")
def cmap():
    return pd.read_csv(MAP, keep_default_na=False, dtype={"canonical_uniprot_position": int})


# --- source preservation --------------------------------------------------------


def test_source_files_preserved_byte_for_byte():
    assert sha256(SOURCE_CSV) == SOURCE_SHA256 and sha256(SOURCE_DOCX) == DOCX_SHA256


def test_source_shape_and_status_counts(source):
    assert list(source.columns) == ch.SOURCE_COLUMNS and len(source) == 434
    assert source["nr_code"].nunique() == 42
    assert source["assignment_status"].value_counts().to_dict() == {
        "structural_homology": 354, "review_required": 53, "published_boundaries": 27}


def test_every_source_row_preserved_verbatim(source, intervals):
    iv = intervals[ch.SOURCE_COLUMNS].astype(str).reset_index(drop=True)
    pd.testing.assert_frame_equal(iv, source.astype(str).reset_index(drop=True))
    assert intervals["source_row_number"].tolist() == list(range(1, 435))


# --- policy -----------------------------------------------------------------------


def test_final_policy_uses_every_professor_label(intervals):
    assert (intervals["analysis_label"] == intervals["canonical_helix"]).all()
    assert (intervals["label_policy"] == ch.FINAL_LABEL_POLICY).all()
    rr = intervals[intervals["assignment_status"] == "review_required"]
    assert len(rr) == 53 and (rr["review_required_flag"].astype(str) == "True").all()  # metadata only


def test_vocabulary_and_expansion():
    assert ch.expand_interval(5, 7) == [5, 6, 7]
    with pytest.raises(ValueError):
        ch.expand_interval(8, 7)
    with pytest.raises(ValueError):
        ch.element_rank("H13")
    with pytest.raises(ValueError):
        ch.status_policy("confirmed")
    assert ch.classify_position_final([]) == ch.FPOS_INTER
    assert ch.classify_position_final([], receptor_mapped=False) == ch.FPOS_NO_MAPPING
    a = {"canonical_helix": "H3"}
    assert ch.classify_position_final([a]) == ch.FPOS_ELEMENT
    assert ch.classify_position_final([a, {"canonical_helix": "H4"}]) == ch.FPOS_MULTI


def test_element_position_sets_equal_professor_union(intervals, cmap):
    for (uid, el), g in intervals.groupby(["project_uniprot_id", "canonical_helix"]):
        exp = set().union(*[set(range(a, b + 1)) for a, b in zip(g["start"], g["end"])])
        m = cmap[(cmap["uniprot_id"] == uid) & cmap["canonical_elements"].str.split("|").apply(lambda xs: el in xs)]
        assert set(m["canonical_uniprot_position"]) == exp, (uid, el)


def test_overlapping_professor_intervals_carry_both_labels(cmap):
    multi = cmap[cmap["position_status"] == ch.FPOS_MULTI]
    got = {(n, p, h) for n, p, h in zip(multi["nr_code"], multi["canonical_uniprot_position"], multi["canonical_elements"])}
    assert got == {("NR1I2", 261, "H3|H4"), ("NR1I3", 344, "H11'|H12")}
    conf = pd.read_csv(CONFLICTS, keep_default_na=False)
    ov = conf[conf["issue_type"] == "PROFESSOR_INTERVALS_OVERLAP"]
    assert len(ov) == 2 and (ov["final_policy_handling"] == "BOTH_PROFESSOR_LABELS_APPLIED_AT_SHARED_POSITIONS").all()


def test_split_element_used_exactly(cmap):
    rorg = cmap[cmap["nr_code"] == "NR1F3"].set_index("canonical_uniprot_position")
    assert (rorg.loc[467:469, "position_status"] == ch.FPOS_INTER).all()
    assert rorg.loc[466, "canonical_elements"] == rorg.loc[470, "canonical_elements"] == "H10/H11"


def test_inter_helix_positions_are_derived_never_loops(cmap):
    inter = cmap[cmap["position_status"] == ch.FPOS_INTER]
    assert (inter["canonical_elements"] == "").all()
    for c in ["canonical_elements", "position_status"]:
        assert not cmap[c].str.contains("loop|beta_sheet", case=False).any()
    labels = {e for s in cmap["canonical_elements"] for e in s.split("|") if e}
    assert labels <= set(ch.ELEMENT_ORDER)


def test_receptor_specific_elements_preserved(source, cmap):
    for name in sorted(ch.RECEPTOR_SPECIFIC_ELEMENTS):
        s = set(source[source["canonical_helix"] == name]["nr_code"])
        m = set(cmap[cmap["canonical_elements"].str.split("|").apply(lambda xs: name in xs)]["nr_code"])
        assert s == m, name


def test_pparg_omega_region_between_h2prime_and_h3(cmap):
    p = cmap[cmap["nr_code"] == "NR1C3"].set_index("canonical_uniprot_position")
    assert (p.loc[279:291, "canonical_elements"] == "H2'").all()
    assert (p.loc[292:303, "position_status"] == ch.FPOS_INTER).all()
    assert (p.loc[292:303, "inter_helix_flank_before"] == "H2'").all()
    assert (p.loc[292:303, "inter_helix_flank_after"] == "H3").all()


# --- identity / TRalpha / span ------------------------------------------------


def test_identity_audit():
    idt = pd.read_csv(IDENTITY, keep_default_na=False)
    assert len(idt) == 48
    assert idt["identity_status"].value_counts().to_dict() == {
        "MATCH": 41, "MATCH_PROFESSOR_ISOFORM_IS_FINAL_REFERENCE": 1, "NOT_IN_SOURCE": 6}
    tra = idt.set_index("nr_code").loc["NR1A1"]
    assert tra["final_reference_accession"] == "P10827-2" and int(tra["reference_sequence_length"]) == 410


def test_tralpha_final_numbering_is_p10827_2(cmap):
    assert ch.reference_accession("P10827") == "P10827-2"
    assert ch.final_crosswalk("P10827-2", "P10827").identical_range is None
    with pytest.raises(ValueError):
        ch.final_crosswalk("P10827", "P10827")
    tra = cmap[cmap["nr_code"] == "NR1A1"]
    assert set(tra["reference_accession"]) == {"P10827-2"}
    h12 = tra[tra["canonical_elements"] == "H12"]["canonical_uniprot_position"]
    assert (h12.min(), h12.max()) == (398, 406)  # professor TRalpha1 H12, used exactly
    iso1 = ch.read_fasta_sequence(SOURCE_DIR / "uniprot_isoform_sequences" / "P10827-2.fasta")
    iso2 = ch.read_fasta_sequence(SOURCE_DIR / "uniprot_isoform_sequences" / "P10827.fasta")
    assert iso1[:370] == iso2[:370] and ch.first_sequence_difference(iso1, iso2) == 371


def test_every_interval_inside_final_analysis_span(intervals):
    assert (intervals["fully_inside_analysis_span"].astype(str) == "True").all()
    span = pd.read_csv(SPAN).set_index("uniprot_id")
    for code, uid in [("NR1C1", "Q07869"), ("NR1C2", "Q03181"), ("NR3C1", "P04150"), ("NR3C2", "P08235"),
                      ("NR3C3", "P06401"), ("NR3C4", "P10275")]:
        g = intervals[intervals["nr_code"] == code]
        assert span.loc[uid, "analysis_span_start"] <= g["start"].min()
        assert span.loc[uid, "analysis_span_end"] >= g["end"].max()
    un = span[~span["professor_mapping"]]  # unmapped receptors keep the UniProt/PROSITE span
    assert ((un["analysis_span_start"] == un["prosite_lbd_start_v1"]) & (un["analysis_span_end"] == un["prosite_lbd_end_v1"])).all()


def test_coverage_42_mapped_6_unavailable():
    cov = pd.read_csv(COVERAGE, keep_default_na=False)
    un = cov[cov["professor_mapping_exists"].astype(str) == "False"]
    assert dict(zip(un["nr_code"], un["common_name"])) == NOT_IN_SOURCE
    assert (un["mapping_status"] == "CANONICAL_MAPPING_UNAVAILABLE").all()
    m = cov[cov["professor_mapping_exists"].astype(str) == "True"]
    assert len(m) == 42 and m["n_review_required"].sum() == 53


def test_unmapped_receptors_have_no_labels(cmap):
    u = cmap[cmap["nr_code"].isin(NOT_IN_SOURCE)]
    assert set(u["position_status"]) == {ch.FPOS_NO_MAPPING} and (u["canonical_elements"] == "").all()


# --- audit history (does not influence labels) ----------------------------------


def test_v3_comparison_is_audit_history_only():
    r = pd.read_csv(V3_RECEPTOR, keep_default_na=False)
    assert (r["note"] == "AUDIT HISTORY ONLY - does not influence any final label").all()
    rec = r.set_index("nr_code")["professor_elements_not_named_in_v3"].str.split()
    assert {"H7", "H8", "H9"} <= set(rec["NR3A1"]) and {"H6", "H7", "H8"} <= set(rec["NR1I1"])


def test_endpoints_agree_with_deposited_helices_under_independent_alignment():
    s = pd.read_csv(ENDPOINTS_BY_RECEPTOR, keep_default_na=False)
    chk = s["n_endpoints_checkable_alignment"].astype(int).sum()
    assert chk >= 700 and s["n_endpoints_exact_alignment"].astype(int).sum() / chk >= 0.99


def test_lbd_window_audit_history_reconciles():
    a = pd.read_csv(LBD_WINDOW, keep_default_na=False)
    assert len(a) == 434
    assert a["window_class"].value_counts().to_dict() == {
        "FULLY_INSIDE": 407, "PARTIALLY_TRUNCATED_N_TERMINAL": 14, "PARTIALLY_TRUNCATED_C_TERMINAL": 7,
        "FULLY_OUTSIDE": 4, "ISOFORM_POSITIONS_NOT_REPRESENTABLE_IN_PROJECT_ACCESSION": 2}
    r = pd.read_csv(LBD_WINDOW_RECEPTORS, keep_default_na=False).set_index("nr_code")
    assert set(r[r["recommendation"] == "B_EXTEND_LBD_TO_CANONICAL_HELIX_BOUNDARY"].index) == {
        "NR1C1", "NR1C2", "NR3C1", "NR3C2", "NR3C3", "NR3C4"}


# --- determinism / frozen -------------------------------------------------------


def test_rebuild_is_byte_identical(tmp_path):
    load_script("97_build_canonical_helix_map.py").build(tmp_path)
    for p in BUILD_OUTPUTS:
        rel = p.relative_to(PROJECT_ROOT)
        assert (tmp_path / rel).read_bytes() == p.read_bytes(), rel


@pytest.mark.skipif(not HAVE_RAW, reason="local mmCIF / UniProt caches not present")
def test_endpoint_audit_rebuild_is_byte_identical(tmp_path):
    load_script("98_audit_canonical_helix_endpoints.py").build(PROJECT_ROOT / "data" / "raw" / "mmcif", tmp_path)
    for p in ENDPOINT_OUTPUTS:
        rel = p.relative_to(PROJECT_ROOT)
        assert (tmp_path / rel).read_bytes() == p.read_bytes(), rel


@pytest.mark.skipif(not HAVE_RAW, reason="local mmCIF / UniProt caches not present")
def test_lbd_window_audit_rebuild_is_byte_identical(tmp_path):
    load_script("100_audit_prebuild_lbd_window.py").build(tmp_path)
    for name in ["lbd_window_interval_audit_v1.csv", "lbd_window_receptor_summary_v1.csv", "lbd_window_helix_summary_v1.csv"]:
        assert (tmp_path / "reports" / "canonical_remap" / name).read_bytes() == (TABLES / name).read_bytes(), name


def test_historical_v3_map_unchanged():
    assert sha256(ANALYSIS / "nr_lbd_structural_segment_map.csv") == V3_MAP_SHA256


@pytest.mark.parametrize("rel", sorted(FROZEN_SHA256))
def test_frozen_foundation_unchanged(rel):
    assert sha256(PROJECT_ROOT / rel) == FROZEN_SHA256[rel]
