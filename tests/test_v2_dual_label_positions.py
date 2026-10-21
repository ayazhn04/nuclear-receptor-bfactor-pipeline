"""
Dual-label safety of the v2 professor map (pre-Phase 2).

Two positions lie inside two professor intervals: PXR (O75469) 261 in H3 and H4, and
CAR (Q14994) 344 in H11' and H12. Each must stay a member of BOTH elements in every
membership/aggregation table, with no other position duplicated anywhere. A joined
string ("H3|H4") is allowed only as a display / one-row lookup value, never as an
element.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils import canonical_helix as ch  # noqa: E402

ANALYSIS_V2 = PROJECT_ROOT / "data" / "rebuild_v2" / "analysis"
SEGMAP = ANALYSIS_V2 / "nr_lbd_structural_segment_map.csv"
INTERVALS = PROJECT_ROOT / "data" / "canonical_remap" / "nr_lbd_canonical_helix_intervals_professor_2026.csv"
KEY = ["uniprot_id", "canonical_uniprot_position"]
DUAL = {("O75469", 261): ("H3", "H4"), ("Q14994", 344): ("H11'", "H12")}

pytestmark = pytest.mark.skipif(not SEGMAP.exists(), reason="v2 foundation not built")


@pytest.fixture(scope="module")
def seg_map():
    return pd.read_csv(SEGMAP, keep_default_na=False)


@pytest.fixture(scope="module")
def memberships(seg_map):
    return ch.position_memberships(seg_map[KEY + ["structural_segment"]].itertuples(index=False))


def test_only_the_two_intended_positions_have_two_rows(seg_map):
    counts = seg_map.groupby(KEY).size()
    multi = counts[counts > 1]
    assert {k: int(v) for k, v in multi.items()} == {k: 2 for k in DUAL}
    for (uid, pos), labels in DUAL.items():
        rows = seg_map[(seg_map["uniprot_id"] == uid) & (seg_map["canonical_uniprot_position"] == pos)]
        assert tuple(rows["structural_segment"]) == labels
        assert set(rows["structural_class"]) == {"HELIX"}


def test_memberships_are_exactly_both_elements(memberships, seg_map):
    multi = {k: v for k, v in memberships.items() if len(v) > 1}
    assert multi == DUAL
    # one membership per HELIX row: nothing lost, nothing added
    assert sum(len(v) for v in memberships.values()) == int((seg_map["structural_class"] == "HELIX").sum())
    assert len(memberships) == seg_map.loc[seg_map["structural_segment"] != "", KEY].drop_duplicates().shape[0]


def test_display_label_is_never_an_element(memberships, seg_map):
    assert ch.display_label(DUAL[("O75469", 261)]) == "H3|H4"
    assert ch.display_label(DUAL[("Q14994", 344)]) == "H11'|H12"
    assert not seg_map["structural_segment"].str.contains(ch.DISPLAY_JOIN, regex=False).any()
    for labels in memberships.values():
        for label in labels:
            ch.element_rank(label)  # every membership is a real professor element


def test_element_sizes_equal_professor_interval_lengths(seg_map):
    iv = pd.read_csv(INTERVALS, keep_default_na=False)
    prof = iv.assign(n=iv["interval_length"].astype(int)).groupby(["project_uniprot_id", "analysis_label"])["n"].sum()
    helix = seg_map[seg_map["structural_class"] == "HELIX"]
    mapped = helix.groupby(["uniprot_id", "structural_segment"]).size()
    prof.index = prof.index.set_names(mapped.index.names)
    assert mapped.sort_index().astype(int).equals(prof.sort_index().astype(int))
    # both dual positions are counted once in each of their two elements
    for (uid, pos), labels in DUAL.items():
        for label in labels:
            members = helix[(helix["uniprot_id"] == uid) & (helix["structural_segment"] == label)]
            assert (members["canonical_uniprot_position"] == pos).sum() == 1


@pytest.mark.parametrize("name,key", [
    ("pdb_segment_flexibility_features.parquet", ["uniprot_id", "pdb_id", "structural_segment"]),
    ("pdb_segment_coverage.parquet", ["uniprot_id", "pdb_id", "structural_segment"]),
    ("receptor_segment_flexibility_profiles.parquet", ["uniprot_id", "structural_segment"]),
])
def test_element_tables_have_no_duplicated_rows(name, key):
    df = pd.read_parquet(ANALYSIS_V2 / name)
    assert not df.duplicated(key).any()
    assert not df["structural_segment"].str.contains(ch.DISPLAY_JOIN, regex=False).any()


def test_both_elements_present_in_profiles():
    prof = pd.read_parquet(ANALYSIS_V2 / "receptor_segment_flexibility_profiles.parquet")
    have = set(zip(prof["uniprot_id"], prof["structural_segment"]))
    for (uid, _), labels in DUAL.items():
        for label in labels:
            assert (uid, label) in have


def test_legacy_lookup_pattern_fails_and_replacement_is_one_row(seg_map, memberships):
    """Scripts 65 and 79 look a position up with Series.get on a non-unique index.
    At a dual position that returns a Series and `if seg:` raises, so both scripts
    must use the membership lookup in Phase 2."""
    legacy = seg_map.set_index(KEY)["structural_segment"]
    with pytest.raises(ValueError):
        bool(legacy.get(("O75469", 261), ""))

    # script 65: one residue-test row per position, joined label for display only
    positions = seg_map[KEY].drop_duplicates()
    annotated = positions.assign(structural_segment=[
        ch.display_label(memberships.get((u, int(p)), ())) for u, p in positions.itertuples(index=False)
    ])
    assert len(annotated) == len(positions)
    assert not annotated.duplicated(KEY).any()
    got = annotated.set_index(KEY)["structural_segment"]
    assert got[("O75469", 261)] == "H3|H4" and got[("Q14994", 344)] == "H11'|H12"

    # script 79: a site residue at a dual position marks EACH of its elements
    site_associated = set()
    for key in DUAL:
        site_associated.update(memberships.get(key, ()))
    assert site_associated == {"H3", "H4", "H11'", "H12"}
