"""
Residue-numbering correction layer v1 tests: deletion/insertion-aware
alignment, the VDR delta165-215 case, isoforms, point mutation vs numbering
offset, withholding of ambiguous residues, value-identity of unaffected
mappings, determinism, rebuild-input regions, count reconciliation, and
protection of the frozen inputs, the professor source and the historical
frozen identity audit.

Synthetic cases need nothing on disk. Recompute checks additionally need the
gitignored mmCIF / UniProt caches and are skipped without them.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils import sequence_numbering as sn  # noqa: E402
from scripts.utils.residue_mapping import build_entity_to_uniprot_map  # noqa: E402

CROSSWALK = PROJECT_ROOT / "data" / "canonical_remap" / "residue_numbering_crosswalk_v1.parquet"
REGIONS = PROJECT_ROOT / "data" / "canonical_remap" / "entity_uniprot_regions_corrected_v1.csv"
STAGE = PROJECT_ROOT / "reports" / "canonical_remap"
INSTANCES = STAGE / "residue_numbering_instance_audit_v1.csv"
SUMMARY = STAGE / "residue_numbering_summary_v1.json"
VDR = STAGE / "vdr_entry_numbering_audit_v1.csv"
IMPACT = STAGE / "numbering_impact_by_receptor_v1.csv"
STAGE_IDENTITY = STAGE / "stage_receptor_identity_audit_v1.csv"
ENTITY = STAGE / "numbering_entity_consistency_v1.csv"
FROZEN_RESIDUES = PROJECT_ROOT / "data" / "processed" / "lbd_residue_bfactors_all.parquet"
HAVE_RAW = (PROJECT_ROOT / "data" / "raw" / "mmcif").is_dir() and (PROJECT_ROOT / "data" / "raw" / "api" / "uniprot").is_dir()

FROZEN_SHA256 = {
    "data/processed/structures.csv": "b9aacc22ea0fcf803b576e883d5b4d0e4b3c86d2b9239a6e3eee385500142246",
    "data/processed/lbd_residue_bfactors_all.parquet": "28cc0f7feff59d27cf6b724becad23391068c6ed74181db3a6b9923750f81a27",
    "data/analysis/normalized_lbd_bfactors_all.parquet": "620fe63a74a10953a1c6b2b2d216d2d0e5290f870d1a686037aa78ea66f9b904",
    "data/analysis/pdb_position_bfactor_profiles.parquet": "6dfa602e05a0fdf41d3ee117f10dbc0ded87aaf0f1e16f76edcd858d8010beee",
    "data/analysis/receptor_position_flexibility_profiles.parquet": "05e3ba6bd2ab76493376e197b401264fddc817bd6aaa5fbde4c32f832aa2b222",
    "data/analysis/nr_lbd_structural_segment_map.csv": "b1b18be7d6b2f6450da6c494cb19ef168a0d0f92ec050295ffa4f3c92d8faa2a",
    "data/manifests/rcsb_sequence_coordinate_alignments.csv": None,  # filled from git below
    "reports/tables/final_receptor_identity_audit.csv": None,
    "data/sources/professor_canonical_helix_table/NR_canonical_helix_boundaries.csv":
        "ae1fb1005e8e4a10d7cf568ac3657e5bf692c62f5eac7572117b558b92e4268c",
    "data/sources/professor_canonical_helix_table/Assignment of conventional LBD helix identities across hNRs.docx":
        "6d12b046a40a856c6ba81cc8004ff567162bfadf5543619f859e3154fec95758",
}
LEGACY_SCANNED_DIRS = ["data/analysis", "data/manifests", "data/processed", "reports/tables"]

# A deterministic 220-residue reference (VDR 118-337) for synthetic cases.
REF = ("DSLRPKLSEEQQRIIAILLDAHHKTYDPTYSDFCQFRPPVRVNDGGGSHPSRPNSRHTPSFSGDSSSSCSDHCITSSDMMDSSSFSNLDLSEEDSDDPSVTLEL"
       "SQLSMLPHLADLVSYSIQKVIGFAKMIPGFRDLTSEDQIVLLKSSAIEVIMLRSNESFTMDDMSWTCGNQDYKYRVSDVTKAGHSLELIEPLIKFQVGLKKLN")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_head_sha256(rel: str) -> str:
    blob = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=PROJECT_ROOT, capture_output=True, check=True).stdout
    return hashlib.sha256(blob).hexdigest()


def run_calls(query: str, target: str):
    res = sn.align_unique(query, target)
    calls, summ = sn.call_residues(query, target, res)
    return [c.target_index for c in calls], [c.reason for c in calls], summ, res


@pytest.fixture(scope="module")
def cw():
    return pd.read_parquet(CROSSWALK)


@pytest.fixture(scope="module")
def inst():
    return pd.read_csv(INSTANCES, keep_default_na=False)


# --- synthetic alignment behaviour ---------------------------------------------


def test_deletion_aware_alignment():
    chain = REF[:50] + REF[101:200]          # 51-residue internal deletion
    pos, _, summ, res = run_calls(chain, REF)
    assert summ["sufficient_evidence"]
    gaps = sn.gaps_in_alignment(res.blocks)
    assert [g["type"] for g in gaps] == ["DELETION_IN_CHAIN"]
    for i, p in enumerate(pos):
        if p is not None:  # never a shifted placement: either the true position or withheld
            assert p == (i if i < 50 else i + 51)
    assert sum(p is None for p in pos) <= 2  # only junction residues may be withheld


def test_insertion_aware_alignment():
    ins = "GSGSGSGSGS"
    chain = REF[:60] + ins + REF[60:150]
    pos, _, _, res = run_calls(chain, REF)
    assert any(g["type"] == "INSERTION_IN_CHAIN" for g in sn.gaps_in_alignment(res.blocks))
    for i, p in enumerate(pos):
        if i < 60:
            assert p in (i, None)
        elif i < 60 + len(ins):
            assert p is None or chain[i] == REF[p]  # inserted residues never forced onto a mismatch
        else:
            assert p in (i - len(ins), None)
    assert sum(p is None for p in pos[60 + len(ins):]) <= 2


def test_point_mutations_are_substitutions_not_offsets():
    chain = list(REF[:150])
    for i in (30, 75, 120):
        chain[i] = "W" if chain[i] != "W" else "A"
    chain = "".join(chain)
    pos, reasons, _, _ = run_calls(chain, REF)
    assert pos == list(range(150))
    assert all(r == "MAPPED" for r in reasons)


def test_ambiguous_junction_residue_is_withheld():
    # residue at the junction equals the residue on both sides of the deletion
    left, right = REF[:40] + "K", "K" + REF[100:160]
    target = left + REF[41:99] + right
    chain = left + right[1:]
    pos, reasons, _, res = run_calls(chain, target)
    assert res.n_optimal >= 2
    assert pos[40] is None and reasons[40] == "WITHHELD_ALTERNATIVE_OPTIMAL_ALIGNMENT"


def test_unrelated_or_divergent_tail_is_not_forced():
    tail = "VTDLRMIGACHASRFLHMKVECPTELFPPLFLEVFEDQEV"          # TRalpha1-specific C-terminus
    chain = REF[:120] + tail
    target = REF[:120] + "EREVQSSILYKGAAAEGRPGGSLGVHPEGQQLLGMHVVQG"  # TRalpha2-like divergent tail
    pos, _, _, _ = run_calls(chain, target)
    assert all(p is None or chain[i] == target[p] and p == i for i, p in enumerate(pos))
    mapped_tail = [p for p in pos[120:] if p is not None]
    assert len(mapped_tail) <= 1


def test_depositor_tag_residues_withheld_and_isoform_rule():
    assert sn.is_non_reference_residue("expression tag")
    assert sn.is_non_reference_residue("Linker")
    assert not sn.is_non_reference_residue("engineered mutation")
    assert sn.isoform_missing_match(232, 235, [(233, 236)])
    assert not sn.isoform_missing_match(192, 192, [(174, 210)])  # overlap alone is not an isoform match


# --- real data -----------------------------------------------------------------


def test_count_reconciliation(cw, inst):
    frozen = pd.read_parquet(FROZEN_RESIDUES, columns=["instance_id", "is_construct_mapped", "ca_record_present"])
    assert frozen["instance_id"].nunique() == len(inst) == 3159
    ynw = frozen[frozen["instance_id"] == "1YNW.D"]
    assert (ynw["is_construct_mapped"] & ynw["ca_record_present"]).sum() == 0  # why the old audit saw 3158
    assert inst.set_index("instance_id").loc["1YNW.D", "numbering_class"] == "CLEAN"
    s = json.loads(SUMMARY.read_text())
    assert s["after_mapped"] + s["after_withheld"] == s["observed_frozen_rows"]
    assert (s["after_identical"] + s["after_substitution_documented"] + s["after_substitution_undocumented"]
            + s["after_nonstandard_residue_unresolved"] + s["after_modified_residue_parent_identical"]) == s["after_mapped"]
    assert s["before_identical"] + s["before_nonidentical"] == s["observed_frozen_rows"]


def test_vdr_delta165_215_correction(cw):
    x = cw[cw["instance_id"] == "1DB1.A"].set_index("label_seq_id")
    assert (x.loc[1:47, "corrected_uniprot_position"].astype(int) == range(118, 165)).all()
    assert pd.isna(x.loc[48, "corrected_uniprot_position"])  # junction residue: matches UniProt 165 and 216
    assert (x.loc[49:259, "corrected_uniprot_position"].astype(int) == range(217, 428)).all()
    assert (x.loc[49:259, "old_canonical_uniprot_position"].astype(int) == range(166, 377)).all()
    v = pd.read_csv(VDR, keep_default_na=False).set_index("pdb_id")
    assert len(v) == 48
    assert v["status"].value_counts().to_dict() == {"CORRECTED_HIGH_CONFIDENCE": 45, "OLD_MAPPING_CORRECT": 3}
    assert set(v[v["status"] == "OLD_MAPPING_CORRECT"].index) == {"3P8X", "8IQN", "8IQT"}
    for p in ("5V39", "4G2I", "3W0Y", "2HAR"):
        assert v.loc[p, "offsets_new_minus_old"].startswith("+51x")
        assert v.loc[p, "uniprot_165_215_in_chain"] == "DELETED"
        assert float(v.loc[p, "corrected_mapping_after_breakpoint_identity"]) >= 0.99


def test_isoform_cases(cw, inst):
    # FINAL policy: TRalpha numbered on the professor accession P10827-2 (TRalpha1)
    tra = cw[cw["uniprot_id"] == "P10827"]
    assert set(tra["uniprot_sequence_accession"]) == {"P10827-2"}
    assert 370 < tra["corrected_uniprot_position"].dropna().max() <= 410
    tail = tra[tra["corrected_uniprot_position"] > 370]
    assert (tail["correction_status"] == "NEWLY_MAPPED").all() and not tail["in_frozen_lbd_table"].any()
    assert (tail["residue_relation"] == "IDENTICAL").mean() > 0.95
    car = inst[inst["uniprot_id"] == "Q14994"]
    assert set(car["numbering_class"]) == {"ISOFORM_DIFFERENCE_OFFSET"}
    assert car["alignment_gaps"].str.contains("depositor-declared 4/4").all()


def test_engineered_mutant_is_not_a_numbering_offset(inst):
    x = inst.set_index("instance_id").loc["7NFB.A"]
    assert x["numbering_class"] == "CLEAN" and int(x["n_corrected"]) == 0
    assert int(x["n_substitutions_documented"]) == 28


def test_no_silent_coordinate_shift(cw, inst):
    corrected = cw[cw["correction_status"] == "CORRECTED"]
    assert (corrected["correction_reason"] != "").all()
    assert corrected["residue_relation"].isin(["IDENTICAL", "SUBSTITUTION"]).all()
    clean = set(inst[inst["numbering_class"] == "CLEAN"]["instance_id"])
    c = cw[cw["instance_id"].isin(clean) & cw["old_canonical_uniprot_position"].notna()]
    assert (c["old_canonical_uniprot_position"] == c["corrected_uniprot_position"]).all()
    assert set(cw["correction_status"]) == {"UNCHANGED", "CORRECTED", "WITHHELD_PREVIOUSLY_MAPPED", "NEWLY_MAPPED"}


def test_ambiguous_rows_stay_unmapped(cw):
    w = cw[cw["correction_status"] == "WITHHELD_PREVIOUSLY_MAPPED"]
    assert w["corrected_uniprot_position"].isna().all() and (w["confidence"] == "WITHHELD").all()
    assert w["call_reason"].str.startswith(("WITHHELD_", "NOT_ALIGNED")).all()
    newly = cw[cw["correction_status"] == "NEWLY_MAPPED"]
    assert not newly["in_frozen_lbd_table"].any()  # nothing newly forced into the frozen LBD frame
    assert (newly["confidence"] == "HIGH").all()


def test_unaffected_mappings_value_identical_to_frozen(cw):
    frozen = pd.read_parquet(FROZEN_RESIDUES, columns=["instance_id", "entity_label_seq_id", "canonical_uniprot_position",
                                                       "is_construct_mapped"])
    frozen = frozen[frozen["is_construct_mapped"]]
    frozen = frozen.assign(label_seq_id=frozen["entity_label_seq_id"].astype("Int64"))
    j = cw[cw["in_frozen_lbd_table"]].merge(frozen, on=["instance_id", "label_seq_id"], how="outer", indicator=True,
                                             validate="one_to_one")
    assert (j["_merge"] == "both").all()  # every frozen construct-mapped row is represented exactly once
    assert (j["old_canonical_uniprot_position"] == j["canonical_uniprot_position"]).all()
    u = j[j["correction_status"] == "UNCHANGED"]
    assert (u["corrected_uniprot_position"] == u["canonical_uniprot_position"]).all()


def test_corrected_regions_reproduce_crosswalk_with_project_builder(cw):
    reg = pd.read_csv(REGIONS)
    assert len(reg) == cw.groupby(["uniprot_id", "polymer_entity_id"]).ngroups
    sample = reg[reg["polymer_entity_id"].isin(["1DB1_1", "8D8I_1", "1XVP_2", "1UPV_1", "2H79_1", "7NFB_1"])]
    for r in sample.itertuples():
        m = build_entity_to_uniprot_map(json.loads(r.aligned_regions_json))
        assert m["status"] == "OK"
        x = cw[(cw["polymer_entity_id"] == r.polymer_entity_id) & (cw["uniprot_id"] == r.uniprot_id)
               & cw["corrected_uniprot_position"].notna()]
        expect = dict(zip(x["label_seq_id"].astype(int), x["corrected_uniprot_position"].astype(int)))
        assert m["mapping"] == expect


def test_entity_instances_agree():
    e = pd.read_csv(ENTITY)
    assert (e["n_label_seq_with_conflicting_corrected_positions"] == 0).all()


@pytest.mark.skipif(not HAVE_RAW, reason="local mmCIF / UniProt caches not present")
def test_corrected_mapping_deterministic(tmp_path, cw, inst):
    spec = importlib.util.spec_from_file_location("b96", PROJECT_ROOT / "scripts" / "96_build_residue_numbering_crosswalk.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    pdbs = ["1DB1", "5V39", "8D8I", "1XVP", "1UPV", "7NFB", "2H79", "1YNW", "6IAR", "5Y44"]
    mod.build(tmp_path, pdbs, workers=1)
    got = pd.read_parquet(tmp_path / mod.CROSSWALK_OUT)
    exp = cw[cw["pdb_id"].isin(pdbs)].reset_index(drop=True)
    pd.testing.assert_frame_equal(got, exp)
    gi = pd.read_csv(tmp_path / mod.INSTANCE_OUT, keep_default_na=False)
    pd.testing.assert_frame_equal(gi, inst[inst["pdb_id"].isin(pdbs)].reset_index(drop=True))


# --- protection ----------------------------------------------------------------


@pytest.mark.parametrize("rel", sorted(FROZEN_SHA256))
def test_frozen_and_source_files_unchanged(rel):
    expected = FROZEN_SHA256[rel] or git_head_sha256(rel)
    assert sha256(PROJECT_ROOT / rel) == expected


def test_no_stage_files_in_legacy_scanned_directories():
    """The historical identity audit (scripts/93) counts every table in these
    directories; the remap stage must not add any (keeps the frozen report and
    its determinism test valid)."""
    out = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all", *LEGACY_SCANNED_DIRS],
                         cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout
    assert out.strip() == ""


def test_stage_identity_audit_clean():
    a = pd.read_csv(STAGE_IDENTITY)
    assert len(a) == 48 and (a["n_inconsistencies"] == 0).all()
    assert (a["stage_files_checked"] > 0).all()


def test_impact_severity_matches_crosswalk(inst):
    imp = pd.read_csv(IMPACT, keep_default_na=False).set_index("nr_code")
    shifted = set(inst[inst["n_corrected"].astype(int) > 0]["nr_code"])
    assert set(imp[imp["severity"] == "SHIFTED_RESIDUES"].index) == shifted
    assert shifted == {"NR1I1", "NR1D1", "NR1I3", "NR1H2", "NR1I2", "NR3A1", "NR2B1", "NR3A2", "NR1H4"}
