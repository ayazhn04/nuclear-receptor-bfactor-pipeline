"""
Phase 1 technical validation of the corrected v2 foundation
(data/rebuild_v2/, reports/rebuild_v2/): corrected residue numbering + final
analysis span + professor canonical helix map.

These tests check that the professor's map and the numbering correction were
IMPLEMENTED correctly; they do not test the professor's biology.
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

from scripts.utils import canonical_helix as ch  # noqa: E402

V2 = PROJECT_ROOT / "data" / "rebuild_v2"
V2R = PROJECT_ROOT / "reports" / "rebuild_v2"
INTERVALS = PROJECT_ROOT / "data" / "canonical_remap" / "nr_lbd_canonical_helix_intervals_professor_2026.csv"
CROSSWALK = PROJECT_ROOT / "data" / "canonical_remap" / "residue_numbering_crosswalk_v1.parquet"
SEGMAP = V2 / "analysis" / "nr_lbd_structural_segment_map.csv"
DETERMINISM_SANDBOX = PROJECT_ROOT / "data" / "interim" / "rebuild_v2_determinism"
LEGACY_SCANNED_DIRS = ["data/analysis", "data/manifests", "data/processed", "reports/tables"]
AA = {"MSE": "MET"}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def res():
    return pd.read_parquet(V2 / "processed" / "lbd_residue_bfactors_all.parquet",
                           columns=["uniprot_id", "pdb_id", "instance_id", "canonical_uniprot_position",
                                    "canonical_residue_name", "entity_label_seq_id", "is_construct_mapped",
                                    "ca_record_present", "label_comp_ids_json"])


@pytest.fixture(scope="module")
def observed(res):
    o = res[res["is_construct_mapped"] & res["ca_record_present"]].copy()
    o["obs"] = o["label_comp_ids_json"].map(lambda s: json.loads(s)[0] if s not in ("", "[]") else "").replace(AA)
    o["label_seq_id"] = o["entity_label_seq_id"].astype("Int64")
    return o


@pytest.fixture(scope="module")
def seg():
    return pd.read_csv(SEGMAP, keep_default_na=False)


@pytest.fixture(scope="module")
def iv():
    return pd.read_csv(INTERVALS, keep_default_na=False)


# 1-2. professor intervals transferred exactly, never truncated -----------------

def test_every_professor_interval_transferred_exactly(iv, seg):
    for (uid, el), g in iv.groupby(["project_uniprot_id", "canonical_helix"]):
        exp = set().union(*[set(range(a, b + 1)) for a, b in zip(g["start"], g["end"])])
        got = set(seg[(seg["uniprot_id"] == uid) & (seg["structural_segment"] == el)]["canonical_uniprot_position"])
        assert got == exp, (uid, el)


def test_no_professor_helix_truncated_by_old_window(iv, res):
    frame = set(map(tuple, res[["uniprot_id", "canonical_uniprot_position"]].drop_duplicates().to_numpy()))
    with_data = set(res["uniprot_id"])
    for t in iv[iv["project_uniprot_id"].isin(with_data)].itertuples():
        assert all((t.project_uniprot_id, p) in frame for p in range(t.start, t.end + 1)), (t.receptor, t.canonical_helix)
    p = pd.read_parquet(V2 / "analysis" / "receptor_segment_flexibility_profiles.parquet").set_index(["uniprot_id", "structural_segment"])
    for key in [("Q07869", "H1"), ("Q07869", "H2"), ("Q03181", "H1"), ("Q03181", "H2"), ("P04150", "H12"),
                ("P08235", "H12"), ("P06401", "H12"), ("P10275", "H12")]:
        assert p.loc[key, "n_pdb_entries_with_segment"] > 0 and p.loc[key, "median_segment_coverage"] == 1.0, key


# 3-7. numbering ---------------------------------------------------------------

def test_every_observed_residue_at_its_sequence_evidence_position(observed):
    cw = pd.read_parquet(CROSSWALK, columns=["instance_id", "label_seq_id", "corrected_uniprot_position"])
    j = observed.merge(cw, on=["instance_id", "label_seq_id"], how="left")
    assert len(j) == len(observed) and (j["corrected_uniprot_position"] == j["canonical_uniprot_position"]).all()


def test_vdr_plus51_shift_gone(observed):
    v = observed[observed["uniprot_id"] == "P11473"]
    assert (v["obs"] == v["canonical_residue_name"]).mean() > 0.999
    x = observed[observed["instance_id"] == "1DB1.A"].set_index("label_seq_id")
    assert int(x.loc[49, "canonical_uniprot_position"]) == 217 and int(x.loc[47, "canonical_uniprot_position"]) == 164
    assert 48 not in x.index  # ambiguous junction residue withheld, not guessed


def test_deletion_and_isoform_constructs(observed):
    # (instance, last label before the gap, offset before, first label after, offset after) -- every observed residue
    for inst, before_end, off_before, after_start, off_after in [
            ("8D8I.A", 46, 276, 48, 377),     # Rev-erbalpha: UniProt 324-424 deleted (281-322 | 425-614)
            ("5EC9.A", 16, 228, 17, 246),     # RXRalpha: 245-262 deleted
            ("1XVP.B", 129, 102, 131, 106)]:  # CAR1 isoform: 232-235 missing
        x = observed[observed["instance_id"] == inst]
        d = x["canonical_uniprot_position"] - x["label_seq_id"].astype(int)
        assert (d[x["label_seq_id"] <= before_end] == off_before).all(), inst
        assert (d[x["label_seq_id"] >= after_start] == off_after).all(), inst
        assert (x["label_seq_id"] >= after_start).sum() > 50, inst


def test_tralpha_numbered_on_p10827_2(observed):
    seq = ch.read_fasta_sequence(PROJECT_ROOT / "data" / "sources" / "professor_canonical_helix_table"
                                 / "uniprot_isoform_sequences" / "P10827-2.fasta")
    t = observed[observed["uniprot_id"] == "P10827"]
    assert (t["canonical_uniprot_position"] > 370).sum() > 0 and t["canonical_uniprot_position"].max() <= 410
    aa1 = {v: k for k, v in {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
                             "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
                             "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}.items()}
    tail = t[t["canonical_uniprot_position"] > 370]
    assert all(aa1[n] == seq[p - 1] for n, p in zip(tail["canonical_residue_name"], tail["canonical_uniprot_position"]))


def test_real_mutations_remain_mutations(observed):
    n = observed[observed["instance_id"] == "7NFB.A"]
    assert (n["obs"] != n["canonical_residue_name"]).sum() == 28


def test_numbering_disagreements_zero(observed):
    cw = pd.read_parquet(CROSSWALK, columns=["instance_id", "label_seq_id", "residue_relation"])
    mis = observed[observed["obs"] != observed["canonical_residue_name"]].merge(cw, on=["instance_id", "label_seq_id"], how="left")
    # every non-identical observed residue is a real substitution or a modified residue whose parent is identical
    assert set(mis["residue_relation"]) <= {"SUBSTITUTION", "IDENTICAL"}


def test_withheld_residues_absent(res):
    cw = pd.read_parquet(CROSSWALK, columns=["instance_id", "label_seq_id", "correction_status"])
    w = cw[cw["correction_status"] == "WITHHELD_PREVIOUSLY_MAPPED"]
    m = res[res["is_construct_mapped"]].assign(label_seq_id=lambda d: d["entity_label_seq_id"].astype("Int64"))
    assert not w.merge(m, on=["instance_id", "label_seq_id"]).shape[0]


# 9. missing stays missing -----------------------------------------------------

def test_missing_values_stay_missing():
    n = pd.read_parquet(V2 / "analysis" / "normalized_lbd_bfactors_all.parquet",
                        columns=["b_factor_z", "selected_ca_b_iso_or_equiv", "row_normalization_status", "normalization_status"])
    assert not (n["selected_ca_b_iso_or_equiv"].isna() & np.isfinite(n["b_factor_z"])).any()
    assert not (n["normalization_status"].ne("NORMALIZED") & np.isfinite(n["b_factor_z"])).any()
    p = pd.read_parquet(V2 / "analysis" / "receptor_position_flexibility_profiles.parquet")
    z = p[p["n_pdb_entries_observed_at_position"] == 0]
    assert z["median_z"].isna().all() and z["q25_z"].isna().all()


def test_new_vdr_dbd_entries_have_no_fabricated_z():
    d = pd.read_csv(V2R / "vdr_dbd_new_entries_v2.csv")
    assert set(d["pdb_id"]) == {"1KB2", "1KB4", "1KB6", "1YNW"} and len(d) == 7
    assert (d["mapped_residues_in_final_vdr_span"] == 1).all() and (d["mapped_positions"].astype(str) == "125").all()
    assert (d["observed_ca_residues"] == 0).all() and (d["normalized_rows_finite_z"] == 0).all()
    assert (d["instance_normalization_status"] == "NOT_NORMALIZABLE_FEWER_THAN_2_USABLE_B").all()


def test_normalization_edge_rules_unchanged_and_applied():
    inst = pd.read_csv(V2 / "analysis" / "instance_bfactor_summary.csv")
    assert set(inst["normalization_status"]) <= {"NORMALIZED", "NOT_NORMALIZABLE_FEWER_THAN_2_USABLE_B", "NOT_NORMALIZABLE_ZERO_SD"}
    assert (inst[inst["n_ca_usable"] < 2]["normalization_status"] == "NOT_NORMALIZABLE_FEWER_THAN_2_USABLE_B").all()
    assert (inst[inst["normalization_status"] == "NORMALIZED"]["n_ca_usable"] >= 2).all()


# 10. no v3 leak -----------------------------------------------------------------

def test_no_v3_label_in_final_map(seg):
    labels = set(seg["structural_segment"]) - {""}
    assert labels <= set(ch.ELEMENT_ORDER)
    assert not seg["structural_segment"].str.contains("loop|beta_sheet").any()
    assert set(seg["structural_class"]) == {"HELIX", "INTER_HELIX_DERIVED", "NO_PROFESSOR_MAPPING"}
    assert (seg[seg["structural_segment"] != ""]["assignment_method"] == "PROFESSOR_CANONICAL_TABLE_FINAL").all()
    p = pd.read_parquet(V2 / "analysis" / "receptor_segment_flexibility_profiles.parquet")
    assert set(p["structural_segment"]) <= set(ch.ELEMENT_ORDER)


# 11. determinism ------------------------------------------------------------------

@pytest.mark.skipif(not (DETERMINISM_SANDBOX / "logs" / "33_build_analysis_foundation_reports.log").exists(),
                    reason="independent determinism rebuild not complete")
def test_independent_rebuild_is_byte_identical():
    exported = json.loads((V2R / "foundation_export_manifest.json").read_text())["exported_files"]
    for rel in exported:
        if not rel.endswith((".csv", ".parquet")):
            continue  # JSON/markdown provenance carries generation timestamps
        sub = rel.replace("data/rebuild_v2/", "data/").replace("reports/rebuild_v2/", "reports/")
        assert sha256(PROJECT_ROOT / rel) == sha256(DETERMINISM_SANDBOX / sub), rel


# 12. frozen history intact ------------------------------------------------------

def test_no_stage_files_in_legacy_directories():
    out = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all", *LEGACY_SCANNED_DIRS],
                         cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout
    assert out.strip() == ""


FROZEN_BASELINE = "18b8fc912cc3da3e1b1e0108a63a24b0e56b289a"  # frozen v1 release (main)


def test_tracked_files_unmodified():
    """Every file tracked at the frozen v1 baseline is unchanged (content, mode, existence)
    in the working tree, committed or not. Files added after the baseline (v2 work and
    later teammate work) are free to change."""
    subprocess.run(["git", "cat-file", "-e", f"{FROZEN_BASELINE}^{{commit}}"], cwd=PROJECT_ROOT, check=True)
    out = subprocess.run(["git", "diff", "--name-status", "--diff-filter=MDT", FROZEN_BASELINE, "--"],
                         cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout
    assert out.strip() == "", f"frozen baseline files changed:\n{out}"


# count reconciliation ----------------------------------------------------------------

def test_count_reconciliation():
    c = pd.read_csv(V2R / "count_reconciliation_v2.csv").set_index("quantity")
    assert c.loc["receptors (project list)", "v2_corrected"] == 48
    assert c.loc["discovery candidates (receptor x PDB x polymer-entity master rows)", "v2_corrected"] == 2072
    assert c.loc["unique PDB entries with a downloaded mmCIF", "change"] == 3
    assert c.loc["coordinate-assessable candidates (master rows with coordinate data)", "change"] == 4
    assert c.loc["normalized residue rows (finite z)", "v2_corrected"] == c.loc["B-factor-usable residue rows", "v2_corrected"]
