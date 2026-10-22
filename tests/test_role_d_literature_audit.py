"""Integrity checks for the all-receptor Role D literature-evidence audit."""
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "reports" / "role_d" / "literature_review"
SITES = {
    "AF2_H12",
    "LIGAND_BINDING_POCKET",
    "COACTIVATOR_GROOVE",
    "DIMERIZATION_INTERFACE",
}


def read_csv(name):
    with (AUDIT / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_all_48_receptors_and_four_site_classes_are_audited_once():
    rows = read_csv("per_receptor_site_evidence_matrix.csv")
    keys = {(row["uniprot_id"], row["site_class"]) for row in rows}

    assert len(rows) == 48 * 4
    assert len(keys) == 48 * 4
    assert len({row["uniprot_id"] for row in rows}) == 48
    assert Counter(row["site_class"] for row in rows) == {site: 48 for site in SITES}


def test_matrix_citations_resolve_and_search_log_has_same_audit_universe():
    matrix = read_csv("per_receptor_site_evidence_matrix.csv")
    bibliography = read_csv("bibliography.csv")
    searches = read_csv("pubmed_systematic_search_log.csv")

    bibliography_keys = {row["citation_key"] for row in bibliography}
    for row in matrix:
        for key in filter(None, row["citation_keys"].split(";")):
            assert key in bibliography_keys

    matrix_keys = {(row["uniprot_id"], row["site_class"]) for row in matrix}
    search_keys = {(row["uniprot_id"], row["site_class"]) for row in searches}
    assert len(searches) == 48 * 4
    assert len(search_keys) == 48 * 4
    assert search_keys == matrix_keys
    assert {row["request_status"] for row in searches} == {"OK"}
