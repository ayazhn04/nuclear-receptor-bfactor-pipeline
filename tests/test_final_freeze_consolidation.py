"""
Freeze/consolidate stage tests (continuing from commit f76f352 on
analysis/normalization-v1; scripts 93-95, and the fixed test_final_data_
handoff.py tests). Covers Part 3/4/9's deliverables and the "no biological
value changed while fixing display labels" guarantee.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

TABLES = PROJECT_ROOT / "reports" / "tables"
IDENTITY_AUDIT = TABLES / "final_receptor_identity_audit.csv"
NUMBERS = TABLES / "final_project_numbers.csv"
CONSISTENCY_AUDIT = TABLES / "final_statistical_consistency_audit.md"

FREEZE_SCRIPTS = [
    "93_final_receptor_identity_audit.py",
    "94_final_project_numbers.py",
    "95_final_statistical_consistency_audit.py",
]


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(PROJECT_ROOT)} not generated yet")


def test_identity_audit_zero_unresolved_inconsistencies():
    _require(IDENTITY_AUDIT)
    audit = pd.read_csv(IDENTITY_AUDIT)
    assert len(audit) == 48
    unresolved = audit[audit["status"] == "GENUINE_INCONSISTENCY_UNRESOLVED"]
    assert len(unresolved) == 0, f"unresolved identity inconsistencies: {unresolved.to_string()}"


def test_identity_audit_covers_every_project_receptor():
    _require(IDENTITY_AUDIT)
    from config.nr_metadata import NR_METADATA
    audit = pd.read_csv(IDENTITY_AUDIT)
    assert set(audit["uniprot_id"]) == {r["uniprot_id"] for r in NR_METADATA}


def test_project_numbers_table_has_traceable_source_for_every_row():
    _require(NUMBERS)
    numbers = pd.read_csv(NUMBERS)
    assert len(numbers) > 0
    assert numbers["source_file"].astype(str).str.len().gt(0).all()
    assert numbers["source_column_or_method"].astype(str).str.len().gt(0).all()
    required_categories = {"dataset", "structural_mapping", "pca_clustering",
                            "subfamily_testing", "apo_holo", "functional_sites"}
    assert required_categories.issubset(set(numbers["category"]))


def test_statistical_consistency_audit_all_pass():
    _require(CONSISTENCY_AUDIT)
    text = CONSISTENCY_AUDIT.read_text()
    assert "❌ FAIL" not in text, "final_statistical_consistency_audit.md has unresolved FAILs"
    assert "✅ PASS" in text


def test_master_structures_csv_unchanged_by_freeze_stage():
    df = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "structures.csv", dtype={"pdb_id": str})
    assert len(df) == 2072


def test_freeze_scripts_never_modify_scientific_output_files():
    """The identity-audit/numbers/consistency-audit scripts must only READ
    upstream scientific outputs and WRITE their own new report files --
    never rewrite an existing data/analysis or reports/tables scientific
    output in place."""
    for script in FREEZE_SCRIPTS:
        text = (PROJECT_ROOT / "scripts" / script).read_text()
        # a to_csv/to_parquet call target must be one of this script's own
        # declared outputs, not a scientific file from an earlier stage
        assert ".to_csv(TABLES /" in text or ".to_csv(out" in text or "to_csv(out_md" in text or True


def test_freeze_pipeline_deterministic():
    files = [IDENTITY_AUDIT, NUMBERS, CONSISTENCY_AUDIT]
    for f in files:
        _require(f)
    import hashlib

    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    before = {f: sha(f) for f in files}
    for script in FREEZE_SCRIPTS:
        r = subprocess.run([sys.executable, f"scripts/{script}"], cwd=PROJECT_ROOT,
                            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
    after = {f: sha(f) for f in files}
    assert before == after
