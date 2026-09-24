"""
Stage 0 infrastructure tests.

These tests verify the *bootstrap* state of the repository: that the
directory scaffold exists and that `config/nr_metadata.py` is a well-formed
NR_METADATA module. They intentionally know nothing about receptor biology,
PDB IDs, or UniProt — that validation belongs to later stages.

Note on `NR_METADATA` emptiness: at Stage 0, this suite asserted
`NR_METADATA == []` (the bootstrap placeholder, before any receptor data
existed). As of Stage 1.5, `config/nr_metadata.py` is intentionally
generated and populated (48 UniProt-validated records) by
`scripts/03_freeze_validated_metadata.py` — see
`tests/test_stage1_5_frozen_metadata.py` for the detailed content checks.
The test below now asserts that transition explicitly, so a regression back
to an empty/placeholder file after Stage 1.5 is caught.
"""

from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent

REQUIRED_DIRECTORIES = [
    "config",
    "data/raw/mmcif",
    "data/raw/api",
    "data/interim/mappings",
    "data/interim/qc",
    "data/interim/ligand_annotations",
    "data/processed",
    "data/manifests",
    "scripts",
    "scripts/utils",
    "notebooks",
    "logs",
    "reports/figures",
    "reports/tables",
    "tests",
]


@pytest.mark.parametrize("rel_dir", REQUIRED_DIRECTORIES)
def test_required_directory_exists(rel_dir: str) -> None:
    assert (PROJECT_ROOT / rel_dir).is_dir(), f"missing required directory: {rel_dir}"


def test_nr_metadata_exists_and_is_a_list() -> None:
    import sys

    sys.path.insert(0, str(PROJECT_ROOT))
    from config.nr_metadata import NR_METADATA

    assert isinstance(NR_METADATA, list)


def test_nr_metadata_reflects_post_stage1_5_state() -> None:
    """Superseded by Stage 1.5 (see module docstring). Historically this
    asserted `NR_METADATA == []`; it now asserts the post-freeze state:
    non-empty, and specifically 48 records (the validated master)."""
    import sys

    sys.path.insert(0, str(PROJECT_ROOT))
    from config.nr_metadata import NR_METADATA

    assert len(NR_METADATA) == 48, (
        "NR_METADATA should contain the 48 Stage-1-validated records after "
        "Stage 1.5 froze the working receptor master — if this is 0, the "
        "freeze script may not have been run; if some other number, "
        "investigate before proceeding"
    )


def test_gitignore_protects_biological_data_dirs() -> None:
    gitignore_text = (PROJECT_ROOT / ".gitignore").read_text()
    for pattern in ("data/raw/mmcif/", "data/raw/api/", "data/interim/", "data/processed/", "logs/"):
        assert pattern in gitignore_text, f".gitignore is missing pattern: {pattern}"


def test_gitignore_does_not_exclude_manifests() -> None:
    ignore_patterns = {
        line.strip()
        for line in (PROJECT_ROOT / ".gitignore").read_text().splitlines()
        if line.strip() and not line.strip().startswith("#")
    }
    assert "data/manifests/" not in ignore_patterns
    assert "data/manifests" not in ignore_patterns
