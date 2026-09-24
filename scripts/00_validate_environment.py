"""
Stage 0 environment validation.

Checks that the `biol363_nr` conda environment and the project directory
scaffold are correctly set up before any data-supply work begins.

This script performs NO network access. It only inspects the local Python
environment, writes/reads a throwaway local Parquet file, and checks the
local project directory structure.
"""

from __future__ import annotations

import importlib
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

REQUIRED_PACKAGES = [
    "pandas",
    "numpy",
    "scipy",
    "requests",
    "Bio",  # biopython's importable package name
    "pyarrow",
    "tqdm",
    "yaml",  # pyyaml's importable package name
    "matplotlib",
    "jupyterlab",
    "pytest",
]

# Human-readable label -> importable module name, for the printed report.
PACKAGE_LABELS = {
    "Bio": "biopython",
    "yaml": "pyyaml",
}

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


def check_python_version() -> None:
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    if sys.version_info < (3, 11):
        raise RuntimeError(
            f"Python 3.11+ is required, found {sys.version_info.major}."
            f"{sys.version_info.minor}"
        )


def check_required_packages() -> None:
    print("\nRequired package versions:")
    for module_name in REQUIRED_PACKAGES:
        label = PACKAGE_LABELS.get(module_name, module_name)
        try:
            module = importlib.import_module(module_name)
        except ImportError as exc:
            raise RuntimeError(f"Required package '{label}' failed to import: {exc}") from exc

        # jupyterlab, pytest, biopython, pyyaml all expose __version__
        # in slightly different ways; fall back gracefully.
        version = getattr(module, "__version__", None)
        if version is None and module_name == "jupyterlab":
            import jupyterlab._version  # type: ignore

            version = jupyterlab._version.__version__
        if version is None:
            try:
                from importlib.metadata import version as pkg_version

                version = pkg_version(label)
            except Exception:
                version = "unknown"
        print(f"  {label:12s} {version}")


def check_parquet_roundtrip() -> None:
    print("\nParquet round-trip check:")
    import pandas as pd

    original = pd.DataFrame(
        {
            "receptor_code": ["NR3C1", "ESR1", "RARA"],
            "example_bfactor": [23.4, 31.8, 19.2],
        }
    )

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir) / "stage0_validation.parquet"
        original.to_parquet(tmp_path)
        reloaded = pd.read_parquet(tmp_path)

        if not original.equals(reloaded):
            raise RuntimeError("Parquet round-trip failed: reloaded DataFrame does not match original")

        if tmp_path.exists():
            tmp_path.unlink()

    print("  wrote, read back, and verified a small DataFrame via Parquet — OK")


def check_project_directories() -> None:
    print("\nProject directory check:")
    missing = []
    for rel_dir in REQUIRED_DIRECTORIES:
        full_path = PROJECT_ROOT / rel_dir
        if not full_path.is_dir():
            missing.append(rel_dir)
    if missing:
        raise RuntimeError(f"Missing required project directories: {missing}")
    print(f"  all {len(REQUIRED_DIRECTORIES)} required directories present — OK")


def check_nr_metadata() -> None:
    """Validate config/nr_metadata.py against the CURRENT project state.

    Historical note: at Stage 0 this check required NR_METADATA to be an
    empty placeholder. As of Stage 1.5, scripts/03_freeze_validated_metadata.py
    intentionally populates it with the 48 UniProt-validated receptor
    records, so this check now validates that populated state instead
    (record count, and uniqueness of uniprot_id / nr_code). It remains
    network-free — it only inspects the already-generated module.
    """
    print("\nconfig/nr_metadata.py check:")
    sys.path.insert(0, str(PROJECT_ROOT))
    try:
        from config.nr_metadata import NR_METADATA
    except ImportError as exc:
        raise RuntimeError(f"config/nr_metadata.py failed to import: {exc}") from exc

    if not isinstance(NR_METADATA, list):
        raise RuntimeError("NR_METADATA must be a list")

    expected_record_count = 48
    if len(NR_METADATA) != expected_record_count:
        raise RuntimeError(
            f"NR_METADATA should contain {expected_record_count} Stage-1.5-validated "
            f"records; found {len(NR_METADATA)}. If this is 0, "
            f"scripts/03_freeze_validated_metadata.py may not have been run yet."
        )

    uniprot_ids = [r.get("uniprot_id") for r in NR_METADATA]
    nr_codes = [r.get("nr_code") for r in NR_METADATA]
    if len(set(uniprot_ids)) != expected_record_count:
        raise RuntimeError(
            f"NR_METADATA uniprot_id values are not all unique "
            f"({len(set(uniprot_ids))} unique of {expected_record_count})"
        )
    if len(set(nr_codes)) != expected_record_count:
        raise RuntimeError(
            f"NR_METADATA nr_code values are not all unique "
            f"({len(set(nr_codes))} unique of {expected_record_count})"
        )

    print(
        f"  NR_METADATA imports successfully: {len(NR_METADATA)} records, "
        f"{len(set(uniprot_ids))} unique uniprot_id, {len(set(nr_codes))} unique nr_code — OK"
    )


def main() -> None:
    check_python_version()
    check_required_packages()
    check_parquet_roundtrip()
    check_project_directories()
    check_nr_metadata()
    print("\nENVIRONMENT VALIDATION: PASS")


if __name__ == "__main__":
    main()
