"""
Fixtures Atomization Scorer tests.
"""

import pytest
from pathlib import Path


# ---------------------------------------------------------------------
# Fixture: directory containing test input files
# ---------------------------------------------------------------------
@pytest.fixture
def fixtures_directory():
    """
    Returns the path to the fixtures directory for test input files.
    """
    return Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------
# Fixture: temporary output directory for test results
# ---------------------------------------------------------------------
@pytest.fixture
def output_dir(tmp_path: Path):
    """
    Returns a temporary output directory for CLI results.
    """
    path = tmp_path / "output"
    path.mkdir(exist_ok=True)
    return path


# ---------------------------------------------------------------------
# Optional helper fixtures for direct access to mini files
# ---------------------------------------------------------------------
@pytest.fixture
def mini_fasta(fixtures_directory: Path) -> Path:
    """Path to the mini FASTA file in fixtures."""
    return fixtures_directory / "mini.fa"


@pytest.fixture
def mini_geese(fixtures_directory: Path) -> Path:
    """Path to the mini GEESE file in fixtures."""
    return fixtures_directory / "mini.geese"
