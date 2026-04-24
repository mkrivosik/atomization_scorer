"""
Tests for the cli.py command-line interface of Atomization Scorer.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
import subprocess
import sys
from pathlib import Path


# --------------------------------------------------------------------------------------
# Test: CLI runs successfully with valid input files
# --------------------------------------------------------------------------------------
def test_cli_valid(mini_fasta: Path, mini_geese: Path, output_dir: Path):
    """CLI should run successfully and produce output with valid input files."""
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "atomization_scorer.cli",
            str(mini_fasta),
            str(mini_geese),
            str(output_dir),
        ],
        capture_output=True,
        text=True
    )

    print(f"Output directory: {output_dir}")
    assert result.returncode == 0
    assert f"Processing files: genomes={mini_fasta}" in result.stderr
    assert f"atomization={mini_geese}" in result.stderr
    assert f"output={output_dir}" in result.stderr
    assert "Overall score result:" in result.stderr
    assert output_dir.exists()


# --------------------------------------------------------------------------------------
# Test: CLI exits if a file is missing
# --------------------------------------------------------------------------------------
def test_cli_missing_file(mini_geese: Path, output_dir: Path):
    """CLI should exit with an error if a required input file is missing."""
    missing_fasta = mini_geese.parent / "missing.fa"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "atomization_scorer.cli",
            str(missing_fasta),
            str(mini_geese),
            str(output_dir),
        ],
        capture_output=True,
        text=True
    )

    assert result.returncode != 0
    assert "file not found" in result.stderr.lower()


# --------------------------------------------------------------------------------------
# Test: CLI exits if file has wrong extension
# --------------------------------------------------------------------------------------
def test_cli_wrong_extension(mini_geese: Path, tmp_path: Path, output_dir: Path):
    """CLI should exit with an error if an input file has the wrong extension."""
    wrong_extension_file = tmp_path / "example.txt"
    wrong_extension_file.write_text(">sequence1\nATGC")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "atomization_scorer.cli",
            str(wrong_extension_file),
            str(mini_geese),
            str(output_dir),
        ],
        capture_output=True,
        text=True
    )

    assert result.returncode != 0
    assert "must have one of the following extensions" in result.stderr.lower()


# --------------------------------------------------------------------------------------
# Test: CLI creates output directory if it does not exist
# --------------------------------------------------------------------------------------
def test_cli_creates_output_dir(mini_fasta: Path, mini_geese: Path, tmp_path: Path):
    """CLI should automatically create the output directory if it is missing."""
    new_output_directory = tmp_path / "new_output"

    assert not new_output_directory.exists()

    subprocess.run(
        [
            sys.executable,
            "-m",
            "atomization_scorer.cli",
            str(mini_fasta),
            str(mini_geese),
            str(new_output_directory),
        ],
        check=False
    )

    assert new_output_directory.exists()
    assert new_output_directory.is_dir()
