"""
Tests for the cli.py command-line interface of Atomization Scorer.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
import json
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
# Test: --quiet suppresses INFO-level log messages
# --------------------------------------------------------------------------------------
def test_cli_quiet_suppresses_info(mini_fasta: Path, mini_geese: Path, output_dir: Path):
    """CLI with --quiet should suppress INFO-level log messages."""
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "atomization_scorer.cli",
            str(mini_fasta),
            str(mini_geese),
            str(output_dir),
            "--quiet",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert "Processing files:" not in result.stderr
    assert "Overall score result:" not in result.stderr


# --------------------------------------------------------------------------------------
# Test: run_parameters.json is written with correct structure
# --------------------------------------------------------------------------------------
def test_cli_run_parameters_json(mini_fasta: Path, mini_geese: Path, output_dir: Path):
    """CLI should write run_parameters.json with the correct structure and input paths."""
    subprocess.run(
        [
            sys.executable,
            "-m",
            "atomization_scorer.cli",
            str(mini_fasta),
            str(mini_geese),
            str(output_dir),
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    run_parameters_file = output_dir / "run_parameters.json"
    assert run_parameters_file.exists()

    with run_parameters_file.open() as file:
        data = json.load(file)

    assert set(data.keys()) == {"timestamp", "inputs", "scoring", "pipeline", "minimap2", "diagnostics"}
    assert data["inputs"]["genomes_file"] == str(mini_fasta.resolve())
    assert data["inputs"]["atomization_file"] == str(mini_geese.resolve())


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
