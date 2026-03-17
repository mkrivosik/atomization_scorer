"""
Tests for align_with_minimap2() function.
"""

import subprocess
from pathlib import Path

import pytest

from atomization_scorer import align_with_minimap2


# --------------------------------------------------------------------------------------
# Helper: create minimal genome FASTA and representatives FASTA
# --------------------------------------------------------------------------------------
def create_minimal_files(tmp_path: Path):
    """Create minimal genome FASTA and representatives FASTA files for testing."""
    genomes = tmp_path / "genomes.fa"
    genomes.write_text(
        ">genome1\nATGCGTACGTAGCTAGCTAG\n"
    )

    representatives = tmp_path / "representatives.fa"
    representatives.write_text(
        ">representative1\nATGCGTACGT\n"
    )

    return genomes, representatives


# --------------------------------------------------------------------------------------
# Test: basic run uses the expected Minimap2 command
# --------------------------------------------------------------------------------------
def test_align_with_minimap2_basic(tmp_path: Path, monkeypatch):
    """align_with_minimap2 should call Minimap2 with the expected default command."""
    genomes, representatives = create_minimal_files(tmp_path)
    output_paf = tmp_path / "alignment.paf"

    calls = []

    # Fake minimap2 run
    def fake_run(*args, **kwargs):
        calls.append(args[0])
        kwargs["stdout"].write("PAF-DATA\n")
        return subprocess.CompletedProcess(args[0], 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    out_path = align_with_minimap2(target=genomes, query=representatives, output_path=output_paf)

    expected_command = [
        "minimap2",
        "-x", "asm20",
        "-c",
        "-p", "0.1",
        str(genomes),
        str(representatives),
    ]

    assert out_path == output_paf
    assert out_path.is_file()
    assert calls == [expected_command]


# --------------------------------------------------------------------------------------
# Test: output directory is created if a nested path is missing
# --------------------------------------------------------------------------------------
def test_align_with_minimap2_creates_output_dir(tmp_path: Path, monkeypatch):
    """align_with_minimap2 should create a missing nested output directory."""
    genomes, representatives = create_minimal_files(tmp_path)
    output_paf = tmp_path / "missing" / "nested" / "alignment.paf"

    assert not output_paf.parent.exists()

    # Fake minimap2 run
    def fake_run(*args, **kwargs):
        kwargs["stdout"].write("PAF-DATA\n")
        return subprocess.CompletedProcess(args[0], 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    out_path = align_with_minimap2(target=genomes, query=representatives, output_path=output_paf)

    assert out_path.is_file()
    assert output_paf.parent.exists()


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if the target or query is missing
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("missing_side", ["target", "query"])
def test_align_with_minimap2_missing_input(tmp_path: Path, monkeypatch, missing_side: str):
    """align_with_minimap2 should raise FileNotFoundError if an input file is missing."""
    target = tmp_path / "target.fa"
    query = tmp_path / "query.fa"

    if missing_side == "target":
        query.write_text(">sequence1\nATGC\n")
    else:
        target.write_text(">sequence1\nATGC\n")

    output_paf = tmp_path / "alignment.paf"
    was_called = False

    # Fake minimap2 run
    def fake_run(*args, **_kwargs):
        nonlocal was_called
        was_called = True
        return subprocess.CompletedProcess(args[0], 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(FileNotFoundError):
        align_with_minimap2(target=target, query=query, output_path=output_paf)

    assert was_called is False


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError if Minimap2 is missing from PATH
# --------------------------------------------------------------------------------------
def test_align_with_minimap2_missing_minimap2_executable(tmp_path: Path, monkeypatch):
    """align_with_minimap2 should raise FileNotFoundError if Minimap2 is not on PATH."""
    genomes, representatives = create_minimal_files(tmp_path)
    output_paf = tmp_path / "alignment.paf"

    # Fake minimap2 run
    def fake_run(*_args, **_kwargs):
        raise FileNotFoundError("minimap2")

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(FileNotFoundError, match="minimap2 executable not found on PATH"):
        align_with_minimap2(target=genomes, query=representatives, output_path=output_paf)


# --------------------------------------------------------------------------------------
# Test: subprocess error if Minimap2 fails
# --------------------------------------------------------------------------------------
def test_align_with_minimap2_failure(tmp_path: Path, monkeypatch):
    """align_with_minimap2 should raise an exception if Minimap2 fails."""
    genomes, representatives = create_minimal_files(tmp_path)
    output_paf = tmp_path / "alignment.paf"

    # Fake minimap2 run
    def fake_run(*args, **_kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=args[0])

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(subprocess.CalledProcessError):
        align_with_minimap2(target=genomes, query=representatives, output_path=output_paf)


# --------------------------------------------------------------------------------------
# Test: missing or empty PAF output is rejected
# --------------------------------------------------------------------------------------
def test_align_with_minimap2_empty_output_raises_value_error(tmp_path: Path, monkeypatch):
    """align_with_minimap2 should raise ValueError if Minimap2 writes no PAF output."""
    genomes, representatives = create_minimal_files(tmp_path)
    output_paf = tmp_path / "alignment.paf"

    # Fake minimap2 run
    def fake_run(*args, **_kwargs):
        return subprocess.CompletedProcess(args[0], 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(ValueError, match="non-empty PAF output file"):
        align_with_minimap2(target=genomes, query=representatives, output_path=output_paf)
