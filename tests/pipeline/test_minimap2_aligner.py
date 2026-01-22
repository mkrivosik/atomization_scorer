"""
Tests for align_with_minimap2() function.
"""

import subprocess
from pathlib import Path

import pytest

from atomization_scorer import align_with_minimap2


# ------------------------------------------------------------------------------
# Helper: create minimal genome FASTA and representatives FASTA
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
# Test: basic run, PAF path returned
# ------------------------------------------------------------------------------
def test_align_with_minimap2_basic(tmp_path: Path, monkeypatch):
    """align_with_minimap2 should call minimap2 and return PAF path."""
    genomes, representatives = create_minimal_files(tmp_path)
    output_paf = tmp_path / "alignment.paf"

    calls = []

    # Fake minimap2 run
    def fake_run(*args, **_kwargs):
        calls.append(args[0])
        output_paf.write_text("PAF-DATA")
        return subprocess.CompletedProcess(args[0], 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    out_path = align_with_minimap2(target=genomes, query=representatives, output_path=output_paf)

    assert out_path.is_file()
    assert calls
    assert "minimap2" in calls[0][0]


# ------------------------------------------------------------------------------
# Test: subprocess error if minimap2 fails
# ------------------------------------------------------------------------------
def test_align_with_minimap2_failure(tmp_path: Path, monkeypatch):
    """align_with_minimap2 should raise an exception if minimap2 fails."""
    genomes, representatives = create_minimal_files(tmp_path)
    output_paf = tmp_path / "alignment.paf"

    # Fake minimap2 run
    def fake_run(*args, **_kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=args[0])

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(subprocess.CalledProcessError):
        align_with_minimap2(target=genomes, query=representatives, output_path=output_paf)


# ------------------------------------------------------------------------------
# Test: output directory is created if missing
# ------------------------------------------------------------------------------
def test_align_with_minimap2_creates_output_dir(tmp_path: Path, output_dir: Path, monkeypatch):
    """align_with_minimap2 should create the output directory if it doesn't exist."""
    genomes, representatives = create_minimal_files(tmp_path)
    output_paf = output_dir / "alignment.paf"

    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 0))

    out_path = align_with_minimap2(target=genomes, query=representatives, output_path=output_paf)
    assert out_path.is_file()
    assert output_paf.parent.exists()


# ------------------------------------------------------------------------------
# Test: raises FileNotFoundError if target is missing
# ------------------------------------------------------------------------------
def test_align_with_minimap2_missing_target(tmp_path: Path):
    """align_with_minimap2 should raise FileNotFoundError if target file is missing."""
    missing_target = tmp_path / "missing_target.fa"
    query = tmp_path / "query.fa"
    query.write_text(">sequence1\nATGC\n")
    output_paf = tmp_path / "alignment.paf"

    with pytest.raises(FileNotFoundError):
        align_with_minimap2(target=missing_target, query=query, output_path=output_paf)


# ------------------------------------------------------------------------------
# Test: raises FileNotFoundError if query is missing
# ------------------------------------------------------------------------------
def test_align_with_minimap2_missing_query(tmp_path: Path):
    """align_with_minimap2 should raise FileNotFoundError if query file is missing."""
    target = tmp_path / "target.fa"
    target.write_text(">sequence1\nATGC\n")
    missing_query = tmp_path / "missing_query.fa"
    output_paf = tmp_path / "alignment.paf"

    with pytest.raises(FileNotFoundError):
        align_with_minimap2(target=target, query=missing_query, output_path=output_paf)
