"""
Tests for the read_fasta() function.
"""

from pathlib import Path

import pytest
from Bio.Seq import Seq

from atomization_scorer import read_fasta


# ---------------------------------------------------------------------
# Test: proper reading of a FASTA file
# ---------------------------------------------------------------------
def test_read_fasta_basic(tmp_path: Path):
    """read_fasta should correctly load sequences from a FASTA file."""
    fasta = tmp_path / "test.fa"
    fasta.write_text(
        ">sequence1\nATG\nCGT\n"
        ">sequence2\nGCT\nAGC\n"
    )

    result = read_fasta(fasta_file=fasta)

    assert isinstance(result, dict)
    assert len(result) == 2
    assert result["sequence1"] == Seq("ATGCGT")
    assert result["sequence2"] == Seq("GCTAGC")


# ---------------------------------------------------------------------
# Test: empty FASTA file
# ---------------------------------------------------------------------
def test_read_fasta_empty_file(tmp_path: Path):
    """read_fasta should return an empty dictionary for an empty FASTA file."""
    fasta = tmp_path / "empty.fa"
    fasta.write_text("")

    result = read_fasta(fasta_file=fasta)

    assert result == {}


# ---------------------------------------------------------------------
# Test: missing FASTA file
# ---------------------------------------------------------------------
def test_read_fasta_missing_file(tmp_path: Path):
    """read_fasta should raise FileNotFoundError if the file does not exist."""
    missing = tmp_path / "missing.fa"

    with pytest.raises(FileNotFoundError):
        read_fasta(fasta_file=missing)
