"""
Tests for the read_fasta() function.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

import pytest
from Bio.Seq import Seq

from atomization_scorer import read_fasta


# --------------------------------------------------------------------------------------
# Test: basic FASTA reading
# --------------------------------------------------------------------------------------
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


# --------------------------------------------------------------------------------------
# Test: empty FASTA file
# --------------------------------------------------------------------------------------
def test_read_fasta_empty_file(tmp_path: Path):
    """read_fasta should return an empty dictionary for an empty FASTA file."""
    fasta = tmp_path / "empty.fa"
    fasta.write_text("")

    result = read_fasta(fasta_file=fasta)

    assert result == {}


# --------------------------------------------------------------------------------------
# Test: FASTA headers are normalized to the Biopython record ID
# --------------------------------------------------------------------------------------
def test_read_fasta_header_normalization(tmp_path: Path):
    """read_fasta should use the normalized FASTA record ID as the dictionary key."""
    fasta = tmp_path / "headers.fa"
    fasta.write_text(
        ">sequence1 full description here\nATGC\n"
        ">sequence2 another description\nGGTT\n"
    )

    result = read_fasta(fasta_file=fasta)

    assert "sequence1" in result
    assert "sequence1 full description here" not in result
    assert result["sequence1"] == Seq("ATGC")
    assert result["sequence2"] == Seq("GGTT")


# --------------------------------------------------------------------------------------
# Test: malformed FASTA input is rejected
# --------------------------------------------------------------------------------------
def test_read_fasta_malformed_file(tmp_path: Path):
    """read_fasta should raise ValueError if the FASTA content is malformed."""
    fasta = tmp_path / "malformed.fa"
    fasta.write_text(">\nATGC\n")

    with pytest.raises(ValueError, match="Malformed FASTA file"):
        read_fasta(fasta_file=fasta)


# --------------------------------------------------------------------------------------
# Test: sequence content is preserved exactly
# --------------------------------------------------------------------------------------
def test_read_fasta_preserves_sequence_content(tmp_path: Path):
    """read_fasta should preserve lowercase bases and ambiguous symbols exactly."""
    fasta = tmp_path / "preserve.fa"
    fasta.write_text(">sequence1\naTgCNn-\n")

    result = read_fasta(fasta_file=fasta)

    assert result["sequence1"] == Seq("aTgCNn-")


# --------------------------------------------------------------------------------------
# Test: missing FASTA file
# --------------------------------------------------------------------------------------
def test_read_fasta_missing_file(tmp_path: Path):
    """read_fasta should raise FileNotFoundError if the FASTA file does not exist."""
    missing = tmp_path / "missing.fa"

    with pytest.raises(FileNotFoundError):
        read_fasta(fasta_file=missing)


# --------------------------------------------------------------------------------------
# Test: duplicate FASTA IDs are rejected
# --------------------------------------------------------------------------------------
def test_read_fasta_duplicate_ids(tmp_path: Path):
    """read_fasta should raise ValueError if the FASTA file contains duplicate IDs."""
    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(
        ">sequence1\nATGC\n"
        ">sequence1\nGGTT\n"
    )

    with pytest.raises(ValueError, match="Duplicate FASTA ID found: sequence1"):
        read_fasta(fasta_file=fasta)
