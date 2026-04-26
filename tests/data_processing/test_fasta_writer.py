"""
Tests for the fasta_writer.py module.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
import pytest

from atomization_scorer import read_fasta, write_fasta


# --------------------------------------------------------------------------------------
# Test: basic FASTA writing
# --------------------------------------------------------------------------------------
def test_write_fasta_basic(tmp_path: Path):
    """write_fasta should write sequences with correct FASTA formatting."""
    output = tmp_path / "out.fa"
    sequences = {
        "sequence1": Seq("ATGCGT"),
        "sequence2": Seq("GCTAGC")
    }

    path = write_fasta(sequences=sequences, output_path=output)

    assert path == output
    assert output.exists()

    records = list(SeqIO.parse(output, "fasta"))
    assert len(records) == 2
    assert records[0].id == "sequence1"
    assert str(records[0].seq) == "ATGCGT"
    assert records[1].id == "sequence2"
    assert str(records[1].seq) == "GCTAGC"


# --------------------------------------------------------------------------------------
# Test: parent directory creation
# --------------------------------------------------------------------------------------
def test_write_fasta_creates_parent_dir(output_dir):
    """write_fasta should create an output directory if it does not exist."""
    nested_directory = output_dir / "nested" / "folder"
    output = nested_directory / "sequences.fa"
    sequences = {"sequence": Seq("AAA")}

    path = write_fasta(sequences=sequences, output_path=output)

    assert output.exists()
    assert path == output


# --------------------------------------------------------------------------------------
# Test: writing an empty FASTA
# --------------------------------------------------------------------------------------
def test_write_fasta_empty_dict(tmp_path: Path):
    """write_fasta should produce an empty file when sequences are empty."""
    output = tmp_path / "empty.fa"
    path = write_fasta(sequences={}, output_path=output)

    assert output.exists()
    assert path == output
    assert output.read_text() == ""


# --------------------------------------------------------------------------------------
# Test: sequence content is preserved for unusual characters
# --------------------------------------------------------------------------------------
def test_write_fasta_preserves_unusual_characters(tmp_path: Path):
    """write_fasta should preserve lowercase bases and unusual sequence symbols exactly."""
    output = tmp_path / "unusual.fa"
    sequences = {"sequence1": Seq("aTgCNn-")}

    write_fasta(sequences=sequences, output_path=output)

    result = read_fasta(fasta_file=output)

    assert result["sequence1"] == Seq("aTgCNn-")


# --------------------------------------------------------------------------------------
# Test: write_fasta output can be read back with read_fasta
# --------------------------------------------------------------------------------------
def test_write_fasta_round_trip_with_read_fasta(tmp_path: Path):
    """write_fasta should produce FASTA that round-trips correctly through read_fasta."""
    output = tmp_path / "round_trip.fa"
    sequences = {
        "sequence1": Seq("ATGCGT"),
        "sequence2": Seq("GCTAGC"),
    }

    write_fasta(sequences=sequences, output_path=output)
    result = read_fasta(fasta_file=output)

    assert result == sequences
def test_write_fasta_wraps_sequences_at_80_characters(tmp_path: Path):
    """write_fasta should wrap FASTA sequence lines to 80 characters."""
    output = tmp_path / "wrapped.fa"
    sequence = Seq("A" * 85)

    write_fasta(sequences={"sequence1": sequence}, output_path=output)

    lines = output.read_text().splitlines()
    assert lines[0] == ">sequence1"
    assert lines[1] == "A" * 80
    assert lines[2] == "A" * 5

    records = list(SeqIO.parse(output, "fasta"))
    assert len(records) == 1
    assert records[0].id == "sequence1"
    assert records[0].seq == sequence


# --------------------------------------------------------------------------------------
# Test: FASTA headers may include spaces
# --------------------------------------------------------------------------------------
def test_write_fasta_allows_spaces_in_header(tmp_path: Path):
    """write_fasta should preserve headers containing spaces."""
    output = tmp_path / "header_with_spaces.fa"
    sequences = {"sequence 1": Seq("ATGC")}

    write_fasta(sequences=sequences, output_path=output)

    records = list(SeqIO.parse(output, "fasta"))
    assert len(records) == 1
    assert records[0].description == "sequence 1"
    assert str(records[0].seq) == "ATGC"


# --------------------------------------------------------------------------------------
# Test: empty FASTA sequences are rejected
# --------------------------------------------------------------------------------------
def test_write_fasta_empty_sequence(tmp_path: Path):
    """write_fasta should raise ValueError if any FASTA sequence is empty."""
    output = tmp_path / "empty_sequence.fa"
    sequences = {"sequence1": Seq("")}

    with pytest.raises(ValueError, match="FASTA sequence must not be empty"):
        write_fasta(sequences=sequences, output_path=output)


# --------------------------------------------------------------------------------------
# Test: unsafe FASTA headers are rejected
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("header", "error_message"),
    [
        ("", "FASTA header must not be empty"),
        ("sequence>1", "FASTA header contains '>'"),
        ("sequence\n1", "FASTA header contains newline characters"),
        ("sequence\r1", "FASTA header contains newline characters"),
    ],
)
def test_write_fasta_header_safety(tmp_path: Path, header: str, error_message: str):
    """write_fasta should raise ValueError if a FASTA header contains unsafe characters."""
    output = tmp_path / "unsafe_header.fa"
    sequences = {header: Seq("ATGC")}

    with pytest.raises(ValueError, match=error_message):
        write_fasta(sequences=sequences, output_path=output)


# --------------------------------------------------------------------------------------
# Test: sequence output is wrapped to 80 characters per line
# --------------------------------------------------------------------------------------
