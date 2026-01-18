"""
Tests for write_fasta() function.
"""

from pathlib import Path
from Bio import SeqIO
from atomization_scorer import write_fasta


# ---------------------------------------------------------------------
# Test: basic FASTA writing
# ---------------------------------------------------------------------
def test_write_fasta_basic(tmp_path: Path):
    """write_fasta should write sequences with correct FASTA formatting."""
    output = tmp_path / "out.fa"
    sequences = {
        "sequence1": "ATGCGT",
        "sequence2": "GCTAGC"
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


# ---------------------------------------------------------------------
# Test: parent directory creation
# ---------------------------------------------------------------------
def test_write_fasta_creates_parent_dir(output_dir):
    """write_fasta should create output directory if it does not exist."""
    nested_directory = output_dir / "nested" / "folder"
    output = nested_directory / "sequences.fa"
    sequences = {"sequence": "AAA"}

    path = write_fasta(sequences=sequences, output_path=output)

    assert output.exists()
    assert path == output


# ---------------------------------------------------------------------
# Test: writing an empty FASTA
# ---------------------------------------------------------------------
def test_write_fasta_empty_dict(tmp_path: Path):
    """write_fasta should produce an empty file when sequences are empty."""
    output = tmp_path / "empty.fa"
    path = write_fasta(sequences={}, output_path=output)

    assert output.exists()
    assert path == output
    assert output.read_text() == ""
