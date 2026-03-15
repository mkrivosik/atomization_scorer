"""
Tests for the extract_representatives() function.
"""

import subprocess
from pathlib import Path

import pytest
from Bio import SeqIO

from atomization_scorer import extract_representatives


# ---------------------------------------------------------------------
# Helper: create minimal genome FASTA and GEESE
# ---------------------------------------------------------------------
def create_minimal_files(tmp_path: Path):
    """Create minimal genome FASTA and GEESE files for testing."""
    fasta_file = tmp_path / "genome.fa"
    fasta_file.write_text(
        ">sequence1\nATGCGTACGTAGCTAGCTAG\n"
        ">sequence2\nGCTAGCTAGCTAGCTAGCTA\n"
    )

    geese_file = tmp_path / "example.geese"
    geese_file.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
        "sequence1\t2\t1\t+\t10\t20\n"
        "sequence2\t3\t2\t+\t0\t10\n"
    )

    return fasta_file, geese_file


# ---------------------------------------------------------------------
# Test: extract representatives with "first" mode
# ---------------------------------------------------------------------
def test_extract_representatives_first(tmp_path: Path):
    """extract_representatives should correctly select the first atom as a representative."""
    fasta_file, geese_file = create_minimal_files(tmp_path)
    output_fasta = tmp_path / "representatives_first.fa"

    out_path = extract_representatives(
        genomes_file=fasta_file,
        atomization_file=geese_file,
        output_path=output_fasta,
        mode="first"
    )
    assert out_path.is_file()

    records = list(SeqIO.parse(out_path, "fasta"))
    assert len(records) == 2
    assert records[0].id == "sequence1|class_1"
    assert str(records[0].seq) == "ATGCGTACGT"
    assert records[1].id == "sequence2|class_2"
    assert str(records[1].seq) == "GCTAGCTAGC"


# ---------------------------------------------------------------------
# Test: extract representatives with "mash" mode
# ---------------------------------------------------------------------
def test_extract_representatives_mash(tmp_path: Path, monkeypatch):
    """extract_representatives should correctly select a representative using mash distances."""
    fasta_file, geese_file = create_minimal_files(tmp_path)
    output_fasta = tmp_path / "representatives_mash.fa"

    # Fake mash output
    def fake_mash(cmd, **_kwargs):
        temporary_fasta_path = cmd[-1]
        headers = []
        with open(temporary_fasta_path) as file:
            for line in file:
                if line.startswith(">"):
                    headers.append(line[1:].strip().split()[0])
        output = []
        for h1 in headers:
            for h2 in headers:
                distance = 0.0 if h1 == h2 else 0.1
                output.append(f"{h1}\t{h2}\t{distance}")
        return "\n".join(output)

    monkeypatch.setattr(subprocess, "check_output", fake_mash)

    out_path = extract_representatives(
        genomes_file=fasta_file,
        atomization_file=geese_file,
        output_path=output_fasta,
        mode="mash"
    )
    assert out_path.is_file()

    records = list(SeqIO.parse(out_path, "fasta"))
    assert len(records) == 2

    expected = {
        "sequence1|class_1": "ATGCGTACGT",
        "sequence2|class_2": "GCTAGCTAGC"
    }

    for record in records:
        assert record.id in expected
        assert str(record.seq) == expected[record.id]


# ---------------------------------------------------------------------
# Test: raises FileNotFoundError
# ---------------------------------------------------------------------
def test_extract_representatives_missing_files(tmp_path: Path):
    """extract_representatives should raise FileNotFoundError if a genome or atomization file is missing."""
    missing_fasta = tmp_path / "missing.fa"
    geese_file = tmp_path / "example.geese"
    geese_file.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "sequence1\t1\t1\t+\t0\t10\n"
    )

    with pytest.raises(FileNotFoundError):
        extract_representatives(
            genomes_file=missing_fasta,
            atomization_file=geese_file,
            output_path=tmp_path / "out.fa"
        )

    fasta_file = tmp_path / "genome.fa"
    fasta_file.write_text(">sequence1\nATGCGTACGTAG\n")
    missing_geese = tmp_path / "missing.geese"

    with pytest.raises(FileNotFoundError):
        extract_representatives(
            genomes_file=fasta_file,
            atomization_file=missing_geese,
            output_path=tmp_path / "out.fa"
        )


# ---------------------------------------------------------------------
# Test: raises ValueError for invalid mode
# ---------------------------------------------------------------------
def test_extract_representatives_invalid_mode(tmp_path: Path):
    """extract_representatives should raise ValueError if the mode is not 'first' or 'mash'."""
    fasta_file, geese_file = create_minimal_files(tmp_path)

    with pytest.raises(ValueError):
        extract_representatives(
            genomes_file=fasta_file,
            atomization_file=geese_file,
            output_path=tmp_path / "out.fa",
            mode="invalid"
        )


def test_extract_representatives_uses_half_open_coordinates(tmp_path: Path):
    """extract_representatives should slice genomes using half-open intervals."""
    fasta_file, geese_file = create_minimal_files(tmp_path)
    output_fasta = tmp_path / "representatives_half_open.fa"

    out_path = extract_representatives(
        genomes_file=fasta_file,
        atomization_file=geese_file,
        output_path=output_fasta,
        mode="first"
    )

    records = {record.id: str(record.seq) for record in SeqIO.parse(out_path, "fasta")}
    assert records["sequence1|class_1"] == "ATGCGTACGT"
