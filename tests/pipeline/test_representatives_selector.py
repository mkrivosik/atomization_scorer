"""
Tests for the extract_representatives() function.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
import subprocess
from pathlib import Path

import pytest
from Bio import SeqIO

from atomization_scorer import extract_representatives


# --------------------------------------------------------------------------------------
# Helper: create minimal genome FASTA and GEESE
# --------------------------------------------------------------------------------------
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
        "sequence1\t1\t1\t+\t0\t5\n"
        "sequence1\t2\t1\t+\t5\t10\n"
        "sequence2\t3\t1\t+\t0\t5\n"
        "sequence2\t4\t2\t+\t5\t15\n"
    )

    return fasta_file, geese_file


# --------------------------------------------------------------------------------------
# Test: extract representatives with "first" mode
# --------------------------------------------------------------------------------------
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
    assert str(records[0].seq) == "ATGCG"
    assert records[1].id == "sequence2|class_2"
    assert str(records[1].seq) == "CTAGCTAGCT"


# --------------------------------------------------------------------------------------
# Test: extract representatives with "mash" mode
# --------------------------------------------------------------------------------------
def test_extract_representatives_mash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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
        expected_headers = {"sequence1_1", "sequence1_2", "sequence2_3"}
        if set(headers) == expected_headers:
            distances = {
                ("sequence1_1", "sequence1_1"): 0.0,
                ("sequence1_1", "sequence1_2"): 0.3,
                ("sequence1_1", "sequence2_3"): 0.4,
                ("sequence1_2", "sequence1_1"): 0.3,
                ("sequence1_2", "sequence1_2"): 0.0,
                ("sequence1_2", "sequence2_3"): 0.1,
                ("sequence2_3", "sequence1_1"): 0.4,
                ("sequence2_3", "sequence1_2"): 0.1,
                ("sequence2_3", "sequence2_3"): 0.0,
            }
        else:
            distances = {(header, header): 0.0 for header in headers}

        return "\n".join(f"{h1}\t{h2}\t{distance}" for (h1, h2), distance in distances.items())

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
        "sequence1|class_1": "TACGT",
        "sequence2|class_2": "CTAGCTAGCT"
    }

    for record in records:
        assert record.id in expected
        assert str(record.seq) == expected[record.id]


# --------------------------------------------------------------------------------------
# Test: raises FileNotFoundError
# --------------------------------------------------------------------------------------
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


# --------------------------------------------------------------------------------------
# Test: raises ValueError for invalid mode
# --------------------------------------------------------------------------------------
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


# --------------------------------------------------------------------------------------
# Test: raises ValueError if atomization references a missing genome
# --------------------------------------------------------------------------------------
def test_extract_representatives_missing_genome_name_raises_value_error(tmp_path: Path):
    """extract_representatives should raise ValueError if atomization references a missing genome."""
    fasta_file = tmp_path / "genome.fa"
    fasta_file.write_text(">sequence1\nATGCGTACGTAG\n")

    geese_file = tmp_path / "example.geese"
    geese_file.write_text(
        "#name\tatom_nr\tclass\tstrand\tstart\tend\n"
        "missing_sequence\t1\t1\t+\t0\t5\n"
    )

    with pytest.raises(ValueError, match="missing_sequence"):
        extract_representatives(
            genomes_file=fasta_file,
            atomization_file=geese_file,
            output_path=tmp_path / "out.fa",
            mode="first"
        )


# --------------------------------------------------------------------------------------
# Test: cleans up temporary file after Mash failure
# --------------------------------------------------------------------------------------
def test_extract_representatives_mash_failure_cleans_up_temporary_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """extract_representatives should delete the temporary Mash FASTA even when Mash fails."""
    fasta_file, geese_file = create_minimal_files(tmp_path)
    temporary_fasta_path: Path | None = None

    def fake_mash(cmd, **_kwargs):
        nonlocal temporary_fasta_path
        temporary_fasta_path = Path(cmd[-1])
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd)

    monkeypatch.setattr(subprocess, "check_output", fake_mash)

    with pytest.raises(subprocess.CalledProcessError):
        extract_representatives(
            genomes_file=fasta_file,
            atomization_file=geese_file,
            output_path=tmp_path / "out.fa",
            mode="mash"
        )

    if temporary_fasta_path is None:
        pytest.fail("Temporary Mash FASTA path was not captured during the test.")

    assert not temporary_fasta_path.exists()


# --------------------------------------------------------------------------------------
# Test: raises ValueError for malformed Mash output
# --------------------------------------------------------------------------------------
def test_extract_representatives_malformed_mash_output_raises_value_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """extract_representatives should raise ValueError for malformed Mash output."""
    fasta_file, geese_file = create_minimal_files(tmp_path)

    def fake_mash(*_args, **_kwargs):
        return "sequence1_1\tunknown_header\t0.1\n"

    monkeypatch.setattr(subprocess, "check_output", fake_mash)

    with pytest.raises(ValueError, match="Malformed Mash output"):
        extract_representatives(
            genomes_file=fasta_file,
            atomization_file=geese_file,
            output_path=tmp_path / "out.fa",
            mode="mash"
        )
