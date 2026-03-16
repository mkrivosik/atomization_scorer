"""
Tests for plot_genome_atomization() function.
"""

from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from PIL import Image

from atomization_scorer import plot_genome_atomization

# --------------------------------------------------------------------------------------
# Helper: create minimal genomes and atomization files
# --------------------------------------------------------------------------------------
def create_minimal_genomes_and_atoms(tmp_path: Path):
    """Create minimal genomes FASTA and predicted/true atomization TSV files."""
    genomes_file = tmp_path / "genomes.fasta"
    records = [
        SeqRecord(Seq("ATGC" * 50), id="genome1"),
        SeqRecord(Seq("ATGC" * 60), id="genome2")
    ]
    SeqIO.write(records, genomes_file, "fasta")

    true_atoms_file = tmp_path / "true_atoms.geese"
    df_true = pd.DataFrame({
        "name": ["genome1", "genome1", "genome2"],
        "class": ["A", "B", "A"],
        "start": [0, 100, 170],
        "end": [50, 150, 190]
    })
    df_true.to_csv(true_atoms_file, sep="\t", index=False)

    # Predicted atoms
    predicted_atoms_file = tmp_path / "predicted_atoms.geese"
    df_pred = pd.DataFrame({
        "name": ["genome1", "genome2", "genome2"],
        "class": ["X", "Y", "Z"],
        "start": [10, 5, 160],
        "end": [60, 25, 195]
    })
    df_pred.to_csv(predicted_atoms_file, sep="\t", index=False)

    return genomes_file, true_atoms_file, predicted_atoms_file

# --------------------------------------------------------------------------------------
# Test: basic visualization generation
# --------------------------------------------------------------------------------------
def test_plot_genome_atomization(tmp_path: Path, output_dir: Path):
    """plot_genome_atomization should create one PNG file per genome in the output directory."""
    genomes_file, true_atoms_file, predicted_atoms_file = create_minimal_genomes_and_atoms(tmp_path)

    plot_genome_atomization(
        genomes_file=genomes_file,
        true_atoms_file=true_atoms_file,
        predicted_atoms_file=predicted_atoms_file,
        output_directory=output_dir
    )

    for genome_name in ["genome1", "genome2"]:
        png_file = output_dir / f"{genome_name}.png"
        assert png_file.is_file(), f"Visualization PNG missing for {genome_name}"

        with Image.open(png_file) as img:
            img.verify()