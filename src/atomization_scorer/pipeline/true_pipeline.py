
"""
true_pipeline.py

Provides the pipeline used to compute true (gold standard) genome atomization alignment.

Functions
---------
compute_true_alignment : Computes true (gold standard) alignment and returns a PSL file path.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path

from atomization_scorer.data_processing import filter_paf, paf_to_geese

from .minimap2_aligner import align_with_minimap2
from .representatives_selector import extract_representatives

# ---------------------------------------------------------------------
# True (Gold Standard) Alignment Pipeline Entry Point
# ---------------------------------------------------------------------

def compute_true_alignment(
    genomes_file: Path,
    atomization_file: Path,
    output_directory: Path,
    mode: str = "mash",
    minimum_similarity: float = 0.95,
    minimum_alignment_length: int = 500
) -> Path:
    """
    Run full true (gold standard) genome atomization pipeline: extract representatives, align genome
    sequences on them with minimap2, filter PAF, and convert to GEESE format.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genomes sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.
    output_directory : Path
        Path to the output directory where results are stored.
    mode : str, optional, default: "mash"
        Representative selection mode ("mash" or "first").
    minimum_similarity : float, optional, default: 0.95
        Minimum similarity for PAF filtering.
    minimum_alignment_length : int, optional, default: 500
        Minimum alignment length for PAF filtering.

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.

    Returns
    -------
    str
        Path to the generated true (gold standard) PSL alignment file.
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")

    output_directory.mkdir(parents=True, exist_ok=True)

    # Extract representatives
    representatives_fasta = output_directory / f"{mode}_representatives.fa"
    extract_representatives(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        output_path=representatives_fasta,
        mode=mode
    )

    # Minimap2 alignment
    paf_file = output_directory / "minimap2_alignments.paf"
    align_with_minimap2(
        query=genomes_file,
        target=representatives_fasta,
        output_path=paf_file
    )

    # Filter PAF
    filtered_paf = output_directory / "minimap2_alignment_filtered.paf"
    filter_paf(
        paf_file=paf_file,
        output_file=filtered_paf,
        minimum_similarity=minimum_similarity,
        minimum_alignment_length=minimum_alignment_length
    )

    # Convert PAF to GEESE
    geese_file = output_directory / "true_atomization.geese"
    paf_to_geese(
        paf_file=filtered_paf,
        output_file=geese_file
    )

    return geese_file
