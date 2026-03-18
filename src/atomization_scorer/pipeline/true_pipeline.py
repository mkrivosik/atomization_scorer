"""
true_pipeline.py

Provides the pipeline used to compute true (gold standard) genome atomization.

Functions
---------
compute_true_alignment : Compute true (gold standard) atomization and return a GEESE file path.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
import logging
from pathlib import Path

from atomization_scorer.data_processing import filter_paf, paf_to_geese
from atomization_scorer.visualization import plot_genome_atomization

from .minimap2_aligner import align_with_minimap2
from .representatives_selector import extract_representatives

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# True (Gold Standard) Alignment Pipeline Entry Point
# --------------------------------------------------------------------------------------

def compute_true_alignment(
    genomes_file: Path,
    atomization_file: Path,
    output_directory: Path,
    mode: str = "mash",
    minimum_similarity: float = 0.95,
    minimum_alignment_length: int = 500
) -> Path:
    """
    Run a full true (gold standard) genome atomization pipeline: extract representatives, align genome
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
    Path
        Generated true (gold standard) GEESE file path.
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")

    output_directory.mkdir(parents=True, exist_ok=True)
    logger.info(
        "Computing true alignment with mode=%s for genomes=%s atomization=%s into output=%s",
        mode,
        genomes_file,
        atomization_file,
        output_directory,
    )

    # Extract representatives
    representatives_fasta = output_directory / f"{mode}_representatives.fa"
    logger.info("Extracting representatives to %s", representatives_fasta)
    extract_representatives(
        genomes_file=genomes_file,
        atomization_file=atomization_file,
        output_path=representatives_fasta,
        mode=mode
    )

    # Minimap2 alignment
    paf_file = output_directory / "minimap2_alignments.paf"
    logger.info("Running minimap2 alignment to %s", paf_file)
    align_with_minimap2(
        query=genomes_file,
        target=representatives_fasta,
        output_path=paf_file
    )

    # Filter PAF
    filtered_paf = output_directory / "minimap2_alignment_filtered.paf"
    logger.info(
        "Filtering PAF into %s with minimum_similarity=%s and minimum_alignment_length=%s",
        filtered_paf,
        minimum_similarity,
        minimum_alignment_length,
    )
    filter_paf(
        paf_file=paf_file,
        output_file=filtered_paf,
        minimum_similarity=minimum_similarity,
        minimum_alignment_length=minimum_alignment_length
    )

    # Convert PAF to GEESE
    geese_file = output_directory / "true_atomization.geese"
    logger.info("Converting filtered PAF to GEESE at %s", geese_file)
    paf_to_geese(
        paf_file=filtered_paf,
        output_file=geese_file
    )

    visualization_directory = output_directory / "atomization_visualization"
    logger.info("Generating visualization into %s", visualization_directory)
    plot_genome_atomization(
        genomes_file=genomes_file,
        true_atoms_file=geese_file,
        predicted_atoms_file=atomization_file,
        output_directory=visualization_directory
    )

    logger.info("True alignment pipeline finished with output %s", geese_file)
    return geese_file
