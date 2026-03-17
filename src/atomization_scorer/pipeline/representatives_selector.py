"""
representatives_selector.py

Utility function for extracting representative sequences for a true (gold standard) genome atomization pipeline.

Functions
---------
extract_representatives : Select class representatives using different strategies ("first", "mash").
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
import logging
import subprocess
import tempfile
from pathlib import Path

from atomization_scorer.data_processing import read_fasta, read_geese, write_fasta

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Representative Extraction
# --------------------------------------------------------------------------------------

def _extract_atom_sequence(
    row,
    genomes,
    class_id,
) -> tuple[str, int, int, int, str]:
    """
    Validate an atom row and return its genome name, atom number, coordinates, and subsequence.

    Parameters
    ----------
    row
        Atomization row containing the genome name, atom number, and coordinates.
    genomes
        Dictionary mapping genome names to genome sequences.
    class_id
        Atomization class identifier used for error reporting.

    Raises
    ------
    ValueError
        Raised if the genome referenced by the atom row is not present in the FASTA input.

    Returns
    -------
    tuple[str, int, int, int, str]
        Tuple containing:
        - genome name
        - atom number
        - start coordinate
        - end coordinate
        - atom subsequence
    """
    sequence_name = row["name"]
    atom_nr = int(row["atom_nr"])
    start = int(row["start"])
    end = int(row["end"])

    # Validate genome reference before slicing the atom subsequence.
    if sequence_name not in genomes:
        raise ValueError(
            f"Genome '{sequence_name}' referenced in atomization was not found in FASTA "
            f"(class={class_id}, atom_nr={atom_nr})."
        )

    subsequence = str(genomes[sequence_name][start:end])
    return sequence_name, atom_nr, start, end, subsequence


def extract_representatives(
    genomes_file: Path,
    atomization_file: Path,
    output_path: Path,
    mode: str = "mash"
) -> Path:
    """
    Extract one representative sequence for each atomization class.

    Parameters
    ----------
    genomes_file : Path
        Input genomes FASTA file containing the genomes sequences.
    atomization_file : Path
        Input GEESE file containing the predicted atomization.
    output_path : Path
        Path to the output FASTA file.
    mode : str, optional, default: "mash"
        Selection mode: "first" or "mash".

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.
    ValueError
        Raised if the mode is not one of: "first" or "mash" or if a genome
        sequence referenced in atomization is missing.
    subprocess.CalledProcessError
        Raised if the "mash" command fails during execution.

    Returns
    -------
    Path
        Generated representative FASTA file path.
    """
    if not genomes_file.is_file():
        raise FileNotFoundError(f"Genomes FASTA file not found: {genomes_file}")
    if not atomization_file.is_file():
        raise FileNotFoundError(f"Atomization file not found: {atomization_file}")
    if mode not in ("first", "mash"):
        raise ValueError("Mode must be one of: 'first', 'mash'")

    genomes = read_fasta(fasta_file=genomes_file)
    df_atoms = read_geese(geese_file=atomization_file)
    representatives = {}

    logger.info(
        "Extracting representatives with mode=%s from genomes=%s atomization=%s into output=%s",
        mode,
        genomes_file,
        atomization_file,
        output_path,
    )

    for class_id, group in df_atoms.groupby("class"):
        row = None
        logger.debug("Selecting representative for class=%s with %s candidate atoms", class_id, len(group))

        if mode == "first":
            row = group.iloc[0]
        elif mode == "mash":
            logger.info("Selecting representative for class=%s using Mash distances", class_id)
            # Create temporary FASTA for mash
            with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fa") as temporary_file:
                distances = {}
                atom_numbers = {}
                rows_by_header = {}
                for _, candidate_row in group.iterrows():
                    row = candidate_row
                    sequence_name, atom_nr, _start, _end, subsequence = _extract_atom_sequence(
                        row=row,
                        genomes=genomes,
                        class_id=class_id,
                    )
                    header = f"{sequence_name}_{atom_nr}"
                    temporary_file.write(f">{header}\n{subsequence}\n")
                    distances[header] = 0.0
                    atom_numbers[header] = atom_nr
                    rows_by_header[header] = row

                temporary_file.flush()
                temporary_path = Path(temporary_file.name)

            try:
                # Run mash
                mash_output = subprocess.check_output(
                    ["mash", "dist", "-i", str(temporary_path), str(temporary_path)],
                    text=True
                )

                # Parse mash output
                for line_number, line in enumerate(mash_output.splitlines(), start=1):
                    if not line.strip():
                        continue

                    fields = line.split("\t")
                    if len(fields) < 3:
                        raise ValueError(
                            f"Malformed Mash output at line {line_number}: expected at least 3 tab-separated fields."
                        )

                    reference, query = fields[0], fields[1]
                    if reference not in distances or query not in distances:
                        raise ValueError(
                            f"Malformed Mash output at line {line_number}: unexpected header '{reference}' or '{query}'."
                        )

                    try:
                        mash_distance = float(fields[2])
                    except ValueError as error:
                        raise ValueError(
                            f"Malformed Mash output at line {line_number}: invalid distance '{fields[2]}'."
                        ) from error

                    if reference != query:
                        distances[reference] += mash_distance
                        distances[query] += mash_distance

                selected_header = min(
                    distances,
                    key=lambda candidate_header: (
                        distances[candidate_header],
                        atom_numbers[candidate_header],
                    ),
                )
                tied_headers = [
                    header
                    for header, distance in distances.items()
                    if distance == distances[selected_header]
                ]
                if len(tied_headers) > 1:
                    logger.debug(
                        "Tie detected for class=%s among headers=%s; selected lowest atom_nr header=%s",
                        class_id,
                        tied_headers,
                        selected_header,
                    )

                row = rows_by_header[selected_header]

            finally:
                # Cleanup temporary file
                if temporary_path.exists():
                    temporary_path.unlink()

        # Extract sequence
        sequence_name, atom_nr, start, end, subsequence = _extract_atom_sequence(
            row=row,
            genomes=genomes,
            class_id=class_id,
        )
        representatives[f"{sequence_name}|class_{class_id}"] = subsequence
        logger.debug(
            "Selected representative for class=%s: genome=%s atom_nr=%s start=%s end=%s",
            class_id,
            sequence_name,
            atom_nr,
            start,
            end,
        )

    output_fasta = write_fasta(sequences=representatives, output_path=output_path)
    logger.info("Representative FASTA saved to %s", output_fasta)
    print(f"Representative FASTA saved to {output_fasta}")
    return output_fasta
