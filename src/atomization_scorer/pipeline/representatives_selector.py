"""
representatives_selector.py

Utility function for extracting representative sequences for true (gold standard) genome atomization pipeline.

Functions
---------
extract_representatives : Select class representatives using different strategies ("first", "mash").
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import subprocess
import tempfile
from pathlib import Path

from atomization_scorer.data_processing import read_fasta, read_geese, write_fasta

# ---------------------------------------------------------------------
# Representative Extraction
# ---------------------------------------------------------------------

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
        Path to the output directory where resulting FASTA file is stored.
    mode : str, optional, default: "mash"
        Selection mode: "first" or "mash".

    Raises
    ------
    FileNotFoundError
        Raised if genomes_file or atomization_file do not exist.
    ValueError
        Raised if mode is not one of: "first" or "mash" or if a genome
        sequence referenced in atomization is missing.
    subprocess.CalledProcessError
        Raised if the "mash" command fails during execution.

    Returns
    -------
    Path
        Path to generated representative FASTA file.
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

    for class_id, group in df_atoms.groupby("class"):
        row = None
        if mode == "first":
            row = group.iloc[0]
        elif mode == "mash":
            # Create temporary FASTA for mash
            with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fa") as temporary_file:
                distances = {}
                for _, r in group.iterrows():
                    sequence_name = r["name"]
                    atom_nr = int(r["atom_nr"])
                    start = int(r["start"])
                    end = int(r["end"])

                    subsequence = genomes[sequence_name][start:end+1]
                    header = f"{sequence_name}_{atom_nr}"
                    temporary_file.write(f">{header}\n{subsequence}\n")
                    distances[header] = 0.0

                temporary_file.flush()
                temporary_path = Path(temporary_file.name)

            try:
                # Run mash
                mash_output = subprocess.check_output(
                    ["mash", "dist", "-i", str(temporary_path), str(temporary_path)],
                    text=True
                )

                # Parse mash output
                for line in mash_output.strip().split("\n"):
                    fields = line.split("\t")
                    reference, query, mash_distance = fields[0], fields[1], float(fields[2])
                    if reference != query:
                        distances[reference] += mash_distance
                        distances[query] += mash_distance

                selected_header = min(distances, key=distances.get)
                selected_atom_nr = int(selected_header.split("_")[-1])
                row = group[group["atom_nr"] == selected_atom_nr].iloc[0]

            finally:
                # Cleanup temporary file
                if temporary_path.exists():
                    temporary_path.unlink()

        # Extract sequence
        sequence_name = row["name"]
        start = int(row["start"])
        end = int(row["end"])

        subsequence = genomes[sequence_name][start:end+1]
        representatives[f"{sequence_name}|class_{class_id}"] = str(subsequence)

    output_fasta = write_fasta(sequences=representatives, output_path=output_path)
    print(f"Representative FASTA saved to {output_fasta}")
    return output_fasta
