"""
paf_to_geese.py

Utility function for converting PAF files to GEESE files.

Functions
---------
paf_to_geese : Converts PAF files to GEESE files.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
from pathlib import Path

# ---------------------------------------------------------------------
# Converting PAF To GEESE
# ---------------------------------------------------------------------

def paf_to_geese(paf_file: Path, output_file: Path) -> Path:
    """
    Converts PAF files to GEESE files.

    Parameters
    ----------
    paf_file : Path
        Path to input a PAF file from minimap2.
    output_file : Path
        File where the GEESE file should be written.

    Raises
    ------
    FileNotFoundError
        Raised if the input PAF file does not exist.

    Returns
    -------
    Path
        The Path to the GEESE file.
    """
    if not paf_file.is_file():
        raise FileNotFoundError(f"PAF file {paf_file} not found.")

    with paf_file.open("r") as paf, output_file.open("w") as geese:
        geese.write("#name\tclass\tstart\tend\n")

        for line in paf:
            fields = line.strip().split("\t")

            query_name = fields[0]
            target_name = fields[5]
            query_start = int(fields[2])
            query_end = int(fields[3])

            # The genome name comes from the query; the class is encoded in the representative header.
            if "|class_" in target_name:
                _, class_id = target_name.split("|class_")

            geese.write(f"{query_name}\t{class_id}\t{query_start}\t{query_end}\n")

    return output_file
