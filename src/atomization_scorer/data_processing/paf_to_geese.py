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
        Path to input PAF file from minimap2.
    output_file : Path
        File where the GEESE file should be written.

    Raises
    ------
    FileNotFoundError
        Raised if the input PAF file does not exist.

    Returns
    -------
    Path
        Path to the GEESE file.
    """
    if not paf_file.is_file():
        raise FileNotFoundError(f"PAF file {paf_file} not found.")

    with paf_file.open("r") as paf, output_file.open("w") as geese:
        geese.write("#name\tclass\tstart\tend\n")

        for line in paf:
            fields = line.strip().split("\t")

            target_name = fields[5]
            target_start = int(fields[7])
            target_end = int(fields[8])

            # Split target_name into name and class
            if "|class_" in target_name:
                name, class_id = target_name.split("|class_")
            else:
                name = target_name
                class_id = "0"

            geese.write(f"{name}\t{class_id}\t{target_start}\t{target_end}\n")

    return output_file
