"""
paf_to_geese.py

Utility function for converting PAF files to GEESE files.

Functions
---------
paf_to_geese : Convert a PAF file to a GEESE file.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
from pathlib import Path

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------
# Converting PAF To GEESE
# --------------------------------------------------------------------------------------
def paf_to_geese(paf_file: Path, output_file: Path) -> Path:
    """
    Convert a PAF file to a GEESE file.

    Parameters
    ----------
    paf_file : Path
        Path to the input PAF file from minimap2.
    output_file : Path
        Path to the output GEESE file.

    Raises
    ------
    FileNotFoundError
        Raised if the input PAF file does not exist.
    ValueError
        Raised if the PAF file contains malformed content.

    Returns
    -------
    Path
        Written GEESE file path.
    """
    if not paf_file.is_file():
        raise FileNotFoundError(f"PAF file {paf_file} not found.")

    log.info("Converting PAF file %s to GEESE file %s", paf_file, output_file)

    converted_rows = 0
    converted_entries = []

    with paf_file.open("r") as paf:
        for line_number, line in enumerate(paf, start=1):
            fields = line.strip().split("\t")

            if len(fields) < 12:
                raise ValueError(
                    f"Malformed PAF row at line {line_number}: "
                    "expected at least 12 tab-separated fields."
                )

            query_name = fields[0]
            target_name = fields[5]

            try:
                query_start = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"Invalid query start at line {line_number}: "
                    "expected an integer in field 3."
                ) from error

            try:
                query_end = int(fields[3])
            except ValueError as error:
                raise ValueError(
                    f"Invalid query end at line {line_number}: "
                    "expected an integer in field 4."
                ) from error

            # The genome name comes from the query; the class is encoded in the representative header.
            if "|class_" not in target_name:
                raise ValueError(
                    f"Missing class tag at line {line_number}: "
                    "target header must contain '|class_'."
                )

            _, class_id = target_name.split("|class_")
            converted_entries.append((query_name, class_id, query_start, query_end))
            converted_rows += 1

    converted_entries.sort(key=lambda entry: (entry[0], entry[2], entry[3], entry[1]))

    with output_file.open("w") as geese:
        geese.write("#name\tatom_nr\tclass\tstart\tend\n")
        for atom_nr, (query_name, class_id, query_start, query_end) in enumerate(converted_entries, start=1):
            geese.write(f"{query_name}\t{atom_nr}\t{class_id}\t{query_start}\t{query_end}\n")

    log.info("GEESE file saved to %s with %s converted rows", output_file, converted_rows)

    return output_file
