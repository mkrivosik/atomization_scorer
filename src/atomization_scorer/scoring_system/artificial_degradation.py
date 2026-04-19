"""
artificial_degradation.py

Provides functionality for generating artificially degraded atomization files.

Functions
---------
degrade_atomization : Randomly change atom classes for a selected fraction of atoms
                      and write the degraded atomization to a GEESE file.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
from __future__ import annotations
import logging
from pathlib import Path
import random

from atomization_scorer.data_processing import read_geese, write_geese

# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------
log = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------
# Artificial Degradation
# --------------------------------------------------------------------------------------
def degrade_atomization(
    atomization_file: Path,
    output_directory: Path,
    degradation_fraction: float,
    random_seed: int | None = None,
) -> Path:
    """
    Randomly change atom classes for a selected fraction of atoms and write the result.

    Parameters
    ----------
    atomization_file : Path
        Input GEESE file containing the atomization to degrade.
    output_directory : Path
        Output directory where the degraded GEESE file is written.
    degradation_fraction : float
        Fraction of atoms whose classes should be changed. Must be between 0.0 and 1.0.
    random_seed : int or None, optional, default=None
        Seed used for reproducible random degradation.

    Raises
    ------
    ValueError
        Raised if degradation_fraction is outside the interval [0.0, 1.0] or the
        atomization contains fewer than two distinct classes.

    Returns
    -------
    Path
        Written GEESE file path containing the degraded atomization.
    """
    if degradation_fraction < 0.0 or degradation_fraction > 1.0:
        raise ValueError("degradation_fraction must be between 0.0 and 1.0.")

    output_file = output_directory / "degraded_atomization.geese"
    log.info(
        "Generating degraded atomization from %s into %s with fraction=%s seed=%s",
        atomization_file,
        output_file,
        degradation_fraction,
        random_seed,
    )

    df_atoms = read_geese(geese_file=atomization_file).copy()
    distinct_classes = sorted(df_atoms["class"].unique().tolist())
    if len(distinct_classes) < 2:
        raise ValueError("Artificial degradation requires at least two distinct atomization classes.")

    n_atoms_to_change = round(len(df_atoms) * degradation_fraction)
    rng = random.Random(random_seed)

    if n_atoms_to_change > 0:
        selected_indices = rng.sample(df_atoms.index.tolist(), k=n_atoms_to_change)
        for index in selected_indices:
            original_class = df_atoms.at[index, "class"]
            candidate_classes = [atom_class for atom_class in distinct_classes if atom_class != original_class]
            df_atoms.at[index, "class"] = rng.choice(candidate_classes)

    log.info(
        "Changed classes for %s of %s atoms in degraded atomization",
        n_atoms_to_change,
        len(df_atoms),
    )
    return write_geese(df_atoms=df_atoms, output_path=output_file)
