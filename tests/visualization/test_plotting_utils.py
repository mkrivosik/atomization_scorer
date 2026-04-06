"""
Tests for interactive visualization helper functions.
"""

import pandas as pd
import pytest

from atomization_scorer.visualization import plotting_utils as pu


# --------------------------------------------------------------------------------------
# Test: HTML is the only supported visualization output format
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("output_format", "expected"),
    [
        ("html", "html"),
        ("HTML", "html"),
    ],
)
def test_normalize_output_format_supported_values(output_format: str, expected: str):
    """normalize_output_format should lowercase supported interactive formats."""
    assert pu.normalize_output_format(output_format) == expected


def test_normalize_output_format_rejects_unsupported_value():
    """normalize_output_format should reject legacy static figure formats."""
    with pytest.raises(ValueError, match="Unsupported output format"):
        pu.normalize_output_format("png")


# --------------------------------------------------------------------------------------
# Test: initial window follows target row count and min/max bounds
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("genome_length", "target_rows", "minimum", "maximum", "expected"),
    [
        (20_000, 20, 10_000, 250_000, 10_000),
        (10_000_000, 20, 10_000, 250_000, 250_000),
        (3_000, 20, 10_000, 250_000, 3_000),
    ],
)
def test_compute_initial_window(
    genome_length: int,
    target_rows: int,
    minimum: int,
    maximum: int,
    expected: int,
):
    """compute_initial_window should clamp the visible window to the configured bounds."""
    assert pu.compute_initial_window(genome_length, target_rows, minimum, maximum) == expected


# --------------------------------------------------------------------------------------
# Test: atom extraction infers genome-global display numbers
# --------------------------------------------------------------------------------------
def test_get_atoms_for_genome_infers_global_display_numbers():
    """get_atoms_for_genome should assign genome-global display numbers in start order."""
    df = pd.DataFrame(
        {
            "name": ["genome1", "genome1", "genome1"],
            "class": ["A", "A", "B"],
            "start": [100, 400, 250],
            "end": [200, 500, 350],
        }
    )

    atoms = pu.get_atoms_for_genome(
        df=df,
        genome_name="genome1",
        genome_length=1_000,
        label="True",
        source="true",
    )

    assert atoms == [
        {
            "genome_name": "genome1",
            "source": "true",
            "class_id": "A",
            "atom_number": 1,
            "atom_id": "A:1",
            "start": 100,
            "end": 200,
            "length": 100,
        },
        {
            "genome_name": "genome1",
            "source": "true",
            "class_id": "B",
            "atom_number": 3,
            "atom_id": "B:3",
            "start": 250,
            "end": 350,
            "length": 100,
        },
        {
            "genome_name": "genome1",
            "source": "true",
            "class_id": "A",
            "atom_number": 2,
            "atom_id": "A:2",
            "start": 400,
            "end": 500,
            "length": 100,
        },
    ]


# --------------------------------------------------------------------------------------
# Test: explicit atom numbers are preserved when present
# --------------------------------------------------------------------------------------
def test_get_atoms_for_genome_preserves_explicit_atom_numbers():
    """get_atoms_for_genome should use an explicit atom_number column when available."""
    df = pd.DataFrame(
        {
            "name": ["genome1", "genome1"],
            "class": ["A", "A"],
            "atom_number": [3, 7],
            "start": [100, 400],
            "end": [200, 500],
        }
    )

    atoms = pu.get_atoms_for_genome(
        df=df,
        genome_name="genome1",
        genome_length=1_000,
        label="Predicted",
        source="predicted",
    )

    assert [atom["atom_number"] for atom in atoms] == [3, 7]


# --------------------------------------------------------------------------------------
# Test: atom_nr is preferred for displayed atom numbers when present in GEESE
# --------------------------------------------------------------------------------------
def test_get_atoms_for_genome_prefers_atom_nr_for_display_numbers():
    """get_atoms_for_genome should display original GEESE atom_nr values when they are present."""
    df = pd.DataFrame(
        {
            "name": ["genome1", "genome1"],
            "class": ["A", "B"],
            "atom_nr": [11, 42],
            "start": [100, 400],
            "end": [200, 500],
        }
    )

    atoms = pu.get_atoms_for_genome(
        df=df,
        genome_name="genome1",
        genome_length=1_000,
        label="Predicted",
        source="predicted",
    )

    assert [atom["atom_number"] for atom in atoms] == [11, 42]


# --------------------------------------------------------------------------------------
# Test: invalid intervals still raise and overlapping intervals only warn
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("df", "expected_message"),
    [
        (
            pd.DataFrame({"name": ["genome1"], "class": ["A"], "start": [-1], "end": [5]}),
            "contains a negative coordinate",
        ),
        (
            pd.DataFrame({"name": ["genome1"], "class": ["A"], "start": [5], "end": [5]}),
            "must satisfy start < end",
        ),
        (
            pd.DataFrame({"name": ["genome1"], "class": ["A"], "start": [0], "end": [2_000]}),
            "ends outside genome length",
        ),
    ],
)
def test_get_atoms_for_genome_rejects_invalid_intervals(df: pd.DataFrame, expected_message: str):
    """get_atoms_for_genome should reject invalid coordinates."""
    with pytest.raises(ValueError, match=expected_message):
        pu.get_atoms_for_genome(
            df=df,
            genome_name="genome1",
            genome_length=1_000,
            label="True",
            source="true",
        )


def test_get_atoms_for_genome_logs_overlapping_intervals(caplog: pytest.LogCaptureFixture):
    """get_atoms_for_genome should warn and continue when intervals overlap."""
    df = pd.DataFrame(
        {
            "name": ["genome1", "genome1"],
            "class": ["A", "B"],
            "start": [100, 150],
            "end": [200, 250],
        }
    )

    with caplog.at_level("WARNING"):
        atoms = pu.get_atoms_for_genome(
            df=df,
            genome_name="genome1",
            genome_length=1_000,
            label="Predicted",
            source="predicted",
        )

    assert len(atoms) == 2
    assert "visualization will continue" in caplog.text


# --------------------------------------------------------------------------------------
# Test: class colors are deterministic and matching pairs use class plus overlap
# --------------------------------------------------------------------------------------
def test_build_class_color_map_is_deterministic():
    """build_class_color_map should assign stable colors from sorted class names."""
    color_map = pu.build_class_color_map(["B", "A", "B"])

    assert set(color_map) == {"A", "B"}
    assert color_map["A"] != color_map["B"]


def test_pair_atoms_matches_same_class_atoms_that_overlap():
    """pair_atoms should match true and predicted atoms when they share a class and overlap in genome coordinates."""
    true_atoms = [
        pu.AtomRecord(
            genome_name="genome1",
            source="true",
            class_id="A",
            atom_number=1,
            atom_id="A:1",
            start=10,
            end=20,
            length=10,
        ),
        pu.AtomRecord(
            genome_name="genome1",
            source="true",
            class_id="A",
            atom_number=2,
            atom_id="A:2",
            start=30,
            end=40,
            length=10,
        ),
    ]
    predicted_atoms = [
        pu.AtomRecord(
            genome_name="genome1",
            source="predicted",
            class_id="A",
            atom_number=1,
            atom_id="A:1",
            start=12,
            end=22,
            length=10,
        ),
        pu.AtomRecord(
            genome_name="genome1",
            source="predicted",
            class_id="A",
            atom_number=2,
            atom_id="A:2",
            start=34,
            end=44,
            length=10,
        ),
        pu.AtomRecord(
            genome_name="genome1",
            source="predicted",
            class_id="B",
            atom_number=3,
            atom_id="B:3",
            start=30,
            end=40,
            length=10,
        ),
    ]

    matched_pairs, unmatched_true, unmatched_predicted = pu.pair_atoms(true_atoms, predicted_atoms)

    assert matched_pairs == [
        (true_atoms[0], predicted_atoms[0]),
        (true_atoms[1], predicted_atoms[1]),
    ]
    assert unmatched_true == []
    assert unmatched_predicted == [predicted_atoms[2]]


def test_pair_atoms_emits_all_same_class_overlaps_for_visualization():
    """pair_atoms should emit every same-class overlap pair so split atoms can connect to one larger partner."""
    true_atoms = [
        pu.AtomRecord(
            genome_name="genome1",
            source="true",
            class_id="A",
            atom_number=1,
            atom_id="A:1",
            start=10,
            end=50,
            length=40,
        ),
    ]
    predicted_atoms = [
        pu.AtomRecord(
            genome_name="genome1",
            source="predicted",
            class_id="A",
            atom_number=1,
            atom_id="A:1",
            start=12,
            end=20,
            length=8,
        ),
        pu.AtomRecord(
            genome_name="genome1",
            source="predicted",
            class_id="A",
            atom_number=2,
            atom_id="A:2",
            start=30,
            end=40,
            length=10,
        ),
    ]

    matched_pairs, unmatched_true, unmatched_predicted = pu.pair_atoms(true_atoms, predicted_atoms)

    assert matched_pairs == [
        (true_atoms[0], predicted_atoms[0]),
        (true_atoms[0], predicted_atoms[1]),
    ]
    assert unmatched_true == []
    assert unmatched_predicted == []


# --------------------------------------------------------------------------------------
# Test: filesystem-safe output names are generated for genome files
# --------------------------------------------------------------------------------------
def test_sanitize_output_stem():
    """sanitize_output_stem should preserve safe characters and hash empty results."""
    assert pu.sanitize_output_stem("genome.1") == "genome.1"
    assert pu.sanitize_output_stem("genome/1") == "genome_1"
    assert pu.sanitize_output_stem("|||").startswith("genome_")
