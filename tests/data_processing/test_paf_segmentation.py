"""
Tests for the paf_segmentation() function.
"""

from pathlib import Path

import pytest

from atomization_scorer import resolve_paf_overlaps, validate_non_overlapping_geese


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def create_paf_file(tmp_path: Path, contents: str, filename: str = "input.paf") -> Path:
    """Create a custom PAF file for segmentation tests."""
    paf_file = tmp_path / filename
    paf_file.write_text(contents)
    return paf_file


# --------------------------------------------------------------------------------------
# Test: keeps best non-overlapping alignments
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_keeps_best_non_overlapping_alignments(tmp_path: Path):
    """resolve_paf_overlaps should retain the best-ranked non-overlapping alignment on a query."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t100\t400\t+\trep|class_2\t300\t0\t300\t290\t300\t60\ttp:A:P\tdv:f:0.0333\n"
        "genome1\t5000\t150\t450\t+\trep|class_1\t300\t0\t300\t300\t300\t60\ttp:A:P\tdv:f:0.0\n"
        "genome1\t5000\t450\t700\t+\trep|class_3\t250\t0\t250\t250\t250\t50\ttp:A:P\tdv:f:0.0\n",
        filename="overlap.paf",
    )
    output_file = tmp_path / "resolved.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 2
    assert "rep|class_1" in lines[0]
    assert "rep|class_3" in lines[1]


# --------------------------------------------------------------------------------------
# Test: allows half-open adjacency
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_allows_half_open_adjacency(tmp_path: Path):
    """resolve_paf_overlaps should allow adjacent half-open query intervals."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t0\t100\t+\trep|class_1\t100\t0\t100\t100\t100\t60\ttp:A:P\tdv:f:0.0\n"
        "genome1\t5000\t100\t200\t+\trep|class_2\t100\t0\t100\t100\t100\t60\ttp:A:P\tdv:f:0.0\n",
        filename="adjacent.paf",
    )
    output_file = tmp_path / "resolved_adjacent.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 2


# --------------------------------------------------------------------------------------
# Test: resolves multiple overlaps within one query
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_keeps_clean_non_overlapping_set_for_multiple_overlaps(tmp_path: Path):
    """resolve_paf_overlaps should keep a clean non-overlapping set when several candidates overlap."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t100\t300\t+\trep|class_1\t200\t0\t200\t200\t200\t60\ttp:A:P\tdv:f:0.0\n"
        "genome1\t5000\t150\t350\t+\trep|class_2\t200\t0\t200\t180\t200\t50\ttp:A:P\tdv:f:0.1\n"
        "genome1\t5000\t280\t420\t+\trep|class_3\t140\t0\t140\t100\t140\t40\ttp:A:S\tdv:f:0.2\n"
        "genome1\t5000\t420\t520\t+\trep|class_4\t100\t0\t100\t100\t100\t60\ttp:A:P\tdv:f:0.0\n",
        filename="multiple_overlap.paf",
    )
    output_file = tmp_path / "resolved_multiple.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 2
    assert "rep|class_1" in lines[0]
    assert "rep|class_4" in lines[1]


# --------------------------------------------------------------------------------------
# Test: resolves overlaps independently per query
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_treats_queries_independently(tmp_path: Path):
    """resolve_paf_overlaps should resolve overlaps independently for different query genomes."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t100\t300\t+\trep|class_1\t200\t0\t200\t200\t200\t60\ttp:A:P\tdv:f:0.0\n"
        "genome1\t5000\t150\t320\t+\trep|class_2\t170\t0\t170\t150\t170\t50\ttp:A:S\tdv:f:0.1\n"
        "genome2\t5000\t100\t300\t+\trep|class_3\t200\t0\t200\t200\t200\t60\ttp:A:P\tdv:f:0.0\n"
        "genome2\t5000\t150\t320\t+\trep|class_4\t170\t0\t170\t150\t170\t50\ttp:A:S\tdv:f:0.1\n",
        filename="independent_queries.paf",
    )
    output_file = tmp_path / "resolved_independent.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 2
    assert lines[0].startswith("genome1\t")
    assert "rep|class_1" in lines[0]
    assert lines[1].startswith("genome2\t")
    assert "rep|class_3" in lines[1]


# --------------------------------------------------------------------------------------
# Test: primary alignment wins over secondary alignment
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_prefers_primary_over_secondary(tmp_path: Path):
    """resolve_paf_overlaps should prefer a primary overlapping alignment over a secondary one."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t100\t300\t+\trep|class_1\t200\t0\t200\t180\t200\t40\ttp:A:S\tdv:f:0.0\n"
        "genome1\t5000\t120\t320\t+\trep|class_2\t200\t0\t200\t180\t200\t40\ttp:A:P\tdv:f:0.0\n",
        filename="primary_secondary.paf",
    )
    output_file = tmp_path / "resolved_primary_secondary.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 1
    assert "rep|class_2" in lines[0]


# --------------------------------------------------------------------------------------
# Test: quality tie-breakers prefer the stronger alignment
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("filename", "contents"),
    [
        (
            "mapq_tie.paf",
            "genome1\t5000\t100\t300\t+\trep|class_1\t200\t0\t200\t190\t200\t20\ttp:A:P\tdv:f:0.01\n"
            "genome1\t5000\t120\t320\t+\trep|class_2\t200\t0\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n",
        ),
        (
            "similarity_tie.paf",
            "genome1\t5000\t100\t300\t+\trep|class_1\t200\t0\t200\t180\t200\t60\ttp:A:P\tdv:f:0.10\n"
            "genome1\t5000\t120\t320\t+\trep|class_2\t200\t0\t200\t180\t200\t60\ttp:A:P\tdv:f:0.01\n",
        ),
        (
            "length_tie.paf",
            "genome1\t5000\t100\t260\t+\trep|class_1\t160\t0\t160\t150\t160\t60\ttp:A:P\tdv:f:0.01\n"
            "genome1\t5000\t120\t320\t+\trep|class_2\t200\t0\t200\t188\t200\t60\ttp:A:P\tdv:f:0.01\n",
        ),
        (
            "matches_tie.paf",
            "genome1\t5000\t100\t300\t+\trep|class_1\t200\t0\t200\t180\t200\t60\ttp:A:P\tdv:f:0.01\n"
            "genome1\t5000\t120\t320\t+\trep|class_2\t200\t0\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n",
        ),
    ],
)
def test_resolve_paf_overlaps_prefers_stronger_quality_tie_breaker(
    tmp_path: Path,
    filename: str,
    contents: str,
):
    """resolve_paf_overlaps should prefer the stronger overlapping alignment across quality tie-breakers."""
    paf_file = create_paf_file(tmp_path, contents, filename=filename)
    output_file = tmp_path / f"resolved_{filename}"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 1
    assert "rep|class_2" in lines[0]


# --------------------------------------------------------------------------------------
# Test: deterministic fallback resolves fully tied candidates
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_uses_deterministic_fallback_for_full_ties(tmp_path: Path):
    """resolve_paf_overlaps should use deterministic ordering when quality fields are fully tied."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t120\t320\t+\trep|class_b\t200\t0\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n"
        "genome1\t5000\t100\t300\t+\trep|class_a\t200\t0\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n",
        filename="full_tie.paf",
    )
    output_file = tmp_path / "resolved_full_tie.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 1
    assert "rep|class_a" in lines[0]


# --------------------------------------------------------------------------------------
# Test: unknown alignment type ranks after known types
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_prefers_known_alignment_type_over_unexpected_type(tmp_path: Path):
    """resolve_paf_overlaps should rank unexpected alignment types after known types."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t100\t300\t+\trep|class_1\t190\t200\t200\t190\t200\t60\ttp:A:X\tdv:f:0.01\n"
        "genome1\t5000\t120\t320\t+\trep|class_2\t190\t200\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n",
        filename="unexpected_type.paf",
    )
    output_file = tmp_path / "resolved_unexpected_type.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 1
    assert "rep|class_2" in lines[0]


# --------------------------------------------------------------------------------------
# Test: missing PAF file is rejected
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_missing_file(tmp_path: Path):
    """resolve_paf_overlaps should raise FileNotFoundError if the PAF file does not exist."""
    missing_file = tmp_path / "missing.paf"

    with pytest.raises(FileNotFoundError):
        resolve_paf_overlaps(paf_file=missing_file, output_file=tmp_path / "resolved.paf")


# --------------------------------------------------------------------------------------
# Test: empty PAF file is handled
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_empty_file(tmp_path: Path):
    """resolve_paf_overlaps should create an empty resolved file when the input PAF is empty."""
    paf_file = create_paf_file(tmp_path, "", filename="empty.paf")
    output_file = tmp_path / "resolved_empty.paf"

    result = resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    assert result == output_file
    assert output_file.is_file()
    assert output_file.read_text() == ""


# --------------------------------------------------------------------------------------
# Test: malformed PAF rows are rejected
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("contents", "expected_message"),
    [
        ("genome1\t5000\t100\n", "expected at least 12 tab-separated fields"),
        (
            "genome1\t5000\tstart\t300\t+\trep|class_1\t200\t0\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n",
            "invalid numeric value in mandatory fields",
        ),
        (
            "genome1\t5000\t300\t300\t+\trep|class_1\t200\t0\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n",
            "query start must be smaller than query end",
        ),
    ],
)
def test_resolve_paf_overlaps_rejects_malformed_rows(
    tmp_path: Path,
    contents: str,
    expected_message: str,
):
    """resolve_paf_overlaps should fail clearly on malformed PAF rows."""
    paf_file = create_paf_file(tmp_path, contents, filename="malformed.paf")

    with pytest.raises(ValueError, match=expected_message):
        resolve_paf_overlaps(paf_file=paf_file, output_file=tmp_path / "resolved.paf")


# --------------------------------------------------------------------------------------
# Test: rejects overlapping half-open GEESE intervals
# --------------------------------------------------------------------------------------
def test_validate_non_overlapping_geese_rejects_overlap(tmp_path: Path):
    """validate_non_overlapping_geese should reject overlapping half-open intervals."""
    geese_file = tmp_path / "true_atomization.geese"
    geese_file.write_text(
        "#name\tclass\tstart\tend\n"
        "genome1\t1\t0\t100\n"
        "genome1\t2\t90\t150\n"
    )

    with pytest.raises(ValueError, match="Overlapping true atoms detected"):
        validate_non_overlapping_geese(geese_file)


# --------------------------------------------------------------------------------------
# Test: allows half-open GEESE adjacency
# --------------------------------------------------------------------------------------
def test_validate_non_overlapping_geese_allows_half_open_adjacency(tmp_path: Path):
    """validate_non_overlapping_geese should allow adjacent half-open intervals."""
    geese_file = tmp_path / "true_atomization.geese"
    geese_file.write_text(
        "#name\tclass\tstart\tend\n"
        "genome1\t1\t0\t100\n"
        "genome1\t2\t100\t150\n"
    )

    validate_non_overlapping_geese(geese_file)
