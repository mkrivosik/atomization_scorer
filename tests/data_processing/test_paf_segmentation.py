"""
Tests for the resolve_paf_overlaps() and validate_non_overlapping_geese() functions.
"""

# --------------------------------------------------------------------------------------
# Imports
# --------------------------------------------------------------------------------------
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


def parse_paf_line(line: str) -> tuple[str, int, int, str]:
    """Return (query_name, query_start, query_end, target_name) from a PAF line."""
    fields = line.split("\t")
    return fields[0], int(fields[2]), int(fields[3]), fields[5]


# --------------------------------------------------------------------------------------
# Test: best-ranked alignment is accepted intact; lower-ranked is trimmed to surviving tail
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_accepts_best_and_trims_lower_ranked_to_non_overlapping_tail(tmp_path: Path):
    """resolve_paf_overlaps should accept the best-ranked alignment intact and trim each
    lower-ranked candidate to its non-overlapping portion."""
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
    assert len(lines) == 3
    _, start0, end0, target0 = parse_paf_line(lines[0])
    _, start1, end1, target1 = parse_paf_line(lines[1])
    _, start2, end2, target2 = parse_paf_line(lines[2])
    assert target0 == "rep|class_1" and start0 == 150 and end0 == 450
    assert target1 == "rep|class_2" and start1 == 100 and end1 == 150
    assert target2 == "rep|class_3" and start2 == 450 and end2 == 700


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
# Test: multi-overlap scenario trims all candidates, producing a non-overlapping covering
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_trims_all_candidates_in_multi_overlap_scenario(tmp_path: Path):
    """resolve_paf_overlaps should trim every overlapping candidate to its surviving portion
    when multiple candidates interact."""
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
    assert len(lines) == 4
    _, start0, end0, target0 = parse_paf_line(lines[0])
    _, start1, end1, target1 = parse_paf_line(lines[1])
    _, start2, end2, target2 = parse_paf_line(lines[2])
    _, start3, end3, target3 = parse_paf_line(lines[3])
    assert target0 == "rep|class_1" and start0 == 100 and end0 == 300
    assert target1 == "rep|class_4" and start1 == 420 and end1 == 520
    assert target2 == "rep|class_2" and start2 == 300 and end2 == 350
    assert target3 == "rep|class_3" and start3 == 350 and end3 == 420


# --------------------------------------------------------------------------------------
# Test: resolves overlaps independently per query genome including their trimmed tails
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_treats_queries_independently(tmp_path: Path):
    """resolve_paf_overlaps should resolve overlaps and trim candidates independently for each
    query genome."""
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
    assert len(lines) == 4
    name0, start0, end0, target0 = parse_paf_line(lines[0])
    name1, start1, end1, target1 = parse_paf_line(lines[1])
    name2, start2, end2, target2 = parse_paf_line(lines[2])
    name3, start3, end3, target3 = parse_paf_line(lines[3])
    assert name0 == "genome1" and target0 == "rep|class_1" and start0 == 100 and end0 == 300
    assert name1 == "genome1" and target1 == "rep|class_2" and start1 == 300 and end1 == 320
    assert name2 == "genome2" and target2 == "rep|class_3" and start2 == 100 and end2 == 300
    assert name3 == "genome2" and target3 == "rep|class_4" and start3 == 300 and end3 == 320


# --------------------------------------------------------------------------------------
# Test: fully contained alignment produces no surviving sub-intervals and is discarded
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_discards_fully_contained_alignment(tmp_path: Path):
    """resolve_paf_overlaps should fully discard a lower-ranked alignment that is completely
    contained within an already-accepted interval."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t0\t400\t+\trep|class_1\t400\t0\t400\t400\t400\t60\ttp:A:P\tdv:f:0.0\n"
        "genome1\t5000\t100\t300\t+\trep|class_2\t200\t0\t200\t200\t200\t40\ttp:A:P\tdv:f:0.0\n",
        filename="contained.paf",
    )
    output_file = tmp_path / "resolved_contained.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 1
    _, start0, end0, target0 = parse_paf_line(lines[0])
    assert target0 == "rep|class_1" and start0 == 0 and end0 == 400


# --------------------------------------------------------------------------------------
# Test: accepted interval in the interior of a candidate splits it into two surviving fragments
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_splits_alignment_into_two_surviving_intervals(tmp_path: Path):
    """resolve_paf_overlaps should split a lower-ranked alignment into two rows when a
    higher-ranked interval covers its interior."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t0\t700\t+\trep|class_2\t700\t0\t700\t650\t700\t40\ttp:A:P\tdv:f:0.01\n"
        "genome1\t5000\t200\t500\t+\trep|class_1\t300\t0\t300\t300\t300\t60\ttp:A:P\tdv:f:0.0\n",
        filename="split.paf",
    )
    output_file = tmp_path / "resolved_split.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 3
    _, start0, end0, target0 = parse_paf_line(lines[0])
    _, start1, end1, target1 = parse_paf_line(lines[1])
    _, start2, end2, target2 = parse_paf_line(lines[2])
    assert target0 == "rep|class_1" and start0 == 200 and end0 == 500
    assert target1 == "rep|class_2" and start1 == 0 and end1 == 200
    assert target2 == "rep|class_2" and start2 == 500 and end2 == 700


# --------------------------------------------------------------------------------------
# Test: minimum_alignment_length suppresses surviving sub-intervals that are too short
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_minimum_alignment_length_discards_short_surviving_fragment(tmp_path: Path):
    """resolve_paf_overlaps should discard a surviving sub-interval whose length is strictly
    below minimum_alignment_length."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t100\t500\t+\trep|class_1\t400\t0\t400\t400\t400\t60\ttp:A:P\tdv:f:0.0\n"
        "genome1\t5000\t480\t520\t+\trep|class_2\t40\t0\t40\t40\t40\t40\ttp:A:P\tdv:f:0.0\n",
        filename="short_fragment.paf",
    )
    output_file = tmp_path / "resolved_short.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file, minimum_alignment_length=30)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 1
    _, start0, end0, target0 = parse_paf_line(lines[0])
    assert target0 == "rep|class_1" and start0 == 100 and end0 == 500


# --------------------------------------------------------------------------------------
# Test: primary alignment wins over secondary; secondary is trimmed
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_prefers_primary_over_secondary(tmp_path: Path):
    """resolve_paf_overlaps should accept a primary alignment ahead of an overlapping secondary
    and trim the secondary to its non-overlapping portion."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t100\t300\t+\trep|class_1\t200\t0\t200\t180\t200\t40\ttp:A:S\tdv:f:0.0\n"
        "genome1\t5000\t120\t320\t+\trep|class_2\t200\t0\t200\t180\t200\t40\ttp:A:P\tdv:f:0.0\n",
        filename="primary_secondary.paf",
    )
    output_file = tmp_path / "resolved_primary_secondary.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 2
    _, start0, end0, target0 = parse_paf_line(lines[0])
    _, start1, end1, target1 = parse_paf_line(lines[1])
    assert target0 == "rep|class_2" and start0 == 120 and end0 == 320
    assert target1 == "rep|class_1" and start1 == 100 and end1 == 120


# --------------------------------------------------------------------------------------
# Test: quality tie-breakers accept the stronger alignment; weaker is trimmed
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
    """resolve_paf_overlaps should accept the stronger overlapping alignment intact and trim the
    weaker one to its non-overlapping tail across all quality tie-breaker dimensions."""
    paf_file = create_paf_file(tmp_path, contents, filename=filename)
    output_file = tmp_path / f"resolved_{filename}"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 2
    _, start0, end0, target0 = parse_paf_line(lines[0])
    _, start1, end1, target1 = parse_paf_line(lines[1])
    assert target0 == "rep|class_2" and start0 == 120 and end0 == 320
    assert target1 == "rep|class_1" and start1 == 100 and end1 == 120


# --------------------------------------------------------------------------------------
# Test: deterministic fallback resolves fully tied candidates; runner-up keeps surviving tail
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_uses_deterministic_fallback_for_full_ties(tmp_path: Path):
    """resolve_paf_overlaps should break full quality ties deterministically by target name and
    trim the runner-up to its surviving tail."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t120\t320\t+\trep|class_b\t200\t0\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n"
        "genome1\t5000\t100\t300\t+\trep|class_a\t200\t0\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n",
        filename="full_tie.paf",
    )
    output_file = tmp_path / "resolved_full_tie.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 2
    _, start0, end0, target0 = parse_paf_line(lines[0])
    _, start1, end1, target1 = parse_paf_line(lines[1])
    assert target0 == "rep|class_a" and start0 == 100 and end0 == 300
    assert target1 == "rep|class_b" and start1 == 300 and end1 == 320


# --------------------------------------------------------------------------------------
# Test: unknown alignment type ranks after known types; lower-ranked is trimmed
# --------------------------------------------------------------------------------------
def test_resolve_paf_overlaps_prefers_known_alignment_type_over_unexpected_type(tmp_path: Path):
    """resolve_paf_overlaps should rank an unexpected alignment type after all known types and
    trim it to its surviving portion."""
    paf_file = create_paf_file(
        tmp_path,
        "genome1\t5000\t100\t300\t+\trep|class_1\t190\t200\t200\t190\t200\t60\ttp:A:X\tdv:f:0.01\n"
        "genome1\t5000\t120\t320\t+\trep|class_2\t190\t200\t200\t190\t200\t60\ttp:A:P\tdv:f:0.01\n",
        filename="unexpected_type.paf",
    )
    output_file = tmp_path / "resolved_unexpected_type.paf"

    resolve_paf_overlaps(paf_file=paf_file, output_file=output_file)

    lines = output_file.read_text().splitlines()
    assert len(lines) == 2
    _, start0, end0, target0 = parse_paf_line(lines[0])
    _, start1, end1, target1 = parse_paf_line(lines[1])
    assert target0 == "rep|class_2" and start0 == 120 and end0 == 320
    assert target1 == "rep|class_1" and start1 == 100 and end1 == 120


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
