"""
Rep-vs-rep minimap2 alignment for large P-P signature class pairs.

For each unique class pair in large_pp_signatures.tsv:
  1. Identify the representative sequence for each class from the filtered PAF.
  2. Extract those sequences from the representatives FASTA.
  3. Run minimap2 representative_a vs representative_b.
  4. Report aligned fraction relative to the shorter representative to decide whether
     the two classes share substantial sequence (class ambiguity) or not.
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO

from tsv_utils import write_tsv

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
STRONG_AMBIGUITY_THRESHOLD = 0.5
MODERATE_SIMILARITY_THRESHOLD = 0.2
SUMMARY_FIELDS = [
    "class_a",
    "class_b",
    "representative_a",
    "representative_b",
    "representative_a_length",
    "representative_b_length",
    "aligned_length",
    "identity",
    "fraction_of_shorter_representative",
    "n_blocks",
    "strands",
    "verdict",
]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
def extract_class_to_representative(paf_path: Path) -> dict[str, str]:
    """Return mapping of class_id -> full representative name from PAF target field."""
    class_to_representative = {}
    with paf_path.open() as file:
        for line in file:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            target_name = fields[5]
            if "|class_" not in target_name:
                continue
            class_id = target_name.split("|class_")[1]
            class_to_representative.setdefault(class_id, target_name)
    return class_to_representative


def read_class_pairs(signatures_path: Path) -> list[tuple[str, str]]:
    """Read unique class pairs from large_pp_signatures.tsv."""
    pairs = []
    seen = set()
    with signatures_path.open() as file:
        for index, line in enumerate(file):
            if index < 2:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            pair = (fields[0].strip(), fields[1].strip())
            if pair not in seen:
                pairs.append(pair)
                seen.add(pair)
    return pairs


def write_single_fasta(path: Path, name: str, sequence: str) -> None:
    """Write one sequence to a FASTA file."""
    with path.open("w") as file:
        file.write(f">{name}\n")
        for start in range(0, len(sequence), 80):
            file.write(sequence[start:start + 80] + "\n")


# ---------------------------------------------------------------------------
# minimap2
# ---------------------------------------------------------------------------
def run_minimap2(query_fasta: Path, target_fasta: Path) -> list[str]:
    """Run minimap2 and return PAF lines."""
    result = subprocess.run(
        ["minimap2", "-x", "asm20", "-c", "-p", "0.1", str(target_fasta), str(query_fasta)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed:\n{result.stderr}")
    return [line for line in result.stdout.splitlines() if line]


# ---------------------------------------------------------------------------
# Result Parsing
# ---------------------------------------------------------------------------
def parse_alignment_result(
    paf_lines: list[str],
    representative_a_length: int,
    representative_b_length: int,
) -> dict:
    """Summarize minimap2 PAF output into alignment statistics."""
    if not paf_lines:
        return {
            "aligned_length": 0,
            "identity": 0.0,
            "fraction_of_shorter_representative": 0.0,
            "n_blocks": 0,
            "strands": "none",
            "verdict": "no_alignment",
        }

    total_matches = 0
    total_aligned = 0
    strands = set()
    for line in paf_lines:
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        total_matches += int(fields[9])
        total_aligned += int(fields[10])
        strands.add(fields[4])

    shorter_representative = min(representative_a_length, representative_b_length)
    fraction = total_aligned / shorter_representative if shorter_representative > 0 else 0.0
    identity = total_matches / total_aligned if total_aligned > 0 else 0.0

    if fraction >= STRONG_AMBIGUITY_THRESHOLD:
        verdict = "class_ambiguity"
    elif fraction >= MODERATE_SIMILARITY_THRESHOLD:
        verdict = "moderate_similarity"
    else:
        verdict = "representatives_distinct"

    if "-" in strands and verdict == "class_ambiguity":
        verdict = "inverted_repeat"

    return {
        "aligned_length": total_aligned,
        "identity": round(identity, 3),
        "fraction_of_shorter_representative": round(fraction, 3),
        "n_blocks": len(paf_lines),
        "strands": ",".join(sorted(strands)),
        "verdict": verdict,
    }


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(
    signatures_path: Path,
    paf_path: Path,
    representatives_fasta: Path,
    output_directory: Path,
) -> None:
    """For each class pair in signatures, run rep-vs-rep minimap2 and write summary."""
    class_pairs = read_class_pairs(signatures_path)
    class_to_representative = extract_class_to_representative(paf_path)
    needed_classes = {cls for pair in class_pairs for cls in pair}
    missing = needed_classes - set(class_to_representative)
    if missing:
        print(f"Warning: no representative found in PAF for classes: {sorted(missing)}")

    needed_representatives = {class_to_representative[cls] for cls in needed_classes if cls in class_to_representative}
    sequences = {
        record.id: str(record.seq)
        for record in SeqIO.parse(representatives_fasta, "fasta")
        if record.id in needed_representatives
    }

    missing_sequences = needed_representatives - set(sequences)
    if missing_sequences:
        print(f"Warning: sequences not found in FASTA for: {sorted(missing_sequences)}")

    output_directory.mkdir(parents=True, exist_ok=True)
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for class_a, class_b in class_pairs:
            representative_a_name = class_to_representative.get(class_a)
            representative_b_name = class_to_representative.get(class_b)

            if not representative_a_name or not representative_b_name:
                print(f"  Skipping {class_a}/{class_b}: representative name missing")
                continue
            if representative_a_name not in sequences or representative_b_name not in sequences:
                print(f"  Skipping {class_a}/{class_b}: sequence not in FASTA")
                continue

            sequence_a = sequences[representative_a_name]
            sequence_b = sequences[representative_b_name]

            fasta_a = tmp_path / f"class_{class_a}.fa"
            fasta_b = tmp_path / f"class_{class_b}.fa"
            write_single_fasta(fasta_a, representative_a_name, sequence_a)
            write_single_fasta(fasta_b, representative_b_name, sequence_b)

            print(f"  Aligning class {class_a} vs {class_b} ...")
            paf_lines = run_minimap2(fasta_a, fasta_b)

            paf_out = output_directory / f"class_{class_a}_vs_{class_b}.paf"
            paf_out.write_text("\n".join(paf_lines) + "\n" if paf_lines else "")

            stats = parse_alignment_result(paf_lines, len(sequence_a), len(sequence_b))
            rows.append({
                "class_a": class_a,
                "class_b": class_b,
                "representative_a": representative_a_name,
                "representative_b": representative_b_name,
                "representative_a_length": len(sequence_a),
                "representative_b_length": len(sequence_b),
                **stats,
            })

    summary_path = output_directory / "rep_vs_rep_summary.tsv"
    rows_sorted = sorted(rows, key=lambda r: -r["fraction_of_shorter_representative"])
    write_tsv(summary_path, rows_sorted, SUMMARY_FIELDS)

    print(f"\nSummary ({len(rows)} pairs): {summary_path}")
    print()
    for row in rows_sorted:
        print(
            f"  {row['class_a']:>6} vs {row['class_b']:>6}  "
            f"fraction={row['fraction_of_shorter_representative']:.2f}  "
            f"verdict={row['verdict']}"
        )


# ---------------------------------------------------------------------------
# Command Line Interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run rep-vs-rep minimap2 for each class pair in large_pp_signatures.tsv."
    )
    parser.add_argument(
        "--signatures",
        type=Path,
        default=Path("reports/large_pp/large_pp_signatures.tsv"),
        help="Path to large_pp_signatures.tsv",
    )
    parser.add_argument(
        "--paf",
        type=Path,
        default=Path("input/minimap2_alignment_filtered.paf"),
        help="Path to filtered minimap2 PAF used to identify representatives",
    )
    parser.add_argument(
        "--representatives",
        type=Path,
        default=Path("input/mash_representatives.fa"),
        help="Path to mash_representatives.fa (default: input/mash_representatives.fa)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/rep_vs_rep"),
        help="Output directory for per-pair PAF files and summary TSV",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Script Entry Point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    run(args.signatures, args.paf, args.representatives, args.output)
