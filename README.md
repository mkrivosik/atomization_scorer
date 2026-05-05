# Atomization Scorer

Atomization Scorer is a tool for evaluating the quality of genome atomization produced
by [GEESE](https://gitlab.ub.uni-bielefeld.de/gi/geese). DNA sequences of related
organisms evolve through local mutations and large-scale events such as duplications and
rearrangements, producing mosaic genomes composed of homologous regions interspersed with
regions of low or no similarity. The IMP algorithm [1] segments DNA into non-overlapping
atomic segments (atoms) and partitions them into classes such that atoms within a class
are mutually similar and atoms from different classes share no significant similarity.
These atoms serve as fine-scale genomic markers for evolutionary analyses. Atomization
Scorer evaluates how accurately a predicted atomization reflects the underlying genomic
structure by constructing a gold-standard atomization from the input sequences via
representative-based alignments and comparing the prediction against it.

The overall quality score $S \in [0, 1]$ is the weighted geometric mean of two metrics:
an alignment score $A$ measuring structural accuracy of predicted atoms against the gold
standard, and a coverage score $C$ measuring what fraction of total genome length is
covered by predicted atoms.

## Installation

```bash
git clone <repository-url>
cd atomization_scorer
pip install -e .
```

The following external tools must be available on `PATH`: `minimap2` [3] (always),
`mash` [2] (for the default and recommended representative selection mode), and `docker` (optional; for Dotter [4]
dot-plot visualization, which also requires an X11 display).

## Usage

```
atomization_scorer GENOMES_FASTA ATOMIZATION_GEESE OUTPUT_DIRECTORY [OPTIONS]
```

The three positional arguments are required. `GENOMES_FASTA` must have extension `.fa`
or `.fasta`; `ATOMIZATION_GEESE` must have extension `.geese`. The output directory is
created automatically if it does not exist. Run `atomization_scorer --help` for the full
option reference.

Key options are organized into four groups:

**Scoring.** `--level {interval,base}` (default: `interval`) selects scoring granularity.
Interval-level scoring treats each predicted atom as a unit and counts it as a true
positive if its Jaccard overlap with any true atom meets `--min-overlap-ratio`
($\theta = 0.8$ by default). Base-level scoring classifies each nucleotide independently
as $\mathrm{TP}$, $\mathrm{FP}$, or $\mathrm{FN}$. `--per-class` computes $A$, $C$,
and $S$ independently for each atomization class. `--alignment-weight` ($w_a = 0.7$ by
default) and `--coverage-weight` ($w_c = 0.3$ by default) must satisfy $w_a + w_c = 1.0$.

**True alignment pipeline.** `--representative-mode {mash,first}` (default: `mash`)
selects the representative per class: `mash` picks the medoid by Mash [2] distance;
`first` picks the first atom in the GEESE file. `--min-similarity` (default: `0.95`) and
`--min-alignment-length` (default: `500`) filter PAF rows before overlap resolution.

**Minimap2.** `--minimap2-preset` (default: `asm20`), `--minimap2-secondary-ratio`
(default: `0.1`), and `--minimap2-no-cigar` control the minimap2 [3] call directly.

**Overlap diagnostics.** `--overlap-diagnostics` enables inspection of overlapping
alignments in the filtered PAF before overlap resolution, writing tabular reports, a JSON
overlap record, and per-anchor FASTA inputs for Dotter [4] dot-plot visualization.
`--skip-dotter` suppresses only the Dotter launch; all diagnostic outputs are still written.

Logging is written to stderr at INFO level by default; `--verbose` enables DEBUG and
`--quiet` suppresses INFO messages.

## Input Files

**FASTA.** Multi-record FASTA (`.fa` or `.fasta`) with one record per sequence. Only the first whitespace-delimited token of each header is used as the sequence name.

**GEESE.** Tab-separated file with required columns `#name` (sequence name), `class`
(integer class label), `start` (0-based, inclusive), and `end` (exclusive). Coordinates
follow half-open interval semantics $[\text{start}, \text{end})$. Additional columns are
preserved. Predicted atoms within the same sequence are non-overlapping by the IMP atom definition.

## Output

```
OUTPUT_DIRECTORY/
├── run_parameters.json                  # all parameters, resolved input paths, timestamp
├── scores.tsv                           # alignment, coverage, overall scores and weights (without --per-class)
├── scores_per_class.tsv                 # per-class scores (with --per-class, replaces scores.tsv)
├── true_atomization/
│   ├── {mode}_representatives.fa        # one representative per predicted class
│   ├── minimap2_alignments.paf          # raw minimap2 output
│   ├── minimap2_alignment_filtered.paf  # after similarity and length filtering
│   ├── minimap2_alignment_resolved.paf  # after overlap resolution
│   └── true_atomization.geese           # gold-standard atomization (validated)
├── metrics/                             # TP/FP/FN status and metric TSVs (filenames vary with --level and --per-class)
├── visualization/
│   ├── atomization_visualization.html   # index linking all per-sequence plots
│   └── genomes_visualization/           # one interactive HTML plot per sequence
└── overlap_diagnostics/                 # only with --overlap-diagnostics
    ├── summary.tsv
    ├── anchor_genome_summary.tsv
    ├── overlaps.json
    └── anchors/<anchor_id>/
        ├── X.fasta
        ├── Y.fasta
        ├── pairs.tsv
        └── dotter.pdf
```

## Method

**Gold-standard construction.** For each predicted atomization class, one representative
sequence is selected. With the default `mash` mode, the medoid is chosen: the atom
minimizing total Mash [2] distance to all other atoms in the class. Mash uses the MinHash
technique to reduce sequences to compact k-mer sketches and estimates genomic distance
from their Jaccard similarity without full pairwise alignment. All sequences are
then aligned to the representative sequences using minimap2 [3] with preset `asm20`,
secondary ratio `0.1`, and CIGAR output enabled. Alignments with similarity below `0.95` or aligned length below `500` bp are discarded. Non-primary alignments are then discarded, and remaining primary alignments that overlap on query-sequence coordinates are resolved into a non-overlapping partition. Each accepted alignment
inherits the class label of the representative it was aligned to, so every position in
the gold-standard atomization is assigned exactly one class. The resolved alignments are
converted to GEESE format and validated for non-overlap to form the gold-standard
atomization.

**Coverage score.** Let $\mathcal{P}$ be the set of predicted atoms and $L = \sum_j |\mathcal{S}_j|$ the total length of all input sequences. The coverage score $C \in [0, 1]$ is:

$$C = \frac{\sum_{p \in \mathcal{P}}\,(\text{end}_p - \text{start}_p)}{L}$$

Atoms are non-overlapping by definition. For class $k$, let $\mathcal{P}_k = \{p \in \mathcal{P} \mid \text{class}(p) = k\}$ denote the predicted atoms of class $k$:

$$C_k = \frac{\sum_{p \in \mathcal{P}_k}(\text{end}_p - \text{start}_p)}{L}, \qquad \sum_k C_k = C$$

**Alignment score.** Let $\mathcal{T}$ be the set of true atoms. At interval level, predicted atom $p \in \mathcal{P}$ is a true positive ($\mathrm{TP}$) if
there exists an as-yet-unmatched $t \in \mathcal{T}$ on the same sequence with
$\text{class}(t) = \text{class}(p)$ such that

$$\frac{|\,[p_s,\,p_e) \cap [t_s,\,t_e)\,|}{|\,[p_s,\,p_e) \cup [t_s,\,t_e)\,|} \geq \theta$$

where $p_s, p_e$ and $t_s, t_e$ are the start and end coordinates of $p$ and $t$
respectively, and $\theta \in [0,1]$ is the overlap threshold (`--min-overlap-ratio`,
default $\theta = 0.8$). Each true atom is consumed by at most one predicted atom.
Predicted atoms with no qualifying match are false positives ($\mathrm{FP}$); unmatched
true atoms are false negatives ($\mathrm{FN}$). The alignment score $A$ is:

$$\text{Precision} = \frac{|\mathrm{TP}|}{|\mathrm{TP}| + |\mathrm{FP}|}, \qquad \text{Recall} = \frac{|\mathrm{TP}|}{|\mathrm{TP}| + |\mathrm{FN}|}$$

$$A = F_1 = \frac{2 \cdot \text{Precision} \cdot \text{Recall}}{\text{Precision} + \text{Recall}}$$

At base level, each nucleotide position is classified independently as $\mathrm{TP}$, $\mathrm{FP}$, or $\mathrm{FN}$ and $F_1$ is computed over nucleotide counts.

**Overall score.** Let $w_a, w_c \in (0,1)$ with $w_a + w_c = 1$ (defaults: $w_a = 0.7$, $w_c = 0.3$). The overall score is:

$$S = A^{w_a} \cdot C^{w_c}$$

The weighted geometric mean ensures $S = 0$ if either $A = 0$ or $C = 0$, and $S = 1$ only if $A = C = 1$. Since $A, C \in [0,1]$, the score $S \in [0,1]$ by construction. Per-class scores use $A_k$, the $F_1$ score restricted to atoms of class $k$, and $C_k$ as defined above:

$$S_k = A_k^{w_a} \cdot C_k^{w_c}$$

## References

[1] Višňovská M, Vinař T, Brejová B. DNA Sequence Segmentation Based on Local
Similarity. *ITAT 2013 Proceedings*, CEUR Workshop Proceedings Vol. 1003, pp. 36–43.
http://ceur-ws.org/Vol-1003

[2] Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM.
Mash: fast genome and metagenome distance estimation using MinHash.
*Genome Biology* (2016) **17**:132.
https://doi.org/10.1186/s13059-016-0997-x

[3] Li H. Minimap2: pairwise alignment for nucleotide sequences.
*Bioinformatics* (2018) **34**(18):3094–3100.
https://doi.org/10.1093/bioinformatics/bty191

[4] Sonnhammer EL, Durbin R. A dot-matrix program with dynamic threshold control suited
for genomic DNA and protein sequence analysis. *Gene* (1995) **167**(1-2):GC1-10.
https://sonnhammer.sbc.su.se/Dotter.html

**Software:**
- GEESE: https://gitlab.ub.uni-bielefeld.de/gi/geese
- Mash: https://github.com/marbl/Mash
- Minimap2: https://minimap2.com/
- Dotter: https://sonnhammer.sbc.su.se/Dotter.html

## License

GPL-3.0.

## Contact

Matej Krivošík  
krivosik7@uniba.sk  
Faculty of Mathematics, Physics, and Informatics, Comenius University, Bratislava

Supervisor: doc. Mgr. Tomáš Vinař, PhD.  
tomas.vinar@fmph.uniba.sk  
Faculty of Mathematics, Physics, and Informatics, Comenius University, Bratislava
