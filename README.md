# RfxCas13 guide selection — BLAST mapping and annotation

Maps guide/target sequences from a CSV onto a reference genome with `blastn`, and
writes the original table back out with reference coordinates, BLAST statistics
and (for a GenBank reference) overlapping gene features appended.

- `src/blast_utils.py` — reusable library: database building, blastn, GenBank features
- `scripts/blast_guides.py` — the command-line entry point
- `tests/` — pytest suite (80 tests; integration tests skip if BLAST+ or the data are absent)

## Requirements

- BLAST+ (`blastn`, `makeblastdb`) on `PATH`
- Python with `pandas`, and `biopython` only if you use a GenBank reference
- `pytest` to run the tests

In this workspace the `pik3-guide` conda env has all of these:

```bash
/home/unimelb.edu.au/rbengtsson/miniconda3/envs/pik3-guide/bin/python -m pytest tests/ -q
```

## Usage

```bash
python scripts/blast_guides.py \
    --input data/mers_guides.csv \
    --seq-column "Target Sequence" \
    --reference data/NC_019843/NC_019843.3.gbk \
    --output results/mers_guides_NC_019843_annotated.csv \
    --all-hits results/mers_guides_NC_019843_annotated_all_hits.csv \
    --id-column Index \
    --threads 4
```

The reference may be FASTA or GenBank; the format is detected from the file's
first line. A GenBank reference turns **annotation mode** on automatically: the
record's sequence is converted to a cached FASTA for BLAST (BLAST+ cannot read
GenBank), and its feature table is used to add the annotation columns below.

### Options

| Option | Meaning |
| --- | --- |
| `--input` | Input CSV. |
| `--seq-column` | **Required.** Column holding the sequences to BLAST. No default, so the choice is always explicit; an unknown name exits 2 and lists the available columns. |
| `--reference` | Reference FASTA or GenBank file. |
| `--output` | Output CSV: the input columns plus the fields documented below. |
| `--all-hits` | Optional second CSV with *every* HSP rather than one row per input sequence. |
| `--id-column` | Column used to label queries in `--all-hits` (default: row number). |
| `--db-dir` | Where to cache the BLAST database (default: `<output dir>/blastdb`). Must not contain whitespace. |
| `--force-makedb` | Rebuild the BLAST database and the GenBank→FASTA conversion even if cached. |
| `--ref-format` | `auto` (default), `fasta` or `genbank`, overriding detection. |
| `--no-annotate` | Skip the annotation columns even for a GenBank reference. |
| `--feature-types` | Comma-separated GenBank feature types to report, or `all` for every type but `source`. Default `gene,CDS,5'UTR,3'UTR`. |
| `--task` | blastn task. Default `blastn-short`, appropriate for the 23–30 nt guides here. |
| `--word-size`, `--evalue`, `--max-target-seqs`, `--threads` | Passed through to `blastn`. Defaults 11, 1e-3, 50, 1. |
| `--min-pident`, `--min-coverage` | Discard HSPs below this % identity / % query coverage before the best hit is chosen. |

## Output fields

Every input row appears in the output exactly once, in the input order, with all
original columns unchanged. Rows whose sequence produced no alignment are kept
and flagged in `blast_status` rather than dropped, so the output can be diffed
row-for-row against the input.

Where a sequence has several alignments, these columns describe the **best** one:
highest `bitscore`, then lowest `evalue`, then longest `aln_length`. Use
`--all-hits` to see the rest.

### Alignment and position

| Field | Meaning |
| --- | --- |
| `ref_accession` | Reference sequence the hit falls on — the FASTA record ID, or the GenBank accession.version. For a multi-record reference this says *which* record. |
| `ref_start` | Start of the alignment on the reference, 1-based inclusive. Always ≤ `ref_end`, on both strands. |
| `ref_end` | End of the alignment on the reference, 1-based inclusive. |
| `strand` | `+` or `-`: the orientation of the hit on the reference. A guide and its target sequence map to the *same* interval on *opposite* strands, so the interval is comparable whichever column you BLAST. |
| `qstart` | First base of the **query** that aligned, 1-based. Greater than 1 means the start of your sequence did not align. |
| `qend` | Last base of the query that aligned. |
| `qlen` | Full length of the query sequence. In `data/mers_guides.csv` this is 23 for RfxCas13d and 30 for PspCas13b. |
| `aln_length` | Length of the alignment, including any gaps. |

Note `ref_start`/`ref_end` are normalised low→high, which differs from raw
`blastn` output: in `-outfmt 6` a minus-strand hit reports `sstart > send`.
Strand is kept in its own column instead.

### Alignment quality

| Field | Meaning |
| --- | --- |
| `pident` | % identity over the aligned region. |
| `qcovhsp` | % of the query covered by this alignment. `pident` 100 with `qcovhsp` 65 means a perfect match of only part of the guide — check both. |
| `mismatch` | Mismatched positions in the alignment. |
| `gapopen` | Number of gap openings. |
| `evalue` | BLAST expect value. Short queries give modest e-values even for perfect matches (~1e-13 for a 30 nt exact match here), so judge quality on `pident` and `qcovhsp`. |
| `bitscore` | BLAST bit score; the primary key for choosing the best hit. |

### Row status

| Field | Meaning |
| --- | --- |
| `n_hits` | Number of HSPs for this sequence that survived `--min-pident`/`--min-coverage`. Greater than 1 means the sequence aligns more than once — within-reference off-target potential, or simply overlapping records. `0` for `no_hit` and `empty_sequence`. |
| `blast_status` | `hit` — an alignment was found and the columns above are populated. `no_hit` — the sequence was BLASTed but nothing passed the e-value threshold and filters. `empty_sequence` — the cell was blank or missing, so nothing was BLASTed. |

### Annotation (GenBank reference only)

These appear only when annotation mode is on. They are derived by interval
overlap between `ref_start`–`ref_end` and the GenBank feature table, so a hit is
annotated with **every** feature it touches, not one winner.

| Field | Meaning |
| --- | --- |
| `n_features` | How many features the hit overlaps by at least 1 bp. |
| `genes` | `/gene` qualifiers of those features, de-duplicated. Empty when no overlapping feature carries one — notably the `5'UTR` and `3'UTR` features in `NC_019843.3`, so read `features` rather than `genes` when UTRs matter. |
| `locus_tags` | `/locus_tag` qualifiers, de-duplicated. |
| `products` | `/product` qualifiers, de-duplicated. This is what distinguishes two same-named features, e.g. `1AB polyprotein` vs `1A polyprotein`. |
| `features` | One `type:name` entry per overlapping feature, e.g. `CDS:S;gene:S`. `name` is the `/gene`, else `/locus_tag`, else `/product`, else the feature type. |
| `feature_types` | The feature type of each entry, e.g. `gene;CDS;3'UTR`. |
| `feature_overlap_bp` | Bases shared between the hit and each feature. Less than `aln_length` means the hit extends past that feature's boundary. |
| `offset_in_feature` | Where the hit starts **inside** each feature, 1-based: `ref_start − feature_start + 1` for a plus-strand feature, `feature_end − ref_end + 1` for a minus-strand one, so it is always measured from the feature's own 5' end. `1` is the feature's first base. **Negative** means the hit begins upstream of the feature — read it with `feature_overlap_bp` to see how much actually lands inside. |

`features`, `feature_types`, `feature_overlap_bp` and `offset_in_feature` are
**parallel `;`-joined lists** in positional order: the *n*th entry of each refers
to the same feature. `genes`, `locus_tags` and `products` are de-duplicated and
therefore **not** parallel to them.

Two things to read carefully:

- A name can repeat where the record holds genuinely distinct features for it.
  `NC_019843.3` has both the frameshifted ORF1ab CDS, `join(279..13433,13433..21514)`,
  and a shorter ORF1a CDS, `279..13454`, so `CDS:orf1ab` appears twice with
  different overlaps; `products` tells them apart.
- `offset_in_feature` is a nucleotide offset along the feature's span, **not** a
  protein position: it does not skip introns and does not apply `codon_start`.
  Do not convert it to a codon number for `orf1ab`, whose offset runs straight
  through the frameshift point.

### `--all-hits` columns

One row per HSP. It carries the same alignment, quality and annotation columns as
above (no `n_hits` or `blast_status`, which are per-input-row), preceded by:

| Field | Meaning |
| --- | --- |
| `input_row` | 0-based position of the source row in the input CSV, for joining back. |
| `query_label` | Value of `--id-column` for that row, else the row number. |
| `sequence` | The cleaned sequence that was BLASTed (uppercased, whitespace removed, `U`→`T`). |

## Notes on behaviour

- **Sequence cleaning.** Values are uppercased, internal whitespace is removed,
  and `U` is mapped to `T`, so RNA-style input works. Blank or missing cells are
  reported as `empty_sequence`.
- **One BLAST process.** All sequences go into a single `blastn` call rather than
  one call per sequence, which is why 18,676 guides take a few seconds. The
  library still exposes `blast_one()` for the single-sequence case.
- **Query IDs are positional internally** (`q0`, `q1`, …), so duplicate or
  whitespace-bearing values in your own ID column cannot collide or corrupt the
  query FASTA.
- **Caching.** The BLAST database and the GenBank→FASTA conversion are rebuilt
  only when missing or older than the source; `--force-makedb` forces both.
- **Whitespace in paths.** BLAST+ parses `-in`/`-db` as space-separated file
  lists, so a reference path containing whitespace is passed via a whitespace-free
  symlink. A `--db-dir` containing whitespace is refused outright.
