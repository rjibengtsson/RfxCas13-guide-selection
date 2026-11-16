#!/usr/bin/env python
"""Map sequences from a delimited table onto a reference genome with blastn.

Reads a CSV or TSV, blasts a user-chosen sequence column against a
user-supplied reference, and writes the original table back out with reference
coordinates and BLAST statistics appended. Every input row appears in the
output; rows without a hit are marked in the ``blast_status`` column. With a
GenBank reference, the overlapping genes are appended too.

``--revcomp-query`` reverse-complements each sequence before blasting, which is
how a table of Cas13 *spacers* is searched: a spacer is the reverse complement
of the RNA it targets, so the target-sense DNA is what the genome contains.

``--outfmt6`` additionally writes every HSP as a standard ``-outfmt 6`` table
with a header line, named by the input's own ID column.

Examples
--------
    # target-sense sequences, CSV in, MERS GenBank reference
    python scripts/blast_guides.py \
        --input data/MERS/mers_guides.csv \
        --seq-column "Target Sequence" \
        --reference data/MERS/NC_019843/NC_019843.3.fasta \
        --output results/MERS/mers_guides_NC_019843.csv \
        --all-hits results/MERS/mers_guides_NC_019843_all_hits.csv \
        --threads 4

    # PspCas13b spacers, TSV in, converted to target DNA first
    python scripts/blast_guides.py \
        --input data/pspcas13b_guides/panSarbeCov_guides.tsv \
        --seq-column Spacer --id-column Name --revcomp-query \
        --reference data/SARS-CoV-2/MT007544.1/MT007544.1.gbk \
        --output results/SARS-CoV-2/panSarbeCov_guides_MT007544.1_annotated.csv \
        --outfmt6 results/SARS-CoV-2/panSarbeCov_guides_MT007544.1_fmt6.tsv \
        --threads 4

Needs pandas and BLAST+ on PATH (e.g. the `pik3-guide` conda env); a GenBank
reference also needs Biopython.
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from blast_utils import (  # noqa: E402
    DEFAULT_EVALUE,
    DEFAULT_FEATURE_TYPES,
    DEFAULT_MAX_TARGET_SEQS,
    DEFAULT_TASK,
    DEFAULT_WORD_SIZE,
    BlastError,
    best_hit_per_query,
    count_hits_per_query,
    detect_reference_format,
    genbank_to_fasta,
    load_features,
    make_blast_db,
    overlapping_features,
    reverse_complement,
    run_blastn,
    write_outfmt6,
)
from table_io import read_table  # noqa: E402

# Columns appended to the input table, in output order.
HIT_COLUMNS = [
    "ref_accession",
    "ref_start",
    "ref_end",
    "strand",
    "pident",
    "aln_length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "qlen",
    "qcovhsp",
    "evalue",
    "bitscore",
]

# Columns appended on top of HIT_COLUMNS when annotation mode is on.
ANNOTATION_COLUMNS = [
    "n_features",
    "genes",
    "locus_tags",
    "products",
    "features",
    "feature_types",
    "feature_overlap_bp",
    "offset_in_feature",
]


def annotate(features, hit) -> dict:
    """Build the annotation columns for one hit from the GenBank feature table.

    ``features``/``feature_types``/``feature_overlap_bp``/``offset_in_feature``
    are parallel ';'-joined lists in positional order; ``genes``/``locus_tags``/
    ``products`` are de-duplicated, since a guide usually overlaps the gene, its
    CDS and a mat_peptide all carrying the same name.

    Two things to read carefully. Entries can repeat a name where the record
    holds genuinely distinct features for it -- a frameshifted ORF1ab yields
    both a join() CDS and a shorter ORF1a CDS, so ``CDS:orf1ab`` appears twice
    with different overlaps; ``products`` tells them apart. And
    ``offset_in_feature`` goes negative when the hit begins upstream of the
    feature, which together with a clipped ``feature_overlap_bp`` flags a guide
    straddling the boundary.
    """
    found = overlapping_features(features, hit.subject_id, hit.ref_start, hit.ref_end)

    def uniq(values):
        seen = {}
        for v in values:
            if v:
                seen[v] = None
        return ";".join(seen)

    return {
        "n_features": len(found),
        "genes": uniq(f.gene for f, _ in found),
        "locus_tags": uniq(f.locus_tag for f, _ in found),
        "products": uniq(f.product for f, _ in found),
        "features": ";".join(f"{f.type}:{f.label}" for f, _ in found),
        "feature_types": ";".join(f.type for f, _ in found),
        "feature_overlap_bp": ";".join(str(bp) for _, bp in found),
        "offset_in_feature": ";".join(
            str(f.offset_of(hit.ref_start, hit.ref_end)) for f, _ in found
        ),
    }


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True, help="input table (CSV or TSV)")
    p.add_argument(
        "--input-sep",
        default="auto",
        help="input delimiter: auto (from the file extension), tab, csv, ';', '|', "
        "whitespace, or a single character (default: auto)",
    )
    p.add_argument("--seq-column", required=True, help="column holding the sequences to blast")
    p.add_argument(
        "--revcomp-query",
        action="store_true",
        help="reverse-complement each sequence before blasting, to turn Cas13 "
        "spacers into the target-sense DNA the genome actually contains",
    )
    p.add_argument("--reference", required=True, help="reference to blast against (FASTA or GenBank)")
    p.add_argument("--output", required=True, help="output CSV (input columns + blast results)")
    p.add_argument("--all-hits", help="optional CSV of every HSP, not just the best per sequence")
    p.add_argument(
        "--outfmt6",
        help="optional headed blastn -outfmt 6 table of every HSP; gene columns "
        "are appended after the canonical fields when annotating",
    )
    p.add_argument("--id-column", help="column used to label queries in --all-hits (default: row number)")
    p.add_argument("--db-dir", help="where to cache the BLAST database (default: <output dir>/blastdb)")
    p.add_argument("--force-makedb", action="store_true", help="rebuild the BLAST database even if cached")

    a = p.add_argument_group("reference format and annotation")
    a.add_argument(
        "--ref-format",
        choices=("auto", "fasta", "genbank"),
        default="auto",
        help="reference format (default: auto-detect from the file's first line)",
    )
    a.add_argument(
        "--no-annotate",
        action="store_true",
        help="do not add feature columns even for a GenBank reference",
    )
    a.add_argument(
        "--feature-types",
        default=",".join(DEFAULT_FEATURE_TYPES),
        help="comma-separated GenBank feature types to report, or 'all' for every "
        f"type but source (default: {','.join(DEFAULT_FEATURE_TYPES)})",
    )

    b = p.add_argument_group("blastn parameters")
    b.add_argument("--task", default=DEFAULT_TASK, help=f"blastn task (default: {DEFAULT_TASK})")
    b.add_argument("--word-size", type=int, default=DEFAULT_WORD_SIZE)
    b.add_argument("--evalue", type=float, default=DEFAULT_EVALUE)
    b.add_argument("--max-target-seqs", type=int, default=DEFAULT_MAX_TARGET_SEQS)
    b.add_argument("--threads", type=int, default=1)

    f = p.add_argument_group("hit filters")
    f.add_argument("--min-pident", type=float, default=0.0, help="discard hits below this %% identity")
    f.add_argument("--min-coverage", type=float, default=0.0, help="discard hits below this %% query coverage")
    return p.parse_args(argv)


def clean_sequence(value) -> str:
    """Uppercase, drop whitespace, and map RNA U to T (blastn has no U)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return "".join(str(value).split()).upper().replace("U", "T")


def main(argv=None) -> int:
    args = parse_args(argv)

    df = read_table(args.input, args.input_sep)
    if args.seq_column not in df.columns:
        print(
            f"error: column {args.seq_column!r} not in {args.input}\n"
            f"available columns: {', '.join(map(repr, df.columns))}",
            file=sys.stderr,
        )
        return 2
    if args.id_column and args.id_column not in df.columns:
        print(f"error: --id-column {args.id_column!r} not in {args.input}", file=sys.stderr)
        return 2

    # Internal query IDs are positional, so duplicate or whitespace-bearing
    # values in the user's own ID column can never collide or break the FASTA.
    sequences = df[args.seq_column].map(clean_sequence)

    # Spacers are blasted as their reverse complement, and the derived sequence
    # is carried into the output so the coordinates can be traced back to it.
    if args.revcomp_query:
        try:
            sequences = sequences.map(lambda s: reverse_complement(s) if s else "")
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        df = df.assign(target_dna=sequences.values)
        print(f"--revcomp-query: blasting the reverse complement of {args.seq_column!r} "
              "(written to the target_dna column)")

    records = [(f"q{i}", seq) for i, seq in enumerate(sequences) if seq]
    n_empty = len(df) - len(records)
    print(f"read {len(df)} rows from {args.input}; {len(records)} sequences to blast"
          + (f" ({n_empty} empty)" if n_empty else ""))

    db_dir = args.db_dir or os.path.join(os.path.dirname(os.path.abspath(args.output)), "blastdb")

    ref_format = args.ref_format if args.ref_format != "auto" else detect_reference_format(args.reference)
    reference_fasta = args.reference
    features = None
    if ref_format == "genbank":
        # BLAST+ cannot read GenBank, so convert; annotation mode follows from it.
        reference_fasta = genbank_to_fasta(
            args.reference,
            os.path.join(db_dir, os.path.basename(args.reference) + ".fasta"),
            force=args.force_makedb,
        )
        print(f"reference format: genbank -> converted to {reference_fasta}")
        if not args.no_annotate:
            types = () if args.feature_types.strip().lower() == "all" else tuple(
                t.strip() for t in args.feature_types.split(",") if t.strip()
            )
            features = load_features(args.reference, types)
            n_feat = sum(len(v) for v in features.values())
            print(f"annotation mode: on ({n_feat} features across {len(features)} record(s))")
            if n_feat == 0:
                print(
                    f"warning: no features of type(s) {args.feature_types} found; "
                    "annotation columns will be empty",
                    file=sys.stderr,
                )
    else:
        print("reference format: fasta (no annotation available)")

    db = make_blast_db(reference_fasta, db_dir, force=args.force_makedb)
    print(f"blast database: {db}")

    hits = run_blastn(
        records,
        db,
        task=args.task,
        word_size=args.word_size,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        threads=args.threads,
    )
    print(f"blastn returned {len(hits)} HSPs")

    if args.min_pident or args.min_coverage:
        kept = [h for h in hits if h.pident >= args.min_pident and h.qcovhsp >= args.min_coverage]
        print(f"{len(kept)} HSPs pass the filters (pident >= {args.min_pident}, "
              f"coverage >= {args.min_coverage})")
        hits = kept

    best = best_hit_per_query(hits)
    counts = count_hits_per_query(hits)

    rows = []
    for i, seq in enumerate(sequences):
        qid = f"q{i}"
        hit = best.get(qid)
        if not seq:
            status = "empty_sequence"
        elif hit is None:
            status = "no_hit"
        else:
            status = "hit"
        row = {"n_hits": counts.get(qid, 0), "blast_status": status}
        if hit is not None:
            rec = hit.as_record()
            row.update({c: rec.get(c) for c in HIT_COLUMNS if c != "ref_accession"})
            row["ref_accession"] = hit.subject_id
            if features is not None:
                row.update(annotate(features, hit))
        rows.append(row)

    columns = HIT_COLUMNS + ["n_hits", "blast_status"]
    if features is not None:
        columns += ANNOTATION_COLUMNS
    out = pd.concat(
        [df.reset_index(drop=True), pd.DataFrame(rows, columns=columns)],
        axis=1,
    )
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    out.to_csv(args.output, index=False)

    summary = out["blast_status"].value_counts().to_dict()
    print(f"wrote {len(out)} rows to {args.output}  {summary}")

    # Row labels, shared by --all-hits and --outfmt6, so both name queries the
    # same way. Falls back to the row number when no ID column is given.
    labels = (
        df[args.id_column].astype(str).tolist() if args.id_column else [str(i) for i in range(len(df))]
    )

    if hits:
        strands = {}
        for hit in hits:
            strands[hit.strand] = strands.get(hit.strand, 0) + 1
        print(f"strand of reported HSPs: {strands}")

    if args.outfmt6:
        # Keyed by id(hit) because one query can yield several HSPs, each with
        # its own overlapping features.
        annotations = (
            {id(h): annotate(features, h) for h in hits} if features is not None else None
        )
        n = write_outfmt6(
            hits,
            args.outfmt6,
            labels={f"q{i}": labels[i] for i in range(len(labels))},
            annotations=annotations,
            annotation_columns=ANNOTATION_COLUMNS if features is not None else (),
        )
        print(f"wrote {n} HSPs to {args.outfmt6} (headed outfmt 6"
              + (" + gene columns)" if features is not None else ")"))

    if args.all_hits:
        all_rows = []
        for hit in hits:
            rec = hit.as_record()
            pos = int(hit.query_id[1:])
            all_rows.append(
                {
                    "input_row": pos,
                    "query_label": labels[pos],
                    "sequence": sequences.iloc[pos],
                    "ref_accession": hit.subject_id,
                    **{c: rec.get(c) for c in HIT_COLUMNS if c != "ref_accession"},
                    **(annotate(features, hit) if features is not None else {}),
                }
            )
        pd.DataFrame(all_rows).to_csv(args.all_hits, index=False)
        print(f"wrote {len(all_rows)} HSPs to {args.all_hits}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (BlastError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
