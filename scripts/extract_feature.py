#!/usr/bin/env python
"""Extract an annotated feature, or an arbitrary interval, to FASTA.

Pulls a named gene (or CDS, UTR, any feature type) out of a GenBank record and
writes it as FASTA, so a single gene can be used as the reference or the design
window for a later step. A minus-strand feature is reverse-complemented, so the
output always reads 5'->3' in the feature's own orientation. A compound
``join(...)`` location has its segments concatenated in positional order.

With ``--coords`` the same extraction works on a plain FASTA reference, for a
region that carries no annotation.

Examples
--------
    # the N gene of SARS-CoV-2 MT007544.1
    python scripts/extract_feature.py \
        --reference data/SARS-CoV-2/MT007544.1/MT007544.1.gbk \
        --gene N \
        --output results/SARS-CoV-2/MT007544.1_N.fasta

    # the same region by coordinate, from the FASTA
    python scripts/extract_feature.py \
        --reference data/SARS-CoV-2/MT007544.1/MT007544.1.fasta \
        --coords 28274-29533 \
        --output results/SARS-CoV-2/MT007544.1_N.fasta

Needs Biopython for a GenBank reference (e.g. the `pik3-guide` conda env).
"""

from __future__ import annotations

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from blast_utils import (  # noqa: E402
    BlastError,
    detect_reference_format,
    load_features,
    require_biopython,
    reverse_complement,
)

LINE_WIDTH = 60  # FASTA sequence wrapping, the usual NCBI width


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--reference", required=True, help="GenBank or FASTA reference")
    p.add_argument("--output", required=True, help="output FASTA")

    w = p.add_argument_group("what to extract (give one)")
    w.add_argument("--gene", help="gene name, as in the /gene qualifier (case-insensitive)")
    w.add_argument("--locus-tag", help="locus tag, as in the /locus_tag qualifier")
    w.add_argument("--product", help="substring of the /product qualifier (case-insensitive)")
    w.add_argument(
        "--coords",
        help="1-based inclusive interval START-END instead of a named feature; "
        "needs no annotation, so it works on a FASTA reference",
    )
    w.add_argument(
        "--strand",
        choices=("+", "-"),
        default="+",
        help="orientation for --coords (default: +); a named feature uses its own",
    )

    o = p.add_argument_group("options")
    o.add_argument(
        "--feature-type",
        default="gene",
        help="feature type to search, or 'any' for all types (default: gene). "
        "For a coding sequence use CDS",
    )
    o.add_argument("--record", help="restrict to this record ID, for a multi-record file")
    o.add_argument(
        "--all",
        action="store_true",
        help="write every match rather than failing when a name is ambiguous",
    )
    o.add_argument(
        "--name",
        help="sequence name for the FASTA header (default: <record>_<feature label>)",
    )
    return p.parse_args(argv)


def parse_coords(spec: str) -> tuple[int, int]:
    """Parse '28274-29533' (also '28274..29533', with optional commas) to ints."""
    cleaned = spec.replace(",", "").replace("..", "-").strip()
    m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", cleaned)
    if not m:
        raise ValueError(f"cannot read --coords {spec!r}; expected START-END, e.g. 28274-29533")
    start, end = int(m.group(1)), int(m.group(2))
    if start < 1:
        raise ValueError(f"--coords start must be 1 or greater, got {start}")
    if end < start:
        raise ValueError(f"--coords end ({end}) is before its start ({start})")
    return start, end


def load_sequences(reference: str, ref_format: str) -> dict[str, str]:
    """Read a GenBank or FASTA reference into ``{record_id: sequence}``."""
    SeqIO = require_biopython()
    sequences: dict[str, str] = {}
    for record in SeqIO.parse(reference, ref_format):
        try:
            seq = str(record.seq)
        except Exception as exc:
            raise BlastError(f"{reference} record {record.id} carries no sequence") from exc
        if not seq:
            raise BlastError(f"{reference} record {record.id} carries no sequence")
        sequences[record.id] = seq.upper()
    if not sequences:
        raise BlastError(f"no records parsed from {reference}")
    return sequences


def slice_parts(seq: str, parts, strand: str) -> str:
    """Concatenate 1-based inclusive *parts* of *seq*, honouring *strand*.

    Segments are joined in ascending positional order first, then the whole
    thing is reverse-complemented for a minus-strand feature -- which is the
    correct order of operations for a spliced location.
    """
    joined = "".join(seq[start - 1 : end] for start, end in parts)
    return reverse_complement(joined) if strand == "-" else joined


def write_fasta(path: str, entries: list[tuple[str, str]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as out:
        for header, seq in entries:
            out.write(f">{header}\n")
            for i in range(0, len(seq), LINE_WIDTH):
                out.write(seq[i : i + LINE_WIDTH] + "\n")


def select_features(features, args) -> list:
    """Features matching --gene/--locus-tag/--product, in positional order."""
    wanted_type = args.feature_type.strip().lower()
    found = []
    for record_id, feats in features.items():
        if args.record and record_id != args.record:
            continue
        for feat in feats:
            if wanted_type != "any" and feat.type.lower() != wanted_type:
                continue
            if args.gene and (feat.gene or "").lower() != args.gene.lower():
                continue
            if args.locus_tag and (feat.locus_tag or "") != args.locus_tag:
                continue
            if args.product and args.product.lower() not in (feat.product or "").lower():
                continue
            found.append(feat)
    return found


def main(argv=None) -> int:
    args = parse_args(argv)

    selectors = [args.gene, args.locus_tag, args.product, args.coords]
    if sum(s is not None for s in selectors) != 1:
        print(
            "error: give exactly one of --gene, --locus-tag, --product or --coords",
            file=sys.stderr,
        )
        return 2

    ref_format = detect_reference_format(args.reference)
    sequences = load_sequences(args.reference, ref_format)

    # --- an explicit interval, which needs no annotation -------------------
    if args.coords:
        start, end = parse_coords(args.coords)
        record_id = args.record or next(iter(sequences))
        if record_id not in sequences:
            print(
                f"error: record {record_id!r} not in {args.reference}; "
                f"available: {', '.join(sequences)}",
                file=sys.stderr,
            )
            return 2
        seq = sequences[record_id]
        if end > len(seq):
            print(
                f"error: --coords end {end} is past the end of {record_id} ({len(seq)} nt)",
                file=sys.stderr,
            )
            return 2
        sub = slice_parts(seq, [(start, end)], args.strand)
        name = args.name or f"{record_id}_{start}-{end}"
        header = f"{name} {record_id}:{start}-{end}({args.strand}) length={len(sub)}"
        write_fasta(args.output, [(header, sub)])
        print(f"wrote {len(sub)} nt to {args.output}  [{record_id}:{start}-{end} {args.strand}]")
        return 0

    # --- a named feature ---------------------------------------------------
    if ref_format != "genbank":
        print(
            f"error: {args.reference} is FASTA, which carries no annotation; "
            "use --coords, or pass a GenBank reference",
            file=sys.stderr,
        )
        return 2

    types = () if args.feature_type.strip().lower() == "any" else (args.feature_type,)
    features = load_features(args.reference, types)
    matches = select_features(features, args)

    if not matches:
        criterion = args.gene or args.locus_tag or args.product
        available = sorted(
            {f.label for feats in features.values() for f in feats if f.label}
        )
        print(
            f"error: no {args.feature_type} feature matching {criterion!r} in {args.reference}\n"
            f"available labels: {', '.join(available) if available else '(none)'}",
            file=sys.stderr,
        )
        return 2

    if len(matches) > 1 and not args.all:
        listed = "; ".join(
            f"{f.record_id} {f.type} {f.start}-{f.end}{f.strand}" for f in matches
        )
        print(
            f"error: {len(matches)} features match; pass --all to write them all, "
            f"or narrow with --feature-type/--record: {listed}",
            file=sys.stderr,
        )
        return 2

    entries = []
    for feat in matches:
        # sequence_parts, not parts: a frameshift join() overlaps by one base,
        # and that base must be kept twice to preserve the reading frame.
        sub = slice_parts(sequences[feat.record_id], feat.sequence_parts, feat.strand)
        name = args.name or f"{feat.record_id}_{feat.label}"
        header = (
            f"{name} {feat.record_id}:{feat.start}-{feat.end}({feat.strand}) "
            f"{feat.type} length={len(sub)}"
        )
        if feat.product:
            header += f" product={feat.product}"
        entries.append((header, sub))
        print(
            f"{feat.type} {feat.label}: {feat.record_id}:{feat.start}-{feat.end}"
            f"{feat.strand}, {len(sub)} nt"
            + (
                f" ({len(feat.sequence_parts)} segments)"
                if len(feat.sequence_parts) > 1
                else ""
            )
        )

    write_fasta(args.output, entries)
    print(f"wrote {len(entries)} sequence(s) to {args.output}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (BlastError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
