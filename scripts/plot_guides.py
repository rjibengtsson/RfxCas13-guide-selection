#!/usr/bin/env python
"""Plot RfxCas13d guides along a sequence, marking the PspCas13b sites.

The normal input is the located table from ``scripts/locate_guides.py``, which
already holds each guide's position, id and overlap class. That script is the
single place position is computed, so the figure and the table cannot disagree.

``--predictions`` remains as a shortcut for plotting a raw TIGER on-target file
without writing an intermediate first: positions are then derived here, by the
same function ``locate_guides.py`` calls. Prefer ``--located`` for anything
whose numbers you intend to quote.

Each guide is coloured by whether it falls fully inside a PspCas13b footprint,
partly over one, or clear of them all. Nothing is filtered and no score
threshold is applied: every guide is drawn.

Examples
--------
    # from the located table (preferred)
    python scripts/plot_guides.py \
        --located results/SARS-CoV-2/rfx_guides_MT007544.1_N_located.csv \
        --sequence results/SARS-CoV-2/MT007544.1_N.fasta \
        --output figures/SARS-CoV-2/MT007544.1_N_guides.html

    # straight from the predictions, deriving positions on the fly
    python scripts/plot_guides.py \
        --predictions results/SARS-CoV-2/on_target.csv \
        --sequence results/SARS-CoV-2/MT007544.1_N.fasta \
        --psp-hits results/SARS-CoV-2/panSarbeCov_guides_MT007544.1_fmt6.tsv \
        --output figures/SARS-CoV-2/MT007544.1_N_guides.html

Needs pandas and plotly (e.g. the `pik3-guide` conda env).
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from table_io import read_table  # noqa: E402
from visualise_scatterlines import (  # noqa: E402
    PREDICTION_COLUMNS,
    build_guides,
    build_psp_guides,
    check_against_sequence,
    drop_offsequence,
    frame_to_guides,
    plot_guides,
    psps_from_frame,
    read_fasta,
    region_start_from_header,
)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--located",
        help="located guides from scripts/locate_guides.py: positions, ids and "
        "overlap classes are read from it rather than recomputed",
    )
    source.add_argument(
        "--predictions",
        help="raw on-target predictions (CSV or TSV); positions are derived here",
    )

    p.add_argument("--sequence", required=True, help="FASTA of the sequence to plot along")
    p.add_argument("--output", required=True, help="output HTML")

    p.add_argument(
        "--psp-hits",
        help="blastn table of PspCas13b guides to mark (headed outfmt 6, as "
        "written by scripts/blast_guides.py --outfmt6). With --located the "
        "footprints are taken from the table itself, so this is only needed to "
        "mark a footprint no guide lies near",
    )
    p.add_argument(
        "--region-start",
        type=int,
        help="1-based position where the plotted sequence begins in the "
        "reference the hits were called against (default: read from the FASTA "
        "header, else 1)",
    )
    p.add_argument(
        "--reference",
        help="keep only hits whose subject is this accession (default: keep all)",
    )
    p.add_argument("--title", help="plot title (default: generated from the sequence name)")
    p.add_argument(
        "--cdn-plotlyjs",
        action="store_true",
        help="link plotly.js from a CDN instead of inlining it: a ~3 MB smaller "
        "file that needs a network connection to open",
    )
    return p.parse_args(argv)


def _resolve_region_start(args: argparse.Namespace, header: str) -> int:
    """The offset putting reference-coordinate hits onto the plotted sequence.

    Getting it wrong silently slides every footprint along the sequence, and
    the plot would still look plausible, so the source is always printed.
    """
    if args.region_start is not None:
        region_start, source = args.region_start, "--region-start"
    else:
        from_header = region_start_from_header(header)
        region_start = from_header if from_header is not None else 1
        source = "FASTA header" if from_header is not None else "default"
    print(f"region start: {region_start:,} (from {source})")
    return region_start


def _finish(args, guides, psps, sequence) -> int:
    """Report the class counts and write the figure. Shared by both inputs."""
    for psp in sorted(psps, key=lambda p: p.start):
        print(f"  {psp.name}: {psp.start:,}-{psp.end:,} ({psp.length} nt)")

    by_class: dict[str, int] = {}
    for guide in guides:
        by_class[guide.classification] = by_class.get(guide.classification, 0) + 1
    print(
        "guides by overlap: "
        + ", ".join(f"{k} {v:,}" for k, v in sorted(by_class.items()))
    )

    written = plot_guides(
        guides,
        psps,
        length=len(sequence),
        output=args.output,
        title=args.title,
        inline_plotlyjs=not args.cdn_plotlyjs,
    )
    if written is None:
        return 1
    print(f"wrote {written}")
    return 0


def main(argv=None) -> int:
    args = parse_args(argv)

    header, sequence = read_fasta(args.sequence)
    name = header.split()[0]
    print(f"sequence: {name}, {len(sequence):,} nt")

    # --- the located table: positions already computed, just read them -----
    if args.located:
        located = read_table(args.located)
        print(f"located guides: {len(located):,} rows from {args.located}")

        psps = []
        if args.psp_hits:
            region_start = _resolve_region_start(args, header)
            psps = build_psp_guides(
                read_table(args.psp_hits),
                region_start=region_start,
                reference=args.reference,
            )
            psps = drop_offsequence(psps, len(sequence))
            print(f"PspCas13b footprints: {len(psps)} from {args.psp_hits}")
        else:
            psps = psps_from_frame(located)
            print(f"PspCas13b footprints: {len(psps)} named in the located table")

        guides = frame_to_guides(located, psps)
        # The table's coordinates were computed against some sequence; this is
        # what proves it was this one.
        check_against_sequence(guides, sequence)
        print("checked: every guide matches the sequence at its stated position")
        return _finish(args, guides, psps, sequence)

    # --- raw predictions: derive the positions here -------------------------
    predictions = read_table(args.predictions)
    missing = [c for c in PREDICTION_COLUMNS if c not in predictions.columns]
    if missing:
        print(
            f"error: {args.predictions} is missing column(s) "
            f"{', '.join(map(repr, missing))}\n"
            f"available columns: {', '.join(map(repr, predictions.columns))}",
            file=sys.stderr,
        )
        return 2
    print(f"predictions: {len(predictions):,} rows from {args.predictions}")

    psps = []
    if args.psp_hits:
        region_start = _resolve_region_start(args, header)
        hits = read_table(args.psp_hits)
        psps = build_psp_guides(hits, region_start=region_start, reference=args.reference)
        psps = drop_offsequence(psps, len(sequence))
        print(f"PspCas13b footprints: {len(psps)} from {args.psp_hits}")

    guides = build_guides(predictions, sequence, psps, transcript=name)
    return _finish(args, guides, psps, sequence)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
