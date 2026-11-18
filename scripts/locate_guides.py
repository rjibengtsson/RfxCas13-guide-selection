#!/usr/bin/env python
"""Locate predicted RfxCas13d guides on a sequence and write them out with ids.

The TIGER on-target output states no coordinates, only the target each guide
was scored against. This script recovers each guide's position by locating that
target in the reference FASTA, gives the guide an id carrying its span, and
writes one row per guide. The result is the intermediate table the rest of the
pipeline joins on, so position is computed once rather than re-derived by every
step that needs it.

Each row carries **two** coordinate systems: ``start``/``end`` on the plotted
sequence (what the id uses) and ``ref_start``/``ref_end`` on the reference the
region was cut from, so a guide can be matched back to the blast output without
redoing the offset arithmetic.

When PspCas13b hits are supplied, each guide is also classified by how it sits
against them -- ``full`` (entirely inside a footprint), ``partial`` (shares a
base but hangs over an edge) or ``none`` -- with the nearest footprint and the
signed gap to it. Nothing is filtered: every predicted guide is written.

A ``.params.txt`` sidecar is written next to the output recording the inputs and
the region offset, so a run can be reproduced.

Example
-------
    python scripts/locate_guides.py \
        --predictions results/SARS-CoV-2/on_target.csv \
        --sequence results/SARS-CoV-2/MT007544.1_N.fasta \
        --psp-hits results/SARS-CoV-2/panSarbeCov_guides_MT007544.1_fmt6.tsv \
        --output results/SARS-CoV-2/rfx_guides_MT007544.1_N_located.csv

Needs pandas (e.g. the `pik3-guide` conda env).
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from table_io import read_table  # noqa: E402
from visualise_scatterlines import (  # noqa: E402
    PREDICTION_COLUMNS,
    build_guides,
    build_psp_guides,
    drop_offsequence,
    filter_psps,
    guides_to_frame,
    read_fasta,
    region_start_from_header,
)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--predictions", required=True, help="on-target predictions (CSV or TSV)")
    p.add_argument("--sequence", required=True, help="FASTA of the sequence to locate against")
    p.add_argument("--output", required=True, help="output CSV of located guides")

    p.add_argument(
        "--psp-hits",
        help="blastn table of PspCas13b guides (headed outfmt 6). Omit to write "
        "positions without the overlap classification",
    )
    p.add_argument(
        "--region-start",
        type=int,
        help="1-based position where the sequence begins in the reference the "
        "hits were called against (default: read from the FASTA header, else 1)",
    )
    p.add_argument(
        "--reference",
        help="keep only hits whose subject is this accession, and record it as "
        "the reference (default: taken from the hits)",
    )
    p.add_argument(
        "--psp-name",
        action="append",
        metavar="NAME",
        help="restrict the footprints to this guide, repeatable. Guides are "
        "then classified against only the named footprints, so a guide "
        "overlapping an excluded one is reported as no overlap (default: use "
        "every footprint in the hits)",
    )
    p.add_argument(
        "--no-params",
        action="store_true",
        help="skip the .params.txt sidecar",
    )
    return p.parse_args(argv)


def write_params(path: str, args: argparse.Namespace, extra: dict) -> str:
    """Record the run's inputs beside its output, so it can be reproduced."""
    params_path = os.path.splitext(path)[0] + ".params.txt"
    with open(params_path, "w") as out:
        out.write(f"# {os.path.basename(path)} -- run parameters\n")
        out.write(f"generated: {datetime.now().astimezone().isoformat(timespec='seconds')}\n")
        out.write(f"script: scripts/locate_guides.py\n")
        for key, value in vars(args).items():
            out.write(f"{key}: {value}\n")
        for key, value in extra.items():
            out.write(f"{key}: {value}\n")
    return params_path


def main(argv=None) -> int:
    args = parse_args(argv)

    header, sequence = read_fasta(args.sequence)
    name = header.split()[0]
    print(f"sequence: {name}, {len(sequence):,} nt")

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

    # The offset matters more than it looks: getting it wrong silently slides
    # every genome coordinate, and the output would still look plausible. So
    # say which source it came from.
    if args.region_start is not None:
        region_start, source = args.region_start, "--region-start"
    else:
        from_header = region_start_from_header(header)
        region_start = from_header if from_header is not None else 1
        source = "FASTA header" if from_header is not None else "default"
    print(f"region start: {region_start:,} (from {source})")

    psps = []
    if args.psp_hits:
        hits = read_table(args.psp_hits)
        psps = build_psp_guides(hits, region_start=region_start, reference=args.reference)
        psps = drop_offsequence(psps, len(sequence))
        print(f"PspCas13b footprints: {len(psps)} from {args.psp_hits}")
        if args.psp_name:
            psps = filter_psps(psps, args.psp_name)
            # Named explicitly, because the restriction decides the overlap
            # classes written to the output, not just what a figure draws.
            print(
                f"  restricted to {len(psps)}: "
                + ", ".join(p.name for p in psps)
            )

    guides = build_guides(predictions, sequence, psps, transcript=name)
    frame = guides_to_frame(
        guides,
        region_start=region_start,
        reference=args.reference,
        length=len(sequence),
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    frame.to_csv(args.output, index=False)

    counts = frame["psp_overlap"].value_counts().to_dict()
    print(
        f"wrote {len(frame):,} guides to {args.output}\n"
        f"  positions {frame['start'].min():,}-{frame['end'].max():,} "
        f"on {name}, {frame['ref_start'].min():,}-{frame['ref_end'].max():,} "
        f"on the reference\n"
        f"  overlap with PspCas13b: "
        + ", ".join(f"{k} {v:,}" for k, v in sorted(counts.items()))
    )
    print(f"  ids like {frame['guide_id'].iloc[0]} .. {frame['guide_id'].iloc[-1]}")

    if not args.no_params:
        params_path = write_params(
            args.output,
            args,
            {
                "resolved_region_start": region_start,
                "region_start_source": source,
                "sequence_name": name,
                "sequence_length": len(sequence),
                "n_predictions": len(predictions),
                "n_psp_footprints": len(psps),
                "n_guides_written": len(frame),
            },
        )
        print(f"  parameters: {params_path}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
