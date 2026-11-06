#!/usr/bin/env python
"""Extract the guides sitting within and around each guide of another Cas type.

Reads the annotated CSV from blast_guides.py and writes one row per
(PspCas13b guide, nearby RfxCas13d guide) pair, carrying both guides' details
and the distance between them. The window may be asymmetric.

Examples
--------
    # every PspCas13b guide, 100 bp either side
    python scripts/extract_neighbours.py \
        --input results/mers_guides_NC_019843_annotated.csv \
        --output results/psp_neighbours_100bp.csv --window 100

    # one guide, with a wider view downstream
    python scripts/extract_neighbours.py \
        --input results/mers_guides_NC_019843_annotated.csv \
        --output results/psp0_neighbours.csv \
        --query-sequence TTAAGTGAATAGCTTGGCTATCTCACTTCC \
        --upstream 100 --downstream 300

Needs pandas (e.g. the `pik3-guide` conda env).
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from guide_neighbours import (  # noqa: E402
    DEFAULT_WINDOW,
    NEIGHBOUR_CAS,
    NO_NEIGHBOUR,
    QUERY_CAS,
    NeighbourDataError,
    find_neighbours,
    load_annotated,
    summarise,
    write_neighbours,
)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, help="annotated CSV from blast_guides.py")
    p.add_argument("--output", required=True, help="output CSV of guide pairs")
    p.add_argument("--summary", help="optional CSV with one row per query guide")

    w = p.add_argument_group("window")
    w.add_argument(
        "--window",
        type=int,
        default=DEFAULT_WINDOW,
        metavar="BP",
        help=f"symmetric window either side of the query guide (default: {DEFAULT_WINDOW})",
    )
    w.add_argument(
        "--upstream",
        type=int,
        metavar="BP",
        help="bases to search below the query guide, overriding --window",
    )
    w.add_argument(
        "--downstream",
        type=int,
        metavar="BP",
        help="bases to search above the query guide, overriding --window",
    )

    s = p.add_argument_group("selection")
    s.add_argument("--query-cas", default=QUERY_CAS, help=f"guides to search around (default: {QUERY_CAS})")
    s.add_argument(
        "--neighbour-cas", default=NEIGHBOUR_CAS, help=f"guides to collect (default: {NEIGHBOUR_CAS})"
    )
    s.add_argument(
        "--query-index",
        help="comma-separated Index values to restrict the query guides to",
    )
    s.add_argument(
        "--query-sequence",
        help="restrict to the query guide with this Guide or Target Sequence",
    )
    s.add_argument(
        "--drop-empty",
        action="store_true",
        help="omit query guides that have no neighbour, instead of marking them",
    )
    return p.parse_args(argv)


def resolve_indexes(df, args) -> list | None:
    """Turn --query-index / --query-sequence into a list of Index values."""
    indexes = None
    if args.query_index:
        indexes = [int(i) for i in args.query_index.split(",") if i.strip()]
    if args.query_sequence:
        wanted = args.query_sequence.strip().upper().replace("U", "T")
        hit = df[
            (df.cas == args.query_cas)
            & (
                (df.guide_sequence.str.upper() == wanted)
                | (df.target_sequence.str.upper() == wanted)
            )
        ]
        if hit.empty:
            raise NeighbourDataError(
                f"no {args.query_cas} guide has the sequence {args.query_sequence!r} "
                "in its Guide Sequence or Target Sequence column"
            )
        found = list(hit.guide_index)
        indexes = found if indexes is None else [i for i in indexes if i in found]
    return indexes


def main(argv=None) -> int:
    args = parse_args(argv)

    upstream = args.window if args.upstream is None else args.upstream
    downstream = args.window if args.downstream is None else args.downstream
    if upstream < 0 or downstream < 0:
        print("error: window sizes must be 0 or more", file=sys.stderr)
        return 2

    df = load_annotated(args.input)
    indexes = resolve_indexes(df, args)

    pairs = find_neighbours(
        df,
        upstream=upstream,
        downstream=downstream,
        query_cas=args.query_cas,
        neighbour_cas=args.neighbour_cas,
        query_indexes=indexes,
        include_empty=not args.drop_empty,
    )
    write_neighbours(pairs, args.output)

    n_queries = pairs.query_guide_index.nunique() if not pairs.empty else 0
    found = pairs[pairs.relative_position != NO_NEIGHBOUR] if not pairs.empty else pairs
    print(
        f"searched {upstream} bp upstream and {downstream} bp downstream of "
        f"{n_queries} {args.query_cas} guide(s)"
    )
    print(
        f"found {len(found)} {args.neighbour_cas} guide(s) "
        f"({found.neighbour_guide_index.nunique() if len(found) else 0} distinct) "
        f"-> {args.output}"
    )

    summary = summarise(pairs)
    if not summary.empty:
        empty = int((summary.n_neighbours == 0).sum())
        if empty:
            print(f"{empty} {args.query_cas} guide(s) had no {args.neighbour_cas} neighbour")
        print()
        print(summary.to_string())
        if args.summary:
            os.makedirs(os.path.dirname(os.path.abspath(args.summary)), exist_ok=True)
            summary.to_csv(args.summary)
            print(f"\nsummary -> {args.summary}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (NeighbourDataError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
