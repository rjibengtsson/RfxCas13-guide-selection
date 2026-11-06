"""Pull out the guides sitting within and around each guide of another Cas type.

Works on the annotated CSV that ``scripts/blast_guides.py`` writes, using its
genome coordinates. For every query guide (by default each PspCas13b guide) it
collects the neighbour guides (by default RfxCas13d) that overlap it or fall
within a given distance, and reports both guides' details side by side.

The window may be asymmetric: *upstream* and *downstream* are counted in genome
coordinates, so upstream means the lower coordinates. Guides overlapping the
query are always included, whatever the window, and are reported as a full or a
partial overlap.
"""

from __future__ import annotations

import os
from typing import Sequence

import numpy as np
import pandas as pd

QUERY_CAS = "PspCas13b"
NEIGHBOUR_CAS = "RfxCas13d"
DEFAULT_WINDOW = 100

FULL_OVERLAP = "full_overlap"
PARTIAL_OVERLAP = "partial_overlap"
UPSTREAM = "upstream"
DOWNSTREAM = "downstream"
NO_NEIGHBOUR = "no_neighbour"

# The CSV's display names mapped to identifiers safe for attribute access;
# 'Index' in particular collides with pandas' own row index.
COLUMN_ALIASES = {
    "Cas": "cas",
    "Index": "guide_index",
    "Guide Sequence": "guide_sequence",
    "Target Sequence": "target_sequence",
    "Guide Score": "guide_score",
}

REQUIRED_COLUMNS = (
    "Cas",
    "Index",
    "Guide Sequence",
    "Target Sequence",
    "Guide Score",
    "ref_start",
    "ref_end",
    "strand",
    "features",
)

# Deliberately no BLAST statistics: identity, e-value and the rest live in the
# annotated CSV and are not what this table is for.
QUERY_FIELDS = (
    "cas",
    "guide_index",
    "gene",
    "guide_sequence",
    "target_sequence",
    "ref_start",
    "ref_end",
)
# No cas or strand: every neighbour is the one neighbour_cas asked for, and the
# guides are blasted as given, so the strand is the same on every row. No gene
# either: a neighbour is a neighbour by position, and its coordinates say where
# it sits -- the annotated CSV still carries its features.
NEIGHBOUR_FIELDS = (
    "guide_index",
    "guide_sequence",
    "target_sequence",
    "guide_score",
    "ref_start",
    "ref_end",
)

OUTPUT_COLUMNS = (
    [f"query_{f}" for f in QUERY_FIELDS]
    + ["n_neighbours", "relative_position"]
    + [f"neighbour_{f}" for f in NEIGHBOUR_FIELDS]
)


class NeighbourDataError(ValueError):
    """Raised when the input CSV cannot support the requested extraction."""


def gene_label(features) -> str:
    """The gene(s) a guide sits in, from the annotation's ``features`` column.

    Takes each overlapping feature's name, de-duplicated and in positional
    order, so a guide in the ORF1ab/S overlap reads ``orf1ab;S``. The name falls
    back to the feature type where a feature carries no /gene qualifier, which is
    how the UTRs appear.
    """
    if pd.isna(features):
        return ""
    names = {}
    for entry in str(features).split(";"):
        f_type, _, name = entry.partition(":")
        names[name or f_type] = None
    return ";".join(names)


def load_annotated(path: str) -> pd.DataFrame:
    """Read an annotated guide CSV, adding a single ``gene`` label per guide."""
    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise NeighbourDataError(
            f"{path} lacks the column(s) {', '.join(missing)}; it must be the "
            "annotated output of blast_guides.py run against a GenBank reference"
        )
    df = df.rename(columns=COLUMN_ALIASES)
    # Unmapped guides have no coordinates, so they cannot take part.
    df = df[df.ref_start.notna() & df.ref_end.notna()].copy()
    df["ref_start"] = df.ref_start.astype(int)
    df["ref_end"] = df.ref_end.astype(int)
    df["gene"] = df.features.map(gene_label)
    return df


def _interval_gap(a_start, a_end, b_start, b_end):
    """Bases separating two 1-based inclusive intervals.

    0 when they overlap or merely touch, 1 when a single base lies between them.
    Broadcasts, so columns against scalars give a whole column of distances.
    """
    return np.maximum(0, np.maximum(b_start - a_end, a_start - b_end) - 1)


def _positions(q_start, q_end, n_start, n_end):
    """Label each neighbour's position relative to the query guide.

    An overlap is split in two, because the two cases mean different things at
    the bench: ``full_overlap`` is one guide's footprint lying wholly inside the
    other's -- for a 23 nt neighbour against a 30 nt query, the neighbour
    targeting entirely within the query's site -- while ``partial_overlap`` is
    an overlap that runs off one end, which at the extreme is a single shared
    base. Containment is tested both ways round, so a neighbour longer than the
    query counts as a full overlap too.

    Guides that do not overlap at all are ``upstream`` or ``downstream``, in
    genome coordinates; a merely abutting guide is one of those, not an overlap.
    """
    overlapping = (n_start <= q_end) & (n_end >= q_start)
    contained = ((n_start >= q_start) & (n_end <= q_end)) | (
        (q_start >= n_start) & (q_end <= n_end)
    )
    return np.select(
        [overlapping & contained, overlapping, n_end < q_start],
        [FULL_OVERLAP, PARTIAL_OVERLAP, UPSTREAM],
        default=DOWNSTREAM,
    )


def find_neighbours(
    df: pd.DataFrame,
    upstream: int = DEFAULT_WINDOW,
    downstream: int = DEFAULT_WINDOW,
    query_cas: str = QUERY_CAS,
    neighbour_cas: str = NEIGHBOUR_CAS,
    query_indexes: Sequence | None = None,
    include_empty: bool = True,
) -> pd.DataFrame:
    """Collect the *neighbour_cas* guides within the window of each query guide.

    A neighbour qualifies when it overlaps the query guide -- wholly
    (``full_overlap``) or in part (``partial_overlap``) -- or lies at most
    *upstream* bases below it, or at most *downstream* bases above it, measured
    between genome intervals on the same reference sequence. Gene boundaries are
    not consulted: a neighbour in an adjoining feature still counts.

    One row per (query guide, neighbour) pair, so a neighbour close to two query
    guides is reported under both. With *include_empty*, a query guide with no
    neighbour gets a single row marked ``no_neighbour`` instead of vanishing.
    Rows run in query order, then by the neighbour's genome position.
    """
    if upstream < 0 or downstream < 0:
        raise NeighbourDataError("upstream and downstream must be 0 or more")

    queries = df[df.cas == query_cas]
    if query_indexes is not None:
        queries = queries[queries.guide_index.isin(list(query_indexes))]
    if queries.empty:
        raise NeighbourDataError(
            f"no {query_cas} guides found"
            + (" for the requested Index values" if query_indexes is not None else "")
            + f"; the file holds: {', '.join(sorted(df.cas.unique()))}"
        )
    pool = df[df.cas == neighbour_cas]

    frames = []
    for q in queries.itertuples(index=False):
        candidates = pool
        if "ref_accession" in pool.columns:
            candidates = pool[pool.ref_accession == q.ref_accession]

        gap = _interval_gap(
            candidates.ref_start.to_numpy(),
            candidates.ref_end.to_numpy(),
            q.ref_start,
            q.ref_end,
        )
        position = _positions(
            q.ref_start,
            q.ref_end,
            candidates.ref_start.to_numpy(),
            candidates.ref_end.to_numpy(),
        )
        # Each side gets its own allowance; an overlap is always kept.
        allowance = np.select(
            [position == UPSTREAM, position == DOWNSTREAM], [upstream, downstream], default=gap
        )
        # The gap decides the window but is not reported: it saturates at 0 for
        # anything overlapping or abutting, which reads as a distance of zero
        # between guides that may share all 23 bases. relative_position says
        # where the neighbour sits, and the coordinates say how far.
        hits = candidates.assign(relative_position=position)[gap <= allowance]
        hits = hits.sort_values(["ref_start", "ref_end"])

        header = {f"query_{f}": getattr(q, f) for f in QUERY_FIELDS}
        header["n_neighbours"] = len(hits)

        if hits.empty:
            if include_empty:
                frames.append(
                    pd.DataFrame([{**header, "relative_position": NO_NEIGHBOUR}])
                )
            continue

        block = hits[list(NEIGHBOUR_FIELDS) + ["relative_position"]].rename(
            columns={f: f"neighbour_{f}" for f in NEIGHBOUR_FIELDS}
        )
        for key, value in header.items():
            block.insert(0, key, value)
        frames.append(block)

    if not frames:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    # reindex, not selection: with no neighbours anywhere, no frame carries the
    # neighbour columns at all.
    pairs = pd.concat(frames, ignore_index=True).reindex(columns=OUTPUT_COLUMNS)
    return _restore_integers(pairs)


# Whole-numbered columns, which the blank rows of a query without neighbours
# would otherwise turn into floats ("212.0" in the CSV).
INTEGER_COLUMNS = (
    "n_neighbours",
    "neighbour_guide_index",
    "neighbour_ref_start",
    "neighbour_ref_end",
)


def _restore_integers(pairs: pd.DataFrame) -> pd.DataFrame:
    for column in INTEGER_COLUMNS:
        try:
            pairs[column] = pairs[column].astype("Int64")
        except (TypeError, ValueError):
            pass  # not whole-numbered in this input; leave it as read
    return pairs


def summarise(pairs: pd.DataFrame) -> pd.DataFrame:
    """One row per query guide: where it sits and how its neighbours fall."""
    if pairs.empty:
        return pd.DataFrame()
    real = pairs[pairs.relative_position != NO_NEIGHBOUR]
    counts = (
        real.groupby("query_guide_index").relative_position.value_counts().unstack(fill_value=0)
        if not real.empty
        else pd.DataFrame()
    )
    summary = (
        pairs.drop_duplicates(subset=["query_guide_index"])
        .set_index("query_guide_index")[
            ["query_gene", "query_ref_start", "query_ref_end", "n_neighbours"]
        ]
        .join(counts)
    )
    for col in (FULL_OVERLAP, PARTIAL_OVERLAP, UPSTREAM, DOWNSTREAM):
        if col not in summary.columns:
            summary[col] = 0
        summary[col] = summary[col].fillna(0).astype(int)
    return summary[
        ["query_gene", "query_ref_start", "query_ref_end", "n_neighbours",
         FULL_OVERLAP, PARTIAL_OVERLAP, UPSTREAM, DOWNSTREAM]
    ]


def write_neighbours(pairs: pd.DataFrame, path: str) -> str:
    """Write the pair table to *path*, creating its directory if needed."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    pairs.to_csv(path, index=False)
    return path
