#!/usr/bin/env python
"""Plot the RfxCas13d guide scores around each PspCas13b guide.

Reads the pair table from extract_neighbours.py plus the index file naming each
PspCas13b guide, and writes one interactive HTML plot per query guide to
results/plots/.

Examples
--------
    python scripts/visualize.py

Needs pandas and plotly (e.g. the `pik3-guide` conda env).
"""

from __future__ import annotations
import os
import re

import pandas as pd
from typing import Iterable, NamedTuple
import plotly.graph_objects as go

# Paths are given relative to the repository root, not to this file, so that
# results/ stays where it is now that the script lives in scripts/.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SURFACE = "#fcfcfb"      # plot background, and the ring that separates stacked points
GRID = "#e8e7e2"
INK = "#52514e"          # axis and tick text
LINE = "#d9d8d3"         # the score profile joining the points, kept recessive

# One class per trace, drawn in this order, so the legend reads from inside the
# query guide outwards. The two overlap classes carry the hue because they are
# what the plot is for; upstream and downstream are context and stay achromatic
# -- which side a point is on is already unmistakable from the shaded band it
# sits beside. Both neutrals clear 3:1 on this surface (3.4 and 7.7), and the
# symbol differs as well as the colour, so the classes do not rely on hue alone.
#
# The hues are slots 1 and 2 of the data-viz reference palette, validated for
# all pairs in light mode: CVD delta E 24.7, normal-vision 33.6, both well over
# the 8 and 15 floors. Four hues cannot clear those floors on an all-pairs
# chart, which is why only the signal classes get one.
CLASS_STYLE = {
    "full_overlap": ("Full overlap", "#2a78d6", "circle"),
    "partial_overlap": ("Partial overlap", "#eb6834", "circle"),
    "upstream": ("Upstream", "#8a8a85", "triangle-left"),
    "downstream": ("Downstream", "#52514e", "triangle-right"),
}


class QueryGuide(NamedTuple):
    index: int             # query_guide_index: how the pair table keys this guide
    experimental_id: str   # MERS_n: how the guide is named at the bench
    gene: str
    start: int
    end: int


class NeighbourGuide(NamedTuple):
    query_index: int
    index: int
    sequence: str
    score: float
    ref_start: int
    ref_end: int
    relative_position: str


def load_experimental_ids(path: str) -> dict[int, str]:
    """Map each PspCas13b guide's query_guide_index to its experimental_id.

    The pair table keys guides by their row Index in the input CSV, which means
    nothing at the bench; this index file carries the MERS_n name each guide is
    actually known by.
    """
    index = pd.read_csv(path)
    missing = [c for c in ("query_guide_index", "experimental_id") if c not in index.columns]
    if missing:
        raise ValueError(f"{path} lacks the column(s) {', '.join(missing)}")
    return {int(row.query_guide_index): str(row.experimental_id) for row in index.itertuples()}


def build_query_guide(df: pd.DataFrame, experimental_ids: dict[int, str]) -> list[QueryGuide]:
    """One QueryGuide per query guide, keyed on its Index rather than its gene.

    A gene can hold several query guides -- orf1ab has three, 16 kb apart -- so
    grouping by query_gene and taking min(start)/max(end) would collapse them
    into a single span covering the gaps between them too. A guide's own
    coordinates are the same on all of its rows, hence the first row of each.
    """
    query_guides = []
    for index, group in df.groupby("query_guide_index"):
        row = group.iloc[0]
        index = int(index)
        if index not in experimental_ids:
            print(f"warning: query guide {index} has no experimental_id in the index file")
        query_guides.append(QueryGuide(
            index=index,
            experimental_id=experimental_ids.get(index, f"guide {index}"),
            gene=row["query_gene"],
            start=int(row["query_ref_start"]),
            end=int(row["query_ref_end"]),
        ))
    return query_guides


def build_neighbour_guides(df: pd.DataFrame) -> list[NeighbourGuide]:
    """Every (query guide, neighbour) pair in the table, carrying the query's Index.

    A neighbour within reach of two query guides appears once per pair, which is
    what the table is for -- so each tuple remembers *which* query guide it was
    found for, and a plot can take only its own.

    Rows marked no_neighbour carry no neighbour at all and are skipped; they
    would otherwise become tuples of NaN.
    """
    neighbour_guides = []
    for _, row in df[df["relative_position"] != "no_neighbour"].iterrows():
        neighbour_guides.append(NeighbourGuide(
            query_index=int(row["query_guide_index"]),
            index=int(row["neighbour_guide_index"]),
            sequence=row["neighbour_guide_sequence"],
            score=row["neighbour_guide_score"],
            ref_start=int(row["neighbour_ref_start"]),
            ref_end=int(row["neighbour_ref_end"]),
            relative_position=row["relative_position"],
        ))
    return neighbour_guides


def file_slug(text: str) -> str:
    """A filename-safe gene label: ``N;3'UTR`` -> ``N_3_UTR``.

    The semicolon in a multi-gene label is a command separator in a shell, and
    the apostrophe in the UTR labels needs quoting, so neither belongs in a path.
    """
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", text)).strip("_")


def plot_query_guide(query: QueryGuide, neighbours: Iterable[NeighbourGuide]) -> None:

    # Selected on the query guide's Index, not on the gene: a neighbour is a
    # neighbour because of where it sits, and it may be annotated to an
    # adjoining feature.
    neighbour = sorted(
        (n for n in neighbours if n.query_index == query.index),
        key=lambda n: n.ref_start,
    )

    # Indexes 10-14 have no RfxCas13d guide within the window, so there is
    # nothing to draw; an empty frame has no columns at all and Plotly cannot
    # find "start".
    
    if not neighbour:
        print(f"skipping {query.experimental_id} ({query.gene}): no neighbour guides")
        return None

    # Convert the NamedTuples into a DataFrame that Plotly can read by column name
    df = pd.DataFrame(
        [
            {
                "start": n.ref_start,
                "end": n.ref_end,
                "guide_score": n.score,
                "relative_position": n.relative_position,
                "guide_index": n.index,
                "guide_sequence": n.sequence,
            }
            for n in neighbour
        ]
    )

    fig = go.Figure()
    # The profile first and in one piece: colouring it by class would break it
    # into disconnected segments, and it is the points that carry the class.
    fig.add_trace(go.Scatter(
        x=df["start"], y=df["guide_score"],
        mode="lines",
        line=dict(color=LINE, width=2),
        hoverinfo="skip",       # the markers' own hover covers every point
        showlegend=False,
        name="score profile",
    ))

    # One trace per class, in CLASS_STYLE order rather than the order the
    # classes happen to appear in, so the legend reads the same on every plot.
    # A class with no guides is left out: an empty legend entry would suggest a
    # category that could have been there at this score, rather than one that
    # cannot occur -- MERS_1 has no full overlap at all, the pool stopping at
    # 13432.
    for position, (label, colour, symbol) in CLASS_STYLE.items():
        points = df[df["relative_position"] == position]
        if points.empty:
            continue
        fig.add_trace(go.Scatter(
            x=points["start"], y=points["guide_score"],
            mode="markers",
            name=label,
            marker=dict(
                color=colour,
                symbol=symbol,
                size=9,
                # A ring in the surface colour keeps touching points separate:
                # at a 23 nt tiling they sit one base apart.
                line=dict(color=SURFACE, width=1.5),
            ),
            customdata=points[["end", "guide_index", "guide_sequence"]],
            hovertemplate=(
                f"<b>{label}</b><br>"
                "score %{y:.3f}<br>"
                "%{x:,d}-%{customdata[0]:,d}<br>"
                "Index %{customdata[1]}<br>"
                "%{customdata[2]}<extra></extra>"
            ),
        ))

    fig.update_layout(
        title=(
            f"RfxCas13d guide scores around PspCas13b {query.experimental_id} "
            f"({query.gene})"

        ),
        xaxis_title="Position",
        yaxis_title="RfxCas13d Guide Score",
        legend_title_text="Position vs query guide",
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK),
        hovermode="closest",
    )
    fig.update_yaxes(gridcolor=GRID, zeroline=False)
    # Genome coordinates, so show them in full: Plotly's default SI formatting
    # turns 22604 into "22.7k", which is useless for finding a guide.
    # Tick spacing is left to Plotly: a fixed dtick of 10 puts 25 six-character
    # labels on a 250 nt window, which collide, and it does not refine as you
    # zoom in.
    fig.update_xaxes(
        tickformat=",d",    # tick labels: 22,604
        hoverformat=",d",   # the x value in the hover box, likewise
        gridcolor=GRID,
    )
    # Neutral, not orange: the band is the query guide's footprint, and a hue
    # here would read as a fifth class next to the orange partial overlaps.
    fig.add_vrect(
        x0=query.start - 0.5, x1=query.end + 0.5,  # whole bases, not their centres
        fillcolor="#8a8a85",
        opacity=0.15,                 # lower = more see-through
        line_width=0,                 # removes the border
        layer="below",                # keeps the band behind your data
        annotation_text=f"PspCas13b {query.experimental_id}",
        annotation_position="top left",
        annotation_font_color=INK,
    )
    name = f"{file_slug(query.experimental_id)}_{file_slug(query.gene)}_guide_scores.html"
    return fig.write_html(
        os.path.join(ROOT, "results", "plots", name)
    )


def main():
    os.makedirs(os.path.join(ROOT, "results", "plots"), exist_ok=True)
    data_file = os.path.join(ROOT, "results", "psp_neighbours_100bp.csv")
    index_file = os.path.join(ROOT, "results", "psp_mers_index.csv")

    df = pd.read_csv(data_file, header=0)
    experimental_ids = load_experimental_ids(index_file)
    query_guides = build_query_guide(df, experimental_ids)
    neighbour_guides = build_neighbour_guides(df)

    for query in query_guides:
        plot_query_guide(query, neighbour_guides)


if __name__ == "__main__":
    main()
