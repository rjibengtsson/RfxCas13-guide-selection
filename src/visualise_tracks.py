"""Stacked guide-position tracks: one lane per gene, two Cas types per lane.

Each guide is drawn as a horizontal bar at its position along a gene. Guides of
one Cas type are packed above a midline and the other below it, so a reader can
see at a glance where the two sets coincide along the sequence.

This answers a different question from ``visualise_scatterlines``. There, every
guide in a tiling is a point and the y-axis carries its score; here a *selected
handful* of guides are intervals, and the y-axis carries nothing but stacking.
Position is the only variable, which is why the bars are packed as tightly as
non-overlap allows (see :func:`pack_intervals`): vertical space is a drawing
artefact, not data, so it is minimised.

Coordinates are 1-based inclusive and **local to the gene**, so a bar's span is
``end - start + 1`` and the x-axis runs from 0 to the gene's length. Positions
arriving in genome coordinates must be shifted first -- see
:func:`shift_to_gene`. Getting that shift wrong slides every bar along the axis
and the figure still looks plausible, so the shift is always reported by the
caller.

Specification: ``docs/guide-track-plot-spec.md``. This is the library half; it
parses no arguments and reads no file by a hard-coded path.
``scripts/plot_guide_tracks.py`` holds the configuration block and the CLI.
"""

from __future__ import annotations

import os
from typing import Iterable, NamedTuple, Sequence

import plotly.graph_objects as go
from plotly.subplots import make_subplots


class Interval(NamedTuple):
    """One guide as a drawable bar, in gene-local 1-based inclusive coordinates.

    ``fields`` holds the already-formatted hover lines as ``(label, value)``
    pairs. Formatting happens before this point on purpose: the two tracks draw
    different field sets from different source tables, and a bar needs to carry
    its own hover text rather than an index into a frame.
    """

    name: str
    start: int
    end: int
    fields: tuple[tuple[str, str], ...] = ()

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class Track(NamedTuple):
    """One Cas type's guides, with how to draw them.

    ``annotate`` prints each guide's name beside its bar -- above it on an
    upper track, below it on a lower one, so a label never sits over the
    midline gap where it could be read as belonging to the other track. Use it
    only where the names are few and short enough not to collide; the hover
    carries the name either way, so this is for a figure meant to be read
    without interaction.
    """

    label: str
    intervals: tuple[Interval, ...]
    color: str
    side: str = "upper"          # "upper" -> y >= 0, "lower" -> y < 0
    annotate: bool = False


def shift_to_gene(start: int, end: int, gene_start: int) -> tuple[int, int]:
    """Move a span from reference coordinates into a gene's own frame.

    *gene_start* is the gene's 1-based start in the reference, so a reference
    position ``p`` lands at ``p - gene_start + 1`` and the gene's first base is
    1. The span is ordered on the way out: a blastn hit on the minus strand
    reports ``sstart > send``, and a bar only needs the extent.
    """
    lo, hi = min(start, end), max(start, end)
    return lo - gene_start + 1, hi - gene_start + 1


def gene_span(genbank_path: str, gene: str) -> tuple[int, int]:
    """``(start, end)`` of *gene* in its record, 1-based inclusive.

    Read from the annotation rather than taken as a parameter so the gene's
    extent and the coordinate shift come from the same place -- a mismatch
    between the two is exactly the error that draws a convincing wrong figure.

    Prefers a ``gene`` feature and falls back to the ``CDS``, which is the same
    span for a single-exon viral gene but need not be in general. Raises if the
    gene is absent, naming what the file does hold.
    """
    from Bio import SeqIO            # imported here: only this function needs it

    found: dict[str, tuple[int, int]] = {}
    names: set[str] = set()
    for record in SeqIO.parse(genbank_path, "genbank"):
        for feature in record.features:
            if feature.type not in ("gene", "CDS"):
                continue
            for name in feature.qualifiers.get("gene", []):
                names.add(name)
                if name == gene and feature.type not in found:
                    # int() drops biopython's ExactPosition wrapper; start is
                    # 0-based in biopython's half-open convention, so +1.
                    found[feature.type] = (
                        int(feature.location.start) + 1,
                        int(feature.location.end),
                    )

    for feature_type in ("gene", "CDS"):
        if feature_type in found:
            return found[feature_type]

    raise ValueError(
        f"{genbank_path} has no gene or CDS named {gene!r}; it holds: "
        + (", ".join(sorted(names)) if names else "no named genes")
    )


def pack_intervals(intervals: Iterable[Interval]) -> list[tuple[Interval, int]]:
    """Assign each interval a 0-based lane, greedy first-fit by start.

    Sorted by start, each interval goes on the lowest-numbered lane whose last
    bar ends before this one begins; otherwise a new lane opens. So guides that
    do not overlap share a lane and overlapping guides stack, which is the
    whole reason the figure needs more than one row per Cas type.

    Bars are **touching, not overlapping**, when ``start == previous end + 1``:
    these are adjacent bases, so they share a lane. Drawn at 10 px they will
    appear to abut, which is accurate -- they do.

    Returned in packing order, each interval with its lane index.
    """
    placed: list[tuple[Interval, int]] = []
    lane_ends: list[int] = []        # last occupied base of each open lane

    for interval in sorted(intervals, key=lambda i: (i.start, i.end)):
        for lane, lane_end in enumerate(lane_ends):
            if interval.start > lane_end:
                lane_ends[lane] = interval.end
                placed.append((interval, lane))
                break
        else:
            lane_ends.append(interval.end)
            placed.append((interval, len(lane_ends) - 1))
    return placed


def lane_to_y(lane: int, side: str) -> int:
    """Lane index to y position: upper lanes 0, 1, 2 up; lower -1, -2 down.

    The asymmetry is deliberate. The midline is y = 0 and belongs to the upper
    track, so the two tracks are never drawn on the same row and the gap
    between them is exactly one unit.
    """
    if side == "upper":
        return lane
    if side == "lower":
        return -(lane + 1)
    raise ValueError(f"side must be 'upper' or 'lower', not {side!r}")


def hover_template(label: str, fields: Sequence[tuple[str, str]]) -> str:
    """A bar's hover box: bold track label, then one line per field.

    Values are baked in rather than passed as ``customdata``: each bar is its
    own trace here, so there is no array to index into and a literal template
    is both simpler and impossible to misalign.

    ``<extra></extra>`` suppresses Plotly's trace-name side box, which would
    otherwise repeat the track label on every hover.
    """
    lines = [f"<b>{label}</b>"]
    lines += [f"{name}: {value}" for name, value in fields]
    return "<br>".join(lines) + "<extra></extra>"


def add_track(
    fig: go.Figure,
    track: Track,
    row: int,
    line_width: int,
    show_legend: bool,
    annotation_font_size: int = 10,
) -> int:
    """Draw one track's bars on *row*. Returns the number of lanes used.

    *show_legend* is set for the first row only, so the legend carries one
    entry per Cas type rather than one per bar, and clicking it toggles that
    Cas type across every gene. ``legendgroup`` is what ties them together.
    """
    placed = pack_intervals(track.intervals)

    for index, (interval, lane) in enumerate(placed):
        y = lane_to_y(lane, track.side)

        if track.annotate:
            # Offset from the bar's own lane, not from the axis, so a label
            # stays with its bar when the track stacks. 0.32 lane-units clears
            # a 10 px bar without reaching the lane above.
            offset = 0.32 if track.side == "upper" else -0.32
            fig.add_annotation(
                text=interval.name,
                x=(interval.start + interval.end) / 2,
                y=y + offset,
                yanchor="bottom" if track.side == "upper" else "top",
                showarrow=False,
                font=dict(size=annotation_font_size, color=track.color),
                row=row,
                col=1,
            )

        fig.add_trace(
            go.Scatter(
                # Half a base either side, so a bar covers whole bases rather
                # than the span between their centres.
                x=[interval.start - 0.5, interval.end + 0.5],
                y=[y, y],
                mode="lines",
                line=dict(color=track.color, width=line_width),
                name=track.label,
                legendgroup=track.label,
                showlegend=show_legend and index == 0,
                hovertemplate=hover_template(track.label, interval.fields),
                hoverlabel=dict(namelength=-1),
            ),
            row=row,
            col=1,
        )

    return max((lane for _, lane in placed), default=-1) + 1


def plot_tracks(
    genes: Sequence[tuple[str, int, Sequence[Track]]],
    output: str,
    title: str,
    legend_title: str = "Cas type",
    x_axis_title: str = "Position (nt)",
    line_width: int = 10,
    y_pad: int = 2,
    vertical_spacing: float = 0.05,
    height_per_gene: int = 320,
    margin: dict | None = None,
    annotation_font_size: int = 10,
    subplot_title_font_size: int = 16,
    inline_plotlyjs: bool = True,
) -> str:
    """One subplot row per gene, each row packing its tracks around a midline.

    *genes* is ``(name, length, tracks)`` per row. *length* is the gene's full
    length, so the x-axis shows the whole gene and a guide's position is read
    against the gene rather than against the extent of the selection.

    Row heights scale with how many lanes a gene needs, plus one unit for the
    midline gap, so a crowded gene is not squeezed to the same height as a
    sparse one.
    """
    if not genes:
        raise ValueError("no genes to plot")

    # Pack once up front: the row heights depend on the lane counts, and
    # make_subplots needs them before any trace is added.
    lane_counts = []
    for _, _, tracks in genes:
        upper = max(
            (
                len(set(lane for _, lane in pack_intervals(t.intervals)))
                for t in tracks
                if t.side == "upper"
            ),
            default=0,
        )
        lower = max(
            (
                len(set(lane for _, lane in pack_intervals(t.intervals)))
                for t in tracks
                if t.side == "lower"
            ),
            default=0,
        )
        lane_counts.append((upper, lower))

    totals = [upper + lower + 1 for upper, lower in lane_counts]   # +1: the gap
    row_heights = [total / sum(totals) for total in totals]

    fig = make_subplots(
        rows=len(genes),
        cols=1,
        shared_xaxes=False,          # coordinates are local to each gene
        subplot_titles=[name for name, _, _ in genes],
        row_heights=row_heights,
        vertical_spacing=vertical_spacing,
    )

    # The subplot titles are the only annotations at this point, so size them
    # now. Doing it after the tracks were added would also resize the guide
    # labels, which carry their own smaller font.
    for annotation in fig.layout.annotations:
        annotation.font.size = subplot_title_font_size

    for row, ((name, length, tracks), (upper, lower)) in enumerate(
        zip(genes, lane_counts), start=1
    ):
        for track in tracks:
            add_track(
                fig,
                track,
                row,
                line_width,
                show_legend=(row == 1),
                annotation_font_size=annotation_font_size,
            )

        # 0 rather than 0.5: the axis is the gene, and starting it at the
        # gene's own origin keeps the two ends of the bar readable.
        fig.update_xaxes(range=[0, length], tickformat=",d", hoverformat=",d",
                         row=row, col=1)
        # The y-axis carries no data, only stacking, so it is unlabelled.
        fig.update_yaxes(
            range=[-lower - 0.5 - y_pad, upper - 0.5 + y_pad],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            title_text=None,
            row=row,
            col=1,
        )

    # Bottom row only: repeating the axis title on every row would say the same
    # thing n times, and the rows share the same unit.
    fig.update_xaxes(title_text=x_axis_title, row=len(genes), col=1)

    fig.update_layout(
        title=title,
        height=height_per_gene * len(genes),
        autosize=True,
        margin=margin or dict(l=30, r=30, t=60, b=20),
        # itemsizing="constant" stops the legend swatch inheriting the 10 px
        # bar width, which would make the two entries different sizes.
        legend=dict(title=legend_title, itemsizing="constant"),
        hovermode="closest",
    )

    directory = os.path.dirname(os.path.abspath(output))
    os.makedirs(directory, exist_ok=True)
    fig.write_html(output, include_plotlyjs=True if inline_plotlyjs else "cdn")
    return output
