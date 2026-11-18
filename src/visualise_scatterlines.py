"""Plot RfxCas13d predicted guides along a target sequence, with the PspCas13b
footprints marked.

Each RfxCas13d guide is classified by how it sits against the PspCas13b guides
on the same sequence:

``full``
    the RfxCas13d guide lies **entirely inside** a PspCas13b footprint.
``partial``
    it shares at least one base with a footprint but hangs over an edge.
``none``
    it shares no base with any footprint.

The RfxCas13d guides here are 23 nt and the PspCas13b guides 30 nt, so "full"
means the shorter guide is contained in the longer -- containment is the only
sense in which a 23 nt guide can fully overlap a 30 nt one. A guide is
classified once, by its **strongest** relationship to any footprint, so this
module emits one record per guide rather than one per (guide, footprint) pair.

Score is on the y-axis and is deliberately **not** encoded in colour as well:
colour carries the overlap state, which is the question the plot exists to
answer. No threshold is applied and nothing is filtered -- every predicted
guide is drawn.

Coordinates are 1-based inclusive, so a guide spans ``end - start + 1`` bases
and ``sequence[start - 1:end]`` is the match. Sequences are 5' -> 3'.
``spacer_seq`` is the guide, never the target.

The predictions file carries no coordinates, only the 23 nt target each guide
was called against, so positions are recovered by locating that target in the
reference sequence -- see :func:`locate_targets`.

This is the library half of the stage: it parses nothing from the command line
and reads no file by a hard-coded path. ``scripts/plot_guides.py`` is the CLI
around it.
"""

from __future__ import annotations

import os
import re
import sys
from typing import Iterable, NamedTuple, Sequence

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


SURFACE = "#fcfcfb"      # plot background, and the ring that separates stacked points
GRID = "#e8e7e2"
INK = "#52514e"          # axis and tick text
LINE = "#d9d8d3"         # the score profile joining the points, kept recessive

# The PspCas13b footprint: a band behind the data and a lane above it. Red is
# the one hue in the figure that is not a guide state, so a band reads as the
# footprint rather than as another class of guide.
PSP_BAND = "#d6453f"
PSP_LANE_FILL = "#c2403b"

# Marker fill by overlap state. Grey is the default and is deliberately
# recessive -- most guides touch no footprint, and they are context rather than
# the subject. Orange and blue are the two states worth finding.
#
# The two hues are slots 1 and 2 of the data-viz reference palette, taken
# as-is because that palette is already validated for colour-vision
# deficiency; picking two hues by eye is exactly what that validation exists to
# prevent. They are also ordered by how much of the guide is involved, so the
# stronger relationship gets the stronger colour.
NONE_FILL = "#b0aea8"     # no overlap with any PspCas13b guide
PARTIAL_FILL = "#eb6834"  # overlaps, but hangs over an edge
FULL_FILL = "#2a78d6"     # lies entirely inside a footprint

# label, outline colour, marker size. Drawn in this order, so the states worth
# seeing land on top of the grey mass rather than under it.
#
# The ring is in SURFACE on every state: at a 1 nt tiling step adjacent points
# touch, and a hairline in the background colour is what keeps them countable.
# The overlapping states are additionally larger, because a fill change alone
# is easy to lose among a thousand points.
CLASS_STYLE = {
    "none": ("No overlap", NONE_FILL, 5),
    "partial": ("Partial overlap", PARTIAL_FILL, 9),
    "full": ("Full overlap", FULL_FILL, 11),
}
CLASS_ORDER = ("none", "partial", "full")

PREDICTION_COLUMNS = (
    "Transcript ID",
    "Target Sequence",
    "Guide Sequence",
    "Guide Score",
)
BLASTN_COLUMNS = ("qseqid", "sseqid", "sstart", "send")

# Lane geometry, as fractions of the figure height. The lane is deliberately
# thin: it is a positional reference for the scatter below, not the subject.
LANE_HEIGHT = 0.03
LANE_GAP = 0.03


class PspGuide(NamedTuple):
    """One PspCas13b guide footprint, in the plotted sequence's coordinates."""

    name: str
    reference: str
    start: int                  # 1-based inclusive, always <= end
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class RfxGuide(NamedTuple):
    """One predicted RfxCas13d guide, classified against the footprints.

    ``psp`` is the footprint the classification refers to: the one this guide
    sits inside, or failing that the nearest. It is None only when the sequence
    carries no footprint at all.
    """

    transcript: str
    target_seq: str             # the 23 nt target, 5' -> 3'
    spacer_seq: str             # the guide, 5' -> 3'
    score: float
    start: int                  # 1-based inclusive
    end: int
    psp: PspGuide | None = None
    classification: str = "none"
    distance: int | None = None  # signed nt gap; negative upstream, 0 on overlap


def read_fasta(path: str) -> tuple[str, str]:
    """``(header, sequence)`` of the first record, sequence uppercased."""
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    break           # a second record: the first is the subject
                header = line[1:].rstrip("\n")
            else:
                chunks.append(line.strip())
    if header is None:
        raise ValueError(f"{path} holds no FASTA record")
    sequence = "".join(chunks).upper()
    if not sequence:
        raise ValueError(f"{path} record {header!r} carries no sequence")
    return header, sequence


def region_start_from_header(header: str) -> int | None:
    """The sub-sequence's 1-based start in its source reference, if stated.

    ``scripts/extract_feature.py`` writes a header like
    ``MT007544.1_N MT007544.1:28274-29533(+) gene length=1260``, so an
    extracted gene carries the offset needed to put genome-coordinate hits onto
    its own coordinates. Returns None when the header says nothing, leaving the
    caller to require it explicitly rather than guess at 1.
    """
    match = re.search(r":(\d+)-(\d+)", header)
    return int(match.group(1)) if match else None


def locate_targets(
    df: pd.DataFrame, sequence: str, column: str = "Target Sequence"
) -> list[int]:
    """1-based start of each row's target within *sequence*.

    The predictions file states no coordinates, so position is recovered by
    matching. Both failure modes are raised rather than worked around, because
    either one would silently put guides at the wrong place:

    * a target that does not occur at all -- the predictions were made against
      a different sequence than the one being plotted;
    * a target occurring more than once -- the position is genuinely ambiguous,
      and picking the first would be a guess.

    The error names the offending rows, since a handful of bad rows in an
    otherwise good file is the likely case.
    """
    starts: list[int] = []
    missing: list[str] = []
    repeated: list[str] = []

    for row, target in enumerate(df[column].astype(str).str.upper()):
        count = sequence.count(target)
        if count == 0:
            missing.append(f"row {row}: {target}")
        elif count > 1:
            repeated.append(f"row {row}: {target} ({count}x)")
        starts.append(sequence.find(target) + 1)

    def summarise(items: list[str]) -> str:
        head = "; ".join(items[:5])
        return head + (f"; ... and {len(items) - 5} more" if len(items) > 5 else "")

    if missing:
        raise ValueError(
            f"{len(missing)} target sequence(s) do not occur in the reference "
            f"({len(sequence)} nt); the predictions may be for a different "
            f"sequence: {summarise(missing)}"
        )
    if repeated:
        raise ValueError(
            f"{len(repeated)} target sequence(s) occur more than once, so their "
            f"position is ambiguous: {summarise(repeated)}"
        )
    return starts


def build_psp_guides(
    df: pd.DataFrame, region_start: int = 1, reference: str | None = None
) -> list[PspGuide]:
    """Footprints from a blastn table, shifted onto the plotted coordinates.

    *region_start* is where the plotted sequence begins in the reference the
    hits were called against, so a hit at genome position g lands at
    ``g - region_start + 1``. It defaults to 1, i.e. the hits are already in the
    plotted sequence's own coordinates.

    One footprint per query, taking its best HSP by bitscore when the column is
    present. A hit falling outside the plotted sequence is dropped with a note
    on stderr rather than silently clipped: it is a real hit, just not on the
    region being drawn.
    """
    if df.empty:
        return []

    if "bitscore" in df.columns:
        df = df.loc[df.groupby("qseqid")["bitscore"].idxmax()]
    else:
        df = df.drop_duplicates("qseqid")

    guides: list[PspGuide] = []
    for _, row in df.iterrows():
        if reference is not None and str(row["sseqid"]) != reference:
            continue
        # blastn reports sstart > send on a minus-strand hit; the span is the
        # same either way and the footprint only needs it ordered.
        lo, hi = int(row["sstart"]), int(row["send"])
        start = min(lo, hi) - region_start + 1
        end = max(lo, hi) - region_start + 1
        guides.append(
            PspGuide(
                name=str(row["qseqid"]),
                reference=str(row["sseqid"]),
                start=start,
                end=end,
            )
        )
    return guides


def drop_offsequence(
    guides: Iterable[PspGuide], length: int
) -> list[PspGuide]:
    """Keep only footprints that fall within ``1..length``, noting any dropped."""
    kept, dropped = [], []
    for guide in guides:
        if guide.end < 1 or guide.start > length:
            dropped.append(guide)
        else:
            kept.append(guide)
    for guide in dropped:
        print(
            f"note: {guide.name} at {guide.start}-{guide.end} falls outside the "
            f"plotted sequence (1-{length}); not drawn",
            file=sys.stderr,
        )
    return kept


def filter_psps(
    guides: Iterable[PspGuide], names: Iterable[str]
) -> list[PspGuide]:
    """Keep only the footprints named in *names*, in their original order.

    Selecting a subset changes what the overlap classes mean -- a guide is
    classified against the footprints it is given, so restricting the set makes
    guides that only touched a dropped footprint fall back to ``none``. That is
    the intended behaviour, and the reason this filter belongs upstream of
    :func:`build_guides` rather than at drawing time: a figure whose colours
    were decided by a footprint it does not draw cannot be read.

    A name that matches nothing is an error, not an empty selection. The likely
    cause is a typo, and the alternative is a figure that is silently missing a
    footprint the caller believes is on it.
    """
    wanted = list(dict.fromkeys(names))     # de-duplicate, keep caller's order
    if not wanted:
        return list(guides)

    guides = list(guides)
    available = {g.name for g in guides}
    unknown = [n for n in wanted if n not in available]
    if unknown:
        raise ValueError(
            "no footprint named "
            + ", ".join(map(repr, unknown))
            + "; available: "
            + ", ".join(sorted(available))
        )
    return [g for g in guides if g.name in set(wanted)]


def classify_overlap(start: int, end: int, psp: PspGuide) -> str:
    """``full``, ``partial`` or ``none`` for one guide against one footprint.

    ``full`` is containment: the guide lies entirely within the footprint. With
    a 23 nt guide and a 30 nt footprint that is the only way a guide can be
    fully covered, and it is the case worth separating -- a contained guide
    targets the same site as the PspCas13b guide, while one hanging over an
    edge only shares part of it.
    """
    if end < psp.start or start > psp.end:
        return "none"
    return "full" if start >= psp.start and end <= psp.end else "partial"


def signed_distance(start: int, end: int, psp: PspGuide) -> int:
    """Gap in nt to *psp*: negative upstream, positive downstream, 0 if they
    share a base."""
    if end < psp.start:
        return end - psp.start
    if start > psp.end:
        return start - psp.end
    return 0


def best_match(
    start: int, end: int, psps: Sequence[PspGuide]
) -> tuple[PspGuide | None, str, int | None]:
    """The footprint that decides a guide's class, with class and signed gap.

    Ranked by the strength of the relationship -- ``full`` beats ``partial``
    beats ``none`` -- and, among equals, by the smaller absolute gap. So a
    guide inside a footprint is reported against that one, and a guide touching
    nothing is reported against its nearest neighbour, which is what the hover
    needs in order to say how far away it is.
    """
    if not psps:
        return None, "none", None

    rank = {"full": 2, "partial": 1, "none": 0}
    best: tuple[PspGuide, str, int] | None = None
    for psp in psps:
        classification = classify_overlap(start, end, psp)
        distance = signed_distance(start, end, psp)
        candidate = (psp, classification, distance)
        if best is None or (
            rank[classification],
            -abs(distance),
        ) > (rank[best[1]], -abs(best[2])):
            best = candidate
    return best


def build_guides(
    df: pd.DataFrame,
    sequence: str,
    psps: Sequence[PspGuide] = (),
    transcript: str | None = None,
) -> list[RfxGuide]:
    """One :class:`RfxGuide` per prediction row, classified against *psps*.

    Nothing is filtered: every row in *df* yields exactly one guide. Positions
    come from :func:`locate_targets`, so the frame and the sequence must agree.
    """
    starts = locate_targets(df, sequence)

    guides: list[RfxGuide] = []
    for (_, row), start in zip(df.iterrows(), starts):
        target = str(row["Target Sequence"]).upper()
        end = start + len(target) - 1
        psp, classification, distance = best_match(start, end, psps)
        guides.append(
            RfxGuide(
                transcript=transcript or str(row["Transcript ID"]),
                target_seq=target,
                spacer_seq=str(row["Guide Sequence"]).upper(),
                score=float(row["Guide Score"]),
                start=start,
                end=end,
                psp=psp,
                classification=classification,
                distance=distance,
            )
        )
    return guides


def pad_width(length: int) -> int:
    """Digits to zero-pad a coordinate to, for a sequence of *length* bases.

    At least 4, so a short sequence still gives ids that line up, and wider for
    anything over 9,999 bases -- the point of the padding is that ids sort
    correctly as text, which breaks the moment one coordinate is wider than
    the rest.
    """
    return max(4, len(str(int(length))))


def guide_id(transcript: str, start: int, end: int, width: int = 4) -> str:
    """``MT007544.1_N_0004-0026``: the guide's span, zero-padded.

    The position is in the id rather than only in its own column so a guide
    stays identifiable once it has been copied into a spreadsheet, an order
    form or a figure label. Padding keeps text sorting in coordinate order.
    """
    return f"{transcript}_{start:0{width}d}-{end:0{width}d}"


def guides_to_frame(
    guides: Iterable[RfxGuide],
    region_start: int = 1,
    reference: str | None = None,
    length: int | None = None,
) -> pd.DataFrame:
    """One row per guide: ids, both coordinate systems, sequences and class.

    Carries **two** coordinate systems on purpose. ``start``/``end`` are on the
    plotted sequence, which is what the ids and the figure use; ``ref_start``/
    ``ref_end`` are on the reference the region was cut from, so a guide can be
    matched back to the earlier blast output without redoing the arithmetic.
    *region_start* is where the plotted sequence begins in that reference.

    Sorted by position, so the file reads along the sequence.
    """
    guides = sorted(guides, key=lambda g: (g.start, g.end))
    if not guides:
        return pd.DataFrame()

    width = pad_width(length if length is not None else max(g.end for g in guides))

    rows = []
    for g in guides:
        rows.append(
            {
                "guide_id": guide_id(g.transcript, g.start, g.end, width),
                "transcript": g.transcript,
                "start": g.start,
                "end": g.end,
                "length": g.end - g.start + 1,
                "reference": reference or (g.psp.reference if g.psp else ""),
                "ref_start": g.start + region_start - 1,
                "ref_end": g.end + region_start - 1,
                "target_sequence": g.target_seq,
                "guide_sequence": g.spacer_seq,
                "guide_score": g.score,
                "psp_overlap": g.classification,
                "nearest_psp": g.psp.name if g.psp else "",
                "psp_start": g.psp.start if g.psp else "",
                "psp_end": g.psp.end if g.psp else "",
                "distance_nt": g.distance if g.distance is not None else "",
            }
        )
    return pd.DataFrame(rows)


#: Columns :func:`guides_to_frame` writes and :func:`frame_to_guides` reads
#: back. The round trip is what lets position be computed once, in
#: ``scripts/locate_guides.py``, and merely read by everything downstream.
LOCATED_COLUMNS = (
    "guide_id",
    "transcript",
    "start",
    "end",
    "target_sequence",
    "guide_sequence",
    "guide_score",
    "psp_overlap",
    "nearest_psp",
    "psp_start",
    "psp_end",
    "distance_nt",
)


def _blank(value) -> bool:
    """True for an empty or missing cell, however pandas chose to render it."""
    return value is None or value == "" or pd.isna(value)


def psps_from_frame(df: pd.DataFrame, reference: str = "") -> list[PspGuide]:
    """Rebuild the footprints named in a located table.

    Only the footprints some guide was classified against can be recovered this
    way -- a footprint with no guide anywhere near it leaves no trace in the
    table. In this pipeline the predictions tile the whole sequence at step 1,
    so every footprint on it is the nearest to something; a caller that needs
    the guarantee should pass the blastn hits instead.
    """
    if df.empty or "nearest_psp" not in df.columns:
        return []

    seen: dict[str, PspGuide] = {}
    for name, start, end in zip(df["nearest_psp"], df["psp_start"], df["psp_end"]):
        if _blank(name) or _blank(start) or _blank(end):
            continue
        name = str(name)
        if name not in seen:
            seen[name] = PspGuide(
                name=name,
                reference=reference,
                start=int(start),
                end=int(end),
            )
    return sorted(seen.values(), key=lambda p: p.start)


def frame_to_guides(df: pd.DataFrame, psps: Sequence[PspGuide] = ()) -> list[RfxGuide]:
    """Read back the guides a located table holds, positions and classes intact.

    The stored ``psp_overlap``/``distance_nt`` are kept rather than recomputed:
    the located file is the single source of position, and silently recomputing
    would let a plot disagree with the table it was drawn from. Pass *psps* to
    attach the footprint objects; without them the names still travel in the
    frame, and are recovered by :func:`psps_from_frame`.
    """
    missing = [c for c in LOCATED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            "located table is missing column(s) "
            + ", ".join(map(repr, missing))
            + "; was it written by scripts/locate_guides.py?"
        )

    by_name = {p.name: p for p in psps}
    guides: list[RfxGuide] = []
    for _, row in df.iterrows():
        name = row["nearest_psp"]
        psp = None if _blank(name) else by_name.get(str(name))
        if psp is None and not _blank(name) and not _blank(row["psp_start"]):
            # The table names a footprint the caller did not supply; take it
            # from the row rather than dropping the association.
            psp = PspGuide(
                name=str(name),
                reference="",
                start=int(row["psp_start"]),
                end=int(row["psp_end"]),
            )
        guides.append(
            RfxGuide(
                transcript=str(row["transcript"]),
                target_seq=str(row["target_sequence"]).upper(),
                spacer_seq=str(row["guide_sequence"]).upper(),
                score=float(row["guide_score"]),
                start=int(row["start"]),
                end=int(row["end"]),
                psp=psp,
                classification=str(row["psp_overlap"]),
                distance=None if _blank(row["distance_nt"]) else int(row["distance_nt"]),
            )
        )
    return guides


def check_against_sequence(guides: Iterable[RfxGuide], sequence: str) -> None:
    """Raise unless every guide's stored target really sits at its coordinates.

    A located table carries positions computed elsewhere, so this is what stops
    one made against a different sequence -- a different isolate, a different
    gene, a re-extracted region -- from being plotted against this one. Without
    it the figure would render perfectly and be wrong.
    """
    bad = []
    for guide in guides:
        if sequence[guide.start - 1 : guide.end] != guide.target_seq:
            bad.append(f"{guide.start}-{guide.end}")
        if len(bad) > 5:
            break
    if bad:
        raise ValueError(
            f"{len(bad)}+ guide(s) do not match the sequence at their stated "
            f"position ({', '.join(bad[:5])}); the located table and the FASTA "
            "disagree -- check they are the same region"
        )


def profile_with_gaps(
    df: pd.DataFrame,
) -> tuple[list[int | None], list[float | None]]:
    """``(xs, ys)`` for the score profile, with ``None`` wherever a gap falls.

    Plotly breaks a line at a ``None``, so the profile is drawn as one trace
    per covered stretch rather than one continuous line. A break is inserted
    between two consecutive guides that **do not overlap**: sorted by start,
    the next guide begins after every earlier guide has ended. Guides that
    merely abut -- ``start == previous end + 1`` -- share no base and so are
    also broken apart.

    Why it matters: the line is a profile of a *tiled* region, so joining
    across a gap would draw a score for positions no guide covers. Over the
    step-1 tiling in this pipeline nothing is dropped and no break appears,
    which is itself worth being able to see.

    The comparison is against the running maximum end, not the previous row's,
    so a guide nested inside an earlier, longer one still counts as covered.
    """
    if df.empty:
        return [], []

    ordered = df.sort_values(["start", "end"])

    xs: list[int | None] = []
    ys: list[float | None] = []
    covered_to: int | None = None
    for start, end, score in zip(ordered["start"], ordered["end"], ordered["score"]):
        if covered_to is not None and start > covered_to:
            xs.append(None)
            ys.append(None)
        xs.append(int(start))
        ys.append(float(score))
        covered_to = int(end) if covered_to is None else max(covered_to, int(end))
    return xs, ys


def plot_guides(
    guides: Iterable[RfxGuide],
    psps: Iterable[PspGuide],
    length: int,
    output: str,
    title: str | None = None,
    inline_plotlyjs: bool = True,
) -> str | None:
    """Score against position for every guide, with the footprints marked.

    *length* is the full length of the plotted sequence, so uncovered stretches
    read as genuinely uncovered rather than being cropped away by autoscaling.
    """
    guides = sorted(guides, key=lambda g: g.start)
    if not guides:
        print("skipping: no guides to plot", file=sys.stderr)
        return None

    psps = sorted(psps, key=lambda p: p.start)
    transcript = guides[0].transcript

    df = pd.DataFrame(
        [
            {
                "start": g.start,
                "end": g.end,
                "score": g.score,
                "classification": g.classification,
                "spacer_seq": g.spacer_seq,
                "target_seq": g.target_seq,
                "psp": g.psp.name if g.psp else "-",
                "distance": g.distance if g.distance is not None else "-",
                # The same id as the intermediate table, so a point picked off
                # the figure can be looked up in the file by name.
                "guide_id": guide_id(
                    g.transcript, g.start, g.end, pad_width(length)
                ),
            }
            for g in guides
        ]
    )

    # One lane above the scatter for the footprints, sharing its x-axis so a
    # footprint sits directly over the guides that target it. Without a
    # footprint there is nothing to put in a lane, so the figure is one panel.
    has_lane = bool(psps)
    scatter_row = 2 if has_lane else 1

    fig = make_subplots(
        rows=scatter_row,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0,       # domains are set explicitly below
    )

    # The profile, not coloured by class: colouring it would break it into
    # disconnected segments, and it is the points that carry the class.
    profile_x, profile_y = profile_with_gaps(df)
    fig.add_trace(go.Scattergl(
        x=profile_x, y=profile_y,
        mode="lines",
        line=dict(color=LINE, width=1),
        connectgaps=False,        # the Nones above are the gaps; do not bridge them
        hoverinfo="skip",         # the markers' own hover covers every point
        showlegend=False,
        name="score profile",
    ), row=scatter_row, col=1)

    # One trace per class, in CLASS_ORDER rather than the order the classes
    # happen to appear, so the legend reads the same on every plot and the
    # overlaps land on top of the grey. A class with no guides is left out, so
    # a plot with no fully-overlapping guide says so by omission.
    for cls in CLASS_ORDER:
        points = df[df["classification"] == cls]
        if points.empty:
            continue
        label, fill, size = CLASS_STYLE[cls]
        fig.add_trace(go.Scattergl(
            x=points["start"], y=points["score"],
            mode="markers",
            name=f"{label} ({len(points):,})",
            marker=dict(
                color=fill,
                symbol="circle",
                size=size,
                line=dict(color=SURFACE, width=0.6 if cls == "none" else 1.4),
            ),
            customdata=points[
                ["end", "spacer_seq", "target_seq", "psp", "distance", "guide_id"]
            ],
            hovertemplate=(
                # The id first: it already states the span, and it is what
                # identifies the guide in the intermediate table.
                "<b>%{customdata[5]}</b><br>"
                f"{label}<br>"
                # Position relative to the plotted sequence, not the genome.
                f"{transcript} %{{x:,d}}-%{{customdata[0]:,d}}<br>"
                "score %{y:.4f}<br>"
                # Ends marked: every sequence in this pipeline is 5' -> 3'.
                "guide 5'-%{customdata[1]}-3'<br>"
                "target 5'-%{customdata[2]}-3'<br>"
                "PspCas13b %{customdata[3]}, gap %{customdata[4]} nt"
                "<extra></extra>"
            ),
        ), row=scatter_row, col=1)

    # The footprint bands, drawn **after** the traces on purpose.
    # `add_vrect` defaults to `exclude_empty_subplots=True`, so a band added to
    # a row that has no trace on it yet is silently dropped -- no error, no
    # shape, just a figure with the footprints missing. Adding them here means
    # the row is populated; `layer="below"` still puts them behind the data, so
    # insertion order costs nothing visually.
    for psp in psps:
        fig.add_vrect(
            x0=psp.start - 0.5, x1=psp.end + 0.5,   # whole bases, not centres
            fillcolor=PSP_BAND,
            opacity=0.13,
            line_width=0,
            layer="below",
            row=scatter_row, col=1,
            exclude_empty_subplots=False,   # belt and braces alongside the ordering
        )

    # The lane, drawn as bars on a one-category axis so each footprint keeps
    # its own hover: a shape would be silent, and the guide's name is most of
    # the point of having the lane.
    if has_lane:
        fig.add_trace(go.Bar(
            x=[p.length for p in psps],
            # `base` is where each bar starts; x is its width. Half a base is
            # subtracted so a bar covers whole bases rather than centres, as
            # the bands do.
            base=[p.start - 0.5 for p in psps],
            y=["PspCas13b"] * len(psps),
            orientation="h",
            marker=dict(color=PSP_LANE_FILL, line=dict(color=SURFACE, width=1)),
            showlegend=False,
            customdata=[(p.name, p.start, p.end, p.length) for p in psps],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "%{customdata[1]:,d}-%{customdata[2]:,d}<br>"
                "%{customdata[3]:,d} nt<extra></extra>"
            ),
        ), row=1, col=1)

        # Each footprint named above its bar. Eight short names across 1.2 kb
        # fit; the hover carries the coordinates.
        for psp in psps:
            fig.add_annotation(
                text=psp.name,
                x=(psp.start + psp.end) / 2, xref="x",
                y=1.0, yref="paper", yanchor="bottom",
                showarrow=False,
                font=dict(color=INK, size=9),
            )

    fig.update_layout(
        # Title only, no subtitle: the per-class counts are already in the
        # legend labels, so a subheading would state them twice.
        title=title or f"RfxCas13d predicted guides along {transcript}",
        # The lane bar is the full height of its lane. Without this it would be
        # drawn at a bar chart's default 80%, leaving a gap that reads as if
        # the annotation stopped short.
        bargap=0,
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK),
        hovermode="closest",
        legend=dict(title_text="Overlap with PspCas13b"),
    )

    # Sequence coordinates, so show them in full: Plotly's default SI
    # formatting turns 1238 into "1.2k", which is useless for finding a guide.
    #
    # The range is applied to every row: the lane shares the scatter's x-axis,
    # but `shared_xaxes` links only the *interaction*, so without an explicit
    # range per row the lane could autoscale to its own features and a
    # footprint would then sit over the wrong guides.
    fig.update_xaxes(
        tickformat=",d",
        hoverformat=",d",
        range=[0.5, length + 0.5],
    )
    fig.update_xaxes(gridcolor=GRID, row=scatter_row, col=1)
    fig.update_xaxes(
        title_text=f"Position in {transcript} (nt)", row=scatter_row, col=1
    )

    fig.update_yaxes(
        gridcolor=GRID, zeroline=False,
        title_text="Predicted guide score",
        row=scatter_row, col=1,
    )
    if has_lane:
        # No grid and no zero line on the lane: a gridline across a
        # one-category axis reads as a feature boundary.
        fig.update_yaxes(
            showgrid=False, zeroline=False, ticks="",
            domain=[1.0 - LANE_HEIGHT, 1.0],
            row=1, col=1,
        )
        fig.update_yaxes(domain=[0.0, 1.0 - LANE_HEIGHT - LANE_GAP],
                         row=scatter_row, col=1)

    # Inlining plotly.js costs ~3 MB but opens without a network connection,
    # which is what an archived or emailed plot needs; "cdn" keeps it small.
    directory = os.path.dirname(os.path.abspath(output))
    os.makedirs(directory, exist_ok=True)
    fig.write_html(output, include_plotlyjs=True if inline_plotlyjs else "cdn")
    return output
