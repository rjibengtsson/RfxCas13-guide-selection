#!/usr/bin/env python
"""Draw selected PspCas13b and RfxCas13d guides as stacked position tracks.

Implements ``docs/guide-track-plot-spec.md`` for the SARS-CoV-2 N gene. The
figure answers one question: where do the two Cas types' guides sit relative to
each other along the gene. PspCas13b guides are packed above the midline and
RfxCas13d below it; a bar's hover carries its identity and sequence.

Two differences from the spec's own configuration, both forced by the inputs:

* the spec expects one ``closest_guides.csv`` pairing each upper guide with a
  nearby lower one, so it de-duplicates per track. Here the two tracks come
  from separate tables and each guide already appears once, so there is nothing
  to de-duplicate;
* the spec cites ``src/entrez.py`` for gene lengths. That module is not in this
  repository, so the GenBank file is parsed with biopython instead
  (``visualise_tracks.gene_span``).

PspCas13b positions arrive in **genome** coordinates from the blastn table and
are shifted into the gene's frame; the RfxCas13d located table is already in
it. The shift is printed, because getting it wrong slides every PspCas13b bar
along the axis and the figure still looks right.

Alongside the figure it writes the selected RfxCas13d guides as a CSV, built
from the same selection, so the table and the figure cannot disagree about
what was chosen. The table's columns are named by ``RFX_EXPORT`` and spell out
which molecule each sequence is -- ``target_dna``, ``spacer_dna``,
``spacer_rna`` -- while the hover box keeps the located table's own names.

Usage
-----
    python scripts/plot_guide_tracks.py
    python scripts/plot_guide_tracks.py --output figures/other.html --cdn-plotlyjs
    python scripts/plot_guide_tracks.py --no-csv        # figure only

Selections and paths live in the CONFIG block below; edit them there rather
than passing them in, so a figure's inputs are recorded in the file that drew
it. Needs pandas, plotly and biopython (the `pik3-guide` conda env).
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import pandas as pd  # noqa: E402

from table_io import read_table  # noqa: E402
from visualise_tracks import (  # noqa: E402
    Interval,
    Track,
    gene_span,
    plot_tracks,
    shift_to_gene,
)


# ============================ CONFIG ========================================

GENBANK_FILE = "data/SARS-CoV-2/MT007544.1/MT007544.1.gbk"
GENE = "N"
OUTPUT_HTML = "figures/SARS-CoV-2/MT007544.1_N_guide_tracks.html"
# The selected RfxCas13d guides as a table, carrying the same columns the
# figure's hover box shows. Written from the same selection that is drawn, so
# the table and the figure cannot disagree about what was chosen.
OUTPUT_RFX_CSV = "results/SARS-CoV-2/rfx_guides_MT007544.1_N_selected.csv"

# ---------- Upper track: PspCas13b ----------
# Sequences and positions live in different files: the TSV names the guides and
# gives their spacers, the blastn table says where each one landed.
PSP_GUIDES = "data/pspcas13b_guides/panSarbeCov_guides.tsv"
PSP_HITS = "results/SARS-CoV-2/panSarbeCov_guides_MT007544.1_fmt6.tsv"
PSP_SELECTION = (
    "Cand_3_N_Sarbecov",
    "N_cand_93",
    "NCP1",
    "Cand_5_N_Sarbecov",
    "A8_NCP",
    "Cand_1_N_Sarbecov",
)
# (hover label, source) per line. "start"/"end" are the gene-local coordinates
# computed here; the rest are columns of PSP_GUIDES.
PSP_HOVER = (
    ("guide_id", "Name"),
    ("start", "start"),
    ("end", "end"),
    ("spacer", "Spacer"),
)
PSP_COLOR = "#ff7f0e"      # orange
# Print each PspCas13b guide's name above its bar. Only six, and the names are
# short, so they fit without colliding; the RfxCas13d ids are long and stay in
# the hover.
PSP_ANNOTATE = True

# ---------- Lower track: RfxCas13d ----------
# Already carries gene-local start/end, written by scripts/locate_guides.py.
RFX_LOCATED = "results/SARS-CoV-2/rfx_guides_MT007544.1_N_located.csv"
RFX_SELECTION = (
    "MT007544.1_N_1001-1023",
    "MT007544.1_N_1074-1096",
    "MT007544.1_N_1185-1207",
    "MT007544.1_N_1201-1223",
    "MT007544.1_N_0161-0183",
    "MT007544.1_N_0160-0182",
    "MT007544.1_N_0394-0416",
    "MT007544.1_N_0070-0092",
    "MT007544.1_N_0463-0485",
)
# Shorten the accession to a Cas-type prefix in the figure only: on the plot
# the accession is already in the title and the gene in the subplot heading, so
# it costs a reader nothing and leaves the long ids legible at 10 px.
#
# The written table keeps the full ids (see RFX_EXPORT), because those are what
# join back to rfx_guides_MT007544.1_N_located.csv. The two therefore disagree
# by design, which is why this rewrite happens here, at the hover, and not in
# the frame both outputs are built from.
RFX_ID_REPLACE = ("MT007544.1_", "rfx_")

RFX_HOVER = (
    ("guide_id", "guide_id"),
    ("start", "start"),
    ("end", "end"),
    ("target_sequence", "target_sequence"),
    ("guide_sequence", "guide_sequence"),
    ("guide_score", "guide_score"),
)
RFX_COLOR = "#1f77b4"      # blue

# Columns of the exported table, as (output name, source column). The names
# differ from the hover box's on purpose: the figure labels fields the way the
# located table does, while the table spells out which molecule each sequence
# is, because that is what gets read when a guide is ordered.
#
# spacer_rna is derived, not read -- see transcribe().
RFX_EXPORT = (
    ("guide_id", "guide_id"),
    ("start", "start"),
    ("end", "end"),
    ("target_dna", "target_sequence"),
    ("spacer_dna", "guide_sequence"),
    ("spacer_rna", None),              # derived from guide_sequence
    ("guide_score", "guide_score"),
)

# A score needs rounding for display; everything else is shown as it is stored.
HOVER_FORMATS = {"guide_score": "{:.4f}"}

# ---------- Layout ----------
TITLE = "Selected PspCas13b and RfxCas13d guides along the SARS-CoV-2 N gene"
LEGEND_TITLE = "Cas type"
X_AXIS_TITLE = "Position in N (nt)"
LINE_WIDTH = 10
Y_PAD = 2
ANNOTATION_FONT_SIZE = 10
SUBPLOT_TITLE_FONT_SIZE = 16
VERTICAL_SPACING = 0.05
HEIGHT_PER_GENE = 320
MARGIN = dict(l=30, r=30, t=60, b=20)

# ============================ END CONFIG ====================================


def format_value(column: str, value) -> str:
    """One hover value as text, rounded if HOVER_FORMATS says so."""
    template = HOVER_FORMATS.get(column)
    if template is None:
        return str(value)
    try:
        return template.format(float(value))
    except (TypeError, ValueError):
        # A malformed cell is shown raw rather than crashing the figure; it is
        # visible in the hover, which is where it would be noticed.
        return str(value)


def require_all(found: set[str], wanted: tuple[str, ...], source: str) -> None:
    """Raise unless every wanted name was found.

    A missing selection is an error, not a smaller figure: the likely cause is
    a typo, and the alternative is a plot quietly short a guide the reader
    believes is on it.
    """
    missing = [name for name in wanted if name not in found]
    if missing:
        raise ValueError(
            f"{len(missing)} selected name(s) not in {source}: "
            + ", ".join(map(repr, missing))
        )


def build_psp_intervals(gene_start: int, gene_length: int) -> list[Interval]:
    """The selected PspCas13b guides, shifted into the gene's coordinates."""
    spacers = read_table(PSP_GUIDES)
    hits = read_table(PSP_HITS)

    by_name = dict(zip(spacers["Name"].astype(str), spacers["Spacer"].astype(str)))
    require_all(set(by_name), PSP_SELECTION, PSP_GUIDES)

    # One row per guide: keep the best HSP by bitscore, as build_psp_guides does
    # for the scatter figure, so the two figures place a guide identically.
    hits = hits.loc[hits.groupby("qseqid")["bitscore"].idxmax()]
    hit_rows = {str(row["qseqid"]): row for _, row in hits.iterrows()}
    require_all(set(hit_rows), PSP_SELECTION, PSP_HITS)

    intervals = []
    for name in PSP_SELECTION:
        row = hit_rows[name]
        start, end = shift_to_gene(
            int(row["sstart"]), int(row["send"]), gene_start
        )
        if end < 1 or start > gene_length:
            raise ValueError(
                f"{name} maps to {start}-{end}, outside {GENE} (1-{gene_length}); "
                "check the GenBank gene and the blastn subject agree"
            )
        values = {"Name": name, "start": start, "end": end, "Spacer": by_name[name]}
        intervals.append(
            Interval(
                name=name,
                start=start,
                end=end,
                fields=tuple(
                    (label, format_value(column, values[column]))
                    for label, column in PSP_HOVER
                ),
            )
        )
    return intervals


def transcribe(dna: str) -> str:
    """A DNA sequence as the RNA of the same strand: T -> U, nothing else.

    The spacer is already the reverse complement of its target, so this is a
    transcription of the strand as written, **not** a reverse complement --
    doing both would give the target's RNA, not the spacer's, and the mistake
    would not be visible by inspection.

    Raises on any character that is not an unambiguous DNA base, because a
    degenerate or masked base silently carried into an ordered oligo is worse
    than a failed run. 5' -> 3' in, 5' -> 3' out.
    """
    sequence = str(dna).strip().upper()
    unexpected = sorted(set(sequence) - set("ACGT"))
    if unexpected:
        raise ValueError(
            f"{dna!r} is not unambiguous DNA; unexpected character(s): "
            + ", ".join(map(repr, unexpected))
        )
    return sequence.replace("T", "U")


def export_rfx(frame):
    """The selected guides as the table to write: renamed and with the RNA.

    Built from *frame*, the same rows the figure draws, so the table cannot
    name a different selection than the figure shows.
    """
    out = {}
    for name, source in RFX_EXPORT:
        if name == "spacer_rna":
            out[name] = frame["guide_sequence"].map(transcribe)
        else:
            out[name] = frame[source]
    return pd.DataFrame(out)


def select_rfx():
    """The selected RfxCas13d rows from the located table, ordered by position.

    Returns a frame holding the columns :data:`RFX_HOVER` names, in that order.
    It feeds both the bars' hover text and, via :func:`export_rfx`, the written
    table, so one selection drives both and neither can name a guide the other
    does not.

    Sorted by start so the table reads along the gene, matching how the bars
    are laid out.
    """
    located = read_table(RFX_LOCATED)
    needed = {column for _, column in RFX_HOVER} | {"guide_id", "start", "end"}
    missing = sorted(needed - set(located.columns))
    if missing:
        raise ValueError(
            f"{RFX_LOCATED} is missing column(s) {', '.join(map(repr, missing))}; "
            "was it written by scripts/locate_guides.py?"
        )

    require_all(set(located["guide_id"].astype(str)), RFX_SELECTION, RFX_LOCATED)

    selected = located[located["guide_id"].astype(str).isin(RFX_SELECTION)]
    columns = [column for _, column in RFX_HOVER]
    return selected.sort_values(["start", "end"])[columns].reset_index(drop=True)


def display_id(guide_id: str) -> str:
    """A guide id as the figure shows it, abbreviated per RFX_ID_REPLACE.

    Only the leading accession is rewritten, and only if it is there, so an id
    that does not carry the expected prefix is shown unchanged rather than
    mangled. The full id stays in the written table.
    """
    prefix, replacement = RFX_ID_REPLACE
    guide_id = str(guide_id)
    if not prefix or not guide_id.startswith(prefix):
        return guide_id
    return replacement + guide_id[len(prefix):]


def build_rfx_intervals(frame) -> list[Interval]:
    """One bar per selected RfxCas13d guide, hover text included.

    The hover's id is abbreviated by :func:`display_id`; ``Interval.name``
    keeps the full id, so anything keyed on identity -- a label, a lookup back
    into the located table -- still uses the real one.
    """
    intervals = []
    for _, row in frame.iterrows():
        values = {column: row[column] for _, column in RFX_HOVER}
        values["guide_id"] = display_id(row["guide_id"])
        intervals.append(
            Interval(
                name=str(row["guide_id"]),
                start=int(row["start"]),
                end=int(row["end"]),
                fields=tuple(
                    (label, format_value(column, values[column]))
                    for label, column in RFX_HOVER
                ),
            )
        )
    return intervals


def write_params(output: str, extra: dict) -> str:
    """Record the run's inputs and selections beside the figure."""
    params_path = os.path.splitext(output)[0] + ".params.txt"
    with open(params_path, "w") as out:
        out.write(f"# {os.path.basename(output)} -- run parameters\n")
        out.write(
            f"generated: {datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        )
        out.write("script: scripts/plot_guide_tracks.py\n")
        out.write("spec: docs/guide-track-plot-spec.md\n")
        out.write(f"genbank: {GENBANK_FILE}\n")
        out.write(f"gene: {GENE}\n")
        out.write(f"psp_guides: {PSP_GUIDES}\n")
        out.write(f"psp_hits: {PSP_HITS}\n")
        out.write(f"rfx_located: {RFX_LOCATED}\n")
        for key, value in extra.items():
            out.write(f"{key}: {value}\n")
        out.write(f"psp_selection ({len(PSP_SELECTION)}):\n")
        for name in PSP_SELECTION:
            out.write(f"  - {name}\n")
        out.write(f"rfx_selection ({len(RFX_SELECTION)}):\n")
        for name in RFX_SELECTION:
            out.write(f"  - {name}\n")
    return params_path


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--output", default=OUTPUT_HTML, help=f"output HTML (default: {OUTPUT_HTML})")
    p.add_argument("--title", default=TITLE, help="plot title")
    p.add_argument(
        "--cdn-plotlyjs",
        action="store_true",
        help="link plotly.js from a CDN instead of inlining it: a ~3 MB smaller "
        "file that needs a network connection to open",
    )
    p.add_argument(
        "--rfx-csv",
        default=OUTPUT_RFX_CSV,
        help=f"table of the selected RfxCas13d guides, with the hover box's "
        f"columns (default: {OUTPUT_RFX_CSV})",
    )
    p.add_argument(
        "--no-csv", action="store_true", help="skip the RfxCas13d table"
    )
    p.add_argument("--no-params", action="store_true", help="skip the .params.txt sidecar")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    gene_start, gene_end = gene_span(GENBANK_FILE, GENE)
    gene_length = gene_end - gene_start + 1
    print(f"gene {GENE}: {gene_start:,}-{gene_end:,} ({gene_length:,} nt) from {GENBANK_FILE}")
    # Stated explicitly: this offset decides where every PspCas13b bar lands.
    print(f"shifting PspCas13b genome coordinates by -{gene_start - 1:,}")

    psp = build_psp_intervals(gene_start, gene_length)
    rfx_frame = select_rfx()
    rfx = build_rfx_intervals(rfx_frame)

    if not args.no_csv:
        os.makedirs(os.path.dirname(os.path.abspath(args.rfx_csv)), exist_ok=True)
        # Full stored precision, not the 4 dp the hover rounds to: the figure is
        # for reading and the table is for computing on.
        exported = export_rfx(rfx_frame)
        exported.to_csv(args.rfx_csv, index=False)
        print(
            f"wrote {len(exported)} selected RfxCas13d guides to {args.rfx_csv}\n"
            f"  columns: {', '.join(exported.columns)}"
        )

    print(f"PspCas13b: {len(psp)} guides")
    for interval in sorted(psp, key=lambda i: i.start):
        print(f"  {interval.name}: {interval.start:,}-{interval.end:,} ({interval.length} nt)")
    print(f"RfxCas13d: {len(rfx)} guides")
    for interval in sorted(rfx, key=lambda i: i.start):
        print(f"  {interval.name}: {interval.start:,}-{interval.end:,} ({interval.length} nt)")

    tracks = (
        # annotate: six short names across 1.2 kb fit above the bars. The
        # RfxCas13d names are long ids and there are nine of them, so those are
        # left to the hover.
        Track(
            label="PspCas13b", intervals=tuple(psp), color=PSP_COLOR,
            side="upper", annotate=PSP_ANNOTATE,
        ),
        Track(label="RfxCas13d", intervals=tuple(rfx), color=RFX_COLOR, side="lower"),
    )

    written = plot_tracks(
        genes=[(GENE, gene_length, tracks)],
        output=args.output,
        title=args.title,
        legend_title=LEGEND_TITLE,
        x_axis_title=X_AXIS_TITLE,
        line_width=LINE_WIDTH,
        y_pad=Y_PAD,
        vertical_spacing=VERTICAL_SPACING,
        height_per_gene=HEIGHT_PER_GENE,
        margin=MARGIN,
        annotation_font_size=ANNOTATION_FONT_SIZE,
        subplot_title_font_size=SUBPLOT_TITLE_FONT_SIZE,
        inline_plotlyjs=not args.cdn_plotlyjs,
    )
    print(f"wrote {written}")

    if not args.no_params:
        params_path = write_params(
            args.output,
            {
                "gene_start": gene_start,
                "gene_end": gene_end,
                "gene_length": gene_length,
                "coordinate_shift": -(gene_start - 1),
                "n_psp_drawn": len(psp),
                "n_rfx_drawn": len(rfx),
                "rfx_csv": "none" if args.no_csv else args.rfx_csv,
            },
        )
        print(f"  parameters: {params_path}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FileNotFoundError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
