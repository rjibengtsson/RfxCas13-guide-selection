# Scatterline plot: specification

A description of a reusable plot type, written so it can be handed to a coding
agent as the brief for a new script. It is domain-neutral: the plot draws
*points along a one-dimensional coordinate*, and knows nothing about the data
that produced them. The worked example at the end fills every variable in with
the values used for the RfxCas13d guide figure in this project.

Reference implementation: `src/visualise_scatterlines.py` (library) and
`scripts/plot_guides.py` (CLI). Output is a self-contained Plotly HTML file.

---

## 1. Figure anatomy

Reading from the bottom up, the figure has four layers:

1. **Main scatter panel.** One point per record, placed at its coordinate on
   the x-axis and its measured value on the y-axis. Points are coloured by a
   categorical field.
2. **Connecting line.** A single recessive line joining the points in
   coordinate order, drawn *under* them. It is one neutral colour and is never
   coloured by category — colouring it would break it into disconnected
   segments and duplicate information the points already carry. The line is
   broken wherever the records stop covering the coordinate axis (§6.2).
3. **Interval bands.** Shaded vertical bands behind the scatter, marking
   regions of interest on the coordinate axis. Low opacity, no border, drawn
   below the data.
4. **Upper tracks (lanes).** One or more thin horizontal lanes stacked above
   the main panel, sharing its x-axis. Each lane draws a set of intervals as
   horizontal bars, so an interval sits directly above the points that fall
   inside it. Lanes are positional references, not the subject of the figure,
   so they are thin.

Bands and lanes usually draw the *same* intervals: the band gives context
behind every point, the lane gives the interval a label and its own hover.
Either can be used without the other.

A figure with no intervals has no lane and no bands, and collapses to a single
panel. The script must handle that case rather than emitting an empty lane.

---

## 2. Variables

Everything below is a parameter of the plotting function or a module-level
constant. Nothing in this list should be hard-coded inside the drawing code.

### 2.1 Colour

| Variable | Type | Purpose |
|---|---|---|
| `SURFACE` | hex colour | Plot and paper background. Also the marker outline colour, so adjacent points stay countable where they touch. |
| `GRID` | hex colour | Gridlines on both axes of the main panel. |
| `INK` | hex colour | Axis titles, tick text, lane annotations. |
| `LINE` | hex colour | The connecting line. Deliberately low-contrast. |
| `BAND_FILL` | hex colour | Interval bands behind the scatter. |
| `BAND_OPACITY` | float | Band opacity. Low enough that points read clearly through it. |
| `LANE_FILL` | hex colour | Bars in the upper track. |
| `CATEGORY_STYLE` | ordered map | The core colour variable — see below. |

`CATEGORY_STYLE` maps each category value to a display label, a fill colour, a
marker size and an outline width:

```
CATEGORY_STYLE = {
    "<category value>": {
        "label":   "<legend text>",
        "fill":    "<hex colour>",
        "size":    <marker size in px>,
        "outline": <outline width in px>,
    },
    ...
}
```

Rules that matter:

- **Order is draw order.** The map is iterated in order, so put the common,
  uninteresting category first and the rare, interesting ones last — later
  traces land on top. The legend then reads identically on every figure
  regardless of which categories happen to be present.
- **A category with no records is omitted**, not drawn empty. Its absence from
  the legend is itself informative.
- **Size encodes importance, not magnitude.** Give the categories you are
  looking for larger markers and a thicker outline; a fill change alone is easy
  to lose among a thousand points.
- **One variable, one channel.** The y-axis already carries the measured value,
  so do not also encode it in colour. Colour carries the category, which is the
  question the figure exists to answer.
- Pick categorical hues from a palette already validated for colour-vision
  deficiency rather than by eye. Reserve the band/lane hue for intervals so a
  band never reads as another category of point.

### 2.2 X axis (the coordinate axis)

| Variable | Type | Purpose |
|---|---|---|
| `x_column` | column name | Record's position on the coordinate axis. |
| `x_title` | string | Axis title, typically templated with the sequence/series name and unit. |
| `x_tickformat` | Plotly format | Tick text format. |
| `x_hoverformat` | Plotly format | Hover text format for the coordinate. |
| `axis_length` | int | Full extent of the coordinate space. The x range is pinned to `[0.5, axis_length + 0.5]` so uncovered stretches stay visible instead of being cropped away by autoscaling. |

### 2.3 Y axis (the measured value)

| Variable | Type | Purpose |
|---|---|---|
| `y_column` | column name | Record's measured value. |
| `y_title` | string | Axis title. |
| `y_zeroline` | bool | Whether to draw the zero line. Usually off. |
| `y_range` | pair or None | Explicit limits; None autoscales. |

### 2.4 Hover box

A **list of field specs**, not a hard-coded template string, so fields can be
added, reordered or dropped without rewriting Plotly syntax:

```
HOVER_FIELDS = [
    {"label": "<prefix text or ''>", "column": "<column>", "format": "<fmt>", "bold": <bool>},
    ...
]
```

The script builds `customdata` from the named columns and assembles the
`hovertemplate` from the list, in order. Conventions:

- First row is the record's **identifier**, in bold, on its own line.
- Second row is the **category label** for the trace being drawn.
- Then position, then measured value, then any payload fields.
- Suppress the trace-name box (`<extra></extra>`).
- A field whose value is missing renders as a placeholder (`-`), never as
  `None` or `nan`.

### 2.5 Upper tracks

A **list of track definitions**, so one lane generalises to several:

```
TRACKS = [
    {
        "name":        "<lane y-axis category label>",
        "intervals":   [(id, start, end), ...],   # 1-based inclusive
        "fill":        "<hex colour>",
        "annotate":    <bool>,    # print each interval's id above its bar
        "hover_fields": [...],    # as in §2.4, for the bars
    },
    ...
]
```

| Variable | Type | Purpose |
|---|---|---|
| `LANE_HEIGHT` | float | Height of one lane as a fraction of figure height. Thin. |
| `LANE_GAP` | float | Gap between the lane stack and the main panel, same units. |
| `annotate_fontsize` | int | Size of the per-interval labels above the bars. |

Annotation only works when the labels are short and the intervals are few
enough not to collide; the hover carries the coordinates either way, so turn
annotation off when the lane gets crowded.

### 2.6 Output

| Variable | Type | Purpose |
|---|---|---|
| `output` | path | HTML file to write. Parent directories are created. |
| `title` | string or None | Figure title; a default is generated from the series name. |
| `inline_plotlyjs` | bool | True inlines plotly.js (~3 MB, opens with no network — right for an archived or emailed figure); False links a CDN copy. |

---

## 3. Data contract

The plotting function takes a tidy table, one row per point, with at least:

- an **identifier** column — stable, and ideally encoding the position so a
  point picked off the figure stays identifiable once copied into a
  spreadsheet;
- the **coordinate** columns `start` and `end` (1-based, inclusive, `start <=
  end`), where `start` is the plotted x position and `end` is used for the span
  in the hover and for the gap logic;
- the **measured value** column;
- the **category** column, whose values are the keys of `CATEGORY_STYLE`;
- any further columns named by `HOVER_FIELDS`.

Intervals for bands and lanes arrive separately as `(id, start, end)` triples
in the same coordinate space.

Nothing is filtered inside the plotting function. Thresholding and selection
are upstream decisions; the figure draws every row it is given, so that what is
absent from the figure is absent from the data.

---

## 4. Layout

- Two rows via `make_subplots(shared_xaxes=True, vertical_spacing=0)`, with the
  lane as row 1 and the scatter as row 2. Domains are set explicitly rather
  than via `vertical_spacing`: lane `[1 - LANE_HEIGHT, 1]`, scatter
  `[0, 1 - LANE_HEIGHT - LANE_GAP]`.
- With no intervals, build a single-row figure and put the scatter in row 1.
- `bargap=0` so a lane bar fills the full height of its lane; the default 80%
  leaves a gap that reads as if the bar stopped short.
- Legend title names the categorical variable.
- Per-category counts belong in the legend labels (`"No overlap (1,180)"`), not
  in a subtitle — a subtitle would state them twice.

---

## 5. Interaction

- `hovermode="closest"`.
- The connecting line sets `hoverinfo="skip"` and `showlegend=False`: the
  markers' own hover covers every point, and a second hover target on the same
  coordinate is noise.
- Lane intervals are drawn as **bars, not shapes**, specifically so each keeps
  its own hover. A shape is silent, and the interval's identity is most of the
  reason the lane exists.

---

## 6. Behaviours that are easy to get wrong

These are the non-obvious rules in the reference implementation. Each one
produced a plausible-looking but wrong figure before it was fixed, so carry
them over.

### 6.1 Half-unit offsets

Intervals cover whole units, not centres. A band or bar spanning positions
`a..b` is drawn from `a - 0.5` to `b + 0.5`. For a lane bar this means
`base = start - 0.5` and `width = end - start + 1`. Without the offset every
interval is drawn half a unit short at each end and appears not to cover the
points at its edges.

### 6.2 Break the line at coverage gaps

Insert `(None, None)` into the line's x/y sequences wherever the records stop
covering the axis, and set `connectgaps=False`. Compare each record's `start`
against the **running maximum `end`** of all earlier records, not the previous
row's `end`, so a record nested inside an earlier, longer one still counts as
covered. Records that merely abut (`start == previous end + 1`) share no unit
and are also broken apart.

Without this the line interpolates a value across positions no record covers,
which is a fabricated reading. Over a fully tiled axis no break appears — and
being able to see that is itself the point.

### 6.3 Add bands after the traces

`add_vrect` defaults to `exclude_empty_subplots=True`, so a band added to a row
that carries no trace yet is **silently dropped** — no error, no shape, just a
figure with the bands missing. Add all traces first, then the bands, and pass
`exclude_empty_subplots=False` as well. `layer="below"` still puts them behind
the data, so the later insertion costs nothing visually.

### 6.4 Pin the x range on every row

`shared_xaxes=True` links the *interaction* between rows, not the ranges. Apply
the explicit range to all rows, or the lane autoscales to its own intervals and
every bar sits over the wrong points — a figure that renders perfectly and is
wrong.

### 6.5 Format coordinates in full

Set `tickformat=",d"` and `hoverformat=",d"`. Plotly's default SI formatting
turns 1238 into "1.2k", which is useless for locating a record.

### 6.6 Outline every marker in the background colour

At a dense tiling, adjacent points touch. A hairline outline in `SURFACE` is
what keeps them countable.

### 6.7 Validate positions against the source before drawing

If coordinates were computed in an earlier stage and read back from a file,
check that each record's stored payload still matches the source at its stated
position, and fail loudly if not. A table computed against a different version
of the source renders a perfectly convincing, entirely wrong figure. Likewise,
coordinates that come from a different reference frame need an explicit offset,
and the script should **print the offset and where it came from**, since
getting it wrong slides every interval along the axis without any visible
symptom.

---

## 7. CLI shape

The library parses no arguments and reads no hard-coded path; a thin CLI wraps
it. The CLI should:

- accept either a **precomputed table** (positions and categories already
  assigned upstream — preferred, because one place owns position and the figure
  cannot disagree with the table) or a **raw input** from which positions are
  derived on the fly as a convenience;
- take the coordinate source (e.g. the reference sequence) separately, and use
  it for the axis length and the validation in §6.7;
- take the intervals from a file, or recover them from the precomputed table
  when they are named there;
- print a short report before writing: the series name and length, the row
  count and its source, each interval with its span, and the per-category
  counts;
- write the figure and print the output path.

Per the project's output rules, the parameters of each run are recorded
alongside its results so the figure can be reproduced.

---

## 8. Worked example — RfxCas13d guide figure

The figure at `figures/SARS-CoV-2/MT007544.1_N_guides.html`, produced by
`scripts/plot_guides.py`. Points are predicted RfxCas13d guides along the
SARS-CoV-2 N gene; intervals are the footprints of previously tested
PspCas13b guides; the category is how each guide sits against those footprints.

**Colour**

```
SURFACE      = "#fcfcfb"
GRID         = "#e8e7e2"
INK          = "#52514e"
LINE         = "#d9d8d3"
BAND_FILL    = "#d6453f"   BAND_OPACITY = 0.13
LANE_FILL    = "#c2403b"

CATEGORY_STYLE = {
    "none":    {"label": "No overlap",      "fill": "#b0aea8", "size":  5, "outline": 0.6},
    "partial": {"label": "Partial overlap", "fill": "#eb6834", "size":  9, "outline": 1.4},
    "full":    {"label": "Full overlap",    "fill": "#2a78d6", "size": 11, "outline": 1.4},
}
```

Grey is the default and is deliberately recessive: most guides touch no
footprint and are context, not subject. Orange and blue are slots 1 and 2 of
the data-viz reference palette, taken as-is, and are ordered so the stronger
relationship gets the stronger colour. Red is the one hue that is not a
category, so a band reads as a footprint rather than as another class of point.

**Axes**

```
x_column = "start"      x_title = "Position in {transcript} (nt)"
x_tickformat = ",d"     x_hoverformat = ",d"
axis_length = len(sequence)          # full gene length, so uncovered stretches show

y_column = "score"      y_title = "Predicted guide score"
y_zeroline = False      y_range = None
```

**Hover fields**

```
HOVER_FIELDS = [
    {"label": "",            "column": "guide_id",    "bold": True},
    {"label": "",            "column": "<category label>"},
    {"label": "{transcript} ", "column": "start-end",  "format": ",d"},
    {"label": "score ",      "column": "score",       "format": ".4f"},
    {"label": "guide 5'-",   "column": "spacer_seq",  "suffix": "-3'"},
    {"label": "target 5'-",  "column": "target_seq",  "suffix": "-3'"},
    {"label": "PspCas13b ",  "column": "psp"},
    {"label": ", gap ",      "column": "distance",    "suffix": " nt"},
]
```

The id leads because it already states the span and is what identifies the
guide in the intermediate table. Strand ends are marked because every sequence
in this pipeline is 5' → 3'. Position is relative to the plotted gene, not the
genome; the genome coordinates live in the table, not the hover.

**Upper track**

```
TRACKS = [{
    "name": "PspCas13b",
    "intervals": [(name, start, end) for each blastn-located PspCas13b guide],
    "fill": LANE_FILL,
    "annotate": True,        # eight short names across 1.2 kb fit
    "hover_fields": [id (bold), "{start}-{end}", "{length} nt"],
}]
LANE_HEIGHT = 0.03
LANE_GAP    = 0.03
```

**Categories, as assigned upstream**

`full` = the 23 nt guide lies entirely inside a 30 nt footprint (containment is
the only way the shorter can be fully covered); `partial` = shares at least one
base but hangs over an edge; `none` = shares no base with any footprint. Each
guide is classified once, by its strongest relationship to any footprint, so
there is one record per guide rather than one per (guide, footprint) pair.

**Caveat on reuse.** The category definitions above are specific to a 23 nt
guide against a 30 nt footprint. For guides and footprints of other lengths,
"full" may no longer mean containment, and the classification rule has to be
restated before the figure means anything.
