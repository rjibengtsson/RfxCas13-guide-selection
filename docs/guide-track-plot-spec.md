> "Write a Python script that produces the plot described in docs/guide_track_plot_spec.md
> using the configuration in it."

# Plot spec: stacked guide-position tracks per gene
---

## 1. Configuration (edit these)

```python
# ---------- Input ----------
INPUT_CSV    = ""  # one row per (upper guide, lower guide) pair
GENBANK_FILE = ""     # source of each gene's length for the x-axis
GENE_COL     = "gene"                                     # column that splits the data into subplots
GENE_ORDER   = None     # e.g. ["5'UTR","orf1ab","S","E","M","N","3'UTR"]; None = order of first appearance
OUTPUT_HTML  = ""

# ---------- Upper tracks (above the midline, y = 0, 1, 2, ...) ----------
UPPER_TRACK = {
    "label":     "PspCas13b",          # legend name + bold first line of the hover box
    "id_col":    "PspCas13b index",    # unique guide ID; rows de-duplicated on this
    "start_col": "qstart",
    "end_col":   "qend",
    "color":     "#ff7f0e",            # orange
    "hover_fields": [                  # (hover label, column(s), format); format None = raw value
        ("Position",   ("qstart", "qend"),      "range"),
        ("Guide seq",  "PspCas13b guide seq",   None),
        # ("Target seq", "PspCas13b target seq", None),
    ],
}

# ---------- Lower tracks (below the midline, y = -1, -2, -3, ...) ----------
LOWER_TRACK = {
    "label":     "RfxCas13d",
    "id_col":    "RfxCas13d index",
    "start_col": "sstart",
    "end_col":   "send",
    "color":     "#1f77b4",            # blue
    "hover_fields": [
        ("Guide idx",   "RfxCas13d index",      None),
        ("Position",    ("sstart", "send"),     "range"),
        ("Guide seq",   "RfxCas13d guide seq",  None),
        ("Target seq",  "RfxCas13d target seq", None),
        ("Guide score", "guide score",          ".3f bold"),  # 3 decimals, wrapped in <b></b>
    ],
}

# ---------- Layout ----------
TITLE            = "Guide position overlap by subgenomic region"
LEGEND_TITLE     = "Cas type"
X_AXIS_TITLE     = "Position (nt)"   # bottom subplot only
LINE_WIDTH       = 10                # px thickness of each guide bar
Y_PAD            = 2                 # empty track-units above the top and below the bottom track
VERTICAL_SPACING = 0.05
HEIGHT_PER_GENE  = 320               # figure height = HEIGHT_PER_GENE * n_genes
MARGIN           = dict(l=30, r=30, t=60, b=20)
```

## 2. Input data

Each row pairs one upper-track guide with one nearby lower-track guide (output of
`bedtools closest` / `window`, filtered to the top 10 per gene by distance then guide score).
Columns in the original data:

`gene, PspCas13b index, qstart, qend, PspCas13b guide seq, PspCas13b target seq,
RfxCas13d index, sstart, send, RfxCas13d guide seq, RfxCas13d target seq, guide score, distance`

Guides repeat across pairs, so **de-duplicate each track separately per gene** on its `id_col`
(keep the first row) before drawing. Select columns by name, not by position.

Gene lengths: parse `GENBANK_FILE` (existing helper `entrez.get_feature_locations` in
`src/entrez.py`), then length = end − start for each feature whose name matches a `GENE_COL` value.

## 3. Figure structure

- Plotly (`plotly.graph_objects` + `make_subplots`), written out with `fig.write_html(OUTPUT_HTML)`.
- **One subplot row per gene**, single column, `shared_xaxes=False`; subplot title = gene name
  (font size 16, centred above the row).
- **Row heights scale with content**: per gene, `n_tracks = n_upper + n_lower + 1` (the +1 is a
  separator gap); `row_heights = [n_tracks_i / sum(n_tracks)]`.
- Each subplot's **x-axis** runs `0` → that gene's length (coordinates are local to the gene).
- Each subplot's **y-axis**: no tick labels, no title, range
  `[-(n_lower) - 0.5 - Y_PAD, n_upper - 0.5 + Y_PAD]`.
- X-axis title on the last row only.
- Figure: `title=TITLE`, `height=HEIGHT_PER_GENE * n_genes`, `autosize=True`, `margin=MARGIN`,
  `legend=dict(title=LEGEND_TITLE, itemsizing="constant")`, default Plotly template.

## 4. Track packing (interval stacking)

Pack each track's guides independently, greedy first-fit:

1. Sort guides by `start_col`.
2. Place each guide on the lowest-numbered track whose last guide ends at or before this guide's
   start (`start >= track_end`); otherwise open a new track.
3. Track index `t` (0-based) → y position:
   - upper track: `y = t` (0, 1, 2, … upward)
   - lower track: `y = -(t + 1)` (−1, −2, … downward)

Overlapping guides therefore stack away from the midline, non-overlapping guides share a row.

## 5. Guide bars (traces)

- Each guide is one `go.Scatter(mode="lines")`: `x=[start, end]`, `y=[y, y]`,
  `line=dict(width=LINE_WIDTH, color=<track color>)`.
- `name` and `legendgroup` = the track's `label`. Set `showlegend=True` for only the **first**
  trace of each track type in the whole figure, so the legend has exactly two entries and
  clicking one toggles that track across all subplots.

## 6. Hover box

Build each trace's `hovertemplate` from its track's `hover_fields`, filling values in Python
(f-string per trace) rather than via `customdata`:

```
<b>{label}</b><br>
{field label}: {value}<br>     # one line per hover_fields entry, in order
<extra></extra>                # suppresses the trace-name side box
```

Format codes:
- `"range"` → `f"{start} - {end}"` for a `(start_col, end_col)` tuple.
- `".3f bold"` → `f"<b>{float(v):.3f}</b>"`.
- `None` → the raw cell value.

So an upper-track hover box shows three lines (bold label, position range, guide sequence) and a
lower-track hover box shows six (bold label, guide index, position range, guide sequence, target
sequence, bold 3-decimal score). Adding or removing a hover line should only mean editing
`hover_fields`.

## 7. Notes / known quirks to fix in the regenerated script

- The original legend label for the lower track was misspelled `RxfCas13d` while the hover box
  said `RfxCas13d`; use one spelling from `UPPER_TRACK["label"]` / `LOWER_TRACK["label"]`.
- The original sliced columns positionally (`iloc[:, 0:6]`, `np.r_[0, 6:12]`), which breaks if
  column order changes — select by the configured column names instead.
- The original hard-coded two Cas types; prefer a list of track dicts (upper/lower, or N tracks)
  driven by the configuration block.
