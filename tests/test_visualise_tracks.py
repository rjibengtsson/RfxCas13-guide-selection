"""Tests for the stacked guide-track figure: packing, coordinates, hover."""

import os
import re

import pytest

from visualise_tracks import (
    Interval,
    Track,
    gene_span,
    hover_template,
    lane_to_y,
    pack_intervals,
    plot_tracks,
    shift_to_gene,
)


def interval(start, end, name=None):
    return Interval(name=name or f"{start}-{end}", start=start, end=end)


class TestShiftToGene:
    def test_the_gene_start_becomes_base_one(self):
        assert shift_to_gene(28274, 28303, 28274) == (1, 30)

    def test_a_known_footprint_lands_where_the_located_table_says(self):
        # N_cand_93 hits MT007544.1 at 29474-29503; the N gene starts at 28274,
        # and scripts/locate_guides.py independently placed it at 1201-1230.
        assert shift_to_gene(29474, 29503, 28274) == (1201, 1230)

    def test_a_minus_strand_hit_is_ordered(self):
        # blastn reports sstart > send on the minus strand; the extent is the
        # same and a bar only needs it ordered.
        assert shift_to_gene(29503, 29474, 28274) == (1201, 1230)

    def test_an_unshifted_span_is_unchanged(self):
        assert shift_to_gene(100, 129, 1) == (100, 129)


class TestPackIntervals:
    def test_non_overlapping_guides_share_a_lane(self):
        packed = pack_intervals([interval(1, 30), interval(100, 129)])
        assert [lane for _, lane in packed] == [0, 0]

    def test_overlapping_guides_stack(self):
        packed = pack_intervals([interval(1, 30), interval(20, 49)])
        assert [lane for _, lane in packed] == [0, 1]

    def test_touching_guides_share_a_lane(self):
        # 31 begins the base after 30 ends: adjacent, not overlapping.
        packed = pack_intervals([interval(1, 30), interval(31, 60)])
        assert [lane for _, lane in packed] == [0, 0]

    def test_a_one_base_overlap_stacks(self):
        packed = pack_intervals([interval(1, 30), interval(30, 59)])
        assert [lane for _, lane in packed] == [0, 1]

    def test_the_lowest_free_lane_is_reused(self):
        # Three mutually overlapping, then one clear of the first: it must drop
        # back to lane 0 rather than opening a fourth lane.
        packed = pack_intervals(
            [interval(1, 30), interval(10, 39), interval(20, 49), interval(50, 79)]
        )
        assert [lane for _, lane in packed] == [0, 1, 2, 0]

    def test_input_order_does_not_matter(self):
        forward = pack_intervals([interval(1, 30), interval(20, 49)])
        reverse = pack_intervals([interval(20, 49), interval(1, 30)])
        assert forward == reverse

    def test_every_interval_is_placed_exactly_once(self):
        intervals = [interval(i, i + 22) for i in (70, 160, 394, 463, 1001, 1074)]
        packed = pack_intervals(intervals)
        assert sorted(i.name for i, _ in packed) == sorted(i.name for i in intervals)

    def test_an_empty_selection_packs_to_nothing(self):
        assert pack_intervals([]) == []

    def test_the_nine_selected_rfx_guides_need_two_lanes(self):
        # MT007544.1_N_0160-0182 and _0161-0183 overlap, as do _1185-1207 and
        # _1201-1223; nothing else does, so two lanes is the expected packing.
        selected = [
            interval(1001, 1023), interval(1074, 1096), interval(1185, 1207),
            interval(1201, 1223), interval(161, 183), interval(160, 182),
            interval(394, 416), interval(70, 92), interval(463, 485),
        ]
        lanes = {lane for _, lane in pack_intervals(selected)}
        assert lanes == {0, 1}


class TestLaneToY:
    def test_upper_lanes_run_up_from_the_midline(self):
        assert [lane_to_y(i, "upper") for i in range(3)] == [0, 1, 2]

    def test_lower_lanes_run_down_from_below_the_midline(self):
        # -1, not 0: the midline belongs to the upper track, so the two never
        # share a row.
        assert [lane_to_y(i, "lower") for i in range(3)] == [-1, -2, -3]

    def test_an_unknown_side_raises(self):
        with pytest.raises(ValueError, match="upper"):
            lane_to_y(0, "middle")


class TestHoverTemplate:
    def test_the_label_is_bold_on_its_own_line(self):
        assert hover_template("PspCas13b", ()).startswith("<b>PspCas13b</b>")

    def test_one_line_per_field_in_order(self):
        template = hover_template(
            "RfxCas13d", (("guide_id", "g1"), ("start", "70"), ("end", "92"))
        )
        assert template.index("guide_id: g1") < template.index("start: 70")
        assert template.index("start: 70") < template.index("end: 92")

    def test_the_trace_name_box_is_suppressed(self):
        assert hover_template("PspCas13b", ()).endswith("<extra></extra>")

    def test_adding_a_field_only_adds_a_line(self):
        fewer = hover_template("X", (("a", "1"),))
        more = hover_template("X", (("a", "1"), ("b", "2")))
        assert more.count("<br>") == fewer.count("<br>") + 1


class TestGeneSpan:
    GBK = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "data",
        "SARS-CoV-2",
        "MT007544.1",
        "MT007544.1.gbk",
    )

    @pytest.mark.skipif(
        not os.path.exists(GBK), reason="MT007544.1.gbk not present"
    )
    def test_the_n_gene_span_matches_the_annotation(self):
        # 28274..29533 in the GenBank file; 1,260 nt, which is what the
        # extracted N FASTA and the located table both assume.
        assert gene_span(self.GBK, "N") == (28274, 29533)

    @pytest.mark.skipif(
        not os.path.exists(GBK), reason="MT007544.1.gbk not present"
    )
    def test_an_absent_gene_raises_and_lists_what_is_there(self):
        with pytest.raises(ValueError, match="no gene or CDS named"):
            gene_span(self.GBK, "NOT_A_GENE")


class TestPlotTracks:
    def test_a_figure_is_written(self, tmp_path):
        tracks = (
            Track("PspCas13b", (interval(160, 189), interval(1201, 1230)), "#ff7f0e", "upper"),
            Track("RfxCas13d", (interval(160, 182), interval(161, 183)), "#1f77b4", "lower"),
        )
        output = str(tmp_path / "tracks.html")
        written = plot_tracks(
            genes=[("N", 1260, tracks)], output=output, title="t",
            inline_plotlyjs=False,
        )
        assert written == output
        assert os.path.getsize(output) > 0

    def test_the_legend_carries_one_entry_per_track(self, tmp_path):
        # Two bars per track, but the legend must show two entries, not four.
        tracks = (
            Track("PspCas13b", (interval(1, 30), interval(100, 129)), "#ff7f0e", "upper"),
            Track("RfxCas13d", (interval(1, 23), interval(100, 122)), "#1f77b4", "lower"),
        )
        output = str(tmp_path / "legend.html")
        plot_tracks(genes=[("N", 1260, tracks)], output=output, title="t",
                    inline_plotlyjs=False)
        with open(output) as handle:
            content = handle.read()
        assert content.count('"showlegend":true') == 2

    def test_annotations_are_off_by_default(self, tmp_path):
        tracks = (Track("PspCas13b", (interval(160, 189, "NCP1"),), "#ff7f0e", "upper"),)
        fig_path = str(tmp_path / "plain.html")
        plot_tracks(genes=[("N", 1260, tracks)], output=fig_path, title="t",
                    inline_plotlyjs=False)
        with open(fig_path) as handle:
            content = handle.read()
        # Only the subplot title is annotated; the guide name is not printed.
        assert content.count('"text":"NCP1"') == 0

    def test_an_annotated_track_prints_each_name(self, tmp_path):
        tracks = (
            Track(
                "PspCas13b",
                (interval(160, 189, "NCP1"), interval(1201, 1230, "N_cand_93")),
                "#ff7f0e", "upper", annotate=True,
            ),
        )
        fig_path = str(tmp_path / "named.html")
        plot_tracks(genes=[("N", 1260, tracks)], output=fig_path, title="t",
                    inline_plotlyjs=False)
        with open(fig_path) as handle:
            content = handle.read()
        assert '"text":"NCP1"' in content
        assert '"text":"N_cand_93"' in content

    def test_a_guide_label_keeps_its_own_font_size(self, tmp_path):
        # The subplot titles are resized to 16; a guide label must not be
        # caught by that pass, or every name is drawn at title size.
        tracks = (
            Track("PspCas13b", (interval(160, 189, "NCP1"),), "#ff7f0e", "upper",
                  annotate=True),
        )
        fig_path = str(tmp_path / "fonts.html")
        plot_tracks(
            genes=[("N", 1260, tracks)], output=fig_path, title="t",
            annotation_font_size=10, subplot_title_font_size=16,
            inline_plotlyjs=False,
        )
        with open(fig_path) as handle:
            content = handle.read()
        # Match the annotation object that holds the label, rather than a
        # window of characters around it: the subplot title is the annotation
        # immediately before, and it legitimately carries size 16.
        objects = re.findall(r"\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}", content)
        label = next(o for o in objects if '"text":"NCP1"' in o)
        assert '"size":10' in label
        assert '"size":16' not in label

        title = next(o for o in objects if '"text":"N"' in o)
        assert '"size":16' in title

    def test_an_upper_label_sits_above_its_bar(self, tmp_path):
        from visualise_tracks import add_track
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        fig = make_subplots(rows=1, cols=1)
        track = Track("P", (interval(10, 39, "g"),), "#000", "upper", annotate=True)
        add_track(fig, track, row=1, line_width=10, show_legend=True)
        annotation = fig.layout.annotations[0]
        assert annotation.y > 0            # the bar sits at y = 0
        assert annotation.yanchor == "bottom"
        assert annotation.x == (10 + 39) / 2

    def test_a_lower_label_sits_below_its_bar(self, tmp_path):
        from visualise_tracks import add_track
        from plotly.subplots import make_subplots

        fig = make_subplots(rows=1, cols=1)
        track = Track("R", (interval(10, 32, "g"),), "#000", "lower", annotate=True)
        add_track(fig, track, row=1, line_width=10, show_legend=True)
        annotation = fig.layout.annotations[0]
        assert annotation.y < -1           # the bar sits at y = -1
        assert annotation.yanchor == "top"

    def test_a_stacked_label_follows_its_own_lane(self, tmp_path):
        from visualise_tracks import add_track
        from plotly.subplots import make_subplots

        fig = make_subplots(rows=1, cols=1)
        # Overlapping, so the second bar is on lane 1 and its label must be
        # offset from lane 1, not from the axis.
        track = Track(
            "P", (interval(10, 39, "a"), interval(20, 49, "b")), "#000", "upper",
            annotate=True,
        )
        add_track(fig, track, row=1, line_width=10, show_legend=True)
        ys = {a.text: a.y for a in fig.layout.annotations}
        assert ys["b"] > ys["a"]
        assert ys["b"] - ys["a"] == pytest.approx(1.0)

    def test_no_genes_raises(self, tmp_path):
        with pytest.raises(ValueError, match="no genes"):
            plot_tracks(genes=[], output=str(tmp_path / "x.html"), title="t")
