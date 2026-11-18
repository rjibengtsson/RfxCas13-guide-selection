"""Tests for the guide/footprint classification and position recovery."""

import os

import pandas as pd
import pytest

from visualise_scatterlines import (
    CLASS_ORDER,
    CLASS_STYLE,
    PspGuide,
    best_match,
    build_guides,
    build_psp_guides,
    LOCATED_COLUMNS,
    check_against_sequence,
    classify_overlap,
    frame_to_guides,
    guide_id,
    guides_to_frame,
    pad_width,
    psps_from_frame,
    drop_offsequence,
    filter_psps,
    locate_targets,
    plot_guides,
    profile_with_gaps,
    read_fasta,
    region_start_from_header,
    signed_distance,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N_FASTA = os.path.join(ROOT, "results", "SARS-CoV-2", "MT007544.1_N.fasta")
ON_TARGET = os.path.join(ROOT, "results", "SARS-CoV-2", "on_target.csv")
PSP_HITS = os.path.join(
    ROOT, "results", "SARS-CoV-2", "panSarbeCov_guides_MT007544.1_fmt6.tsv"
)

needs_real_data = pytest.mark.skipif(
    not (os.path.isfile(N_FASTA) and os.path.isfile(ON_TARGET)),
    reason="N gene / predictions not generated",
)

# A 30 nt footprint at 100-129, as the real ones are.
PSP = PspGuide(name="A8_NCP", reference="MT007544.1", start=100, end=129)


class TestClassifyOverlap:
    def test_a_contained_guide_is_full(self):
        # 23 nt inside the 30 nt footprint.
        assert classify_overlap(105, 127, PSP) == "full"

    def test_exactly_filling_the_footprint_is_full(self):
        assert classify_overlap(100, 129, PSP) == "full"

    def test_sharing_the_start_edge_is_full(self):
        assert classify_overlap(100, 122, PSP) == "full"

    def test_sharing_the_end_edge_is_full(self):
        assert classify_overlap(107, 129, PSP) == "full"

    def test_hanging_over_the_start_is_partial(self):
        assert classify_overlap(99, 121, PSP) == "partial"

    def test_hanging_over_the_end_is_partial(self):
        assert classify_overlap(108, 130, PSP) == "partial"

    def test_a_single_shared_base_at_the_start_is_partial(self):
        assert classify_overlap(78, 100, PSP) == "partial"

    def test_a_single_shared_base_at_the_end_is_partial(self):
        assert classify_overlap(129, 151, PSP) == "partial"

    def test_a_guide_spanning_the_whole_footprint_is_partial(self):
        # Longer than the footprint, so it is not contained by it.
        assert classify_overlap(90, 140, PSP) == "partial"

    def test_abutting_without_sharing_a_base_is_none(self):
        assert classify_overlap(77, 99, PSP) == "none"
        assert classify_overlap(130, 152, PSP) == "none"

    def test_clear_of_the_footprint_is_none(self):
        assert classify_overlap(1, 23, PSP) == "none"
        assert classify_overlap(500, 522, PSP) == "none"


class TestSignedDistance:
    def test_upstream_is_negative(self):
        assert signed_distance(70, 92, PSP) == -8

    def test_downstream_is_positive(self):
        assert signed_distance(137, 159, PSP) == 8

    def test_an_overlap_is_zero(self):
        assert signed_distance(105, 127, PSP) == 0

    def test_abutting_is_one(self):
        # One base of gap either side, so neither is 0.
        assert signed_distance(77, 99, PSP) == -1
        assert signed_distance(130, 152, PSP) == 1


class TestBestMatch:
    OTHER = PspGuide("NCP1", "MT007544.1", 200, 229)

    def test_no_footprints_gives_none(self):
        psp, cls, dist = best_match(10, 32, [])
        assert (psp, cls, dist) == (None, "none", None)

    def test_containment_beats_a_nearer_non_overlap(self):
        # Inside PSP, but closer in absolute terms to nothing else; the
        # contained relationship must win regardless of ordering.
        psp, cls, _ = best_match(105, 127, [self.OTHER, PSP])
        assert (psp.name, cls) == ("A8_NCP", "full")

    def test_a_partial_overlap_beats_a_non_overlap(self):
        psp, cls, _ = best_match(108, 130, [self.OTHER, PSP])
        assert (psp.name, cls) == ("A8_NCP", "partial")

    def test_an_overlap_is_preferred_over_a_clear_neighbour(self):
        # 180-202 overlaps OTHER (200-229) at 200-202, so that is partial.
        psp, cls, _ = best_match(180, 202, [PSP, self.OTHER])
        assert (psp.name, cls) == ("NCP1", "partial")

    def test_with_no_overlap_the_nearest_is_reported(self):
        # 150-172 sits between the two: 21 nt past PSP (100-129) and 28 nt
        # short of OTHER (200-229), so PSP is nearer.
        psp, cls, dist = best_match(150, 172, [PSP, self.OTHER])
        assert (psp.name, cls, dist) == ("A8_NCP", "none", 21)

    def test_the_nearest_of_two_distant_footprints(self):
        psp, cls, dist = best_match(135, 157, [PSP, self.OTHER])
        assert cls == "none"
        assert psp.name == "A8_NCP"
        assert dist == 6


class TestLocateTargets:
    SEQ = "ATGGGGCCCCTTTTAAAAGGGGCCCC"

    def test_a_target_is_found(self):
        df = pd.DataFrame({"Target Sequence": ["ATGGGG"]})
        assert locate_targets(df, self.SEQ) == [1]

    def test_positions_are_one_based(self):
        df = pd.DataFrame({"Target Sequence": ["TGGGG"]})
        assert locate_targets(df, self.SEQ) == [2]

    def test_lowercase_input_still_matches(self):
        df = pd.DataFrame({"Target Sequence": ["atgggg"]})
        assert locate_targets(df, self.SEQ) == [1]

    def test_an_absent_target_raises(self):
        df = pd.DataFrame({"Target Sequence": ["ATGGGG", "TTTTTTTTTTTT"]})
        with pytest.raises(ValueError, match="do not occur in the reference"):
            locate_targets(df, self.SEQ)

    def test_a_repeated_target_raises_rather_than_guessing(self):
        # GGGGCCCC occurs twice; picking the first would place it arbitrarily.
        df = pd.DataFrame({"Target Sequence": ["GGGGCCCC"]})
        with pytest.raises(ValueError, match="more than once"):
            locate_targets(df, self.SEQ)

    def test_the_error_names_the_offending_rows(self):
        df = pd.DataFrame({"Target Sequence": ["ATGGGG", "ZZZZZZ"]})
        with pytest.raises(ValueError, match="row 1"):
            locate_targets(df, self.SEQ)


class TestRegionStartFromHeader:
    def test_an_extract_feature_header(self):
        header = "MT007544.1_N MT007544.1:28274-29533(+) gene length=1260"
        assert region_start_from_header(header) == 28274

    def test_a_plain_header_gives_none(self):
        assert region_start_from_header("MT007544.1 complete genome") is None


class TestBuildPspGuides:
    def test_coordinates_are_shifted_onto_the_region(self):
        # A hit at genome 28341 in a region starting at 28274 lands at 68.
        df = pd.DataFrame(
            {"qseqid": ["A8_NCP"], "sseqid": ["MT007544.1"],
             "sstart": [28341], "send": [28370]}
        )
        psp = build_psp_guides(df, region_start=28274)[0]
        assert (psp.start, psp.end, psp.length) == (68, 97, 30)

    def test_the_default_offset_leaves_coordinates_alone(self):
        df = pd.DataFrame(
            {"qseqid": ["g"], "sseqid": ["r"], "sstart": [10], "send": [39]}
        )
        assert build_psp_guides(df)[0].start == 10

    def test_a_minus_strand_hit_is_ordered(self):
        df = pd.DataFrame(
            {"qseqid": ["g"], "sseqid": ["r"], "sstart": [39], "send": [10]}
        )
        psp = build_psp_guides(df)[0]
        assert (psp.start, psp.end) == (10, 39)

    def test_the_best_hsp_per_query_is_kept(self):
        df = pd.DataFrame(
            {
                "qseqid": ["g", "g"],
                "sseqid": ["r", "r"],
                "sstart": [10, 500],
                "send": [39, 529],
                "bitscore": [20.0, 60.0],
            }
        )
        guides = build_psp_guides(df)
        assert len(guides) == 1
        assert guides[0].start == 500

    def test_a_reference_filter_drops_other_subjects(self):
        df = pd.DataFrame(
            {"qseqid": ["a", "b"], "sseqid": ["keep", "drop"],
             "sstart": [1, 1], "send": [30, 30]}
        )
        guides = build_psp_guides(df, reference="keep")
        assert [g.name for g in guides] == ["a"]

    def test_an_empty_frame_gives_no_guides(self):
        assert build_psp_guides(pd.DataFrame()) == []


class TestDropOffsequence:
    def test_a_footprint_past_the_end_is_dropped(self, capsys):
        inside = PspGuide("in", "r", 10, 39)
        outside = PspGuide("out", "r", 2000, 2029)
        kept = drop_offsequence([inside, outside], length=1260)
        assert [g.name for g in kept] == ["in"]
        assert "not drawn" in capsys.readouterr().err

    def test_a_footprint_before_the_start_is_dropped(self):
        kept = drop_offsequence([PspGuide("neg", "r", -40, -11)], length=100)
        assert kept == []

    def test_a_footprint_straddling_the_edge_is_kept(self):
        # Partly on the region: still worth drawing, clipped by the axis.
        kept = drop_offsequence([PspGuide("edge", "r", 1250, 1279)], length=1260)
        assert len(kept) == 1


class TestFilterPsps:
    ALL = [
        PspGuide("A8_NCP", "r", 68, 97),
        PspGuide("NCP1", "r", 160, 189),
        PspGuide("Cand_3_N_Sarbecov", "r", 1074, 1103),
        PspGuide("N_cand_93", "r", 1201, 1230),
    ]

    def test_only_the_named_footprints_are_kept(self):
        kept = filter_psps(self.ALL, ["NCP1", "N_cand_93"])
        assert [g.name for g in kept] == ["NCP1", "N_cand_93"]

    def test_the_original_order_is_preserved(self):
        # Named out of order; the result still reads along the sequence.
        kept = filter_psps(self.ALL, ["N_cand_93", "A8_NCP"])
        assert [g.name for g in kept] == ["A8_NCP", "N_cand_93"]

    def test_an_unknown_name_raises(self):
        # A typo must not quietly yield a figure missing a footprint.
        with pytest.raises(ValueError, match="NCP2"):
            filter_psps(self.ALL, ["NCP1", "NCP2"])

    def test_the_error_lists_what_is_available(self):
        with pytest.raises(ValueError, match="A8_NCP"):
            filter_psps(self.ALL, ["nope"])

    def test_no_names_keeps_everything(self):
        assert filter_psps(self.ALL, []) == self.ALL

    def test_a_repeated_name_is_not_duplicated(self):
        kept = filter_psps(self.ALL, ["NCP1", "NCP1"])
        assert [g.name for g in kept] == ["NCP1"]


class TestProfileWithGaps:
    def test_an_empty_frame_is_empty(self):
        assert profile_with_gaps(pd.DataFrame()) == ([], [])

    def test_overlapping_guides_are_not_broken(self):
        df = pd.DataFrame({"start": [1, 2], "end": [23, 24], "score": [0.1, 0.2]})
        xs, ys = profile_with_gaps(df)
        assert xs == [1, 2] and None not in ys

    def test_a_gap_inserts_a_break(self):
        df = pd.DataFrame({"start": [1, 500], "end": [23, 522], "score": [0.1, 0.2]})
        xs, ys = profile_with_gaps(df)
        assert xs == [1, None, 500]
        assert ys[1] is None

    def test_abutting_guides_are_broken(self):
        # start == previous end + 1: they share no base.
        df = pd.DataFrame({"start": [1, 24], "end": [23, 46], "score": [0.1, 0.2]})
        assert profile_with_gaps(df)[0] == [1, None, 24]

    def test_a_nested_guide_does_not_create_a_gap(self):
        df = pd.DataFrame({"start": [1, 5], "end": [100, 27], "score": [0.1, 0.2]})
        assert None not in profile_with_gaps(df)[0]


class TestBuildGuides:
    SEQ = "ATG" + "ACGT" * 60   # 243 nt, no repeated 23-mers at step 1
    DF = pd.DataFrame(
        {
            "Transcript ID": ["t", "t"],
            "Target Sequence": ["ATGACGTACGTACGTACGTACGT", "TGACGTACGTACGTACGTACGTA"],
            "Guide Sequence": ["ACGTACGTACGTACGTACGTCAT", "TACGTACGTACGTACGTACGTCA"],
            "Guide Score": [0.5, 0.9],
        }
    )

    def test_every_row_yields_one_guide_and_nothing_is_filtered(self):
        guides = build_guides(self.DF, self.SEQ)
        assert len(guides) == len(self.DF)

    def test_spans_are_derived_from_the_target_length(self):
        guide = build_guides(self.DF, self.SEQ)[0]
        assert (guide.start, guide.end) == (1, 23)
        assert guide.end - guide.start + 1 == len(guide.target_seq)

    def test_with_no_footprints_everything_is_unclassified(self):
        guides = build_guides(self.DF, self.SEQ)
        assert {g.classification for g in guides} == {"none"}
        assert all(g.psp is None for g in guides)

    def test_a_footprint_classifies_the_guides_it_covers(self):
        psp = PspGuide("p", "t", 1, 30)
        guides = build_guides(self.DF, self.SEQ, [psp])
        # Guide 1 is 1-23, inside 1-30; guide 2 is 2-24, also inside.
        assert [g.classification for g in guides] == ["full", "full"]

    def test_the_score_is_carried_through_unchanged(self):
        guides = build_guides(self.DF, self.SEQ)
        assert [g.score for g in guides] == [0.5, 0.9]


class TestGuideIds:
    def test_the_span_is_zero_padded(self):
        assert guide_id("MT007544.1_N", 4, 26, 4) == "MT007544.1_N_0004-0026"

    def test_ids_sort_in_coordinate_order_as_text(self):
        # The reason for padding: unpadded, "1216" sorts before "4".
        ids = [guide_id("t", s, s + 22, 4) for s in (4, 876, 1216)]
        assert sorted(ids) == ids

    def test_the_width_covers_the_sequence_length(self):
        assert pad_width(1260) == 4
        assert pad_width(29893) == 5
        assert pad_width(100) == 4      # never narrower than 4

    def test_a_wider_sequence_widens_every_id(self):
        width = pad_width(29893)
        assert guide_id("g", 4, 26, width) == "g_00004-00026"


class TestGuidesToFrame:
    SEQ = "ATG" + "ACGT" * 60
    DF = pd.DataFrame(
        {
            "Transcript ID": ["t"],
            "Target Sequence": ["ATGACGTACGTACGTACGTACGT"],
            "Guide Sequence": ["ACGTACGTACGTACGTACGTCAT"],
            "Guide Score": [0.5],
        }
    )

    def test_an_empty_input_gives_an_empty_frame(self):
        assert guides_to_frame([]).empty

    def test_both_coordinate_systems_are_written(self):
        guides = build_guides(self.DF, self.SEQ, transcript="MT007544.1_N")
        frame = guides_to_frame(guides, region_start=28274, length=1260)
        row = frame.iloc[0]
        assert (row["start"], row["end"]) == (1, 23)
        # 1 on the region is 28274 on the reference.
        assert (row["ref_start"], row["ref_end"]) == (28274, 28296)

    def test_the_default_offset_leaves_the_reference_coordinates_equal(self):
        guides = build_guides(self.DF, self.SEQ)
        row = guides_to_frame(guides).iloc[0]
        assert row["ref_start"] == row["start"]

    def test_the_id_matches_the_region_coordinates(self):
        guides = build_guides(self.DF, self.SEQ, transcript="N")
        row = guides_to_frame(guides, region_start=28274, length=1260).iloc[0]
        assert row["guide_id"] == f"N_{row['start']:04d}-{row['end']:04d}"

    def test_rows_are_sorted_by_position(self):
        df = pd.DataFrame(
            {
                "Transcript ID": ["t", "t"],
                "Target Sequence": [
                    "TGACGTACGTACGTACGTACGTA",      # starts at 2
                    "ATGACGTACGTACGTACGTACGT",      # starts at 1
                ],
                "Guide Sequence": ["a" * 23, "b" * 23],
                "Guide Score": [0.1, 0.2],
            }
        )
        frame = guides_to_frame(build_guides(df, self.SEQ))
        assert list(frame["start"]) == [1, 2]

    def test_nothing_is_filtered(self):
        frame = guides_to_frame(build_guides(self.DF, self.SEQ))
        assert len(frame) == len(self.DF)

    def test_the_overlap_columns_are_blank_without_footprints(self):
        row = guides_to_frame(build_guides(self.DF, self.SEQ)).iloc[0]
        assert row["psp_overlap"] == "none"
        assert row["nearest_psp"] == "" and row["distance_nt"] == ""

    def test_the_overlap_columns_are_filled_with_a_footprint(self):
        psp = PspGuide("p", "MT007544.1", 1, 30)
        row = guides_to_frame(build_guides(self.DF, self.SEQ, [psp])).iloc[0]
        assert row["psp_overlap"] == "full"
        assert (row["nearest_psp"], row["psp_start"], row["psp_end"]) == ("p", 1, 30)
        assert row["distance_nt"] == 0


class TestLocatedRoundTrip:
    """guides_to_frame -> frame_to_guides must preserve every guide exactly."""

    SEQ = "ATG" + "ACGT" * 60
    DF = pd.DataFrame(
        {
            "Transcript ID": ["t", "t"],
            "Target Sequence": ["ATGACGTACGTACGTACGTACGT", "TGACGTACGTACGTACGTACGTA"],
            "Guide Sequence": ["ACGTACGTACGTACGTACGTCAT", "TACGTACGTACGTACGTACGTCA"],
            "Guide Score": [0.5, 0.9],
        }
    )
    PSP = PspGuide("p", "MT007544.1", 1, 30)

    def round_trip(self, psps=()):
        guides = build_guides(self.DF, self.SEQ, psps, transcript="t")
        frame = guides_to_frame(guides, region_start=28274, length=1260)
        return guides, frame, frame_to_guides(frame, psps)

    def test_the_frame_holds_every_declared_column(self):
        _, frame, _ = self.round_trip()
        assert set(LOCATED_COLUMNS) <= set(frame.columns)

    def test_positions_survive_the_round_trip(self):
        before, _, after = self.round_trip()
        assert [(g.start, g.end) for g in after] == [(g.start, g.end) for g in before]

    def test_sequences_and_scores_survive(self):
        before, _, after = self.round_trip()
        assert [g.spacer_seq for g in after] == [g.spacer_seq for g in before]
        assert [g.target_seq for g in after] == [g.target_seq for g in before]
        assert [g.score for g in after] == [g.score for g in before]

    def test_the_classification_survives(self):
        before, _, after = self.round_trip([self.PSP])
        assert [g.classification for g in after] == [g.classification for g in before]
        assert {g.classification for g in after} == {"full"}

    def test_the_distance_survives(self):
        before, _, after = self.round_trip([self.PSP])
        assert [g.distance for g in after] == [g.distance for g in before]

    def test_the_footprint_is_reattached(self):
        _, _, after = self.round_trip([self.PSP])
        assert all(g.psp is not None and g.psp.name == "p" for g in after)

    def test_a_guide_with_no_footprint_round_trips_as_none(self):
        _, _, after = self.round_trip()
        assert all(g.psp is None and g.distance is None for g in after)

    def test_the_stored_class_is_trusted_not_recomputed(self):
        # A table whose stored class disagrees with the geometry must come back
        # as stored: the located file is the single source of truth, and
        # silently recomputing would let a figure contradict it.
        _, frame, _ = self.round_trip([self.PSP])
        frame.loc[0, "psp_overlap"] = "partial"
        assert frame_to_guides(frame)[0].classification == "partial"

    def test_a_missing_column_is_refused(self):
        _, frame, _ = self.round_trip()
        with pytest.raises(ValueError, match="missing column"):
            frame_to_guides(frame.drop(columns=["start"]))


class TestPspsFromFrame:
    def test_footprints_are_recovered_and_deduplicated(self):
        frame = pd.DataFrame(
            {
                "nearest_psp": ["a", "a", "b"],
                "psp_start": [10, 10, 100],
                "psp_end": [39, 39, 129],
            }
        )
        psps = psps_from_frame(frame)
        assert [(p.name, p.start, p.end) for p in psps] == [("a", 10, 39), ("b", 100, 129)]

    def test_blank_rows_are_skipped(self):
        frame = pd.DataFrame(
            {"nearest_psp": ["", None], "psp_start": ["", None], "psp_end": ["", None]}
        )
        assert psps_from_frame(frame) == []

    def test_an_empty_frame_gives_nothing(self):
        assert psps_from_frame(pd.DataFrame()) == []

    def test_they_come_back_sorted_by_position(self):
        frame = pd.DataFrame(
            {"nearest_psp": ["late", "early"], "psp_start": [100, 10],
             "psp_end": [129, 39]}
        )
        assert [p.name for p in psps_from_frame(frame)] == ["early", "late"]


class TestCheckAgainstSequence:
    SEQ = "ATG" + "ACGT" * 60

    def test_matching_guides_pass(self):
        guides = build_guides(TestLocatedRoundTrip.DF, self.SEQ)
        check_against_sequence(guides, self.SEQ)   # must not raise

    def test_a_shifted_coordinate_is_caught(self):
        guides = build_guides(TestLocatedRoundTrip.DF, self.SEQ)
        shifted = [g._replace(start=g.start + 7, end=g.end + 7) for g in guides]
        with pytest.raises(ValueError, match="do not match the sequence"):
            check_against_sequence(shifted, self.SEQ)

    def test_a_different_sequence_is_caught(self):
        # The real failure this guards: a located table made against another
        # region, plotted against this one.
        guides = build_guides(TestLocatedRoundTrip.DF, self.SEQ)
        with pytest.raises(ValueError, match="disagree"):
            check_against_sequence(guides, "T" * len(self.SEQ))


class TestStyleTable:
    def test_every_class_has_a_style(self):
        assert set(CLASS_ORDER) == set(CLASS_STYLE)

    def test_the_overlap_states_are_drawn_after_the_grey(self):
        # Draw order is read from CLASS_ORDER, so the overlaps must come last
        # or they would be buried under a thousand grey points.
        assert CLASS_ORDER.index("none") == 0
        assert CLASS_ORDER[-1] == "full"

    def test_the_three_fills_are_distinct(self):
        fills = {fill for _, fill, _ in CLASS_STYLE.values()}
        assert len(fills) == 3


@needs_real_data
class TestAgainstTheRealData:
    def test_every_prediction_is_placed(self):
        _, sequence = read_fasta(N_FASTA)
        df = pd.read_csv(ON_TARGET)
        starts = locate_targets(df, sequence)
        assert len(starts) == len(df)
        assert min(starts) >= 1
        assert max(starts) + 22 <= len(sequence)

    def test_the_tiling_is_step_one(self):
        _, sequence = read_fasta(N_FASTA)
        starts = locate_targets(pd.read_csv(ON_TARGET), sequence)
        assert all(b - a == 1 for a, b in zip(starts, starts[1:]))

    @pytest.mark.skipif(not os.path.isfile(PSP_HITS), reason="psp hits not generated")
    def test_the_figure_is_written(self, tmp_path):
        header, sequence = read_fasta(N_FASTA)
        psps = build_psp_guides(
            pd.read_csv(PSP_HITS, sep="\t"),
            region_start=region_start_from_header(header),
        )
        guides = build_guides(pd.read_csv(ON_TARGET), sequence, psps)
        out = tmp_path / "plot.html"
        assert plot_guides(
            guides, psps, length=len(sequence), output=str(out),
            inline_plotlyjs=False,
        ) == str(out)
        assert out.stat().st_size > 10_000

    @pytest.mark.skipif(not os.path.isfile(PSP_HITS), reason="psp hits not generated")
    def test_a_band_is_drawn_for_every_footprint(self, tmp_path, monkeypatch):
        """Regression: add_vrect silently drops bands on an empty subplot.

        `add_vrect(row=...)` defaults to exclude_empty_subplots=True, so bands
        added before the row had any trace vanished with no error -- the figure
        rendered fine, just with no footprints marked.
        """
        import plotly.graph_objects as go

        import visualise_scatterlines as module

        header, sequence = read_fasta(N_FASTA)
        psps = build_psp_guides(
            pd.read_csv(PSP_HITS, sep="\t"),
            region_start=region_start_from_header(header),
        )
        guides = build_guides(pd.read_csv(ON_TARGET), sequence, psps)

        captured = {}
        monkeypatch.setattr(
            go.Figure, "write_html", lambda self, *a, **k: captured.setdefault("fig", self)
        )
        module.plot_guides(guides, psps, len(sequence), str(tmp_path / "p.html"))

        rects = [s for s in captured["fig"].layout.shapes if s.type == "rect"]
        assert len(rects) == len(psps) == 8
        # Each band covers whole bases, so it is half a base wider either side.
        spans = sorted((r.x0, r.x1) for r in rects)
        assert spans[0] == (psps[0].start - 0.5, psps[0].end + 0.5)

    @pytest.mark.skipif(not os.path.isfile(PSP_HITS), reason="psp hits not generated")
    def test_the_footprints_land_inside_the_gene(self):
        header, sequence = read_fasta(N_FASTA)
        psps = build_psp_guides(
            pd.read_csv(PSP_HITS, sep="\t"),
            region_start=region_start_from_header(header),
        )
        assert len(psps) == 8
        assert all(1 <= p.start and p.end <= len(sequence) for p in psps)
        assert all(p.length == 30 for p in psps)
