"""Tests for src/guide_neighbours.py and scripts/extract_neighbours.py."""

import os

import pandas as pd
import pytest

from extract_neighbours import main
from guide_neighbours import (
    DEFAULT_WINDOW,
    DOWNSTREAM,
    NEIGHBOUR_CAS,
    FULL_OVERLAP,
    NO_NEIGHBOUR,
    OUTPUT_COLUMNS,
    PARTIAL_OVERLAP,
    QUERY_CAS,
    UPSTREAM,
    NeighbourDataError,
    _interval_gap,
    _positions,
    find_neighbours,
    gene_label,
    load_annotated,
    summarise,
    write_neighbours,
)

# Query guides are 30 nt, neighbours 23 nt, as in the real data. Q1 sits at
# 1000-1029 and the neighbour gaps from it are, in order: N6 177, N1 77, N2 7,
# N3 overlapping, N4 10, N5 170.
ROWS = [
    # Cas, Index, Guide Seq, Target Seq, Score, ref_start, ref_end, strand, features
    (QUERY_CAS, 1, "Q1g", "Q1t", "Unknown", 1000, 1029, "+", "gene:A;CDS:A"),
    (QUERY_CAS, 2, "Q2g", "Q2t", "Unknown", 1100, 1129, "+", "gene:B"),
    (QUERY_CAS, 3, "Q3g", "Q3t", "Unknown", 5000, 5029, "+", "gene:C"),
    (NEIGHBOUR_CAS, 11, "N1g", "N1t", "0.11", 900, 922, "+", "gene:A"),
    (NEIGHBOUR_CAS, 12, "N2g", "N2t", "0.22", 970, 992, "+", "gene:A;CDS:A"),
    (NEIGHBOUR_CAS, 13, "N3g", "N3t", "0.33", 1010, 1032, "-", "gene:A;gene:B"),
    (NEIGHBOUR_CAS, 14, "N4g", "N4t", "0.44", 1040, 1062, "+", "gene:B"),
    (NEIGHBOUR_CAS, 15, "N5g", "N5t", "0.55", 1200, 1222, "+", "3'UTR:3'UTR"),
    (NEIGHBOUR_CAS, 16, "N6g", "N6t", "0.66", 800, 822, "+", "gene:A"),
]
COLUMNS = [
    "Cas", "Index", "Guide Sequence", "Target Sequence", "Guide Score",
    "ref_start", "ref_end", "strand", "features",
]


@pytest.fixture
def annotated_csv(tmp_path):
    path = tmp_path / "annotated.csv"
    pd.DataFrame(ROWS, columns=COLUMNS).to_csv(path, index=False)
    return str(path)


@pytest.fixture
def df(annotated_csv):
    return load_annotated(annotated_csv)


def gaps(pairs):
    """The bases between each pair's two guides, from the reported coordinates.

    The table no longer carries the gap, so a test that cares about the window
    measures it the same way find_neighbours does.
    """
    return _interval_gap(
        pairs.neighbour_ref_start, pairs.neighbour_ref_end,
        pairs.query_ref_start, pairs.query_ref_end,
    )


def indexes(pairs, query=1):
    """Neighbour Index values reported for one query guide, in row order."""
    sub = pairs[(pairs.query_guide_index == query) & (pairs.relative_position != NO_NEIGHBOUR)]
    return list(sub.neighbour_guide_index)


class TestGeneLabel:
    def test_takes_the_gene_name(self):
        assert gene_label("gene:S;CDS:S") == "S"

    def test_deduplicates_repeated_names(self):
        assert gene_label("CDS:orf1ab;CDS:orf1ab;gene:orf1ab") == "orf1ab"

    def test_reports_both_genes_of_an_overlap(self):
        assert gene_label("CDS:orf1ab;gene:orf1ab;CDS:S;gene:S") == "orf1ab;S"

    def test_falls_back_to_the_feature_type_for_a_utr(self):
        assert gene_label("5'UTR:5'UTR") == "5'UTR"
        assert gene_label("CDS:N;gene:N;3'UTR:3'UTR") == "N;3'UTR"

    def test_missing_annotation_is_empty(self):
        assert gene_label(None) == ""
        assert gene_label(float("nan")) == ""


class TestLoadAnnotated:
    def test_columns_are_aliased_to_safe_names(self, df):
        assert {"cas", "guide_index", "guide_sequence", "target_sequence", "guide_score"} <= set(
            df.columns
        )

    def test_a_gene_column_is_added(self, df):
        assert df.set_index("guide_index").loc[13, "gene"] == "A;B"

    def test_coordinates_become_integers(self, df):
        assert df.ref_start.dtype.kind == "i" and df.ref_end.dtype.kind == "i"

    def test_unmapped_guides_are_dropped(self, tmp_path):
        """A no_hit row has no coordinates, so it cannot take part."""
        rows = ROWS + [(NEIGHBOUR_CAS, 99, "Xg", "Xt", "0.9", None, None, None, None)]
        path = tmp_path / "with_unmapped.csv"
        pd.DataFrame(rows, columns=COLUMNS).to_csv(path, index=False)
        assert 99 not in set(load_annotated(str(path)).guide_index)

    def test_a_csv_without_annotation_is_refused(self, tmp_path):
        path = tmp_path / "plain.csv"
        pd.DataFrame({"Cas": [QUERY_CAS], "Index": [1]}).to_csv(path, index=False)
        with pytest.raises(NeighbourDataError, match="features|blast_guides"):
            load_annotated(str(path))


class TestIntervalGap:
    def test_overlapping_and_touching_are_zero(self):
        assert _interval_gap(10, 20, 15, 25) == 0
        assert _interval_gap(10, 20, 21, 30) == 0

    def test_the_gap_counts_the_bases_between(self):
        assert _interval_gap(10, 20, 26, 30) == 5
        assert _interval_gap(26, 30, 10, 20) == 5

    def test_it_is_symmetric(self):
        assert _interval_gap(10, 20, 100, 110) == _interval_gap(100, 110, 10, 20)


class TestPositions:
    def test_labels_each_side(self):
        assert _positions(1000, 1029, 1010, 1032) == PARTIAL_OVERLAP
        assert _positions(1000, 1029, 900, 922) == UPSTREAM
        assert _positions(1000, 1029, 1040, 1062) == DOWNSTREAM

    def test_a_touching_neighbour_is_not_overlapping(self):
        assert _positions(1000, 1029, 1030, 1052) == DOWNSTREAM

    def test_a_neighbour_inside_the_query_overlaps_fully(self):
        """A 23 nt neighbour wholly within the query's 30 nt site."""
        assert _positions(1000, 1029, 1003, 1025) == FULL_OVERLAP
        assert _positions(1000, 1029, 1000, 1029) == FULL_OVERLAP

    def test_a_neighbour_containing_the_query_overlaps_fully(self):
        """Containment counts whichever guide is the longer one."""
        assert _positions(1010, 1020, 1000, 1030) == FULL_OVERLAP

    def test_an_overlap_running_off_an_end_is_partial(self):
        assert _positions(1000, 1029, 980, 1002) == PARTIAL_OVERLAP
        assert _positions(1000, 1029, 1028, 1050) == PARTIAL_OVERLAP

    def test_a_single_shared_base_is_a_partial_overlap(self):
        assert _positions(1000, 1029, 978, 1000) == PARTIAL_OVERLAP
        assert _positions(1000, 1029, 1029, 1051) == PARTIAL_OVERLAP


class TestFindNeighbours:
    def test_a_symmetric_window_takes_both_sides(self, df):
        pairs = find_neighbours(df, upstream=100, downstream=100)
        assert indexes(pairs) == [11, 12, 13, 14]  # 77, 7, overlap, 10

    def test_rows_are_ordered_by_genome_position(self, df):
        pairs = find_neighbours(df, upstream=200, downstream=200)
        sub = pairs[pairs.query_guide_index == 1]
        assert list(sub.neighbour_ref_start) == sorted(sub.neighbour_ref_start)

    def test_upstream_only(self, df):
        pairs = find_neighbours(df, upstream=100, downstream=0)
        assert indexes(pairs) == [11, 12, 13]

    def test_downstream_only(self, df):
        pairs = find_neighbours(df, upstream=0, downstream=100)
        assert indexes(pairs) == [13, 14]

    def test_the_two_sides_are_independent(self, df):
        """77 bp upstream and 10 bp downstream: each allowance applies alone."""
        pairs = find_neighbours(df, upstream=77, downstream=9)
        assert indexes(pairs) == [11, 12, 13]
        pairs = find_neighbours(df, upstream=76, downstream=10)
        assert indexes(pairs) == [12, 13, 14]

    def test_an_overlapping_neighbour_survives_a_zero_window(self, df):
        pairs = find_neighbours(df, upstream=0, downstream=0)
        assert indexes(pairs) == [13]
        assert pairs[pairs.query_guide_index == 1].relative_position.tolist() == [PARTIAL_OVERLAP]

    def test_the_gap_is_not_a_column(self, df):
        """The window is measured in gap bases, but a gap of 0 covers both a
        23 bp overlap and a mere abutment, so it is not reported. See
        TestIntervalGap for the measure itself.
        """
        pairs = find_neighbours(df, upstream=100, downstream=100)
        assert "distance_bp" not in pairs.columns

    def test_whole_numbered_columns_stay_integers(self, df):
        """A query with no neighbours must not turn coordinates into floats."""
        pairs = find_neighbours(df, upstream=100, downstream=100)
        for column in ("n_neighbours", "neighbour_guide_index",
                       "neighbour_ref_start", "neighbour_ref_end"):
            assert str(pairs[column].dtype) == "Int64", column

    def test_integers_survive_the_round_trip_to_csv(self, df, tmp_path):
        path = write_neighbours(find_neighbours(df, upstream=100, downstream=100),
                                str(tmp_path / "out.csv"))
        text = open(path).read()
        assert ".0," not in text

    def test_a_wide_window_takes_everything(self, df):
        pairs = find_neighbours(df, upstream=200, downstream=200)
        assert indexes(pairs) == [16, 11, 12, 13, 14, 15]

    def test_a_neighbour_near_two_queries_is_reported_under_both(self, df):
        pairs = find_neighbours(df, upstream=100, downstream=100)
        assert 14 in indexes(pairs, query=1) and 14 in indexes(pairs, query=2)

    def test_each_query_gets_its_own_position(self, df):
        """N4 sits above Q1 and below Q2, so the one neighbour is labelled from
        each query guide's own point of view."""
        pairs = find_neighbours(df, upstream=100, downstream=100)
        shared = pairs[pairs.neighbour_guide_index == 14].set_index("query_guide_index")
        assert shared.loc[1, "relative_position"] == DOWNSTREAM
        assert shared.loc[2, "relative_position"] == UPSTREAM

    def test_a_query_with_no_neighbour_is_marked_not_dropped(self, df):
        pairs = find_neighbours(df, upstream=100, downstream=100)
        q3 = pairs[pairs.query_guide_index == 3]
        assert len(q3) == 1
        assert q3.relative_position.iloc[0] == NO_NEIGHBOUR
        assert q3.n_neighbours.iloc[0] == 0
        assert pd.isna(q3.neighbour_guide_index.iloc[0])

    def test_empty_queries_can_be_dropped(self, df):
        pairs = find_neighbours(df, upstream=100, downstream=100, include_empty=False)
        assert 3 not in set(pairs.query_guide_index)

    def test_n_neighbours_matches_the_rows_present(self, df):
        pairs = find_neighbours(df, upstream=100, downstream=100)
        for q, sub in pairs.groupby("query_guide_index"):
            real = sub[sub.relative_position != NO_NEIGHBOUR]
            assert sub.n_neighbours.iloc[0] == len(real), q

    def test_the_query_set_can_be_restricted(self, df):
        pairs = find_neighbours(df, query_indexes=[1])
        assert set(pairs.query_guide_index) == {1}

    def test_both_guides_details_are_reported(self, df):
        pairs = find_neighbours(df, upstream=100, downstream=100)
        row = pairs[pairs.neighbour_guide_index == 12].iloc[0]
        assert row.query_guide_sequence == "Q1g" and row.query_target_sequence == "Q1t"
        assert row.query_gene == "A" and row.query_ref_start == 1000
        assert row.neighbour_guide_sequence == "N2g" and row.neighbour_target_sequence == "N2t"
        assert float(row.neighbour_guide_score) == 0.22

    def test_the_columns_are_exactly_the_documented_set(self, df):
        assert list(find_neighbours(df).columns) == OUTPUT_COLUMNS

    def test_no_blast_statistics_are_carried_over(self, df):
        """The annotated CSV's stats stay behind; this table is not for them."""
        cols = " ".join(find_neighbours(df).columns)
        for stat in ("pident", "evalue", "bitscore", "qcovhsp", "mismatch", "aln_length"):
            assert stat not in cols

    def test_an_unknown_query_cas_is_an_error_listing_what_exists(self, df):
        with pytest.raises(NeighbourDataError, match=QUERY_CAS):
            find_neighbours(df, query_cas="NoSuchCas")

    def test_an_unmatched_query_index_is_an_error(self, df):
        with pytest.raises(NeighbourDataError, match="Index"):
            find_neighbours(df, query_indexes=[999])

    def test_a_negative_window_is_an_error(self, df):
        with pytest.raises(NeighbourDataError, match="0 or more"):
            find_neighbours(df, upstream=-1)

    def test_no_neighbour_cas_leaves_every_query_empty(self, df):
        """No frame carries the neighbour columns, so the table must still be
        built with the full column set rather than raising."""
        pairs = find_neighbours(df, neighbour_cas="NoSuchCas")
        assert (pairs.relative_position == NO_NEIGHBOUR).all()
        assert len(pairs) == 3
        assert list(pairs.columns) == OUTPUT_COLUMNS
        assert pairs.neighbour_guide_index.isna().all()

    def test_the_default_window_is_100(self):
        assert DEFAULT_WINDOW == 100


class TestSummarise:
    def test_one_row_per_query_with_the_position_breakdown(self, df):
        summary = summarise(find_neighbours(df, upstream=100, downstream=100))
        assert list(summary.index) == [1, 2, 3]
        assert summary.loc[1, PARTIAL_OVERLAP] == 1  # N3 runs off the query's end
        assert summary.loc[1, FULL_OVERLAP] == 0
        assert summary.loc[1, UPSTREAM] == 2
        assert summary.loc[1, DOWNSTREAM] == 1
        assert summary.loc[1, "n_neighbours"] == 4

    def test_an_empty_query_reads_as_zero(self, df):
        summary = summarise(find_neighbours(df, upstream=100, downstream=100))
        assert summary.loc[3, "n_neighbours"] == 0
        assert summary.loc[
            3, [FULL_OVERLAP, PARTIAL_OVERLAP, UPSTREAM, DOWNSTREAM]
        ].tolist() == [0, 0, 0, 0]

    def test_the_query_gene_and_position_are_carried(self, df):
        summary = summarise(find_neighbours(df))
        assert summary.loc[1, "query_gene"] == "A"
        assert summary.loc[1, "query_ref_start"] == 1000

    def test_an_empty_table_summarises_to_nothing(self):
        assert summarise(pd.DataFrame()).empty


class TestWriteNeighbours:
    def test_it_creates_the_directory(self, df, tmp_path):
        path = write_neighbours(find_neighbours(df), str(tmp_path / "new" / "out.csv"))
        assert os.path.isfile(path)

    def test_the_file_round_trips(self, df, tmp_path):
        pairs = find_neighbours(df, upstream=100, downstream=100)
        path = write_neighbours(pairs, str(tmp_path / "out.csv"))
        assert list(pd.read_csv(path).columns) == OUTPUT_COLUMNS
        assert len(pd.read_csv(path)) == len(pairs)


class TestCli:
    def test_a_symmetric_window(self, annotated_csv, tmp_path, capsys):
        out = tmp_path / "out.csv"
        assert main(["--input", annotated_csv, "--output", str(out), "--window", "100"]) == 0
        printed = capsys.readouterr().out
        assert "100 bp upstream and 100 bp downstream" in printed
        assert indexes(pd.read_csv(out)) == [11, 12, 13, 14]

    def test_an_asymmetric_window(self, annotated_csv, tmp_path, capsys):
        out = tmp_path / "out.csv"
        assert main(["--input", annotated_csv, "--output", str(out),
                     "--upstream", "0", "--downstream", "100"]) == 0
        assert "0 bp upstream and 100 bp downstream" in capsys.readouterr().out
        assert indexes(pd.read_csv(out)) == [13, 14]

    def test_one_side_overrides_the_symmetric_window(self, annotated_csv, tmp_path):
        out = tmp_path / "out.csv"
        main(["--input", annotated_csv, "--output", str(out),
              "--window", "100", "--upstream", "0"])
        assert indexes(pd.read_csv(out)) == [13, 14]

    def test_a_query_can_be_named_by_sequence(self, annotated_csv, tmp_path):
        out = tmp_path / "out.csv"
        assert main(["--input", annotated_csv, "--output", str(out),
                     "--query-sequence", "Q1t"]) == 0
        assert set(pd.read_csv(out).query_guide_index) == {1}

    def test_the_sequence_lookup_accepts_the_guide_strand_and_rna(self, annotated_csv, tmp_path):
        out = tmp_path / "out.csv"
        main(["--input", annotated_csv, "--output", str(out), "--query-sequence", "q1g"])
        assert set(pd.read_csv(out).query_guide_index) == {1}

    def test_an_unknown_sequence_exits_one(self, annotated_csv, tmp_path, capsys):
        with pytest.raises(NeighbourDataError, match="no PspCas13b guide"):
            main(["--input", annotated_csv, "--output", str(tmp_path / "out.csv"),
                  "--query-sequence", "ACGTACGT"])

    def test_a_query_can_be_named_by_index(self, annotated_csv, tmp_path):
        out = tmp_path / "out.csv"
        main(["--input", annotated_csv, "--output", str(out), "--query-index", "1,2"])
        assert set(pd.read_csv(out).query_guide_index) == {1, 2}

    def test_drop_empty_removes_the_barren_query(self, annotated_csv, tmp_path):
        out = tmp_path / "out.csv"
        main(["--input", annotated_csv, "--output", str(out), "--drop-empty"])
        assert 3 not in set(pd.read_csv(out).query_guide_index)

    def test_a_summary_file_can_be_written(self, annotated_csv, tmp_path):
        out, summary = tmp_path / "out.csv", tmp_path / "summary.csv"
        main(["--input", annotated_csv, "--output", str(out), "--summary", str(summary)])
        assert os.path.isfile(summary)
        assert "n_neighbours" in pd.read_csv(summary).columns

    def test_a_negative_window_exits_two(self, annotated_csv, tmp_path, capsys):
        code = main(["--input", annotated_csv, "--output", str(tmp_path / "o.csv"),
                     "--window", "-5"])
        assert code == 2
        assert "0 or more" in capsys.readouterr().err

    def test_the_empty_query_count_is_reported(self, annotated_csv, tmp_path, capsys):
        main(["--input", annotated_csv, "--output", str(tmp_path / "o.csv")])
        assert "had no RfxCas13d neighbour" in capsys.readouterr().out


REAL_CSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "results",
    "MERS",
    "mers_guides_NC_019843_annotated.csv",
)


@pytest.mark.skipif(not os.path.isfile(REAL_CSV), reason="annotated results not generated")
class TestAgainstRealResults:
    @pytest.fixture(scope="class")
    @classmethod
    def real(cls):
        return load_annotated(REAL_CSV)

    def test_all_seventeen_psp_guides_are_reported(self, real):
        pairs = find_neighbours(real, upstream=100, downstream=100)
        assert pairs.query_guide_index.nunique() == 17

    def test_the_named_guide_has_131_neighbours(self, real):
        """TTAAGTGAATAGCTTGGCTATCTCACTTCC is Index 0, at genome 4-33."""
        pairs = find_neighbours(real, upstream=100, downstream=100, query_indexes=[0])
        assert pairs.query_target_sequence.iloc[0] == "TTAAGTGAATAGCTTGGCTATCTCACTTCC"
        assert (pairs.query_ref_start.iloc[0], pairs.query_ref_end.iloc[0]) == (4, 33)
        assert len(pairs) == 131
        # 8 of the 23 nt neighbours sit wholly inside the 30 nt query site
        assert (pairs.relative_position == FULL_OVERLAP).sum() == 8
        assert (pairs.relative_position == PARTIAL_OVERLAP).sum() == 22
        assert (pairs.relative_position == DOWNSTREAM).sum() == 101

    def test_distances_stay_inside_the_window(self, real):
        pairs = find_neighbours(real, upstream=100, downstream=100, include_empty=False)
        assert gaps(pairs).max() <= 100

    def test_an_asymmetric_window_is_honoured(self, real):
        pairs = find_neighbours(real, upstream=0, downstream=100, include_empty=False)
        gap = gaps(pairs)
        # Only a directly abutting guide qualifies upstream: a 0 bp allowance
        # still admits a 0 bp gap.
        assert (gap[pairs.relative_position == UPSTREAM] == 0).all()
        downstream_gaps = gap[pairs.relative_position == DOWNSTREAM]
        assert (downstream_gaps > 0).any()
        assert downstream_gaps.max() <= 100

    def test_six_psp_guides_have_no_neighbour_within_100bp(self, real):
        """RfxCas13d coverage has real gaps: no guides at all between 13432 and
        21459, 25495 and 27593, or 27816 and 29104."""
        summary = summarise(find_neighbours(real, upstream=100, downstream=100))
        assert sorted(summary[summary.n_neighbours == 0].index) == [3, 10, 11, 12, 13, 14]

    def test_widening_the_window_reaches_the_m_gene_guides(self, real):
        """Psp Index 10's nearest RfxCas13d guide is 105 bp away."""
        summary = summarise(find_neighbours(real, upstream=150, downstream=150))
        assert summary.loc[10, "n_neighbours"] > 0
        assert summary.loc[3, "n_neighbours"] == 0  # still 3,779 bp from any guide

    def test_the_saturated_count_matches_the_guide_tiling(self, real):
        """Where RfxCas13d guides tile every base, a 30 nt query overlaps
        30 + 23 - 1 = 52 of them, with 101 more at gaps 0-100 on each side.

        Of those 52, the 30 - 23 + 1 = 8 starting inside the query's first 8
        bases fit wholly within it; the other 44 run off one end or the other.
        """
        summary = summarise(find_neighbours(real, upstream=100, downstream=100))
        saturated = summary[summary.n_neighbours == 254]
        assert not saturated.empty
        assert (saturated[FULL_OVERLAP] == 8).all()
        assert (saturated[PARTIAL_OVERLAP] == 44).all()
        assert (saturated[UPSTREAM] == 101).all()
        assert (saturated[DOWNSTREAM] == 101).all()

    def test_the_query_guide_gene_is_reported(self, real):
        pairs = find_neighbours(real, upstream=100, downstream=100, include_empty=False)
        assert (pairs.query_gene.str.len() > 0).all()
        assert "N;3'UTR" in set(pairs.query_gene)  # Index 15 spans both
