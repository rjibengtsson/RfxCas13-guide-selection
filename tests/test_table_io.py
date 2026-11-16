"""Tests for the delimiter handling in src/table_io.py."""

import os

import pytest

from table_io import read_table, resolve_separator, sniff_separator

GUIDES_TSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "pspcas13b_guides",
    "panSarbeCov_guides.tsv",
)


class TestResolveSeparator:
    @pytest.mark.parametrize(
        "spec,expected",
        [
            ("auto", None),
            ("tab", "\t"),
            ("tsv", "\t"),
            ("\\t", "\t"),
            ("csv", ","),
            ("comma", ","),
            (",", ","),
            (";", ";"),
            ("pipe", "|"),
            ("TAB", "\t"),
            ("  tab  ", "\t"),
        ],
    )
    def test_known_spellings(self, spec, expected):
        assert resolve_separator(spec) == expected

    def test_none_means_auto(self):
        assert resolve_separator(None) is None

    def test_a_single_character_is_taken_literally(self):
        assert resolve_separator("~") == "~"

    def test_an_unknown_word_is_rejected(self):
        with pytest.raises(ValueError, match="unrecognised"):
            resolve_separator("tabulated")


class TestSniffSeparator:
    @pytest.mark.parametrize(
        "name,expected",
        [("g.tsv", "\t"), ("g.csv", ","), ("g.tab", "\t"), ("g.txt", "\t")],
    )
    def test_extension_decides(self, tmp_path, name, expected):
        p = tmp_path / name
        p.write_text("Name,Spacer\na,ACGT\n")  # content disagrees on purpose
        assert sniff_separator(str(p)) == expected

    def test_unknown_extension_falls_back_to_the_header(self, tmp_path):
        p = tmp_path / "guides.dat"
        p.write_text("Name\tSpacer\na\tACGT\n")
        assert sniff_separator(str(p)) == "\t"

    def test_a_single_column_file_defaults_to_comma(self, tmp_path):
        p = tmp_path / "guides.dat"
        p.write_text("Spacer\nACGT\n")
        assert sniff_separator(str(p)) == ","

    def test_leading_blank_lines_are_skipped(self, tmp_path):
        p = tmp_path / "guides.dat"
        p.write_text("\n\nName;Spacer\na;ACGT\n")
        assert sniff_separator(str(p)) == ";"

    def test_an_empty_file_is_an_error(self, tmp_path):
        p = tmp_path / "guides.dat"
        p.write_text("")
        with pytest.raises(ValueError, match="empty"):
            sniff_separator(str(p))


class TestReadTable:
    def test_a_tsv_is_split_into_columns(self, tmp_path):
        p = tmp_path / "g.tsv"
        p.write_text("Name\tSpacer\nguide1\tACGT\nguide2\tTTTT\n")
        df = read_table(str(p))
        assert list(df.columns) == ["Name", "Spacer"]
        assert df["Spacer"].tolist() == ["ACGT", "TTTT"]

    def test_an_explicit_separator_overrides_the_extension(self, tmp_path):
        # Tab-separated content in a .csv file: 'auto' would get this wrong.
        p = tmp_path / "g.csv"
        p.write_text("Name\tSpacer\nguide1\tACGT\n")
        assert read_table(str(p), "tab")["Spacer"].tolist() == ["ACGT"]

    def test_a_missing_file_is_an_error(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_table(str(tmp_path / "absent.tsv"))

    def test_a_single_column_parse_warns(self, tmp_path, capsys):
        p = tmp_path / "g.csv"
        p.write_text("Name\tSpacer\nguide1\tACGT\n")  # wrong delimiter for .csv
        read_table(str(p))
        assert "single column" in capsys.readouterr().out

    @pytest.mark.skipif(not os.path.isfile(GUIDES_TSV), reason="guide TSV not present")
    def test_the_real_guide_tsv_reads_as_two_columns(self):
        df = read_table(GUIDES_TSV)
        assert list(df.columns) == ["Name", "Spacer"]
        assert len(df) == 8
        # Every spacer is a 30 nt PspCas13b spacer.
        assert set(df["Spacer"].str.len()) == {30}
