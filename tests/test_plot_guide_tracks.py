"""Tests for the guide-track CLI's derived columns."""

import pandas as pd
import pytest

from plot_guide_tracks import (
    RFX_EXPORT,
    build_rfx_intervals,
    display_id,
    export_rfx,
    transcribe,
)


class TestDisplayId:
    def test_the_accession_becomes_the_cas_prefix(self):
        assert display_id("MT007544.1_N_0070-0092") == "rfx_N_0070-0092"

    def test_the_gene_and_coordinates_are_kept(self):
        assert display_id("MT007544.1_N_1201-1223") == "rfx_N_1201-1223"

    def test_there_is_no_double_underscore(self):
        assert "__" not in display_id("MT007544.1_N_0070-0092")

    def test_an_id_without_the_prefix_is_unchanged(self):
        # A different accession must be shown as it is, not mangled.
        assert display_id("NC_019843.3_N_0070-0092") == "NC_019843.3_N_0070-0092"

    def test_an_already_shortened_id_is_unchanged(self):
        assert display_id("rfx_N_0070-0092") == "rfx_N_0070-0092"

    def test_the_prefix_is_only_replaced_at_the_start(self):
        # The accession appearing later in a name is not a prefix.
        assert display_id("x_MT007544.1_N_1-23") == "x_MT007544.1_N_1-23"


class TestTranscribe:
    def test_thymine_becomes_uracil(self):
        assert transcribe("ACGT") == "ACGU"

    def test_a_real_spacer_is_transcribed(self):
        # MT007544.1_N_0070-0092's spacer.
        assert transcribe("TCTCCATTCTGGTTACTGCCAGT") == "UCUCCAUUCUGGUUACUGCCAGU"

    def test_the_sequence_is_not_reverse_complemented(self):
        # The spacer is already the reverse complement of its target, so
        # transcription must leave the strand alone. Reverse complementing as
        # well would give the target's RNA and look just as plausible.
        assert transcribe("AAAAC") == "AAAAC".replace("T", "U")
        assert transcribe("ACGT")[0] == "A"

    def test_length_is_unchanged(self):
        spacer = "TGCAATTGTTTGGAGAAATCATC"
        assert len(transcribe(spacer)) == len(spacer)

    def test_no_thymine_is_left(self):
        assert "T" not in transcribe("TTTTTTT")

    def test_lowercase_is_accepted_and_upper_cased(self):
        assert transcribe("acgt") == "ACGU"

    def test_surrounding_whitespace_is_ignored(self):
        assert transcribe(" ACGT\n") == "ACGU"

    def test_a_degenerate_base_raises(self):
        # N would travel into an ordered oligo unnoticed.
        with pytest.raises(ValueError, match="unambiguous"):
            transcribe("ACGTN")

    def test_uracil_in_the_input_raises(self):
        # Already-RNA input means the column was built from the wrong source.
        with pytest.raises(ValueError, match="'U'"):
            transcribe("ACGU")

    def test_an_empty_sequence_gives_an_empty_result(self):
        assert transcribe("") == ""


class TestExportRfx:
    FRAME = pd.DataFrame(
        {
            "guide_id": ["MT007544.1_N_0070-0092"],
            "start": [70],
            "end": [92],
            "target_sequence": ["ACTGGCAGTAACCAGAATGGAGA"],
            "guide_sequence": ["TCTCCATTCTGGTTACTGCCAGT"],
            "guide_score": [0.9869713634409136],
        }
    )

    def test_the_columns_are_renamed_in_order(self):
        out = export_rfx(self.FRAME)
        assert list(out.columns) == [name for name, _ in RFX_EXPORT]
        assert list(out.columns) == [
            "guide_id", "start", "end",
            "target_dna", "spacer_dna", "spacer_rna", "guide_score",
        ]

    def test_the_renamed_columns_carry_the_original_values(self):
        out = export_rfx(self.FRAME)
        assert out["target_dna"][0] == self.FRAME["target_sequence"][0]
        assert out["spacer_dna"][0] == self.FRAME["guide_sequence"][0]

    def test_the_rna_column_is_the_transcribed_spacer(self):
        out = export_rfx(self.FRAME)
        assert out["spacer_rna"][0] == transcribe(self.FRAME["guide_sequence"][0])

    def test_the_rna_is_derived_from_the_spacer_not_the_target(self):
        # The distinction that matters: spacer_rna must pair with target_dna,
        # not repeat it.
        out = export_rfx(self.FRAME)
        assert out["spacer_rna"][0] != transcribe(self.FRAME["target_sequence"][0])

    def test_the_score_keeps_full_precision(self):
        out = export_rfx(self.FRAME)
        assert out["guide_score"][0] == self.FRAME["guide_score"][0]

    def test_the_row_count_is_unchanged(self):
        doubled = pd.concat([self.FRAME, self.FRAME], ignore_index=True)
        assert len(export_rfx(doubled)) == 2

    def test_the_exported_id_is_not_abbreviated(self):
        # The table keeps the full id so it joins back to the located table;
        # only the figure's hover is shortened.
        out = export_rfx(self.FRAME)
        assert out["guide_id"][0] == "MT007544.1_N_0070-0092"


class TestBuildRfxIntervals:
    FRAME = TestExportRfx.FRAME

    def test_the_hover_shows_the_abbreviated_id(self):
        fields = dict(build_rfx_intervals(self.FRAME)[0].fields)
        assert fields["guide_id"] == "rfx_N_0070-0092"

    def test_the_interval_keeps_the_full_id(self):
        # Identity is still the real id; only the displayed text is shortened.
        assert build_rfx_intervals(self.FRAME)[0].name == "MT007544.1_N_0070-0092"

    def test_the_coordinates_come_through_unchanged(self):
        bar = build_rfx_intervals(self.FRAME)[0]
        assert (bar.start, bar.end) == (70, 92)

    def test_the_score_is_rounded_for_display(self):
        fields = dict(build_rfx_intervals(self.FRAME)[0].fields)
        assert fields["guide_score"] == "0.9870"
