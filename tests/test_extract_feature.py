"""Tests for scripts/extract_feature.py."""

import os

import pytest

from conftest import GENOME_FASTA, GENOME_GBK, needs_genome
from blast_utils import reverse_complement
from extract_feature import main, parse_coords, slice_parts

SARS_GBK = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "SARS-CoV-2",
    "MT007544.1",
    "MT007544.1.gbk",
)
SARS_FASTA = SARS_GBK.replace(".gbk", ".fasta")

needs_sars = pytest.mark.skipif(
    not os.path.isfile(SARS_GBK), reason="MT007544.1 not downloaded"
)


def read_fasta(path):
    """Return (header, sequence) of the first record, sequence unwrapped."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:].rstrip("\n")
            else:
                seq.append(line.strip())
    return header, "".join(seq)


def read_all(path):
    """Return [(header, sequence), ...] for every record in the file."""
    entries, seq = [], []
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq)))
                header, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq)))
    return entries


class TestParseCoords:
    @pytest.mark.parametrize(
        "spec", ["28274-29533", "28274..29533", "28,274-29,533", " 28274 - 29533 "]
    )
    def test_accepted_spellings(self, spec):
        assert parse_coords(spec) == (28274, 29533)

    def test_a_single_base(self):
        assert parse_coords("5-5") == (5, 5)

    @pytest.mark.parametrize("spec", ["", "abc", "100", "100-", "-100", "1e3-2e3"])
    def test_malformed(self, spec):
        with pytest.raises(ValueError, match="cannot read"):
            parse_coords(spec)

    def test_a_reversed_interval_is_rejected(self):
        # Silently swapping would hide a typo and give a plausible wrong answer.
        with pytest.raises(ValueError, match="before its start"):
            parse_coords("29533-28274")

    def test_a_zero_start_is_rejected(self):
        with pytest.raises(ValueError, match="1 or greater"):
            parse_coords("0-100")


class TestSliceParts:
    SEQ = "AAAACCCCGGGGTTTT"  # 1-based: A 1-4, C 5-8, G 9-12, T 13-16

    def test_a_single_segment_on_the_plus_strand(self):
        assert slice_parts(self.SEQ, [(5, 8)], "+") == "CCCC"

    def test_a_minus_strand_segment_is_reverse_complemented(self):
        assert slice_parts(self.SEQ, [(1, 4)], "-") == "TTTT"

    def test_segments_are_joined_in_positional_order(self):
        assert slice_parts(self.SEQ, [(1, 2), (9, 10)], "+") == "AAGG"

    def test_a_spliced_minus_strand_feature_joins_then_revcomps(self):
        # Order of operations matters: join ascending, then revcomp the whole.
        assert slice_parts(self.SEQ, [(1, 2), (9, 10)], "-") == reverse_complement("AAGG")

    def test_the_whole_sequence(self):
        assert slice_parts(self.SEQ, [(1, len(self.SEQ))], "+") == self.SEQ


@needs_sars
class TestExtractTheNGene:
    @pytest.fixture(scope="class")
    @classmethod
    def genome(cls):
        return read_fasta(SARS_FASTA)[1]

    def test_the_n_gene_matches_the_genome_slice(self, tmp_path, genome):
        out = tmp_path / "N.fasta"
        assert main(["--reference", SARS_GBK, "--gene", "N", "--output", str(out)]) == 0
        header, seq = read_fasta(str(out))
        # N is 28274..29533 on the plus strand in this record.
        assert seq == genome[28274 - 1 : 29533]
        assert len(seq) == 1260
        assert header.startswith("MT007544.1_N ")
        assert "28274-29533(+)" in header

    def test_the_n_gene_is_a_complete_reading_frame(self, tmp_path):
        # A sanity check on the coordinates themselves, not on the code.
        out = tmp_path / "N.fasta"
        main(["--reference", SARS_GBK, "--feature-type", "CDS", "--gene", "N",
              "--output", str(out)])
        seq = read_fasta(str(out))[1]
        assert len(seq) % 3 == 0
        assert seq.startswith("ATG")
        assert seq[-3:] in {"TAA", "TAG", "TGA"}

    def test_coords_reproduce_the_named_feature(self, tmp_path):
        by_name = tmp_path / "a.fasta"
        by_coords = tmp_path / "b.fasta"
        main(["--reference", SARS_GBK, "--gene", "N", "--output", str(by_name)])
        main(["--reference", SARS_FASTA, "--coords", "28274-29533", "--output", str(by_coords)])
        assert read_fasta(str(by_name))[1] == read_fasta(str(by_coords))[1]

    def test_the_sequence_is_wrapped(self, tmp_path):
        out = tmp_path / "N.fasta"
        main(["--reference", SARS_GBK, "--gene", "N", "--output", str(out)])
        lines = [l for l in out.read_text().splitlines() if not l.startswith(">")]
        assert all(len(l) <= 60 for l in lines)
        assert len(lines) > 1

    def test_a_custom_name_is_used(self, tmp_path):
        out = tmp_path / "N.fasta"
        main(["--reference", SARS_GBK, "--gene", "N", "--name", "my_N", "--output", str(out)])
        assert read_fasta(str(out))[0].startswith("my_N ")

    def test_gene_matching_is_case_insensitive(self, tmp_path):
        out = tmp_path / "N.fasta"
        assert main(["--reference", SARS_GBK, "--gene", "n", "--output", str(out)]) == 0

    def test_a_product_substring_selects_the_cds(self, tmp_path):
        out = tmp_path / "N.fasta"
        code = main(["--reference", SARS_GBK, "--feature-type", "CDS",
                     "--product", "nucleocapsid", "--output", str(out)])
        assert code == 0
        assert len(read_fasta(str(out))[1]) == 1260

    def test_the_output_directory_is_created(self, tmp_path):
        out = tmp_path / "new" / "dir" / "N.fasta"
        main(["--reference", SARS_GBK, "--gene", "N", "--output", str(out)])
        assert out.is_file()


@needs_sars
class TestErrors:
    def test_an_unknown_gene_lists_what_is_available(self, tmp_path, capsys):
        out = tmp_path / "x.fasta"
        code = main(["--reference", SARS_GBK, "--gene", "NOPE", "--output", str(out)])
        assert code == 2
        err = capsys.readouterr().err
        assert "no gene feature matching" in err
        assert "available labels" in err and "N" in err
        assert not out.exists()

    def test_no_selector_is_an_error(self, tmp_path, capsys):
        code = main(["--reference", SARS_GBK, "--output", str(tmp_path / "x.fasta")])
        assert code == 2
        assert "exactly one of" in capsys.readouterr().err

    def test_two_selectors_are_an_error(self, tmp_path, capsys):
        code = main(["--reference", SARS_GBK, "--gene", "N", "--coords", "1-10",
                     "--output", str(tmp_path / "x.fasta")])
        assert code == 2
        assert "exactly one of" in capsys.readouterr().err

    def test_a_named_feature_against_a_fasta_reference_is_refused(self, tmp_path, capsys):
        # FASTA has no annotation, so this must say so rather than fail obscurely.
        code = main(["--reference", SARS_FASTA, "--gene", "N",
                     "--output", str(tmp_path / "x.fasta")])
        assert code == 2
        assert "carries no annotation" in capsys.readouterr().err

    def test_coords_past_the_end_are_refused(self, tmp_path, capsys):
        code = main(["--reference", SARS_FASTA, "--coords", "29800-30500",
                     "--output", str(tmp_path / "x.fasta")])
        assert code == 2
        assert "past the end" in capsys.readouterr().err

    def test_an_ambiguous_match_is_refused_without_all(self, tmp_path, capsys):
        # 'any' matches both the gene and the CDS for N.
        code = main(["--reference", SARS_GBK, "--feature-type", "any", "--gene", "N",
                     "--output", str(tmp_path / "x.fasta")])
        assert code == 2
        assert "features match" in capsys.readouterr().err

    def test_all_writes_every_match(self, tmp_path):
        out = tmp_path / "x.fasta"
        code = main(["--reference", SARS_GBK, "--feature-type", "any", "--gene", "N",
                     "--all", "--output", str(out)])
        assert code == 0
        entries = read_all(str(out))
        assert len(entries) == 2  # gene + CDS
        assert {len(s) for _, s in entries} == {1260}

    def test_a_missing_reference_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            main(["--reference", str(tmp_path / "absent.gbk"), "--gene", "N",
                  "--output", str(tmp_path / "x.fasta")])


@needs_genome
class TestAgainstTheMersRecord:
    """The same code path on a record with a compound (join) location."""

    def test_the_spliced_orf1ab_cds_is_a_whole_number_of_codons(self, tmp_path):
        out = tmp_path / "orf1ab.fasta"
        code = main(["--reference", GENOME_GBK, "--feature-type", "CDS",
                     "--gene", "orf1ab", "--all", "--output", str(out)])
        assert code == 0
        for header, seq in read_all(str(out)):
            assert len(seq) % 3 == 0, header
            assert seq.startswith("ATG"), header

    def test_the_frameshift_base_is_kept_twice(self, tmp_path):
        """Regression: merging the join() segments drops the slippage base.

        The ORF1ab CDS is join(a..b, b..c) -- the two segments share base b,
        which the ribosome reads twice. Extracting from the *merged* parts
        yields one base too few and shifts the frame of everything downstream.
        """
        from blast_utils import load_features

        feats = load_features(GENOME_GBK, ("CDS",))
        orf1ab = [
            f
            for f in feats["NC_019843.3"]
            if f.gene == "orf1ab" and len(f.sequence_parts) > 1
        ]
        assert orf1ab, "expected a compound ORF1ab CDS in the MERS record"
        feat = orf1ab[0]
        # The unmerged segments overlap, so they span one more base than the
        # merged ones.
        raw_span = sum(e - s + 1 for s, e in feat.sequence_parts)
        merged_span = sum(e - s + 1 for s, e in feat.parts)
        assert raw_span == merged_span + 1
        assert raw_span % 3 == 0 and merged_span % 3 != 0

    def test_a_plain_interval_from_the_mers_fasta(self, tmp_path):
        out = tmp_path / "region.fasta"
        assert main(["--reference", GENOME_FASTA, "--coords", "1-100",
                     "--output", str(out)]) == 0
        assert len(read_fasta(str(out))[1]) == 100
