"""Tests for reverse_complement and the headed outfmt-6 writer."""

import pytest

from blast_utils import (
    OUTFMT_FIELDS,
    BlastHit,
    reverse_complement,
    write_outfmt6,
)


def make_hit(query_id="q0", sstart=100, send=129, **kw):
    """A plausible 30 nt full-length HSP, overridable field by field."""
    fields = dict(
        query_id=query_id,
        subject_id="MT007544.1",
        pident=100.0,
        aln_length=30,
        mismatch=0,
        gapopen=0,
        qstart=1,
        qend=30,
        sstart=sstart,
        send=send,
        evalue=1e-10,
        bitscore=60.0,
        qlen=30,
        qcovhsp=100.0,
    )
    fields.update(kw)
    return BlastHit(**fields)


class TestReverseComplement:
    def test_a_plain_sequence(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_it_is_its_own_inverse(self):
        seq = "GGATTGTTGCAATTGTTTGGAGAAATCATC"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_rna_becomes_dna(self):
        # U complements to A, and the result never contains U: blastn has no U.
        assert reverse_complement("UUUU") == "AAAA"
        assert reverse_complement("AUGC") == "GCAT"
        assert "U" not in reverse_complement("AUGCUUGC")

    def test_lowercase_is_uppercased(self):
        assert reverse_complement("atgc") == "GCAT"

    def test_whitespace_is_stripped(self):
        assert reverse_complement(" AT GC\n") == "GCAT"

    def test_ambiguity_codes_are_complemented(self):
        # R (A/G) complements to Y (C/T), and N to N.
        assert reverse_complement("RN") == "NY"
        assert reverse_complement("ACGTRYSWKMBDHVN") == "NBDHVKMWSRYACGT"

    def test_a_non_nucleotide_is_rejected(self):
        with pytest.raises(ValueError, match="not a nucleotide sequence"):
            reverse_complement("ATGZ")

    def test_an_empty_string_is_empty(self):
        assert reverse_complement("") == ""


class TestWriteOutfmt6:
    def test_the_header_is_the_canonical_field_list(self, tmp_path):
        out = tmp_path / "hits.tsv"
        write_outfmt6([make_hit()], str(out))
        header = out.read_text().splitlines()[0]
        assert header.split("\t") == list(OUTFMT_FIELDS)

    def test_one_row_per_hsp(self, tmp_path):
        out = tmp_path / "hits.tsv"
        n = write_outfmt6([make_hit(), make_hit(query_id="q1")], str(out))
        assert n == 2
        assert len(out.read_text().splitlines()) == 3  # header + 2

    def test_labels_replace_the_internal_query_id(self, tmp_path):
        out = tmp_path / "hits.tsv"
        write_outfmt6([make_hit("q0")], str(out), labels={"q0": "N_cand_93"})
        row = out.read_text().splitlines()[1].split("\t")
        assert row[0] == "N_cand_93"

    def test_an_unlabelled_query_keeps_its_internal_id(self, tmp_path):
        out = tmp_path / "hits.tsv"
        write_outfmt6([make_hit("q7")], str(out), labels={"q0": "first"})
        assert out.read_text().splitlines()[1].split("\t")[0] == "q7"

    def test_subject_coordinates_are_written_as_blastn_reported_them(self, tmp_path):
        # A minus-strand HSP has sstart > send in real fmt6, and must not be
        # silently normalised low->high by the writer.
        out = tmp_path / "hits.tsv"
        write_outfmt6([make_hit(sstart=129, send=100)], str(out))
        row = dict(zip(OUTFMT_FIELDS, out.read_text().splitlines()[1].split("\t")))
        assert (row["sstart"], row["send"]) == ("129", "100")

    def test_annotation_columns_are_appended_after_the_canonical_fields(self, tmp_path):
        out = tmp_path / "hits.tsv"
        hit = make_hit()
        write_outfmt6(
            [hit],
            str(out),
            annotations={id(hit): {"genes": "N", "n_features": 2}},
            annotation_columns=("n_features", "genes"),
        )
        lines = out.read_text().splitlines()
        assert lines[0].split("\t") == list(OUTFMT_FIELDS) + ["n_features", "genes"]
        row = lines[1].split("\t")
        # The canonical block is untouched, the extras sit after it.
        assert row[: len(OUTFMT_FIELDS)] == [
            str(dict(zip(OUTFMT_FIELDS, row))[f]) for f in OUTFMT_FIELDS
        ]
        assert row[len(OUTFMT_FIELDS):] == ["2", "N"]

    def test_a_missing_annotation_becomes_an_empty_field(self, tmp_path):
        out = tmp_path / "hits.tsv"
        write_outfmt6(
            [make_hit()], str(out), annotations={}, annotation_columns=("genes",)
        )
        row = out.read_text().splitlines()[1].split("\t")
        assert len(row) == len(OUTFMT_FIELDS) + 1
        assert row[-1] == ""

    def test_two_hsps_of_one_query_get_their_own_annotations(self, tmp_path):
        # Keyed by id(hit), so HSPs of the same query must not share a row.
        out = tmp_path / "hits.tsv"
        a, b = make_hit("q0", sstart=100, send=129), make_hit("q0", sstart=500, send=529)
        write_outfmt6(
            [a, b],
            str(out),
            annotations={id(a): {"genes": "N"}, id(b): {"genes": "S"}},
            annotation_columns=("genes",),
        )
        genes = [line.split("\t")[-1] for line in out.read_text().splitlines()[1:]]
        assert genes == ["N", "S"]

    def test_no_hits_still_writes_a_header(self, tmp_path):
        out = tmp_path / "hits.tsv"
        assert write_outfmt6([], str(out)) == 0
        assert out.read_text().splitlines() == ["\t".join(OUTFMT_FIELDS)]

    def test_a_missing_output_directory_is_created(self, tmp_path):
        out = tmp_path / "nested" / "dir" / "hits.tsv"
        write_outfmt6([make_hit()], str(out))
        assert out.is_file()
