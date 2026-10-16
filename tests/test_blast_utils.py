"""Unit tests for src/blast_utils.py: parsing, coordinates and GenBank features."""

import os

import pytest

from blast_utils import (
    DEFAULT_FEATURE_TYPES,
    OUTFMT_FIELDS,
    BlastError,
    BlastHit,
    Feature,
    _merge_intervals,
    best_hit_per_query,
    count_hits_per_query,
    detect_reference_format,
    genbank_to_fasta,
    load_features,
    overlapping_features,
)
from conftest import needs_blast, needs_gbk, needs_genome


def make_row(**over):
    """A tab-separated -outfmt 6 row, overridable per field."""
    values = {
        "qseqid": "q0",
        "sseqid": "NC_019843.3",
        "pident": "100.000",
        "length": "23",
        "mismatch": "0",
        "gapopen": "0",
        "qstart": "1",
        "qend": "23",
        "sstart": "100",
        "send": "122",
        "evalue": "1e-05",
        "bitscore": "46.1",
        "qlen": "23",
        "qcovhsp": "100",
    }
    values.update({k: str(v) for k, v in over.items()})
    return "\t".join(values[f] for f in OUTFMT_FIELDS)


class TestBlastHit:
    def test_parses_every_field(self):
        hit = BlastHit.from_row(make_row())
        assert (hit.query_id, hit.subject_id) == ("q0", "NC_019843.3")
        assert (hit.pident, hit.aln_length, hit.bitscore) == (100.0, 23, 46.1)
        assert hit.evalue == 1e-05

    def test_plus_strand_keeps_order(self):
        hit = BlastHit.from_row(make_row(sstart=100, send=122))
        assert (hit.strand, hit.ref_start, hit.ref_end) == ("+", 100, 122)

    def test_minus_strand_is_normalised_low_to_high(self):
        hit = BlastHit.from_row(make_row(sstart=122, send=100))
        assert (hit.strand, hit.ref_start, hit.ref_end) == ("-", 100, 122)

    def test_as_record_exposes_derived_fields(self):
        rec = BlastHit.from_row(make_row(sstart=122, send=100)).as_record()
        assert rec["strand"] == "-" and rec["ref_start"] == 100 and rec["ref_end"] == 122

    def test_wrong_column_count_is_an_error(self):
        with pytest.raises(BlastError, match="blastn columns"):
            BlastHit.from_row("q0\tNC_019843.3\t100.0")


class TestHitSelection:
    def test_best_hit_prefers_higher_bitscore(self):
        hits = [
            BlastHit.from_row(make_row(bitscore=30.2, length=15, qcovhsp=65)),
            BlastHit.from_row(make_row(bitscore=46.1, length=23)),
        ]
        assert best_hit_per_query(hits)["q0"].bitscore == 46.1

    def test_ties_on_bitscore_break_on_evalue(self):
        hits = [
            BlastHit.from_row(make_row(bitscore=46.1, evalue=1e-2)),
            BlastHit.from_row(make_row(bitscore=46.1, evalue=1e-9)),
        ]
        assert best_hit_per_query(hits)["q0"].evalue == 1e-9

    def test_keeps_queries_separate(self):
        hits = [
            BlastHit.from_row(make_row(qseqid="q0", bitscore=46.1)),
            BlastHit.from_row(make_row(qseqid="q1", bitscore=30.2)),
        ]
        best = best_hit_per_query(hits)
        assert set(best) == {"q0", "q1"}

    def test_counts_hits_per_query(self):
        hits = [BlastHit.from_row(make_row(qseqid=q)) for q in ("q0", "q0", "q1")]
        assert count_hits_per_query(hits) == {"q0": 2, "q1": 1}

    def test_empty_input(self):
        assert best_hit_per_query([]) == {} and count_hits_per_query([]) == {}


class TestDetectReferenceFormat:
    def test_detects_fasta_and_genbank_by_content(self, tmp_path, minimal_gbk):
        fasta = tmp_path / "x.fasta"
        fasta.write_text(">a\nACGT\n")
        assert detect_reference_format(str(fasta)) == "fasta"
        assert detect_reference_format(minimal_gbk) == "genbank"

    def test_content_beats_a_misleading_extension(self, tmp_path):
        """A FASTA named .gbk is still a FASTA."""
        path = tmp_path / "mislabelled.gbk"
        path.write_text(">a\nACGT\n")
        assert detect_reference_format(str(path)) == "fasta"

    def test_leading_blank_lines_are_skipped(self, tmp_path):
        path = tmp_path / "x.fa"
        path.write_text("\n\n>a\nACGT\n")
        assert detect_reference_format(str(path)) == "fasta"

    def test_extension_is_the_fallback_when_content_is_inconclusive(self, tmp_path):
        for name, expected in (("empty.gbk", "genbank"), ("empty.fna", "fasta")):
            path = tmp_path / name
            path.write_text("")
            assert detect_reference_format(str(path)) == expected

    def test_unrecognisable_file_raises_and_names_the_override(self, tmp_path):
        path = tmp_path / "guides.csv"
        path.write_text("Cas,Index\nRfxCas13d,0\n")
        with pytest.raises(ValueError, match="--ref-format"):
            detect_reference_format(str(path))

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            detect_reference_format(str(tmp_path / "nope.fasta"))


class TestMergeIntervals:
    def test_touching_intervals_merge(self):
        assert _merge_intervals([(1, 10), (11, 20)]) == [(1, 20)]

    def test_overlapping_intervals_merge(self):
        """The ORF1ab frameshift shares one base between its two arms."""
        assert _merge_intervals([(279, 13433), (13433, 21514)]) == [(279, 21514)]

    def test_disjoint_intervals_are_kept_apart(self):
        assert _merge_intervals([(21, 30), (1, 10)]) == [(1, 10), (21, 30)]


class TestFeature:
    def feature(self, **over):
        kw = dict(
            record_id="TEST",
            type="gene",
            start=100,
            end=200,
            strand="+",
            parts=((100, 200),),
        )
        kw.update(over)
        return Feature(**kw)

    def test_overlap_counts_shared_bases(self):
        assert self.feature().overlap_bp(150, 160) == 11

    def test_overlap_is_clipped_at_the_boundary(self):
        assert self.feature().overlap_bp(195, 220) == 6

    def test_no_overlap_is_zero(self):
        assert self.feature().overlap_bp(201, 230) == 0

    def test_single_base_overlap_counts(self):
        assert self.feature().overlap_bp(200, 230) == 1

    def test_overlap_skips_the_gap_in_a_compound_location(self):
        feat = self.feature(start=1, end=30, parts=((1, 10), (21, 30)))
        assert feat.overlap_bp(5, 25) == 11  # 5-10 plus 21-25, not 11-20

    def test_offset_is_one_based_from_the_feature_start(self):
        assert self.feature().offset_of(100, 122) == 1
        assert self.feature().offset_of(150, 172) == 51

    def test_offset_goes_negative_upstream_of_the_feature(self):
        """Flags a guide straddling the 5' boundary, as two PspCas13b guides do."""
        assert self.feature().offset_of(90, 112) == -9

    def test_offset_is_measured_from_the_other_end_on_the_minus_strand(self):
        assert self.feature(strand="-").offset_of(150, 160) == 41

    def test_label_falls_back_through_the_qualifiers(self):
        assert self.feature(gene="orf1ab").label == "orf1ab"
        assert self.feature(locus_tag="A_01").label == "A_01"
        assert self.feature(product="alpha protein").label == "alpha protein"
        assert self.feature(type="5'UTR").label == "5'UTR"


class TestLoadFeatures:
    def test_reads_qualifiers_and_strand(self, minimal_gbk):
        feats = load_features(minimal_gbk, ("gene", "CDS"))
        assert set(feats) == {"TEST.1"}
        by_label = {(f.type, f.label): f for f in feats["TEST.1"]}
        alpha = by_label[("gene", "alpha")]
        assert (alpha.start, alpha.end, alpha.strand) == (1, 30, "+")
        beta = by_label[("gene", "beta")]
        assert (beta.start, beta.end, beta.strand, beta.locus_tag) == (41, 60, "-", "B_01")
        cds = by_label[("CDS", "alpha")]
        assert cds.product == "alpha protein" and cds.locus_tag == "A_01"

    def test_compound_location_keeps_its_parts(self, minimal_gbk):
        cds = [f for f in load_features(minimal_gbk, ("CDS",))["TEST.1"]][0]
        assert cds.parts == ((1, 10), (21, 30))
        assert (cds.start, cds.end) == (1, 30)

    def test_feature_types_filter_is_honoured(self, minimal_gbk):
        assert [f.type for f in load_features(minimal_gbk, ("CDS",))["TEST.1"]] == ["CDS"]

    def test_empty_types_keeps_everything_but_source(self, minimal_gbk):
        types = [f.type for f in load_features(minimal_gbk, ())["TEST.1"]]
        assert "source" not in types
        assert sorted(types) == ["CDS", "gene", "gene", "mat_peptide"]

    def test_mat_peptide_is_excluded_by_default(self, minimal_gbk):
        types = {f.type for f in load_features(minimal_gbk, DEFAULT_FEATURE_TYPES)["TEST.1"]}
        assert "mat_peptide" not in types

    def test_features_come_back_in_positional_order(self, minimal_gbk):
        starts = [f.start for f in load_features(minimal_gbk, ())["TEST.1"]]
        assert starts == sorted(starts)


class TestOverlappingFeatures:
    def test_returns_each_overlapping_feature_with_its_overlap(self, minimal_gbk):
        feats = load_features(minimal_gbk, ())
        found = overlapping_features(feats, "TEST.1", 5, 8)
        assert {(f.type, bp) for f, bp in found} == {("gene", 4), ("CDS", 4), ("mat_peptide", 4)}

    def test_reports_nothing_in_an_intergenic_gap(self, minimal_gbk):
        assert overlapping_features(load_features(minimal_gbk, ()), "TEST.1", 33, 38) == []

    def test_unknown_record_id_is_empty_not_an_error(self, minimal_gbk):
        assert overlapping_features(load_features(minimal_gbk, ()), "ABSENT", 1, 10) == []

    def test_results_are_in_positional_order(self, minimal_gbk):
        found = overlapping_features(load_features(minimal_gbk, ()), "TEST.1", 1, 60)
        assert [f.start for f, _ in found] == sorted(f.start for f, _ in found)


class TestGenbankToFasta:
    def test_writes_one_record_named_by_accession(self, minimal_gbk, tmp_path):
        out = genbank_to_fasta(minimal_gbk, str(tmp_path / "out.fasta"))
        text = open(out).read()
        assert text.startswith(">TEST.1\n")
        assert len(text.split("\n")[1]) == 60

    def test_conversion_is_cached_but_force_rebuilds(self, minimal_gbk, tmp_path):
        out = genbank_to_fasta(minimal_gbk, str(tmp_path / "out.fasta"))
        os.utime(out, (0, 0))
        genbank_to_fasta(minimal_gbk, out)  # stale: mtime older than the gbk
        assert os.path.getmtime(out) > 0

    def test_record_without_a_sequence_is_rejected(self, tmp_path):
        path = tmp_path / "noseq.gbk"
        path.write_text(
            "LOCUS       EMPTY                     10 bp    DNA     linear   UNA 01-JAN-2024\n"
            "ACCESSION   EMPTY\nVERSION     EMPTY.1\n"
            "FEATURES             Location/Qualifiers\n     source          1..10\n"
            "ORIGIN\n//\n"
        )
        with pytest.raises(BlastError, match="no sequence|ORIGIN"):
            genbank_to_fasta(str(path), str(tmp_path / "out.fasta"))


@needs_blast
@needs_genome
class TestBlastIntegration:
    """Known-answer tests: pull a sequence out of the genome, blast it back."""

    def test_a_slice_of_the_genome_maps_to_where_it_came_from(self, genome_db, genome_seq):
        from blast_utils import blast_one

        start, end = 22604, 22633  # a real guide locus from the S gene
        hits = blast_one("probe", genome_seq[start - 1 : end], genome_db)
        best = best_hit_per_query(hits)["probe"]
        assert (best.ref_start, best.ref_end, best.strand) == (start, end, "+")
        assert best.pident == 100.0 and best.qcovhsp == 100.0

    def test_the_reverse_complement_maps_to_the_same_interval(self, genome_db, genome_seq):
        from blast_utils import blast_one

        start, end = 22604, 22633
        rc = genome_seq[start - 1 : end][::-1].translate(str.maketrans("ACGT", "TGCA"))
        best = best_hit_per_query(blast_one("probe", rc, genome_db))["probe"]
        assert (best.ref_start, best.ref_end, best.strand) == (start, end, "-")

    def test_one_call_handles_many_queries_independently(self, genome_db, genome_seq):
        from blast_utils import run_blastn

        loci = {"a": (4, 33), "b": (13441, 13470), "c": (27722, 27751)}
        records = [(k, genome_seq[s - 1 : e]) for k, (s, e) in loci.items()]
        best = best_hit_per_query(run_blastn(records, genome_db))
        assert {k: (h.ref_start, h.ref_end) for k, h in best.items()} == loci

    def test_a_sequence_absent_from_the_genome_returns_no_hits(self, genome_db):
        from blast_utils import blast_one

        assert blast_one("probe", "ACGT" * 8, genome_db) == []

    def test_blank_and_empty_inputs_skip_blast_entirely(self, genome_db):
        from blast_utils import blast_one, run_blastn

        assert run_blastn([], genome_db) == []
        assert blast_one("probe", "", genome_db) == []

    def test_duplicate_query_ids_are_rejected(self, genome_db):
        from blast_utils import run_blastn

        with pytest.raises(ValueError, match="unique"):
            run_blastn([("q0", "ACGTACGTACGTACGTACGTACG"), ("q0", "ACGT")], genome_db)

    def test_database_is_cached_across_calls(self, tmp_path):
        from blast_utils import make_blast_db
        from conftest import GENOME_FASTA

        db_dir = str(tmp_path / "db")
        first = make_blast_db(GENOME_FASTA, db_dir)
        mtime = os.path.getmtime(first + ".nin")
        assert make_blast_db(GENOME_FASTA, db_dir) == first
        assert os.path.getmtime(first + ".nin") == mtime  # not rebuilt

    def test_whitespace_in_the_db_dir_is_refused_with_advice(self, tmp_path):
        from blast_utils import make_blast_db
        from conftest import GENOME_FASTA

        with pytest.raises(BlastError, match="--db-dir"):
            make_blast_db(GENOME_FASTA, str(tmp_path / "has space"))

    def test_whitespace_in_the_reference_path_is_handled(self, tmp_path, genome_seq):
        """BLAST+ reads -in as a file list, so a spaced path needs the symlink route."""
        import shutil as sh

        from blast_utils import blast_one, make_blast_db
        from conftest import GENOME_FASTA

        spaced = tmp_path / "ref dir"
        spaced.mkdir()
        sh.copy(GENOME_FASTA, spaced / "NC_019843.3.fasta")
        db = make_blast_db(str(spaced / "NC_019843.3.fasta"), str(tmp_path / "db"))
        best = best_hit_per_query(blast_one("probe", genome_seq[22603:22633], db))["probe"]
        assert best.ref_start == 22604


@needs_blast
@needs_gbk
class TestGenbankReferenceIntegration:
    def test_the_real_record_parses_with_its_frameshift_cds(self):
        from conftest import GENOME_GBK

        feats = load_features(GENOME_GBK, ("gene", "CDS", "5'UTR", "3'UTR"))
        assert list(feats) == ["NC_019843.3"]
        orf1ab = [f for f in feats["NC_019843.3"] if f.type == "CDS" and f.gene == "orf1ab"]
        # join(279..13433,13433..21514) merges to one span; ORF1a is the second CDS.
        assert (21514, (279, 21514)) in [(f.end, (f.start, f.end)) for f in orf1ab]
        assert any(f.end == 13454 for f in orf1ab)

    def test_converted_fasta_matches_the_plain_fasta_sequence(self, tmp_path, genome_seq):
        from conftest import GENOME_GBK

        out = genbank_to_fasta(GENOME_GBK, str(tmp_path / "from_gbk.fasta"))
        converted = "".join(l.strip() for l in open(out) if not l.startswith(">")).upper()
        assert converted == genome_seq

    def test_known_loci_annotate_to_the_expected_genes(self):
        from conftest import GENOME_GBK

        feats = load_features(GENOME_GBK, ("gene", "CDS", "5'UTR", "3'UTR"))
        cases = {
            (4, 33): {"5'UTR"},
            (22604, 22633): {"S"},
            (27722, 27751): {"E"},
            (21477, 21499): {"orf1ab", "S"},  # the frameshift arm overlaps S
        }
        for (start, end), expected in cases.items():
            found = overlapping_features(feats, "NC_019843.3", start, end)
            genes = {f.gene for f, _ in found if f.gene} or {f.type for f, _ in found}
            assert genes == expected, (start, end)
