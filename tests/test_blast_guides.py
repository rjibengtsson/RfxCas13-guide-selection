"""End-to-end tests for scripts/blast_guides.py."""

import os

import pandas as pd
import pytest

from blast_guides import ANNOTATION_COLUMNS, HIT_COLUMNS, annotate, clean_sequence, main
from blast_utils import BlastHit, load_features
from conftest import GENOME_FASTA, GENOME_GBK, needs_blast, needs_gbk, needs_genome

# Real guide loci, taken from the head of data/mers_guides.csv.
LOCI = {0: (4, 33), 1: (1410, 1439), 4: (22604, 22633), 7: (27722, 27751)}
REGIONS = {0: "5'UTR", 1: "ORF1ab", 4: "S", 7: "E"}


class TestCleanSequence:
    def test_uppercases(self):
        assert clean_sequence("acgt") == "ACGT"

    def test_strips_internal_and_trailing_whitespace(self):
        assert clean_sequence(" ac gt\n") == "ACGT"

    def test_maps_rna_u_to_t(self):
        assert clean_sequence("ACGU") == "ACGT"

    def test_missing_values_become_empty(self):
        assert clean_sequence(None) == ""
        assert clean_sequence(float("nan")) == ""
        assert clean_sequence("") == ""


@needs_gbk
class TestAnnotate:
    def hit(self, start, end):
        row = "\t".join(
            ["q0", "NC_019843.3", "100.0", str(end - start + 1), "0", "0", "1",
             str(end - start + 1), str(start), str(end), "1e-5", "46.1",
             str(end - start + 1), "100"]
        )
        return BlastHit.from_row(row)

    @pytest.fixture
    def feats(self):
        return load_features(GENOME_GBK, ("gene", "CDS", "5'UTR", "3'UTR"))

    def test_parallel_columns_have_matching_lengths(self, feats):
        ann = annotate(feats, self.hit(22604, 22633))
        lengths = {
            len(ann[c].split(";"))
            for c in ("features", "feature_types", "feature_overlap_bp", "offset_in_feature")
        }
        assert lengths == {ann["n_features"]}

    def test_name_columns_are_deduplicated(self, feats):
        """A guide hits gene S and CDS S; 'genes' should say S once."""
        ann = annotate(feats, self.hit(22604, 22633))
        assert ann["genes"] == "S" and ann["n_features"] == 2

    def test_a_guide_in_two_genes_reports_both(self, feats):
        ann = annotate(feats, self.hit(21477, 21499))
        assert ann["genes"] == "orf1ab;S"

    def test_utr_features_have_no_gene_qualifier_but_do_appear(self, feats):
        ann = annotate(feats, self.hit(4, 33))
        assert ann["genes"] == "" and "5'UTR" in ann["features"]

    def test_a_clipped_guide_shows_partial_overlap_and_negative_offset(self, feats):
        ann = annotate(feats, self.hit(27576, 27605))  # E starts at 27590
        assert ann["feature_overlap_bp"] == "16;16"
        assert ann["offset_in_feature"] == "-13;-13"

    def test_all_annotation_columns_are_produced(self, feats):
        assert set(annotate(feats, self.hit(22604, 22633))) == set(ANNOTATION_COLUMNS)


@pytest.fixture
def guides_csv(tmp_path):
    """A small guides CSV whose expected coordinates are known."""
    path = tmp_path / "guides.csv"
    pd.DataFrame(
        {
            "Index": list(LOCI),
            "Region": [REGIONS[i] for i in LOCI],
            "Target Sequence": [None] * len(LOCI),
        }
    ).to_csv(path, index=False)
    return str(path)


def fill_sequences(csv_path, genome_seq, mutate=None):
    """Populate the sequence column from the genome at each known locus."""
    df = pd.read_csv(csv_path)
    df["Target Sequence"] = [genome_seq[s - 1 : e] for s, e in (LOCI[i] for i in df.Index)]
    if mutate:
        mutate(df)
    df.to_csv(csv_path, index=False)
    return df


def run_cli(*args):
    code = main([str(a) for a in args])
    assert code == 0, f"CLI exited {code}"


@needs_blast
@needs_genome
class TestCliWithFasta:
    def test_output_keeps_every_input_row_and_column(self, guides_csv, genome_seq, tmp_path):
        original = fill_sequences(guides_csv, genome_seq)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_FASTA, "--output", out)
        df = pd.read_csv(out)
        assert len(df) == len(original)
        assert list(original.columns) == list(df.columns)[: len(original.columns)]
        assert set(HIT_COLUMNS + ["n_hits", "blast_status"]) <= set(df.columns)

    def test_coordinates_match_the_known_loci(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_FASTA, "--output", out)
        df = pd.read_csv(out).set_index("Index")
        for idx, (start, end) in LOCI.items():
            assert (df.loc[idx, "ref_start"], df.loc[idx, "ref_end"]) == (start, end)
        assert (df.blast_status == "hit").all()
        assert (df.strand == "+").all()

    def test_no_annotation_columns_for_a_fasta_reference(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_FASTA, "--output", out)
        assert not set(ANNOTATION_COLUMNS) & set(pd.read_csv(out).columns)

    def test_rows_are_never_dropped_and_are_labelled(self, guides_csv, genome_seq, tmp_path):
        def mutate(df):
            df.loc[df.Index == 0, "Target Sequence"] = ""           # empty
            df.loc[df.Index == 1, "Target Sequence"] = "ACGT" * 8   # absent from genome
            df.loc[df.Index == 4, "Target Sequence"] = df.loc[
                df.Index == 4, "Target Sequence"
            ].str.lower().str.replace("T", "U")                      # lowercase RNA

        fill_sequences(guides_csv, genome_seq, mutate)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_FASTA, "--output", out)
        df = pd.read_csv(out).set_index("Index")
        assert len(df) == len(LOCI)
        assert df.loc[0, "blast_status"] == "empty_sequence"
        assert df.loc[1, "blast_status"] == "no_hit"
        assert df.loc[4, "blast_status"] == "hit"       # lowercase + U still maps
        assert df.loc[4, "ref_start"] == LOCI[4][0]
        assert df.loc[[0, 1], "n_hits"].tolist() == [0, 0]

    def test_guide_strand_input_gives_the_same_interval(self, guides_csv, genome_seq, tmp_path):
        """Blasting the reverse complement flips strand but not the interval."""
        df = pd.read_csv(guides_csv)
        df["Target Sequence"] = [
            genome_seq[s - 1 : e][::-1].translate(str.maketrans("ACGT", "TGCA"))
            for s, e in (LOCI[i] for i in df.Index)
        ]
        df.to_csv(guides_csv, index=False)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_FASTA, "--output", out)
        res = pd.read_csv(out).set_index("Index")
        assert (res.strand == "-").all()
        for idx, (start, end) in LOCI.items():
            assert (res.loc[idx, "ref_start"], res.loc[idx, "ref_end"]) == (start, end)

    def test_all_hits_file_lists_every_hsp(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        out, allh = tmp_path / "out.csv", tmp_path / "all.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_FASTA, "--output", out,
                "--all-hits", allh, "--id-column", "Index")
        hits = pd.read_csv(allh)
        assert len(hits) == pd.read_csv(out).n_hits.sum()
        assert set(hits.query_label) == set(LOCI)

    def test_min_pident_filter_can_discard_every_hit(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_FASTA, "--output", out, "--min-pident", "100.1")
        assert (pd.read_csv(out).blast_status == "no_hit").all()

    def test_a_bad_column_name_exits_two(self, guides_csv, genome_seq, tmp_path, capsys):
        fill_sequences(guides_csv, genome_seq)
        code = main(["--input", guides_csv, "--seq-column", "Nope",
                     "--reference", GENOME_FASTA, "--output", str(tmp_path / "out.csv")])
        assert code == 2
        assert "available columns" in capsys.readouterr().err

    def test_a_bad_id_column_exits_two(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        assert main(["--input", guides_csv, "--seq-column", "Target Sequence",
                     "--id-column", "Nope", "--reference", GENOME_FASTA,
                     "--output", str(tmp_path / "out.csv")]) == 2


@needs_blast
@needs_gbk
@needs_genome
class TestCliWithGenbank:
    def test_annotation_turns_on_automatically(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_GBK, "--output", out)
        df = pd.read_csv(out).set_index("Index")
        assert set(ANNOTATION_COLUMNS) <= set(df.columns)
        assert df.loc[4, "genes"] == "S" and df.loc[7, "genes"] == "E"
        assert (df.n_features > 0).all()

    def test_no_annotate_suppresses_the_columns(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_GBK, "--output", out, "--no-annotate")
        assert not set(ANNOTATION_COLUMNS) & set(pd.read_csv(out).columns)

    def test_genbank_coordinates_equal_the_fasta_coordinates(self, guides_csv, genome_seq, tmp_path):
        """The GenBank -> FASTA conversion must not shift anything."""
        fill_sequences(guides_csv, genome_seq)
        cols = ["ref_start", "ref_end", "strand"]
        results = []
        for name, ref in (("fa", GENOME_FASTA), ("gb", GENOME_GBK)):
            out = tmp_path / f"{name}.csv"
            run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                    "--reference", ref, "--output", out)
            results.append(pd.read_csv(out)[cols])
        assert results[0].equals(results[1])

    def test_feature_types_all_adds_mat_peptide(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        out = tmp_path / "out.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_GBK, "--output", out, "--feature-types", "all")
        df = pd.read_csv(out).set_index("Index")
        assert "mat_peptide" in df.loc[1, "feature_types"]

    def test_ref_format_override_is_respected(self, guides_csv, genome_seq, tmp_path):
        """Forcing 'fasta' on a GenBank file must not silently annotate."""
        fill_sequences(guides_csv, genome_seq)
        out = tmp_path / "out.csv"
        with pytest.raises(Exception):
            run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                    "--reference", GENOME_GBK, "--output", out, "--ref-format", "fasta")

    def test_annotated_all_hits_carries_the_feature_columns(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        out, allh = tmp_path / "out.csv", tmp_path / "all.csv"
        run_cli("--input", guides_csv, "--seq-column", "Target Sequence",
                "--reference", GENOME_GBK, "--output", out,
                "--all-hits", allh, "--id-column", "Index")
        assert set(ANNOTATION_COLUMNS) <= set(pd.read_csv(allh).columns)

    def test_a_csv_as_reference_is_refused_by_detection(self, guides_csv, genome_seq, tmp_path):
        fill_sequences(guides_csv, genome_seq)
        with pytest.raises(ValueError, match="--ref-format"):
            main(["--input", guides_csv, "--seq-column", "Target Sequence",
                  "--reference", guides_csv, "--output", str(tmp_path / "out.csv")])


@needs_blast
@needs_gbk
def test_the_real_guide_file_runs_end_to_end(tmp_path):
    """Smoke test over all 18k guides, if the real input is present."""
    from conftest import GUIDES_CSV

    if not os.path.isfile(GUIDES_CSV):
        pytest.skip(f"missing {GUIDES_CSV}")
    out = tmp_path / "out.csv"
    run_cli("--input", GUIDES_CSV, "--seq-column", "Target Sequence",
            "--reference", GENOME_GBK, "--output", out, "--threads", "4")
    df = pd.read_csv(out)
    assert len(df) == len(pd.read_csv(GUIDES_CSV))
    assert (df.blast_status == "hit").all()
    assert (df.n_features > 0).all()
    # Region labels and annotated genes must agree where both name a gene. The
    # Region must be among the genes listed, not necessarily first: guides in the
    # ORF1ab/S overlap are labelled Region=S but annotate as "orf1ab;S".
    named = df[df.genes.notna() & df.Region.isin(["S", "E", "M", "N"])]
    listed = [r.Region in r.genes.split(";") for r in named.itertuples()]
    assert all(listed), f"{listed.count(False)} rows disagree"
