"""Shared fixtures: put src/ and scripts/ on the path, and locate the real data."""

import os
import shutil
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for sub in ("src", "scripts"):
    path = os.path.join(ROOT, sub)
    if path not in sys.path:
        sys.path.insert(0, path)

GENOME_FASTA = os.path.join(ROOT, "data", "NC_019843", "NC_019843.3.fasta")
GENOME_GBK = os.path.join(ROOT, "data", "NC_019843", "NC_019843.3.gbk")
GUIDES_CSV = os.path.join(ROOT, "data", "mers_guides.csv")

# Integration tests need BLAST+ and the MERS reference; unit tests need neither.
needs_blast = pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="BLAST+ not on PATH",
)
needs_genome = pytest.mark.skipif(
    not os.path.isfile(GENOME_FASTA), reason=f"missing {GENOME_FASTA}"
)
needs_gbk = pytest.mark.skipif(not os.path.isfile(GENOME_GBK), reason=f"missing {GENOME_GBK}")


MINIMAL_GBK = """LOCUS       TEST                      60 bp    DNA     linear   UNA 01-JAN-2024
DEFINITION  synthetic record for tests.
ACCESSION   TEST
VERSION     TEST.1
FEATURES             Location/Qualifiers
     source          1..60
                     /organism="synthetic construct"
     gene            1..30
                     /gene="alpha"
     CDS             join(1..10,21..30)
                     /gene="alpha"
                     /locus_tag="A_01"
                     /product="alpha protein"
     gene            complement(41..60)
                     /gene="beta"
                     /locus_tag="B_01"
     mat_peptide     5..15
                     /gene="alpha"
                     /product="alpha peptide"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""


@pytest.fixture
def minimal_gbk(tmp_path):
    """A small GenBank file exercising join(), complement() and mat_peptide."""
    path = tmp_path / "test.gbk"
    path.write_text(MINIMAL_GBK)
    return str(path)


@pytest.fixture
def genome_seq():
    """The MERS genome as one uppercase string, for building known-answer queries."""
    if not os.path.isfile(GENOME_FASTA):
        pytest.skip(f"missing {GENOME_FASTA}")
    with open(GENOME_FASTA) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">")).upper()


@pytest.fixture
def genome_db(tmp_path):
    """A BLAST database built from the MERS genome FASTA."""
    from blast_utils import make_blast_db

    if shutil.which("makeblastdb") is None:
        pytest.skip("BLAST+ not on PATH")
    if not os.path.isfile(GENOME_FASTA):
        pytest.skip(f"missing {GENOME_FASTA}")
    return make_blast_db(GENOME_FASTA, str(tmp_path / "db"))
