"""Reusable blastn helpers, plus GenBank feature annotation.

The core routine is :func:`run_blastn`, which blasts many sequences in a single
blastn process. :func:`blast_one` is a thin wrapper for the single-sequence case.

:func:`detect_reference_format` sniffs whether a reference is FASTA or GenBank.
For GenBank references, :func:`genbank_to_fasta` produces the FASTA that BLAST+
needs and :func:`load_features` reads the feature table, so hits can be
annotated with :func:`overlapping_features`.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, asdict
from typing import Iterable, Sequence

# Explicit -outfmt 6 field list, in order. Keep in sync with BlastHit.from_row.
OUTFMT_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "qcovhsp",
)

# Defaults tuned for ~30 nt guide RNAs; every one is overridable.
DEFAULT_TASK = "blastn-short"
DEFAULT_WORD_SIZE = 11
DEFAULT_EVALUE = 1e-3
DEFAULT_MAX_TARGET_SEQS = 50


class BlastError(RuntimeError):
    """Raised when a BLAST+ executable is missing or exits non-zero."""


@dataclass(frozen=True)
class BlastHit:
    """One HSP, with coordinates normalised to always run low -> high."""

    query_id: str
    subject_id: str
    pident: float
    aln_length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    qcovhsp: float

    @classmethod
    def from_row(cls, row: str) -> "BlastHit":
        f = row.rstrip("\n").split("\t")
        if len(f) != len(OUTFMT_FIELDS):
            raise BlastError(
                f"expected {len(OUTFMT_FIELDS)} blastn columns, got {len(f)}: {row!r}"
            )
        return cls(
            query_id=f[0],
            subject_id=f[1],
            pident=float(f[2]),
            aln_length=int(f[3]),
            mismatch=int(f[4]),
            gapopen=int(f[5]),
            qstart=int(f[6]),
            qend=int(f[7]),
            sstart=int(f[8]),
            send=int(f[9]),
            evalue=float(f[10]),
            bitscore=float(f[11]),
            qlen=int(f[12]),
            qcovhsp=float(f[13]),
        )

    @property
    def strand(self) -> str:
        """'+' or '-', from the orientation of the subject coordinates."""
        return "-" if self.sstart > self.send else "+"

    @property
    def ref_start(self) -> int:
        """1-based start on the reference, always <= ref_end."""
        return min(self.sstart, self.send)

    @property
    def ref_end(self) -> int:
        return max(self.sstart, self.send)

    def as_record(self) -> dict:
        rec = asdict(self)
        rec.update(strand=self.strand, ref_start=self.ref_start, ref_end=self.ref_end)
        return rec


def _require(exe: str) -> str:
    path = shutil.which(exe)
    if path is None:
        raise BlastError(f"{exe} not found on PATH; install BLAST+ (e.g. conda install -c bioconda blast)")
    return path


def make_blast_db(fasta: str, db_dir: str, force: bool = False) -> str:
    """Build a nucleotide BLAST database for *fasta* inside *db_dir*.

    Returns the database prefix. Rebuilds only when the database is missing,
    older than the FASTA, or *force* is set.
    """
    if not os.path.isfile(fasta):
        raise FileNotFoundError(f"reference FASTA not found: {fasta}")
    if any(c.isspace() for c in db_dir):
        raise BlastError(
            f"BLAST cannot read a database path containing whitespace: {db_dir!r}; "
            "pass a whitespace-free --db-dir"
        )
    os.makedirs(db_dir, exist_ok=True)
    db_prefix = os.path.join(db_dir, "_".join(os.path.basename(fasta).split()))

    index = db_prefix + ".nin"
    if not force and os.path.exists(index) and os.path.getmtime(index) >= os.path.getmtime(fasta):
        return db_prefix

    # BLAST+ parses -in as a space-separated list of files, so a path containing
    # whitespace has to reach it through a whitespace-free symlink.
    fasta_in = os.path.abspath(fasta)
    if any(c.isspace() for c in fasta_in):
        fasta_in = db_prefix + ".source.fasta"
        if os.path.islink(fasta_in) or os.path.exists(fasta_in):
            os.unlink(fasta_in)
        os.symlink(os.path.abspath(fasta), fasta_in)

    proc = subprocess.run(
        [_require("makeblastdb"), "-in", fasta_in, "-dbtype", "nucl", "-out", db_prefix],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise BlastError(f"makeblastdb failed:\n{proc.stderr.strip()}")
    return db_prefix


def _write_fasta(records: Sequence[tuple[str, str]], handle) -> int:
    written = 0
    for seq_id, seq in records:
        if not seq:
            continue
        handle.write(f">{seq_id}\n{seq}\n")
        written += 1
    handle.flush()
    return written


def run_blastn(
    records: Iterable[tuple[str, str]],
    db: str,
    task: str = DEFAULT_TASK,
    word_size: int = DEFAULT_WORD_SIZE,
    evalue: float = DEFAULT_EVALUE,
    max_target_seqs: int = DEFAULT_MAX_TARGET_SEQS,
    threads: int = 1,
    dust: str = "no",
    extra_args: Sequence[str] | None = None,
) -> list[BlastHit]:
    """Blast every ``(query_id, sequence)`` in *records* against *db*.

    All sequences go into one blastn process, so query IDs must be unique and
    free of whitespace. Returns every HSP reported, in blastn's output order.
    """
    records = [(str(i), s) for i, s in records]
    if not records:
        return []

    ids = [i for i, _ in records]
    if len(set(ids)) != len(ids):
        raise ValueError("query IDs passed to run_blastn must be unique")

    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
        query_path = fh.name
        n = _write_fasta(records, fh)
    try:
        if n == 0:
            return []
        cmd = [
            _require("blastn"),
            "-query", query_path,
            "-db", db,
            "-task", task,
            "-word_size", str(word_size),
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_target_seqs),
            "-num_threads", str(threads),
            "-dust", dust,
            "-outfmt", "6 " + " ".join(OUTFMT_FIELDS),
        ]
        if extra_args:
            cmd.extend(extra_args)

        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise BlastError(f"blastn failed:\n{proc.stderr.strip()}")
        return [BlastHit.from_row(line) for line in proc.stdout.splitlines() if line.strip()]
    finally:
        os.unlink(query_path)


def blast_one(seq_id: str, sequence: str, db: str, **kwargs) -> list[BlastHit]:
    """Blast a single sequence. Convenience wrapper over :func:`run_blastn`."""
    return run_blastn([(seq_id, sequence)], db, **kwargs)


def best_hit_per_query(hits: Iterable[BlastHit]) -> dict[str, BlastHit]:
    """Pick one hit per query: highest bitscore, then lowest e-value, then longest."""
    best: dict[str, BlastHit] = {}
    for hit in hits:
        current = best.get(hit.query_id)
        if current is None or (
            (hit.bitscore, -hit.evalue, hit.aln_length)
            > (current.bitscore, -current.evalue, current.aln_length)
        ):
            best[hit.query_id] = hit
    return best


def count_hits_per_query(hits: Iterable[BlastHit]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for hit in hits:
        counts[hit.query_id] = counts.get(hit.query_id, 0) + 1
    return counts


# --------------------------------------------------------------------------
# Reference format detection and GenBank annotation
# --------------------------------------------------------------------------

GENBANK_SUFFIXES = (".gb", ".gbk", ".gbf", ".gbff", ".genbank")
FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".seq")

# mat_peptide is deliberately excluded: in a coronavirus record it adds the full
# nsp list to nearly every ORF1ab row. Opt in via --feature-types.
DEFAULT_FEATURE_TYPES = ("gene", "CDS", "5'UTR", "3'UTR")


def detect_reference_format(path: str) -> str:
    """Return 'genbank' or 'fasta' by sniffing the first non-blank line.

    Content wins over the file name; the suffix is only consulted when the
    content is inconclusive.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"reference not found: {path}")

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("LOCUS"):
                return "genbank"
            break

    suffix = os.path.splitext(path)[1].lower()
    if suffix in GENBANK_SUFFIXES:
        return "genbank"
    if suffix in FASTA_SUFFIXES:
        return "fasta"
    raise ValueError(
        f"cannot tell whether {path} is FASTA or GenBank: it starts with neither "
        "'>' nor 'LOCUS'. Pass --ref-format to say which it is."
    )


def _require_biopython():
    try:
        from Bio import SeqIO  # noqa: PLC0415
    except ImportError as exc:  # pragma: no cover
        raise BlastError(
            "reading a GenBank reference needs Biopython (conda install -c conda-forge biopython)"
        ) from exc
    return SeqIO


def genbank_to_fasta(gbk_path: str, out_fasta: str, force: bool = False) -> str:
    """Write the sequences of a GenBank file to *out_fasta* for BLAST+.

    Records are named by their accession/version, so blastn's subject IDs match
    the keys returned by :func:`load_features`. Skips the conversion when the
    FASTA is already newer than the GenBank file.
    """
    if (
        not force
        and os.path.exists(out_fasta)
        and os.path.getmtime(out_fasta) >= os.path.getmtime(gbk_path)
    ):
        return out_fasta

    SeqIO = _require_biopython()
    os.makedirs(os.path.dirname(os.path.abspath(out_fasta)), exist_ok=True)
    n = 0
    with open(out_fasta, "w") as out:
        for record in SeqIO.parse(gbk_path, "genbank"):
            # A record with no ORIGIN block parses fine but has undefined
            # contents, which only raises on access.
            try:
                seq = str(record.seq)
            except Exception as exc:
                raise BlastError(
                    f"{gbk_path} record {record.id} carries no sequence; a GenBank "
                    "reference must include its ORIGIN block"
                ) from exc
            if not seq or set(seq.upper()) <= {"N"}:
                raise BlastError(
                    f"{gbk_path} record {record.id} carries no sequence; a GenBank "
                    "reference must include its ORIGIN block"
                )
            out.write(f">{record.id}\n{seq}\n")
            n += 1
    if n == 0:
        raise BlastError(f"no records parsed from {gbk_path}")
    return out_fasta


def _merge_intervals(parts: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge touching/overlapping 1-based inclusive intervals."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(parts):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


@dataclass(frozen=True)
class Feature:
    """A GenBank feature, with 1-based inclusive coordinates.

    ``parts`` holds the merged segments of a compound (``join(...)``) location,
    such as the ribosomal-frameshift ORF1ab CDS; ``start``/``end`` span them all.
    """

    record_id: str
    type: str
    start: int
    end: int
    strand: str
    parts: tuple[tuple[int, int], ...]
    gene: str = ""
    locus_tag: str = ""
    product: str = ""

    @property
    def label(self) -> str:
        return self.gene or self.locus_tag or self.product or self.type

    def overlap_bp(self, start: int, end: int) -> int:
        """Bases shared with the 1-based inclusive interval [start, end]."""
        return sum(
            max(0, min(end, p_end) - max(start, p_start) + 1) for p_start, p_end in self.parts
        )

    def offset_of(self, start: int, end: int) -> int:
        """1-based position of the interval within this feature, respecting strand.

        Measured along the feature's span, so it is approximate for a compound
        location with introns.
        """
        return start - self.start + 1 if self.strand != "-" else self.end - end + 1


def load_features(
    gbk_path: str, feature_types: Sequence[str] = DEFAULT_FEATURE_TYPES
) -> dict[str, list[Feature]]:
    """Read *gbk_path* into ``{record_id: [Feature, ...]}``, sorted by position.

    Pass ``feature_types=()`` to keep every type except ``source``.
    """
    SeqIO = _require_biopython()
    wanted = {t.lower() for t in feature_types}
    features: dict[str, list[Feature]] = {}

    for record in SeqIO.parse(gbk_path, "genbank"):
        collected: list[Feature] = []
        for feat in record.features:
            if feat.type.lower() == "source":
                continue
            if wanted and feat.type.lower() not in wanted:
                continue
            parts = _merge_intervals(
                [(int(p.start) + 1, int(p.end)) for p in feat.location.parts]
            )
            if not parts:
                continue
            q = feat.qualifiers
            collected.append(
                Feature(
                    record_id=record.id,
                    type=feat.type,
                    start=parts[0][0],
                    end=parts[-1][1],
                    strand="-" if feat.location.strand == -1 else "+",
                    parts=tuple(parts),
                    gene=q.get("gene", [""])[0],
                    locus_tag=q.get("locus_tag", [""])[0],
                    product=q.get("product", [""])[0],
                )
            )
        collected.sort(key=lambda f: (f.start, f.end, f.type, f.label))
        features[record.id] = collected

    return features


def overlapping_features(
    features: dict[str, list[Feature]], record_id: str, start: int, end: int
) -> list[tuple[Feature, int]]:
    """Features on *record_id* sharing at least one base with [start, end].

    Returns ``(feature, overlap_bp)`` in positional order. A linear scan: GenBank
    records hold few enough features that an interval index would not pay off.
    """
    found = []
    for feat in features.get(record_id, ()):
        if feat.start > end:
            break
        overlap = feat.overlap_bp(start, end)
        if overlap > 0:
            found.append((feat, overlap))
    return found
