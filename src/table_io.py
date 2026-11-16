"""Delimited-table reading, kept separate so the CLI never hardcodes a format.

The guide tables arriving from collaborators are a mix of comma- and
tab-separated files, so the delimiter is inferred from the extension rather
than assumed. :func:`read_table` is the only entry point the scripts need.
"""

from __future__ import annotations

import os

import pandas as pd

# Extension -> delimiter. Anything not listed falls back to sniffing the
# header line, so an unusual name still works without a flag.
SEPARATORS = {
    ".csv": ",",
    ".tsv": "\t",
    ".tab": "\t",
    ".txt": "\t",
}

# Spellings accepted for --input-sep, so the user never has to pass a raw tab.
SEP_ALIASES = {
    "auto": None,
    "comma": ",",
    "csv": ",",
    ",": ",",
    "tab": "\t",
    "tsv": "\t",
    "\t": "\t",
    "\\t": "\t",
    "semicolon": ";",
    ";": ";",
    "pipe": "|",
    "|": "|",
    "whitespace": r"\s+",
}


def resolve_separator(spec: str | None) -> str | None:
    """Turn an --input-sep spelling into a delimiter, or None to auto-detect."""
    if spec is None:
        return None
    key = spec.strip().lower()
    if key in SEP_ALIASES:
        return SEP_ALIASES[key]
    # A single literal character is taken as the delimiter itself.
    if len(spec) == 1:
        return spec
    raise ValueError(
        f"unrecognised --input-sep {spec!r}; use one of "
        f"{', '.join(sorted(k for k in SEP_ALIASES if k))} or a single character"
    )


def sniff_separator(path: str) -> str:
    """Infer the delimiter of *path* from its extension, then its header line.

    Extension first because it is the author's stated intent. When the
    extension says nothing, the header line decides: whichever candidate
    delimiter splits it into the most fields wins. Ties go to comma, so a
    single-column file keeps pandas' default behaviour.
    """
    suffix = os.path.splitext(path)[1].lower()
    if suffix in SEPARATORS:
        return SEPARATORS[suffix]

    if not os.path.isfile(path):
        raise FileNotFoundError(f"input table not found: {path}")

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.strip():
                header = line.rstrip("\n")
                break
        else:
            raise ValueError(f"{path} is empty; nothing to read")

    # Comma last in the candidate order so it wins ties via >= below.
    best, best_n = ",", 1
    for candidate in ("\t", ";", "|", ","):
        n = len(header.split(candidate))
        if n > best_n:
            best, best_n = candidate, n
    return best


def read_table(path: str, sep: str | None = None) -> pd.DataFrame:
    """Read a delimited table into a DataFrame.

    *sep* may be an --input-sep spelling ('tab', 'csv', ','...) or None to
    infer it. Values are left exactly as pandas parses them; sequence cleaning
    is the caller's job.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"input table not found: {path}")

    delimiter = resolve_separator(sep)
    if delimiter is None:
        delimiter = sniff_separator(path)

    # engine='python' because a regex delimiter (whitespace) needs it, and the
    # guide tables are far too small for the C parser's speed to matter.
    df = pd.read_csv(path, sep=delimiter, engine="python")
    if df.shape[1] == 1:
        # Almost always a delimiter mismatch rather than a genuine one-column
        # file, and it would otherwise surface as a confusing missing-column
        # error further down.
        print(
            f"warning: {path} parsed as a single column "
            f"({df.columns[0]!r}); check the delimiter"
        )
    return df
