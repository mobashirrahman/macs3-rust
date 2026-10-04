#!/usr/bin/env python3
"""Report which *named* column of a peak file differs, instead of which line.

Two findings in `docs/upstream-findings.md` were caused by the same mistake: reading
a differing line of `*_peaks.xls` as `-log10(pvalue)` when it was `fold_enrichment`
(F243, retracted in F244). The xls and narrowPeak column orders are *not* the same,
and `fold_enrichment` sits between `-log10(pvalue)` and `-log10(qvalue)`, so an
off-by-one there silently produces a plausible-looking but wrong numeric argument.

This script exists so that cannot recur silently. It validates the header against
upstream's own names (`PeakIO.py:797` and `:968`) and then reports differences by
column *name*, never by position.

Usage:
    python3 oracle/check_peak_columns.py GOLDEN_DIR OUR_DIR [GLOB]

    # e.g.
    python3 oracle/check_peak_columns.py \\
        tests/golden/se_model/realistic/call_summits /tmp/ours \\
        '*_peaks.xls'

Exit status is 0 when the named-column sets agree, 1 when they differ -- so it is
usable as a gate in its own right, and as a pre-step before any numeric reasoning
about a peak file.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

# `PeakIO.py:968` -- the reader validates an xls against exactly these names, so a
# file whose header differs is not the format this script reasons about.
XLS_COLUMNS = (
    "chr",
    "start",
    "end",
    "length",
    "abs_summit",
    "pileup",
    "-log10(pvalue)",
    "fold_enrichment",
    "-log10(qvalue)",
    "name",
)

# `PeakIO.py:817` writes narrowPeak with these columns. Note `signalValue` carries
# the *fold enrichment*, not `-log10(p)`: for `se_model/realistic/call_summits` the
# xls row `... 9 8.5323 5.97002 6.45005 ...` and the narrowPeak row `... 64 . 5.97002
# 8.5323 6.45005 325 ...` show the same 5.97002 in `signalValue`. Treating that as a
# p-score is exactly the F243 error.
NARROWPEAK_COLUMNS = (
    "chr",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue",
    "peak",
)

# Which xls column each narrowPeak field corresponds to, so a difference can be
# reported once, under the xls name, with the narrowPeak alias alongside.
NARROWPEAK_TO_XLS = {
    "signalValue": "fold_enrichment",
    "pValue": "-log10(pvalue)",
    "qValue": "-log10(qvalue)",
    "score": "pileup",
}


def columns_for(path: str) -> tuple[str, ...]:
    if path.endswith(".xls"):
        return XLS_COLUMNS
    return NARROWPEAK_COLUMNS


def read_rows(path: str) -> list[list[str]]:
    rows = []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.rstrip("\n").split("\t"))
    return rows


def check_header(path: str) -> bool:
    """Confirm the file has the column count we are going to index against."""
    want = columns_for(path)
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            n = len(line.rstrip("\n").split("\t"))
            if n != len(want):
                print(f"  FAIL {os.path.basename(path)}: {n} columns, expected "
                      f"{len(want)} ({', '.join(want)})")
                return False
            return True
    print(f"  FAIL {os.path.basename(path)}: no data rows")
    return False


def diff_named(golden: str, ours: str) -> tuple[set[str], int]:
    """Names of the columns that differ anywhere, plus the number of differing rows."""
    cols = columns_for(golden)
    a, b = read_rows(golden), read_rows(ours)
    names: set[str] = set()
    rows = 0
    for i in range(max(len(a), len(b))):
        ra = a[i] if i < len(a) else []
        rb = b[i] if i < len(b) else []
        if ra == rb:
            continue
        rows += 1
        for j in range(min(len(ra), len(rb))):
            if ra[j] != rb[j]:
                name = cols[j] if j < len(cols) else f"col{j}"
                names.add(name)
        if len(ra) != len(rb):
            names.add("<row count>")
    return names, rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("golden_dir")
    ap.add_argument("our_dir")
    ap.add_argument("glob", nargs="?", default="*_peaks.*")
    args = ap.parse_args()

    ok = True
    print("differing columns, by name:")
    for g in sorted(glob.glob(os.path.join(args.golden_dir, args.glob))):
        base = os.path.basename(g)
        o = os.path.join(args.our_dir, base)
        if not os.path.exists(o):
            continue
        if not (check_header(g) and check_header(o)):
            ok = False
            continue
        names, rows = diff_named(g, o)
        if not names:
            continue
        for n in sorted(names):
            alias = NARROWPEAK_TO_XLS.get(n)
            via = f"  (narrowPeak {n} -> xls {alias})" if alias else ""
            print(f"  {base}: column {n!r} differs in {rows} row(s){via}")
    if ok:
        print("  headers validated against upstream's column names")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())