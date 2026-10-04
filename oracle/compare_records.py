#!/usr/bin/env python
"""Record-level differential: upstream MACS3 parsers vs the Rust port.

Compares, for every fixture input, the records each side parsed. Both dumps are
sorted on (chrom, pos, strand) as tuples so the comparison is order-independent
but still position-exact: a plain text sort would order "100" before "49" and
report differences that do not exist.

Usage:
    compare_records.py [--dump-bin PATH] [FIXTURE_ROOT]
"""

import argparse
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORACLE_DUMP = os.path.join(REPO, "oracle", "dump_oracle_records.py")
RUST_DUMP = os.path.join(REPO, "target", "debug", "macs-io-dump")

# extension -> format name understood by both dump programs
FORMATS = {
    ".bed": "bed",
    ".bed.gz": "bed",
    ".bedpe": "bedpe",
    ".frag": "frag",
    ".bedGraph": "bedgraph",
    ".bedgraph": "bedgraph",
}


def fmt_for(path):
    for ext, fmt in sorted(FORMATS.items(), key=lambda kv: -len(kv[0])):
        if path.endswith(ext):
            return fmt
    return None


def sort_key(line):
    """Sort each column so integers compare numerically.

    Sorting on raw text orders "100" before "49" and would report differences that
    do not exist. It also fails to tie-break on the trailing count column of a
    FRAG record, which leaves the relative order of equal-coordinate records up to
    the sort and shows up as a pile of false mismatches. Comparing every numeric
    column as a number fixes both.
    """
    parts = line.rstrip("\n").split("\t")
    key = []
    for part in parts:
        # (0, int) sorts before (1, str), keeping the key totally ordered
        key.append((0, int(part), "") if part.lstrip("+-").isdigit() else (1, 0, part))
    return tuple(key)


def load(text):
    return sorted((l for l in text.splitlines() if l), key=sort_key)


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr.strip()


def compare_one(path, rust_bin):
    fmt = fmt_for(path)
    if fmt is None:
        return None

    oc, oout, oerr = run([sys.executable, ORACLE_DUMP, fmt, path])
    rc, rout, rerr = run([rust_bin, fmt, path])

    if oc != 0:
        return ("oracle-error", f"{fmt}: upstream failed: {oerr.splitlines()[-1] if oerr else '?'}")
    if rc != 0:
        return ("rust-error", f"{fmt}: {rerr.splitlines()[-1] if rerr else '?'}")

    a, b = load(oout), load(rout)
    if a == b:
        return ("ok", f"{fmt}: {len(a)} records identical")

    # localise: first differing index, and how many records differ overall
    n = max(len(a), len(b))
    first = next((i for i in range(n) if (a[i] if i < len(a) else None) != (b[i] if i < len(b) else None)), None)
    diff = sum(1 for i in range(n) if (a[i] if i < len(a) else None) != (b[i] if i < len(b) else None))
    ctx = ""
    if first is not None and first < len(a) and first < len(b):
        ctx = f"\n    oracle: {a[first]}\n    rust:   {b[first]}"
    return (
        "mismatch",
        f"{fmt}: {diff} of {n} records differ (oracle {len(a)}, rust {len(b)});"
        f" first at index {first}{ctx}",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-bin", default=RUST_DUMP)
    ap.add_argument("root", nargs="?", default=os.path.join(REPO, "tests", "fixtures"))
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    targets = []
    for dirpath, _dirnames, filenames in os.walk(args.root):
        for name in sorted(filenames):
            p = os.path.join(dirpath, name)
            if fmt_for(p):
                targets.append(p)
    targets.sort()

    if not targets:
        sys.exit(f"no parsable inputs under {args.root}")

    tally = {}
    for path in targets:
        result = compare_one(path, args.dump_bin)
        if result is None:
            continue
        status, msg = result
        tally[status] = tally.get(status, 0) + 1
        if status != "ok" or args.verbose:
            rel = os.path.relpath(path, args.root)
            print(f"[{status:>11}] {rel}: {msg}")

    ok = tally.get("ok", 0)
    bad = sum(v for k, v in tally.items() if k != "ok")
    print(f"\n{ok} files identical, {bad} files differing, {len(targets)} inputs total")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
