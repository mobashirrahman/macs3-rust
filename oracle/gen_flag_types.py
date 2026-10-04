#!/usr/bin/env python3
"""Regenerate `oracle/flag_types.tsv` from the pinned oracle's argparse.

The port's flag surface is derived from `oracle/flag_matrix.tsv` (which records
which flags exist, their arity and their choices) but that matrix carries no
*type* column. Type matters for accept/reject parity: argparse rejects
`--extsize notanumber` with exit 2, and a port that treats the value as a string
silently falls back to a default instead.

Rather than infer types from defaults -- `--verbose` defaults to `2` and is an
int, `--gsize` defaults to `hs` and is a string, `--pvalue` defaults to `1e-5`
and is a float -- they are read straight out of the `add_argument(...)` calls in
the oracle's own `bin/macs3`.

Types are scoped per subcommand. `-g`, `-q`, `-G` and `-w` are typed differently
in different subcommands: `callpeak -g` is a genome-size *string* while
`bdgdiff -g` is an int max-gap, so a single global flag->type map would reject
valid invocations.

Usage:
    python3 oracle/gen_flag_types.py [--check]

`--check` verifies the committed file matches, for CI.
"""
import argparse
import pathlib
import re
import sys

HEADER = """\
# GENERATED from the pinned oracle's argparse -- MACS3/bin/macs3 at
# c5443190e3edfeb301cc94acf450e2b2c026a223 (see oracle/ENV.lock).
# Regenerate with:  python3 oracle/gen_flag_types.py
#
# Lists only flags declared `type=int` or `type=float`. argparse rejects a
# malformed value for one of these with exit 2, so the port must too; a flag
# absent from this table is a string and accepts any value.
#
# Scoped per subcommand deliberately: `-g`, `-q`, `-G` and `-w` are typed
# differently in different subcommands (`callpeak -g` is a genome-size string,
# `bdgdiff -g` is an int max-gap), so a global flag->type map would be wrong.
subcommand\tflag\ttype"""


def oracle_bin() -> pathlib.Path:
    """Locate the oracle's `bin/macs3` via `oracle/ENV.lock`."""
    repo = pathlib.Path(__file__).resolve().parent.parent
    lock = repo / "oracle" / "ENV.lock"
    for line in lock.read_text().splitlines():
        if line.startswith("MACS3_PATH="):
            init = pathlib.Path(line.split("=", 1)[1])
            # MACS3_PATH is <oracle>/MACS3/__init__.py; bin/macs3 is <oracle>/bin/macs3
            cand = init.parent.parent / "bin" / "macs3"
            if cand.exists():
                return cand
    sys.exit("oracle/ENV.lock does not record a usable MACS3_PATH")


def regions(src: str):
    """Yield `(subcommand, source_text)` for each subparser's block.

    Most parsers are created inline as `subparsers.add_parser("name", ...)`;
    `pileup`'s is built by a helper that takes the name as a defaulted parameter,
    so the literal is recovered from that helper's signature.
    """
    marks = []
    for m in re.finditer(r'add_parser\(\s*(?:"([a-z0-9]+)"|name)', src):
        nm = m.group(1)
        if m.group(0).rstrip().endswith("name"):
            head = src[: m.start()]
            defm = None
            for d in re.finditer(r"def\s+add_\w*parser\s*\((.*?)\):", head, re.S):
                defm = d
            d = (
                re.search(r'name\s*=\s*"([a-z0-9]+)"', defm.group(1))
                if defm
                else None
            )
            nm = d.group(1) if d else None
        marks.append((m.start(), nm))
    for k, (pos, nm) in enumerate(marks):
        end = marks[k + 1][0] if k + 1 < len(marks) else len(src)
        if nm:
            yield nm, src[pos:end]


def table(src: str) -> dict:
    out = {}
    for nm, chunk in regions(src):
        # One segment per add_argument( call. Upstream writes every option spelling
        # on the same line as `add_argument(`, so that line carries the flags and the
        # remainder of the segment carries dest/type/default/help.
        idx = [m.start() for m in re.finditer(r"add_argument\(", chunk)]
        idx.append(len(chunk))
        for k in range(len(idx) - 1):
            seg = chunk[idx[k] : idx[k + 1]]
            tm = re.search(r"(?:^|[\s,(])type\s*=\s*(int|float)\b", seg)
            if not tm:
                continue
            first_line = seg.split("\n", 1)[0]
            for f in re.findall(r'"(-{1,2}[A-Za-z0-9][A-Za-z0-9-]*)"', first_line):
                out[(nm, f)] = tm.group(1)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    src = oracle_bin().read_text()
    t = table(src)
    body = "\n".join(f"{s}\t{f}\t{ty}" for (s, f), ty in sorted(t.items()))
    text = HEADER + "\n" + body + "\n"

    dest = pathlib.Path(__file__).resolve().parent / "flag_types.tsv"
    if a.check:
        cur = dest.read_text() if dest.exists() else ""
        if cur != text:
            print(f"{dest} is stale; re-run oracle/gen_flag_types.py", file=sys.stderr)
            return 1
        print(f"flag_types.tsv up to date ({len(t)} typed flags)")
        return 0
    dest.write_text(text)
    print(f"wrote {dest} ({len(t)} typed flags across "
          f"{len({s for s, _ in t})} subcommands)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
