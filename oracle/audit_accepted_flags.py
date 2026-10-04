#!/usr/bin/env python3
"""Audit: which accepted boolean flags have no observable effect?

`--cutoff-analysis` was accepted, exited 0, and wrote nothing while upstream wrote a
report (F195). Nothing in the golden corpus caught it, because `run_oracle.py` never
passes that flag -- the corpus only exercises flags someone thought to record.

This sweeps every `store_true` flag in the auto-derived matrix, runs the subcommand with
and without it on a small fixture, and reports flags whose output is **byte-identical**.
That is a candidate list, not a verdict: `--verbose` legitimately changes only stderr,
and some flags need inputs a single fixture does not have. The point is that "this flag
does nothing" becomes a checked, reviewable list instead of an unexamined assumption.

Usage:
    python3 oracle/audit_accepted_flags.py            # all subcommands
    python3 oracle/audit_accepted_flags.py callpeak   # one
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "target/release/macs3-rs")
MATRIX = os.path.join(REPO, "oracle", "flag_matrix.tsv")

# A minimal valid invocation per subcommand. Only commands that actually run on the
# fixtures below are swept; the rest report as "no base invocation".
SE_FIX = "tests/fixtures/se_basic/gauss_two_peaks"
PE_FIX = "tests/fixtures/pe_basic/gauss_fragments"
FRAG_FIX = "tests/fixtures/frag_basic/barcode_fragments"
BAM_FIX = "tests/fixtures/bam/reads.bam"
BEDG_FIX = "tests/fixtures/se_basic/gauss_two_peaks"

BASE: dict[str, list[str]] = {
    "callpeak": [
        "-t", f"{SE_FIX}/treat.bed", "-c", f"{SE_FIX}/ctrl.bed",
        "-f", "BED", "--nomodel", "--extsize", "200", "-g", "48200",
    ],
    # No `-B` in the base: it is one of the flags under test, and including it made the
    # sweep compare "-B" with "-B".
    "pileup": [
        "-i", f"{SE_FIX}/treat.bed", "-f", "BED", "--extsize", "200", "-o", "pileup.bdg",
    ],
    # `-o` is needed or the base writes nothing (filterdup defaults to stdout), and a
    # base that writes nothing makes every flag look like a no-op.
    "filterdup": [
        "-i", f"{SE_FIX}/treat.bed", "-f", "BED", "--keep-dup", "1", "--tsize", "200",
        "-o", "filtered.bed",
    ],
    # Only `callpeak`, `bdgpeakcall`, `filterdup`, `pileup` and `hmmratac` have
    # `store_true` flags in the matrix, so no other subcommand needs a base here.
    # (`predictd`/`randsample` have store-valued flags; their equivalence is covered by
    # oracle/check_predictd_randsample.sh instead.)
    # `-g` is bdgpeakcall's *maxgap*, not a genome size; there is no genome-size flag.
    "bdgpeakcall": ["-i", "@BDG@", "-c", "5.0", "-l", "200", "-g", "30",
                    "-o", "bdg_peaks.narrowPeak"],
    "bdgbroadcall": ["-i", "@BDG@", "-c", "10", "-l", "200", "-g", "30",
                     "-o", "bdg_peaks.broadPeak"],
}


def ensure_bedgraph(tmp: str) -> str:
    """A bedGraph built from the SE fixture by `macs3-rs pileup`, for the bedGraph
    commands' base invocations. Generating it from our own tool keeps the audit
    self-contained; it only needs *a* well-formed bedGraph to feed `bdgpeakcall`."""
    p = os.path.join(tmp, "in.bdg")
    if os.path.exists(p):
        return p
    subprocess.run(
        [BIN, "pileup", "-i", f"{SE_FIX}/treat.bed", "-f", "BED", "-B",
         "--extsize", "200", "-o", p, "--outdir", tmp],
        capture_output=True, text=True, timeout=600, check=False,
    )
    return p


def load_matrix() -> dict[str, list[tuple[str, str]]]:
    """subcommand -> [(flag, action)] for store_true flags."""
    with open(MATRIX) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    head = rows[0]
    ic, iflag, iaction = head.index("subcommand"), head.index("flag"), head.index("action")
    out: dict[str, list[tuple[str, str]]] = {}
    for r in rows[1:]:
        if len(r) <= max(ic, iflag, iaction):
            continue
        if r[iaction] == "store_true":
            out.setdefault(r[ic], []).append((r[iflag], r[iaction]))
    return out


# Flags whose value names an output path. Each is redirected into the per-run scratch
# directory, because several commands honour `-o` and ignore `--outdir` -- and a base run
# whose output lands outside the scanned tree makes every flag look like a no-op.
OUT_FLAGS = {
    "-o", "--ofile", "--outputfile", "--outdir", "--o-prefix", "--o",
    "--output", "--bed", "--faidx", "--model", "--barcode",
}


def run(args, outdir):
    os.makedirs(outdir, exist_ok=True)
    argv = []
    i = 0
    while i < len(args):
        a = args[i]
        if a in OUT_FLAGS and i + 1 < len(args):
            argv += [a, os.path.join(outdir, os.path.basename(args[i + 1]))]
            i += 2
            continue
        argv.append(a)
        i += 1
    p = subprocess.run(
        [BIN, *argv, "--outdir", outdir],
        capture_output=True,
        text=True,
        timeout=600,
    )
    p.stagedump_tail = p.stderr.strip().splitlines()[-2:]
    # Recursive: `--o-prefix` and friends create subdirectories.
    files = {}
    for dirpath, _dirs, names in os.walk(outdir):
        for f in sorted(names):
            path = os.path.join(dirpath, f)
            rel = os.path.relpath(path, outdir)
            with open(path, "rb") as fh:
                files[rel] = fh.read()
    return p.returncode, files, " | ".join(p.stagedump_tail)


def sweep(cmd: str, flags: list[str]) -> list[dict]:
    base = BASE.get(cmd)
    if base is None:
        return [{"flag": None, "note": "no base invocation defined"}]
    tmp = tempfile.mkdtemp(prefix="audit-")
    try:
        bdg = ensure_bedgraph(tmp)
        base = [bdg if a == "@BDG@" else a for a in base]
        rc0, files0, err0 = run([cmd, *base], os.path.join(tmp, "base"))
        if rc0 != 0:
            return [{"flag": None, "note": f"base invocation exited {rc0}: {err0}"}]
        out = []
        for flag, _ in flags:
            d = os.path.join(tmp, "f" + flag.strip("-").replace("-", "_"))
            rc, files, err = run([cmd, *base, flag], d)
            same_files = set(files) == set(files0)
            same_bytes = same_files and files == files0
            out.append(
                {
                    "flag": flag,
                    "rc": rc,
                    "new_files": sorted(set(files) - set(files0)),
                    "lost_files": sorted(set(files0) - set(files)),
                    "identical": bool(same_bytes),
                }
            )
        return out
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("subcommand", nargs="*")
    ap.add_argument("--json", help="write the full report here")
    args = ap.parse_args()

    if not os.path.exists(BIN):
        print(f"audit: {BIN} not built", file=sys.stderr)
        return 2

    matrix = load_matrix()
    cmds = args.subcommand or sorted(matrix)
    report = {}
    suspects = 0
    for cmd in cmds:
        flags = matrix.get(cmd, [])
        if not flags:
            continue
        rows = sweep(cmd, flags)
        report[cmd] = rows
        note = next((r["note"] for r in rows if r.get("note")), None)
        print(f"\n{cmd}: {len(flags)} boolean flag(s)")
        if note:
            print(f"  {note}")
            continue
        for r in rows:
            if r["identical"]:
                suspects += 1
                print(f"  IDENTICAL  {r['flag']}   (no observable effect)")
            else:
                why = []
                if r["new_files"]:
                    why.append("new " + ",".join(r["new_files"]))
                if r["lost_files"]:
                    why.append("lost " + ",".join(r["lost_files"]))
                if not why:
                    why.append("contents differ")
                print(f"  effective  {r['flag']}  ({'; '.join(why)})")
    print(f"\n{suspects} flag(s) with no observable effect")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())