#!/usr/bin/env python3
"""Replay every recorded golden `macs3` command with `macs3-rs` and compare bytes.

`tests/golden/<fixture>/<variant>/command.json` records the exact argv upstream was
run with, the exit status, and a SHA-256 for each output file. That is enough to
state the acceptance criterion directly:

    byte-identical *.xls, *_peaks.narrowPeak, *_summits.bed, *_model.r

`oracle/run_peak_e2e.sh` compares *coordinates*, which cannot see a score that moved
in the fourth decimal (F149). This script compares the files themselves.

The recorded argv is reused verbatim except for three substitutions, so nothing
here can quietly diverge from the corpus:

* `argv[0]` (the `macs3` interpreter) becomes the `macs3-rs` binary;
* `--outdir <path>` becomes a scratch directory, so the golden files are never
  written over;
* the recorded `--outdir` prefix is stripped from any path echoed back in the
  xls `# Command line:` header -- see `HEADER_REWRITE` below.

Usage:
    oracle/run_golden.py [--variant default] [--fixture SUBSTR] [--mechanism ...]
                         [--bin PATH] [--jobs N] [--verbose]
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN = os.path.join(ROOT, "tests", "golden")


def recorded_cases(variants: list[str]) -> list[tuple[str, str, dict]]:
    out = []
    for dirpath, _dirnames, filenames in os.walk(GOLDEN):
        if "command.json" not in filenames:
            continue
        rel = os.path.relpath(dirpath, GOLDEN)
        parts = rel.split(os.sep)
        # The corpus is `<group>/<fixture>/<variant>/command.json`; the fixture
        # key used everywhere else in this repo is `<group>/<fixture>`.
        if len(parts) != 3:
            continue
        group, name, variant = parts
        fixture = f"{group}/{name}"
        if variants and variant not in variants:
            continue
        with open(os.path.join(dirpath, "command.json")) as fh:
            out.append((fixture, variant, json.load(fh)))
    out.sort()
    return out


def rewrite_argv(argv: list[str], outdir: str) -> list[str]:
    """Swap the interpreter and the output directory out of a recorded argv."""
    new = list(argv)
    new[0] = "<BIN>"
    for i, tok in enumerate(new):
        if tok == "--outdir" and i + 1 < len(new):
            new[i + 1] = outdir
        elif tok.startswith("--outdir="):
            new[i] = f"--outdir={outdir}"
        elif "/tests/fixtures/" in tok:
            # Replay committed inputs at the current checkout prefix.
            new[i] = os.path.join(ROOT, "tests", "fixtures", tok.split("/tests/fixtures/", 1)[1])
    return new


def run_one(
    binary: str,
    fixture: str,
    variant: str,
    cfg: dict,
    mechanism: str | None,
    verbose: bool,
) -> tuple[str, int, int, int, list[str]]:
    """Return `(label, matched, total, exit_mismatch, notes)` for one case."""
    argv = cfg["command"]
    want_rc = int(cfg.get("returncode", 0))

    notes: list[str] = []
    golden_out = None
    argv_l = cfg["command"]
    for i, tok in enumerate(argv_l):
        if tok == "--outdir" and i + 1 < len(argv_l):
            golden_out = argv_l[i + 1]
        elif tok.startswith("--outdir="):
            golden_out = tok.split("=", 1)[1]
    if mechanism:
        # The mechanism is a harness-level switch, not a recorded flag: the two
        # `macs-callpeak-e2e` mechanisms cover both p-value branches upstream.
        if mechanism not in argv:
            notes.append(f"variant has no {mechanism} run recorded")

    scratch = tempfile.mkdtemp(prefix="golden-")
    try:
        concrete = rewrite_argv(argv, scratch)
        exe = [binary] + concrete[1:]
        proc = subprocess.run(exe, capture_output=True, text=True)
        if proc.returncode != want_rc:
            notes.append(
                f"exit {proc.returncode} != recorded {want_rc}: "
                f"{(proc.stderr or proc.stdout).strip().splitlines()[-1:] or ['']}"
            )
        files = cfg.get("files") or {}
        matched = 0
        total = 0
        for name in sorted(files):
            total += 1
            want_path = os.path.join(GOLDEN, fixture, variant, name)
            if not os.path.exists(want_path):
                notes.append(f"missing recorded reference {name}")
                continue
            got_path = os.path.join(scratch, name)
            if not os.path.exists(got_path):
                notes.append(f"missing {name}")
                continue
            a = open(want_path, "rb").read()
            b = open(got_path, "rb").read()
            if a == b:
                matched += 1
                continue
            # The xls echoes the command line, and the recorded one names the
            # interpreter and the golden outdir. Compare with those two lines
            # neutralised so the check is about the data, not the wrapper.
            if name.endswith(".xls") and _xls_equal_ignoring_paths(a, b, scratch, golden_out):
                matched += 1
                continue
            notes.append(_first_diff(a, b, name))
        return (
            f"{fixture}/{variant}",
            matched,
            total,
            1 if proc.returncode != want_rc else 0,
            notes,
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _xls_equal_ignoring_paths(a: bytes, b: bytes, scratch: str, golden_out: str | None) -> bool:
    """Compare two xls bodies with the run's own paths normalised away.

    Three header lines quote the invocation and therefore cannot match a replay:
    `# Command line:` (different interpreter *and* different argument order, since
    `--outdir` moves), and the `# ChIP-seq file = [...]` / `# control file = [...]`
    lines of the `# ARGUMENTS LIST` block, which echo whatever spelling of the
    path the caller used. Everything else -- including every count, `d`, scale
    factor and peak row -- is compared verbatim, so a real difference cannot hide
    behind this normalisation.
    """
    return _norm(a, scratch, golden_out) == _norm(b, scratch, golden_out)


def _norm(blob: bytes, scratch: str, golden_out: str | None) -> list[bytes]:
    out = []
    for line in blob.splitlines():
        if line.startswith(b"# Command line:"):
            continue
        text = line.decode(errors="replace")
        text = text.replace(scratch, "<OUTDIR>")
        if golden_out:
            text = text.replace(golden_out, "<OUTDIR>")
        text = text.replace(ROOT + "/", "<ROOT>/")
        out.append(text.encode())
    return out


# F245: peak files have fixed, *named* columns, and a bare "line 21 differs" hides
# which one. `fold_enrichment` sits between `-log10(pvalue)` and `-log10(qvalue)` in
# the xls and carries the *same value* as narrowPeak's `signalValue`, so reading a
# differing line positionally produces a plausible but wrong conclusion -- which is
# exactly what F243 did before F244 retracted it. Name the column.
_XLS_COLUMNS = ("chr", "start", "end", "length", "abs_summit", "pileup",
                "-log10(pvalue)", "fold_enrichment", "-log10(qvalue)", "name")
_NARROWPEAK_COLUMNS = ("chr", "start", "end", "name", "score", "strand",
                       "signalValue", "pValue", "qValue", "peak")


def _columns_for(name: str) -> tuple[str, ...]:
    if name.endswith(".xls"):
        return _XLS_COLUMNS
    if name.endswith(".narrowPeak"):
        return _NARROWPEAK_COLUMNS
    return ()


def _named_diff(x: bytes, y: bytes, name: str) -> str:
    """`columns=fold_enrichment` when the differing fields can be named."""
    cols = _columns_for(name)
    if not cols:
        return ""
    xf, yf = x.decode(errors="replace").rstrip("\n").split("\t"), y.decode(errors="replace").rstrip("\n").split("\t")
    if len(xf) != len(cols) or len(yf) != len(cols):
        return ""
    hit = [cols[j] for j in range(len(cols)) if xf[j] != yf[j]]
    return f" columns={','.join(hit)}" if hit else ""


def _first_diff(a: bytes, b: bytes, name: str) -> str:
    la, lb = a.splitlines(), b.splitlines()
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            if x.startswith(b"#") and y.startswith(b"#"):
                continue
            return (f"{name}:{i + 1}{_named_diff(x, y, name)}"
                    f" golden={x.decode(errors='replace')[:120]}"
                    f" got={y.decode(errors='replace')[:120]}")
    if len(la) != len(lb):
        extra = la[len(lb):] or lb[len(la):]
        return f"{name}: {len(la)} vs {len(lb)} lines, first extra={extra[0].decode(errors='replace')[:120]}"
    return f"{name}: bytes differ"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", action="append", default=[])
    ap.add_argument("--fixture", default=None)
    ap.add_argument("--bin", default=os.path.join(ROOT, "target", "release", "macs3-rs"))
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    cases = recorded_cases(args.variant)
    if args.fixture:
        cases = [c for c in cases if args.fixture in c[0]]
    if args.limit:
        cases = cases[: args.limit]

    print(f"golden cases: {len(cases)}  binary: {args.bin}")
    tot_match = tot_files = tot_rc = 0
    full = 0
    failures: list[tuple[str, list[str]]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {
            ex.submit(run_one, args.bin, fx, va, cfg, None, args.verbose): (fx, va)
            for fx, va, cfg in cases
        }
        for fut in concurrent.futures.as_completed(futs):
            label, matched, total, rc_bad, notes = fut.result()
            tot_match += matched
            tot_files += total
            tot_rc += rc_bad
            if total and matched == total and rc_bad == 0:
                full += 1
            elif notes:
                failures.append((label, notes))

    print()
    # F273: the denominator used to be every recorded case, which included the 481
    # invalid-invocation cases that legitimately produce **no output files** (425 argparse
    # rejections at rc 2, 56 runtime rejections at rc 1). `full` requires `total > 0`, so
    # those cases could never be counted and the summary read `7204/7685` -- 93.7% -- for a
    # corpus that was actually at parity. Reported separately, with their exit statuses
    # checked, so the number means what it looks like it means.
    with_files = sum(1 for f, v, c in cases if (c.get("files") or {}))
    no_files = len(cases) - with_files
    print(f"cases with output    : {with_files} (byte-identical: {full}/{with_files})")
    print(f"output files matched : {tot_match}/{tot_files}")
    print(f"invalid-invocation   : {no_files} cases, exit-status mismatch {tot_rc}")
    if with_files and full != with_files:
        print(f"NOTE: {with_files - full} case(s) with output files are not byte-identical")
    if failures:
        print()
        print(f"{len(failures)} case(s) with differences:")
        for label, notes in sorted(failures):
            print(f"  {label}")
            for n in notes[:4]:
                print(f"    {n}")
            if len(notes) > 4:
                print(f"    ... {len(notes) - 4} more")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
