#!/usr/bin/env python3
"""Check the release criteria that are about *behaviour*, not output bytes.

Two of the acceptance criteria are not observable in the golden outputs:

* **Identical accept/reject behaviour on invalid invocations.** Upstream's argparse
  exits `2` on a usage error; a runtime error exits `1`. Anything else -- including
  `0`, and including a panic's `101` -- is a divergence.
* **Errors raised before any output file is created.** Upstream validates names,
  genome size and control-file readability *before* opening the writer, so a
  rejected run leaves no partial output behind. A run that creates
  ``foo_peaks.xls`` and then fails is observably different, because the file exists.

Both are cheap to check and expensive to notice later, so they get their own script
and their own CI job.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

BIN = os.environ.get("MACS3_RS", "target/release/macs3-rs")

SUBCOMMANDS = [
    "callpeak", "bdgpeakcall", "bdgbroadcall", "bdgcmp", "bdgopt", "cmbreps",
    "bdgdiff", "filterdup", "predictd", "pileup", "randsample", "refinepeak",
    "callvar", "hmmratac",
]


def run(args, cwd=None):
    return subprocess.run(
        [os.path.abspath(BIN), *args], capture_output=True, text=True, cwd=cwd
    )


def outputs_of(d):
    """Every file `macs3-rs` may have written into `d`."""
    return sorted(os.listdir(d)) if os.path.isdir(d) else []


class Failure(Exception):
    pass


def expect(condition, message):
    if not condition:
        raise Failure(message)


def check_usage_errors(tmp):
    """Usage errors exit 2, and the exit code is never a panic's 101."""
    cases = [
        (["--no-such-flag"], "unknown global flag"),
        (["callpeak"], "missing required -t/-g"),
        (["not_a_subcommand"], "unknown subcommand"),
        (["callpeak", "--help-there"], "unknown option"),
        (["pileup"], "missing -i/-f"),
    ]
    for args, what in cases:
        p = run(args, cwd=tmp)
        expect(
            p.returncode == 2,
            f"usage error '{what}' ({' '.join(args)}) exited {p.returncode}, expected 2",
        )
    # every subcommand exists and answers --help with 0
    for cmd in SUBCOMMANDS:
        p = run([cmd, "--help"], cwd=tmp)
        expect(p.returncode == 0, f"`{cmd} --help` exited {p.returncode}, expected 0")


def check_no_panics(tmp):
    """Malformed and hostile input must never produce a Rust panic (exit 101)."""
    bed = os.path.join(tmp, "bad.bed")
    with open(bed, "wb") as fh:
        fh.write(
            b"chr1\n"
            b"chr1\tnot_a_number\t200\tr\t0\t+\n"
            b"chr1\t100\t200\tr\t0\t*\n"
            b"chr1\t100\t200\tr\t0\t+\n"
            b"\n"
            b"\t\t\t\n"
            b"chr1\t99999999999999999999\t1\tr\t0\t+\n"
        )
    for args in (
        ["callpeak", "-n", "p", "-g", "1000", "-t", bed, "-f", "BED",
         "--nomodel", "--extsize", "50"],
        ["filterdup", "-i", bed, "-f", "BED", "-o", os.path.join(tmp, "o.bed")],
        ["predictd", "-i", bed, "-m", "50", "50"],
    ):
        p = run(args, cwd=tmp)
        expect(
            p.returncode not in (101, -6, -11),
            f"{' '.join(args[:2])} crashed with exit {p.returncode}",
        )
        expect(
            "panicked at" not in p.stderr and "RUST_BACKTRACE" not in p.stderr,
            f"{' '.join(args[:2])} panicked:\n{p.stderr[-800:]}",
        )


def check_errors_precede_output(tmp):
    """A rejected run must leave no output file behind."""
    bed = os.path.join(tmp, "ok.bed")
    with open(bed, "wb") as fh:
        for i in range(200):
            fh.write(b"chr1\t%d\t%d\tr%d\t0\t+\n" % (100 + i * 20, 150 + i * 20, i))

    # a control file that does not exist is a runtime error (exit 1), not a panic,
    # and must be detected before the writer opens anything
    outdir = os.path.join(tmp, "out_missing_ctrl")
    p = run(["callpeak", "-n", "x", "-g", "1000", "-t", bed,
             "-c", os.path.join(tmp, "no_such_control.bed"), "-f", "BED",
             "--nomodel", "--extsize", "50", "--outdir", outdir], cwd=tmp)
    expect(p.returncode == 1, f"missing control exited {p.returncode}, expected 1")
    expect(not outputs_of(outdir),
           f"missing-control run left {outputs_of(outdir)} behind; errors must "
           f"precede any output file")

    # a treatment file with chromosome names disagreeing with the control is also a
    # runtime rejection (F27, `check_names`), and equally must precede any output
    ctrl = os.path.join(tmp, "other.bed")
    with open(ctrl, "wb") as fh:
        for i in range(200):
            fh.write(b"chrZ\t%d\t%d\tc%d\t0\t+\n" % (100 + i * 20, 150 + i * 20, i))
    outdir2 = os.path.join(tmp, "out_name_mismatch")
    p = run(["callpeak", "-n", "x", "-g", "1000", "-t", bed, "-c", ctrl,
             "-f", "BED", "--nomodel", "--extsize", "50", "--outdir", outdir2],
            cwd=tmp)
    expect(p.returncode == 1, f"chrom mismatch exited {p.returncode}, expected 1")
    expect(not outputs_of(outdir2),
           f"chrom-mismatch run left {outputs_of(outdir2)} behind")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", help="write a machine-readable report here")
    args = ap.parse_args()

    if not os.path.exists(BIN):
        print(f"check_exit_contract: {BIN} not built", file=sys.stderr)
        return 2

    tmp = tempfile.mkdtemp(prefix="exit-contract-")
    checks = {
        "usage_errors_exit_2": check_usage_errors,
        "no_panics_on_malformed_input": check_no_panics,
        "errors_precede_output_files": check_errors_precede_output,
    }
    report = {}
    failed = 0
    try:
        for name, fn in checks.items():
            try:
                fn(tmp)
                report[name] = "pass"
                print(f"pass  {name}")
            except Failure as exc:
                report[name] = f"fail: {exc}"
                failed += 1
                print(f"FAIL  {name}: {exc}", file=sys.stderr)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
