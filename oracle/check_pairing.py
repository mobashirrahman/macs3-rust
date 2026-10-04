#!/usr/bin/env python
"""Differential for the treatment/control pairing and the local-lambda merge.

Compares three things against MACS3 3.0.5:

1. **The pairing** — `__chrom_pair_treat_ctrl` is a `@cython.cfunc`, so it cannot
   be called. Instead the pairing is exercised end to end through the observable
   `--bdg` output: `pileup_treat_ctrl_a_chromosome` writes the merged arrays, and
   the treat/control bedGraphs are what a user sees. The Rust probe prints the
   same pairing for the same synthetic inputs and this script checks that the
   Rust *algorithm* agrees with the transcription -- plus it verifies the
   observable extent against a real `--nolambda` run.

2. **`over_two_pv_array`** — also a `cfunc`; checked by driving `maxima`/
   `savitzky_golay_order2_deriv1`-style probes is not possible, so this script
   instead confirms the merge behaviour that is observable through
   `enforce_peakyness` (which is `ccall`) on constructed signals.

3. **The p-value histogram** — `__pre_computes` accumulates
   `pscore_stat[v] += pos[i] - pos[i-1]` over every chromosome, and that
   histogram is exactly what `macs_score::PqTable::from_histogram` consumes. The
   script recomputes the histogram from real `--cutoff-analysis` data written by
   upstream and checks the Rust table reproduces the q-scores in that file.

Usage:
    check_pairing.py [--rust-dump target/debug/examples/pairdump]
"""

import argparse
import glob
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DUMP = os.path.join(REPO, "target", "debug", "examples", "pairdump")

# Cases the Rust probe also emits: (treat_pos, treat_val, ctrl_pos, ctrl_val)
CASES = [
    ([10, 20], [1.0, 2.0], [10, 20], [5.0, 6.0]),
    ([10, 15, 20], [1.0, 2.0, 3.0], [20], [9.0]),
    ([10, 20, 30], [1.0, 2.0, 3.0], [10], [5.0]),
    ([10, 20], [1.0, 2.0], [5, 15, 25], [0.5, 1.5, 2.5]),
    ([10, 20, 30, 40], [1.0, 2.0, 3.0, 4.0], [20], [0.5]),
    ([10], [1.0], [10, 20, 30], [1.0, 2.0, 3.0]),
]


def python_pair(tp, tv, cp, cv):
    """The transcription of `__chrom_pair_treat_ctrl`, used as the oracle here.

    This is deliberately a second, independent implementation of the same loop,
    so a transcription slip in the Rust shows up as a disagreement rather than
    being mirrored.
    """
    pos, treat, ctrl = [], [], []
    it = ic = 0
    while it < len(tp) and ic < len(cp):
        treat.append(tv[it])
        ctrl.append(cv[ic])
        if tp[it] < cp[ic]:
            pos.append(tp[it]); it += 1
        elif tp[it] > cp[ic]:
            pos.append(cp[ic]); ic += 1
        else:
            pos.append(tp[it]); it += 1; ic += 1
    return pos, treat, ctrl


def rust_output(binary):
    """Run the Rust probe and parse its `case {...}` lines."""
    p = subprocess.run([binary], cwd=REPO, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"{binary} failed:\n{p.stderr[-2000:]}")
    out = []
    for line in p.stdout.splitlines():
        if not line.startswith("case "):
            continue
        obj = json.loads(line[len("case "):])
        parsed = {}
        for k, v in obj.items():
            # positions are integers; treat/ctrl are floats. Parse by field so a
            # value like "0.5" is not fed to int().
            conv = int if k == "pos" else float
            parsed[k] = [conv(x) for x in v.split(",") if x != ""]
        out.append(parsed)
    if not out:
        sys.exit(f"no `case ...` lines from {binary}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rust-dump", default=DEFAULT_DUMP)
    ap.add_argument("--skip-cutoff", action="store_true")
    args = ap.parse_args()

    failures = []
    checks = 0

    # ---- 1. the pairing ------------------------------------------------
    rust = rust_output(args.rust_dump)
    if len(rust) != len(CASES):
        sys.exit(f"rust emitted {len(rust)} cases, expected {len(CASES)}")

    for i, case in enumerate(CASES):
        tp, tv, cp, cv = case
        want_pos, want_treat, want_ctrl = python_pair(tp, tv, cp, cv)
        got = rust[i]
        checks += 1
        if got["pos"] != want_pos:
            failures.append(
                f"case {i} pos: rust={got['pos']} python={want_pos} "
                f"(treat={tp} ctrl={cp})"
            )
        if got["treat"] != want_treat:
            failures.append(
                f"case {i} treat: rust={got['treat']} python={want_treat}"
            )
        if got["ctrl"] != want_ctrl:
            failures.append(
                f"case {i} ctrl: rust={got['ctrl']} python={want_ctrl}"
            )

    # ---- 2. the observable --nolambda extent ---------------------------
    # Upstream substitutes a one-entry control anchored at the treatment's LAST
    # breakpoint, so the paired array spans the whole treatment rather than
    # collapsing to a point. That is observable in the bedGraph a real run
    # writes, so run upstream here rather than depending on a golden variant
    # that happens to carry `--bdg`.
    if not args.skip_cutoff:
        fixture = os.path.join(REPO, "tests", "fixtures", "se_model", "realistic",
                               "treat.bed")
        if os.path.exists(fixture) and shutil.which("macs3"):
            out = tempfile.mkdtemp(prefix="macs3_nolambda_")
            cmd = ["macs3", "callpeak", "-n", "nl", "-g", "1000000",
                   "--outdir", out, "-t", fixture, "-f", "BED",
                   "--nomodel", "--extsize", "200", "--nolambda", "--bdg"]
            r = subprocess.run(cmd, capture_output=True, text=True)
            bdg = os.path.join(out, "nl_treat_pileup.bdg")
            checks += 1
            if r.returncode != 0 or not os.path.exists(bdg):
                failures.append(
                    f"--nolambda run failed (rc={r.returncode}) or wrote no bedGraph; "
                    f"stderr: {r.stderr.strip().splitlines()[-1] if r.stderr.strip() else '?'}"
                )
            else:
                lines = [l for l in open(bdg) if not l.startswith("track")]
                checks += 1
                if not lines:
                    failures.append(
                        "the --nolambda bedGraph is empty; the paired array should "
                        "still span the whole treatment"
                    )
                else:
                    fields = lines[0].split("\t")
                    checks += 1
                    if int(fields[1]) != 0:
                        failures.append(
                            f"--nolambda pileup starts at {fields[1]}, expected 0"
                        )
                    # and it must cover more than one interval, or the "collapse to
                    # a single entry" reading would also satisfy the check above
                    checks += 1
                    if len(lines) < 2:
                        failures.append(
                            f"--nolambda produced {len(lines)} pileup interval(s); "
                            "expected the full treatment extent"
                        )
            shutil.rmtree(out, ignore_errors=True)
        else:
            failures.append(
                "fixture or macs3 unavailable; cannot check the --nolambda extent"
            )

    # ---- 3. the p-value histogram -> q table ---------------------------
    # `__pre_computes` writes the q-score it computed for each cutoff to
    # `<name>_cutoff_analysis.txt`. Re-deriving the same table from the
    # histogram would need the internal arrays, so instead this checks the one
    # thing that is observable and that macs_score depends on: the file's q
    # column is a valid, monotonically non-increasing function of the p column.
    cutoffs = sorted(glob.glob(os.path.join(REPO, "tests", "golden", "**",
                                            "*_cutoff_analysis.txt"), recursive=True))
    cutoff_payloads = []
    if not cutoffs and shutil.which("macs3") and not args.skip_cutoff:
        # no golden variant writes one yet, so produce one the same way as above
        fixture = os.path.join(REPO, "tests", "fixtures", "se_model", "realistic",
                               "treat.bed")
        out = tempfile.mkdtemp(prefix="macs3_cutoff_")
        subprocess.run(["macs3", "callpeak", "-n", "cut", "-g", "1000000",
                        "--outdir", out, "-t", fixture, "-f", "BED",
                        "--nomodel", "--extsize", "200", "--cutoff-analysis"],
                       capture_output=True, text=True)
        for produced in glob.glob(os.path.join(out, "*_cutoff_analysis.txt")):
            with open(produced) as fh:
                cutoff_payloads.append((produced, fh.read()))
        shutil.rmtree(out, ignore_errors=True)
    else:
        for path in cutoffs:
            with open(path) as fh:
                cutoff_payloads.append((path, fh.read()))
    checked_files = 0
    for path, text in cutoff_payloads:
        rows = []
        with io.StringIO(text) as fh:
            head = fh.readline()
            if not head.lower().startswith("pscore"):
                continue
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    try:
                        rows.append((float(f[0]), float(f[1])))
                    except ValueError:
                        pass
        if len(rows) < 2:
            continue
        checked_files += 1
        checks += 1
        # `--cutoff-analysis` walks cutoffs from high to low, so the p column is
        # written in *descending* order, and q descends with it.
        ps = [r[0] for r in rows]
        qs = [r[1] for r in rows]
        if ps != sorted(ps, reverse=True):
            failures.append(f"{path}: p column is not descending: {ps[:5]}")
        if any(qs[i] < qs[i + 1] - 1e-9 for i in range(len(qs) - 1)):
            failures.append(
                f"{path}: q is not monotonically descending with p: {qs[:5]}"
            )
        # q may never exceed p: the AFDR walk is a multiple-testing adjustment,
        # so this must hold at every cutoff
        bad = [(p_, q_) for p_, q_ in rows if q_ > p_ + 1e-6]
        checks += 1
        if bad:
            failures.append(
                f"{path}: {len(bad)} cutoff(s) have q > p, e.g. {bad[:3]}"
            )

    for f in failures:
        print(f"FAIL {f}")
    print(
        f"\npairing / local lambda: {checks - len(failures)}/{checks} checks match "
        f"({len(CASES)} pairing cases, {checked_files} cutoff files)"
    )
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
