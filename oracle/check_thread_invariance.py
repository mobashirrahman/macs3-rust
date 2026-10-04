#!/usr/bin/env python3
"""Acceptance gate: byte-identical output at 1 thread and at 32 threads.

The porting plan requires "Identical output at --threads 1 and --threads 32". The pinned oracle
(MACS3 3.0.5, commit c5443190) has **no `--threads` flag** -- `macs3 callpeak --threads 4` fails with
`error: unrecognized arguments: --threads 4`, and the auto-derived `oracle/flag_matrix.tsv` contains
no such row. Accepting the flag would therefore *break* drop-in compatibility and the
"identical accept/reject behaviour on invalid invocations" criterion, so the criterion is read as
what it can only mean: the number of worker threads must not be observable in the output.

Our knob is the `MACS3_RS_THREADS` environment variable (`configure_pool` in
`crates/macs-cli/src/commands/callpeak.rs`), which is why this is an env var and not a flag.

For each recorded case this runs the command twice -- once pinned to a single worker, once to 32 --
into separate scratch directories and compares every output file byte-for-byte. The `.xls` files
echo their own command line and outdir, so those are compared with the same normalisation
`run_golden.py` uses; without it every case would "differ" on the wrapper and the gate would be
vacuous.

Rayon is only used for per-chromosome work, so the interesting failures are chromosome-ordering
races, shared-cache mutation, and non-deterministic accumulation order -- all of which this catches
only if the corpus actually contains multi-chromosome inputs, which it does.

Usage:
    oracle/check_thread_invariance.py [--jobs N] [--limit N] [--threads 32] [--verbose]
Exit status is 0 only when every compared file is identical.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import run_golden as rg  # noqa: E402


def run_once(binary: str, argv: list[str], threads: int) -> tuple[str, int]:
    scratch = tempfile.mkdtemp(prefix=f"tinv{threads}-")
    env = dict(os.environ, MACS3_RS_THREADS=str(threads))
    try:
        concrete = rg.rewrite_argv(argv, scratch)
        proc = subprocess.run(
            [binary] + concrete[1:], capture_output=True, text=True, env=env
        )
        return scratch, proc.returncode
    except Exception:
        shutil.rmtree(scratch, ignore_errors=True)
        raise


def check_case(binary: str, fixture: str, variant: str, threads: int) -> list[str]:
    """Return a list of human-readable differences for one case."""
    cpath = os.path.join(rg.GOLDEN, fixture, variant, "command.json")
    cfg = json.load(open(cpath))
    argv = cfg["command"]
    want_rc = cfg.get("returncode", 0)
    files = sorted((cfg.get("files") or {}).keys())

    dir_a, rc_a = run_once(binary, argv, 1)
    dir_b, rc_b = run_once(binary, argv, threads)
    notes: list[str] = []
    try:
        if rc_a != rc_b:
            notes.append(f"exit status differs: 1 thread -> {rc_a}, {threads} thread(s) -> {rc_b}")
        if rc_a != want_rc:
            notes.append(f"exit {rc_a} != recorded {want_rc} (1 thread)")
        if rc_b != want_rc:
            notes.append(f"exit {rc_b} != recorded {want_rc} ({threads} threads)")
        for name in files:
            pa, pb = os.path.join(dir_a, name), os.path.join(dir_b, name)
            if not os.path.exists(pa) and not os.path.exists(pb):
                continue
            if not os.path.exists(pa):
                notes.append(f"{name}: produced only at {threads} threads")
                continue
            if not os.path.exists(pb):
                notes.append(f"{name}: produced only at 1 thread")
                continue
            a, b = open(pa, "rb").read(), open(pb, "rb").read()
            if a == b:
                continue
            if name.endswith(".xls") and rg._xls_equal_ignoring_paths(a, b, dir_a, dir_b):
                continue
            notes.append(f"{name}: {rg._first_diff(a, b, name)}")
    finally:
        shutil.rmtree(dir_a, ignore_errors=True)
        shutil.rmtree(dir_b, ignore_errors=True)
    return notes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=os.path.join(rg.ROOT, "target", "release", "macs3-rs"))
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--limit", type=int, default=0, help="0 = every recorded case")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    # `command.json` lives at GOLDEN/<group...>/<fixture>/<variant>/command.json, and the
    # grouping depth is not fixed (`sweep/...`, `se_model/...`, but also `bam/...`), so walk
    # for the marker rather than assuming a depth. An assumed two-level walk silently
    # compared zero cases and still printed a clean summary -- which is exactly the
    # failure mode this gate exists to prevent.
    cases: list[tuple[str, str]] = []
    for dirpath, _dirnames, filenames in os.walk(rg.GOLDEN):
        if "command.json" not in filenames:
            continue
        rel = os.path.relpath(dirpath, rg.GOLDEN)
        parts = rel.split(os.sep)
        if len(parts) < 2:
            continue
        cases.append(("/".join(parts[:-1]), parts[-1]))
    cases.sort()
    # Deterministic stride sampling keeps CI time bounded without biasing toward the
    # alphabetically-first fixtures, which are the small hand-written ones.
    if args.limit and args.limit < len(cases):
        stride = len(cases) / args.limit
        cases = [cases[int(i * stride)] for i in range(args.limit)]

    if not os.path.exists(args.binary):
        print(f"binary not found: {args.binary}", file=sys.stderr)
        return 2

    if not cases:
        print(f"no recorded cases found under {rg.GOLDEN}", file=sys.stderr)
        return 2

    from concurrent.futures import ThreadPoolExecutor

    bad: list[tuple[str, list[str]]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {
            pool.submit(check_case, args.binary, f, v, args.threads): (f, v) for f, v in cases
        }
        done = 0
        for fut, key in futs.items():
            notes = fut.result()
            done += 1
            if notes:
                bad.append((f"{key[0]}/{key[1]}", notes))
            if args.verbose and done % 200 == 0:
                print(f"  {done}/{len(cases)} ...", file=sys.stderr)

    print(f"cases compared        : {len(cases)}")
    if not cases:
        return 2
    print(f"thread counts         : 1 vs {args.threads}")
    print(f"cases differing       : {len(bad)}")
    for label, notes in bad[:20]:
        print(f"  {label}")
        for n in notes[:4]:
            print(f"    {n}")
    if len(bad) > 20:
        print(f"  ... and {len(bad) - 20} more")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
