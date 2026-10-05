#!/usr/bin/env python
"""Differential for the `--cutoff-analysis` pre-computation.

`__pre_computes` is a `@cython.cdef`, so it cannot be called directly. This
script therefore checks `macs_peaks`'s cutoff accounting against an **independent
transcription** of the same loop, written from
`MACS3/Signal/CallPeakUnit.py` and deliberately structured differently (dense
loop over indices rather than cursors) so that a transcription slip in the Rust
shows up as a disagreement rather than being mirrored.

It also pins the two numeric facts that decide the answer and that are easy to get
backwards:

* the cutoff ladder is `numpy.float64` (from `np.arange`), so `f32 score > f64
  cutoff` promotes the **score** to f64 -- `f32(0.3)` clears the `0.3` cutoff
  (F36);
* `above_cutoff` uses `pos_array[above_cutoff - 1]`, so the first run starts at 0
  (F5).

Usage:
    check_precomputes.py [--rust-dump target/debug/examples/precompdump]
"""

import argparse
import os
import subprocess
import sys

# NumPy lives in the oracle virtualenv, so this re-execs under the provisioned
# interpreter rather than whichever `python3` happens to be first on PATH.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oracle_env import ensure_oracle_python  # noqa: E402

ensure_oracle_python()

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DUMP = os.path.join(REPO, "target", "debug", "examples", "precompdump")

MAX_GAP = 50
MIN_LENGTH = 200


def ladder():
    return [round(x, 5) for x in sorted(list(np.arange(0.3, 10.0, 0.3)), reverse=True)]


def transcript_counts(pos, score, max_gap=MAX_GAP, min_length=MIN_LENGTH):
    """An independent transcription of `__pre_computes`'s per-cutoff loop.

    Written as a plain indexed loop rather than the upstream cursor walk, and
    using numpy for the threshold, so it exercises a different route to the same
    numbers.
    """
    cutoffs = ladder()
    pos = np.asarray(pos, dtype="i4")
    score = np.asarray(score, dtype="f4")
    npeaks = []
    lengths = []
    for cutoff in cutoffs:
        total_l = 0
        total_p = 0
        # numpy value-based comparison against an f64 scalar promotes the array
        above = np.nonzero(score > cutoff)[0]
        if above.size == 0:
            npeaks.append(0)
            lengths.append(0)
            continue
        endpos = pos[above]
        startpos = pos[above - 1] if above[0] != 0 else np.where(above == 0, 0, pos[above - 1])
        # the first chunk's start is forced to 0
        if above[0] == 0:
            startpos = startpos.copy()
            startpos[0] = 0

        region = [(int(startpos[0]), int(endpos[0]))]
        lastp = int(endpos[0])
        for i in range(1, len(startpos)):
            ts = int(startpos[i])
            te = int(endpos[i])
            tl = ts - lastp
            if tl <= max_gap:
                region.append((ts, te))
            else:
                plen = region[-1][1] - region[0][0]
                if plen >= min_length:
                    total_l += plen
                    total_p += 1
                region = [(ts, te)]
            lastp = te
        if region:
            plen = region[-1][1] - region[0][0]
            if plen >= min_length:
                total_l += plen
                total_p += 1
        npeaks.append(total_p)
        lengths.append(total_l)
    return npeaks, lengths


CASES = {
    "single_run": ([0, 500], [0.0, 5.0]),
    "two_merged": ([0, 200, 400], [0.0, 5.0, 5.0]),
    "two_split": ([0, 200, 900, 1100], [0.0, 5.0, 0.0, 5.0]),
    "below_min": ([0, 100], [0.0, 5.0]),
    "exactly_min": ([0, 200], [0.0, 5.0]),
    "empty": ([], []),
    "flat_at_03": ([0, 1000, 2000], [0.3, 0.3, 0.3]),
    "flat_at_99": ([0, 1000, 2000], [9.9, 9.9, 9.9]),
    "all_above": ([0, 500, 1000, 1500], [5.0, 6.0, 7.0, 8.0]),
    "mismatch": ([0, 200], [5.0]),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rust-dump", default=DEFAULT_DUMP)
    args = ap.parse_args()

    p = subprocess.run([args.rust_dump], cwd=REPO, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"{args.rust_dump} failed:\n{p.stderr[-2000:]}")
    lines = [l for l in p.stdout.splitlines() if l.strip()]
    if not lines:
        sys.exit(f"{args.rust_dump} printed nothing")

    rust_ladder = [float(x) for x in lines[0].split(",")]
    want_ladder = ladder()
    failures = []
    checks = 1

    if rust_ladder != want_ladder:
        i = next(
            (j for j, (a, b) in enumerate(zip(want_ladder, rust_ladder)) if a != b),
            min(len(want_ladder), len(rust_ladder)),
        )
        failures.append(
            f"cutoff ladder differs at index {i}: "
            f"numpy={want_ladder[i:i+3]} rust={rust_ladder[i:i+3]} "
            f"(len {len(want_ladder)} vs {len(rust_ladder)})"
        )

    rust = {}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        name, field, vals = parts
        rust.setdefault(name, {})[field] = [int(x) for x in vals.split(",")]

    for name, (pos, score) in CASES.items():
        if name not in rust:
            failures.append(f"{name}: missing from the Rust dump")
            continue
        want_n, want_l = transcript_counts(pos, score)
        checks += 2
        if rust[name]["npeaks"] != want_n:
            failures.append(
                f"{name} npeaks: rust={rust[name]['npeaks']}\n"
                f"          transcription={want_n}"
            )
        if rust[name]["length"] != want_l:
            failures.append(
                f"{name} length: rust={rust[name]['length']}\n"
                f"           transcription={want_l}"
            )

    # The two numeric facts, asserted directly against numpy so a future edit to
    # the Rust comparison cannot quietly change them.
    checks += 1
    if not (np.float64(0.3) < np.float32(0.3).astype(np.float64)):
        failures.append("f32(0.3) is expected to widen above the f64 cutoff 0.3")
    checks += 1
    if np.float32(9.9).astype(np.float64) > 9.9:
        failures.append("f32(9.9) is expected to fall below the f64 cutoff 9.9")

    for f in failures:
        print(f"FAIL {f}")
    print(f"\npre_computes cutoff accounting: {checks - len(failures)}/{checks} checks match")
    print(f"  ({len(CASES)} cases x (npeaks, length), plus the ladder and 2 promotion facts)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
