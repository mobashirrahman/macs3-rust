#!/usr/bin/env python3
"""Generate golden parity vectors for `macs-pileup` from the compiled MACS3.

Calls `MACS3.Signal.PileupV2.pileup_from_PN_shifted` and
`pileup_from_LR_as_list` and records the resulting `pv` arrays verbatim: the
breakpoint positions and the `float32` values.

The Rust side (`crates/macs-pileup/tests/parity_vectors.rs`) replays the same
calls and requires identical breakpoint positions and identical `f32` values.
Breakpoints are compared *exactly* including the ones our canonical form
coalesces away, so the comparison also verifies that coalescing is
value-preserving rather than merely plausible.

Usage:
    /path/to/macs3-venv/bin/python oracle/gen_pileup_vectors.py \\
        > crates/macs-pileup/tests/pileup_vectors.tsv
"""

import os
import struct
import sys

import numpy as np

from MACS3.Signal.PileupV2 import pileup_from_PN_shifted, pileup_from_LR_as_list

sys.stderr.write("generating pileup vectors from the compiled MACS3\n")

rows = []


def f32s(values):
    return ["0x%08x" % struct.unpack(">I", struct.pack(">f", np.float32(v)))[0]
            for v in values]


def emit_macs3(fn, args, pv):
    p, v = pv
    rows.append(f"{fn}\t{args}\t"
                f"{','.join(str(int(x)) for x in p)}\t{','.join(f32s(v))}")


def fnv(*parts):
    """A deterministic, order-independent hash of the fixture parameters."""
    h = 0xCBF29CE484222325
    for part in parts:
        for b in str(part).encode():
            h ^= b
            h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return str(h)


# ---------------------------------------------------------------- single-end
# The parameter grid deliberately includes every branch of the clip/sort/sweep
# and the F5 question: pairs of reads whose inferred fragments abut exactly.
SE_CASES = []
for rlength in (100, 1000, 100_000):
    for five_shift, three_shift in (
        (0, 200),        # directional, no shift
        (10, 210),       # directional, shift 10
        (0, 1),          # minimal extension: every read abuts its neighbour
        (100, 100),      # bidirectional, even d
        (100, 101),      # bidirectional, odd d
        (-50, 250),      # negative five_shift: fragment start before the read
        (0, 0),          # zero-length fragments
        (5, 0),          # zero extension at the 3' end
    ):
        for scale in (1.0, 0.5, 0.0, 2.5):
            for baseline in (0.0, 0.5, 1.0):
                SE_CASES.append((rlength, five_shift, three_shift, scale, baseline))


def make_positions(rng, rlength, n):
    return np.sort(np.array([rng.integers(0, rlength) for _ in range(n)], dtype="i4"))


progress = 0
for (rlength, five_shift, three_shift, scale, baseline) in SE_CASES:
    for n_plus, n_minus in ((1, 0), (0, 1), (3, 2), (10, 10), (50, 50), (200, 0)):
        # a deterministic LCG so the corpus is regenerable
        seed = fnv(rlength, five_shift, three_shift, scale, baseline, n_plus, n_minus)
        rng = np.random.default_rng(abs(hash(seed)) % (2**32))
        plus = make_positions(rng, rlength, n_plus)
        minus = make_positions(rng, rlength, n_minus)
        if len(plus) == 0 and len(minus) == 0:
            continue
        try:
            pv = pileup_from_PN_shifted(plus, minus, five_shift, three_shift,
                                        rlength, scale, baseline)
        except Exception as exc:  # noqa: BLE001
            rows.append(f"pn_error\t{seed}\t{type(exc).__name__}\t")
            continue
        args = (f"{rlength}\t{five_shift}\t{three_shift}\t{scale!r}\t{baseline!r}\t"
                f"{','.join(map(str, plus.tolist()))}\t{','.join(map(str, minus.tolist()))}")
        emit_macs3("pn", args, pv)
        progress += 1
        if progress % 500 == 0:
            sys.stderr.write(f"  pn {progress} cases\n")
            sys.stderr.flush()

# --- the F5 question, isolated -------------------------------------------------
# Two reads exactly `d` apart so one inferred fragment ends where the next
# begins, and a third read in between so the depth at the junction is 2.
sys.stderr.write("F5 junction cases\n")
for rlength in (1000, 100_000):
    for d in (10, 50, 200, 1000):
        for (a, b, c) in ((0, d, 2 * d), (5, 5 + d, 5 + 2 * d), (0, d, 3 * d + 1)):
            plus = np.array(sorted({a, b, c} & set(range(rlength))), dtype="i4")
            if len(plus) == 0:
                continue
            pv = pileup_from_PN_shifted(plus, np.array([], dtype="i4"), 0, d, rlength, 1.0, 0.0)
            args = (f"{rlength}\t0\t{d}\t1.0\t0.0\t{','.join(map(str, plus.tolist()))}\t")
            emit_macs3("pn_junction", args, pv)

# ------------------------------------------------------------------ fragments
sys.stderr.write("fragment cases\n")
for rlength in (100, 1000, 100_000):
    for n in (1, 5, 100, 1000):
        for scale in (1.0, 0.25):
            seed = fnv("lr", rlength, n, scale)
            rng = np.random.default_rng(abs(hash(seed)) % (2**32))
            l = np.sort(rng.integers(0, rlength, size=n)).astype("i4")
            width = rng.integers(1, max(2, rlength // 4), size=n).astype("i4")
            r = np.minimum(l + width, rlength).astype("i4")
            lr = np.zeros(n, dtype=[("l", "i4"), ("r", "i4")])
            lr["l"] = l
            lr["r"] = r
            keep = r > l
            lr = lr[keep]
            if len(lr) == 0:
                continue
            pv = pileup_from_LR_as_list(lr, scale, 0.0)
            args = (f"{rlength}\t{scale!r}\t"
                    f"{','.join(map(str, lr['l'].tolist()))}\t"
                    f"{','.join(map(str, lr['r'].tolist()))}")
            emit_macs3("lr", args, pv)

sys.stderr.write(f"done: {len(rows)} pileup vectors\n")
sys.stdout.write("\n".join(rows) + "\n")
