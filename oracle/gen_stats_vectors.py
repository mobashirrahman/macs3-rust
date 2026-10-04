#!/usr/bin/env python3
"""Generate golden parity vectors for `macs-stats` from the upstream transcription.

Writes `tests/golden/stats_vectors.tsv`, one row per evaluated call, with the
result as its exact IEEE-754 hex bit pattern. The Rust side
(`crates/macs-stats/tests/parity_vectors.rs`) replays the same calls and
requires **bit equality** — not a tolerance. That is how a deliberate upstream
inaccuracy like the log10 Poisson early-exit stays pinned.

Usage:
    python3 oracle/gen_stats_vectors.py > crates/macs-stats/tests/stats_vectors.tsv
"""

import sys
import struct
import math

import prob_reference as P


def bits(x):
    if x != x:  # NaN
        return "nan"
    if x == float("inf"):
        return "inf"
    if x == float("-inf"):
        return "-inf"
    return "0x%016x" % struct.unpack(">Q", struct.pack(">d", float(x)))[0]


def f32bits(x):
    import numpy as np
    return "0x%08x" % struct.unpack(">I", struct.pack(">f", np.float32(x)))[0]


out = []


def emit(fn, args, value, kind="f64"):
    b = f32bits(value) if kind == "f32" else bits(value)
    out.append("\t".join([fn, args, b]))


# ---------------------------------------------------------------- Poisson CDF
# The grid deliberately spans every branch of upstream's dispatch:
#   lambda <= 700  -> linear small-lambda path
#   lambda >  700  -> rescaled large-lambda path
#   k == 0, k small, k near lambda, k >> lambda, k = 0-length edge
LAMS = [
    "1e-8", "0.001", "0.1", "0.5", "1.0", "2.0", "3.7", "5.0", "10.0", "50.0",
    "100.0", "299.0", "300.0", "500.0", "699.0", "699.9", "700.0", "700.1",
    "701.0", "1000.0", "5000.0", "50000.0",
]
KS = [0, 1, 2, 5, 10, 37, 100, 200, 500, 1000]

for lam_s in LAMS:
    lam = float(lam_s)
    for k in KS:
        for lower in (0, 1):
            for log10 in (0, 1):
                try:
                    v = P.poisson_cdf(k, lam, bool(lower), bool(log10))
                except AssertionError:
                    continue
                except Exception:
                    continue
                emit(
                    "poisson_cdf",
                    f"{k}\t{lam_s}\t{lower}\t{log10}",
                    v,
                )

# ------------------------------------------------------- Poisson inverse CDF
for lam in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 300.0, 700.0, 739.0]:
    for cdf in [0.0, 1e-12, 0.001, 0.1, 0.36787944117144233, 0.5, 0.9, 0.99, 0.999, 1.0]:
        for maximum in [100, 1000]:
            v = P.poisson_cdf_inv(cdf, lam, maximum)
            emit("poisson_cdf_inv", f"{lam}\t{cdf!r}\t{maximum}", v, kind="i32")

# -------------------------------------------------------------- Poisson PDF
for lam in [0.5, 1.0, 5.0, 20.0, 100.0]:
    for k in [0, 1, 5, 20, 50, 100, 150]:
        v = math.exp(-1 * lam) * lam ** k
        f = 1.0
        for i in range(2, k + 1):
            f *= i
        v = v / f
        emit("poisson_pdf", f"{k}\t{lam!r}", v)

# ----------------------------------------------------------------- Binomial
for a in [1, 2, 5, 10, 50, 100, 1000, 20000]:
    for b in [1e-8, 1e-6, 0.001, 0.1, 0.5, 0.9, 0.999]:
        for x in sorted({0, 1, a // 4, a // 2, a - 1, a}):
            if x < 0 or x > a:
                continue
            emit("binomial_pdf", f"{x}\t{a}\t{b!r}", P.binomial_pdf(x, a, b))
            emit("binomial_cdf_f", f"{x}\t{a}\t{b!r}", P._binomial_cdf_f(x, a, b))
            emit("binomial_cdf_r", f"{x}\t{a}\t{b!r}", P._binomial_cdf_r(x, a, b))
            emit("binomial_sf_f", f"{x}\t{a}\t{b!r}", P.binomial_sf(x, a, b, True))

# ----------------------------------------------- MACS --keep-dup auto curve
# This is the function callpeak actually makes:
#   max_dup = binomial_cdf_inv(1 - p, N_total, 1 / effective_genome_size)
print("# binomial_cdf_inv is the --keep-dup auto curve", file=sys.stderr)
for gsize in [hs for hs in [293128983, 265278350, 1000000, 100000]]:
    b = 1.0 / gsize
    for total in [1000, 10000, 100000, 1000000]:
        for keep_prob in [0.01, 0.1, 0.5]:
            v = P.binomial_cdf_inv(1.0 - keep_prob, total, b)
            emit("binomial_cdf_inv", f"{1.0 - keep_prob!r}\t{total}\t{b!r}", v, kind="i64")

# ------------------------------------------------------------- pduplication
for n_obs in [10, 100, 1000, 10000, 100000]:
    for p in [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]:
        pmf = [p, p, 2 * p, p]
        v = P.pduplication(pmf, n_obs)
        emit("pduplication", f"{n_obs}\t{p!r}", v, kind="f32")

sys.stdout.write("\n".join(out) + "\n")
