#!/usr/bin/env python3
"""Generate golden parity vectors for `macs-score` from the compiled MACS3.

Drives the real `MACS3.Signal.ScoreTrack.ScoreTrackII` through
`compute_pvalue` + `make_pq_table` and records, per fixture:

* the `-log10` p-score track (as `f32` bit patterns at each breakpoint), and
* the whole p-to-q table as `(p-score bits, q-score bits)` pairs.

The Rust side (`crates/macs-score/tests/parity_vectors.rs`) replays the same
fixtures and requires **bit equality** on both.

Bit equality is the right bar, and it is not achievable by accident: the p-score
lattice is 1e-5 wide (`F1`), the AFDR rank accumulates in `f32` (so it is
quantised above 2^24 base pairs), and `-0.0`/`0.0` share a histogram bucket in
Python's dict but not in a bits-keyed map.

Usage:
    /path/to/macs3-venv/bin/python oracle/gen_score_vectors.py \\
        > crates/macs-score/tests/score_vectors.tsv
"""

import os
import struct
import sys

import numpy as np

from MACS3.Signal.ScoreTrack import ScoreTrackII

sys.stderr.write("generating score vectors from the compiled MACS3\n")

ROWS = []


def f32bits(x):
    return "0x%08x" % struct.unpack(">I", struct.pack(">f", np.float32(x)))[0]


def f64bits(x):
    return "0x%016x" % struct.unpack(">Q", struct.pack(">d", float(x)))[0]


def emit(name, args, payload):
    ROWS.append(f"{name}\t{args}\t{payload}")


def build(cond1, cond2, pseudocount, chrom_len, seed):
    """Run one fixture through the upstream ScoreTrackII.

    `cond1` / `cond2` are lists of (position, treatment, lambda) triples.
    Returns (pvalue_stat, pqtable).
    """
    st = ScoreTrackII(0.0, 0.0, pseudocount)
    # the class stores one (positions, v1, v2, v3) tuple per chromosome; we
    # build the arrays directly, which is what build_chromosome() would do
    n = len(cond1)
    pos = np.zeros(n, dtype="i4")
    v1 = np.zeros(n, dtype="f4")
    v2 = np.zeros(n, dtype="f4")
    for i, (p, t, c) in enumerate(cond1):
        pos[i] = p
        v1[i] = t
        v2[i] = c
    st.data = {"chrT": (pos, v1, v2, np.zeros(n, dtype="f4"))}
    st.datalength = {"chrT": n}
    st.compute_pvalue()
    table = st.make_pq_table()
    return st, table


def collect(st):
    pos, _v1, _v2, v = st.data["chrT"]
    return pos.tolist(), v.tolist()


# ---------------------------------------------------------------- fixtures
# The grid spans the interesting regimes of the AFDR walk:
#   N small (few scored bases) and N large (rank above 2^24 in f32)
#   a single p-score, many distinct p-scores
#   q crossing zero (the walk stops early)
#   treatment/lambda pairs that produce p-scores on both sides of the lattice
PSEUDOCOUNTS = [1.0, 0.5, 2.0]

# deterministic LCG so the corpus is regenerable
class LCG:
    def __init__(self, seed):
        self.s = seed & 0xFFFFFFFFFFFFFFFF

    def next(self, n):
        out = []
        for _ in range(n):
            self.s = (self.s * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
            out.append((self.s >> 33) % n)
        return out


FIXTURES = []

# 1. tiny, hand-checkable
FIXTURES.append(("tiny", 1.0, [(10, 0.0, 0.5), (20, 1.0, 0.5), (30, 2.0, 1.0), (40, 5.0, 0.25)]))

# 2. a single constant score, sweeping N so the walk crosses q = 0
for n_bases in [1, 2, 10, 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000]:
    FIXTURES.append((
        f"const_{n_bases}", 1.0,
        [(0, 3.0, 1.0)] * n_bases,
    ))

# 3. two-level: a strong population and a background population
for total in [1_000, 100_000, 5_000_000]:
    for frac in [0.01, 0.1, 0.5, 0.9]:
        rng = LCG(total * 7919 + int(frac * 1000) + 13)
        picks = rng.next(total)
        rows = []
        for i in range(total):
            depth = 8.0 if (i % 1000) < int(frac * 1000) else 1.0
            rows.append((i, depth, 0.5))
        FIXTURES.append((f"bimodal_{total}_{frac}", 1.0, rows))

# 4. a ladder of distinct p-scores, so the walk has many steps
for n_levels in [2, 5, 20, 100, 500]:
    rng = LCG(n_levels * 104729 + 1)
    rows = []
    for i in range(n_levels):
        treat = float(i)
        lam = 0.5
        rows.append((i * 3, treat, lam))
    FIXTURES.append((f"ladder_{n_levels}", 1.0, rows))

# 5. pseudocount variations on a fixed signal
base = [(i, float(i % 7), 0.5) for i in range(500)]
for pc in PSEUDOCOUNTS:
    FIXTURES.append((f"pc_{pc}", pc, base))

# 6. lambda-heavy: a spread of local lambdas
rows = []
for i in range(1000):
    rows.append((i, float(i % 11), (i % 7) * 0.5))
FIXTURES.append(("lambda_spread", 1.0, rows))

# 7. treatment-heavy: very deep pileups
rows = [(i, float(50 + (i % 500)), 0.25) for i in range(2000)]
FIXTURES.append(("deep", 1.0, rows))


# ---------------------------------------------------------------- run them
sys.stderr.write(f"{len(FIXTURES)} fixtures\n")
for name, pc, rows in FIXTURES:
    try:
        st, table = build(rows, None, pc, 1 << 30, 0)
    except Exception as exc:  # noqa: BLE001
        emit(f"{name}_error", f"{pc!r}", type(exc).__name__)
        sys.stderr.write(f"  {name}: {type(exc).__name__}: {exc}\n")
        continue

    pos, v = collect(st)
    sys.stderr.write(f"  {name}: {len(pos)} breakpoints, {len(table)} table entries\n")

    # the p-score track: positions plus the f32 score at each
    payload = ",".join(str(p) for p in pos) + "\t" + ",".join(f32bits(x) for x in v)
    emit(f"pscore_track\t{name}\t{pc!r}", f"{len(pos)}", payload)

    # the whole p-to-q table
    tpayload = ",".join(f"{f32bits(k)}:{f32bits(val)}" for k, val in sorted(
        table.items(), key=lambda kv: -kv[0]))
    emit(f"pq_table\t{name}\t{pc!r}", f"{len(table)}", tpayload)

sys.stderr.write(f"done: {len(ROWS)} score vectors\n")
sys.stdout.write("\n".join(ROWS) + "\n")
