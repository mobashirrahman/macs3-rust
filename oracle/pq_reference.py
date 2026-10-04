#!/usr/bin/env python3
"""Faithful transcription of MACS3 3.0.5's p-score and p-to-q pipeline.

Transcribed from `MACS3/Signal/ScoreTrack.py:430-544`
(`compute_pvalue`, `make_pq_table`) with `MACS3.Signal.Prob.poisson_cdf` — the
**compiled** one — supplying the p-scores, so the only thing being transcribed is
the histogram construction and the AFDR walk.

`ScoreTrackII` is a Cython `cdef class` whose attributes and whose `get_pscore`
are not reachable from Python, so the table cannot be driven directly. This
transcription is the bridge, and it is validated two ways:

1. `crates/macs-score/tests/parity_vectors.rs` checks the Rust `PqTable` against
   this transcription bit-for-bit over an exhaustive set of histograms.
2. `tests/golden/**/**_cutoff_analysis.txt` is a **direct readout of the real
   compiled MACS3's p-to-q table** (it prints `self.pqtable[cutoff]` for 33
   cutoffs). Those files are captured by `oracle/gen_pq_vectors.py` and are the
   end-to-end check once the callpeak pipeline exists (gate G9).

## Numeric details that decide whether this matches

* `k: cython.float` — the AFDR rank accumulates in **float32**, so above 2^24
  its spacing exceeds 1 and the rank is quantised.
* `f: cython.float` — `-log10(N)` is computed in f64 and then **truncated to
  float32** on assignment.
* `q: cython.float` — `v + (log10(k) + f)` is evaluated in f64 (`v` and `f` are
  f32 widened, `log10` is a Python function) and then **truncated to float32**.
* `pre_q: cython.float`, initialised to `2147483647`, which f32 cannot represent
  and which rounds to 2^31.
* The histogram is keyed by Python `float`, so `-0.0` and `0.0` **share a bucket**
  (`hash(-0.0) == hash(0.0)` and `-0.0 == 0.0`).
* The walk breaks at the first `q <= 0`; that p-score and every lower one are then
  written as `0` by the following loop.
"""

import math
import struct

import numpy as np
from MACS3.Signal.Prob import poisson_cdf

__all__ = [
    "pscore",
    "pvalue_stat",
    "make_pq_table",
    "run_pipeline",
    "f32",
]


def f32(x):
    """Round a Python float to float32, as a C `cython.float` assignment does."""
    return float(np.float32(x))


def f32bits(x):
    return "0x%08x" % struct.unpack(">I", struct.pack(">f", np.float32(x)))[0]


def pscore(observed, expectation, pseudocount=1.0):
    """`get_pscore(int(treat + pseudocount), lambda + pseudocount)`.

    `expectation` is a C `float`, so the p-score is cached and returned at
    float32 precision.
    """
    o = int(observed + pseudocount)
    e = f32(expectation + pseudocount)
    return f32(-1.0 * poisson_cdf(o, e, False, True))


def pvalue_stat(treat, control, pseudocount=1.0):
    """Build upstream's `pvalue_stat`: p-score -> number of base pairs.

    `treat` and `control` are lists of (position, value) with the **right-endpoint**
    convention: `value` applies to `[prev_position, position)`, and
    `prev_position` starts at 0. See `F5` in `docs/upstream-findings.md`.

    Chromosomes are summed in ascending name order, matching
    `for chrom in sorted(self.data.keys())`.
    """
    stat = {}
    for _chrom in sorted(treat.keys()):
        prev = 0
        trows = sorted(treat[_chrom], key=lambda r: r[0])
        crows = dict(control[_chrom])
        for pos, tval in trows:
            v = pscore(tval, crows[pos], pseudocount)
            length = pos - prev
            try:
                stat[v] += length
            except KeyError:
                stat[v] = length
            prev = pos
    return stat


def make_pq_table(pvaluestat):
    """`ScoreTrackII.make_pq_table`."""
    if not pvaluestat:
        return {}

    n = 0
    for v in pvaluestat.values():
        n += v

    k = f32(1.0)
    f = f32(-1.0 * math.log10(n))
    pre_q = f32(2147483647.0)

    table = {}
    unique_values = sorted(pvaluestat.keys(), reverse=True)
    i = 0
    for i, v in enumerate(unique_values):
        ln = pvaluestat[v]
        # f64 arithmetic, truncated to float32 on assignment
        q = f32(v + (math.log10(k) + f))
        if q > pre_q:
            q = pre_q
        if q <= 0:
            break
        table[v] = q
        pre_q = q
        k = f32(k + ln)

    for j in range(i, len(unique_values)):
        table[unique_values[j]] = 0.0
    return table


def run_pipeline(treat, control, pseudocount=1.0):
    """`pvalue_stat` + `make_pq_table`, i.e. the whole p/q preparation."""
    stat = pvalue_stat(treat, control, pseudocount)
    return stat, make_pq_table(stat)
