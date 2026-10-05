#!/usr/bin/env python3
"""Emit the upstream `get_pscore` reference grid for the Rust differential.

`get_pscore(k, lam) = -1 * poisson_cdf(k, lam, lower=False, log10=True)` -- the
upper tail of a Poisson, in log10 space. Writes a TSV of (k, lam, value) that
`crates/macs-score/tests/pscore_vs_upstream.rs` asserts against.

Usage: oracle/gen_pscore_grid.py > crates/macs-score/tests/pscore_grid.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oracle_env import ensure_oracle_python, require_src  # noqa: E402

ensure_oracle_python()
sys.path.insert(0, require_src())

from MACS3.Signal.Prob import poisson_cdf  # noqa: E402

print("# k\tlambda\tupstream_get_pscore")
for k in [1, 2, 3, 5, 8, 9, 10, 12, 15, 20, 25, 30, 40]:
    for lam in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0]:
        print("%d\t%.1f\t%.9g" % (k, lam, -1.0 * poisson_cdf(k, lam, False, True)))
