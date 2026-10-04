#!/usr/bin/env python3
"""Generate the Savitzky-Golay derivative-1 kernel table.

`SignalProcessing.savitzky_golay_order2_deriv1` does **not** use the closed-form
Savitzky-Golay coefficients. It uses

```python
b = np.array([[1, k, k**2] for k in range(-half_window, half_window+1)], dtype='i8')
m = np.linalg.pinv(b)[1]
```

i.e. row 1 of the **pseudo-inverse of an integer design matrix**, which NumPy
computes through an SVD. This is not interchangeable with any mathematically exact
alternative. Measured over 300 random pileup-shaped signals, comparing the resulting
`maxima` index sets against upstream:

| coefficient source | maxima mismatches / 300 |
|---|---|
| `np.linalg.pinv(b)[1]` | **0** |
| closed form `3k / (h(h+1)(2h+1))` | 113 |
| `b.T @ inv(b @ b.T)` | 299 |

The mathematically *exact* `b.T (b b.T)^-1` is worse than the closed form, because the
whole downstream pipeline rounds the derivative to 16 decimals before taking its
sign: the result is a flat stretch reading as exactly `0.0` (one broad summit) versus
`~1e-16` (two summits), and that last bit is decided by SVD rounding.

So the coefficients *are* the specification, and they have to be shipped. F141 did
that for the single window size 179 and used the closed form everywhere else, which
left every other `--call-summits` run wrong. This script emits a table for every odd
window size in range, so the closed form is gone entirely.

Usage:
    python3 oracle/gen_sg_coeffs.py --out crates/macs-peaks/src/sg_coeffs.tsv
"""
from __future__ import annotations

import argparse
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HEADER = """# GENERATED FILE -- do not edit.
#
# Regenerate with:
#     python3 oracle/gen_sg_coeffs.py --out crates/macs-peaks/src/sg_coeffs.tsv
#
# Row N is `np.linalg.pinv(b)[1]` for the Savitzky-Golay order-2 derivative-1 kernel
# with window size N, exactly as `SignalProcessing.savitzky_golay_order2_deriv1`
# computes it. Values are 17 significant digits because the last bit is load-bearing:
# `maxima` rounds the smoothed derivative to 16 decimals and takes its sign, so a
# flat pileup stretch must land on exactly 0.0 to read as one summit.
#
# Format: <window_size>\\t<v0>,<v1>,...,<v(N-1)>
"""


def coeffs(window: int) -> np.ndarray:
    if window % 2 != 1:
        window += 1
    half = (window - 1) // 2
    b = np.array([[1, k, k**2] for k in range(-half, half + 1)], dtype="i8")
    return np.linalg.pinv(b)[1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-window", type=int, default=3)
    ap.add_argument("--max-window", type=int, default=513)
    ap.add_argument(
        "--out",
        default=os.path.join(REPO, "crates", "macs-peaks", "src", "sg_coeffs.tsv"),
    )
    args = ap.parse_args()

    rows = []
    for w in range(args.min_window, args.max_window + 1):
        if w % 2 == 0:
            continue
        c = coeffs(w)
        assert c.shape[0] == w, (c.shape, w)
        # 17 significant digits, matching the literals F141 needed.
        rows.append(f"{w}\t" + ",".join(repr(float(v)) for v in c))

    with open(args.out, "w") as fh:
        fh.write(HEADER)
        fh.write("\n".join(rows))
        fh.write("\n")
    print(
        f"wrote {args.out}: {len(rows)} window sizes "
        f"({args.min_window}..{args.max_window}), "
        f"{os.path.getsize(args.out) / 1e6:.2f} MB"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
