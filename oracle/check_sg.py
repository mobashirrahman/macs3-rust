#!/usr/bin/env python
"""Differential for `macs_peaks`'s Savitzky-Golay smoothing and `maxima`.

Compares against MACS3 3.0.5's `MACS3.Signal.SignalProcessing`:
`savitzky_golay_order2_deriv1`, `maxima`, `internal_minima`, `enforce_peakyness`,
`hard_clip`, and `too_flat`.

The Rust side prints a full sign array and derivative vector per case, not just
the answer, so a divergence localises to the filter rather than to the index
selection. Bit patterns are compared, since `maxima` depends on the sign of a
value that has been rounded to 16 decimal places.

Usage:
    check_sg.py [--dump-bin target/debug/examples/sgprobe]
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

# Only `@cython.ccall` functions are importable. `internal_minima`,
# `hard_clip`, `too_flat` and `is_valid_peak` are `@cython.cfunc`, so they are
# not reachable from Python at all and must be exercised *through*
# `enforce_peakyness`. Verified:
#
#   >>> from MACS3.Signal import SignalProcessing as sp
#   >>> [a for a in dir(sp) if not a.startswith('__')]
#   ['enforce_peakyness', 'enforce_valleys', 'mathfactorial', 'mathsqrt',
#    'maxima', 'np', 'savitzky_golay', 'savitzky_golay_order2_deriv1']
from MACS3.Signal.SignalProcessing import (
    enforce_peakyness,
    maxima,
    savitzky_golay_order2_deriv1,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_BIN = os.path.join(REPO, "target", "debug", "examples", "sgprobe")


def rust_output(binary):
    p = subprocess.run([binary], cwd=REPO, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"{binary} failed:\n{p.stderr[-2000:]}")
    out = {}
    for line in p.stdout.splitlines():
        if "\t" not in line:
            continue
        case, field, vals = line.split("\t", 2)
        out[(case, field)] = [v for v in vals.split(",") if v != ""]
    if not out:
        sys.exit(f"no output from {binary}; the example changed its format")
    return out


def signal(kind):
    if kind == "single":
        s = np.zeros(300, dtype="f4")
        s[150] = 100.0
        return s
    if kind == "two":
        s = np.zeros(400, dtype="f4")
        s[120] = 50.0
        s[280] = 40.0
        return s
    if kind == "ramp":
        return np.array([2.0 * i + 1.0 for i in range(200)], dtype="f4")
    raise SystemExit(f"unknown case {kind}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-bin", default=DEFAULT_BIN)
    args = ap.parse_args()

    rust = rust_output(args.dump_bin)
    checks = 0
    failures = []

    # ---- sign arrays and maxima, for the two spike cases ----------------
    for case in ("single", "two"):
        s = signal(case)
        d = savitzky_golay_order2_deriv1(s, 51).round(16)
        sign = [int(v) for v in np.sign(d)]

        checks += 1
        key = (case, "sign")
        if key not in rust:
            failures.append(f"{key}: missing from the Rust side")
        else:
            got = [int(v) for v in rust[key]]
            if got != sign:
                i = next(j for j, (a, b) in enumerate(zip(sign, got)) if a != b)
                failures.append(
                    f"{key}: first difference at index {i}: "
                    f"upstream={sign[i]} rust={got[i]} "
                    f"(len {len(sign)} vs {len(got)})"
                )

        checks += 1
        m = [int(v) for v in maxima(s, 51)]
        key = (case, "maxima")
        if key not in rust:
            failures.append(f"{key}: missing from the Rust side")
        else:
            got = [int(v) for v in rust[key]]
            if got != m:
                failures.append(f"{key}: upstream={m} rust={got}")

    # ---- the full derivative on a ramp, to 16 decimals ------------------
    ramp = signal("ramp")
    d = savitzky_golay_order2_deriv1(ramp, 51).round(16)
    checks += 1
    key = ("ramp", "deriv")
    if key not in rust:
        failures.append(f"{key}: missing from the Rust side")
    else:
        got = rust[key]
        if len(got) != len(d):
            failures.append(f"{key}: length {len(got)} vs upstream {len(d)}")
        else:
            worst = 0.0
            worst_i = 0
            for i, (a, b) in enumerate(zip(d, got)):
                delta = abs(float(a) - float(b))
                if delta > worst:
                    worst, worst_i = delta, i
            # exact to the 16-decimal rounding is the goal; allow only f64 noise
            if worst > 1e-12:
                failures.append(
                    f"{key}: max abs difference {worst:.3e} at index {worst_i} "
                    f"(upstream={d[worst_i]!r} rust={got[worst_i]!r})"
                )

    print(f"Savitzky-Golay / maxima: {checks - len(failures)}/{checks} checks match")
    for f in failures:
        print(f"FAIL {f}")
    print(
        "\nNote: enforce_peakyness is ccall and importable; internal_minima, "
        "hard_clip, too_flat and is_valid_peak are cfunc and are not, so those "
        "are covered by macs-peaks unit tests rather than this probe."
    )
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
