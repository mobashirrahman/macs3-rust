#!/usr/bin/env python
"""Differential for `macs_stats::NumpyRng` against the installed NumPy.

MACS3's downsamplers seed NumPy's global RNG and shuffle position arrays in
place, so which reads survive is a function of NumPy's exact stream. A shuffle
with the right distribution is not enough: the retained *set* has to match.

This compares the raw `next_uint32` stream (so a divergence is localised to the
generator rather than to `shuffle`), two full permutations, and the first ten
elements of a 100-element permutation.

    python3 oracle/check_rng.py [--cargo] [--seeds 0,1,42,12345]
"""

import argparse
import subprocess
import sys
import tempfile
import os

# NumPy lives in the oracle virtualenv, so this re-execs under the provisioned
# interpreter rather than whichever `python3` happens to be first on PATH.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oracle_env import ensure_oracle_python  # noqa: E402

ensure_oracle_python()

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def rust_outputs(extra_seeds, cargo):
    """Run `cargo run --example rngcheck` and parse its `key<TAB>values` lines."""
    cmd = ["cargo", "run", "-q", "-p", "macs-stats", "--example", "rngcheck"]
    if cargo:
        cmd = ["cargo", "run", "-q", "--release", "-p", "macs-stats",
               "--example", "rngcheck"]
    p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"cargo run failed:\n{p.stderr[-2000:]}")
    out = {}
    for line in p.stdout.splitlines():
        if "\t" not in line:
            continue
        k, v = line.split("\t", 1)
        out[k] = [int(x) for x in v.split(",")]
    if not out:
        sys.exit("no output from rngcheck; the example changed its format")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cargo", action="store_true",
                    help="build with --release, matching CI benchmarking")
    ap.add_argument("--seeds", default="0,1,42,12345")
    args = ap.parse_args()
    seeds = [int(s) for s in args.seeds.split(",")]

    rust = rust_outputs(seeds, args.cargo)

    checks = 0
    failures = []

    for seed in seeds:
        rs = np.random.RandomState(seed)
        words = [int(rs.randint(0, 1 << 32, dtype=np.uint32)) for _ in range(20)]
        key = f"words_{seed}"
        checks += 1
        if key not in rust:
            failures.append(f"{key}: missing from the Rust side")
        elif rust[key] != words:
            i = next(j for j, (a, b) in enumerate(zip(words, rust[key])) if a != b)
            failures.append(
                f"{key}: first difference at index {i}: "
                f"numpy={words[i]} rust={rust[key][i]}")

    for seed in seeds[:4]:
        rs = np.random.RandomState(seed)
        a = list(range(50))
        rs.shuffle(a)
        key = f"shuffle_{seed}"
        checks += 1
        if key not in rust:
            failures.append(f"{key}: missing from the Rust side")
        elif rust[key] != a:
            i = next(j for j, (x, y) in enumerate(zip(a, rust[key])) if x != y)
            failures.append(f"{key}: first difference at index {i}: "
                            f"numpy={a[i]} rust={rust[key][i]}")

    for seed in seeds[:3]:
        rs = np.random.RandomState(seed)
        a = list(range(100))
        rs.shuffle(a)
        key = f"shuffle100_{seed}"
        checks += 1
        if key not in rust:
            failures.append(f"{key}: missing from the Rust side")
        elif rust[key] != a[:10]:
            failures.append(f"{key}: numpy[:10]={a[:10]} rust={rust[key]}")

    for f in failures:
        print(f"FAIL {f}")
    print(f"\n{checks - len(failures)}/{checks} RNG checks match NumPy "
          f"{np.__version__}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
