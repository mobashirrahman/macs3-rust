#!/usr/bin/env python3
"""Generate golden vectors for `macs-score` (the p-score and the p-to-q table).

Two sources, written to one TSV:

**`table` rows** — an explicit `(p-score, base-pair count)` histogram, and the
p-to-q table `MACS3.Signal.ScoreTrack.make_pq_table` produces from it, computed
by `oracle/pq_reference.py` (which takes p-scores from the *compiled*
`poisson_cdf`, and transcribes only the histogram construction and the AFDR
walk). `ScoreTrackII` is a Cython `cdef class` so its table cannot be driven
directly; this is the bridge, and the Rust test requires **bit equality** with it.

**`cutoff` rows** — a *direct* readout of the real compiled MACS3's p-to-q table,
from `--cutoff-analysis`, which prints `self.pqtable[cutoff]` for 33 cutoffs.
These are captured now as a committed artifact and used as the end-to-end check
once the callpeak pipeline exists (gate G9); they are not a substitute for the
`table` rows.

Usage:
    /path/to/macs3-venv/bin/python oracle/gen_pq_vectors.py \\
        > crates/macs-score/tests/pq_vectors.tsv
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pq_reference as P  # noqa: E402

ROWS = []


def emit_table(name, hist):
    # upstream's `pvalue_stat` is keyed by C `float` values (the p-score track is
    # `f32`), so the keys are f32 numbers widened to f64. Feeding raw f64
    # literals would compare a different value than upstream ever sees, and the
    # divergence shows up as a 1-ULP difference in q.
    stat = {P.f32(k): v for k, v in hist}
    table = P.make_pq_table(stat)
    total = sum(stat.values())
    # histogram as p-score:count, sorted by p-score bits for a stable encoding
    h = ",".join(f"{P.f32bits(k)}:{v}" for k, v in sorted(stat.items()))
    t = ",".join(f"{P.f32bits(k)}:{P.f32bits(v)}" for k, v in sorted(table.items()))
    ROWS.append(f"table\t{name}\t{total}\t{len(stat)}\t{len(table)}\t{h}\t{t}")


# --------------------------------------------------------------- histograms
# The grid targets every regime the AFDR walk has:
HISTOGRAMS = []


def add(name, pairs):
    HISTOGRAMS.append((name, pairs))


# 1. a single bucket: q = v - log10(N)
for total in (1, 2, 10, 1000, 1_000_000):
    add(f"single_{total}", [(5.0, total)])

# 2. two levels, the classic strong/weak split
add("two_900_100", [(5.0, 900), (2.0, 100)])
add("two_1_999", [(10.0, 1), (1.0, 999)])
add("two_5000_1", [(20.0, 5000), (1.0, 1)])

# 3. ladders of distinct p-scores, so the walk has many steps and the
#    monotone cap actually binds
for levels in (2, 3, 5, 10, 20, 50):
    add(f"ladder_{levels}", [(0.5 * i, 10 * i + 1) for i in range(1, levels + 1)])

# 4. N large enough that the float32 rank is quantised above 2^24
add("f32_rank_2p24", [(4.0, (1 << 24) + 7), (3.0, 3), (2.0, 1)])
add("f32_rank_2p30", [(6.0, (1 << 30) - 5), (1.0, 1)])
add("f32_rank_many", [(12.0, (1 << 25)), (8.0, 999), (4.0, 17)])

# 5. the walk must stop: f = -log10(N) is very negative, so q <= 0 almost at once
add("cut_early", [(1.0, 1), (0.5, 1), (0.00001, 2_000_000_000)])
add("cut_at_first", [(0.5, 4_000_000_000)])

# 6. -0.0 and 0.0 must share a bucket, as they do in a Python dict
add("signed_zero", [(-0.0, 500), (0.0, 500), (1.0, 10)])
add("signed_zero_only", [(-0.0, 1_000_000)])

# 7. a realistic peak-caller shape: most bases at the background score, a tail
# of stronger scores, all on the 1e-5 lattice
rng = np.random.default_rng(20240930)
realistic = {}
for i in range(4000):
    p = float(rng.choice([0.00001, 0.00002, 0.1, 0.5, 1.0, 2.0, 5.0, 12.34567]))
    realistic[p] = realistic.get(p, 0) + int(rng.integers(1, 5000))
add("realistic", list(realistic.items()))

# 8. a histogram where the monotone cap binds repeatedly (descending p with
#    growing counts can push q above the previous value)
add("cap_binds", [(9.0, 1), (8.9, 10), (8.8, 100), (8.7, 1000), (8.6, 10_000)])

sys.stderr.write(f"{len(HISTOGRAMS)} histograms\n")
for name, pairs in HISTOGRAMS:
    emit_table(name, pairs)


# ---------------------------------------------- real MACS3 cutoff analysis
# `--cutoff-analysis` prints the real table: `pscore qscore npeaks lpeaks ...`
FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "fixtures")
macs3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                     "venv", "bin", "macs3")
if not os.path.exists(macs3):
    macs3 = os.environ.get("MACS3_BIN", "macs3")


def gsize_of(fdir):
    total = 0
    with open(os.path.join(fdir, "genome.txt")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            total += int(line.split()[1])
    return total


def cutoff_pairs(fdir):
    manifest = {}
    for line in open(os.path.join(fdir, "manifest.tsv")):
        parts = line.rstrip("\n").split("\t")
        manifest[parts[0]] = parts[1] if len(parts) > 1 else ""
    mode = manifest.get("mode", "se")
    fmt = {"se": "BED", "pe": "BEDPE", "frag": "FRAG"}[mode]
    treat = os.path.join(fdir, f"treat.{'bed' if mode == 'se' else ('bedpe' if mode == 'pe' else 'frag')}")
    ctrl = os.path.join(fdir, f"ctrl.{'bed' if mode == 'se' else ('bedpe' if mode == 'pe' else 'frag')}")
    base = ["-t", treat]
    if os.path.exists(ctrl):
        base += ["-c", ctrl]
    base += ["-f", fmt]
    if manifest.get("model_capable", "1") != "1":
        base += ["--nomodel", "--extsize", "200"]

    with tempfile.TemporaryDirectory() as td:
        cmd = [macs3, "callpeak", "-n", "pq", "-g", str(gsize_of(fdir)),
               "--outdir", td, "--cutoff-analysis"] + base
        try:
            subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"  {fdir}: {exc}\n")
            return None
        path = os.path.join(td, "pq_cutoff_analysis.txt")
        if not os.path.exists(path):
            return None
        pairs = []
        for line in open(path):
            f = line.split()
            if len(f) < 2 or f[0] == "pscore":
                continue
            pairs.append((float(f[0]), float(f[1])))
        return pairs


if os.path.isdir(FIXTURES):
    found = 0
    for group in sorted(os.listdir(FIXTURES)):
        gdir = os.path.join(FIXTURES, group)
        if not os.path.isdir(gdir):
            continue
        for name in sorted(os.listdir(gdir)):
            fdir = os.path.join(gdir, name)
            if not os.path.isfile(os.path.join(fdir, "manifest.tsv")):
                continue
            pairs = cutoff_pairs(fdir)
            if not pairs:
                continue
            found += 1
            payload = ",".join(f"{P.f32bits(p)}:{P.f32bits(q)}" for p, q in pairs)
            ROWS.append(f"cutoff\t{group}/{name}\t{len(pairs)}\t0\t0\t{payload}\t")
            sys.stderr.write(f"  cutoff {group}/{name}: {len(pairs)} pairs\n")
            sys.stderr.flush()
    sys.stderr.write(f"{found} fixtures contributed real cutoff-analysis tables\n")

sys.stderr.write(f"done: {len(ROWS)} vectors\n")
sys.stdout.write("\n".join(ROWS) + "\n")
