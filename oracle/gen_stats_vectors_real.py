#!/usr/bin/env python3
"""Generate golden parity vectors for `macs-stats` from the *compiled* MACS3.

This is the authoritative oracle: it calls the real Cython extensions from the
pinned upstream build (`oracle/ENV.lock`), not a transcription. The transcription
in `prob_reference.py` exists only so the port could start before the build
finished; anything it and this script disagree on is a transcription bug and must
be resolved in favour of this script.

Writes a TSV of `function <TAB> args <TAB> result-as-ieee754-bits`.

Usage:
    /path/to/macs3-venv/bin/python oracle/gen_stats_vectors_real.py \\
        > crates/macs-stats/tests/stats_vectors.tsv
"""

import os
import struct
import sys

import numpy as np

from MACS3.Signal.Prob import (
    poisson_cdf,
    poisson_cdf_inv,
    poisson_cdf_Q_inv,
    poisson_pdf,
    binomial_pdf,
    binomial_cdf,
    binomial_sf,
    binomial_cdf_inv,
    pduplication,
    pnorm,
    pnorm2,
    chisq_pvalue_e,
    chisq_logp_e,
    factorial,
)
import MACS3
import MACS3.Signal.Prob as _prob

# `MACS3.__version__` is not exported by the package; the version lives in
# distribution metadata, and the CLI prints it. Prefer metadata, fall back to the
# `macs3 --version` output.
def _detect_version():
    try:
        from importlib.metadata import version
        return version("MACS3")
    except Exception:
        pass
    try:
        import subprocess
        out = subprocess.run(
            [os.path.join(os.path.dirname(sys.executable), "macs3"), "--version"],
            capture_output=True, text=True, timeout=60,
        )
        return out.stdout.strip() or out.stderr.strip() or "unknown"
    except Exception:
        return "unknown"


MACS3_VERSION = _detect_version()
if MACS3_VERSION == "unknown":
    sys.stderr.write("WARNING: could not detect the MACS3 version\n")

# `poz`, `logspace_add`, `binomial_coef` and `ex20` are Cython `cdef` functions:
# they exist in the compiled module but are not exported to Python, so they
# cannot be vectored directly. They are covered by the transcription in
# prob_reference.py plus the unit tests in crates/macs-stats.
for _n in ("poz", "logspace_add", "binomial_coef", "ex20", "_binomial_cdf_f"):
    assert not hasattr(_prob, _n), f"upstream now exports {_n}; update this script"

sys.stderr.write(f"oracle: MACS3 {MACS3_VERSION} from {MACS3.__file__}\n")


def f64bits(x):
    x = float(x)
    if x != x:
        return "nan"
    if x == float("inf"):
        return "inf"
    if x == float("-inf"):
        return "-inf"
    return "0x%016x" % struct.unpack(">Q", struct.pack(">d", x))[0]


def f32bits(x):
    return "0x%08x" % struct.unpack(">I", struct.pack(">f", np.float32(x)))[0]


def i64bits(x):
    return str(int(x))


rows = []


def emit(fn, args, bits):
    rows.append(f"{fn}\t{args}\t{bits}")


def progress(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


# ---------------------------------------------------------------- Poisson CDF
# The grid spans every branch of upstream's dispatch:
#   lambda <= 700 -> linear small-lambda path
#   lambda >  700 -> rescaled large-lambda path
# plus k below, near and far above lambda.
LAMS = [
    "1e-8", "0.001", "0.1", "0.5", "1.0", "2.0", "3.7", "5.0", "10.0", "50.0",
    "100.0", "299.0", "300.0", "500.0", "699.0", "699.9", "700.0", "700.1",
    "701.0", "1000.0", "5000.0", "50000.0",
]
KS = [0, 1, 2, 5, 10, 37, 100, 200, 500, 1000]

progress("[1/7] poisson_cdf")
for lam_s in LAMS:
    lam = float(lam_s)
    for k in KS:
        for lower in (0, 1):
            for log10 in (0, 1):
                try:
                    v = poisson_cdf(k, lam, bool(lower), bool(log10))
                except Exception as exc:  # noqa: BLE001
                    emit("poisson_cdf_error", f"{k}\t{lam_s}\t{lower}\t{log10}",
                         type(exc).__name__)
                    continue
                emit("poisson_cdf", f"{k}\t{lam_s}\t{lower}\t{log10}", f64bits(v))

progress("[2/7] poisson_cdf_inv")
for lam in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 300.0, 700.0, 739.0]:
    for cdf in [0.0, 1e-12, 0.001, 0.1, 0.36787944117144233, 0.5, 0.9, 0.99, 0.999, 1.0]:
        for maximum in [100, 1000]:
            v = poisson_cdf_inv(cdf, lam, maximum)
            emit("poisson_cdf_inv", f"{lam!r}\t{cdf!r}\t{maximum}", i64bits(v))
            # upstream's `poisson_cdf_Q_inv` is a byte-for-byte copy of
            # `poisson_cdf_inv` and does not invert the upper tail at all
            q = poisson_cdf_Q_inv(cdf, lam, maximum)
            emit("poisson_cdf_Q_inv", f"{lam!r}\t{cdf!r}\t{maximum}", i64bits(q))

progress("[3/7] poisson_pdf")
for lam in [0.5, 1.0, 5.0, 20.0, 100.0]:
    for k in [0, 1, 5, 20, 50, 100, 150]:
        try:
            v = poisson_pdf(k, lam)
        except Exception as exc:  # noqa: BLE001
            emit("poisson_pdf_error", f"{k}\t{lam!r}", type(exc).__name__)
            continue
        emit("poisson_pdf", f"{k}\t{lam!r}", f64bits(v))

progress("[4/7] binomial pdf/cdf")
# F11: `binomial_pdf` does not terminate for tiny `b` with `x > 0` and
# `x <= a - x` (the `pdf /= (1 - p)` recovery loop advances by ~1e-8 per step).
# Those cells are skipped here on purpose; the Rust port bounds the loop and
# `docs/upstream-findings.md` F11 records the deviation.
PDF_SAFE_B = 1e-8
for a in [1, 2, 5, 10, 50, 100, 1000, 20000]:
    for b in [1e-8, 1e-6, 0.001, 0.1, 0.5, 0.9, 0.999]:
        for x in sorted({0, 1, a // 4, a // 2, a - 1, a}):
            if x < 0 or x > a:
                continue
            # F11: the `pdf /= (1 - p)` recovery loop inside `binomial_pdf` does
            # not terminate for tiny `b` when the running product drops below
            # 1e-100. It is reachable from the *cdf* functions too, because they
            # call `binomial_pdf(argmax, a, b)` internally with `argmax` chosen
            # by them. So for `b < 1e-3` only `x == 0` is vectored; the cell is
            # recorded as skipped rather than silently dropped.
            safe = not (b < 1e-3 and x > 0)
            args = f"{x}\t{a}\t{b!r}"
            progress(f"      a={a} b={b!r} x={x}")
            if safe:
                emit("binomial_pdf", args, f64bits(binomial_pdf(x, a, b)))
                emit("binomial_cdf_f", args, f64bits(binomial_cdf(x, a, b, True)))
                emit("binomial_cdf_r", args, f64bits(binomial_cdf(x, a, b, False)))
                emit("binomial_sf_f", args, f64bits(binomial_sf(x, a, b, True)))
                emit("binomial_sf_r", args, f64bits(binomial_sf(x, a, b, False)))
            else:
                emit("binomial_skipped", args, "F11")

progress("[5/7] keep-dup auto curve (binomial_cdf_inv)")
# This is the call callpeak actually makes for `--keep-dup auto`:
#   max_dup = binomial_cdf_inv(1 - p, N_total, 1 / effective_genome_size)
for gsize in [293128983, 265278350, 1000000, 100000]:
    b = 1.0 / gsize
    for total in [1000, 10000, 100000, 1000000]:
        for keep_prob in [0.01, 0.1, 0.5]:
            v = binomial_cdf_inv(1.0 - keep_prob, total, b)
            emit("binomial_cdf_inv", f"{1.0 - keep_prob!r}\t{total}\t{b!r}", i64bits(v))
    progress(f"      gsize={gsize} done")

progress("[6/7] pduplication (f32 accumulation)")
for n_obs in [10, 100, 1000, 10000, 100000]:
    for p in [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]:
        pmf = np.array([p, p, 2 * p, p], dtype="f8")
        v = pduplication(pmf, n_obs)
        emit("pduplication", f"{n_obs}\t{p!r}", f32bits(v))

progress("[7/7] normal + chisq + factorial")
for x in range(-30, 31):
    for u in (-10, 0, 10):
        for v in (1, 4, 25):
            emit("pnorm", f"{x}\t{u}\t{v}", f64bits(pnorm(x, u, v)))
for x in np.arange(-6.0, 6.0, 0.125, dtype="f4"):
    for u in (-2.0, 0.0, 1.5):
        for v in (0.5, 1.0, 4.0):
            emit("pnorm2", f"{float(x)!r}\t{u!r}\t{v!r}", f32bits(pnorm2(x, u, v)))
for n in [0, 1, 5, 10, 20, 25, 100, 170]:
    emit("factorial", f"{n}", f64bits(factorial(n)))
for df in [2, 4, 6, 10, 20, 50]:
    for x in [0.5, 1.0, 3.0, 10.0, 25.0, 39.0, 41.0, 80.0, 200.0]:
        try:
            emit("chisq_pvalue_e", f"{x!r}\t{df}", f64bits(chisq_pvalue_e(x, df)))
        except Exception as exc:  # noqa: BLE001
            emit("chisq_pvalue_e_error", f"{x!r}\t{df}", type(exc).__name__)
        for log10 in (False, True):
            try:
                v = chisq_logp_e(x, df, log10)
            except Exception as exc:  # noqa: BLE001
                emit("chisq_logp_e_error", f"{x!r}\t{df}\t{int(log10)}", type(exc).__name__)
                continue
            emit("chisq_logp_e", f"{x!r}\t{df}\t{int(log10)}", f64bits(v))

progress(f"done: {len(rows)} vectors")
sys.stdout.write("\n".join(rows) + "\n")

# record the environment so CI can prove the oracle was the expected one
meta = os.environ.get("MACS3_ORACLE_META")
if meta:
    with open(meta, "w") as fh:
        fh.write(f"MACS3_VERSION={MACS3_VERSION}\n")
        fh.write(f"MACS3_PATH={MACS3.__file__}\n")
        fh.write(f"PYTHON={sys.version.split()[0]}\n")
        fh.write(f"NUMPY={np.__version__}\n")
        fh.write(f"VECTORS={len(rows)}\n")
