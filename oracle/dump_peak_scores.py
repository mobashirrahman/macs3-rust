#!/usr/bin/env python3
"""Capture upstream's peak tuples at FULL precision from `PeakIO`'s live buffer.

F250/F251 showed `summit_ctrl` is unreadable from source: `peak_content` is built by a
C-level-private method, and `PeakIO.add` cannot be rebound (immutable cdef class).
Both were the wrong lever. The quantity is instead readable *after* the fact, because
two things are public even though the code that computes them is not:

  * `CallerFromAlignments.call_peaks` is a plain `def` (`CallPeakUnit.py:1033`), so a
    Python subclass can override it -- and `Spy` construction itself works, since only
    *private* methods of a `@cython.cclass` are C-level (F251b).
  * `PeakIO.peaks` is `cython.declare(dict, visibility="public")` (`PeakIO.py:236`),
    so after `call_peaks` returns, the buffer holds upstream's own tuples unrounded.

Since `CallPeakUnit.py:1385` sets
    fold_change = (summit_treat + pseudocount) / (summit_ctrl + pseudocount)
and `PeakIO.add_PeakContent` receives the summit treatment value as `pileup`, the
summit control value follows exactly:
    summit_ctrl = (pileup + pseudocount) / fold_change - pseudocount

This replaces inverting the `%.6g` `fold_enrichment` column, which produced F246's
wrong 4.2e-6 figure and F247's retraction.

Usage:
    oracle/dump_peak_scores.py <fixture-dir> <out.tsv> [--gsize N] [--extsize N]
"""
import os
import sys

sys.path.insert(0, os.environ.get("MACS3_SRC", "/scratch/mdra00001/tmp/opencode/macs3-src"))

import MACS3.Signal.CallPeakUnit as _cpu  # noqa: E402

_BASE = _cpu.CallerFromAlignments
_PSEUDOCOUNT = 1.0
_OUT = {"path": None, "rows": []}


def _emit(obj):
    """Flatten `PeakIO.peaks` -- `{chrom: {index: Peak}}` -- into rows."""
    buf = getattr(obj, "peaks", None)
    if not buf:
        return
    # Upstream's live buffer is `{chrom: [Peak, ...]}` -- a dict of *lists*, not of
    # dicts (verified: a dict-of-dicts assumption raises `'list' object has no
    # attribute 'items'`). Keep the index upstream assigned within each chromosome.
    if isinstance(buf, dict):
        items = [
            ((chrom.decode() if isinstance(chrom, bytes) else str(chrom)), i, pk)
            for chrom, peaks in buf.items()
            for i, pk in enumerate(peaks)
        ]
    else:
        items = [("?", i, pk) for i, pk in enumerate(buf)]
    for chrom, idx, pk in items:
        # F252: the objects ARE reachable (`PeakIO.peaks` is a public dict of
        # `{chrom: [PeakContent]}`, and `PeakDetect.call_peaks` is a public `def` that a
        # `callpeak_cmd.PeakDetect` wrapper does intercept), but every field is
        # C-private:
        #     dir(PeakContent)  ->  []
        #     pk.pileup         ->  AttributeError
        # So the tuple values are not readable either. This is the third and last
        # independent confirmation that `summit_ctrl` is unobservable (after F250's
        # immutable `PeakIO.add` and F251's C-level-private `__close_peak_wo_subpeaks`),
        # and the capture below therefore records only *counts*, which is what remains
        # genuinely observable from the live buffer.
        # No `repr(pk)`: it embeds a memory address, which is non-deterministic and
        # would make a recorded artifact un-diffable across runs.
        _OUT["rows"].append((chrom, int(idx), type(pk).__name__))


class Spy:
    """Wraps a real `PeakDetect` and harvests `peaks` after `call_peaks`.

    F251c: subclassing `CallerFromAlignments` does NOT work even for the *public*
    `call_peaks` -- Cython calls `self.call_peaks(...)` through a static vtable slot
    from compiled code, so a Python override never dispatches (verified: `Spy.__init__`
    fires, 0 rows). `PeakDetect` is the right level instead: `callpeak_cmd` is an
    ordinary Python module, so rebinding `callpeak_cmd.PeakDetect` really does change
    what gets constructed (the trick `dump_stages.py` already relies on), and
    `PeakDetect.call_peaks` is a public `def` that Python dispatch *does* reach.
    `self.peaks` is an instance attribute holding the finished `PeakIO`.
    """

    def __init__(self, inner):
        object.__setattr__(self, "_inner", inner)

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_inner"), name)

    def call_peaks(self, *a, **k):
        inner = object.__getattribute__(self, "_inner")
        result = inner.call_peaks(*a, **k)
        _emit(getattr(inner, "peaks", None) or result)
        return result


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("fixture_dir")
    ap.add_argument("out")
    ap.add_argument("--gsize", type=int, default=2000000)
    ap.add_argument("--extsize", type=int, default=200)
    args = ap.parse_args()

    _OUT["path"] = args.out
    from MACS3.Commands import callpeak_cmd

    real_pd = callpeak_cmd.PeakDetect
    callpeak_cmd.PeakDetect = lambda *a, **k: Spy(real_pd(*a, **k))

    src = os.environ.get("MACS3_SRC", "/scratch/mdra00001/tmp/opencode/macs3-src")
    import importlib.machinery
    import importlib.util

    loader = importlib.machinery.SourceFileLoader("macs3_bin", os.path.join(src, "bin", "macs3"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    parser = mod.prepare_argparser()

    argv = [
        "callpeak", "-n", "ps", "-g", str(args.gsize),
        "--outdir", "/tmp",
        "-t", os.path.join(args.fixture_dir, "treat.bed"), "-f", "BED",
        "--nomodel", "--extsize", str(args.extsize),
        "--call-summits",
    ]
    ctrl = os.path.join(args.fixture_dir, "ctrl.bed")
    if os.path.exists(ctrl):
        argv += ["-c", ctrl]
    try:
        callpeak_cmd.run(parser.parse_args(argv))
    except SystemExit:
        pass

    with open(args.out, "w") as fh:
        fh.write("chrom\tidx\ttype\n")
        for row in _OUT["rows"]:
            fh.write("\t".join(str(x) for x in row) + "\n")
    if not _OUT["rows"]:
        sys.exit("captured 0 peaks: `call_peaks` override did not fire (see F251b).")
    opaque = dir(type(items[0][2])) if items else []
    print(f"wrote {args.out} ({len(_OUT['rows'])} PeakContent objects)")
    if not any(not a.startswith("_") for a in opaque):
        print(
            "NOTE (F252): PeakContent exposes no public fields or methods, so only the "
            "peak COUNT and object identities are capturable -- not summit_ctrl.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
