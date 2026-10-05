#!/usr/bin/env python3
"""Capture upstream's `peak_content` tuples -- the summit argmax's actual input.

Every port-side hypothesis about the summit residual (F46-F59g) has been
eliminated by measurement: the treatment pileup, control lambda, p-score,
p->q walk and histogram are all verified against upstream's own numbers on a
real fixture. What has never been captured is the thing the summit is computed
*from*: the per-entry `(tstart, tend, ttreat_p, tctrl_p, ti)` tuples that
`__close_peak_wo_subpeaks` sweeps with

    if not summit_value or summit_value < tscore:
        tsummit = [(tend + tstart) // 2,]
    midindex = (len(tsummit) + 1) // 2 - 1

`CallerFromAlignments` is a `@cython.cclass`, and it *is* subclassable, but
`__close_peak_wo_subpeaks` is NOT interceptable (see F251 for the mechanism), so
this harness cannot currently produce the rows. It is kept because the subclassing
and import-ordering machinery below is correct and reusable, and because it now
fails loudly instead of silently writing a header-only file.

Writes one TSV row per chunk: region_index, peak_index, tstart, tend, ttreat_p,
tctrl_p, ti.

Usage:
    oracle/dump_peak_content.py <fixture-dir> <out.tsv> [--gsize N]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oracle_env import require_src  # noqa: E402

sys.path.insert(0, require_src())

import MACS3.Signal.CallPeakUnit as _cpu  # noqa: E402

_BASE = _cpu.CallerFromAlignments
_OUT = {"path": None, "rows": [], "region": 0, "peak": 0}

# The private method's real name, after Python name mangling.
_CLOSE = "_CallerFromAlignments__close_peak_wo_subpeaks"


def _dump(peak_content, *args, **kwargs):
    chrom = args[2] if len(args) > 2 else kwargs.get("chrom", b"?")
    for (tstart, tend, ttreat, tctrl, ti) in peak_content:
        _OUT["rows"].append(
            (
                int(_OUT["region"]),
                int(_OUT["peak"]),
                int(tstart),
                int(tend),
                float(ttreat),
                float(tctrl),
                int(ti),
                chrom.decode() if isinstance(chrom, bytes) else str(chrom),
            )
        )
    _OUT["peak"] += 1
    return getattr(_BASE, _CLOSE)(self, peak_content, *args, **kwargs)


class Spy(_BASE):
    """`CallerFromAlignments` that records every `peak_content` it closes."""

    def __init__(self, *a, **k):
        if os.environ.get("MACS3_RS_SPY_DEBUG"):
            sys.stderr.write("!!! Spy.__init__ CALLED\n")
        super().__init__(*a, **k)
        _OUT["region"] += 1
        _OUT["peak"] = 0


setattr(Spy, _CLOSE, _dump)


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("fixture_dir")
    ap.add_argument("out")
    ap.add_argument("--gsize", type=int, default=2000000)
    ap.add_argument("--extsize", type=int, default=200)
    args = ap.parse_args()

    _cpu.CallerFromAlignments = Spy
    import MACS3.Signal.PeakDetect  # noqa: F401  (cached import must see Spy)

    from MACS3.Commands import callpeak_cmd

    import importlib.machinery
    import importlib.util

    _OUT["path"] = args.out
    # F251: this used to rebind `callpeak_cmd.CallerFromAlignments`, which is a no-op
    # -- `callpeak_cmd` never mentions the name, so `Spy` was never constructed and the
    # capture silently wrote 0 rows. The real construction sites are
    #   MACS3/Signal/PeakDetect.py:234 (with control) and :350 (without),
    # both written as `from MACS3.Signal.CallPeakUnit import CallerFromAlignments` at
    # *import* time, so the name is bound into `PeakDetect`'s own namespace.
    # Rebinding it there intercepts the real object without touching upstream's type.
    # F251b: `MACS3.Signal.PeakDetect` is a compiled `.so`, and Cython caches a
    # module-level `from X import Y` in a *C* global during module init. Rebinding
    # `PeakDetect.CallerFromAlignments` at runtime is therefore invisible to the
    # compiled code -- verified: the attribute reads back as `Spy` yet the real
    # caller is still built. The only interception point that works is rebinding
    # `CallPeakUnit.CallerFromAlignments` *before* `PeakDetect` is first imported,
    # so its cached import captures `Spy`.
    _cpu.CallerFromAlignments = Spy
    import MACS3.Signal.PeakDetect  # noqa: F401  (must follow the patch)

    src = require_src()
    loader = importlib.machinery.SourceFileLoader("macs3_bin", os.path.join(src, "bin", "macs3"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    parser = mod.prepare_argparser()
    treat = os.path.join(args.fixture_dir, "treat.bed")
    ctrl = os.path.join(args.fixture_dir, "ctrl.bed")
    argv = [
        "callpeak", "-n", "pc", "-g", str(args.gsize),
        "--outdir", "/tmp", "-t", treat, "-f", "BED",
        "--nomodel", "--extsize", str(args.extsize),
    ]
    if os.path.exists(ctrl):
        argv += ["-c", ctrl]
    try:
        callpeak_cmd.run(parser.parse_args(argv))
    except SystemExit:
        pass

    with open(args.out, "w") as fh:
        fh.write("region\tpeak\ttstart\ttend\tttreat\ttctrl\tti\tchrom\n")
        for row in _OUT["rows"]:
            fh.write("\t".join(str(x) for x in row) + "\n")
    if not _OUT["rows"]:
        sys.exit(
            "captured 0 peak_content chunks. This is NOT a rebinding bug -- the "
            "Spy subclass is constructed (see F251b, which fixed the import order). "
            "`__close_peak_wo_subpeaks` is a *private* method of a `@cython.cclass`, "
            "and Cython gives those C-level names: "
            "    getattr(CallerFromAlignments, '_CallerFromAlignments__close_peak_wo_subpeaks', None)\n"
            "    -> None\n"
            "so it is neither readable nor overridable, and subclass dispatch never "
            "reaches it. Together with F250 (immutable `PeakIO.add`) this means "
            "`peak_content` / `summit_ctrl` have no observational handle inside the "
            "oracle. Do not silently accept the empty capture."
        )
    print(f"wrote {args.out} ({len(_OUT['rows'])} chunks, {_OUT['region']} regions)")


if __name__ == "__main__":
    main()
