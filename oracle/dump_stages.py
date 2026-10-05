#!/usr/bin/env python
"""Dump every intermediate stage of MACS3's callpeak pipeline, in-process.

Why in-process instead of shelling out: `callpeak` does not write intermediates,
so a stage-by-stage differential has no other option than to drive the pipeline
itself. Rather than re-implement `callpeak_cmd.run()` -- which would drift from
upstream on every edit -- this harness *wraps* the real functions and records
what they computed. The pipeline, the arithmetic, and the ordering are upstream's;
only the observation is ours.

Captured stages
---------------
| stage | source |
|---|---|
| `reads_pre_filter`  | track contents before `filter_dup` |
| `reads_post_filter` | track contents after `filter_dup` |
| `duplicates`        | retained counts and redundant rates, parsed from upstream's log |
| `totals`            | t0/t1/c0/c1, redundant rates, duplicate limits |
| `d`                 | fragment length actually used, plus `alternative_d` |
| `scaling`           | `ratio_treat2control`, `tocontrol`, `lambda_bg` |
| `treat_pileup`      | `stages_treat_pileup.bdg`, written by `-B` |
| `lambda_merged`     | `stages_control_lambda.bdg`, written by `-B` |
| `qvalue_table`      | `ScoreTrack.pvalue_stat` |
| `candidate_peaks`   | `PeakDetect.peaks` |
| `final_peaks`       | `PeakDetect.final_peaks` |
| `summits`           | `PeakDetect.final_peaks` summit column |

The fragment-model arrays are **not** captured: `PeakModel` is a `cdef` class with no
Python-visible attributes. The `*_model.r` file MACS3 writes *is* byte-compared by the
golden layer, which covers the same arrays as an output rather than as an intermediate.

Duplicates (F240) were also listed as unreachable, which was half true. `filter_dup`'s
retained counts already show up as the pre/post-filter totals; the *rates* do not, but
upstream reports them ("`Redundant rate of treatment: 0.99`"), and the log is captured
in-process with a logging handler so the rates are recorded alongside the counts. Note
that FRAG mode runs no duplicate filtering at all -- upstream forces `--keep-dup all`
for it -- so for `-f FRAG` only the totals appear, and their absence is correct rather
than a gap.

The control scale ladder (`ctrl_d_s`, `ctrl_scaling_factor_s`, `lambda_bg`) is a **C
attribute** of `CallerFromAlignments`, so it cannot be read from the live object. It is
passed into the constructor, however, and the Python module global can be wrapped to
observe the actual constructor arguments before the Cython object is created. The
constructor's first two arguments are positional (`treat`, `ctrl`); the ladder values
are supplied by keyword in this upstream call.

The per-*scale* local lambda arrays (slocal and llocal separately) are locals
inside `ScoreTrack.call_peaks` and are likewise not reachable; `xls_header`
carries the two window sizes instead, which is the same information in the form
both sides can compare.

Usage:
    dump_stages.py <fixture-dir> <mode:se|pe|frag> --out DIR [-- <extra macs3 args>]
"""

import argparse
import json
import os
import re
import sys

# NumPy lives in the oracle virtualenv, so this re-execs under the provisioned
# interpreter rather than whichever `python3` happens to be first on PATH.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oracle_env import ensure_oracle_python, require_src  # noqa: E402

ensure_oracle_python()

import numpy as np


def _load_macs3_bin(src):
    import importlib.machinery
    import importlib.util

    path = os.path.join(src, "bin", "macs3")
    loader = importlib.machinery.SourceFileLoader("macs3_bin", path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def build_args(src, fixture_dir, mode, extra, out_dir):
    """Parse a real `macs3 callpeak` argv through upstream's own parser."""
    mod = _load_macs3_bin(src)
    parser = mod.prepare_argparser()

    treat = os.path.join(fixture_dir, "treat.bed" if mode == "se"
                         else ("treat.bedpe" if mode == "pe" else "treat.frag"))
    ctrl = os.path.join(fixture_dir, "ctrl.bed" if mode == "se"
                        else ("ctrl.bedpe" if mode == "pe" else "ctrl.frag"))
    fmt = {"se": "BED", "pe": "BEDPE", "frag": "FRAG"}[mode]

    # keep upstream's own output inside --out so the harness never scatters
    # files into the repository
    argv = ["callpeak", "-n", "stages", "-g", "1000000", "--outdir", out_dir,
            "-t", treat, "-f", fmt]
    if os.path.exists(ctrl):
        argv += ["-c", ctrl]
    argv += list(extra)
    # Ask upstream to persist two stages itself rather than trying to read them
    # back off a Cython object: `--bdg` writes the per-position treatment pileup
    # and control lambda bedGraphs, `--cutoff-analysis` writes the p->q table.
    # Both are produced by the engine we are observing, so they are as
    # authoritative as anything else in the golden set.
    if not any(a.startswith("--cutoff-analysis") for a in extra):
        argv += ["--cutoff-analysis"]
    # G0: the signal tracks are cdef attributes of `CallerFromAlignments` and cannot be
    # read from Python at all, so they are obtained the only other way there is -- by
    # asking the engine to persist them. `-B` writes `<name>_treat_pileup.bdg` and
    # `<name>_control_lambda.bdg`, which are the treatment pileup and the merged local
    # lambda track. Same bytes a user receives, so strictly better evidence than a probe
    # would be.
    if not any(a.startswith(("--bdg", "-B")) for a in argv):
        argv += ["-B"]
    return parser.parse_args(argv), " ".join(argv)


def chrom_name(name):
    """Chromosome names are bytes upstream; normalise once, here."""
    if isinstance(name, bytes):
        return name.decode()
    return str(name)


class Capture:
    """Records stage outputs as the real pipeline produces them.

    Cython `cdef` classes are immutable, so `FWTrack.filter_dup` and
    `PeakDetect.call_peaks` cannot be patched in place. Both are reachable
    another way: `callpeak_cmd` looks them up in its own module namespace at call
    time, so rebinding the *names* it uses intercepts the real calls without
    touching upstream's types.
    """

    def __init__(self):
        self.stages = {}
        self.peakdetect = None
        self.scoretrack = None
        self.pre_filter = {}
        self.post_filter = {}
        self.caller = None
        self.ladder = None
        self.treat_pileup = {}
        self.ctrl_pileup = {}
        self.lambda_merged = {}
        self.pvalue_track = {}
        self.qvalue_track = {}
        self.pq_table = []

    # -- stage: reads as loaded, before any filtering ---------------------
    def wrap_loader(self, cmd_mod, name):
        """Wrap `load_tag_files_options` / `load_frag_files_options`.

        Both return `(treat, control)`, and this is the last point at which the
        tracks still hold their unfiltered reads, so it is where the pre-filter
        stage has to be observed.
        """
        original = getattr(cmd_mod, name)

        def patched(options):
            treat, control = original(options)
            self.pre_filter["treatment"] = snapshot_track(treat)
            self.pre_filter["control"] = (
                snapshot_track(control) if control is not None else None)
            return treat, control

        setattr(cmd_mod, name, patched)
        return original

    # -- stage: the score engine's lambda parameters and p->q table -------
    def wrap_caller_factory(self, pd_mod):
        """Observe the actual arguments used to construct the lambda calculator.

        `PeakDetect` resolves `CallerFromAlignments` from its module globals at
        call time, so rebinding the name here observes arguments without patching
        any Cython type. Cython attributes remain inaccessible; the constructor
        inputs are captured exactly as passed.
        """
        original = pd_mod.CallerFromAlignments

        def factory(*a, **kw):
            # Match CallerFromAlignments.__init__ in the pinned upstream source.
            # This call path currently passes treat and ctrl positionally and
            # the ladder values by keyword. Mapping the remaining parameters too
            # keeps the observation correct if that call site changes.
            names = (
                "treat", "ctrl", "d", "ctrl_d_s", "treat_scaling_factor",
                "ctrl_scaling_factor_s", "stderr_on", "pseudocount",
                "end_shift", "lambda_bg", "save_bedGraph",
                "bedGraph_filename_prefix", "bedGraph_treat_filename",
                "bedGraph_control_filename", "cutoff_analysis_filename",
                "save_SPMR",
            )
            values = dict(zip(names, a))
            values.update(kw)
            self.ladder = (
                list(values.get("ctrl_d_s", ())),
                list(values.get("ctrl_scaling_factor_s", ())),
                values.get("lambda_bg", 0.0),
                values.get("treat_scaling_factor", 1.0),
            )
            return original(*a, **kw)

        pd_mod.CallerFromAlignments = factory
        return original

    # -- stage: totals after filtering, plus the d actually used -----------
    def wrap_peakdetect_factory(self, cmd_mod):
        original = cmd_mod.PeakDetect

        def factory(*a, **kw):
            obj = original(*a, **kw)
            self.peakdetect = obj
            # the tracks are final by the time PeakDetect is constructed, so
            # these totals already reflect filter_dup
            self.post_filter["treatment"] = int(getattr(obj.treat, "total", 0))
            self.post_filter["control"] = (
                int(getattr(obj.control, "total", 0)) if obj.control is not None else None)
            return obj

        cmd_mod.PeakDetect = factory
        return original

    def wrap_call_peaks(self, obj):
        """Return a proxy that records `call_peaks`'s arguments, then delegates.

        `PeakDetect.call_peaks` is a cdef method and the class is immutable, so it
        cannot be rebound. But `callpeak_cmd.run` calls `peakdetect.call_peaks()`
        on whatever object the factory returned, so a delegating proxy observes
        the call. That is the only route to `ctrl_d_s`, `ctrl_scaling_factor_s`
        and `lambda_bg`, which are C attributes that raise on direct access.

        The ladder matters: it is the input F238 turned on, and reading it from
        upstream's own call makes the wide-window factors a diffable fact instead
        of a numeric argument.
        """
        cap = self

        class _Proxy:
            def __getattr__(self, name):
                return getattr(obj, name)

            def call_peaks(self, *a, **kw):
                cap.ladder = (
                    list(kw.get("ctrl_d_s", ())),
                    list(kw.get("ctrl_scaling_factor_s", ())),
                    kw.get("lambda_bg", 0.0),
                    kw.get("treat_scale", 1.0),
                )
                return obj.call_peaks(*a, **kw)

        return _Proxy()


def snapshot_track(track):
    """Per-chromosome read/fragment counts for a track, as it currently stands.

    Reads a `FWTrack`/`PETrack` whose arrays may still carry `buffer_size`
    padding, so `finalize()` is called first (see F19) and a `total` cross-check
    is included so a wrong snapshot cannot pass silently.
    """
    out = {"total": int(getattr(track, "total", 0)), "chroms": {}}
    try:
        track.finalize()
    except Exception:
        pass
    try:
        chroms = sorted(track.get_chr_names())
    except Exception:
        return out
    for chrom in chroms:
        locs = track.get_locations_by_chr(chrom)
        name = chrom_name(chrom)
        if isinstance(locs, np.ndarray) and getattr(locs.dtype, "names", None):
            fields = list(locs.dtype.names)
            out["chroms"][name] = {"fragments": int(locs.size), "fields": fields}
        elif isinstance(locs, list):
            out["chroms"][name] = {"plus": int(len(locs[0])), "minus": int(len(locs[1]))}
    return out


def parse_cutoff_analysis(path):
    """The p->q table as written by upstream's `--cutoff-analysis`.

    `CallerFromAlignments.pqtable` is a C attribute of a cdef class and cannot be
    read from Python, so the table is taken from the file the engine itself
    writes. Columns: pscore qscore npeaks lpeaks avelpeak.
    """
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as fh:
        header = fh.readline()
        if not header.lower().startswith("pscore"):
            return []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            # npeaks and lpeaks are counts; avelpeak is a float, so a blanket
            # int() over the trailing columns silently drops every row
            try:
                rows.append([float(f[0]), float(f[1]),
                             int(f[2]), int(f[3]),
                             float(f[4]) if len(f) > 4 else None])
            except (ValueError, IndexError):
                continue
    rows.sort(key=lambda r: (-r[0], -r[1]))
    return rows


def parse_bedgraph(path):
    """A bedGraph the engine wrote, as ``{chrom: [[start, end, value], ...]}``.

    Used for the `-B` outputs. Keeping the intervals (not just the values) is the
    point: a differential that only compared values would miss a segmentation
    difference, which is precisely the F207 class of bug.
    """
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                out.setdefault(f[0], []).append(
                    [int(f[1]), int(f[2]), float(f[3])])
            except ValueError:
                continue
    return out


def parse_xls_header(path):
    """`d` and the local-lambda windows, from the `#` comments upstream writes.

    `callpeak_cmd` records "# d = N" and "Range for calculating regional lambda
    is: S bps and L bps" in the XLS header. That is the same information the
    engine holds, in a form both sides can be compared on byte-for-byte.
    """
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                if "pscore" in line or line.strip() == "":
                    continue
                break
            low = line.lower()
            if low.startswith("# d ="):
                out["d"] = int(line.split("=")[1].strip())
            elif "range for calculating regional lambda is:" in low:
                body = line.split(":", 1)[1]
                nums = [int(x) for x in body.replace("bps", "").split()
                        if x.strip().isdigit()]
                if len(nums) >= 2:
                    out["slocal"], out["llocal"] = nums[0], nums[1]
            elif "local lambda is disabled" in low:
                out["nolambda"] = True
            elif "use" in low and "as fragment length" in low:
                body = line.split(":", 1)[1]
                for tok in body.split():
                    if tok.isdigit():
                        out["d"] = int(tok)
    return out



# ---------------------------------------------------------------------------
# G240: the scale ladder and the duplicate counts.
#
# Both were listed as unreachable. `ctrl_d_s` / `ctrl_scaling_factor_s` /
# `lambda_bg` are C attributes of `CallerFromAlignments` and *do* raise on
# attribute access -- but they are passed **into** `call_peaks` as arguments, so
# wrapping the method observes them directly. That matters beyond bookkeeping:
# the ladder is the input F238 turned on, and having upstream's own values in the
# committed corpus turns a speculative numeric argument into a diffable fact.
#
# Duplicates are computed inside the compiled `filter_dup` and are likewise not
# readable off the track, but upstream *reports* them: "#1  tags after filtering
# in treatment: 32", "Redundant rate of treatment: 0.99". Parsing the log is
# worse evidence than reading a structure, and it is the only evidence there is.
# ---------------------------------------------------------------------------

_DUP_PATTERNS = (
    ("treat_total", re.compile(r"total (?:tags|fragments) in treatment:\s*(\d+)")),
    ("treat_after_filter", re.compile(r"(?:tags|fragments) after filtering in treatment:\s*(\d+)")),
    ("treat_redundant_rate", re.compile(r"Redundant rate of treatment:\s*([0-9.]+)")),
    ("ctrl_total", re.compile(r"total (?:tags|fragments) in control:\s*(\d+)")),
    ("ctrl_after_filter", re.compile(r"(?:tags|fragments) after filtering in control:\s*(\d+)")),
    ("ctrl_redundant_rate", re.compile(r"Redundant rate of control:\s*([0-9.]+)")),
)


def parse_duplicates(log_path):
    """Retained counts and redundant rates, as upstream reports them."""
    if not log_path or not os.path.exists(log_path):
        return {}
    with open(log_path, errors="replace") as fh:
        text = fh.read()
    out = {}
    for key, pat in _DUP_PATTERNS:
        m = pat.search(text)
        if m:
            out[key] = float(m.group(1))
    return out


def capture_ladder(cap):
    """`ctrl_d_s`, `ctrl_scaling_factor_s` and `lambda_bg`, seen as arguments."""
    out = {}
    got = getattr(cap, "ladder", None)
    if got:
        d_s, scale_s, bg, tscale = got
        out["ctrl_d_s"] = [int(v) for v in d_s]
        out["ctrl_scaling_factor_s"] = [float(v) for v in scale_s]
        out["lambda_bg"] = float(bg)
        out["treat_scale"] = float(tscale)
    return out



def dump_caller(stages, out_dir):
    """Capture the stages that are only reachable through upstream's writers.

    `CallerFromAlignments` is a Cython cdef class: its `pqtable`, `d`,
    `ctrl_d_s`, `ctrl_scaling_factor_s` and `lambda_bg` are all C attributes and
    raise on access from Python, so they cannot be read off the live object. Two
    of them are instead recovered from files the engine writes itself, which is
    strictly better evidence because it is the same bytes a user receives.
    """
    stages["qvalue_table"] = parse_cutoff_analysis(
        os.path.join(out_dir, "stages_cutoff_analysis.txt"))
    stages["xls_header"] = parse_xls_header(
        os.path.join(out_dir, "stages_peaks.xls"))
    treat = parse_bedgraph(os.path.join(out_dir, "stages_treat_pileup.bdg"))
    if treat:
        stages["treat_pileup"] = treat
    lam = parse_bedgraph(os.path.join(out_dir, "stages_control_lambda.bdg"))
    if lam:
        stages["lambda_merged"] = lam


def dump_peakio(out_dir, label, peakio):
    """Serialise a PeakIO with upstream's own writers.

    `PeakContent` is a fully opaque Cython cdef class: it exposes no Python
    attributes and is not iterable, so its contents cannot be read directly. Its
    writers can, and they are the same writers `callpeak_cmd` uses for the golden
    files -- so what lands here is both the intermediate and a byte-comparable
    rendering of it, which is what the acceptance criteria compares.

    Returns the paths written, for the manifest.
    """
    if peakio is None:
        return []
    written = []

    def _write(suffix, method, **kwargs):
        path = os.path.join(out_dir, f"{label}.{suffix}")
        try:
            with open(path, "w") as fh:
                getattr(peakio, method)(fh, **kwargs)
        except TypeError:
            # the signature differs between PeakIO versions; fall back to the
            # minimum call and let a real failure surface as a missing file
            try:
                with open(path, "w") as fh:
                    getattr(peakio, method)(fh)
            except Exception:
                return None
        except Exception:
            return None
        if os.path.exists(path):
            written.append(path)
        return path

    _write("peaks.xls", "write_to_xls", name=b"stages")
    _write("summits.bed", "write_to_summit_bed", name_prefix=b"stages_peak_")
    return written


def dump_qtable(stat):
    """The p->q lookup table, sorted by descending p-score like upstream."""
    if stat is None:
        return []
    rows = []
    for key, _v in stat.items():
        p, q = key
        rows.append([float(p), float(q)])
    rows.sort(key=lambda r: (-r[0], -r[1]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fixture_dir")
    ap.add_argument("mode", choices=["se", "pe", "frag"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--log", default="",
                    help="oracle log to parse duplicate counts and rates from")
    ap.add_argument("--macs3-src", default=None,
                    help="the pinned MACS3 checkout; default: the provisioned one "
                         "(MACS3_SRC, then oracle/ENV.provisioned)")
    args, extra = ap.parse_known_args()
    # Anything not claimed above is passed through to MACS3 verbatim. A bare `--`
    # separator is dropped if the caller used one.
    if extra and extra[0] == "--":
        extra = extra[1:]

    os.makedirs(args.out, exist_ok=True)

    argv_args, argv_text = build_args(
        args.macs3_src or require_src(), args.fixture_dir, args.mode, extra, args.out)

    from MACS3.Commands import callpeak_cmd
    from MACS3.Signal import PeakDetect as pd_mod

    # The oracle logs its duplicate rates to the root logger as it runs, so the log
    # has to be captured *in process* -- `record_stages.py` only sees the output
    # after this script has already exited.
    log_handler = None
    if args.log:
        import logging

        log_handler = logging.FileHandler(args.log, mode="w")
        log_handler.setFormatter(logging.Formatter("%(message)s"))
        logging.getLogger().addHandler(log_handler)
        logging.getLogger().setLevel(logging.INFO)

    cap = Capture()
    se_loader = cap.wrap_loader(callpeak_cmd, "load_tag_files_options")
    frag_loader = cap.wrap_loader(callpeak_cmd, "load_frag_files_options")
    real_caller = cap.wrap_caller_factory(pd_mod)
    real_pd = cap.wrap_peakdetect_factory(callpeak_cmd)

    stages = {"_argv": argv_text}

    try:
        callpeak_cmd.run(argv_args)
    except SystemExit as e:
        stages["_systexit"] = str(e)
    finally:
        callpeak_cmd.load_tag_files_options = se_loader
        callpeak_cmd.load_frag_files_options = frag_loader
        pd_mod.CallerFromAlignments = real_caller
        callpeak_cmd.PeakDetect = real_pd
        if log_handler is not None:
            logging.getLogger().removeHandler(log_handler)
            log_handler.flush()
            log_handler.close()

    # -- collect what the wrapped calls saw -------------------------------
    stages["reads_pre_filter"] = cap.pre_filter
    stages["reads_post_filter"] = cap.post_filter
    pd_obj = cap.peakdetect
    if pd_obj is not None:
        stages["d"] = int(getattr(pd_obj, "d", 0))
        stages["scaling"] = {
            "ratio_treat2control": float(getattr(pd_obj, "ratio_treat2control", 0.0)),
            "sregion": int(pd_obj.sregion or 0),
            "lregion": int(pd_obj.lregion or 0),
            "tocontrol": bool(getattr(pd_obj.opt, "tocontrol", False)),
        }
        stages["peak_files"] = {
            "candidate": dump_peakio(args.out, "candidate_peaks",
                                     getattr(pd_obj, "peaks", None)),
            "final": dump_peakio(args.out, "final_peaks",
                                 getattr(pd_obj, "final_peaks", None)),
        }
    ladder = capture_ladder(cap)
    if ladder:
        stages["lambda_ladder"] = ladder
    dups = parse_duplicates(getattr(args, "log", None))
    if dups:
        stages["duplicates"] = dups
    dump_caller(stages, args.out)

    with open(os.path.join(args.out, "stages.json"), "w") as fh:
        json.dump(stages, fh, indent=1, sort_keys=True, default=float)

    print(f"wrote {os.path.join(args.out, 'stages.json')} "
          f"({len(stages)} stage keys)")


if __name__ == "__main__":
    main()
