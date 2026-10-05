#!/usr/bin/env python3
"""Run the pinned MACS3 over the fixture corpus and capture the golden outputs.

This is the authority for the differential tests. For every fixture it runs the
same command line that the Rust CLI will be given, with the outputs MACS3
supports for that fixture's mode, and records:

* every output file verbatim, plus its SHA-256
* the command line, MACS3's version, and the effective genome size used
* the `*_peaks.narrowPeak`, `*_peaks.xls`, `*_summits.bed` and bedGraph contents
  in a normalised form (`manifest.json` + a `parsed/` directory) so the Rust
  comparator does not have to re-parse BED itself
* stderr, so a fixture that makes upstream warn or fail is recorded as such
  rather than silently skipped

It also probes the option matrix (q-value cutoffs, `--call-summits`, broad
calling, `--keep-dup`, `--nolambda`, `--nomodel`) for the fixtures that support
it, because most of the compatibility surface is in the options, not in the data.

# The recording is path-independent

A run names the directory it was made in: `command.json` stores the argv, and every
`*_peaks.xls` header echoes the input paths. Committed verbatim, the corpus would
only replay from the one absolute path the recording machine had, so **everything
recorded here is written with the neutral tokens of `oracle/run_golden.py`** --
`<ROOT>` for this checkout, `macs3` for the interpreter -- in the argv, in the log
tails, and in the recorded bytes themselves. `run_golden.py` expands the token back
to whatever checkout it is running in; nothing downstream has to know where the
recording happened. `oracle/relocate_golden.py` did this once for the corpus that
was recorded before the rule existed, and `--check` is what keeps it true.

Usage:
    /path/to/macs3-venv/bin/python oracle/run_oracle.py \\
        --fixtures tests/fixtures --golden tests/golden
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_fixtures import GENOMES  # noqa: E402
from run_golden import neutralise, recording_substitutions  # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def neutralise_file(path, subs):
    """Rewrite one recorded output in place; True when its bytes changed."""
    with open(path, "rb") as fh:
        old = fh.read()
    new = neutralise(old, subs)
    if new == old:
        return False
    with open(path, "wb") as fh:
        fh.write(new)
    return True


def read_manifest(fixture_dir):
    m = {}
    path = os.path.join(fixture_dir, "manifest.tsv")
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            m[parts[0]] = parts[1] if len(parts) > 1 else ""
            if len(parts) > 2:
                m[parts[0]] = parts[1:]
    return m


def effective_gsize(fixture_dir, override=None):
    if override:
        return override
    total = 0
    with open(os.path.join(fixture_dir, "genome.txt")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            total += int(line.split()[1])
    return total


def macs3_version(macs3_bin):
    try:
        out = subprocess.run([macs3_bin, "--version"], capture_output=True, text=True, timeout=120)
        return (out.stdout + out.stderr).strip().splitlines()[-1]
    except Exception as exc:  # noqa: BLE001
        return f"unknown ({exc})"


def run_one(macs3_bin, fixture_dir, out_dir, opts, gsize, subs, timeout=900):
    os.makedirs(out_dir, exist_ok=True)
    name = os.path.basename(opts.get("name", "run"))
    cmd = [macs3_bin, "callpeak", "-n", name, "-g", str(gsize), "--outdir", out_dir]
    cmd += opts["args"]
    started = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        rc, out, err = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        rc, out, err = -1, "", "TIMEOUT"
    elapsed = time.time() - started

    # Neutralise the recorded bytes *before* hashing them: a recorded digest is a
    # claim about the committed bytes, so it has to describe the tokenised file.
    files = {}
    for fn in sorted(os.listdir(out_dir)):
        p = os.path.join(out_dir, fn)
        if os.path.isfile(p) and fn != "command.json":
            neutralise_file(p, subs)
            files[fn] = sha256(p)

    record = {
        "fixture": opts.get("fixture", ""),
        "variant": opts.get("variant", ""),
        "mode": opts.get("mode", "se"),
        "command": [neutralise(t.encode(), subs).decode() for t in cmd],
        "returncode": rc,
        "stdout_tail": neutralise(out[-4000:].encode(), subs).decode(errors="replace"),
        "stderr_tail": neutralise(err[-4000:].encode(), subs).decode(errors="replace"),
        "wall_seconds": round(elapsed, 3),
        "files": files,
    }
    with open(os.path.join(out_dir, "command.json"), "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
    return record


def model_capable(fixture_dir):
    """Whether this fixture can build a fragment model (>= 100 paired peaks)."""
    m = read_manifest(fixture_dir)
    v = m.get("model_capable", ["1"])
    if isinstance(v, list):
        v = v[0]
    return str(v) == "1"


def variants_for(fixture_dir, mode, gsize):
    """The option matrix for one fixture.

    Kept small but chosen to cover the compatibility surface: q-value cutoffs,
    narrow vs broad, summits, duplicate handling, and the two "skip a stage"
    flags that change the whole code path.
    """
    treat = os.path.join(fixture_dir, "treat.bed" if mode == "se" else
                         ("treat.bedpe" if mode == "pe" else "treat.frag"))
    ctrl = os.path.join(fixture_dir, "ctrl.bed" if mode == "se" else
                        ("ctrl.bedpe" if mode == "pe" else "ctrl.frag"))
    fmt = {"se": "BED", "pe": "BEDPE", "frag": "FRAG"}[mode]
    have_ctrl = os.path.exists(ctrl)

    base = ["-t", treat]
    if have_ctrl:
        base += ["-c", ctrl]
    base += ["-f", fmt]
    if not model_capable(fixture_dir):
        # MACS3 exits 1 when it cannot find 100 paired +/- strand peaks, so a
        # fixture that is too small for the model must declare --nomodel. This is
        # part of the compatibility contract, not a workaround.
        base += ["--nomodel", "--extsize", "200"]

    out = []

    def v(variant, extra, **kw):
        d = dict(kw)
        d["args"] = base + extra
        d["variant"] = variant
        d["mode"] = mode
        out.append(d)

    v("default", [])
    v("q001", ["-q", "0.01"])
    v("q05", ["-q", "0.5"])
    v("call_summits", ["--call-summits"])
    v("broad", ["--broad", "--broad-cutoff", "0.1", "-q", "0.05"])
    if have_ctrl:
        v("keepdup1", ["--keep-dup", "1"])
        v("keepdup_all", ["--keep-dup", "all"])
        v("keepdup_auto", ["--keep-dup", "auto"])
        v("nolambda", ["--nolambda"])
        v("slocal_500_llocal_2000", ["--slocal", "500", "--llocal", "2000"])
    v("nomodel_extsize", ["--nomodel", "--extsize", "200"])
    v("nomodel_shift", ["--nomodel", "--extsize", "200", "--shift", "-100"])
    v("shift_only", ["--shift", "-100", "--extsize", "200"])
    v("bw300", ["--bw", "300"])
    # upstream's --mfold is nargs=2 with default [5, 50] (bin/macs3:253), not
    # three values; `--mfold 2 4 8` is a usage error (exit 2) and that behaviour
    # is itself part of the CLI contract
    v("mfold_3_20", ["--mfold", "3", "20"])
    v("mfold_bad_arity", ["--mfold", "2", "4", "8"])
    v("spmr", ["--SPMR"])
    v("scale_to_large", ["--scale-to", "large"])
    v("B", ["-B"])
    v("gsize_numeric", [])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--macs3", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "venv", "bin", "macs3"))
    ap.add_argument("--fixtures", default="tests/fixtures")
    ap.add_argument("--golden", default="tests/golden")
    ap.add_argument("--only", default=None, help="fixture directory name filter")
    ap.add_argument("--variants", default=None, help="comma-separated variant filter")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    macs3_bin = os.path.abspath(args.macs3)
    if not os.path.exists(macs3_bin):
        sys.stderr.write(f"macs3 not found at {macs3_bin}\n")
        return 2

    fixtures_root = os.path.abspath(args.fixtures)
    golden_root = os.path.abspath(args.golden)
    os.makedirs(golden_root, exist_ok=True)

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # `<ROOT>` stands for the checkout, so a recording is only path-independent if the
    # paths it prints really are under it. Say so rather than writing a corpus that
    # only compares clean by accident.
    for label, path in (("--fixtures", fixtures_root), ("--golden", golden_root)):
        if os.path.commonpath([repo_root, path]) != repo_root:
            sys.stderr.write(
                f"run_oracle: {label} ({path}) is outside the checkout ({repo_root}); "
                f"recorded paths outside it cannot be written as <ROOT>\n"
            )
            return 2
    subs = recording_substitutions(repo_root, interpreter=macs3_bin,
                                   macs3_src=os.environ.get("MACS3_SRC"))

    version = macs3_version(macs3_bin)
    sys.stderr.write(f"oracle: {version}\n")

    summary = {"macs3": version, "fixtures": {}}
    total = 0
    variant_filter = set(args.variants.split(",")) if args.variants else None

    for group in sorted(os.listdir(fixtures_root)):
        gdir = os.path.join(fixtures_root, group)
        if not os.path.isdir(gdir):
            continue
        for name in sorted(os.listdir(gdir)):
            fdir = os.path.join(gdir, name)
            if not os.path.isfile(os.path.join(fdir, "manifest.tsv")):
                continue
            if args.only and args.only not in f"{group}/{name}":
                continue
            manifest = read_manifest(fdir)
            mode = manifest.get("mode", ["se"])[0] if isinstance(manifest.get("mode"), list) \
                else manifest.get("mode", "se")
            gsize = effective_gsize(fdir)
            key = f"{group}/{name}"
            summary["fixtures"][key] = {"mode": mode, "gsize": gsize, "variants": {}}
            for v in variants_for(fdir, mode, gsize):
                if variant_filter and v["variant"] not in variant_filter:
                    continue
                out_dir = os.path.join(golden_root, group, name, v["variant"])
                if os.path.isdir(out_dir):
                    shutil.rmtree(out_dir)
                v["fixture"] = key
                v["name"] = f"{name}_{v['variant']}"
                rec = run_one(macs3_bin, fdir, out_dir, v, gsize, subs)
                summary["fixtures"][key]["variants"][v["variant"]] = {
                    "returncode": rec["returncode"],
                    "files": rec["files"],
                    "wall_seconds": rec["wall_seconds"],
                }
                total += 1
                sys.stderr.write(
                    f"  {key}/{v['variant']}: rc={rec['returncode']} "
                    f"{len(rec['files'])} files {rec['wall_seconds']}s\n")
                sys.stderr.flush()
                if args.limit and total >= args.limit:
                    break
            if args.limit and total >= args.limit:
                break
        if args.limit and total >= args.limit:
            break

    with open(os.path.join(golden_root, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    sys.stderr.write(f"wrote {total} golden runs to {golden_root}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
