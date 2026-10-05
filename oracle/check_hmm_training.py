#!/usr/bin/env python3
"""Fresh differential check for HMMRATAC self-training and decoded regions.

Runs the pinned MACS3 source tree and this checkout on the same paired-end
input. It checks both Gaussian and Poisson training, training artifacts, model
parameters, and a second inference pass using the newly written models.

Example:
    python3 oracle/check_hmm_training.py --input /tmp/atr/yeast_500k_SRR1822137.bedpe.gz

Build the Rust binary first with ``cargo build --release -p macs-cli``. ``--input``,
``--oracle-src``, ``--venv``, and ``--binary`` can be supplied for local layouts.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oracle_env import require_python, require_src  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MODEL_ATOL = 1e-8
TRAINING_ATOL = 5e-12
MIN_ACCESSIBLE_BASE_JACCARD = 0.98


def oracle_source(explicit):
    """The pinned checkout, or exit saying how to get one. See oracle_env.py."""
    return require_src(explicit or None)


def run(command, cwd=None, env=None):
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True)
    if result.returncode:
        print("command failed:", " ".join(command), file=sys.stderr)
        print(result.stdout, file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"command exited {result.returncode}")
    return result


def parse_training_data(path):
    rows = []
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) != 6:
                raise AssertionError(f"{path}:{line_no}: expected six columns")
            rows.append((f[0], int(f[1]), tuple(map(float, f[2:]))))
    return rows


def compare_training_data(a_path, b_path):
    a, b = parse_training_data(a_path), parse_training_data(b_path)
    if len(a) != len(b):
        raise AssertionError(f"training row counts differ: {len(a)} vs {len(b)}")
    max_abs = 0.0
    for i, (x, y) in enumerate(zip(a, b)):
        if x[:2] != y[:2]:
            raise AssertionError(f"training coordinate differs at row {i}: {x[:2]} vs {y[:2]}")
        for xv, yv in zip(x[2], y[2]):
            max_abs = max(max_abs, abs(xv - yv))
    if max_abs > TRAINING_ATOL:
        raise AssertionError(f"training values differ by {max_abs:.3g} (limit {TRAINING_ATOL})")
    return len(a), max_abs


def flatten_numbers(value):
    if isinstance(value, list):
        out = []
        for item in value:
            out.extend(flatten_numbers(item))
        return out
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return [float(value)]
    return []


def compare_models(a_path, b_path):
    with open(a_path) as fh:
        a = json.load(fh)
    with open(b_path) as fh:
        b = json.load(fh)
    for key in ("hmm_type", "n_features", "hmm_binsize", "i_open_region",
                "i_background_region", "i_nucleosomal_region"):
        if a.get(key) != b.get(key):
            raise AssertionError(f"model field {key} differs: {a.get(key)!r} vs {b.get(key)!r}")
    fields = ("startprob", "transmat", "means", "covars") if a["hmm_type"] == "gaussian" else (
        "startprob", "transmat", "lambdas")
    maximum = 0.0
    for key in fields:
        x, y = flatten_numbers(a[key]), flatten_numbers(b[key])
        if len(x) != len(y):
            raise AssertionError(f"model field {key} shape differs: {len(x)} vs {len(y)} values")
        maximum = max(maximum, *(abs(u - v) for u, v in zip(x, y)))
    if not math.isfinite(maximum) or maximum > MODEL_ATOL:
        raise AssertionError(f"model parameter difference {maximum:.3g} exceeds {MODEL_ATOL}")
    return maximum


def merge_intervals(path):
    by_chrom = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            by_chrom.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    merged = {}
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        out = []
        for start, end in intervals:
            if out and start <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], end))
            else:
                out.append((start, end))
        merged[chrom] = out
    return merged


def accessible_base_jaccard(a_path, b_path):
    a, b = merge_intervals(a_path), merge_intervals(b_path)
    intersection = union = 0
    for chrom in set(a) | set(b):
        av, bv = a.get(chrom, []), b.get(chrom, [])
        total_a = sum(e - s for s, e in av)
        total_b = sum(e - s for s, e in bv)
        overlap = 0
        i = j = 0
        while i < len(av) and j < len(bv):
            overlap += max(0, min(av[i][1], bv[j][1]) - max(av[i][0], bv[j][0]))
            if av[i][1] <= bv[j][1]:
                i += 1
            else:
                j += 1
        intersection += overlap
        union += total_a + total_b - overlap
    return (intersection / union if union else 1.0), sum(map(len, a.values())), sum(map(len, b.values()))


def compare_bytes(a, b, label):
    if open(a, "rb").read() != open(b, "rb").read():
        raise AssertionError(f"{label} differs byte-for-byte")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=None)
    ap.add_argument("--oracle-src", default=None)
    ap.add_argument("--venv", default=None)
    ap.add_argument("--binary", default=os.path.join(ROOT, "target/release/macs3-rs"))
    ap.add_argument("--gaussian-model", help="optional external Gaussian model for full-decoding checks")
    ap.add_argument("--poisson-model", help="optional external Poisson model for full-decoding checks")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    src = oracle_source(args.oracle_src)
    args.input = args.input or os.environ.get("MACS3_HMM_INPUT")
    if not args.input:
        candidates = ("/tmp/atr/yeast_500k_SRR1822137.bedpe.gz",
                      os.path.join(src, "test/yeast_500k_SRR1822137.bedpe.gz"))
        args.input = next((p for p in candidates if os.path.exists(p)), candidates[-1])
    oracle_python = require_python(
        os.path.join(args.venv, "bin", "python") if args.venv else None
    )
    for path in (args.input, oracle_python, args.binary):
        if not os.path.exists(path):
            print(f"missing required path: {path}", file=sys.stderr)
            return 2

    work = args.outdir or tempfile.mkdtemp(prefix="hmm-training-")
    os.makedirs(work, exist_ok=True)
    oracle_env = dict(os.environ, PYTHONPATH=src)
    summary = []
    for hmm_type in ("gaussian", "poisson"):
        print(f"== {hmm_type} self-training ==", file=sys.stderr)
        odir = os.path.join(work, f"oracle-{hmm_type}")
        rdir = os.path.join(work, f"rust-{hmm_type}")
        os.makedirs(odir, exist_ok=True)
        os.makedirs(rdir, exist_ok=True)
        prefix = f"{hmm_type}"
        common = ["hmmratac", "-i", args.input, "-f", "BEDPE", "--hmm-type", hmm_type,
                  "--modelonly", "--save-training-data", "--save-digested", "-n", prefix]
        run([oracle_python, "bin/macs3", *common, "--outdir", odir], cwd=src, env=oracle_env)
        run([args.binary, *common, "--outdir", rdir])
        for suffix in ("training_regions.bed", "training_lengths.txt", "digested_short.bdg",
                       "digested_mono.bdg", "digested_di.bdg", "digested_tri.bdg",
                       "cutoff_analysis.tsv"):
            compare_bytes(os.path.join(odir, f"{prefix}_{suffix}"),
                          os.path.join(rdir, f"{prefix}_{suffix}"), suffix)
        rows, row_delta = compare_training_data(
            os.path.join(odir, f"{prefix}_training_data.txt"),
            os.path.join(rdir, f"{prefix}_training_data.txt"))
        model_delta = compare_models(os.path.join(odir, f"{prefix}_model.json"),
                                     os.path.join(rdir, f"{prefix}_model.json"))

        print(f"== {hmm_type} full decode ==", file=sys.stderr)
        odir_full = os.path.join(work, f"oracle-{hmm_type}-decode")
        rdir_full = os.path.join(work, f"rust-{hmm_type}-decode")
        os.makedirs(odir_full, exist_ok=True)
        os.makedirs(rdir_full, exist_ok=True)
        supplied_model = getattr(args, f"{hmm_type}_model")
        omodel = supplied_model or os.path.join(odir, f"{prefix}_model.json")
        rmodel = supplied_model or os.path.join(rdir, f"{prefix}_model.json")
        decode_common = ["hmmratac", "-i", args.input, "-f", "BEDPE", "--hmm-type", hmm_type,
                         "--save-states", "--save-likelihoods", "--model", omodel, "-n", prefix]
        run([oracle_python, "bin/macs3", *decode_common, "--outdir", odir_full], cwd=src, env=oracle_env)
        rust_common = decode_common.copy()
        rust_common[rust_common.index(omodel)] = rmodel
        run([args.binary, *rust_common, "--outdir", rdir_full])
        stem = f"{prefix}_accessible_regions.narrowPeak"
        jac, acount, bcount = accessible_base_jaccard(os.path.join(odir_full, stem),
                                                       os.path.join(rdir_full, stem))
        if jac < MIN_ACCESSIBLE_BASE_JACCARD:
            raise AssertionError(f"{hmm_type} accessible base Jaccard {jac:.5f} below {MIN_ACCESSIBLE_BASE_JACCARD}")
        compare_bytes(os.path.join(odir_full, f"{prefix}_states.bed"),
                      os.path.join(rdir_full, f"{prefix}_states.bed"), "states BED")
        summary.append((hmm_type, rows, row_delta, model_delta, jac, acount, bcount))

    # The upstream command-line suite also exercises the counted, barcode-filtered
    # FRAG input path. It uses PETrackII's MT19937 sampler, unlike BAMPE/BEDPE's
    # PETrackI path, and is the regression case for per-count downsampling.
    frag_input = os.path.join(src, "test", "test.fragments.tsv.gz")
    barcode_file = os.path.join(src, "test", "barcodes.txt")
    for path in (frag_input, barcode_file):
        if not os.path.exists(path):
            raise FileNotFoundError(f"missing pinned scATAC fixture: {path}")
    print("== scATAC Poisson self-training ==", file=sys.stderr)
    odir = os.path.join(work, "oracle-scatac")
    rdir = os.path.join(work, "rust-scatac")
    os.makedirs(odir, exist_ok=True)
    os.makedirs(rdir, exist_ok=True)
    prefix = "scatac"
    common = ["hmmratac", "-i", frag_input, "-f", "FRAG", "--barcodes", barcode_file,
              "--hmm-type", "poisson", "--modelonly", "--save-training-data",
              "--save-digested", "-n", prefix]
    run([oracle_python, "bin/macs3", *common, "--outdir", odir], cwd=src, env=oracle_env)
    run([args.binary, *common, "--outdir", rdir])
    for suffix in ("training_regions.bed", "training_lengths.txt", "digested_short.bdg",
                   "digested_mono.bdg", "digested_di.bdg", "digested_tri.bdg",
                   "cutoff_analysis.tsv"):
        compare_bytes(os.path.join(odir, f"{prefix}_{suffix}"),
                      os.path.join(rdir, f"{prefix}_{suffix}"), f"scATAC {suffix}")
    rows, row_delta = compare_training_data(
        os.path.join(odir, f"{prefix}_training_data.txt"),
        os.path.join(rdir, f"{prefix}_training_data.txt"))
    model_delta = compare_models(os.path.join(odir, f"{prefix}_model.json"),
                                 os.path.join(rdir, f"{prefix}_model.json"))

    print("== scATAC Poisson full decode ==", file=sys.stderr)
    odir_full = os.path.join(work, "oracle-scatac-decode")
    rdir_full = os.path.join(work, "rust-scatac-decode")
    os.makedirs(odir_full, exist_ok=True)
    os.makedirs(rdir_full, exist_ok=True)
    omodel = os.path.join(odir, f"{prefix}_model.json")
    rmodel = os.path.join(rdir, f"{prefix}_model.json")
    decode_common = ["hmmratac", "-i", frag_input, "-f", "FRAG", "--barcodes", barcode_file,
                     "--hmm-type", "poisson", "--save-states", "--save-likelihoods",
                     "--model", omodel, "-n", prefix]
    run([oracle_python, "bin/macs3", *decode_common, "--outdir", odir_full],
        cwd=src, env=oracle_env)
    rust_common = decode_common.copy()
    rust_common[rust_common.index(omodel)] = rmodel
    run([args.binary, *rust_common, "--outdir", rdir_full])
    stem = f"{prefix}_accessible_regions.narrowPeak"
    jac, acount, bcount = accessible_base_jaccard(os.path.join(odir_full, stem),
                                                   os.path.join(rdir_full, stem))
    if jac < MIN_ACCESSIBLE_BASE_JACCARD:
        raise AssertionError(f"scATAC accessible base Jaccard {jac:.5f} below {MIN_ACCESSIBLE_BASE_JACCARD}")
    compare_bytes(os.path.join(odir_full, f"{prefix}_states.bed"),
                  os.path.join(rdir_full, f"{prefix}_states.bed"), "scATAC states BED")
    summary.append(("scATAC-Poisson", rows, row_delta, model_delta, jac, acount, bcount))

    print("type      bins  training-max-abs  model-max-abs  accessible-base-J  peaks-oracle/rust")
    for typ, rows, row_delta, model_delta, jac, a, b in summary:
        print(f"{typ:<8} {rows:>6} {row_delta:>17.3g} {model_delta:>14.3g} {jac:>18.6f} {a}/{b}")
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, RuntimeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
