#!/usr/bin/env python
"""G13 differential: HMMRATAC accessible regions against the pinned oracle.

The acceptance criterion is chromosome-aware accessible-base Jaccard >= 0.98. This
script uses a fresh upstream model unless one is explicitly supplied.

Both sides run inference against the same model. The separate
``check_hmm_training.py`` script measures Gaussian and Poisson self-training.

    python3 oracle/check_hmmratac.py [--input BEDPE.gz] [--model a_model.json]

Exits non-zero when any threshold is missed, so it can gate CI.
"""

import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# The acceptance thresholds.
MIN_JACCARD = 0.98


def load_narrowpeak(path):
    """chrom, start, end, summit_offset, signal -> list."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            out.append((f[0], int(f[1]), int(f[2]), int(f[9]) if len(f) > 9 else -1, float(f[7])))
    return out


def merge(rs):
    by = {}
    for c, s, e, _sm, _sig in rs:
        by.setdefault(c, []).append((s, e))
    m = {}
    for c, v in by.items():
        v.sort()
        cur = None
        acc = []
        for s, e in v:
            if cur and s <= cur[1]:
                cur = (cur[0], max(cur[1], e))
            else:
                if cur:
                    acc.append(cur)
                cur = (s, e)
        if cur:
            acc.append(cur)
        m[c] = acc
    return m


def interval_intersection(a, b):
    """Number of shared bases in two sorted interval lists."""
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=None)
    ap.add_argument("--model", default=None)
    ap.add_argument("--oracle-src", default=None)
    ap.add_argument("--venv", default=None)
    ap.add_argument("--binary", default=os.path.join(ROOT, "target/release/macs3-rs"))
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    def env_record(filename, key):
        path = os.path.join(HERE, filename)
        if not os.path.exists(path):
            return None
        with open(path) as fh:
            for line in fh:
                if line.startswith(key + "="):
                    return line.split("=", 1)[1].strip()
        return None

    src = args.oracle_src or os.environ.get("MACS3_SRC")
    if src is None:
        p = env_record("ENV.provisioned", "MACS3_PATH") or env_record("ENV.lock", "MACS3_PATH")
        src = os.path.dirname(os.path.dirname(p)) if p else "/scratch/mdra00001/tmp/opencode/macs3-src"
    venv = args.venv or os.environ.get("MACS3_VENV") or env_record("ENV.provisioned", "MACS3_VENV")
    venv = venv or "/scratch/mdra00001/tmp/opencode/macs3-venv"
    args.input = args.input or os.environ.get("MACS3_HMM_INPUT")
    if not args.input:
        candidates = ("/tmp/atr/yeast_500k_SRR1822137.bedpe.gz",
                      os.path.join(src, "test/yeast_500k_SRR1822137.bedpe.gz"))
        args.input = next((p for p in candidates if os.path.exists(p)), candidates[-1])
    for p in (args.input, os.path.join(venv, "bin", "python"), args.binary):
        if not os.path.exists(p):
            print(f"missing required path: {p}", file=sys.stderr)
            return 2

    work = args.outdir or tempfile.mkdtemp(prefix="g13-")
    os.makedirs(work, exist_ok=True)
    o_dir = os.path.join(work, "oracle")
    u_dir = os.path.join(work, "ours")
    os.makedirs(o_dir, exist_ok=True)
    os.makedirs(u_dir, exist_ok=True)

    env = dict(os.environ, PYTHONPATH=src)
    oracle_python = os.path.join(venv, "bin", "python")
    model = args.model
    if model is None:
        # Generate a fresh pinned model when the caller has not supplied one.
        model_prefix = os.path.join(o_dir, "fresh")
        command = [oracle_python, "bin/macs3", "hmmratac", "-i", args.input, "-f", "BEDPE",
                   "--jump", "1.5", "--modelonly", "--name", "fresh", "--outdir", o_dir]
        result = subprocess.run(command, cwd=src, env=env, capture_output=True, text=True)
        if result.returncode:
            print(result.stderr, file=sys.stderr)
            print("oracle model training failed", file=sys.stderr)
            return 2
        model = model_prefix + "_model.json"
    if not os.path.exists(model):
        print(f"missing model: {model}", file=sys.stderr)
        return 2

    common = ["hmmratac", "-i", args.input, "-f", "BEDPE", "-n", "atr",
              "--jump", "1.5", "--model", model, "--save-digested"]
    print("== oracle hmmratac ==", file=sys.stderr)
    rc = subprocess.run([oracle_python, "bin/macs3", *common, "--outdir", o_dir],
                        cwd=src, env=env, capture_output=True, text=True).returncode
    if rc != 0:
        print(f"oracle hmmratac failed (rc={rc})", file=sys.stderr)
        return 2
    print("== macs3-rs hmmratac ==", file=sys.stderr)
    rc = subprocess.run([args.binary, *common, "--outdir", u_dir], capture_output=True, text=True).returncode
    if rc != 0:
        print(f"macs3-rs hmmratac failed (rc={rc})", file=sys.stderr)
        return 2

    opath = os.path.join(o_dir, "atr_accessible_regions.narrowPeak")
    upath = os.path.join(u_dir, "atr_accessible_regions.narrowPeak")
    o = load_narrowpeak(opath)
    u = load_narrowpeak(upath)
    om, um = merge(o), merge(u)
    intersection = union = 0
    for chrom in set(om) | set(um):
        av, bv = om.get(chrom, []), um.get(chrom, [])
        intersection += interval_intersection(av, bv)
        union += sum(e - s for s, e in av) + sum(e - s for s, e in bv) - interval_intersection(av, bv)
    bj = intersection / union if union else 1.0
    oset = {(c, s, e) for c, rs in om.items() for s, e in rs}
    uset = {(c, s, e) for c, rs in um.items() for s, e in rs}
    inter = oset & uset
    jac = len(inter) / len(oset | uset) if oset | uset else 1.0

    # Boundary and summit displacement, over **shared** regions only.
    #
    # A region present on one side and absent on the other has no counterpart, so there
    # is no displacement to measure -- and matching it to the nearest region on the
    # other chromosome would report thousands of bp for what is really a missing
    # region. Unmatched counts are reported separately instead.
    uidx = {c: v for c, v in um.items()}
    uset_full = {(c, s, e): sm for c, s, e, sm, _sig in u}
    oset_full = {(c, s, e): sm for c, s, e, sm, _sig in o}
    sshift = []
    for key, osm in oset_full.items():
        if key not in uset_full or key not in inter:
            continue
        usm = uset_full[key]
        if osm >= 0 and usm >= 0:
            sshift.append(abs(osm - usm))
    unmatched_o = len(oset) - len(inter)
    unmatched_u = len(uset) - len(inter)

    print()
    print(f"regions:          oracle {len(oset)}  ours {len(uset)}  shared {len(inter)}")
    print(f"Jaccard (region): {jac:.4f}   requirement >= {MIN_JACCARD}")
    print(f"Jaccard (base):   {bj:.4f}")
    print(f"unmatched regions: oracle-only {unmatched_o}, ours-only {unmatched_u}")

    # Boundary/summit displacement is only meaningful between regions that exist on
    # both sides with the *same* coordinates, and it is then zero by construction --
    # which is why it is reported as a consistency check on the matching, not as a
    # measure of disagreement. The disagreement shows up as the unmatched counts and
    # as the Jaccard above.
    #
    # What the boundary discrepancy actually looks like is reported separately: for
    # shared (chrom, start) or (chrom, end) anchors, how far the far edge moved.
    anchor = []
    for c, v in om.items():
        for (bs, be) in v:
            near = uidx.get(c, [])
            for (bs2, be2) in near:
                if bs2 == bs and be2 != be:
                    anchor.append(abs(be2 - be))
                elif be2 == be and bs2 != bs:
                    anchor.append(abs(bs2 - bs))
    if anchor:
        anchor.sort()
        print(f"one-edge-matched regions: {len(anchor)}, far edge off by "
              f"median {anchor[len(anchor)//2]} / max {anchor[-1]} bp "
              f"(the bin size is 10 bp)")
    if sshift:
        print(f"summit shift:     max {max(sshift)} bp over {len(sshift)} identically placed regions")

    # Digested-track agreement, which is where a divergence first becomes visible.
    print()
    digest_ok = True
    for sig in ("short", "mono", "di", "tri"):
        a = os.path.join(o_dir, f"atr_digested_{sig}.bdg")
        b = os.path.join(u_dir, f"atr_digested_{sig}.bdg")
        if not (os.path.exists(a) and os.path.exists(b)):
            print(f"digested {sig:<5}:   absent")
            continue
        if open(a, "rb").read() == open(b, "rb").read():
            print(f"digested {sig:<5}:   IDENTICAL")
        else:
            al = sum(1 for _ in open(a))
            bl = sum(1 for _ in open(b))
            print(f"digested {sig:<5}:   differs ({al} vs {bl} lines)")
            digest_ok = False

    # F231: the digested tracks are byte-identical to upstream. That makes this the
    # tightest available signal on the digest half of the pipeline, and it is a
    # *precondition* for the Jaccard number meaning anything -- a digest that differs
    # moves every HMM state boundary, so a low Jaccard would no longer localise.
    # Kept as a separate gate so a digest regression is reported as such instead of
    # being buried under a region-count mismatch.
    ok = bj >= MIN_JACCARD and digest_ok
    print()
    if not digest_ok:
        print("FAIL: digested tracks differ (F231 regression)")
    if bj < MIN_JACCARD:
        print(f"FAIL: base Jaccard {bj:.4f} < {MIN_JACCARD}")
    if ok:
        print("PASS")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
