#!/usr/bin/env python3
"""Benchmark macs3-rs against pinned MACS3 on upstream's own real data.

This is the source of the README's performance table. Inputs are the files MACS3
ships for its speed test (`test/CTCF_12878_5M.bed.gz`, 5.0 M ChIP reads, and
`test/Input_12878_5M.bed.gz`, 5.0 M control reads) plus the chr22 paired-end BAMs.
The bedGraph commands read the treatment pileup, control lambda and p-score tracks
that MACS3 itself writes for those reads, so both sides parse identical input.

Each workload runs `--reps` times per implementation, alternating sides. Wall clock
and peak RSS come from `/usr/bin/time -v` (the kernel's high-water mark); the
median of each is reported. On the first repetition every output file is also
byte-compared between the two sides, ignoring only the `# Command line:` header
of `*.xls`, which names the executable and output directory.

Usage:
    scripts/bench_real.py --macs3 .oracle/venv/bin/macs3 --src .oracle/macs3-src \\
        --work /tmp/macs3-rs-bench [--reps 3] [--json out.json] [--only NAME]
"""
from __future__ import annotations

import argparse
import filecmp
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def workloads(t: Path, prep: Path):
    chip, ctrl = t / "CTCF_12878_5M.bed.gz", t / "Input_12878_5M.bed.gz"
    pe_chip, pe_ctrl = t / "CTCF_PE_ChIP_chr22_50k.bam", t / "CTCF_PE_CTRL_chr22_50k.bam"
    treat, lam, ppois = prep / "prep_treat_pileup.bdg", prep / "prep_control_lambda.bdg", prep / "prep_ppois.bdg"
    se = ["callpeak", "-t", chip, "-g", "hs", "--nomodel", "--extsize", "200", "-n", "x"]
    return [
        ("callpeak, single-end -B", se + ["-c", ctrl, "-B", "--outdir", "{out}"]),
        ("callpeak, single-end --SPMR", se + ["-c", ctrl, "--SPMR", "--outdir", "{out}"]),
        ("callpeak, single-end --broad", se + ["-c", ctrl, "--broad", "--outdir", "{out}"]),
        ("callpeak, no control", se + ["--outdir", "{out}"]),
        ("callpeak, paired-end BAM", ["callpeak", "-f", "BAMPE", "-t", pe_chip, "-c", pe_ctrl,
                                      "-g", "52000000", "-n", "x", "--outdir", "{out}"]),
        ("pileup", ["pileup", "-f", "BED", "-i", chip, "--extsize", "200", "-o", "x.bdg", "--outdir", "{out}"]),
        ("filterdup", ["filterdup", "-f", "BED", "-i", chip, "-g", "hs", "-o", "x.bed", "--outdir", "{out}"]),
        ("randsample", ["randsample", "-f", "BED", "-i", chip, "-n", "2500000", "--seed", "42",
                        "-o", "x.bed", "--outdir", "{out}"]),
        ("bdgpeakcall", ["bdgpeakcall", "-i", ppois, "-c", "5", "-o", "x.narrowPeak", "--outdir", "{out}"]),
        ("bdgopt -m p2q", ["bdgopt", "-i", ppois, "-m", "p2q", "-o", "x.bdg", "--outdir", "{out}"]),
        ("cmbreps -m max", ["cmbreps", "-i", treat, lam, "-m", "max", "-o", "x.bdg", "--outdir", "{out}"]),
        ("bdgcmp -m ppois", ["bdgcmp", "-t", treat, "-c", lam, "-m", "ppois", "-o", "x.bdg", "--outdir", "{out}"]),
    ]


def timed(binary: str, argv: list, out: Path) -> tuple[float, float]:
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    cmd = ["/usr/bin/time", "-v", binary] + [str(a).replace("{out}", str(out)) for a in argv]
    p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.exit(f"FAILED ({p.returncode}): {' '.join(cmd)}\n{p.stderr[-2000:]}")
    m = re.search(r"Elapsed \(wall clock\) time.*?: (?:(\d+):)?(\d+):([\d.]+)", p.stderr)
    secs = int(m.group(1) or 0) * 3600 + int(m.group(2)) * 60 + float(m.group(3))
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", p.stderr).group(1)) / 1024
    return secs, rss


def same(a: Path, b: Path) -> bool:
    if a.suffix != ".xls":
        return filecmp.cmp(a, b, shallow=False)
    strip = lambda p: [l for l in p.read_text().splitlines() if not l.startswith("# Command line:")]
    return strip(a) == strip(b)


def compare(up: Path, rs: Path) -> tuple[int, int]:
    names = sorted(p.name for p in up.iterdir())
    ok = sum(1 for n in names if (rs / n).exists() and same(up / n, rs / n))
    return ok, len(names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--macs3", required=True)
    ap.add_argument("--src", required=True, help="pinned MACS3 source tree")
    ap.add_argument("--rs", default=str(ROOT / "target/release/macs3-rs"))
    ap.add_argument("--work", required=True)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--json")
    ap.add_argument("--only")
    a = ap.parse_args()
    os.environ.setdefault("OPENBLAS_CORETYPE", "Haswell")
    work, t = Path(a.work), Path(a.src).resolve() / "test"
    macs3, rs = str(Path(a.macs3).absolute()), str(Path(a.rs).absolute())
    prep = work / "prep"
    if not (prep / "prep_ppois.bdg").exists():
        prep.mkdir(parents=True, exist_ok=True)
        subprocess.run([macs3, "callpeak", "-t", t / "CTCF_12878_5M.bed.gz", "-c", t / "Input_12878_5M.bed.gz",
                        "-g", "hs", "--nomodel", "--extsize", "200", "-B", "-n", "prep", "--outdir", prep],
                       check=True, capture_output=True)
        subprocess.run([macs3, "bdgcmp", "-t", prep / "prep_treat_pileup.bdg", "-c", prep / "prep_control_lambda.bdg",
                        "-m", "ppois", "-o", "prep_ppois.bdg", "--outdir", prep], check=True, capture_output=True)
    rows = []
    for name, argv in workloads(t, prep):
        if a.only and a.only not in name:
            continue
        slug = re.sub(r"\W+", "_", name)
        up_t, up_m, rs_t, rs_m, cmp_ = [], [], [], [], None
        for rep in range(a.reps):
            s, m = timed(macs3, argv, work / slug / "up"); up_t.append(s); up_m.append(m)
            s, m = timed(rs, argv, work / slug / "rs"); rs_t.append(s); rs_m.append(m)
            if rep == 0:
                cmp_ = compare(work / slug / "up", work / slug / "rs")
        med = statistics.median
        row = dict(name=name, macs3_s=med(up_t), macs3_mb=med(up_m), rs_s=med(rs_t), rs_mb=med(rs_m),
                   identical=cmp_[0], files=cmp_[1], runs=dict(up_t=up_t, up_m=up_m, rs_t=rs_t, rs_m=rs_m))
        rows.append(row)
        print(f"{name:30s} macs3 {row['macs3_s']:6.2f} s {row['macs3_mb']:5.0f} MB | rs {row['rs_s']:6.2f} s "
              f"{row['rs_mb']:5.0f} MB | {row['macs3_s'] / row['rs_s']:.2f}x {100 * row['rs_mb'] / row['macs3_mb']:.0f}% "
              f"| files identical {cmp_[0]}/{cmp_[1]}", flush=True)
    if a.json:
        Path(a.json).write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
