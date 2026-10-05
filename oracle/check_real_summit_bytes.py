#!/usr/bin/env python3
"""Compare complete summit-related callpeak outputs with the pinned oracle.

Runs the source-tree CTCF chr22 single-end, BEDPE, and BAMPE examples through
both MACS3 and macs3-rs. Every generated summit-related file is compared byte
for byte. The only normalization is the XLS ``# Command line:`` header, which
necessarily names different executables and scratch output directories.

Requires ``MACS3_SRC`` and ``MACS3_ORACLE_BIN`` (or ``--source`` and
``--oracle``). The Rust binary defaults to ``target/release/macs3-rs``.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
COMMON_OUTPUTS = (
    "target_control_lambda.bdg",
    "target_peaks.narrowPeak",
    "target_peaks.xls",
    "target_summits.bed",
    "target_treat_pileup.bdg",
)
OUTPUTS = {
    # Single-end callpeak writes the fitted fragment model in addition to the
    # peak/summit outputs. Keep it in the same fresh differential gate.
    "SE": (*COMMON_OUTPUTS, "target_model.r"),
    "BEDPE": COMMON_OUTPUTS,
    "BAMPE": COMMON_OUTPUTS,
}
CASES = (
    (
        "SE",
        "BED",
        "CTCF_SE_ChIP_chr22_50k.bed.gz",
        "CTCF_SE_CTRL_chr22_50k.bed.gz",
        ("--d-min", "15"),
    ),
    (
        "BEDPE",
        "BEDPE",
        "CTCF_PE_ChIP_chr22_50k.bedpe.gz",
        "CTCF_PE_CTRL_chr22_50k.bedpe.gz",
        (),
    ),
    (
        "BAMPE",
        "BAMPE",
        "CTCF_PE_ChIP_chr22_50k.bam",
        "CTCF_PE_CTRL_chr22_50k.bam",
        (),
    ),
)


def read_file(path: Path) -> bytes:
    return path.read_bytes()


def normalized_xls(path: Path) -> bytes:
    return b"\n".join(
        line for line in read_file(path).splitlines()
        if not line.startswith(b"# Command line:")
    )


def run(binary: str, case: tuple, source_test: Path, outdir: Path) -> None:
    _label, fmt, treatment, control, extra = case
    outdir.mkdir(parents=True)
    argv = [
        binary,
        "callpeak",
        "-t",
        str(source_test / treatment),
        "-c",
        str(source_test / control),
        "-f",
        fmt,
        "-g",
        "52000000",
        "-B",
        "--call-summits",
        *extra,
        "--outdir",
        str(outdir),
        "-n",
        "target",
    ]
    proc = subprocess.run(argv, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(
            f"{binary} failed for {case[0]} (exit {proc.returncode}):\n"
            f"{proc.stderr[-4000:]}\n{proc.stdout[-1000:]}"
        )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", default=os.environ.get("MACS3_SRC"))
    ap.add_argument("--oracle", default=os.environ.get("MACS3_ORACLE_BIN"))
    ap.add_argument("--rust", default=str(ROOT / "target/release/macs3-rs"))
    args = ap.parse_args()
    if not args.source or not args.oracle:
        ap.error("set MACS3_SRC and MACS3_ORACLE_BIN or pass --source and --oracle")
    source_test = Path(args.source).resolve() / "test"
    if not source_test.is_dir():
        ap.error(f"source test directory not found: {source_test}")
    for case in CASES:
        for name in (case[2], case[3]):
            if not (source_test / name).is_file():
                ap.error(f"missing fixture: {source_test / name}")
    for binary in (args.oracle, args.rust):
        if not shutil.which(binary) and not Path(binary).is_file():
            ap.error(f"executable not found: {binary}")

    failures = 0
    with tempfile.TemporaryDirectory(prefix="macs3-summit-bytes-") as tmp:
        base = Path(tmp)
        for case in CASES:
            label = case[0]
            oracle_out = base / label / "oracle"
            rust_out = base / label / "rust"
            run(args.oracle, case, source_test, oracle_out)
            run(args.rust, case, source_test, rust_out)
            missing = [
                f"{side}/{name}"
                for side, directory in (("oracle", oracle_out), ("rust", rust_out))
                for name in OUTPUTS[label]
                if not (directory / name).is_file()
            ]
            if missing:
                print(f"{label}: missing output(s): {', '.join(missing)}")
                failures += 1
                continue
            mismatches = []
            for name in OUTPUTS[label]:
                a, b = oracle_out / name, rust_out / name
                left = normalized_xls(a) if name.endswith(".xls") else read_file(a)
                right = normalized_xls(b) if name.endswith(".xls") else read_file(b)
                if left != right:
                    mismatches.append(name)
            if mismatches:
                print(f"{label}: byte mismatches: {', '.join(mismatches)}")
                failures += 1
            else:
                print(
                    f"{label}: all {len(OUTPUTS[label])} summit-related outputs match byte for byte"
                )

    if failures:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
