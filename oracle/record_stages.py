#!/usr/bin/env python3
"""Record upstream's intermediate stages for a set of fixtures.

The porting plan's method step 1 requires the oracle to freeze "final outputs AND
intermediates": reads, duplicates, `d`, model arrays, treatment pileup, control
pileup, per-scale lambda, merged lambda, p-score track, q-score table,
candidate/merged intervals, summits.

`dump_stages.py` can observe all of those, but it only ever ran ad hoc and nothing
was kept, so there was no recorded stage corpus to diff against -- which is exactly
what the stage-by-stage test layer (L3) needs. This script drives it over a fixture
list and commits the result.

`dump_stages.py` drives upstream's real functions in-process, so what lands here is
upstream's arithmetic and ordering; only the observation is ours. Nothing is
patched into the oracle, so `verify_oracle_clean.sh` still passes.

Usage:
    python3 oracle/record_stages.py --out tests/stages [--limit N]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES = os.path.join(REPO, "tests", "fixtures")

# (fixture group, fixture name, mode)
SELECTION = [
    ("se_basic", "gauss_two_peaks", "se"),
    ("pe_basic", "gauss_fragments", "pe"),
    ("pe_basic", "atac_short", "pe"),
    ("pe_basic", "nucleosome_ladder", "pe"),
    ("frag_basic", "barcode_fragments", "frag"),
]

# The extra `macs3` arguments each mode needs. `--gsize` is mandatory upstream.
ARGS = {
    "se": ["--nomodel", "--extsize", "200"],
    "pe": [],
    "frag": [],
}


def manifest_mode(fixture_dir: str) -> str | None:
    path = os.path.join(fixture_dir, "manifest.tsv")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] == "mode":
                return parts[1]
    return None


def genome_size(fixture_dir: str) -> str:
    """The numeric effective genome size for this fixture.

    `genome.txt` is a *table* of per-chromosome lengths with a `# sum=` comment, not a
    value upstream can resolve -- passing it verbatim makes `callpeak` reject the run,
    and `dump_stages.py` then records a stage dump full of empty payloads with
    `_systexit: "1"`. The sum is what the golden corpus uses, and
    `manifest.tsv`'s `effective_gsize` column agrees.
    """
    p = os.path.join(fixture_dir, "genome.txt")
    total = 0
    if os.path.exists(p):
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    total += int(parts[1])
    return str(total) if total else "2000000"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(REPO, "tests", "stages"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    env = dict(os.environ)
    src = os.environ.get("MACS3_SRC", "/scratch/mdra00001/tmp/opencode/macs3-src")
    env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")

    made, skipped = [], []
    for group, name, mode in SELECTION:
        fixture_dir = os.path.join(FIXTURES, group, name)
        if not os.path.isdir(fixture_dir):
            skipped.append(f"{group}/{name}: no such fixture")
            continue
        declared = manifest_mode(fixture_dir)
        if declared and declared != mode:
            skipped.append(f"{group}/{name}: manifest says mode={declared}, asked {mode}")
            continue
        out = os.path.join(args.out, group, name, mode)
        if args.dry_run:
            print(f"would record {group}/{name} [{mode}] -> {out}")
            made.append(out)
            continue
        os.makedirs(out, exist_ok=True)
        # upstream's own log is the only source for the duplicate *rates*; the
        # retained counts are already visible as the pre/post filter totals.
        log_path = os.path.join(out, "oracle.log")
        cmd = [
            sys.executable,
            os.path.join(REPO, "oracle", "dump_stages.py"),
            fixture_dir,
            mode,
            "--out",
            out,
            "--log",
            log_path,
            "--",
            "--gsize",
            genome_size(fixture_dir),
            *ARGS[mode],
            "--cutoff-analysis",
        ]
        p = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=1800)
        with open(log_path, "w") as fh:
            fh.write(p.stderr)
        ok = False
        sj = os.path.join(out, "stages.json")
        if os.path.exists(sj):
            try:
                with open(sj) as fh:
                    blob = json.load(fh)
                # `dump_stages.py` writes the file even when the pipeline raised, and
                # records the exception under `_systexit`. A dump full of empty
                # payloads with a non-null `_systexit` is a failed run, not data.
                if str(blob.get("_systexit", "0")) != "0":
                    print(
                        f"  {group}/{name}: pipeline exited {blob['_systexit']}",
                        file=sys.stderr,
                    )
                else:
                    ok = True
            except (OSError, ValueError) as e:
                print(f"  {group}/{name}: unreadable stages.json: {e}", file=sys.stderr)
        (made if ok else skipped).append(
            f"{group}/{name} [{mode}]: {'ok' if ok else 'FAILED rc=%s' % p.returncode}"
        )
        if not ok:
            print(p.stderr[-600:], file=sys.stderr)

    for m in made:
        print(f"recorded  {m}")
    for s in skipped:
        print(f"skipped   {s}", file=sys.stderr)
    if args.limit:
        made = made[: args.limit]
    print(f"\n{len(made)} recorded, {len(skipped)} skipped")
    return 0 if made else 1


if __name__ == "__main__":
    raise SystemExit(main())
