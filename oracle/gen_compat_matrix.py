#!/usr/bin/env python3
"""Generate the compatibility matrix from the recorded corpus.

The release definition requires a *CI-generated, never hand-edited* matrix. This
script is its only writer: it reads `tests/golden/**/command.json` -- the recorded
invocations, which are the authority -- replays each against the built `macs3-rs`,
and emits a per-(command, variant) table of byte-identical file counts.

Two design choices worth stating, because they are what keep the matrix honest:

* **It replays; it does not trust a summary.** `tests/golden/summary.json` is a
  cached rollup. If the matrix were built from it, a stale cache would produce a
  confident, wrong table.
* **Absence is reported as such.** A `(command, variant)` pair with zero recorded
  invocations is emitted with `0` and an explicit "not covered" marker, so a gap in
  the matrix reads as a gap instead of silently shrinking.
* **It replays in this checkout, from a recording that names none.** The recorded
  argv spells the checkout root `<ROOT>` (`oracle/relocate_golden.py`), and
  `run_golden.run_one` expands it, so the matrix is the same set of numbers in any
  directory. The counters below come from that replay; they are not a lookup.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
import run_golden

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.environ.get("MACS3_RS", os.path.join(REPO, "target/release/macs3-rs"))

HEADER = """<!-- GENERATED FILE -- do not edit.
Regenerate with:  python3 oracle/gen_compat_matrix.py --out docs/compatibility-matrix.md
Source of truth: tests/golden/**/command.json (the recorded invocations)
-->"""


def recorded():
    """Every recorded invocation, as (group, variant, argv, files)."""
    root = os.path.join(REPO, "tests", "golden")
    out = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if "command.json" not in filenames:
            continue
        with open(os.path.join(dirpath, "command.json")) as fh:
            rec = json.load(fh)
        argv = rec.get("command") or []
        if len(argv) < 2:
            continue
        sub = argv[1]
        if sub.startswith("-"):
            continue
        rel = os.path.relpath(dirpath, root)
        parts = rel.split(os.sep)
        variant = rec.get("variant") or (parts[-1] if len(parts) > 2 else "default")
        group = os.sep.join(parts[:-1]) or parts[0]
        out.append(
            {
                "subcommand": sub,
                "variant": variant,
                "fixture": rec.get("fixture", group),
                "argv": argv,
                "golden_dir": dirpath,
                "returncode": rec.get("returncode"),
                "fixture_path": group,
                "recording": rec,
            }
        )
    return out


def verify(rec):
    """Use the same complete byte and exit-status comparator as the golden gate."""
    _, same, total, exit_mismatch, notes = run_golden.run_one(
        BIN, rec["fixture_path"], rec["variant"], rec["recording"], None, False
    )
    missing = sum(note.startswith("missing ") for note in notes)
    return rec, same, total, missing, exit_mismatch


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(REPO, "docs/compatibility-matrix.md"))
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--list-only", action="store_true",
                    help="report coverage from the recordings without replaying")
    args = ap.parse_args()

    if not os.path.exists(BIN):
        print(f"gen_compat_matrix: {BIN} not built", file=sys.stderr)
        return 2

    recs = recorded()
    if not recs:
        print("gen_compat_matrix: no recorded invocations found", file=sys.stderr)
        return 2

    cells = {}
    exit_mismatches = 0
    if args.list_only:
        for r in recs:
            k = (r["subcommand"], r["variant"])
            cells.setdefault(k, [0, 0, 0, 0])
            cells[k][0] += 1
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for rec, same, compared, missing, exit_mismatch in pool.map(verify, recs):
                exit_mismatches += exit_mismatch
                k = (rec["subcommand"], rec["variant"])
                cell = cells.setdefault(k, [0, 0, 0, 0])
                cell[0] += 1
                cell[1] += same
                cell[2] += compared
                cell[3] += missing

    subs = sorted({k[0] for k in cells})
    variants = sorted({k[1] for k in cells})

    lines = [HEADER, "",
             "# Compatibility matrix", "",
             f"Generated from **{len(recs)} recorded invocations** by "
             "`oracle/gen_compat_matrix.py`. Each cell is "
             "`byte-identical files / recorded files` across the fixtures for that "
             "(subcommand, option variant). A dash means the pair is not covered by "
             "the recorded corpus, which is a gap in the matrix, not a pass.", "",
             "| subcommand | " + " | ".join(variants) + " |",
             "|---|" + "---|" * len(variants)]
    for s in subs:
        row = [f"`{s}`"]
        for v in variants:
            c = cells.get((s, v))
            if c is None:
                row.append("-- not covered --")
            elif args.list_only:
                row.append(f"{c[0]} invocations")
            else:
                row.append(f"{c[1]}/{c[2]}")
        lines.append("| " + " | ".join(row) + " |")

    totals = [0, 0, 0]
    for c in cells.values():
        totals[0] += c[1]
        totals[1] += c[2]
        totals[2] += c[3]
    lines += ["", f"**Totals: {totals[0]}/{totals[1]} compared output files "
                  f"byte-identical** across {len(cells)} (subcommand, variant) pairs. "
                  f"{totals[2]} recorded files were not produced by the replay "
                  f"and remain in the denominator. Exit-status mismatches: {exit_mismatches}."]
    if not args.list_only and totals[1] and totals[0] != totals[1]:
        lines += ["", "> Not all recorded files are byte-identical yet; see "
                      "`docs/upstream-findings.md` for the open findings."]
    lines.append("")

    with open(args.out, "w") as fh:
        fh.write("\n".join(lines))
    print(f"wrote {args.out}: {len(subs)} subcommands, {len(variants)} variants, "
          f"{totals[0]}/{totals[1]} files identical "
          f"({totals[2]} not produced by the replay)")
    return int(not args.list_only and (totals[0] != totals[1] or exit_mismatches > 0))


if __name__ == "__main__":
    sys.exit(main())
