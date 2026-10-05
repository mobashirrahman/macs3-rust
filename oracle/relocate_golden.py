#!/usr/bin/env python3
"""Rewrite the recorded corpora so they name no checkout.

`tests/golden/**` and `tests/stages/**` are *recorded data*: 7,685 recorded MACS3
invocations plus the recorded intermediate stages. Every `command.json` and every
`*_peaks.xls` header in them names the directory the corpus was recorded in --
`<root>/tests/fixtures/...`, `<root>/tests/golden/...`, and the interpreter it was
driven with. Committed that way, the corpus only replays from the one absolute path
the recording machine had: `oracle/run_golden.py` can neutralise *this* checkout's
root, so every other checkout -- a GitHub runner, a second clone -- reports 0/7204
cases byte-identical, and the cause is a path, not a difference in the output.

This script is the one-time rewrite that removes that coupling, and it is committed
so the rewrite is reproducible and auditable rather than a pile of hand edits:

    <recorded root>            -> <ROOT>        the checkout root
    <recorded interpreter>     -> macs3         argv[0]; never expanded again
    <recorded MACS3 checkout>  -> <MACS3_SRC>   the pinned upstream source tree
    <author scratchpad>        -> <TMP>         nothing else is expected here
    <author home>              -> <HOME>        last-resort catch-all

Every replacement is a pure byte substitution of a long, specific ASCII prefix; all
other bytes of every file are copied verbatim. `command.json`'s per-file SHA-256 map
is then recomputed for each file whose bytes moved, because a recorded digest is a
claim about those bytes and a stale one would be a lie the corpus cannot detect.

`run_golden.py` expands the token back to the checkout it is running in when it
replays a command and maps this checkout's root back to the token when it compares;
`run_oracle.py` and `record_stages.py` write the token when they record. This script
only had to run once; `--check` is what keeps it true.

Idempotent: the tokens contain none of the recorded prefixes, so a second run
rewrites nothing and reports zero.

Usage:
    oracle/relocate_golden.py            # rewrite (a no-op once already relocated)
    oracle/relocate_golden.py --check    # fail if a recorded path came back
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from run_golden import (  # noqa: E402
    HOME_TOKEN,
    MACS3_SRC_TOKEN,
    MACS3_TOKEN,
    ROOT_TOKEN,
    TMP_TOKEN,
    neutralise,
    substitution_table,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN = os.path.join(REPO, "tests", "golden")
STAGES = os.path.join(REPO, "tests", "stages")

# The absolute prefixes baked into the committed recordings. These are constants *of
# the corpus*, not of this checkout: the whole point is that the rewrite does not
# depend on where the repository happens to be.
AUTHOR_HOME = "/scratch/mdra00001"
AUTHOR_ROOT = AUTHOR_HOME + "/MACS3-rust"
AUTHOR_INTERPRETER = AUTHOR_HOME + "/tmp/opencode/macs3-venv/bin/macs3"
AUTHOR_SRC = AUTHOR_HOME + "/tmp/opencode/macs3-src"
AUTHOR_TMP = AUTHOR_HOME + "/tmp/opencode"

# Longest needle first, so `<MACS3_SRC>` wins over `<TMP>` and `<ROOT>` over `<HOME>`.
# `substitution_table` sorts on length regardless; this order only documents intent.
RECORDED_SUBSTITUTIONS = substitution_table(
    [
        (AUTHOR_INTERPRETER, MACS3_TOKEN),
        (AUTHOR_SRC, MACS3_SRC_TOKEN),
        (AUTHOR_ROOT, ROOT_TOKEN),
        (AUTHOR_TMP, TMP_TOKEN),
        (AUTHOR_HOME, HOME_TOKEN),
    ]
)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rewrite(path: str) -> bool:
    """Apply `RECORDED_SUBSTITUTIONS` to `path` in place; True when the bytes moved."""
    with open(path, "rb") as fh:
        old = fh.read()
    new = neutralise(old, RECORDED_SUBSTITUTIONS)
    if new == old:
        return False
    with open(path, "wb") as fh:
        fh.write(new)
    return True


def walk_files(root: str):
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in sorted(filenames):
            # `*.log` is git-ignored: a local run's log, never part of the corpus.
            if name.endswith(".log"):
                continue
            yield os.path.join(dirpath, name)


def case_dirs(root: str):
    """Every recorded run directory under `root`, as path components."""
    for dirpath, _dirnames, filenames in os.walk(root):
        if "command.json" in filenames:
            yield os.path.relpath(dirpath, root).split(os.sep)


def recorded_digests(case_dir: str) -> dict:
    with open(os.path.join(case_dir, "command.json")) as fh:
        return json.load(fh).get("files") or {}


def patch_digests(json_path: str, digests: dict[bytes, bytes]) -> int:
    """Write the recomputed digests back into a `command.json`/`summary.json`.

    Patched as a byte substitution of the quoted digest rather than by re-dumping the
    JSON, so the rewrite touches nothing but the digests that actually changed and the
    diff stays auditable. Two files with identical bytes share a digest before and
    after the rewrite -- same bytes in, same bytes out -- so replacing every occurrence
    of a stale digest is unambiguous.

    Returns the number of digests that were stale.
    """
    if not digests:
        return 0
    with open(json_path, "rb") as fh:
        blob = fh.read()
    for old, new in digests.items():
        blob = blob.replace(b'"' + old + b'"', b'"' + new + b'"')
    with open(json_path, "wb") as fh:
        fh.write(blob)
    return len(digests)


def relocate() -> int:
    """Rewrite both corpora, then recompute every digest the rewrite invalidated."""
    stale: dict[bytes, bytes] = {}
    cases = run_files = stage_files = 0

    # A run directory is rewritten as a unit: the outputs first, then the record that
    # describes them. `summary.json` is handled after the walk because its digests
    # repeat the ones the walk recomputed.
    for rel in case_dirs(GOLDEN):
        case_dir = os.path.join(GOLDEN, *rel)
        cases += 1
        recorded = recorded_digests(case_dir)
        digests: dict[bytes, bytes] = {}
        for path in walk_files(case_dir):
            run_files += int(rewrite(path))
            name = os.path.basename(path)
            if name in recorded:
                fresh = sha256(path)
                if fresh != recorded[name]:
                    digests[recorded[name].encode()] = fresh.encode()
        for old, new in digests.items():
            stale.setdefault(old, new)
        patch_digests(os.path.join(case_dir, "command.json"), digests)

    for path in walk_files(STAGES):
        stage_files += int(rewrite(path))

    summary = os.path.join(GOLDEN, "summary.json")
    summary_stale = patch_digests(summary, stale) if os.path.exists(summary) else 0

    print(f"recorded run directories : {cases}")
    print(f"golden files rewritten   : {run_files}")
    print(f"stage files rewritten    : {stage_files}")
    print(f"sha-256 digests recomputed: {len(stale)}"
          f"{f' (+{summary_stale} cached in summary.json)' if summary_stale else ''}")
    return 0


def check() -> int:
    """The postcondition, as a gate: no recorded prefix, no stale digest."""
    named = []
    files = 0
    live: set[str] = set()
    for root in (GOLDEN, STAGES):
        for path in walk_files(root):
            files += 1
            with open(path, "rb") as fh:
                blob = fh.read()
            if AUTHOR_HOME.encode() in blob:
                named.append(os.path.relpath(path, REPO))
    if named:
        print(f"{len(named)} file(s) still name a recorded path, e.g.:", file=sys.stderr)
        for p in named[:5]:
            print(f"  {p}", file=sys.stderr)
        return 1

    mismatched = 0
    for rel in case_dirs(GOLDEN):
        case_dir = os.path.join(GOLDEN, *rel)
        for name, want in recorded_digests(case_dir).items():
            target = os.path.join(case_dir, name)
            if not os.path.isfile(target):
                continue
            got = sha256(target)
            live.add(got)
            if got != want:
                mismatched += 1
    if mismatched:
        print(f"{mismatched} recorded digest(s) do not match their file", file=sys.stderr)
        return 1

    # `summary.json` is a cached rollup that repeats those digests. Nothing reads it
    # for a verdict, but a cache that disagrees with the files beside it is a trap.
    summary = os.path.join(GOLDEN, "summary.json")
    if os.path.exists(summary):
        with open(summary) as fh:
            for group in json.load(fh).get("fixtures", {}).values():
                for variant in group.get("variants", {}).values():
                    for digest in (variant.get("files") or {}).values():
                        if digest not in live:
                            mismatched += 1
        if mismatched:
            print(f"{mismatched} digest(s) do not match any recorded file", file=sys.stderr)
            return 1

    print(f"files scanned            : {files}")
    print(f"tokens                   : {ROOT_TOKEN} {MACS3_TOKEN} {MACS3_SRC_TOKEN} "
          f"{TMP_TOKEN} {HOME_TOKEN}")
    print("relocated corpus is clean")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="verify the postcondition instead of rewriting")
    args = ap.parse_args()
    return check() if args.check else relocate()


if __name__ == "__main__":
    sys.exit(main())