#!/usr/bin/env python3
"""Require every documented upstream finding to reference Rust source or reviewed debt."""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FINDINGS = ROOT / "docs/upstream-findings.md"
DEBT = ROOT / "docs/finding-reference-debt.txt"


def finding_ids(text: str) -> set[str]:
    return set(re.findall(r"^## (F\d+)\b", text, re.MULTILINE))


def debt_ids(text: str) -> set[str]:
    return {line.strip() for line in text.splitlines()
            if re.fullmatch(r"F\d+", line.strip())}


def main() -> int:
    findings = finding_ids(FINDINGS.read_text())
    if not findings:
        print("no documented findings were parsed", file=sys.stderr)
        return 1
    debt = debt_ids(DEBT.read_text())
    source = "\n".join(path.read_text(errors="replace")
                      for path in (ROOT / "crates").rglob("*.rs"))
    referenced = findings & set(re.findall(r"\bF\d+\b", source))
    missing = findings - referenced
    stale = debt - missing
    unreviewed = missing - debt
    print(f"{len(findings)} findings: {len(referenced)} referenced in Rust, "
          f"{len(missing)} reviewed historical debt")
    if stale:
        print("remove resolved IDs from docs/finding-reference-debt.txt: " +
              ", ".join(sorted(stale, key=lambda value: int(value[1:]))),
              file=sys.stderr)
    if unreviewed:
        print("new findings need a Rust reference or reviewed debt entry: " +
              ", ".join(sorted(unreviewed, key=lambda value: int(value[1:]))),
              file=sys.stderr)
    unknown = debt - findings
    if unknown:
        print("debt IDs are not documented findings: " +
              ", ".join(sorted(unknown, key=lambda value: int(value[1:]))),
              file=sys.stderr)
    return int(bool(stale or unreviewed or unknown))


if __name__ == "__main__":
    raise SystemExit(main())
