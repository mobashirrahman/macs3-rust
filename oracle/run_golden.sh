#!/bin/bash
# Golden byte-parity gate: replay every recorded `macs3` command with `macs3-rs`
# and compare the output files byte-for-byte.
#
# This is the acceptance criterion "byte-identical *.xls, *_peaks.narrowPeak,
# *_summits.bed, *_model.r across the whole corpus and flag matrix". The existing
# `run_peak_e2e.sh` compares *coordinates*, which cannot see a score that moved in
# the fourth decimal; `run_golden.py` compares the files.
#
# Usage:
#   oracle/run_golden.sh                       # every recorded case
#   oracle/run_golden.sh --variant default     # one flag-matrix variant
#   oracle/run_golden.sh --fixture sweep       # fixtures matching a substring
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIN="$ROOT/target/release/macs3-rs"
if [ ! -x "$BIN" ]; then
    echo "building macs3-rs (release)..."
    cargo build --release --bin macs3-rs || exit 1
fi

python3 "$ROOT/oracle/run_golden.py" --bin "$BIN" "$@"