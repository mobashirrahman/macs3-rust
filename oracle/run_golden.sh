#!/bin/bash
# Golden byte-parity gate: replay every recorded `macs3` command with `macs3-rs`
# and compare the output files byte-for-byte.
#
# This is the acceptance criterion "byte-identical *.xls, *_peaks.narrowPeak,
# *_summits.bed, *_model.r across the whole corpus and flag matrix". The existing
# `run_peak_e2e.sh` compares *coordinates*, which cannot see a score that moved in
# the fourth decimal; `run_golden.py` compares the files.
#
# The corpus is committed with its paths rewritten to `<ROOT>` and friends, so the
# gate is the same set of numbers from any checkout -- a runner, a second clone.
# `relocate_golden.py --check` runs *after* the replay and is a second, independent
# condition on the same corpus: no recorded file names a checkout, and every recorded
# SHA-256 still describes the file beside it. It is a second gate rather than a
# preflight because the byte-parity verdict is this script's headline, and a corpus
# that has regressed should be reported as the failing cases it is.
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
rc=$?

echo
echo "== recorded corpus is path-independent =="
python3 "$ROOT/oracle/relocate_golden.py" --check || rc=1

exit "$rc"