#!/usr/bin/env bash
# G0 hermeticity check: the pinned oracle tree must be byte-identical to the
# upstream commit recorded in ENV.lock, with no local instrumentation.
#
# Option 1a of the porting plan instrumented MACS3 in place to capture
# intermediates. That is no longer needed -- `dump_stages.py` drives the real
# pipeline and records what it computed, so it can observe everything the
# instrumentation exposed without touching the oracle. This script guards
# against that instrumentation creeping back in.
set -uo pipefail
# resolve relative to this script so it works from any cwd (CI runs it from the
# workspace root, `cargo test` from a crate dir)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${MACS3_SRC:-/scratch/mdra00001/tmp/opencode/macs3-src}"
LOCK="$HERE/ENV.lock"
fail=0

want=$(grep '^MACS3_COMMIT=' "$LOCK" | cut -d= -f2)
have=$(git -C "$SRC" rev-parse HEAD 2>/dev/null)
if [ "$want" != "$have" ]; then
  echo "FAIL commit: want $want, have $have"; fail=1
else
  echo "ok   commit $have"
fi

# any tracked modification to the source tree
dirty=$(git -C "$SRC" status --porcelain)
if [ -n "$dirty" ]; then
  echo "FAIL oracle tree is modified:"; echo "$dirty" | sed 's/^/     /'; fail=1
else
  echo "ok   oracle tree clean"
fi

# Any instrumentation hooks left behind. Only the `.py` sources matter -- the
# `.c` files beside them are Cython build artifacts, are gitignored, and are not
# imported at runtime. They are still stale copies of the instrumented build, so
# they are reported separately rather than counted as a failure; the loaded
# modules are the `.so` files, which are checked byte-for-byte below.
if grep -rlq --include='*.py' 'MACS3RS_' "$SRC/MACS3" 2>/dev/null; then
  echo "FAIL instrumentation hooks present in:"; grep -rl --include='*.py' 'MACS3RS_' "$SRC/MACS3" | sed 's/^/     /'; fail=1
else
  echo "ok   no instrumentation hooks in Python sources"
fi
stale=$(grep -rl 'MACS3RS_' "$SRC/MACS3" --include='*.c' 2>/dev/null | wc -l)
if [ "$stale" -gt 0 ]; then
  echo "note $stale stale Cython .c artifact(s) mention the hooks (gitignored, not loaded)"
fi

# the built extensions must match the backups taken before instrumentation
for m in CallPeakUnit Pileup PileupV2; do
  so="$SRC/MACS3/Signal/$m.cpython-312-x86_64-linux-gnu.so"
  bk="/tmp/oracle_so_backup/$m.cpython-312-x86_64-linux-gnu.so"
  if [ -f "$so" ] && [ -f "$bk" ] && cmp -s "$so" "$bk"; then
    echo "ok   $m.so matches pre-instrumentation build"
  else
    echo "FAIL $m.so does not match its backup"; fail=1
  fi
done

exit $fail
