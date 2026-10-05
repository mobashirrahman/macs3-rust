#!/usr/bin/env bash
# Verify the pinned checkout is pristine and its loaded extension binaries match
# the hashes captured immediately after provisioning.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCK="$HERE/ENV.lock"
PROVISIONED="$HERE/ENV.provisioned"
fail=0

env_value() {
  local key="$1" file="$2"
  awk -F= -v key="$key" '$1 == key || $1 == "# " key { sub(/^[^=]*=/, ""); print; exit }' "$file" 2>/dev/null
}

SRC="${MACS3_SRC:-}"
if [ -z "$SRC" ]; then
  SRC="$(env_value MACS3_SRC "$PROVISIONED")"
fi
if [ -z "$SRC" ]; then
  recorded="$(env_value MACS3_PATH "$LOCK")"
  SRC="${recorded%/MACS3/__init__.py}"
fi
if [ -z "$SRC" ] || [ ! -d "$SRC/.git" ]; then
  echo "FAIL oracle source not found; set MACS3_SRC or run oracle/provision_oracle.sh" >&2
  exit 1
fi

want="$(env_value MACS3_COMMIT "$LOCK")"
have="$(git -C "$SRC" rev-parse HEAD 2>/dev/null)"
if [ "$want" != "$have" ]; then
  echo "FAIL commit: want $want, have $have"
  fail=1
else
  echo "ok   commit $have"
fi

dirty="$(git -C "$SRC" status --porcelain)"
if [ -n "$dirty" ]; then
  echo "FAIL oracle tracked/untracked tree is modified:"
  echo "$dirty" | sed 's/^/     /'
  fail=1
else
  echo "ok   oracle tree clean"
fi

if grep -rlq --include='*.py' 'MACS3RS_' "$SRC/MACS3" 2>/dev/null; then
  echo "FAIL instrumentation hooks present in Python sources:"
  grep -rl --include='*.py' 'MACS3RS_' "$SRC/MACS3" | sed 's/^/     /'
  fail=1
else
  echo "ok   no instrumentation hooks in Python sources"
fi

MANIFEST="${MACS3_BUILD_MANIFEST:-$(env_value MACS3_BUILD_MANIFEST "$PROVISIONED")}"
if [ -n "$MANIFEST" ] && [ -f "$MANIFEST" ]; then
  manifest_commit="$(env_value MACS3_COMMIT "$MANIFEST")"
  if [ "$manifest_commit" != "$want" ]; then
    echo "FAIL extension manifest commit: $manifest_commit != $want"
    fail=1
  fi
  for key in MACS3_VERSION NUMPY SCIPY SKLEARN HMMLEARN CYTHON_BUILD OPENBLAS_CORETYPE; do
    locked="$(env_value "$key" "$LOCK")"
    recorded="$(env_value "$key" "$MANIFEST")"
    if [ -n "$locked" ] && [ "$locked" != "$recorded" ]; then
      echo "FAIL extension manifest $key: $recorded != lock $locked"
      fail=1
    fi
  done
  lock_python="$(env_value PYTHON "$LOCK" | cut -d. -f1,2)"
  manifest_python="$(env_value PYTHON "$MANIFEST" | cut -d. -f1,2)"
  if [ -n "$lock_python" ] && [ "$lock_python" != "$manifest_python" ]; then
    echo "FAIL extension manifest Python: $manifest_python != lock $lock_python"
    fail=1
  fi
  expected_paths="$(awk '!/^#/ && NF >= 2 {print $2}' "$MANIFEST" | sort)"
  actual_paths="$(cd "$SRC" && find MACS3 -type f -name '*.so' -printf '%p\n' | sort)"
  if [ "$expected_paths" != "$actual_paths" ]; then
    echo "FAIL compiled extension set differs from pristine manifest"
    diff -u <(printf '%s\n' "$expected_paths") <(printf '%s\n' "$actual_paths") || true
    fail=1
  else
    while read -r hash relpath; do
      [[ "$hash" == \#* || -z "$hash" ]] && continue
      actual="$(sha256sum "$SRC/$relpath" | awk '{print $1}')"
      if [ "$actual" != "$hash" ]; then
        echo "FAIL extension hash: $relpath"
        fail=1
      fi
    done < "$MANIFEST"
    [ "$fail" -ne 0 ] || echo "ok   compiled extensions match pristine build manifest"
  fi
else
  echo "note no generated build manifest; using any available local pre-instrumentation backups"
fi

# Preserve the developer's older pre-instrumentation comparison where backups
# exist, while allowing CI checkouts with no such /tmp files to use the manifest.
BACKUP_DIR="${ORACLE_SO_BACKUP_DIR:-/tmp/oracle_so_backup}"
backup_matches=0
for module in CallPeakUnit Pileup PileupV2; do
  so="$(find "$SRC/MACS3/Signal" -maxdepth 1 -type f -name "$module.*.so" -print -quit 2>/dev/null)"
  if [ -n "$so" ] && [ -f "$BACKUP_DIR/$(basename "$so")" ]; then
    if cmp -s "$so" "$BACKUP_DIR/$(basename "$so")"; then
      echo "ok   $module.so matches local pre-instrumentation backup"
      backup_matches=$((backup_matches + 1))
    else
      echo "FAIL $module.so differs from local pre-instrumentation backup"
      fail=1
    fi
  fi
done

if [ ! -f "${MANIFEST:-}" ] && [ "$backup_matches" -eq 0 ]; then
  echo "FAIL no extension build manifest or local pre-instrumentation backup is available"
  fail=1
fi

exit "$fail"
