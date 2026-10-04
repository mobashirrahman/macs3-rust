#!/usr/bin/env bash
# Provision the pinned MACS3 3.0.5 oracle for the L3/L4 differential layers.
#
# The oracle is intentionally *not* vendored: it is a git checkout plus a Python
# environment, and the release definition requires the shipped tool to be free of
# any Python dependency. Vendoring it would also make it easy to "fix" the oracle,
# which would silently invalidate every recorded golden. So CI provisions it from
# the pinned commit recorded in `oracle/ENV.lock` and asserts the checkout is
# byte-clean afterwards (`oracle/verify_oracle_clean.sh`).
#
# Usage:
#   bash oracle/provision_oracle.sh [--dest DIR] [--check-only]
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCK="$REPO_ROOT/oracle/ENV.lock"
DEST="${MACS3_ORACLE_DIR:-$REPO_ROOT/.oracle/macs3-src}"
CHECK_ONLY=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dest) DEST="$2"; shift 2 ;;
    --check-only) CHECK_ONLY=1; shift ;;
    -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
    *) echo "provision_oracle: unknown argument $1" >&2; exit 2 ;;
  esac
done

[ -f "$LOCK" ] || { echo "provision_oracle: missing $LOCK" >&2; exit 1; }
# shellcheck disable=SC1090
COMMIT="$(awk -F= '/^MACS3_COMMIT=/{print $2}' "$LOCK")"
VERSION="$(awk -F= '/^MACS3_VERSION=/{print $2}' "$LOCK")"
URL="https://github.com/macs3-project/MACS3.git"

echo "pinning MACS3 $VERSION at $COMMIT"

if [ "$CHECK_ONLY" = "1" ]; then
  # Fall back to wherever the recorded lock file says the oracle already lives, so a
  # developer who provisioned it by hand can validate the pin without re-cloning.
  if [ ! -d "$DEST/.git" ]; then
    recorded="$(awk -F= '/^MACS3_PATH=/{print $2}' "$LOCK")"
    recorded="${recorded%/MACS3/__init__.py}"
    [ -d "$recorded/.git" ] && DEST="$recorded"
  fi
  [ -d "$DEST/.git" ] || { echo "not provisioned: $DEST" >&2; exit 1; }
  got="$(git -C "$DEST" rev-parse HEAD)"
  [ "$got" = "$COMMIT" ] || { echo "wrong commit: $got != $COMMIT" >&2; exit 1; }
  echo "oracle already provisioned at $DEST ($got)"
  exit 0
fi

if [ ! -d "$DEST/.git" ]; then
  mkdir -p "$(dirname "$DEST")"
  # `--no-checkout` keeps the working tree clean of the pinned branch name, so the
  # checkout below is the only thing that ever touches it
  git clone --quiet --filter=blob:none --no-checkout "$URL" "$DEST"
fi

git -C "$DEST" fetch --quiet --depth 1 origin "$COMMIT" 2>/dev/null \
  || git -C "$DEST" fetch --quiet origin
git -C "$DEST" checkout --quiet --force "$COMMIT"

got="$(git -C "$DEST" rev-parse HEAD)"
if [ "$got" != "$COMMIT" ]; then
  echo "provision_oracle: checked out $got, expected $COMMIT" >&2
  exit 1
fi

PY="${PYTHON:-python3}"
"$PY" -m pip install --quiet --disable-pip-version-check \
  -r "$DEST/requirements.txt" 2>/dev/null \
  || "$PY" -m pip install --quiet --disable-pip-version-check cython numpy scipy

MACS3_SRC="$DEST" bash "$REPO_ROOT/oracle/verify_oracle_clean.sh"

# Record what we actually got, so a drifting dependency is visible in the log rather
# than inferred later from a byte mismatch.
{
  echo "MACS3_VERSION=$VERSION"
  echo "MACS3_COMMIT=$COMMIT"
  echo "MACS3_PATH=$DEST/MACS3/__init__.py"
  echo "PYTHON=$("$PY" -c 'import platform;print(platform.python_version())')"
  echo "NUMPY=$("$PY" -c 'import numpy;print(numpy.__version__)')"
  echo "PROVISIONED_BY=oracle/provision_oracle.sh"
} > "$REPO_ROOT/oracle/ENV.provisioned"

echo "oracle ready: $DEST"
echo "export MACS3_ORACLE_DIR=$DEST"
