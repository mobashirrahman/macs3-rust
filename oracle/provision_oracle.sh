#!/usr/bin/env bash
# Build a private, pinned MACS3 oracle for differential and golden tests.
#
# Usage: bash oracle/provision_oracle.sh [--dest DIR] [--venv DIR] [--check-only]
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCK="$REPO_ROOT/oracle/ENV.lock"
DEFAULT_ROOT="$REPO_ROOT/.oracle"
DEST="${MACS3_ORACLE_DIR:-$DEFAULT_ROOT/macs3-src}"
VENV="${MACS3_ORACLE_VENV:-$DEFAULT_ROOT/venv}"
CHECK_ONLY=0
DEST_EXPLICIT=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dest) DEST="$2"; DEST_EXPLICIT=1; shift 2 ;;
    --venv) VENV="$2"; shift 2 ;;
    --check-only) CHECK_ONLY=1; shift ;;
    -h|--help) sed -n '2,5p' "$0"; exit 0 ;;
    *) echo "provision_oracle: unknown argument $1" >&2; exit 2 ;;
  esac
done

[ -f "$LOCK" ] || { echo "provision_oracle: missing $LOCK" >&2; exit 1; }
lock_value() { awk -F= -v key="$1" '$1 == key { sub(/^[^=]*=/, ""); print; exit }' "$LOCK"; }
COMMIT="$(lock_value MACS3_COMMIT)"
VERSION="$(lock_value MACS3_VERSION)"
PYTHON_MINOR="$(lock_value PYTHON | cut -d. -f1,2)"
NUMPY="$(lock_value NUMPY)"
SCIPY="$(lock_value SCIPY)"
SKLEARN="$(lock_value SKLEARN)"
HMMLEARN="$(lock_value HMMLEARN)"
CYTHON_BUILD="$(lock_value CYTHON_BUILD)"
export OPENBLAS_CORETYPE="$(lock_value OPENBLAS_CORETYPE)"
OPENBLAS_CORETYPE="${OPENBLAS_CORETYPE:-Haswell}"
URL="https://github.com/macs3-project/MACS.git"

# Values are repeated as fallbacks for older local ENV.lock files. New locks
# record them explicitly so CI and local oracle runs use the same dependency set.
PYTHON_MINOR="${PYTHON_MINOR:-3.12}"
NUMPY="${NUMPY:-2.5.3}"
SCIPY="${SCIPY:-1.18.1}"
SKLEARN="${SKLEARN:-1.9.1}"
HMMLEARN="${HMMLEARN:-0.3.3}"
CYTHON_BUILD="${CYTHON_BUILD:-3.3.0}"
MANIFEST="$(dirname "$DEST")/oracle-build-manifest.sha256"
PROVISIONED="$REPO_ROOT/oracle/ENV.provisioned"

if [ -z "$COMMIT" ] || [ -z "$VERSION" ]; then
  echo "provision_oracle: ENV.lock must record MACS3_COMMIT and MACS3_VERSION" >&2
  exit 1
fi

if [ "$CHECK_ONLY" = "1" ]; then
  if [ "$DEST_EXPLICIT" = "0" ] && [ -n "${MACS3_SRC:-}" ] && [ ! -d "$DEST/.git" ]; then
    DEST="$MACS3_SRC"
  fi
  if [ ! -d "$DEST/.git" ]; then
    provisioned_src="$(awk -F= '$1 == "MACS3_SRC" {sub(/^[^=]*=/, ""); print; exit}' "$PROVISIONED" 2>/dev/null || true)"
    recorded="$(lock_value MACS3_PATH)"
    recorded="${recorded%/MACS3/__init__.py}"
    if [ -n "$provisioned_src" ] && [ -d "$provisioned_src/.git" ]; then
      DEST="$provisioned_src"
    elif [ -d "$recorded/.git" ]; then
      DEST="$recorded"
    fi
  fi
  recorded_manifest="$(awk -F= '$1 == "MACS3_BUILD_MANIFEST" {sub(/^[^=]*=/, ""); print; exit}' "$PROVISIONED" 2>/dev/null || true)"
  MANIFEST="${MACS3_BUILD_MANIFEST:-${recorded_manifest:-$(dirname "$DEST")/oracle-build-manifest.sha256}}"
  [ -d "$DEST/.git" ] || { echo "not provisioned: $DEST" >&2; exit 1; }
  MACS3_SRC="$DEST" MACS3_BUILD_MANIFEST="$MANIFEST" bash "$REPO_ROOT/oracle/verify_oracle_clean.sh"
  exit $?
fi

echo "pinning MACS3 $VERSION at $COMMIT"
if [ ! -e "$DEST" ]; then
  mkdir -p "$(dirname "$DEST")"
  git clone --quiet --filter=blob:none --no-checkout "$URL" "$DEST"
  git -C "$DEST" fetch --quiet --depth 1 origin "$COMMIT"
  git -C "$DEST" checkout --quiet "$COMMIT"
elif [ ! -d "$DEST/.git" ]; then
  echo "provision_oracle: refusing to replace existing non-git path: $DEST" >&2
  exit 1
fi

got="$(git -C "$DEST" rev-parse HEAD)"
if [ "$got" != "$COMMIT" ]; then
  echo "provision_oracle: existing source is at $got, expected $COMMIT; preserving it" >&2
  exit 1
fi
dirty="$(git -C "$DEST" status --porcelain)"
if [ -n "$dirty" ]; then
  echo "provision_oracle: refusing to build from a modified source tree:" >&2
  echo "$dirty" >&2
  exit 1
fi

git -C "$DEST" submodule update --init --depth 1 MACS3/fermi-lite/lib

PY="${PYTHON:-python3}"
if [ ! -x "$VENV/bin/python" ]; then
  mkdir -p "$(dirname "$VENV")"
  "$PY" -m venv "$VENV"
fi
VENV_PY="$VENV/bin/python"
actual_minor="$("$VENV_PY" -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')"
if [ "$actual_minor" != "$PYTHON_MINOR" ]; then
  echo "provision_oracle: existing venv uses Python $actual_minor; expected $PYTHON_MINOR (preserving it)" >&2
  exit 1
fi

CORE_LOCK="$(dirname "$DEST")/oracle-core-constraints.txt"
cat > "$CORE_LOCK" <<EOF
numpy==$NUMPY
scipy==$SCIPY
scikit-learn==$SKLEARN
hmmlearn==$HMMLEARN
EOF
# MACS3's requirements contain Cython<3.1, while the pristine generated C files
# at this commit identify Cython 3.3.0. Install the other declared runtime/build
# requirements with exact numerical constraints, then install that compiler.
REQ_NO_CYTHON="$(dirname "$DEST")/oracle-requirements.txt"
sed '/^[[:space:]]*[Cc]ython[[:space:]]*[<=>~!]/d' "$DEST/requirements.txt" > "$REQ_NO_CYTHON"
"$VENV/bin/pip" install --quiet --disable-pip-version-check \
  --constraint "$CORE_LOCK" -r "$REQ_NO_CYTHON"
"$VENV/bin/pip" install --quiet --disable-pip-version-check --no-deps \
  "Cython==$CYTHON_BUILD"
"$VENV/bin/pip" install --quiet --disable-pip-version-check \
  --no-deps --no-build-isolation --editable "$DEST"

# Hash every compiled extension immediately after a clean build. This manifest
# sits outside the oracle checkout and lets verification detect later edits to
# ignored build artifacts as well as tracked source changes.
mkdir -p "$(dirname "$MANIFEST")"
{
  echo "# MACS3_COMMIT=$COMMIT"
  echo "# MACS3_VERSION=$VERSION"
  echo "# PYTHON=$actual_minor"
  echo "# NUMPY=$NUMPY"
  echo "# SCIPY=$SCIPY"
  echo "# SKLEARN=$SKLEARN"
  echo "# HMMLEARN=$HMMLEARN"
  echo "# CYTHON_BUILD=$CYTHON_BUILD"
  echo "# OPENBLAS_CORETYPE=$OPENBLAS_CORETYPE"
  (cd "$DEST" && find MACS3 -type f -name '*.so' -print0 | sort -z | xargs -0 -r sha256sum)
} > "$MANIFEST"
[ -s "$MANIFEST" ] || { echo "provision_oracle: no compiled extensions found" >&2; exit 1; }

MACS3_SRC="$DEST" MACS3_BUILD_MANIFEST="$MANIFEST" \
  ORACLE_SO_BACKUP_DIR="${ORACLE_SO_BACKUP_DIR:-/tmp/oracle_so_backup}" \
  bash "$REPO_ROOT/oracle/verify_oracle_clean.sh"

cat > "$PROVISIONED" <<EOF
MACS3_VERSION=$VERSION
MACS3_COMMIT=$COMMIT
MACS3_SRC=$DEST
MACS3_PATH=$DEST/MACS3/__init__.py
MACS3_VENV=$VENV
MACS3_BUILD_MANIFEST=$MANIFEST
PYTHON=$actual_minor
NUMPY=$NUMPY
SCIPY=$SCIPY
SKLEARN=$SKLEARN
HMMLEARN=$HMMLEARN
CYTHON_BUILD=$CYTHON_BUILD
OPENBLAS_CORETYPE=$OPENBLAS_CORETYPE
PROVISIONED_BY=oracle/provision_oracle.sh
EOF

if [ -n "${GITHUB_ENV:-}" ]; then
  {
    echo "MACS3_SRC=$DEST"
    echo "MACS3_PATH=$DEST/MACS3/__init__.py"
    echo "MACS3_VENV=$VENV"
    echo "MACS3_BUILD_MANIFEST=$MANIFEST"
    echo "MACS3_ORACLE_DIR=$DEST"
    # The entry point and the interpreter that can import it. Later steps and the
    # crate tests read these rather than looking for a `macs3` on PATH or assuming a
    # checkout location: the venv's console script carries an absolute shebang, so it
    # only works on the machine that built it.
    echo "MACS3_ORACLE_BIN=$VENV/bin/macs3"
    echo "MACS3_ORACLE_PYTHON=$VENV/bin/python"
    echo "PYTHON=$VENV/bin/python"
    echo "OPENBLAS_CORETYPE=$OPENBLAS_CORETYPE"
    echo "PATH=$VENV/bin:$PATH"
  } >> "$GITHUB_ENV"
fi

echo "oracle ready: $DEST (Python $actual_minor, NumPy $NUMPY, SciPy $SCIPY, sklearn $SKLEARN, hmmlearn $HMMLEARN)"
