#!/bin/bash
# Build the pinned MACS3 oracle from scratch.
#
# Reproduces the reference implementation used by every golden vector and every
# differential fixture. Pinned to MACS3 3.0.5 at commit c544319; the commit SHA
# is verified after checkout so a silent upstream change cannot slip through.
#
#   ./scripts/build_oracle.sh /path/to/venv
#
# Notes:
#   * the `MACS3/fermi-lite/lib` git submodule is required, because
#     `fermi-lite/ksw.c` includes `lib/x86/sse2.h` from it and the build fails
#     with a confusing "fatal error" without it
#   * MACS3 requires Python >= 3.12 (setup.py:25)

set -euo pipefail

VENV="${1:?usage: build_oracle.sh <venv-dir>}"
MACS3_REPO="https://github.com/macs3-project/MACS.git"
MACS3_COMMIT="c544319"          # v3.0.5
MACS3_TAG="v3.0.5"
SRC="$(dirname "$VENV")/macs3-src"

echo "==> creating venv at $VENV"
python3 -m venv "$VENV"
"$VENV/bin/pip" install -q --upgrade pip setuptools wheel

echo "==> cloning MACS3"
rm -rf "$SRC"
git clone --quiet "$MACS3_REPO" "$SRC"
git -C "$SRC" checkout --quiet "$MACS3_COMMIT"
ACTUAL="$(git -C "$SRC" rev-parse --short HEAD)"
if [ "$ACTUAL" != "$MACS3_COMMIT" ]; then
  echo "FATAL: expected $MACS3_COMMIT, got $ACTUAL" >&2
  exit 1
fi
echo "    pinned at $ACTUAL (tag: $(git -C "$SRC" describe --tags 2>/dev/null || echo '?'))"

echo "==> fetching the simde submodule (needed by fermi-lite/ksw.c)"
git -C "$SRC" submodule update --init --depth 1 MACS3/fermi-lite/lib

echo "==> installing MACS3 (compiles ~25 Cython extensions, several minutes)"
"$VENV/bin/pip" install -e "$SRC" > "$VENV/../macs3-build.log" 2>&1 || {
  echo "FATAL: build failed, see $VENV/../macs3-build.log" >&2
  tail -30 "$VENV/../macs3-build.log" >&2
  exit 1
}

echo "==> verifying"
"$VENV/bin/python" - <<'PY'
import MACS3, importlib.metadata as md
from MACS3.Signal.Prob import poisson_cdf, binomial_cdf_inv
print("MACS3 distribution version:", md.version("MACS3"))
print("MACS3 path:", MACS3.__file__)
# a spot check that the compiled extensions are actually working
v = poisson_cdf(30, 10.0, False, True)
assert v == -7.09779, f"unexpected p-score {v}"
print("poisson_cdf(30, 10.0, upper, log10) =", v, "(expected -7.09779)")
print("binomial_cdf_inv(0.99, 1000, 1/293128983) =", binomial_cdf_inv(0.99, 1000, 1/293128983))
PY
"$VENV/bin/macs3" --version

cat > "$(dirname "$VENV")/oracle-env.sh" <<EOF
# source this to use the oracle
export MACS3_VENV="$VENV"
export PATH="$VENV/bin:\$PATH"
EOF
echo "==> done; write $(dirname "$VENV")/oracle-env.sh to activate"
