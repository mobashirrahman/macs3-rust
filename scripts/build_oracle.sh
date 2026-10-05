#!/usr/bin/env bash
# Build the pinned MACS3 oracle without deleting or resetting an existing source
# checkout or virtual environment.
#
#   ./scripts/build_oracle.sh /path/to/venv [source-dir]
set -euo pipefail

VENV="${1:?usage: build_oracle.sh <venv-dir> [source-dir]}"
SRC="${2:-$(dirname "$VENV")/macs3-src}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [ -e "$SRC" ] && [ ! -d "$SRC/.git" ]; then
  echo "FATAL: refusing to replace existing non-git source path: $SRC" >&2
  exit 1
fi

bash "$REPO_ROOT/oracle/provision_oracle.sh" --dest "$SRC" --venv "$VENV"

cat > "$(dirname "$VENV")/oracle-env.sh" <<EOF
# source this file to activate the pinned MACS3 oracle
export MACS3_SRC="$SRC"
export MACS3_PATH="$SRC/MACS3/__init__.py"
export MACS3_VENV="$VENV"
export PATH="$VENV/bin:\$PATH"
EOF

"$VENV/bin/macs3" --version
echo "==> oracle ready; source $(dirname "$VENV")/oracle-env.sh to activate"
