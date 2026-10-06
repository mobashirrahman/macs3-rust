#!/usr/bin/env bash
# Build the release tarball for `macs3-rs`.
#
# The artifact is what makes the project usable without a Rust toolchain, so it
# is built by a script rather than by hand: the same command on any checkout
# produces the same tarball contents, and the release notes can quote the commit
# and SHA-256 it came from.
#
# Usage:
#   scripts/package_release.sh [--out DIR] [--version X.Y.Z]
#
# Requires the pinned C toolchain for the vendored fermi-lite assembler, which
# `callvar` calls through a five-function FFI. Everything else is pure Rust.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT="$ROOT/dist"
VERSION=""

while [ $# -gt 0 ]; do
    case "$1" in
        --out) OUT="$2"; shift 2 ;;
        --version) VERSION="$2"; shift 2 ;;
        *) echo "usage: $0 [--out DIR] [--version X.Y.Z]" >&2; exit 2 ;;
    esac
done

# Default to the workspace version, and refuse to package a version the source
# does not agree with: a tarball named v0.1.0 built from a tree saying 0.2.0 is
# exactly the kind of mismatch this project exists to catch elsewhere.
if [ -z "$VERSION" ]; then
    VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)"
fi
WS_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)"
if [ "$VERSION" != "$WS_VERSION" ]; then
    echo "refusing to package: asked for $VERSION but the workspace is $WS_VERSION" >&2
    exit 2
fi
# The CLI inherits the version rather than restating it, so confirm that is still
# true -- a literal here would be a second source of truth that could disagree.
if ! grep -q '^version\.workspace = true$' crates/macs-cli/Cargo.toml; then
    echo "refusing to package: crates/macs-cli no longer inherits the workspace version" >&2
    exit 2
fi

echo "building macs3-rs $VERSION (release)..."
cargo build --release --bin macs3-rs

SRC_COMMIT="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
NAME="macs3-rs-v${VERSION}"
STAGE="$OUT/$NAME"
rm -rf "$STAGE"
mkdir -p "$STAGE"

# Debug info is enabled in the release profile so CI crash reports are useful;
# it is not useful in a distributed artifact and triples its size.
cp target/release/macs3-rs "$STAGE/macs3-rs"
strip --strip-debug "$STAGE/macs3-rs" 2>/dev/null || true
chmod +x "$STAGE/macs3-rs"

cp README.md LICENSE "$STAGE/"

# Record provenance inside the tarball. A reader can check that the binary they
# downloaded corresponds to a commit in the repository.
cat > "$STAGE/BUILD.txt" <<EOF
macs3-rs $VERSION
source commit : $SRC_COMMIT
compatibility: MACS3 3.0.5 (commit c5443190e3edfeb301cc94acf450e2b2c026a223)
target        : $(rustc -vV | sed -n 's/^host: //p')
built         : $(date -u +%Y-%m-%dT%H:%M:%SZ)

Requires glibc. Build from source for any other target:
  cargo install --path crates/macs-cli
EOF

TARBALL="$OUT/$NAME-x86_64-linux-gnu.tar.gz"
tar -czf "$TARBALL" -C "$OUT" "$NAME"
rm -rf "$STAGE"

echo
echo "artifact : $TARBALL"
echo "size     : $(du -h "$TARBALL" | cut -f1)"
echo "sha256   : $(sha256sum "$TARBALL" | cut -d' ' -f1)"
echo "commit   : $SRC_COMMIT"
echo
echo "smoke test:"
tar -xzf "$TARBALL" -C "$OUT"
"$OUT/$NAME/macs3-rs" --version
rm -rf "$OUT/$NAME"