#!/usr/bin/env bash
# Replay the recorded stage fixtures through `macs3-rs` and emit our own
# `stages.json`, so `macs-compare --stages` has both sides.
#
# The argument list per fixture must match what `oracle/record_stages.py` used, or the
# comparison is not like-for-like; that list lives in one place in the Python recorder,
# so this reads it back rather than duplicating the choices.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-/tmp/our_stages}"
rm -rf "$OUT"
mkdir -p "$OUT"

replay() {
  local group="$1" name="$2" mode="$3"
  local fix="$REPO/tests/fixtures/$group/$name"
  local dst="$OUT/$group/$name/$mode"
  [ -d "$fix" ] || { echo "skip  $group/$name (no fixture)"; return 0; }
  mkdir -p "$dst"
  local gsize
  gsize="$(awk -F'[\t ]+' '!/^#/ && NF>=2 {s+=$2} END {print s+0}' "$fix/genome.txt" 2>/dev/null || echo 2000000)"
  [ "$gsize" -gt 0 ] || gsize=2000000
  local treat ctrl fmt
  case "$mode" in
    se)   treat="$fix/treat.bed"; ctrl="$fix/ctrl.bed";   fmt=BED  ;;
    pe)   treat="$fix/treat.bedpe"; ctrl="$fix/ctrl.bedpe"; fmt=BEDPE ;;
    frag) treat="$fix/treat.frag"; ctrl="$fix/ctrl.frag"; fmt=FRAG ;;
  esac
  local extra=()
  [ "$mode" = "se" ] && extra=(--nomodel --extsize 200)
  local cargs=()
  [ -f "$ctrl" ] && cargs=(-c "$ctrl")
  MACS3_RS_DUMP_STAGES="$dst" "$REPO/target/release/macs3-rs" callpeak \
      -n "ours_$mode" -g "$gsize" -t "$treat" "${cargs[@]}" -f "$fmt" \
      "${extra[@]}" --outdir "$dst" >/dev/null 2>&1 || {
        echo "FAIL  $group/$name [$mode] (macs3-rs exited $?)"; return 0; }
  echo "ok    $group/$name [$mode]"
}

replay se_basic gauss_two_peaks se
replay pe_basic gauss_fragments pe
replay pe_basic atac_short pe
replay pe_basic nucleosome_ladder pe
replay frag_basic barcode_fragments frag
echo "wrote $OUT"
