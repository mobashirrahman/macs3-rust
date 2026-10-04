#!/usr/bin/env bash
# Compare our `*_cutoff_analysis.txt` against upstream's recorded files, byte for byte.
#
# The feature is not only a report: the ladder's cutoffs are seeded into the AFDR
# histogram before the q-table is built, so a wrong ladder changes every q-score as
# well as the file. A format-only test would miss that; this runs the whole pipeline.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R="$REPO/target/release/macs3-rs"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail=0; ok=0

run() {
  local group="$1" name="$2" mode="$3" fmt="$4" treat="$5" ctrl="$6"
  local fix="$REPO/tests/fixtures/$group/$name"
  local gold="$REPO/tests/stages/$group/$name/$mode/stages_cutoff_analysis.txt"
  [ -f "$gold" ] || { echo "skip  $group/$name (no recorded file)"; return 0; }
  local gsize
  gsize="$(awk -F'[\t ]+' '!/^#/ && NF>=2 {s+=$2} END {print s+0}' "$fix/genome.txt")"
  [ "$gsize" -gt 0 ] || gsize=2000000
  local extra=(); [ "$mode" = "se" ] && extra=(--nomodel --extsize 200)
  local cargs=(); [ -f "$ctrl" ] && cargs=(-c "$ctrl")
  rm -rf "$WORK/o"
  "$R" callpeak -n ours -g "$gsize" -t "$treat" "${cargs[@]}" -f "$fmt" \
      "${extra[@]}" --cutoff-analysis --outdir "$WORK/o" >/dev/null 2>&1
  local ours="$WORK/o/ours_cutoff_analysis.txt"
  if [ ! -f "$ours" ]; then
    echo "FAIL  $group/$name [$mode]: no file written"; fail=$((fail+1)); return 0
  fi
  if diff -q "$gold" "$ours" >/dev/null; then
    echo "ok    $group/$name [$mode]"; ok=$((ok+1))
  else
    echo "FAIL  $group/$name [$mode]"; diff "$gold" "$ours" | head -8; fail=$((fail+1))
  fi
}

F="$REPO/tests/fixtures"
run se_basic gauss_two_peaks      se   BED   "$F/se_basic/gauss_two_peaks/treat.bed"            "$F/se_basic/gauss_two_peaks/ctrl.bed"
run pe_basic gauss_fragments      pe   BEDPE "$F/pe_basic/gauss_fragments/treat.bedpe"          "$F/pe_basic/gauss_fragments/ctrl.bedpe"
run pe_basic atac_short           pe   BEDPE "$F/pe_basic/atac_short/treat.bedpe"              "$F/pe_basic/atac_short/ctrl.bedpe"
run pe_basic nucleosome_ladder    pe   BEDPE "$F/pe_basic/nucleosome_ladder/treat.bedpe"       "$F/pe_basic/nucleosome_ladder/ctrl.bedpe"
run frag_basic barcode_fragments  frag FRAG  "$F/frag_basic/barcode_fragments/treat.frag"        "$F/frag_basic/barcode_fragments/ctrl.frag"

echo
echo "$ok identical, $fail differing"
[ "$fail" -eq 0 ]
