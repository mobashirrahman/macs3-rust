#!/usr/bin/env bash
# Differential check for the BAM reader against the pinned MACS3 oracle.
#
# For each region query it runs both implementations and diffs the outputs.
# The oracle side is `oracle/dump_bam_region.py`, which calls upstream's own
# `BAMaccessor.get_reads_in_region`; the Rust side is `macs-io-dump bam`, which
# renders the same columns. Byte equality here means the flag/MAPQ filter, the
# `rightmost` CIGAR rule, duplicate suppression and record order all agree.
#
# Usage: oracle/run_bam_e2e.sh [--regenerate]
set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BAM="$ROOT/tests/fixtures/bam/reads.bam"
GOLDEN="$ROOT/tests/golden/bam"
PY="${MACS3_ORACLE_PYTHON:-/scratch/mdra00001/tmp/opencode/macs3-venv/bin/python}"
DUMP="$ROOT/target/debug/macs-io-dump"

REGENERATE=0
[ "${1:-}" = "--regenerate" ] && REGENERATE=1

if [ ! -f "$BAM" ]; then
  echo "fixture missing: $BAM"
  echo "run: $PY $ROOT/oracle/make_bam_fixtures.py $ROOT/tests/fixtures/bam"
  exit 1
fi

cargo build -p macs-io --bin macs-io-dump 2>/dev/null || {
  echo "build failed"
  exit 1
}

# name:chrom:left:right:maxdup
QUERIES=(
  "chr1_1000_1200_d1:chr1:1000:1200:1"
  "chr1_1000_1200_d2:chr1:1000:1200:2"
  "chr1_1300_1500_d1:chr1:1300:1500:1"
  "chr1_1990_2060_d1:chr1:1990:2060:1"
  "chr1_1990_2060_d4:chr1:1990:2060:4"
  "chr2_3000_3120_d1:chr2:3000:3120:1"
  "chr1_full:chr1:0:20000:1"
  "chr2_full:chr2:0:8000:1"
)

mkdir -p "$GOLDEN"
tmp_o="$(mktemp)"; tmp_r="$(mktemp)"
trap 'rm -f "$tmp_o" "$tmp_r"' EXIT

pass=0; fail=0
printf '%-24s %10s  %s\n' "QUERY" "READS" "RESULT"
for q in "${QUERIES[@]}"; do
  name="${q%%:*}"; rest="${q#*:}"
  c="${rest%%:*}"; rest="${rest#*:}"
  l="${rest%%:*}"; rest="${rest#*:}"
  r="${rest%%:*}"; d="${rest##*:}"

  if ! "$PY" "$ROOT/oracle/dump_bam_region.py" "$BAM" "$c" "$l" "$r" "$d" > "$tmp_o" 2>/dev/null; then
    printf '%-24s %10s  %s\n' "$name" "-" "ORACLE ERROR"
    fail=$((fail+1)); continue
  fi
  if ! "$DUMP" bam "$BAM" "$c" "$l" "$r" "$d" > "$tmp_r" 2>&1; then
    printf '%-24s %10s  %s\n' "$name" "-" "RUST ERROR: $(head -1 "$tmp_r")"
    fail=$((fail+1)); continue
  fi

  n="$(awk -F'\t' '/^### n_reads/{print $2}' "$tmp_o")"
  if diff -q "$tmp_o" "$tmp_r" > /dev/null; then
    printf '%-24s %10s  %s\n' "$name" "$n" "identical"
    pass=$((pass+1))
    [ "$REGENERATE" = "1" ] && cp "$tmp_o" "$GOLDEN/$name.tsv"
  else
    printf '%-24s %10s  %s\n' "$name" "$n" "DIFFERS"
    diff "$tmp_o" "$tmp_r" | head -4 | sed 's/^/    /'
    fail=$((fail+1))
  fi
done

echo
echo "TOTAL (bam region)   $pass passed, $fail failed"
[ "$fail" -eq 0 ]