#!/usr/bin/env bash
# Differential for the `callvar` variant-calling kernel against pinned MACS3 3.0.5.
#
# This is the only oracle that can exercise `PosReadsInfo`/`PeakVariants`/
# `RACollection`. All three are `cython.cclass`, so their methods are not reachable
# from Python and there is no way to call them one function at a time -- the only
# surface is a full `callvar` run. Hence this compares whole VCF bodies rather than
# individual statistics.
#
# Inputs are upstream's own test data (`test/CTCF_PE_ChIP_chr22_50k.bam` etc.) read
# straight out of the pinned oracle tree, so the fixture never has to be committed.
#
#   bash oracle/check_callvar.sh [--keep] [--regen] [--fermi | --fermi-on]
#
#     --keep   leave the work directory in place and print its path
#     --regen  overwrite crates/macs-callvar/tests/data/callvar_variants.golden
#
# The golden holds the 16 non-`#` records. The header is deliberately not compared
# here: it embeds the run date, and it is already pinned byte-for-byte by
# crates/macs-callvar/tests/header_parity.rs.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCK="$ROOT/oracle/ENV.lock"
[ ! -f "$ROOT/oracle/ENV.provisioned" ] || LOCK="$ROOT/oracle/ENV.provisioned"
GOLDEN="$ROOT/crates/macs-callvar/tests/data/callvar_variants.golden"
GOLDEN_NOFERMI="$ROOT/crates/macs-callvar/tests/data/callvar_variants_nofermi.golden"
BIN="$ROOT/target/release/macs3-rs"

KEEP=0
REGEN=0
# Default to no assembly; the two switches also check automatic and forced
# assembly against fresh upstream records through the bundled fermi-lite bridge.
FERMI_MODE=off
for a in "$@"; do
    case "$a" in
        --keep) KEEP=1 ;;
        --regen) REGEN=1 ;;
        --fermi) FERMI_MODE=auto ;;
        --fermi-on) FERMI_MODE=on ;;
        *) echo "usage: $0 [--keep] [--regen] [--fermi | --fermi-on]" >&2; exit 2 ;;
    esac
done

MACS3_INIT="${MACS3_PATH:-$(awk -F= '/^MACS3_PATH=/{print $2}' "$LOCK")}"
if [ -z "$MACS3_INIT" ]; then
    echo "no MACS3_PATH in $LOCK" >&2
    exit 1
fi
ORACLE_SRC="$(dirname "$(dirname "$MACS3_INIT")")"

if [ ! -x "$BIN" ]; then
    echo "building macs3-rs ..."
    (cd "$ROOT" && cargo build --release --workspace --bin macs3-rs) || exit 1
fi

WORK="$(mktemp -d)"
cleanup() { [ "$KEEP" -eq 1 ] || rm -rf "$WORK"; }
trap cleanup EXIT

PEAKS="$ORACLE_SRC/test/callvar_testing.narrowPeak"
TBAM="$ORACLE_SRC/test/CTCF_PE_ChIP_chr22_50k.bam"
CBAM="$ORACLE_SRC/test/CTCF_PE_CTRL_chr22_50k.bam"
for f in "$PEAKS" "$TBAM" "$CBAM"; do
    [ -f "$f" ] || { echo "missing oracle test input: $f" >&2; exit 1; }
done

echo "== oracle callvar =="
( cd "$ORACLE_SRC" && PYTHONPATH="$ORACLE_SRC" python3 bin/macs3 callvar \
    -F "$FERMI_MODE" -b "$PEAKS" -t "$TBAM" -c "$CBAM" -o "$WORK/oracle.vcf" ) > "$WORK/oracle.log" 2>&1
ORACLE_RC=$?
if [ "$ORACLE_RC" -ne 0 ]; then
    echo "oracle callvar failed (rc=$ORACLE_RC)" >&2
    tail -20 "$WORK/oracle.log" >&2
    exit 1
fi

echo "== macs3-rs callvar =="
"$BIN" callvar -F "$FERMI_MODE" -b "$PEAKS" -t "$TBAM" -c "$CBAM" -o "$WORK/ours.vcf" > "$WORK/ours.log" 2>&1
OURS_RC=$?

grep -v '^#' "$WORK/oracle.vcf" > "$WORK/oracle.body"
[ -f "$WORK/ours.vcf" ] && grep -v '^#' "$WORK/ours.vcf" > "$WORK/ours.body" || : > "$WORK/ours.body"

if [ "$REGEN" -eq 1 ]; then
    case "$FERMI_MODE" in
        off) REGEN_GOLDEN="$GOLDEN_NOFERMI" ;;
        auto) REGEN_GOLDEN="$GOLDEN" ;;
        on) echo "no recorded golden target for forced assembly" >&2; exit 2 ;;
    esac
    cp "$WORK/oracle.body" "$REGEN_GOLDEN"
    echo "regenerated $REGEN_GOLDEN ($(wc -l < "$REGEN_GOLDEN") records)"
    exit 0
fi

if [ "$OURS_RC" -ne 0 ] || [ ! -s "$WORK/ours.body" ]; then
    echo "FAIL: macs3-rs callvar produced no variant records (rc=$OURS_RC)"
    tail -5 "$WORK/ours.log" >&2
    exit 1
fi

# When the default (`off`) run is being checked, also verify the *other* mode.
#
# These are genuinely different answers -- upstream re-calls indel and
# reference-biased peaks from the assembly rather than adding to them, and `off`
# gives 22 records against `auto`'s 16 -- so `--fermi auto` is compared against its own
# golden instead of against `off`. Both must be exact.
if [ "$FERMI_MODE" != "auto" ]; then
    echo
    echo "== cross-check: --fermi auto against its own golden =="
    AUTO="$WORK/auto.vcf"
    "$BIN" callvar -F auto -b "$PEAKS" -t "$TBAM" -c "$CBAM" -o "$AUTO" \
        > "$WORK/auto.log" 2>&1
    ARC=$?
    if [ "$ARC" -ne 0 ] || [ ! -e "$AUTO" ]; then
        echo "FAIL: --fermi auto failed (rc=$ARC); it is implemented and must succeed"
        tail -3 "$WORK/auto.log" >&2
        exit 1
    fi
    grep -v '^#' "$AUTO" > "$WORK/auto.body"
    sort "$WORK/auto.body" -o "$WORK/auto.body.sorted"
    grep -v '^#' "$GOLDEN" > "$WORK/golden.auto"
    sort "$WORK/golden.auto" -o "$WORK/golden.auto.sorted"
    AUTO_BAD=0
    while IFS= read -r line; do
        want_n=$(grep -cxF "$line" "$WORK/golden.auto.sorted" || true)
        got_n=$(grep -cxF "$line" "$WORK/auto.body.sorted" || true)
        [ "$got_n" -eq "$want_n" ] || { AUTO_BAD=$((AUTO_BAD + 1)); printf 'DIFF --fermi auto: %s\n' "$(echo "$line" | cut -f1,2,4,5)"; }
    done < "$WORK/golden.auto"
    AUTO_EXTRA=$(comm -13 <(sort -u "$WORK/golden.auto") <(sort -u "$WORK/auto.body") | wc -l)
    printf -- '--fermi auto: %s/%s records identical\n' \
        "$(( $(wc -l < "$WORK/auto.body.sorted") - AUTO_BAD - AUTO_EXTRA ))" \
        "$(wc -l < "$WORK/golden.auto.sorted")"
    [ "$AUTO_BAD" -eq 0 ] && [ "$AUTO_EXTRA" -eq 0 ] || exit 1
fi

# Compare as a *multiset*: `callvar_testing.narrowPeak` holds two byte-identical
# peaks, and upstream writes every variant inside both of them -- three records
# legitimately appear twice. A membership test would pass a port that emitted them
# once, so counts are part of the comparison.
sort "$WORK/ours.body" -o "$WORK/ours.body.sorted"
sort "$WORK/oracle.body" -o "$WORK/oracle.body.sorted"

TOTAL=0
BAD=0
while IFS= read -r line; do
    TOTAL=$((TOTAL + 1))
    want_n=$(grep -cxF "$line" "$WORK/oracle.body.sorted" || true)
    got_n=$(grep -cxF "$line" "$WORK/ours.body.sorted" || true)
    label="$(echo "$line" | cut -f1,2,4,5,7)"
    if [ "$got_n" -eq "$want_n" ]; then
        printf 'ok   %s x%s\n' "$label" "$want_n"
    else
        BAD=$((BAD + 1))
        printf 'DIFF %s\n     oracle x%s, ours x%s\n     %s\n' "$label" "$want_n" "$got_n" "$line"
    fi
done < "$WORK/oracle.body"

EXTRA=$(comm -13 <(sort -u "$WORK/oracle.body") <(sort -u "$WORK/ours.body") | wc -l)
if [ "$EXTRA" -gt 0 ]; then
    echo "FAIL: $EXTRA record(s) we emit that the oracle does not:"
    comm -13 <(sort -u "$WORK/oracle.body") <(sort -u "$WORK/ours.body") | cut -f1,2,4,5,7
    BAD=$((BAD + EXTRA))
fi

echo
echo "$((TOTAL - BAD))/$TOTAL records identical"
[ "$BAD" -eq 0 ] || exit 1
