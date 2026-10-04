#!/usr/bin/env bash
# Peak-coordinate differential: the gate above `oracle/run_e2e.sh`.
#
# run_e2e.sh checks the *signals* (treatment pileup, control local lambda) against
# upstream's own `--bdg`. This script checks the *output*: peak count, peak
# boundaries and summit coordinates against the golden `*.xls`.
#
# Both score mechanisms are run, because neither is at parity yet (F48):
#
#   poisson  -- macs_score::pscore_track, which is upstream's actual mechanism
#   subtract -- a pointwise subtraction, known-wrong (F46) but it reproduces peak
#               boundaries on more fixtures today
#
# Reporting both is the point. A single number would hide the divergence.
#
#   oracle/run_peak_e2e.sh                  # the default fixture list
#   oracle/run_peak_e2e.sh --all            # every single-end fixture with a golden XLS
#   oracle/run_peak_e2e.sh se_model/realistic
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/target/release/macs-callpeak-e2e"

declare -a FIXTURES=()
VARIANT=default
if [ "${1:-}" = "--variant" ]; then
    VARIANT="${2:-default}"; shift 2
fi
if [ "${1:-}" = "--all" ]; then
    while IFS= read -r fx; do
        # SE fixtures have `treat.bed`; PE fixtures have `treat.bedpe`. Both are
        # driven (F89); only FRAG fixtures are still skipped, since `--format
        # FRAG` needs the weighted-fragment path.
        [ -f "$ROOT/tests/fixtures/$fx/treat.bed" ] || [ -f "$ROOT/tests/fixtures/$fx/treat.bedpe" ] || continue
        [ -f "$ROOT/tests/fixtures/$fx/treat.frag" ] && continue
        FIXTURES+=("$fx")
    done < <(find "$ROOT/tests/golden" -mindepth 4 -maxdepth 4 -name '*_peaks.xls' \
                 | sed "s|^$ROOT/tests/golden/||" | cut -d/ -f1-2 | sort -u)
elif [ -n "${1:-}" ]; then
    FIXTURES=("$1")
else
    FIXTURES=(
        se_model/realistic
        se_model/sharp_spikes
        se_model/spikes_only
        se_basic/gauss_two_peaks
        se_shapes/bias
        se_dup/dup_rate_0.25
    )
fi

echo "building harness..."
(cd "$ROOT" && cargo build -q --release --bin macs-callpeak-e2e) || exit 1

for mech in poisson subtract; do
    echo
    echo "=== mechanism: $mech ==="
    printf '%-34s %5s %5s %7s %7s %9s\n' fixture rust upstream bnd_ok summits worst_bnd
    t_rust=0 t_up=0 t_bnd=0 t_summ=0 t_worst=0 t_fix=0 t_full=0 t_skipped=0
    for fx in "${FIXTURES[@]}"; do
        dir="$ROOT/tests/fixtures/$fx"
        xls=$(ls "$ROOT/tests/golden/$fx/$VARIANT/"*_peaks.xls 2>/dev/null | head -1)
        [ -d "$dir" ] && [ -n "$xls" ] || continue
        # fixtures upstream rejects have no peaks worth comparing
        golden_cfg="$ROOT/tests/golden/$fx/$VARIANT/command.json"
        golden="$golden_cfg"
        if [ -f "$golden" ]; then
            rc=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['returncode'])" "$golden")
            [ "$rc" = "0" ] || { printf '%-34s %s\n' "$fx" "skipped (upstream rc=$rc)"; t_skipped=$((t_skipped+1)); continue; }
        fi

        # F88: replay the fixture's ACTUAL upstream flags, not a hardcoded set.
        # The effective genome size feeds the p-value -> q-value map, so forcing
        # `--gsize 2000000` changed where the q cutoff is crossed: 6 `gonechrom_*`
        # fixtures were generated with `-g 12000` and this port, run at 2000000,
        # called each peak 1-8 bp wider on the left. Every other flag comes from
        # `command.json` too, so a fixture generated with different --slocal /
        # --llocal / --broad cannot silently diverge either.
        mapfile -t ARGS < <(python3 -c '
import json,shlex,sys
d=json.load(open(sys.argv[1]))
argv=shlex.split(d["command"].join(["",""])) if False else d["command"]
keep=[]
i=0
while i < len(argv):
    a=argv[i]
    # Only flags the differential harness actually implements are replayed. The
    # harness is always `--nomodel` (it takes --extsize/--slocal/--llocal
    # directly), and the upstream model path is a separate gate.
    if a == "-f":
        keep += ["--format", argv[i+1]]; i+=2; continue
    if a in ("--call-summits","--broad","--nolambda"):
        keep.append(a); i+=1; continue          # boolean flags take no value
    if a in ("--broad-cutoff","--max-gap","--min-length","--pvalue"):
        keep += [a, argv[i+1]]; i+=2; continue          # boolean flags take no value
    if a in ("-g","--gsize","--extsize","--slocal","--llocal","--qvalue"):
        # normalise upstream short forms to the long form the harness parses
        long = {"-g": "--gsize"}.get(a, a)
        keep += [long, argv[i+1]]; i+=2; continue
    if a.startswith("--gsize=") or a.startswith("--extsize=") or a.startswith("--slocal=") \
       or a.startswith("--llocal=") or a.startswith("--qvalue=") or a.startswith("--pvalue=") \
       or a.startswith("--max-gap=") or a.startswith("--min-length="):
        keep.append(a); i+=1; continue
    i+=1
print("\n".join(keep))' "$golden_cfg")
        out=$(CALLPEAK_SCORE=$mech "$BIN" "$dir" "${ARGS[@]}" --xls "$xls" 2>/dev/null)
        nr=$(echo "$out" | sed -n '1s/[^0-9]*\([0-9]*\).*/\1/p')
        nu=$(echo "$out" | sed -n '2s/[^0-9]*\([0-9]*\).*/\1/p')
        # count peaks whose boundaries agree, and the largest boundary gap
        # F87: on full agreement the harness prints a single `PEAK MATCH` line and
        # no per-peak lines, so the per-peak regex below finds nothing and reports
        # `0 0 0` -- a vacuous zero that looks like a failure in the totals. Count
        # `PEAK MATCH` as full parity for the whole fixture.
        if echo "$out" | grep -q "CALLPEAK PEAK MATCH"; then
            printf '%-34s %5s %5s %7s %7s %9s\n' "$fx" "$nr" "$nu" "$nu" "$nu" 0
            t_fix=$((t_fix+1))
            t_rust=$((t_rust+nr)); t_up=$((t_up+nu))
            t_bnd=$((t_bnd+nu)); t_summ=$((t_summ+nu))
            t_full=$((t_full+1))
            continue
        fi
        stats=$(echo "$out" | python3 -c '
import re,sys
bnd_ok=summ_ok=worst=0
for line in sys.stdin:
    m=re.match(r"\s*peak \d+: rust (\S+):(\d+)-(\d+) summit (\d+) "
               r"vs upstream (\S+):(\d+)-(\d+) summit (\d+)", line)
    if not m: continue
    if (m.group(1),m.group(2),m.group(3))==(m.group(5),m.group(6),m.group(7)): bnd_ok+=1
    else: worst=max(worst,abs(int(m.group(2))-int(m.group(6))),abs(int(m.group(3))-int(m.group(7))))
    if m.group(4)==m.group(8): summ_ok+=1
print(bnd_ok,summ_ok,worst)')
        read -r bok sok worst <<<"$stats"
        printf '%-34s %5s %5s %7s %7s %9s\n' "$fx" "$nr" "$nu" "$bok" "$sok" "$worst"
        t_fix=$((t_fix+1))
        t_rust=$((t_rust+nr)); t_up=$((t_up+nu))
        t_bnd=$((t_bnd+bok)); t_summ=$((t_summ+sok))
        [ "$worst" -gt "$t_worst" ] && t_worst=$worst
        [ "$nr" -gt 0 ] && [ "$nr" = "$bok" ] && [ "$nr" = "$sok" ] && t_full=$((t_full+1))
    done
    echo
    printf '%-34s %5s %5s %7s %7s %9s\n' "TOTAL ($mech/$VARIANT)" "$t_rust" "$t_up" "$t_bnd" "$t_summ" "$t_worst"
    echo "  fixtures compared            : $t_fix  (skipped upstream-reject: $t_skipped)"
    echo "  fixtures with full parity    : $t_full"
done

cat <<'EOF'

bnd_ok    = peaks whose start/end match upstream exactly
summits   = peaks whose summit matches upstream exactly
worst_bnd = largest start/end disagreement, in bp

Both mechanisms are reported so neither can hide behind the other. `poisson` is
upstream's SE default; `subtract` is the other branch (F48). See
docs/upstream-findings.md.
EOF
