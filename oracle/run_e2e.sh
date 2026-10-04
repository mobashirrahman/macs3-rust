#!/usr/bin/env bash
# Stage-by-stage differential for the real pileup pipeline.
#
# For each fixture: regenerate upstream's intermediate stages in-process via
# oracle/dump_stages.py, then run the Rust pipeline and compare the treatment
# pileup and control local-lambda bedGraphs that upstream's own `--bdg` writes.
#
#   oracle/run_e2e.sh                 # default fixture list
#   oracle/run_e2e.sh se_shapes/bias  # just one
#   oracle/run_e2e.sh --all           # every single-end fixture
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV="/scratch/mdra00001/tmp/opencode/macs3-venv"
WORK="${E2E_WORK:-/tmp/e2e_work}"
MACS3_SRC="${MACS3_SRC:-/scratch/mdra00001/tmp/opencode/macs3-src}"
# Must match the golden runs, which fix -g at 2e6 for every fixture.
GSIZE="${E2E_GSIZE:-2000000}"

# shellcheck disable=SC1091
source "$VENV/bin/activate"

declare -a FIXTURES=()
if [ "${1:-}" = "--all" ]; then
    while IFS= read -r d; do
        rel="${d#"$ROOT/tests/fixtures/"}"
        [ -f "$d/treat.bed" ] || continue
        [ -f "$d/ctrl.bed" ] || continue
        [ -f "$d/pe.bed" ] && continue
        FIXTURES+=("$rel")
    done < <(find "$ROOT/tests/fixtures" -mindepth 1 -maxdepth 2 -type d | sort)
elif [ -n "${1:-}" ]; then
    FIXTURES=("$1")
else
    FIXTURES=(
        se_model/realistic
        se_model/sharp_spikes
        se_model/spikes_only
        se_basic/gauss_two_peaks
        se_edge/disjoint_chromosomes
        se_shapes/bias
        se_dup/dup_rate_0.25
        se_edge/very_shallow
    )
fi

echo "building harness..."
(cd "$ROOT" && cargo build -q --release --bin macs-pipeline-e2e) || exit 1
BIN="$ROOT/target/release/macs-pipeline-e2e"

pass=0; fail=0; skip=0; failed=()
for fx in "${FIXTURES[@]}"; do
    dir="$ROOT/tests/fixtures/$fx"
    [ -d "$dir" ] || { echo "SKIP $fx (missing)"; continue; }

    # Some fixtures exist precisely to exercise upstream's rejection paths, and
    # upstream exits non-zero on them. There is no pileup to compare in that
    # case -- those are covered by the error-parity gate instead.
    golden="$ROOT/tests/golden/$fx/default/command.json"
    if [ -f "$golden" ]; then
        rc=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['returncode'])" "$golden")
        if [ "$rc" != "0" ]; then
            echo "SKIP $fx (golden rc=$rc -- upstream rejects this fixture)"
            skip=$((skip + 1)); continue
        fi
    fi

    out="$WORK/$(echo "$fx" | tr '/' '_')"
    rm -rf "$out"; mkdir -p "$out"

    # The golden runs use `-g <gsize> --nomodel --extsize 200`, so the toy
    # 1 Mb genomes never attempt model building (which needs >=100 paired peaks).
    if ! python3 "$ROOT/oracle/dump_stages.py" "$dir" se --out "$out" \
            --macs3-src "$MACS3_SRC" -- -g "$GSIZE" --nomodel --extsize 200 --bdg \
            >"$out/stages.log" 2>&1; then
        echo "FAIL $fx  (upstream stage dump failed; see $out/stages.log)"
        fail=$((fail + 1)); failed+=("$fx"); continue
    fi

    if "$BIN" "$dir" \
            --extsize 200 --gsize "$GSIZE" --slocal 1000 --llocal 10000 \
            --treat-pileup "$out/stages_treat_pileup.bdg" \
            --control-lambda "$out/stages_control_lambda.bdg" \
            >"$out/e2e.log" 2>&1; then
        t=$(grep -o 'max abs diff [0-9.e+-]*' "$out/e2e.log" | head -1)
        c=$(grep -o 'max abs diff [0-9.e+-]*' "$out/e2e.log" | tail -1)
        printf 'PASS %-34s treat %s | ctrl %s\n' "$fx" "${t#max abs diff }" "${c#max abs diff }"
        pass=$((pass + 1))
    else
        n=$(grep -o 'total mismatches: [0-9]*' "$out/e2e.log" | tail -1)
        printf 'FAIL %-34s %s\n' "$fx" "${n#total mismatches: }"
        fail=$((fail + 1)); failed+=("$fx")
    fi
done

echo
echo "e2e pileup: $pass passed, $fail failed, $skip skipped (upstream rejects) (of ${#FIXTURES[@]})"
if [ ${#failed[@]} -gt 0 ]; then
    printf 'open: %s\n' "${failed[@]}"
    exit 1
fi
