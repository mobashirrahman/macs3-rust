#!/usr/bin/env bash
# Byte-compare `predictd` and `randsample` against the pinned oracle.
#
# These two commands are not in the golden corpus (which is callpeak-only), so their
# equivalence was previously unchecked even though both had accept/reject divergences
# from upstream:
#
#   F199  `randsample` required `-s/--tsize`. Upstream auto-detects the tag size from
#         the first <=10 usable lines (`Parser.tsize()`), so every invocation without
#         `-s` was rejected here and accepted there.
#   F200  `predictd` exited 1 when no model could be fitted. Upstream catches
#         `NotEnoughPairsException`, warns, and exits 0 with no output file.
#   F201  `*_model.r` was not byte-identical: `smooth` must scale each window element
#         by 1/n before summing, not sum and divide once.
#
# Each case is a *pair* of runs; a divergence in exit status counts as a failure just
# as a byte difference does.
set -uo pipefail
export OPENBLAS_CORETYPE="${OPENBLAS_CORETYPE:-Haswell}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OURS="$ROOT/target/release/macs3-rs"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

pass=0
fail=0

# Deterministic tag-size / model fixtures.
cat > "$WORK/se.bed" <<'EOF'
chrM	52	203	0	0	+
EOF
for i in $(seq 1 40); do
	printf 'chrM\t%d\t%d\t0\t0\t%s\n' "$((i * 97))" "$((i * 97 + 151))" \
		"$([ $((i % 2)) -eq 0 ] && echo + || echo -)"
done >> "$WORK/se.bed"

# 150 sites of paired +/- clusters: enough for PeakModel to fit (>= 100 paired peaks).
python3 - "$WORK/model.bed" <<'PYEOF'
import random, sys
random.seed(5)
with open(sys.argv[1], 'w') as f:
    for site in range(150):
        x = 1000 + site * 3000
        for _ in range(300):
            f.write(f"chr1\t{x + random.randint(0, 150)}\t{x + 251 + random.randint(0, 150)}\t0\t0\t+\n")
        for _ in range(300):
            f.write(f"chr1\t{x + 200 + random.randint(0, 150)}\t{x + 451 + random.randint(0, 150)}\t0\t0\t-\n")
PYEOF

awk '{print $1"\t"$2"\t"$3}' \
	"$ROOT/tests/fixtures/frag_basic/barcode_fragments/treat.frag" > "$WORK/pe.bedpe"

# run_case <name> <subcommand> <args...>
run_case() {
	local name="$1" sub="$2"
	shift 2
	local u="$WORK/up_$name" r="$WORK/rs_$name"
	mkdir -p "$u" "$r"

	local ua=() ra=()
	local a
	for a in "$@"; do
		ua+=("${a//@OUT@/$u}")
		ra+=("${a//@OUT@/$r}")
	done

	PYTHONPATH="$ORACLE_SRC" timeout 900 "$ORACLE_PY" "$ORACLE_SRC/bin/macs3" \
		"$sub" "${ua[@]}" >"$WORK/$name.up.log" 2>&1
	local urc=$?
	"$OURS" "$sub" "${ra[@]}" >"$WORK/$name.rs.log" 2>&1
	local rrc=$?

	local detail=""
	if [ "$urc" != "$rrc" ]; then
		detail="exit $urc vs $rrc"
	fi

	# Compare the produced trees, ignoring the directory name itself.
	if [ ! -d "$u" ] || [ ! -d "$r" ]; then
		detail="$detail (missing outdir)"
	else
		local diffout
		diffout="$(diff -r "$u" "$r" 2>&1)"
		if [ -n "$diffout" ]; then
			detail="$detail output differs: $(echo "$diffout" | head -3 | tr '\n' ' ')"
		fi
	fi

	if [ -z "$detail" ]; then
		pass=$((pass + 1))
		printf '  ok       %s\n' "$name"
	else
		fail=$((fail + 1))
		printf '  FAIL     %s -- %s\n' "$name" "$detail"
	fi
}

# The pinned oracle root, derived from ENV.lock's MACS3_PATH so the checker follows
# the lockfile rather than a path baked into the script.
LOCK="$ROOT/oracle/ENV.lock"
[ ! -f "$ROOT/oracle/ENV.provisioned" ] || LOCK="$ROOT/oracle/ENV.provisioned"
ORACLE_SRC="${MACS3_SRC:-$(dirname "$(dirname "$(grep '^MACS3_PATH=' "$LOCK" | cut -d= -f2-)")")}"
# Run upstream under the pinned interpreter when one is provisioned: `*_model.r` is
# only byte-identical against the pinned NumPy, and a bare `python3` is whatever the
# host has first on PATH.
ORACLE_PY="python3"
VENV="$(grep '^MACS3_VENV=' "$LOCK" | cut -d= -f2-)"
[ -z "$VENV" ] || [ ! -x "$VENV/bin/python" ] || ORACLE_PY="$VENV/bin/python"

echo "predictd / randsample oracle comparison"

# --- predictd ---------------------------------------------------------------

# A dataset that cannot yield 100 paired peaks: upstream warns and exits 0 with no
# output file at all. This is the F200 case.
run_case predictd_no_model predictd \
	-i "$WORK/se.bed" -g 500000 --outdir @OUT@
# A dataset that does fit: exercises the *_model.r byte-comparison (F201).
run_case predictd_model predictd \
	-i "$WORK/model.bed" -g 500000 --outdir @OUT@
# Paired-end mode never fits a model at all; it prints the mean insertion length.
run_case predictd_pe predictd \
	-i "$WORK/pe.bedpe" -f BEDPE --outdir @OUT@

# --- randsample -------------------------------------------------------------

# No -s: upstream detects the tag size (F199). With a fixed seed the sample is
# reproducible; without one NumPy draws from global entropy and can never match.
run_case randsample_autotsize randsample \
	-i "$WORK/se.bed" -f BED -p 50 --seed 42 -o out.bed --outdir @OUT@
run_case randsample_explicit_tsize randsample \
	-i "$WORK/se.bed" -f BED -p 50 -s 200 --seed 42 -o out.bed --outdir @OUT@
run_case randsample_by_number randsample \
	-i "$WORK/se.bed" -f BED -n 20 --seed 7 -o out.bed --outdir @OUT@

echo
echo "$pass identical, $fail differing"
[ "$fail" -eq 0 ]
