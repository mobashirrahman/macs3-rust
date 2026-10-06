# Upstream MACS3 3.0.5 findings

Reverse-engineering notes recorded while porting. Every entry cites the upstream
line it comes from and states how `macs3-rs` treats it. These are the reason the
port cannot be a "translation": upstream has several behaviours that are
mathematically wrong, one of which is a non-termination, and all of them are part
of the observable output.

Upstream reference: `macs3-project/MACS` @ `c544319` (v3.0.5).

---

## F1. `log10_poisson_cdf_Q_large_lambda` is correct, but only to five decimals

`MACS3/Signal/Prob.py:499` (Cython `log10_poisson_cdf_Q_large_lambda`).

The function sums the upper tail in log space and stops as soon as an
individual term changes the residue by less than `1e-5`:

```
while True:
    m += 1
    logy = logx + ln_lbd - log(m)
    pre_residue = residue
    residue = logspace_add(pre_residue, logy)
    if fabs(pre_residue - residue) < 1e-5:
        break
```

**I originally recorded this as an upstream bug and that was wrong.** The
`1e-5` increment threshold looks alarming because the terms *grow* while
`m < lambda`, but the loop necessarily walks past the mode, after which they
decay geometrically and the threshold is reached only once the tail is
exhausted. Checked against an independent direct summation over all 220
`(k, lambda)` combinations in the golden vector grid:

```
max | oracle p-score - true p-score | = 2.0e-5
```

which is exactly the `round(..., 5)` quantisation. So upstream is right.

What *is* load-bearing here:

* The result is rounded to 5 decimals, so **every p-score MACS produces is a
  multiple of `1e-5`**, and the smallest non-zero p-score is `1e-5`. A different
  rounding would change peak boundaries everywhere.
* The sign convention is inverted relative to the docstring — see F4.
* The cost is `O(lambda)` iterations, not `O(k)` — see F2.

**Consequence for the port:** p-scores must be **bit-identical**, and the golden
vector file `crates/macs-stats/tests/stats_vectors.tsv` enforces exactly that
across 4115 vectors.

---

## F2. `log10_poisson_cdf_Q_large_lambda` iterates O(lambda) times

`MACS3/Signal/ScoreTrack.py:69-88` and `MACS3/Signal/CallPeakUnit.py:55-73`.

```python
memcpy(&y_bits, &y, 4)
key = ((longlong)(uint)x << 32) | y_bits
score = -1 * poisson_cdf(x, y, False, True)
```

`y` (the local lambda) is a C `float`, so the cache key is the exact `f32` bit
pattern. Two lambdas that differ only below `f32` precision are the *same* key.

`macs3-rs` therefore:

* carries local lambda as `f32` on every MACS-compatible path, and
* keys its score cache on `(observed, lambda.to_bits())`.

This is why `macs-rle::BitKey` is bit-based and order-preserving rather than a
`Hash` on the float: a numeric-keyed map would merge `+0.0` with `-0.0` and
turn a valid lookup into a miss.

---

## F4. `poisson_cdf` in log mode returns `log10 p`, not `-log10 p`

Same function. The docstring says `-log10(p)`; the code computes
`(residue - lambda) / log(10)`, which is `log10 p <= 0`. Every call site negates.

`macs3-rs` reproduces the **code**. [`macs_stats::poisson_cdf`] is documented
accordingly so nobody "fixes" the sign later.

---

## F5. The `pv` array is right-endpoint indexed

`MACS3/Signal/PileupV2.py:344` (`_pileup_sorted_unit_as_list`) and its
consumer `_write_pv_to_bedGraph` (`:446`).

`pv` is a pair of arrays `(positions, values)`, and the convention is
**`values[i]` is the value on `[positions[i-1], positions[i])`**, with
`positions[-1] := 0`. This is fixed by the writer:

```python
pre: cython.int = 0
for i in range(l_pv):
    pos = p_ptr[i]
    value = v_ptr[i] * scale_factor
    if value < baseline_value:
        value = baseline_value
    fprintf(fh, b"%s\t%d\t%d\t%.5f\n", chrom_char, pre, pos, value)
    pre = pos
```

Two consequences that are easy to get wrong:

* the **first** entry describes the region *before* the first read, and
* nothing at all is represented after the last breakpoint — the format has no
  "persisting tail".

Reading it the other way round (`values[i]` on `[positions[i], positions[i+1])`)
shifts the entire signal by one interval, which no value-at-a-point spot check
would catch. `SignalTrack::from_breakpoints` / `to_breakpoints` implement the
right-endpoint convention, and
`crates/macs-pileup/tests/parity_vectors.rs` checks all 1776 oracle vectors
under it.

### F5a. Coincident start/end events are correct

The same sweep handles a fragment ending on the base the next one begins with

```python
else:
    i_s += 1
    i_e += 1   # both advance, pileup unchanged
```

**I initially recorded this as a serious upstream bug and that was wrong.** The
net depth delta really is zero, so the value is correct; all that changes is
that no breakpoint is emitted there, because the value does not change. With a
fixed extension `d` it fires for every pair of reads exactly `d` apart, and
with `--shift` for every pair `d +/- 2 * shift` apart, so the breakpoint list
can be much coarser than a naive RLE would be.

Worked example, reads at 5, 55, 105 with `d = 50` (fragments `[5,55)`,
`[55,105)`, `[105,155)`, true depth 1 across `[5,155)`):

| | MACS 3.0.5 | naive RLE |
|---|---|---|
| positions | `5, 155` | `5, 55, 105, 155` |
| values | `0.0, 1.0` | `0.0, 1.0, 1.0, 1.0` |
| value function | identical | identical |

Verified against the compiled oracle (`tests/pileup_vectors.tsv`, the
`pn_junction` rows) and pinned by the test
`coincident_start_and_end_emits_no_breakpoint`.

## F6. `over_two_pv_array` truncates the merge at the shorter breakpoint list

`MACS3/Signal/PileupV2.py:950-1032`.

```python
while i1 < l1 and i2 < l2:
    ...
```

The loop stops as soon as *either* input's breakpoint list is exhausted. The
final value emitted therefore persists over the region where the shorter track
has no breakpoints but the longer one does.

This is load-bearing for the local lambda. `pileup_a_chromosome_c`
(`MACS3/Signal/FixWidthTrack.py:775-852`) folds the `d`, `slocal` and `llocal`
window pileups together with `func="max"`, and the `d` window track is strictly
shorter than the `llocal` one. So `lambda_local` is only defined where the
shortest window track has breakpoints, and MACS's `lambda_bg` floor exists to
cover the rest.

`macs-rle::SignalTrack::zip` documents and reproduces this; `pad_to` is provided
for callers that need the true pointwise max instead.

---

## F7. `pduplication` accumulates in `float32`

`MACS3/Signal/Prob.py:688-698`. `sf` is `float32_t` and each `f64` summand is
truncated back to `f32` on every iteration, and the return value is `float32`.

`macs_stats::pduplication` does the same widening/narrowing dance explicitly.

---

## F8. `chisq_pvalue_e` loses all precision for `x > 40`

`MACS3/Signal/Prob.py:148-192`, via `ex20` which forces `exp(t)` to exactly `0`
for `t < -20`.

* `df == 2`: the answer is literally `ex20(-x/2)`, i.e. exactly `0.0` for
  `x > 40`, where the true value is `exp(-x/2)`.
* `df > 2`, `x/2 > 20`: the rescaled log-space loop starts from a truncated seed,
  losing several percent (e.g. `df=4, x=80`: upstream `0.0`, true `1.74e-16`).

Only `callvar` consumes these. `macs3-rs` reproduces the loss, because the VCF
`QUAL` field is defined by upstream's number.

---

## F9. `poisson_cdf_inv` is unreachable for `cdf < exp(-lambda)`

`MACS3/Signal/Prob.py:553-584`. `sumold` starts at `P(X = 0) = exp(-lambda)`, so
a `cdf` below that can never be bracketed and the function returns `maximum`
(1000 by default).

`--keep-dup auto` calls `binomial_cdf_inv`, not this function, so the pathology
does not reach the CLI. It is still reproduced and tested.

---

## F10. `poisson_pdf` overflows for moderate `k`

`MACS3/Signal/Prob.py:622-631` evaluates `exp(-a) * a**k / k!` directly, which is
`inf / inf = NaN` once `a**k` and `k!` both overflow — `k = 1000, a = 10` is
already past that. Upstream would raise `OverflowError`; MACS never calls this
function (p-scores go through `poisson_cdf`).

`macs_stats::poisson_pdf` returns `0.0` instead so it stays total, and
`poisson_pdf_log` is provided for the stable value. Neither is on a
MACS-compatible path.

---

## F11. `binomial_pdf` non-terminates for tiny `b` with `x <= a - x`

`MACS3/Signal/Prob.py:829-880`.

```python
for q in range(1, mn + 1):
    pdf *= (a - q + 1) * p / (mn - q + 1)
    if pdf < 1e-100:
        while pdf < 1e-3:
            pdf /= 1 - p
            t -= 1
```

The recovery loop divides by `1 - p`. When `p` is tiny, `1 - p` is
`~1.0 + p`, so the division *increases* `pdf` by about one part in `1/p` per
iteration. Lifting `pdf` from `1e-100` to `1e-3` therefore takes ~`1e97 / p`
iterations.

Reproduced with the compiled oracle:

| call | result |
|---|---|
| `binomial_pdf(25, 50, 1e-8)` | does not terminate |
| `binomial_pdf(25, 100, 1e-8)` | does not terminate |
| `binomial_pdf(0, 50, 1e-8)` | `0.99999999...` (no hang) |

MACS's own pipeline never hits this: `pduplication` and the `--keep-dup auto`
curve only ever evaluate `x in {0, 1, 2}` at tiny `b`, where the inner loop is
short and `pdf` never approaches `1e-100`. It is only reachable from ad-hoc
analysis of the function.

`macs3-rs` therefore bounds the recovery loop (returning the value it has once
the bound is hit) so the function is total. This is a **deliberate deviation**
with no observable effect: upstream produces no output for these inputs, so no
peak can differ.

The golden vector grid excludes these cases and says so.

---

## F12. `binomial_cdf` upper tail can also non-terminate

`MACS3/Signal/Prob.py:702-756` (`_binomial_cdf_r`). The upward loop

```python
i = argmax
while True:
    pdf *= (a - i) * b / (1 - b) / (i + 1)
    if pdf == 0.0:
        break
```

is only terminated by `pdf` reaching *exactly* `0.0`. For `i > a` the
multiplier `(a - i)` is negative, so `pdf` flips sign and grows in magnitude
instead of underflowing. The loop has no bound on `i`.

In practice `pdf` underflows to `+0.0` before `i` passes `a` for every parameter
combination MACS uses, so this is latent rather than active — but the
`if x < argmax` branch's *downward* loop
(`for i in range(argmax - 1, x, -1)`, i.e. indices `argmax-1 .. x+1`) is
transcribed exactly, because getting its endpoint wrong changes the sum.

---

## F14. The log-space Poisson path costs O(k + lambda) per call

`MACS3/Signal/Prob.py:499`. Two separate costs:

* an `O(k)` prefix sum for `ln(k!)` before the walk begins, and
* an `O(lambda)` walk, because the loop must pass the mode at `m ~ lambda`
  before the `1e-5` increment criterion can trigger (see F1/F2).

MACS never evaluates this with a large `k` — the observed argument is a read
depth — but a caller that feeds a garbage depth pays `O(k)` for it before the
work is even started. A depth of `2^31` takes about 30 s in Rust and far longer
in Python. `macs_stats` is a direct transcription and keeps both costs, because
short-circuiting either would change where the walk stops and therefore the
p-score (F1).

`macs_score::PScoreCache` is what makes this affordable in practice, and it
exists for exactly this reason: upstream caches on `(observed, f32 bits of
expectation)`, so each distinct depth/lambda pair is evaluated once.

---

## F15. `check_names` raises `TypeError` instead of reporting no common chromosomes

`MACS3/Commands/callpeak_cmd.py:31-42`.

```python
def check_names(treat, control, error_stream):
    tchrnames = set(treat.get_chr_names())
    cchrnames = set(control.get_chr_names())
    commonnames = tchrnames.intersection(cchrnames)
    if len(commonnames) == 0:
        error_stream("No common chromosome names can be found ...")
        error_stream("Chromosome names in treatment: %s" % ",".join(sorted(tchrnames)))
        ...
```

`get_chr_names()` returns **`bytes`**, so `",".join(...)` raises

```
TypeError: sequence item 0: expected str instance, bytes found
```

before the intended message is ever printed. So when treatment and control
share no chromosome names, upstream crashes with a traceback and exit code 1
instead of explaining the problem.

Verified against the compiled oracle on the fixture
`se_edge/disjoint_chromosomes` and `se_edge/one_sided_no_reads`.

`macs3-rs` behaviour: **fail cleanly with `MacsError::NoCommonChromosomes`,
exit non-zero, and create no output files.** The exception is a crash, not a
contract, so it is not reproduced; the *condition* is. It is recorded here
because the golden corpus contains the upstream crash, and CI asserts that the
Rust side fails the same way rather than diverging silently.

---

## F16. `--mfold` takes two values, not three

`bin/macs3:253`:

```python
group_bimodal.add_argument("-m", "--mfold", dest="mfold", type=int,
                           default=[5, 50], nargs=2, ...)
```

`--mfold 2 4 8` is an argparse usage error (exit code 2), *not* a
"low < mid < high" triple as the older documentation suggests. `-m/--mfold` is
also reused by `predictd` with the same `nargs=2`.

`macs3-rs` therefore defines `--mfold` with exactly two integer values and
rejects a third with a usage error, matching the exit code.

---

## F17. `PeakModel` cannot be trained on a synthetic miniature genome

`MACS3/Signal/PeakModel.py:98-138`.

```python
self.peaksize = 2 * self.bw
self.min_tags = int(round(float(self.treatment.total) *
                          self.lmfold * self.peaksize / self.gz / 2))
self.max_tags = int(round(float(self.treatment.total) *
                          self.umfold * self.peaksize / self.gz / 2))
```

and the strand peaks are found on tags that have **already been
de-duplicated** (the default `--keep-dup 1` keeps one tag per position per
strand, and `self.treatment.total` is the post-filter count).

Three consequences, all found by running the real MACS3 against generated
fixtures:

1. **`total` is a function of the number of *distinct positions*, not the read
   count.** A spike of 900 reads over 120 distinct positions collapses to 240
   tags, which moves `min_tags` by three orders of magnitude.
2. **A peak's depth must lie in `[min_tags, max_tags]`, and both scale with
   `total`.** With `S` sites of depth `D` on a genome of size `gz`:
   `total ~ S*D` and the requirement `D > total*lmfold*peaksize/gz/2` gives
   `gz > S * lmfold * peaksize / 2`. For `S = 150`, `lmfold = 5` and
   `peaksize = 600` that is `gz > 225 kb`, and MACS3 additionally demands
   `S >= 100`. The 48 kb `mini` genome is therefore structurally incapable of
   producing a model, no matter how the reads are arranged.
3. Even once the arithmetic works out, the per-site depth must land in a narrow
   band, and MACS3 still reports `Total number of paired peaks: 0` for every
   arrangement tried here (spike widths 25-400 bp, 10-12 sites per contig,
   12-1500 reads per site, 20 contigs, `gz = 2 Mb`).

**Consequence for the port:** the fragment-model gate (G8) is run against
MACS3's own CTCF test dataset — the corpus its documentation uses to explain the
model — rather than a synthetic fixture. The `se_model` fixture group is
generated and carried through every *other* stage with `--nomodel`, and is
explicitly marked `model_capable 0` in its manifest so the runner and the Rust
side agree.

This is also a useful compatibility fact in its own right: `macs3-rs` must
reproduce the `NotEnoughPairsException` and the exit code, because users hit it
constantly on small inputs.

---

## F18. The lowest-ranked p-score bucket always gets q = 0

`MACS3/Signal/ScoreTrack.py:525-543` (`make_pq_table`).

```python
for i in range(len(unique_values)):
    v = unique_values[i]
    ...
    if q <= 0:
        q = 0
        break
    pvalue2qvalue[v] = q
    pre_q = q
    k += ln
# bottom rank pscores all have qscores 0
for j in range(i, len(unique_values)):
    pvalue2qvalue[unique_values[j]] = 0
```

The second loop is written as if it only ran after a `break`, but when the
first loop runs to completion `i == len(unique_values) - 1`. So the loop always
overwrites **the lowest p-score in the histogram with `0`**, discarding the
positive `q` the AFDR walk just computed for it.

Confirmed against the compiled oracle: for the fixture `tiny/two_contigs` the
`--cutoff-analysis` table ends

```
0.60   0.09
0.30   0.00
```

and the histogram has exactly one p-score below 0.60. A hand-computed table
gives `q(0.30) = 0.30 + (log10(k) - log10(N)) > 0`, so the `0.00` is this
overwrite and not the `q <= 0` break.

**Consequence for the port:** the *lowest* p-score present in the genome
therefore always has q-value 0, which is harmless in practice (a q of 0 is the
best possible score and that p-score is by definition the least significant) but
is **not** what the AFDR formula says. `macs_score::PqTable` reproduces it, and
`pq_tables_lowest_bucket_is_forced_to_zero` pins it.

A more defensible reading of the comment would be `range(i + 1, ...)` when the
loop completed without a break, but the intent is not recoverable and the code is
the contract.

## F19 — track arrays are zero-padded until `finalize()`

`FixWidthTrack.add_loc` allocates each chromosome's per-strand numpy array at
`buffer_size` (100000) and tracks the fill count in a parallel pointer dict. The
arrays are resized down to that count only in `finalize()`. `PairedEndTrack` does
the same with a structured `[('l','i4'),('r','i4')]` array.

Consequence: reading `get_locations_by_chr()` straight after `build_fwtrack()`
returns 200000 entries per chromosome, all but the real ones zero. Any port that
skips `finalize()` — or any harness that forgets it — sees plausible-looking zeros
rather than an error. `oracle/dump_oracle_records.py` calls `finalize()` for this
reason, and `SingleEndStore::finalize` is the Rust counterpart.

## F20 — the FRAG barcode is parsed but not stored

`FragParser.pe_parse_line` returns `(chrom, left, right, barcode, count)` with the
barcode in **column 4** and the count in **column 5** — the reverse of the order
the prose in `MACS3/Commands/callpeak_cmd.py` implies, and the single easiest
thing to get wrong in this format. `build_petrack` keeps only the count
(`[('l','i4'),('r','i4'),('c','u2')]`) and discards the barcode string, so a
differential can verify coordinates and counts but has no oracle value for
barcodes. `macs-io` parses the barcode and unit-tests it; the differential compares
the other three fields.

The count is stored in an `unsigned short`. **It truncates; it does not clamp.**
Measured with `FragParser.build_petrack` on MACS3 3.0.5, not inferred:

| file column | stored |
|---|---|
| `70000` | `4464` (= 70000 & 0xFFFF) |
| `-5` | `65531` (two's complement) |
| `65535` | `65535` |

So the `except OverflowError: thiscount = 65535` branch in `pe_parse_line` never
fires -- Cython emits a plain truncating C cast rather than a range check -- and
upstream's own "The count in this line is over 65535, and will be capped at
65535" warning is **dead code**. Clamping here, which is the obvious reading of
the source, would disagree with upstream on every count above 65535.

Pinned by the `frag_counts/out_of_range` fixture, whose counts include 70000, -5,
0, 65535 and 65536, and covered by `crates/macs-io/src/formats.rs`. See F22.

## F21 — `atoi`/`atof` are lenient, and `split()` collapses whitespace but `split(b"\t")` does not

Upstream calls libc `atoi`/`atof` (`Parser.py:26`, `BedGraphIO.py:28`). Unparseable
input yields `0`, trailing garbage is ignored, and overflow saturates at the `int`
extremes. bedGraph uses Python's argument-less `split()`, which collapses runs of
whitespace; BED/BEDPE/FRAG use `split(b'\t')`, which emits an empty field between
adjacent tabs. Rust's `split` with a closure predicate behaves like the latter, so
the whitespace-collapsing case has to be written by hand — see `split_line`.

Separately, `BEDParser.fw_parse_line` reads the **strand column to choose the
coordinate**: `+` reads take `fields[1]` (the start) and `-` reads take
`fields[2]` (the *end*, one base past the last covered base). A six-column BED is
therefore two 5' ends at different coordinates, not one fragment.

## F22 — the FRAG "capped at 65535" branch is dead code

`FragParser.pe_parse_line` wraps the count assignment in
`except OverflowError: thiscount = 65535` and warns that the count "will be capped
at 65535". That branch is unreachable: Cython assigns the C `int` result of `atoi`
to a `cython.ushort` with a truncating cast, not a checked conversion, so no
`OverflowError` is raised for either an over-large or a negative value.

Consequences for a port:

* a count of 70000 is stored as 4464, not 65535;
* a count of -5 is stored as 65531, not 0;
* the only correct transcription is a two's-complement truncation to 16 bits.

This is the one place in the parsers where reading the source and reading the
binary disagree. It was settled by measurement (`oracle/dump_oracle_records.py`
against a fixture carrying the offending counts), not by reading the code, and it
would have been a silent disagreement -- `macs-io` clamped on the first attempt
and matched on all 52 existing fixtures because none of them had a count above 7.

## F23 — MACS3's help text documents flags that do not exist

`macs3 callpeak --help` contains `--exsize`, `--extension` and `--shiftsize`;
`macs3 randsample --help` contains `--num` and `--percent`. None of these are
option strings — they appear only inside `help=` prose. The real ones are
`-e/--extsize`, `-s/--shift`, `-n/--number`, `-p/--percentage`.

Consequence for tooling: the CLI surface must be derived by walking argparse's
action tree (`bin/macs3:prepare_argparser`), not by scraping `--help`. A scraper
would add flags upstream rejects, and the port would accept invocations that must
exit 2. `oracle/gen_flag_matrix.py` walks the tree; `oracle/flag_matrix.tsv` is
the generated 315-row result (14 subcommands).

## F24 — the score engine's state is invisible from Python

`PeakDetect` does not use `ScoreTrack` at all: `self.scoretrack` is initialised
to `None` and never assigned. The real engine is `CallerFromAlignments`
(`MACS3/Signal/CallPeakUnit.py`), and it is a Cython `cdef` class whose
`pqtable`, `d`, `ctrl_d_s`, `ctrl_scaling_factor_s`, `lambda_bg` and
`pvalue_length`/`pvalue_npeaks` are all C attributes that raise on access from
Python. `PeakContent` is worse: it exposes no Python attributes at all and is not
iterable, so peak rows cannot be read back.

Both facts push stage capture the same way — intercept the module globals that
`callpeak_cmd` and `PeakDetect` resolve at call time (`load_tag_files_options`,
`PeakDetect`, `CallerFromAlignments`), and read anything genuinely private off the
files upstream itself writes. `oracle/dump_stages.py` does this. The per-position
**per-scale** lambda arrays remain unreachable: they are locals inside
`CallerFromAlignments.call_peaks`, and capturing them would need either an
upstream patch or a faithful re-derivation. They are the one item on the G0
intermediate list still outstanding.

## F25 — control is duplicate-filtered with *treatment's* limit

`callpeak_cmd.py` computes `control_max_dup_tags` from the control's own tag count
when `--keep-dup auto`, logs it, and then calls

```python
control.filter_dup(treatment_max_dup_tags)
```

using the **treatment** value. `control_max_dup_tags` is dead. So under
`--keep-dup auto` the control's duplicate limit is whichever value the treatment
computed, not its own. This is invisible whenever both samples are the same size,
and only shows up as an asymmetry in control retention between an auto-threshold
run and a manually pinned run.

Reproduce with a treatment far larger than its control and `--keep-dup auto`; the
control's retained count should follow the treatment's binomial threshold. Port
must reproduce the bug, not the intent. Not yet empirically confirmed — flagged
for a fixture that separates the two sample sizes.

## F26 — the fixture generator was not reproducible (our bug, not upstream's)

Not an upstream finding, but it invalidates a golden corpus if unfixed, so it is
recorded next to the ones that are.

Two fixtures were seeded from `hash()` of a value containing a `str`:

```python
rng = random.Random(hash(shape) & 0xFFFF)                       # se_shapes
rng = random.Random(abs(hash((name, label, chrom, i))) & 0xFFFF) # strand choice
```

CPython randomises string hashing per process unless `PYTHONHASHSEED` is set, so
**every regeneration produced different BED strand assignments**. The corpus looked
stable only because nobody diffed two generations against each other. Replaced with
`stable_seed()`, a SHA-256 digest of the parts.

The lesson generalises: a golden corpus is only evidence if regenerating it is a
no-op, so `gen_fixtures.py` must be verified by generating twice and by generating
under two different `PYTHONHASHSEED` values. Both are now clean across 425
fixtures.

## F27 — no-common-chromosome reports an internal error instead of a diagnosis

`callpeak_cmd.check_names` is meant to explain *why* treatment and control share no
chromosome:

```python
error_stream("No common chromosome names can be found from treatment and control!")
error_stream("Please make sure that the treatment and control alignment files "
             "were generated by using the same genome assembly!")
error_stream("Chromosome names in treatment: %s" % ",".join(sorted(tchrnames)))
error_stream("Chromosome names in control: %s" % ",".join(sorted(cchrnames)))
sys.exit()
```

`get_chr_names()` returns a set of **`bytes`**, so the last two lines raise
`TypeError: sequence item 0: expected str instance, bytes found`. The user never
sees the diagnostic; they see a traceback and exit code 1. The first two lines do
print, so the useful part appears and then the crash follows it.

Measured: 38 of the corpus's runs hit this (any fixture whose treatment and control
share no chromosome, e.g. `se_edge/one_sided_no_reads`).

The port must reproduce exit code 1 and the two printed lines. Whether to also
reproduce the traceback is a judgement call — **it is the only place where upstream
emits a Python traceback on a normal, non-malformed input**, and the objective's
"zero panics on malformed input" clause clearly does not ask for one. Recorded here
as an explicit compatibility decision to confirm with the user.

## F28 — ZeroDivisionError when control has reads but no usable control pileup

`PeakDetect.__call_peaks_w_control` computes

```python
control_sum = self.control.total * self.d
self.ratio_treat2control = float(treat_sum)/control_sum
```

With `--shift -100` on `se_edge/contig_edges`, the control's reads are all shifted
past the contig boundary and the control pileup ends up empty, so `control_sum` is
0 and the division raises. Exit code 1, no output files.

18 corpus runs hit this. Unlike F27 this one is a genuine unhandled division rather
than a formatting bug, so the faithful behaviour is simply "exit 1, create nothing".
The port should detect the zero-control case explicitly rather than rely on a
floating-point division to produce a non-finite value, since the two diverge on
whether a partial output is written first.

## F29 — `get_chr_names()` returns an unordered `set`, so chromosome processing order is undefined

`FixWidthTrack.get_chr_names` is annotated `-> set` and returns

```python
return set(sorted(self.locations.keys()))
```

The `sorted()` is thrown away by `set()`. Every consumer then iterates that set
directly — `sort`, `filter_dup`, `sample_percent`, `sample_num` — so the order in
which chromosomes are visited, and therefore the order in which each chromosome's
array consumes draws from the shared NumPy RNG stream, depends on CPython's set
iteration order, which for `bytes` keys is **randomised per process by
`PYTHONHASHSEED`**.

Measured on five chromosomes, three runs, same input:

```
PYTHONHASHSEED=1  ->  chrA, chr2, chrB, chrC, chr1
PYTHONHASHSEED=2  ->  chr2, chr1, chrB, chrC, chrA
PYTHONHASHSEED=3  ->  chr1, chrC, chrB, chr2, chrA
```

Consequences:

* `filter_dup` is unaffected, because each chromosome's result is independent and
  `total` is a sum. The 420 duplicate-filtering vectors are stable, as measured.
* `sample_percent` / `sample_num` **are** affected: chromosome visit order changes
  which array is shuffled first, so a different *subset* survives. The retained
  count stays correct; the retained positions do not.
* Therefore upstream's own downsampled output is **not reproducible across
  processes** for inputs with more than one chromosome.

The port resolves this by iterating chromosomes in **sorted byte order**, which is
the order the author evidently intended and is the only stable choice. The retained
*count* then matches upstream for every input, and the retained *set* matches
whenever upstream happens to run under a hash seed whose set order is sorted.

This is a genuine incompatibility that cannot be closed by better code on this
side: two runs of `macs3 --down-sample` can legitimately disagree with each other.
Recorded so it is not mistaken for a port bug. If exact subset parity is required
for downsampling, the corpus must either pin `PYTHONHASHSEED` when generating
goldens or restrict sampling cases to single-chromosome inputs; both are worse than
declaring the sorted-order choice, which is what upstream's `sorted()` was reaching
for.

## F30 — `sample_num` above the track size zero-pads instead of clamping

`FWTrack.sample_num` converts a target count into a fraction and delegates:

```python
percent = cython.cast(cython.float, samplesize) / self.total
self.sample_percent(percent, seed)
```

and `sample_percent` does

```python
num = cython.cast(cython.int, round(self.locations[k][0].shape[0] * percent, 5))
np.random.shuffle(self.locations[k][0])
self.locations[k][0].resize(num, refcheck=False)
```

When `samplesize > total`, `percent > 1`, so `num` exceeds the array length.
NumPy's `ndarray.resize` to a *larger* size does not clamp — it **zero-fills**:

```
before total: 37 ; sample_num(1000) -> after total: 999
plus: [0, 0, 0, 0, 0, 0, 0, 0] ...
minus: [0, 0, 0, 0, 0, 0, 0, 0] ...
numpy.resize([1,2,3] -> 6) = [1 2 3 0 0 0]
```

So asking for more reads than the track holds yields a track of fabricated
positions at 0, and `total` reports `samplesize`-ish rather than the real count.
These zero positions are then pileup'd as if they were reads.

Measured directly on MACS3 3.0.5. The port reproduces it: `sample_num` past the end
zero-pads rather than saturating. `percent == 1.0` and `percent > 1.0` are therefore
different behaviours, and a "sensible" clamp would silently change both the retained
count and every downstream coordinate.

## F31 — `sample_num`'s percentage is a C `float`, and the f32 rounding is observable

`FWTrack.sample_num` declares

```python
def sample_num(self, samplesize: cython.ulong, seed: cython.int = -1):
    percent: cython.float
    percent = cython.cast(cython.float, samplesize) / self.total
    self.sample_percent(percent, seed)
```

`cython.float` is a C `float`, so the quotient is computed in f64 and then
**truncated to f32** before `sample_percent` uses it. Because `sample_percent` then
computes `int(round(len * percent, 5))` — a *truncating* cast (see
`macs-track`'s `retained_count`) — the f32 error pushes the product just under an
integer and the result loses one.

Measured on MACS3 3.0.5 with 9 plus and 6 minus positions and `sample_num(1000)`:

```
f64: 1000/15            = 66.66666666666667 -> 600.0, 400.0 -> 600 + 400 = 1000
f32: float32(1000/15)   = 66.66666412353516
     9 * that           = 599.9999771118164 -> trunc 599
     6 * that           = 399.99998474121094 -> trunc 399   -> total 998
upstream reports: 998
```

Two arithmetic details had to be right to land on 998:

1. the quotient passes through `f32`;
2. `round(x, 5)` is a round-half-to-even *decimal* rounding and the subsequent
   `cython.int` cast **truncates** — it is not a rounding step.

Getting (1) wrong gives 1000. Getting (2) wrong as a rounding step gives 600+400 as
well, but for the wrong reason: on `sample_percent`'s own tie cases, `0.5` must
become 0.

## F32 emulation sites so far

Collected because the objective asks for a documented list. Each is a place where
upstream's C types make an f64 computation observable in the output:

| site | C type | effect |
|---|---|---|
| `FWTrack.sample_num` percent | `cython.float` | F31, one position lost per strand |
| `ScoreTrack.add_chromosome` score column | `dtype="float32"` | p/q scores are f32 |
| `PileupV2` return arrays | `cython.float` | pileup values are f32 |
| `PeakModel` fragment values | `cython.float` | model `d` is f32 |
| `CallPeakUnit` p/q arrays | `dtype="f4"` | local lambda and p/q are f32 |
| `fix_nomodel_shift` `end_shift` | `cython.float` | shift arithmetic is f32 |

## F32 — `maxima` reports the index of the sign *fall*, and a symmetric spike gives two indices

`SignalProcessing.maxima` smooths the signal, takes the sign of the derivative, and
keeps where that sign decreases:

```text
smoothed = savitzky_golay_order2_deriv1(signal, window_size).round(16)
sign     = np.sign(smoothed)
diff     = np.diff(sign)
m        = np.where(diff <= -1)[0]
```

`np.diff` yields `sign[i+1] - sign[i]`, so `diff <= -1` is a **fall** in the sign. A
symmetric spike has its derivative rising, crossing zero, then falling, so it
contributes *two* adjacent indices — one at the zero crossing and one after it.

Measured on MACS3 3.0.5:

```
single spike at 150, window 51 -> maxima: [149, 150]
spikes at 120 and 280          -> maxima: [119, 120, 279, 280]
```

The natural but wrong transcription is `sign[i] - sign[i+1] <= -1`, which mirrors
the condition and reports the *rising* edges. On the same inputs that yields
`[124, 175]` and `[94, 145, 254, 305]` — the padding boundaries, not the peaks.
Both summits then map to no above-cutoff chunk and are dropped, so `--call-summits`
silently produces no summits at all rather than an obviously wrong answer.

`enforce_peakyness` then removes candidates that are narrower than 50 bp, have
fewer than 6 distinct values, or sit in the 10 bp padding. Note the filter does
*not* narrow a spike: with a single maximum there are no interior minima, so
`enforce_peakyness` returns its input unchanged. Narrowness is rejected by
`maxima` and by the caller, not by this filter.

## F33 — the Savitzky-Golay edge padding is asymmetric, and getting it wrong shifts right-edge maxima

Upstream pads with

```python
firstvals = signal[0] - np.abs(signal[1:half_window+1][::-1] - signal[0])
lastvals  = signal[-1] + np.abs(signal[-half_window-1:-1][::-1] - signal[-1])
```

The two slices are mirrored but walked in **opposite directions**: `signal[1:half+1]`
reversed runs `half, half-1, ..., 1`, while `signal[-half-1:-1]` reversed runs
`n-2, n-3, ..., n-half-1`. So the first value appended on the right uses
`signal[n-2]`, not `signal[n-half-1]`.

Writing both loops in the same direction mirrors the right padding and changes the
derivative in the last `half_window` samples. On a 200-sample ramp with
`window_size = 51` that is a maximum absolute difference of **0.53** at index 184,
which is exactly the width of a summit that sits near a contig end.

The coefficients themselves are `np.linalg.pinv([[1, k, k**2]])[1]`, which for a
quadratic fit equals the closed form `3k / (h(h+1)(2h+1))` over offsets `-h..=h`.
Upstream passes `m[::-1]` to `np.convolve`, which flips internally, so the two
reversals cancel and the coefficients are applied in ascending-offset order.
Getting that backwards flips the derivative's sign, which moves every detected
maximum to the wrong side of the peak.

## F34 — `--nolambda` spans the whole treatment, not a single point

`pileup_treat_ctrl_a_chromosome` substitutes the control pileup when local lambda is
disabled:

```python
if self.no_lambda_flag:
    ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg,], dtype="f4")]
```

That is a **one-entry** control array anchored at the treatment's *last*
breakpoint, carrying only the global lambda. It is tempting to read that as "the
paired output is one entry", because the merge then sees a length-1 input.

It is not. `__chrom_pair_treat_ctrl` advances on whichever end is smaller:

```python
while it < lt and ic < lc:
    ...
    if t_p_ptr[0] < c_p_ptr[0]:   # advance treatment, control head is carried
```

The control cursor only moves when it hits that final anchor position, by which
point the treatment cursor has already walked the whole chromosome. So for a
treatment of `n` runs the paired array has `n` entries and the control value is
`lambda_bg` throughout.

Verified against MACS3 3.0.5: `callpeak --nolambda --bdg` on the
`se_model/realistic` fixture writes **376** pileup intervals spanning
`chr20 0..93003` — the full treatment extent, not one point.

The port implements the merge, not the shortcut. A "helpful" reading that returned
a single paired entry would produce one peak candidate for the entire chromosome.

## F35 — both the pairing and the local-lambda merge drop the longer tail

Two separate cursor merges, both `while i_a < l_a and i_b < l_b` with no tail
drain:

* `__chrom_pair_treat_ctrl` pairs the treatment and control pileups;
* `over_two_pv_array` combines the per-scale local lambda tracks.

In both, when one array is exhausted the remainder of the other is **discarded**.
So the paired arrays are shorter than either input, and combined lambda stops
where the *shorter* scale stopped.

Because a wider window produces more breakpoints, the `llocal` track usually has
the most entries and therefore usually does *not* truncate; but when the control
track is sparse on a chromosome, the `d`-scale track can be the shorter one and
the merged lambda — and hence every peak boundary on that chromosome — stops early.

Reproduced rather than "fixed" by draining the longer array. Whether upstream
intends this is unknowable from the source, but the golden corpora are the
authority and they record the truncation.

Measured on synthetic cases in `oracle/check_pairing.py`: treatment at
`[10, 20, 30]` against control at `[10]` pairs to exactly **one** entry, and the
`[20, 30]` treatment runs are gone.

## F36 — the cutoff comparison promotes the f32 score array to f64, and that changes the answer

`__pre_computes` (and `__chrom_call_peak_using_certain_criteria`) threshold with

```python
above_cutoff = np.nonzero(score_array > cutoff)[0]
```

`score_array` is `f4`. `cutoff` comes from `tmplist`, which is built as

```python
tmplist = [round(x, 5) for x in sorted(list(np.arange(0.3, 10.0, 0.3)), reverse=True)]
```

and `np.arange` yields `numpy.float64`, so every ladder element is an f64 **scalar**,
not a Python float.

NumPy's promotion rule depends on exactly that distinction:

| comparison | result | why |
|---|---|---|
| `f32 array > np.float64(0.3)` | **True** | an f64 *scalar* promotes the array to f64; `f64(float32(0.3))` is `0.30000001192092896` |
| `f32 array > 0.3` (Python float) | False | a Python float is a *weak* scalar and is cast down to f32, giving equality |
| `f32 array > 9.9` (`np.float64`) | False | `f64(float32(9.9))` is `9.899999618530273`, below `9.9` |

So a p-score written as `0.3f32` **does** clear the `0.3` cutoff, while one written
as `9.9f32` does **not** clear the `9.9` cutoff. Both look like "equal" values and
behave differently.

The port compares `f64::from(score) > cutoff`, widening the score, which is what
upstream does. The tempting alternative — narrowing the ladder to `f32` so the
comparison stays in one lane — silently drops the low cutoffs, and since the
cutoff ladder is the *only* place those values are used, the effect is confined to
`--cutoff-analysis` runs. That makes it very easy to believe the tests pass.

Note this is distinct from the `pvalue`/`qvalue` *scores*, which genuinely are f32
upstream (see the f32 site table). Only the cutoff ladder is f64.

## The `--cutoff-analysis` ladder is a fixed 33-value list

`sorted(np.arange(0.3, 10.0, 0.3), reverse=True)` rounded to 5 decimals: 9.9 down
to 0.3 in 0.3 steps. It does not depend on the data, so every input walks the same
cutoffs. Real `*_cutoff_analysis.txt` files have fewer than 33 rows because rows
with `npeaks == 0` are skipped when writing.

## F37 — `enforce_peakyness` is what rejects `--call-summits` output on noisy backgrounds

Not a port bug — upstream behaves the same way — but it is a surprising
interaction worth writing down, because it makes `--call-summits` look broken on
signals that are perfectly reasonable to generate.

`--call-summits` finds maxima by smoothing the pileup with a Savitzky-Golay window
of `min_length` (i.e. the fragment size `d`) and taking sign changes of the first
derivative. A **smooth** signal therefore behaves as expected: measured on two broad
triangular humps over a rising background, `maxima` returns 3 candidates and
`enforce_peakyness` keeps 2, giving 2 summits.

A **sawtooth** background of the same amplitude does not work at all:

```
smooth background:      maxima -> 3 candidates, enforce_peakyness keeps 2
sawtooth background:    maxima -> 67 candidates, enforce_peakyness keeps 0
```

The reason is the interaction between the two filters. A sawtooth makes the smoothed
derivative cross zero every few samples, so `maxima` returns a candidate every ~50
bases. `enforce_peakyness` then computes, for each candidate, the valley between it
and the *adjacent* minima — which are only a few bases away — so every candidate is
bounded by a near-degenerate interval and fails `is_valid_peak`'s 50 bp width
requirement. All 67 are discarded and `--call-summits` emits no summits at all.

Two independent thresholds therefore have to be satisfied for a summit to survive:
the smoothed derivative must have a genuine sign change, **and** the region between
the surrounding minima must be at least 50 bp wide with at least 6 distinct values.
Real ChIP/ATAC pileup profiles are smooth, so this only bites on synthetic fixtures
— but it will silently produce zero summits rather than an error.

## F38 — the zero-padding of F19 becomes a *real read* when `filter_dup` runs before `finalize`

F19 established that `FWTrack.build_fwtrack` returns arrays pre-sized to
`buffer_size`, zero-filled, with `total` left at 0. Measured on
`tests/fixtures/se_shapes/bias/treat.bed`:

```
build_fwtrack : plus len=100000  minus len=100000  total=0
                plus distinct=177  minus distinct=2      <- 0 and 200
after finalize: plus len=3526   minus len=3474   total=7000
                minus values=[200]                     <- the 0 was padding
```

So before `finalize`, position `0` appears as a real value on both strands — but it
is padding. `filter_dup` sorts the array, which puts every padded zero at the front
as one enormous run at position 0, and then truncates runs to `maxnum`. That run
collapses to a single entry, so **deduplication turns the padding into one genuine
read at position 0 on each strand**.

Calling `filter_dup(1)` directly on a freshly built (unfinalized) track therefore
reports **179** retained reads where the same track, finalized first, reports
**177** — and upstream's own `callpeak` log prints 177.

Consequences for a port:

* this crate never zero-pads, so the hazard cannot arise here — but it means a
  differential that drives `filter_dup` on a raw `build_fwtrack` result is testing
  something upstream's real pipeline never does, and will disagree by exactly one
  read per strand;
* any future BAM/SAM path that sizes its arrays up front **must** trim before
  deduplicating, or it will silently gain a phantom read at 0;
* the phantom read sits at position 0 and is extended `d` bases 3'-ward, so it
  contributes a constant `+1` to the treatment pileup over `[0, d)`.

Verified with `MACS3.Signal.FixWidthTrack` / `MACS3.IO.Parser.BEDParser` on
MACS3 3.0.5, and cross-checked against the `tags after filtering` line of the
recorded golden run's stdout.

## F39 — `tocontrol` scales *treatment* down, so it is false when the control is larger

`PeakDetect.__call_peaks_w_control`:

```python
if self.opt.tocontrol:
    lambda_bg   = float(control_sum) / self.gsize
    treat_scale = 1 / self.ratio_treat2control
else:
    lambda_bg   = float(treat_sum) / self.gsize
    treat_scale = 1.0
```

The name is a trap: `tocontrol` means "scale **to** the control's depth", which
happens only when the *treatment* is the larger sample. `callpeak_cmd` with the
default `--to-small` sets it from

```python
if t1 > c1:  options.tocontrol = True    # treatment is bigger -> scale treatment down
else:        options.tocontrol = False   # control is bigger    -> scale control down
```

So for a control-dominant library `tocontrol` is `False`, `treat_scale` is `1.0`, and
`lambda_bg` comes from `treat_sum`. Reading the flag backwards scales the wrong
sample *and* takes lambda_bg from the wrong total, which perturbs every control
lambda value and therefore every peak boundary.

The three control scaling factors are then, with `ratio = treat_sum / control_sum`:

| scale | factor when `not tocontrol` | when `tocontrol` |
|---|---|---|
| `d` | `ratio` | `1.0` |
| `slocal` | `d / slocal * ratio` | `d / slocal` |
| `llocal` (only if `llocal > slocal`) | `d / llocal * ratio` | `d / llocal` |

Note the local-window factors are the region size **into** `d`, and the `llocal`
scale is omitted entirely unless `llocal > slocal` — it is not clamped or merged.

## F43 — the q-score track is `-log10(q)`, and chunk starts are the previous end **plus one**

Wiring the first end-to-end peak comparison (`macs-callpeak-e2e`) turned up two
conventions that are easy to get backwards.

**`-log10(q)`, not `q`.** `qscore_track` feeds `-log10(q)` into the track, so
significance is `-log10(q) >= -log10(qvalue)`. Comparing the track value against
`qvalue` directly is a category error -- with the default `qvalue = 0.05` the
threshold is `1.301`, and every non-zero track value clears `0.05`, so the whole
chromosome becomes one significant region.

**Chunk start is `previous run end + 1`.** The signal arrays are end-indexed
(F5): a run ending at `e` covers through `e` inclusive, so the following run
begins at `e + 1`. Taking the previous end as the start places every peak exactly
one base too early. On `se_model/realistic` this was the *only* start/end error:
all 12 peaks then had byte-identical `start` and `end`.

## Open: summit coordinates, and one merge boundary

With starts and ends fixed, `se_model/realistic` gives **12 peaks against
upstream's 12**, with `start`/`end` identical on 11 of 12. Two real gaps remain:

* **Summits are consistently early.** Every summit this port reports is lower
  than upstream's, by 18 to 66 bases (7709 vs 7732, 15488 vs 15506, 23037 vs
  23056, 30945 vs 30992, 46168 vs 46195). A one-sided, always-negative offset
  points at the candidate/midpoint stage of summit selection rather than at the
  signal arrays, since the signal arrays are now confirmed identical -- if they
  were not, the *ends* would move too. The `lower median of candidate midpoints`
  tie-break (see the `macs-peaks` crate docs) is the prime suspect.
* **One merge boundary.** Peak 4 starts at 38105 where upstream starts at 38462,
  357 bases earlier. The peak to either side agrees, so this is a single
  `max_gap` merge being applied where upstream does not merge, not a systematic
  offset.

Both are recorded as open rather than guessed at. Summit coordinates are an
explicit acceptance criterion, so this is the next thing to close.

## F44 — chunks are contiguous **position-level** above-cutoff spans

Upstream builds `peak_content` from the position-level arrays, not from runs:

```python
for i in range(len(dq_track)):
    q = -1.0 * dq_track[i]
    if q >= self.pqtable.cutoff and d_pileup_d[i] > 0:
        if new_peak:
            new_peak = False
            chunks.append([i, i, d_pileup_d[i], ctrl_pileup_d[i]])
        else:
            chunks[-1][1] = i
            chunks[-1][2] = d_pileup_d[i]   # overwritten every position
            chunks[-1][3] = ctrl_pileup_d[i]
    else:
        new_peak = True
```

Three things follow, and the port now implements all three:

* one chunk spans a whole contiguous above-cutoff stretch;
* `treat`/`ctrl` are the values at the chunk's **last** position, continuously
  overwritten -- not a per-run value;
* `d_pileup_d[i] > 0` is a filter in its own right, so a zero-pileup position
  ends the chunk even when its q-value clears the cutoff.

## F45 — regions are found by `np.nonzero(score > cutoff)`, and start is the *preceding* position

Upstream does not walk a q track to segment regions. In
`__chrom_call_peak_using_certain_criteria`:

```python
score_array_s = [self.__cal_qscore(treat_array, ctrl_array)]   # 'q' by default
above_cutoff = np.nonzero(apply_multiple_cutoffs(score_array_s, score_cutoff_s))[0]
above_cutoff_index_array = np.arange(pos_array.shape[0], dtype="i4")[above_cutoff]
above_cutoff_endpos      = pos_array[above_cutoff]
above_cutoff_startpos    = pos_array[above_cutoff - 1]
```

Three details, all of which this port now matches:

* the test is a **strict `>`** on the `-log10(q)` track, not `>=`;
* `above_cutoff_startpos` is `pos_array[above_cutoff - 1]` -- the position
  *before* the first above-cutoff position. Combined with F5 this is where the
  chunk `start` comes from, and it is why `previous run end + 1` is the right
  construction;
* `apply_multiple_cutoffs` *counts* how many scores beat their cutoffs, so with
  the default single `q` score the condition is exactly `qscore > cutoff`.

`__cal_qscore` also reveals a numerical contract worth stating explicitly:

```python
s_ptr[0] = self.pqtable[get_pscore(cython.cast(cython.int, a1_ptr[0]), a2_ptr[0])]
```

the **treatment value is cast to `int`** before the p-score is computed, and the
table stores `-log10(q)` -- the worst p-scores map to 0 and the best to the
largest value, which is why `score > cutoff` is the correct orientation.

With F43, F44 and F45 in place, `se_model/realistic` reaches **12 peaks against
upstream's 12**, and `oracle/run_peak_e2e.sh` measures **10 of 12** with
byte-identical `start` and `end` (worst disagreement 357 bp, the single `max_gap`
merge boundary). An earlier claim in this file that all 12 matched was wrong: it
came from reading a filtered subset of the diff output, not the whole table.

## Open: summit coordinates only

The sole remaining difference on `se_model/realistic` is `abs_summit`, and it is
one-sided: this port's summit is always **lower** (earlier) than upstream's, by
8 to 95 bases (61722 vs 61743, 69317 vs 69412, 76895 vs 76903, 84597 vs 84670).

`__close_peak_wo_subpeaks` has been compared line by line and matches, including
the lower-median tie-break:

```python
if not summit_value or summit_value < tscore:
    tsummit = [(tend + tstart) // 2,]; tsummit_index = [i,]; summit_value = tscore
elif summit_value == tscore:
    tsummit.append((tend + tstart) // 2); tsummit_index.append(i)
midindex = (len(tsummit) + 1) // 2 - 1
```

A strictly-lower result means the *set* of tied maxima differs, not that the
selection rule does: the summit is `(tend + tstart) // 2` of the chunk with the
largest `ttreat_p`, and `ttreat_p` is the treatment value at that chunk's last
position. Since peak start/end and the p->q tables now agree exactly, the
remaining suspect is chunk *granularity* inside a peak -- specifically that
upstream's `ttreat_p` comes from `d_pileup_d[pos_array[i]]` at the last
above-cutoff **index**, while this port takes the treatment value at the chunk's
last **run end**. Those coincide only when the above-cutoff run ends on a
treatment-pileup boundary.

Next step: compare, per peak, the chunk list and each chunk's `ttreat_p` between
`__chrom_call_peak_using_certain_criteria` and the port, which needs the region
contents rather than just the XLS.

## F46 — the p-score is a Poisson tail probability, not `treat - ctrl`

Upstream's `ScoreTrack` does not subtract the two tracks. Per base it evaluates

```python
v[i] = get_pscore(cython.cast(cython.int, (p[i] + self.pseudocount)),
                  c[i] + self.pseudocount)
```

and `get_pscore` is `-1 * poisson_cdf(k, lam, ...)` -- a Poisson survival
probability, negated. The treatment side is **truncated to an integer after the
f32 add**; the lambda side stays `f32`. `macs-stats/src/poisson.rs` and
`macs-score/src/pscore.rs` already implement exactly this, and
`macs_score::pscore_track(chrom, treat, lambda, cache, pseudocount)` computes the
whole track in one pass.

So the harness in `macs-callpeak-e2e` was using the **wrong score**: a pointwise
subtraction (`over_two_pv_array_track`). The two agree closely enough that peak
boundaries matched on `se_model/realistic` (12/12), which is why it went unnoticed,
but they are not the same function and they must be.

Switching to `pscore_track` regressed the boundaries, which localises the
remaining error:

```text
subtraction   peak 7: 61430-61883   (matches upstream exactly)
pscore_track peak 7: 61453-61820   (~150 bp narrower on both flanks)
```

Narrower on *both* flanks, symmetrically, is the signature of a score that is
systematically too pessimistic rather than shifted or mis-scaled. Prime suspects,
in order: the pseudocount convention (`ScoreTrack` applies it to both terms
before the Poisson call, and `pseudocounted_inputs` does too, so this is more
likely the *orientation* of `get_pscore`), the lambda actually passed in
(`pscore_track` intersects the two tracks' supports, and the control lambda
carries a `lambda_bg` baseline that may not belong in the p-score), or the
int-truncation of the treatment term.

The subtraction is therefore **kept in the harness for now**, because it
reproduces upstream's peak boundaries exactly on the reference fixture, with an
explicit `NOTE(F46)` at the call site so the substitution is not mistaken for the
intended design. This is a known-wrong mechanism with a known-correct
replacement, and closing F46 properly is a prerequisite for summit parity,
because the summit depends on `ttreat_p` comparisons that a slightly different
score track perturbs.

### F46a — `get_pscore` is the *upper* tail in log10 space

```python
def poisson_cdf(n: uint32_t, lam: float64_t, lower: bool = False, log10: bool = False)
# lower=False -> upper tail CDF;  log10=True -> computed in log space
val = -1.0 * poisson_cdf(x, y, False, True)
```

So `get_pscore(k, lam)` is `-log10 P(X >= k)` with `X ~ Poisson(lam)` -- a
positive, larger-is-more-significant `-log10` p-value, which is the orientation
the q table and the `score > cutoff` test (F45) both rely on. Recording it
explicitly because the two extra arguments are positional booleans and the
`lower`/`log10` order is easy to transpose, which would silently return
`log10 P(X <= k)` -- a *negative* score, and therefore no peaks above a positive
cutoff at all.

### F46b — narrowing is symmetric, which points away from the treatment term

`pscore_track` replaced the subtraction and the peaks came out narrower on *both*
flanks by a similar amount. Since the treatment-side int truncation can only move
positions near a pileup step, and since the intersection of the two track
supports is what `pscore_track` walks, the two things that can shrink a peak
symmetrically are (a) the lambda actually supplied being systematically too high,
or (b) the p-score track being shorter than the treatment track so the outer
flanks are simply not scored.

**(b) is ruled out by measurement.** On `se_model/realistic` the tracks are:

```text
treat = 0..9223372036854775807   (clipped at zero, otherwise unbounded)
ctrl  = 0..100100
subtraction p-score track = 0..93003
```

`pscore_track` walks the *intersection* of the two supports, so its extent would
be `0..100100` -- **longer** than the subtraction track's `0..93003`. A track that
reaches further cannot make peaks narrower on both flanks, so the symmetric
narrowing is not a support problem and the cause lies in the score values
themselves: this port's Poisson evaluation yields systematically smaller
`-log10 P(X >= k)` than upstream's at the flanks, where the control lambda is
small. The next check is to compare `pscore` against upstream's `poisson_cdf` at
the specific flank positions, rather than comparing whole tracks.

## F47 — `pscore` is bit-exact; the F46 error is in track assembly, not the math

The obvious suspect for F46 was `macs_score::pscore`. It is not wrong.
`crates/macs-score/tests/pscore_vs_upstream.rs` checks it against a grid
generated straight from the oracle (`oracle/gen_pscore_grid.py`,
`crates/macs-score/tests/pscore_grid.tsv`):

```text
117 vectors, worst abs diff 0e0 at k=0 lam=0
```

Zero, not merely within a tolerance: `pscore(k, lam)` reproduces
`-1 * poisson_cdf(k, lam, lower=False, log10=True)` exactly over
`k` in 1..40 and `lam` in 0.5..20, which covers the whole operating range of the
reference fixture (treat_scale is 1.0 there, so `k` is a small integer and
`lam` is a small fraction).

So the Poisson evaluation, the upper-tail orientation, the log10 space and the
f32 rounding are all correct. F46's narrowing must come from *what is fed in* or
*how the track is assembled*. The remaining candidate is now specific: the merged
control lambda.

`over_two_pv_array_track` walks the **intersection** of the two inputs, so
merging per-scale control lambdas of different extents silently truncates to the
shortest one. Upstream merges with `d_pileup_track_from_control_or_treat`, which
keeps the longer track's tail where the other has nothing. A control lambda that
is too short leaves the outer flanks unpaired, and since `pscore_track` walks the
intersection of treatment and lambda, that is precisely the mechanism that makes
peaks narrower on **both** flanks.

Measured extents are consistent with this:

```text
treat              = 0..9223372036854775807
merged control     = 0..100100          <-- ends well short of the treatment
subtraction p-score= 0..93003           <-- intersection with the shorter input
```

The next check is to confirm that the merged control loses tail that
`pointwise_max` should have kept -- comparing each per-scale track's extent
before and after the merge. That is a concrete, mechanical test rather than a
hypothesis about score magnitudes.

### F47a — correction: upstream's control merge *also* truncates to the intersection

F47 proposed that `over_two_pv_array_track` walking the intersection of its
inputs was the defect. It is not: upstream does exactly the same.

```python
# MACS3/Signal/Pileup.py:835
while i1 < l1 and i2 < l2:
    ret_v_ptr[0] = f(a1_v_ptr[0], a2_v_ptr[0])
    ...
return [ret_pos, ret_v]
```

`FWTrack.pileup_a_chromosome_c` folds the per-scale control lambdas together with
that function and `func="max"`, so the merged control is, like this port's, the
pointwise max over the **intersection** of the scales' extents. The merged
control ending at `100100` rather than running to the treatment's end is upstream
behaviour, not a port defect.

So F47's mechanism is falsified. Combined with F47 (pscore bit-exact) that
leaves exactly one untested difference in the p-score path: the **inputs**.
`pscore_track` is handed `treat` and the merged control `lambda` directly, whereas
upstream's `ScoreTrack` is handed the **paired** arrays `d_pileup_d` /
`ctrl_d_pileup_d` from `__chrom_pair_treat_ctrl`, which are position/value pairs
truncated at the paired extent (`F35`). Everything downstream -- the Poisson
evaluation, the p->q table, the cutoff test -- is common to both. The next check
is therefore to feed `pscore_track` the paired arrays rather than the raw tracks,
which is a wiring change in the harness, not a numerical one.

### F47b — the score mechanism is **exonerated** for the summit error

Both mechanisms are now selectable in `macs-callpeak-e2e` via `CALLPEAK_SCORE`
(`poisson` uses `pscore_track` over the paired/truncated extent, per F47a;
`subtract` uses the pointwise subtraction). Results on `se_model/realistic`:

| mechanism | boundaries | summits |
|---|---|---|
| `subtract` | **all 12 exact** | 7709/7732, 15488/15506, 23037/23056, 30945/30992 -- always low |
| `poisson`  | too narrow (7609-8025 vs 7572-8039) | 7829/7732 (**high**), 15490/15506, 23037/23056, 30947/30992 |

Two things follow.

**Poisson is upstream's real mechanism** (F46), and the summit error *changes
sign* under it -- peak 0 goes from 23 bp early to 97 bp late. A sign flip on a
sensitivity like this is what a correct-but-slightly-mispared input should do,
so the Poisson path is the right target, not the subtraction.

**But it cannot be the cause of the summit error.** Peaks 2 and 3 produce
*byte-identical summits* under both mechanisms (23037 and 30945) and are still
19 and 47 bases low against upstream. Whatever is wrong there is independent of
how the score track is produced, which rules the score mechanism in or out as the
culprit for those peaks.

That redirects the investigation to the **chunk `ttreat_p`** assignment, i.e. the
one input to `__close_peak_wo_subpeaks` that has not yet been pinned:

```python
if not summit_value or summit_value < tscore:
    tsummit = [(tend + tstart) // 2,]; summit_value = tscore
```

`tscore` is `ttreat_p`, the treatment value at the chunk's last above-cutoff
position. Upstream reads `d_pileup_d[pos_array[i]]` at that index; this port
takes the treatment value at the chunk's last *run end*. These agree only when
the above-cutoff stretch happens to end on a treatment-pileup boundary -- which,
given the summit is always low, is evidently not the common case.

`CALLPEAK_SCORE` is left in place so both paths stay comparable while that is
pinned down; `subtract` remains the default because it is the only one that
reproduces peak boundaries exactly today, and is labelled known-wrong at the call
site (F46) so it cannot be mistaken for the intended design.

## F47c — `tscore` is the max treatment **within** a run, not the value at its end

Reading the chunk loop in `__chrom_call_peak_using_certain_criteria` settles F47b:

```python
for i in range(1, above_cutoff_startpos.shape[0]):
    ts = acs_ptr[0]; te = ace_ptr[0]; ti = acia_ptr[0]
    tp = treat_array_ptr[ti]; cp = ctrl_array_ptr[ti]
    tl = ts - lastp
    if tl <= max_gap:
        peak_content.append((ts, te, tp, cp, ti))
    else:
        self.__close_peak_wo_subpeaks(peak_content, peaks, ...)
        peak_content = [(ts, te, tp, cp, ti),]
    lastp = te
```

`peak_content` accumulates **one entry per above-cutoff position**. Several
consecutive entries share the same run's `(ts, te)` -- because `ts` and `te` come
from `above_cutoff_startpos`/`above_cutoff_endpos`, which are constant across a
run -- while `ti` advances position by position, so `tp` is that position's own
treatment value.

So `__close_peak_wo_subpeaks` sweeps `tscore` over *positions*, and because
tied entries share a midpoint, the argmax is the run containing the **highest
treatment value anywhere in the run**. This port was instead assigning each chunk
the treatment value at the run's *last* position, which is why every summit came
out low: the peak of the run is somewhere inside it, not at its tail.

Fixed by taking `max` of the treatment (and control) track over each run's
`[start, end]` span. Effect on `se_model/realistic`, `subtract` mechanism:

| peak | summit before | summit after | upstream |
|---|---|---|---|
| 0 | 7709 (-23) | **7726 (-6)** | 7732 |
| 1 | 15488 (-18) | 15488 (-18) | 15506 |
| 2 | 23037 (-19) | 23037 (-19) | 23056 |
| 3 | 30945 (-47) | **30962 (-30)** | 30992 |

Boundaries remain exact on 10 of 12. Two of four summits improved sharply,
which confirms the mechanism; the residual 6-30 bp is now small and no longer
uniformly tied to run geometry, so the remaining candidates are the
`max_gap` region split and the `int` truncation of the treatment term, both of
which need the region contents rather than the XLS to check.

Note the Poisson mechanism also improved here (peak 0: 7829 -> 7724), so the two
are converging on the same answer, which is further evidence both are close to
right.

## F48 — neither score mechanism is at parity across the corpus; the default is the honest one

F47b exonerated the score mechanism for the summit error using one fixture. That
does **not** generalise. Running `macs-callpeak-e2e` over four SE fixtures, peak
counts and mismatched-peak counts:

| fixture | `subtract` (rust/upstream) | `poisson` (rust/upstream) |
|---|---|---|
| `se_model/realistic` | 12/12 | 12/12 |
| `se_model/spikes_only` | 12/12 | **0/12** |
| `se_model/sharp_spikes` | **4/1** | 1/1 |
| `se_basic/gauss_two_peaks` | **1/0** | 0/0 |

`subtract` **overcalls** on `sharp_spikes` and `gauss_two_peaks`; `poisson`
**undercalls** `spikes_only` badly enough to return nothing at all (its q-score
maximum is 1.032 against a `-log10(0.05)` = 1.301 cutoff, so no position clears
it). Each mechanism is right where the other is wrong, which is the signature of
two different score functions rather than one function with a tuning bug.

Correcting an earlier statement in this file: the harness default is `poisson`,
not `subtract` (F47b claimed otherwise). The default is deliberately left as
`poisson`, because that is upstream's actual mechanism, and switching it to
`subtract` merely to make more boundaries match would hide the divergence
instead of recording it. Both are selectable via `CALLPEAK_SCORE`.

Where boundaries *do* agree, they agree exactly -- `spikes_only` peak 0 is
`7409-8105` on both sides under `subtract` -- so the pileup, lambda and region
segmentation chain is sound and the divergence is confined to the score function
and its q-value consequences. Summits remain off by up to ~174 bp on
`spikes_only`, well outside the <=10 bp criterion.

This is the clearest statement of remaining G9 work: the score mechanism must be
made faithful, and until it is, per-fixture parity cannot be claimed.

## F49 — a proper peak gate, and a correction

`oracle/run_peak_e2e.sh` promotes the peak comparison from an ad-hoc binary
invocation to a repeatable gate, reporting peak count, boundary agreement,
summit agreement and the worst boundary gap -- for **both** score mechanisms, so
the divergence cannot be hidden behind whichever one happens to look better.

```
=== mechanism: poisson ===                         === mechanism: subtract ===
fixture              rust  up   bnd  summ  worst   rust  up   bnd  summ  worst
se_model/realistic    12   12     0     1    262     12   12    10     0     357
se_model/sharp_spikes  1    1     0     0      6      4    1     0     0  11490
se_model/spikes_only   0   12     0     0      0     12   12    10     0     122
se_basic/gauss_two_peaks 0  0     0     0      0      1    0     0     0       0
se_shapes/bias         0    0     0     0      0      1    0     0     0       0
se_dup/dup_rate_0.25   0    0     0     0      0      1    0     0     0       0
```

**Correction.** Two earlier claims in this file that peak boundaries matched
"all 12" on `se_model/realistic` were wrong. The full table shows **10 of 12**,
with the worst disagreement at 357 bp. The original claim came from reading a
grep-filtered subset of the diff output -- the same class of mistake as the
`spikes_only` false positive in F48, and the reason the gate now parses the whole
diff rather than counting filtered lines. Three of the four "0 upstream peaks"
rows are fixtures where upstream's golden XLS contains no data rows at all, so
they are vacuous and must not be read as passes.

The table also makes the F48 picture quantitative: `poisson` finds 1 peak on
`sharp_spikes` where `subtract` finds 4 against upstream's 1, and `poisson`
finds none on `spikes_only` where `subtract` finds 12 against upstream's 12 with
10 exact boundaries. `subtract` is clearly closer to upstream on the fixtures
where upstream produces peaks, which is worth stating plainly even though it is
the mechanism F46 showed to be wrong -- "wrong mechanism, closer answer" is a
real and useful diagnostic, not a reason to ship it.

## F50 — measured shape of the p-score track on `spikes_only`

`CALLPEAK_HIST=<path>` dumps the p-score track as `pscore_bits<TAB>bases` per run,
so the track can be inspected without going through the q-value walk. On
`se_model/spikes_only` with the `poisson` mechanism:

```text
rust: 162 distinct p-score values, 92,744 scored bases
      top buckets: 6.00 (6 bases), 5.82 (21), 5.77 (3), 5.68 (8), 5.63 (3)
```

Upstream's `--cutoff-analysis` for the same fixture starts at the ladder bucket
`pscore 7.50 -> qscore 5.37 -> npeaks 1 -> lpeaks 200`, and reports 3 peaks at
`6.90`, 13 at `0.90`.

Since `pscore` itself is bit-exact (F47), and both the treatment pileup and the
merged control lambda already match upstream's `--bdg` to 5e-6 (F42), the
difference has to come from *what is fed in*. Ruled out by measurement:

* **`lambda_bg` is not the cause.** `spikes_only` is control-dominant (treat 168,
  control 1989), so `tocontrol` is false and `lambda_bg = treat_sum/gsize =
  168*200/2e6 = 0.0168` -- negligible against a local lambda that is itself
  scaled by `ratio = 0.0844`. So the baseline cannot suppress p-scores.
* **The AFDR walk is not the only suspect.** Even with a correct walk, a track
  whose highest value is 6.00 cannot produce the q-score 5.37 that upstream
  reports at its top bucket, so the *track* differs first.

Careful correction on method: upstream's `lpeaks` column is **segmented peak
length, not a base count**, so it cannot be diffed directly against this dump.
Earlier in this investigation that comparison was attempted and it was invalid.
The ladder bucket values (`pscore`, `qscore`) are comparable; the length columns
are not.

Remaining candidate, and the concrete next step: upstream's `ScoreTrack` is
handed the **paired** arrays `d_pileup_d` / `ctrl_d_pileup_d`, whose treatment
side is `d_pileup_d[i]` indexed by *array position*. This port passes a
`SignalTrack` restricted to `0..min(treat.end, ctrl.end)`. If the paired arrays
carry a different treatment value than the raw scaled pileup at the same
coordinate -- which the F43/F45 end-indexing conventions make plausible -- the
p-score track shifts even though both `SignalTrack`s are individually correct.
Checking that requires the paired arrays themselves, which are only reachable
from inside the compiled `CallPeakUnit` (not patchable, per F42).

## F51 — `ScoreTrack` is fed by a two-pointer walk that can pair a position with a control value from *ahead* of it

`make_ScoreTrackII_for_macs` (`BedGraph.py:1170`) does **not** zip the treatment
and control tracks pointwise. It walks two independent run pointers:

```python
p1 = p1n(); v1 = v1n()          # treatment run
p2 = p2n(); v2 = v2n()          # control run
while True:
    if p1 < p2:
        retadd(chrom, p1, v1, v2)     # treatment pos p1, control value v2 ...
        p1 = p1n(); v1 = v1n()         # ... but only p1 advances
    ...
    elif p1 > p2:
        retadd(chrom, p2, v1, v2)
        p2 = p2n(); v2 = v2n()         # only p2 advances
    else:
        retadd(chrom, p1, v1, v2); both advance
```

Each `retadd` supplies **one** position with a treatment value and a control
value, and advances only the pointer that was behind. So when `p1 < p2`, the
position `p1` is paired with `v2` -- the value of the control run whose start
`p2` is *beyond* `p1`. The control value attached to a treatment position can
therefore come from a region of the control track that does not yet cover that
coordinate.

This is a genuine upstream quirk, not a rounding artefact, and it is invisible
to any correct pointwise implementation: `pscore_track` here zips treatment and
lambda, which is precisely the thing this walk is not. It also explains the
whole F46/F48/F50 pattern without needing any further error:

* `pscore` itself is bit-exact (F47), because the function is fine;
* both input tracks match upstream's `--bdg` to 5e-6 (F42), because both are
  fine;
* but the *pairing* of a treatment value with a control value differs, so the
  p-score track differs, so the q-value walk differs, so peaks differ;
* and the effect is fixture-dependent -- it only bites where treatment and
  control runs interleave, which is why `subtract` (which ignores the pairing)
  sometimes lands closer to upstream than the faithful `poisson` does.

Fixing this means reproducing `retadd` exactly: a two-pointer walk over the two
run lists that emits one `(position, treat, ctrl)` per step, advancing only the
lagging pointer, with upstream's handling of the trailing run. That is the
concrete next step, and it is a wiring change in `macs-peaks`, not a numerical
one -- the p-score arithmetic is already proven exact.

### F51a — the `retadd` pairing is faithful, and it changes nothing: hypothesis not supported

F51 was transcribed rather than just described. `macs-callpeak-e2e` now builds the
p-score track by reproducing upstream's lagging-pointer walk exactly -- two run
pointers, `retadd` at `min(p1, p2)`, advancing only the pointer that was behind,
both on equality, stopping on whichever input is exhausted, with
`pseudocounted_inputs` for the int cast.

Result, via `oracle/run_peak_e2e.sh`:

```text
=== mechanism: poisson ===   (identical to the pointwise-zip build)
fixture                             rust upstream  bnd_ok summits worst_bnd
se_model/realistic                    12    12       0       1       262
se_model/sharp_spikes                  1     1       0       0         6
se_model/spikes_only                   0    12       0       0         0
se_basic/gauss_two_peaks               0     0       0       0         0
se_shapes/bias                         0     0       0       0         0
se_dup/dup_rate_0.25                   0     0       0       0         0
```

Byte-for-byte the same table as before. So the lagging-pointer walk is a real
upstream quirk but it is **not** what makes this port's p-score track differ:
on these fixtures it yields the same pairing as a correct pointwise zip. F51's
explanation of the F46/F48/F50 pattern was wrong, and it should not be carried
forward as the cause.

The transcription is kept rather than reverted. It is a faithful reading of
upstream, it costs nothing, and it removes a whole class of "did you check the
pairing?" doubt from the next person. `macs_score::pseudocounted_inputs` is now
exported for it.

So the p-score path now has *three* independently verified-correct components --
`pscore` (bit-exact, F47), the control merge (identical intersection semantics,
F47a), and the pairing walk (just measured identical) -- and the p-score track
still disagrees with upstream. The divergence is therefore upstream of all
three: in the **bedGraph tracks themselves** that feed `make_ScoreTrackII_for_macs`,
not in what happens to them afterwards. The `--bdg` files written by `callpeak`
are a *different* representation from the in-memory `bedGraphTrackI` fed to the
score track, and matching the written file does not imply the in-memory one
matches. That is the next thing to check, and it needs the in-memory tracks,
reachable only from inside the compiled `CallPeakUnit`.

### F52 — the p-score track and histogram construction are verified faithful

The last unverified step in the p-score path is now checked against source.
`ScoreTrackII.add` stores only `endpos`, `chip` and `control` -- it computes no
score -- so the score is produced later, in `track2bedgraph`:

```python
for chrom in sorted(self.data.keys()):
    prev_pos = 0
    pos = self.data[chrom][0]; p = self.data[chrom][1]
    c  = self.data[chrom][2]; v  = self.data[chrom][3]
    ln = self.datalength[chrom]
    for i in range(ln):
        v[i] = get_pscore(cython.cast(cython.int, (p[i] + self.pseudocount)),
                          c[i] + self.pseudocount)
        tmp_l = pos[i] - prev_pos
        self.pvalue_stat[v[i]] += tmp_l
        prev_pos = pos[i]
```

This confirms, against source rather than assumption:

* the pseudocount (default **1.0**) is added to **both** terms, and the
  **treatment** term is then truncated to `int` -- exactly
  `pseudocounted_inputs(treat, ctrl, 1.0)`;
* chromosomes are walked in `sorted(self.data.keys())` order, matching the
  documented requirement on `QScoreSink`;
* the histogram is keyed by the f32 score bit-pattern with value
  `pos[i] - prev_pos`, i.e. **end-indexed run lengths from 0**, which is what
  `CALLPEAK_HIST` dumps.

So the p-score pipeline now has four independently verified-correct stages --
input pairing (F51a), per-position score (F47), histogram accumulation (this
finding), and q-value walk (G7, 23 tables) -- and the track still disagrees with
upstream. That is now a strong result in its own right: the divergence is
confined to the **values of the in-memory `bedGraphTrackI` objects** built from
the pileups, which the `--bdg` files do not expose and which live only inside
the compiled `CallPeakUnit`.

Practical consequence for the next step: capturing those objects requires
patching the upstream oracle source (`BedGraph.py` / `CallPeakUnit.py` are Cython
sources present in the pinned checkout, so `CallPeakUnit.c` can be regenerated),
not a runtime hook. That is a one-off harness change, and it is the only way to
close this without guessing.

## F53 — upstream's `ScoreTrack` is constructible in pure Python, and my track is run-for-run identical in shape

F52 concluded the divergence was confined to values in the in-memory
`bedGraphTrackI` objects and that capturing them needed a Cython rebuild. That
conclusion was too pessimistic: `bedGraphTrackI` and `ScoreTrackII` are ordinary
importable classes, so upstream's score track can be built directly from the
**golden `--bdg` files** with no rebuild and no runtime hook:

```python
from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.Signal.ScoreTrack import ScoreTrackII
t1 = bedGraphTrackI();  [t1.add_loc(chrom, start, end, value) for each bdg run]
t2 = bedGraphTrackI();  ...
ret = t1.make_ScoreTrackII_for_macs(t2, 1.0, 1.0)
pos, treat, ctrl, val = ret.get_data_by_chr(b"chr20")
```

On `se_model/spikes_only`:

```text
upstream ScoreTrack rows = 4783
rust     p-score track    = 4783 runs   (CALLPEAK_HIST, header excluded)
```

**Exact run-for-run agreement.** So the p-score track's *shape* -- how many
positions upstream scores, and therefore the exact length that feeds the AFDR
denominator -- is already right in this port. That eliminates the denominator
hypothesis I raised in F50 outright, and it means the q-value walk cannot be
salvaged or broken by a length mismatch.

Two caveats, both worth carrying forward:

* `ScoreTrackII` has no `compute_pvalue`; that method lives on `ScoreTrackI`
  (used by `CallPeakUnit`). So `val` comes back all-zero from `ScoreTrackII`
  alone, and the per-position p-scores still cannot be read this way. The run
  *positions* are readable and are what was compared here.
* the `treat`/`ctrl` values read back through this path were 0 at the positions
  sampled, so the value comparison is not yet meaningful -- `add_loc`'s argument
  order for an end-indexed track (F5) needs checking before the values can be
  trusted.

Net effect: F50's denominator theory is dead, F47a/F51a/F52 stand, and the
remaining work is a **value-level** comparison on a track whose structure already
agrees -- which needs the `add_loc` convention pinned, not another hypothesis
about mechanisms.

## F54 — upstream's `--bdg` and its score track are built from **different arrays**

The value-level check now says something decisive. Evaluating upstream's own
`get_pscore` on the values in upstream's own golden `--bdg` files, at the position
of maximum treatment on `se_model/spikes_only`:

```text
pos=15611   treat=8.0000   ctrl=0.2770   ->  pscore 5.99933   (upstream's formula)
rust p-score track at the same position  ->  5.99932          (agree to f32)
```

So this port's p-score track **is** correct given the `--bdg` values, and every
stage downstream of the inputs is now verified (F47, F47a, F51a, F52, F53).

But that creates a contradiction that localises the bug to the *inputs*:

* the maximum treatment anywhere in upstream's `--bdg` for this fixture is **8.0**,
  and `get_pscore(9, 1.277) = 5.999`;
* upstream's `--cutoff-analysis` nonetheless reports a top ladder bucket of
  **`pscore 7.50 -> qscore 5.37`**, which requires roughly `treat ~ 11-12`, or a
  far smaller local lambda, than anything the `--bdg` contains.

A p-score of 7.50 is simply unreachable from `treat = 8`. Therefore the arrays
feeding `make_ScoreTrackII_for_macs` are **not** the arrays written to `--bdg`.

This resolves the whole F46 -> F54 investigation in one step. `callpeak --bdg`
writes `d_pileup_d` and `ctrl_d_pileup_d` from the *paired* arrays used for peak
calling, while `compute_pvalue` feeds `ScoreTrack` from a different in-memory
`bedGraphTrackI` pair -- most plausibly the same pileups **before** the treatment
`scaling factor` and/or **before** the control `lambda_bg` baseline is applied.
For `se_model/spikes_only` `treat_scale` is 1.0, so the scaling is not the
difference there; the `lambda_bg` contribution or the pre-pairing extent is.

Concretely, this means the `--bdg` differential (`oracle/run_e2e.sh`, green at
97/97) validates the *peak-calling* signals, and must **not** be read as
validating the score track. They are two different quantities. Any future work
should compare the score track against upstream's in-memory pair -- now
reachable in pure Python per F53, by building both `bedGraphTrackI` objects the
way `compute_pvalue` does rather than from the written files.

## F55 — `callpeak` never uses `ScoreTrack`; its score track has **no pseudocount**

The root cause of the entire G9 divergence, and it had been one assumption deep
the whole time.

`ScoreTrack` appears in exactly four files:

```text
Signal/ScoreTrack.py
Signal/BedGraph.py          (make_ScoreTrackII_for_macs)
Commands/bdgcmp_cmd.py
Commands/bdgdiff_cmd.py
```

**`callpeak` does not use it.** `CallPeakUnit` builds the score track with its own
`__cal_qscore`, and that function applies **no pseudocount**:

```python
s_ptr[0] = self.pqtable[get_pscore(cython.cast(cython.int, a1_ptr[0]), a2_ptr[0])]
```

The treatment is truncated to an integer and the control is used as-is. The
pseudocount in `macs_score::pseudocounted_inputs` is transcribed from
`ScoreTrack`'s `track2bedgraph` (F52) -- the *other* code path.

Measured on `se_model/spikes_only` at the peak position:

```text
get_pscore(8, 0.277)          = 10.686   <- __cal_qscore, what callpeak does
get_pscore(int(8+1), 0.277+1)  =  5.999   <- ScoreTrack, what this port did
upstream cutoff-analysis top bucket: pscore 7.50
```

Only the first can reach the 7.50 ladder bucket, which is exactly what F54 deduced
from the numbers ("requires roughly treat ~ 11-12") without identifying the
mechanism.

**Effect** (`oracle/run_peak_e2e.sh`, `poisson` mechanism -- upstream's real one):

| fixture | before | after |
|---|---|---|
| `se_model/realistic` | 12 peaks, 0 exact bnd, worst 262 bp | 12 peaks, worst gap **0 bp** |
| `se_model/sharp_spikes` | 1 peak, worst 6 bp | 1 peak, worst gap **0 bp** |
| `se_model/spikes_only` | **0 peaks** | 12 peaks, worst gap **0 bp** |

Peak counts now match upstream exactly on every fixture where upstream calls
peaks, and **every peak boundary is byte-identical**. Summit offsets on
`se_model/realistic` fell from 8-95 bp to **3-21 bp** (6, 11, 19, 12, 3, 14, 14,
11, 8, 21, 7) -- still all low, so a residual tie-break or `ttreat_p` detail
remains, but the score path itself is now correct.

This retires F46, F47b, F50, F51a and F54 in one step, and it is a warning worth
keeping: F52 verified the pseudocount faithfully against a function that
`callpeak` does not call. Faithfulness to the wrong reference proves nothing.

## F56 — after F55, peak counts and boundaries are exact across the whole single-end corpus

`oracle/run_peak_e2e.sh --all` (mechanism `poisson`, i.e. upstream's real one),
restricted to single-end fixtures the harness can actually read. Every fixture on
which **upstream calls any peak**:

```text
fixture                          rust  upstream  bnd_ok  worst_bnd
onechrom/single_contig               1         1       1           0
se_edge/disjoint_chromosomes         1         1       1           0
se_edge/very_shallow                 1         1       1           0
se_model/realistic                  12        12      11           0
se_model/sharp_spikes                1         1       1           0
se_model/spikes_only               12        12      11           0
```

* **peak counts match upstream exactly on every one**;
* **worst boundary disagreement is 0 bp everywhere** -- every peak's `start` and
  `end` is byte-identical to upstream's.

The remaining rows report `0 0`, which is **vacuous**: upstream's golden XLS has no
data rows for those fixtures (no significant peaks at `qvalue 0.05`), so they are
not passes and must not be counted as such.

Also fixed this turn: `run_peak_e2e.sh --all` was discovering fixtures three ways
wrong at once -- wrong `find` depth (variants live one level deeper), no
deduplication of variant subdirectories, and no filter excluding paired-end and
FRAG fixtures, which this single-end harness can only report as a meaningless
`0 vs 0`. All three are fixed; the sweep now lists each fixture exactly once and
only where the comparison is meaningful.

With F55 in place the `subtract` mechanism is now clearly dominated -- 10 bnd_ok
and a 357 bp worst gap on `realistic`, against 11 and 0 for `poisson`. It is kept
only as a regression witness for F55 and should not be revisited.

**Open:** summit coordinates remain 3-21 bp low on `se_model/realistic` (6, 11,
19, 12, 3, 14, 14, 11, 8, 21, 7). With the signal arrays and boundaries now
exact, the residue is confined to the `tscore` argmax in
`__close_peak_wo_subpeaks` -- specifically how `ttreat_p` is attached to each
`peak_content` entry and the `(tend + tstart) // 2` midpoint that follows. Next
step: compare the per-entry `peak_content` tuples, which are now reachable because
the score track no longer has to be guessed.

## F57 — chunk granularity was already right: a third negative result on the summit

Upstream's `peak_content` holds one entry per above-cutoff **run**; a new peak
*region* is started only by `ts - lastp > max_gap`, never by run adjacency. This
port was merging contiguous above-cutoff runs into a single `Chunk`, which widens
it -- and since the summit is `(tend + tstart) // 2`, widening pushes the summit
earlier. That is precisely the one-sided low bias observed on every peak, so it
was the obvious remaining suspect.

Removing the merge (one `Chunk` per above-cutoff run, matching upstream's
granularity exactly) produced **no change at all**:

```text
before: 7726 15495 23037 30980 38607 46181 61729 69401 76895 84649 92358
after:  7726 15495 23037 30980 38607 46181 61729 69401 76895 84649 92358
upstream: 7732 15506 23056 30992 38610 46195 61743 69412 76903 84670 92365
```

The merge condition (`last.end == lo`) was evidently already almost never
firing, because the q-track's run boundaries coincide with the above-cutoff
boundaries on these fixtures.

The change is **kept** -- it is the faithful granularity and costs nothing -- but
it is recorded as a *null result*, not a fix. Third consecutive negative on the
summit (F47a, F51a, F57), all recorded rather than dropped.

Where that leaves the summit, stated precisely:

* treatment pileup, control lambda, the q-track shape (F53), `pscore` (F47),
  chunking (F57) and peak `start`/`end` (F56) are all confirmed exact;
* the residue is a one-sided 3-21 bp low bias, i.e. upstream's chosen run sits
  slightly *later* than this port's within an identical region;
* since `tsummit` is built on `summit_value < tscore` with ties appending and the
  **lower median** taken, a strictly-low result means this port finds **more
  entries tied at the maximum** than upstream does -- a tie-breaking or
  `ttreat_p`-equality difference, not a wrong run.

That makes the next check specific and small: count how many runs tie at the
maximum `ttreat_p` per peak on both sides. If this port ties more runs, the
f32 comparison of `tscore` is where to look; if it ties fewer, the `int`
truncation of `ttreat_p` is.

## F58 — `ttreat_p` = max over the run is confirmed better, but not sufficient

With F55 in place, the two candidate definitions of `ttreat_p` were measured
directly on `se_model/realistic`, all else identical:

| `ttreat_p` | mean summit offset | peaks compared |
|---|---|---|
| value at the run's **end** position | **-18.17 bp** | 12 |
| **max over the run** (F47c) | **-11.45 bp** | 11 |

Max-over-run is the better of the two, so F47c's direction survives the F55 fix --
but neither reaches zero, and the gap between them (6.7 bp) is smaller than the
residual (11.45 bp). So the choice of `ttreat_p` is *not* the remaining error;
something upstream of it still is.

Note the peak count differs between the two runs (11 vs 12 compared), which means
the `ttreat_p` definition also changes which chunks win the argmax enough to move
a peak boundary. Both variants keep `worst_bnd = 0`, so this is confined to
summits.

Given the residual is now a uniform ~11 bp low bias across a peak whose extent is
exactly right, and that upstream's `(tend + tstart) // 2` is an integer midpoint,
the leading remaining candidates are (a) `ts`/`te` for the winning chunk being the
`above_cutoff_startpos`/`above_cutoff_endpos` pair (`pos_array[above_cutoff-1]` /
`pos_array[above_cutoff]`) rather than the run's own `[start, end]` -- these differ
by one position at the left edge, which shifts a midpoint by at most 1 and so
cannot alone explain 11 bp; or (b) the argmax landing on a different run because
upstream's `ttreat_p` is `treat_array[ti]` at a *position* rather than a maximum,
so upstream's winning entry can be an earlier position of a later run.

(b) is the better explanation and is directly testable by emitting this port's
per-chunk `(start, end, ttreat_p, midpoint)` for one peak on both the max-over-run
and per-position definitions, and checking which one yields upstream's 7732 for
peak 0.

## F59 — the summit's tie set is right; the chunk *boundaries* carry an extra split

`CALLPEAK_CHUNKS=1` dumps this port's per-chunk `(start, end, ttreat_p, mid)`
for peak 0 of `se_model/realistic`. The peak spans 7572-8039, and the maximum
`ttreat_p` is **9.0**, attained by exactly four chunks:

```text
chunk   7703-7716  ttreat=9.000  mid=7709
chunk   7717-7735  ttreat=9.000  mid=7726
chunk   7833-7839  ttreat=9.000  mid=7836
chunk   7840-7852  ttreat=9.000  mid=7846
```

`tsummit = [7709, 7726, 7836, 7846]`, `midindex = (4+1)//2 - 1 = 1`, so
`summit_pos = 7726` -- exactly what this port reports.

**Upstream reports 7732, which is not any of these midpoints.** The only chunk
with midpoint 7732 in this neighbourhood is one spanning `7717-7747`. Upstream
therefore has a *single* chunk where this port has two, split at `7736`:

```text
upstream:  chunk 7717-7747                      mid=7732   <= MAX, selected
rust:      chunk 7717-7735  ttreat=9.000        mid=7726   <= MAX, selected
rust:      chunk 7736-7746  ttreat=8.000        mid=7741
```

So the tie set and the selection rule are both correct -- the error is that this
port's q-track has an **extra run boundary at 7735/7736 that upstream does not
have**, which splits the winning run and shifts the midpoint of the selected
chunk 14 bp earlier.

Where the extra boundary comes from: chunk boundaries are derived from the q-track
run ends, and the q-track comes from the `retadd` two-pointer walk (F51a), which
emits a run at **every treatment or control run boundary** -- the union of both.
Upstream derives `te` from `pos_array[above_cutoff]`, i.e. the positions of the
**paired** arrays, which have coarser breakpoints: `__chrom_pair_treat_ctrl`
(F35) truncates at the paired extent and carries only positions where both
signals are defined.

This is consistent with F53's run-count match being exact on `spikes_only` while
still differing in *where* some boundaries fall here, and it explains the residual
uniformly: every selected chunk is split at a spurious boundary, so every summit
lands low, by up to half the width of the over-split run (14 bp here, 3-21 bp
across the corpus -- F56/F58).

**The fix, now precisely stated:** chunk boundaries must come from the paired
positions, not from the union-walk run ends -- emit a run only where *both* the
treatment and the control track change.

### F59a -- an attempted implementation of F59 was wrong, and is reverted

A first attempt kept emitting at `min(p1, p2)` but substituted the previously
emitted position when the two pointers disagreed. That is not the same thing: it
does not remove the spurious boundary, it *stretches* the preceding run's end
over it, corrupting the run lengths and therefore the histogram. Measured:

```text
mean summit offset +23119.67 over 3 peaks     (was -11.45 over 11-12)
```

Catastrophically worse, and the peak count collapsed -- the tell that the run
structure itself was broken. Reverted; the tree is back to `mean -11.45,
worst_bnd 0` and 333 passing.

The correct implementation is a three-way walk over the two run lists:

```text
while i1 < len(p1s) and i2 < len(p2s):
    if p1 == p2:  emit run ending at p1 with (v1, v2); advance both
    elif p1 < p2:  advance i1 only  -- do not emit
    else:          advance i2 only  -- do not emit
```

so the emitted positions are exactly the **coincident** run boundaries, and each
emitted run carries the `(v1, v2)` in effect at its end. Note this also changes
the run *lengths*, so the p-score histogram is recomputed -- which is the point,
since the current lengths are known to be wrong (F59).

Two details to get right: flush a final partial run when either pointer is
exhausted, and give each emitted run the extent
`(previous coincident boundary, this one]` rather than `(previous emitted, this one]`.

### F59b -- F59 implemented correctly: coincident boundaries yield **zero peaks**

F59a's three-way walk was implemented properly this time and measured via
`CALLPEAK_BOUNDARY=paired|union` (default `union`, the working path):

```text
paired   ->  rust peaks: 0, upstream peaks: 12      (no summit comparisons)
union    ->  rust peaks: 12, upstream peaks: 12, mean summit offset -11.45
```

So the coincidence reading is **wrong**: treatment and control run ends almost
never coincide, so emitting only at coincident positions leaves an almost-empty
score track and no peaks at all. F59a's +23119 result was a different bug (broken
run lengths); this one is a genuine refutation.

That does not rescue F59's observation, though -- the underlying measurement
stands and is still unexplained:

```text
upstream peak 0: chunk 7717-7747  mid=7732   (selected)
rust      peak 0: chunk 7717-7735  mid=7726   (selected)
                    chunk 7736-7746  mid=7741
```

The winning chunk really is split in two here, and the split is real rather than
an artefact of the coincidence reading. What remains to explain is where 7736
comes from: this port emits it because it is `min(p1, p2)` for some pair of run
ends, and upstream evidently has no position there. Since coincident-only is not
the answer, the next thing to check is whether the treatment and control
**run ends themselves** differ -- i.e. whether the merged control's runs are
carrying boundaries the paired arrays would not have, rather than whether the
walk combines them wrongly. That is a question about `ctrl_d_pileup_d`'s run
structure, not about how two lists are merged.

Both boundary modes remain reachable by environment variable so the measurement
is reproducible; `union` is the default because it is the only one that
reproduces upstream's peak counts and boundaries.

### F59c -- the split boundary comes from the **treatment** pileup, which `--bdg` cannot see

`CALLPEAK_RUNS=1` dumps the raw run structure around the split, and the answer is
unambiguous:

```text
RUN   7702-7716   treat v=8.0000
RUN   7716-7735   treat v=9.0000
RUN   7735-7746   treat v=8.0000
          (control track has NO runs in this range at all)
```

So the boundary at `7735` that splits the winning chunk comes from the
**treatment** pileup's own run structure -- not from the control, and not from the
way the two lists are merged. That eliminates the merge as the cause (consistent
with F59b) and localises it upstream, in the treatment track's granularity.

And that exposes a real gap in the gate that has been green this whole time.
`oracle/run_e2e.sh` compares the `--bdg` files, and
`__write_bedGraph_for_a_chromosome` only emits a line when the value changes by
more than `1e-5` (F42). A run boundary whose two sides differ by less than that
is **invisible in the file**. So this port's treatment track can carry a spurious
run end at `7735` -- splitting the chunk that wins the summit argmax -- while
still matching `--bdg` perfectly, because the values on either side round to the
same five decimals.

That means the 97/97 signal gate is weaker than it looks: it certifies the
*values*, not the *run boundaries*, and the summit depends on the boundaries.
Most likely source: an `f32` rounding difference between this port's scaled
treatment pileup and upstream's, large enough to split a run but far below the
`1e-5` render threshold.

The fix is not to keep tweaking the summit logic -- it is to compare **run
structure**, not rendered values. Two options:

1. dump upstream's `d_pileup_d` run ends directly (needs the compiled
   `CallPeakUnit`, F42) and diff the boundary sets; or
2. tighten the oracle harness to write the bedGraph **without** the `1e-5`
   coalescing, so every run end becomes visible and `run_e2e.sh` starts
   comparing structure rather than a 5-decimal rendering.

Option 2 is the better gate and does not depend on reaching into compiled code;
`dump_stages.py` already writes those files and controls what they contain.
Until then, treat "signal parity 97/97" as value parity only, not structural
parity.

### F59c correction -- the boundary at 7735 IS visible in `--bdg`; F59c is wrong

F59c claimed the spurious boundary is invisible in the coalesced `--bdg` file and
that the 97/97 gate therefore certifies only values, not run structure. Checked
directly against the file:

```text
$ awk -F'\t' '$1=="chr20" && $2>=7690 && $2<=7770' stages_treat_pileup.bdg
7702-7716 8.00000
7716-7735 9.00000
7735-7746 8.00000
7746-7771 7.00000
```

The boundary at **7735 is present in the file**, and this port's run structure
matches it exactly (`7702-7716 / 7716-7735 / 7735-7746 / 7746-7771`). So the two
differ by 8.0 and 7.0 here -- far above the `1e-5` coalescing threshold -- and the
run structure *is* visible and *is* being compared by `run_e2e.sh`.

**F59c's gate-gap conclusion is withdrawn.** The `--bdg` comparison does certify
run structure at this position, and the withdrawal matters: it was about to send
the next person to change the harness for no reason.

What the corrected picture gives instead is sharper than the original claim. Both
sides agree the treatment pileup has runs `7716-7735` (9.0) and `7735-7746` (8.0).
Upstream nevertheless selects a **single** chunk `7717-7747` spanning both, while
this port selects two. Since the chunk boundaries come from
`above_cutoff_startpos = pos_array[above_cutoff - 1]` and
`above_cutoff_endpos = pos_array[above_cutoff]` (F45), and since both sides agree
the treatment runs break at 7735, the only way upstream gets one chunk is for
**7735 not to be in its `above_cutoff` set** -- i.e. upstream's q-value at 7735
falls *below* the cutoff where this port's clears it.

So the residual is a single-position q-value threshold crossing at 7735, not a
granularity problem. That is consistent with everything else measured: the
treatment/control values agree, `pscore` is bit-exact (F47), the score track's
shape agrees (F53), and a value that is only marginally above the cutoff can sit
on either side of it after the q-value walk.

The check is now specific and cheap: print this port's q-value at 7735 against the
cutoff `-log10(0.05) = 1.30103`, and see how far it is from crossing. If it is
within ~1e-4, this is ordinary f32 threshold sensitivity at the `>` comparison
(F45) and the fix is to match upstream's comparison precision exactly rather than
to change any structural logic.

### F59d -- the q value at 7735 is 7.10, not marginal: threshold sensitivity refuted too

F59c's correction suggested the residue might be ordinary f32 threshold
sensitivity at the `>` comparison. Measured directly (`CALLPEAK_Q=1`):

```text
Q   7702-7716   v=7.09660673  above=true
Q   7716-7735   v=7.99770498  above=true
Q   7735-7746   v=7.09660673  above=true      <- the disputed position
Q   7746-7771   v=6.04201889  above=true
Q cutoff = 1.30103004
```

The q value at 7735 is **7.10 against a cutoff of 1.30** -- five times over, not
within `1e-4` of it. So it is not threshold sensitivity, and no amount of care in
the comparison operator will move it. **Refuted.**

That makes the requirement sharp and unusual. For upstream's `above_cutoff` set
to skip 7735 -- which it must, since `peak_content` pairs *consecutive*
`above_cutoff` entries into `(ts, te)` and upstream's selected entry spans
`7716..7746`, skipping 7735 -- upstream's q value there has to be **below 1.30**
where this port's is **7.10**. Meanwhile the treatment and control values at that
position agree exactly (F59c correction: `7735-7746 v=8.00000` in both), and
`pscore` is bit-exact (F47).

Those three facts cannot all be true of a single shared input. Since `pscore` is
bit-exact and the inputs agree, the p-score at 7735 must agree, so the
discrepancy is in the **p-score -> q-value walk** for this track, not in the
score. The walk's output for a given p-score depends on the whole histogram
(G7), and G7's 23 tables were verified on synthetic inputs -- not on this
fixture's real histogram. So the first thing to check is whether the p-score
histogram this fixture produces is one of the verified cases, and if not, feed it
through `check_precomputes.py`.

That is a concrete, bounded next step and it is the first thing in this whole
investigation that the existing 333-test suite does *not* already cover.

### F59e -- a `(pscore, qscore)` dump I built is misaligned; do not trust it

The obvious way to test F59d's hypothesis (the p->q walk on a real histogram) is
to dump `(pscore, qscore)` pairs and compare with the `pscore`/`qscore` columns of
upstream's `--cutoff-analysis`. A first version of that dump produced a dramatic
and completely wrong result:

```text
mine      p=10.84893  q=0.09170        (top p-scores, ~zero q)
upstream  p= 9.30     q=7.04           (top bucket)
```

which would mean the walk assigns near-zero q to the most significant p-scores --
an inversion. It is **not** real. `CALLPEAK_Q`, which reads the q-track directly
and is not subject to any pairing assumption, shows this port's q-track values in
the same neighbourhood as `7.09660673` and `7.99770498` (F59d) -- i.e. healthy
high q-values where they should be.

So the dump's pairing is wrong, most likely because `pv.get(i)` and
`qt.runs()[i]` are not in one-to-one correspondence for this track (the q-track
can carry a different number of spans than the p-track it was built from, and a
run-by-run zip silently misaligns rather than failing). The diagnostic was
removed rather than left in the tree, and it is recorded here because the failure
mode is dangerous: it produced a confident, plausible, wrong inversion.

This is the third time in this investigation a diagnostic has manufactured a
finding that measurement contradicted -- the grep-filtered "12/12" (F49), the
vacuous 0-peak "match" (F48), and F59c's gate-gap claim (withdrawn). The pattern
is worth stating as a rule: **every claim about the q-track must be cross-checked
against a second, independent readout before it is written up.** For this one the
cross-check is `CALLPEAK_Q`, which reads the track directly.

The underlying question from F59d is still open and still worth asking -- G7's 23
tables are synthetic, and a real-histogram check is genuinely missing from the
suite. But it needs a dump whose pairing is verified (assert that p-track and
q-track span counts agree before zipping) before its output means anything.

### F59f -- the p->q walk is sound; but this port's top p-score is 12.06 where upstream's is under 9.6

`CALLPEAK_PQ` rebuilt with the pairing F59e prescribed -- each q-span matched to
the p-span **containing its start position** via `value_at`, never by index.
First finding: span counts legitimately differ, so the index zip really was the
bug:

```text
NOTE 0: p-spans=4763  q-spans=732 (coalesced)
```

`qscore_track` merges the p-score track, because many distinct p-scores share a
q-score and unmapped ones become 0. 4763 -> 732 is expected, and any index-based
zip over that is meaningless. That fully explains F59e's bogus inversion.

With the pairing fixed, the walk is healthy -- high p-score does map to high
q-score, so no inversion:

```text
mine:      p=11.86533  q=8.37252
           p=11.88633  q=8.37252
           p=12.05798  q=8.37252
upstream:  p= 8.70     q=6.55
           p= 9.00     q=6.84
           p= 9.30     q=7.04        (ladder tops out at 9.9)
```

**But a real difference is now visible for the first time.** This port's maximum
p-score is **12.06**; upstream's is **below 9.6** (the 9.6 and 9.9 ladder buckets
have zero length). Since `pscore` is bit-exact (F47), the inputs differ: at the
same treatment and control values this port must be evaluating a smaller lambda,
or a larger treatment.

That is consistent with everything, and it is finally a *magnitude* difference
rather than a structural one. It also explains the summit residue without any
fiddling with `ttreat_p` or the tie-break: if the whole p-score track is shifted
upward, then at the winning chunk the argmax can land on a different position
than upstream's, because the *relative* ordering of near-tied `ttreat_p`-derived
scores changes -- even though the treatment values themselves agree.

Next step is now specific and cheap: at the peak, compare this port's p-score with
`get_pscore` computed from the **same** treatment and control values read out of
the `--bdg` files, exactly as F54 did for one position. If they disagree, the
lambda fed to the walk is the culprit; if they agree, the p->q walk is shifting
them.

`CALLPEAK_PQ` is kept in the tree (env-gated, no effect on normal runs) because
it is now a correct instrument, with the position-based pairing and the span-count
report that make the pairing auditable.

### F59g -- F59f's magnitude difference is an artefact of the cutoff ladder capping at 9.9

F59f reported this port's max p-score as 12.06 against upstream's "< 9.6" and
called it a real difference. It is not. The check F59f prescribed settles it:
sweeping `get_pscore` over peak 0's region using the values read straight out of
upstream's own `--bdg` files,

```text
peak0 region 7717-8039, 322 positions scanned
max get_pscore from --bdg values: 10.84889 at pos=7717 (treat=9.000 ctrl=0.38574)
this port's p-score at the same position: 10.84893
```

**They agree to f32 precision.** So there is no magnitude difference in the p-score
at all.

The apparent one came from the `--cutoff-analysis` table itself: its ladder is
`np.arange(0.3, 10.0, 0.3)`, which **tops out at 9.9**. Any p-score above 9.9 --
and the peak values are ~10.85 -- cannot appear in that table at all. Reading
"upstream's max is under 9.6" off the ladder was reading a property of the ladder,
not of the data. `F59f` is corrected accordingly.

So the p-score and p->q paths are confirmed correct against upstream's own
numbers, on a real fixture, at the position that matters. Together with F47
(`pscore` bit-exact), F52/F53/F54/F55 (mechanism, histogram, pairing, no
pseudocount), that closes off the entire score path as a source of the summit
residue.

**Which leaves the summit unexplained.** Every input and every transform is now
verified against upstream on real data at the disputed position, and the summit is
still a mean -11.45 bp low (F58). The measurements that remain unverified are
exactly the ones that were never made: upstream's own `peak_content` tuples and
`tsummit` list for this peak. That requires the in-memory arrays inside the
compiled `CallPeakUnit` (F42), and F53 showed the pure-Python `bedGraphTrackI` /
`ScoreTrackII` route can read run *positions* but not the per-entry `ttreat_p`
that the summit argmax runs on.

So the honest position after this turn: the score path is exonerated by direct
measurement, and the summit needs upstream's `peak_content` itself -- which is an
oracle-side capture task, not another port-side hypothesis.

## F60 — `peak_content` cannot be captured from Python; the private methods are not exposed

F59g named the remaining step as an oracle-side capture of upstream's
`peak_content`. `oracle/dump_peak_content.py` was written to do it and **it runs
upstream end to end but captures 0 chunks**. Two facts explain that, and together
they close off the route:

```text
>>> [n for n in dir(CallerFromAlignments) if 'close' in n or 'peak' in n]
['call_broadpeaks', 'call_peaks']
```

* `CallerFromAlignments` **is** subclassable, and `callpeak_cmd` does resolve
  `CallerFromAlignments` in its own namespace at call time -- so the interception
  trick from `dump_stages.py` does apply, and the subclass is accepted;
* but `__close_peak_wo_subpeaks`, `__chrom_call_treat_ctrl_a_chromosome` and
  `__pre_computes` are `cdef`-private and simply are **not attributes of the
  class**. There is nothing to override. The same holds for `CallPeakUnit` (F42).

So `peak_content` is only reachable from inside the compiled extension. Two ways
forward, both real work rather than more hypothesis:

1. **Recompile the oracle with an instrumented `CallPeakUnit`.** The Cython
   sources are in the pinned checkout (`Signal/CallPeakUnit.py` plus the `.c`
   it was built from), so a one-line dump inside `__close_peak_wo_subpeaks` and a
   rebuild is mechanical. This is the direct route and gives exactly the
   `(tstart, tend, ttreat_p, tctrl_p, ti)` tuples the summit argmax runs on.
2. **Reimplement `__chrom_call_treat_ctrl_a_chromosome` in the harness** as a
   pure-Python oracle over `bedGraphTrackI` (which F53 showed is constructible in
   pure Python), then diff its `peak_content` against this port's chunks. No
   rebuild, but the harness would be re-deriving upstream's logic rather than
   observing it, so it can confirm the port only as far as the harness itself is
   trustworthy.

Option 1 is the honest instrument. Until one of these exists, the summit residual
(mean -11.45 bp, F58) has been narrowed to a single unobserved quantity and
nothing further can be settled from the outside -- every other input and
transform is verified (F47, F52, F53, F55, F59g).

`oracle/dump_peak_content.py` is kept, with the subclass wiring intact, so that
once the oracle is rebuilt it works immediately; it currently writes a header and
no rows, which is a visible failure rather than a silent one.

### F60a -- F60 option 1 is blocked: Cython is not in the hermetic environment

The `Signal/*.py` files really are Cython sources (they use `cython.cast`,
`cython.pointer`, `@cython.cfunc`), so instrumenting `__close_peak_wo_subpeaks`
and rebuilding is mechanically the right instrument. It cannot be done without
breaking the oracle's hermeticity, though:

```text
>>> import Cython
ModuleNotFoundError: No module named 'cython'
>>> ls MACS3/Signal/CallPeakUnit.pyx
No such file or directory
```

* Cython is **not installed** in the pinned oracle venv (`oracle/ENV.lock`);
* there is no `.pyx` -- the `.py` is the source and `.c`/`.so` were generated
  from it at build time.

So instrumenting the oracle means adding a build-time dependency to a venv that
is pinned by `oracle/ENV.lock`. That is a real tradeoff, not an oversight, and
it is the user's call:

* **option 1a** -- install Cython into the oracle venv and rebuild the two
  extensions (`CallPeakUnit`, `Pileup`). Gets an authoritative `peak_content`
  dump. Cost: the venv stops matching `ENV.lock`, so `ENV.lock` must be
  re-pinned and the oracle rebuilt from source, which weakens "frozen MACS3
  3.0.5 at commit c544319" to "frozen plus one local instrumentation patch".
* **option 1b** -- build the instrumented extension in a *separate* venv or
  prefix, leaving the pinned oracle untouched, and use it only for this
  diagnostic. Preserves hermeticity of the reference; costs a second build and
  the guarantee that the instrumented build behaves identically to the pinned
  one (it should, but it is no longer the same binary).
* **option 2** -- reimplement `__chrom_call_treat_ctrl_a_chromosome` as a
  pure-Python oracle over `bedGraphTrackI` (constructible in Python, F53) and
  diff its `peak_content` against this port's chunks. No rebuild at all. But it
  re-derives upstream's logic rather than observing it, so it can only confirm
  the port as far as the harness itself is right -- which is precisely the kind
  of circularity that made F52 (pseudocount verified against a function
  `callpeak` never calls) a dead end.

I have not chosen. This is a genuine scope/compatibility decision: it trades
oracle hermeticity against the authority of the observation, and it is the
user's call rather than mine. Until it is made, the summit residual (mean
-11.45 bp, F58) stays where F59g left it -- narrowed to a single unobserved
quantity, with every other input and transform verified.

## F61 — `peak_content` captured from upstream; the chunk **boundaries** differ, not the formula

F60a option 1a, done. Cython 3.3.0 installed into the oracle venv (the built
`.so` files were backed up first), a gated dump added at the top of
`__close_peak_wo_subpeaks`, and `MACS3/Signal/CallPeakUnit.py` rebuilt with
`setup.py build_ext --inplace`. **915 `peak_content` tuples** captured for
`se_model/realistic`.

Peak 0, `(tstart, tend)`:

```text
rust:    7572-7590  7591-7608  7609-7635  7636-7652  7653-7657 ...
upstream 7571-7590  7590-7596  7596-7604  7604-7608  7608-7635 ...
```

and around the summit:

```text
upstream  (7702,7705) (7705,7716) (7716,7720) (7720,7730) (7730,7733)
          (7733,7735) (7735,7742) (7742,7746) (7746,7770)
rust      (7703,7716) (7716,7735) (7735,7746)
```

**The formula is right; the boundaries are wrong.** Two things follow immediately.

**First, the summit arithmetic is confirmed as `(tend + tstart + 1) // 2`.** No raw
`(tend + tstart) // 2` equals upstream's 7732, but `(7730 + 7733 + 1) // 2 = 7732`
does -- i.e. the same `+1` that turns `tstart 7571` into the XLS start `7572` also
enters the midpoint. This port already computes `start = tstart + 1` and
`(start + end) / 2`, so this side is already faithful.

**Second -- and this is the actual bug -- this port's chunks are far too coarse
near the peak.** Upstream breaks `7716-7735` into `(7716,7720) (7720,7730)
(7730,7733) (7733,7735)`, four chunks; this port has one. Upstream's chunk
boundaries come from consecutive `above_cutoff` *indices* in the paired
position array -- `ts = pos_array[ti-1]`, `te = pos_array[ti]` -- whereas this
port derives chunk boundaries from the **q-track's RLE runs**, and the q-track
coalesces (F59f: 4763 p-spans become 732 q-spans). Coalescing is exactly what
merges upstream's four chunks into one, which shifts the selected chunk's
midpoint and pulls the summit low.

So the fix is precise and local: chunk boundaries must be taken from the
**p-score track's run ends** (the un-coalesced signal), not from the coalesced
q-track. Every other input and formula in the path is now confirmed correct by
direct observation of upstream's own data.

Caveat to carry: the oracle venv no longer matches `oracle/ENV.lock` (Cython
3.3.0 added, one extension rebuilt from a locally patched source). The pinned
binaries were backed up to `/tmp/oracle_so_backup/` and
`MACS3/Signal/CallPeakUnit.py` to `/tmp/CallPeakUnit.py.orig`. `ENV.lock` should
be re-pinned deliberately rather than silently, and the instrumented build used
only for diagnostics -- `MACS3RS_PC` gates it to a no-op otherwise.

### F61a -- switching chunk boundaries to the p-track regressed; reverted

F61's stated fix was implemented: iterate the **p-score track's runs** for chunk
boundaries instead of the coalesced q-track, and test significance with the
q-track's value at that span. It is a large regression:

```text
before (q-track boundaries):  realistic 12/12, sharp_spikes 1/1, spikes_only 12/12
after  (p-track boundaries):  realistic  1/12, sharp_spikes 0/1,  spikes_only  1/12
```

Reverted; the tree is back to 12/12, 1/1, 12/12 with `worst_bnd 0`.

The obvious flaw is visible on inspection and was not caught before running:
upstream's `above_cutoff` is a selection over indices of the **paired** position
array, so its chunk list is neither the p-track's runs nor the q-track's runs --
it is the subsequence of the paired array's entries whose score clears the
cutoff. Taking *all* p-track runs as candidate boundaries ignores that
selection, so the significance test has to be applied per p-track run against the
q value in effect, and the resulting chunk set must then be grouped by
`ts - lastp > max_gap`. Doing that naively collapsed 12 regions into 1, which
means the run grouping or the significance predicate is wrong -- most likely the
q-value lookup, since the q-track is coalesced and `value_at` on it does not
reproduce the per-index score the selection is made on.

The corrected route is now available because F61's capture exists: rather than
reconstructing the chunk set, **read it**. `/tmp/pc.tsv` already holds upstream's
`peak_content` for `se_model/realistic`; the right next step is to have
`macs-callpeak-e2e` read a `peak_content` capture and use it as the chunk list
directly, which tests the summit stage in isolation from everything upstream of
it. That turns a construction problem into a comparison problem, and it is
decidable today.

Note for whoever picks this up: the capture is available for any fixture by
setting `MACS3RS_PC=<file>` and running `oracle/dump_peak_content.py`.

## F62 -- feeding upstream's own `peak_content` in: the summit formula was off by one

F61a's isolation test is built: `macs-callpeak-e2e --peak-content <file>` replaces
this port's chunk list with upstream's captured tuples, so the summit stage is
tested on its own. Two things came out of it.

**Bug 1: the summit midpoint is `(tend + tstart + 1) // 2`, and this port had
`(tend + tstart) // 2`.** With upstream's chunks and this port's formula, summits
were off by exactly 1 bp. Fixing it took the isolation run from 9 mismatched
peaks to 5, with **9 of 12 summits byte-exact**. Note F61's arithmetic was
consistent with `+1` but ambiguous for it: `(7730+7733+1)//2` and `(7730+7733+2)//2`
both give 7732 because the sum is odd, so the even-sum case -- which is where the
off-by-one actually shows -- could not be distinguished from F61's data alone.
This is the case that needed the summit stage fed directly.

**Bug 2: the remaining three are the `ttreat_p` argmax, now isolated.** Feeding
upstream's chunks still leaves peaks 7, 8 and 10 wrong (61730 vs 61743, 69410 vs
69412, 84666 vs 84670). Since the chunk *list* is now upstream's own, the only
remaining freedom is which chunk wins the argmax -- i.e. the per-entry
`ttreat_p`, which this port computes as the **max treatment over the chunk's
span** (F47c) while upstream uses `treat_array[ti]` at the chunk's own index. The
capture already records `ti` in its fifth column, so the next step is to take
`treat` from the treatment track at exactly that index rather than a span maximum,
and confirm the three peaks land.

That is now a one-line change with a measurement attached, instead of the
reconstruction problem F61a ran into.

Current state of the summit stage, measured:

| run | summits exact |
|---|---|
| this port's own chunks, before F62 | 0 of 12 |
| upstream's chunks, before F62 | 3 of 12 (off by 1) |
| upstream's chunks, after F62 | **9 of 12** |

### F62a -- `ttreat_p` at the chunk's end is *worse* than the span max; F62's hint was wrong

F62 predicted the last three peaks were the `ttreat_p` argmax, and suggested
taking `treat` from `treat_array[ti]` rather than the span max. Since
`pos_array[ti] == tend`, that means the treatment value at the chunk's **end**
position. Measured on `se_model/realistic` with upstream's chunks fed in:

```text
span max (F47c, kept):   5 peak mismatches   (9 of 12 summits exact)
value at tend (F62a):   12 peak mismatches   (0 of 12 summits exact)
```

So the prediction is refuted, and strongly: value-at-end is worse than *no*
progress at all, and 9 of 12 -> 0 of 12. Reverted; the span max stays.

That is the second time a correct-looking derivation about `ti` led somewhere
wrong (F61's own arithmetic had the same flavour), and the pattern is worth
naming: `pos_array[ti] == tend` holds for the **position** array, but
`treat_array[ti]` indexes the **paired treatment** array, which is a different
array with its own run structure. Reading a position out of one array and
assuming the matching value lives in another is exactly the mistake that made
F52 "verify" a pseudocount `callpeak` never applies.

Where the three peaks actually stand, now that the chunk list and the midpoint
are both confirmed: the chunk list is upstream's own, `(tend + tstart + 1) // 2`
is confirmed, and the argmax still differs on peaks 7, 8, 10. So `ttreat_p` is
genuinely one thing this port computes and upstream computes another, but the
span max is closer than either alternative tried so far. The discriminating
measurement is cheap and was available all along: dump the actual `ttreat_p`
values from upstream (they are already in the capture's fifth column as `ti`; a
one-line change to the instrumentation writes `_e[2]`, the value itself) and
compare per entry. That is a second oracle rebuild, not a port change.

Interim state, measured and reproducible via
`macs-callpeak-e2e --peak-content /tmp/pc.tsv`: 9 of 12 summits byte-exact, 5
mismatches, peak counts and boundaries exact.

## F63 -- the summit computation is now fully characterised and reproduced exactly

The instrumentation was extended to write `ttreat_p` (`_e[2]`) and `tctrl_p`
(`_e[3]`) alongside the positions, and `CallPeakUnit` rebuilt again. Running
upstream's own summit argmax over its own captured data for peak 0 of
`se_model/realistic`:

```text
peak0 chunks: 67   max ttreat_p = 9.0   tied: 6
their (tstart + tend + 1) // 2: [7718, 7725, 7732, 7734, 7841, 7848]
midindex = (6 + 1) // 2 - 1 = 2  ->  7732
upstream summit for peak 0     =  7732      EXACT MATCH
```

So every element of the summit computation is confirmed against upstream's own
data:

| step | status |
|---|---|
| chunk list | upstream's, from the paired position array's `above_cutoff` subsequence |
| `tscore` | `treat_array[ti]` -- the paired treatment value at the chunk's own index |
| argmax | strictly-increasing scan, ties appended (`summit_value < tscore`, then `==`) |
| tie-break | **lower** median, `(len + 1) // 2 - 1` -- *not* `div_ceil` (F62) |
| midpoint | `(tstart + tend + 1) // 2` (F62) |

And the first captured `ttreat_p` values show what `ti` actually carries:

```text
(7571, 7590) ti=779  ttreat=3.0  tctrl=0.38573580980300903
(7590, 7596) ti=780  ttreat=4.0  tctrl=0.38573580980300903
(7596, 7604) ti=781  ttreat=4.0  tctrl=0.38573580980300903
```

Note the values are **not** monotone across the chunk and are far smaller than
the treatment pileup's run values later in the peak, confirming `ttreat_p` is a
value at one index of the *paired treatment* array, not a maximum over a span
(F47c) and not the value at the chunk's end (F62a). Both of those were measured
and refuted on this same data.

**What this settles.** With upstream's chunk list *and* upstream's `ttreat_p`
values, this port's summit arithmetic reproduces upstream's answer exactly. So
the entire remaining gap is upstream of the summit stage: this port must derive
the same chunk list from the paired position array. It is a construction
problem, fully specified, with the expected output now known for three fixtures.

**The measured reference for the next attempt** is simply: peak 0 of
`se_model/realistic` must yield 67 chunks with max `ttreat_p` 9.0 tied across 6
chunks at `7718, 7725, 7732, 7734, 7841, 7848`. A candidate implementation can be
checked against that in one run, instead of against an XLS whose summits were
previously wrong for unlocated reasons.

`oracle/dump_peak_content.py` plus `MACS3RS_PC=<file>` produces this capture for
any fixture; the format is now
`min_length, smoothlen, tstart, tend, ti, ttreat_p, tctrl_p`.

## F64 -- with upstream's own chunks AND `ttreat_p`, 7 of 12 summits are exact; two distinct residual causes

`macs-callpeak-e2e --peak-content` now takes `ttreat_p` straight from the
capture's sixth column, so the summit stage runs on upstream's own values.
Result on `se_model/realistic`: **12 peaks, 7 summits byte-exact, 5 mismatching**
(peaks 3, 5, 7, 8, 10).

Replaying the same data in Python over the capture, region by region:

```text
peak 30652-31142: chunks=80  max=9.0  tied=9  midindex->30991  (upstream 30992)
peak 45944-46527: chunks=80  max=8.0  tied=5  midindex->46194  (upstream 46195)
peak 61429-61884: chunks=61  max=9.0  tied=1  midindex->61743  (upstream 61743)  OK
```

Two **separate** causes fall out, which is why this is worth recording rather
than folding into F63:

**1. Rust's region grouping does not reproduce the capture's region.** Peak 7 has
a *single* tied maximum, so its summit is unambiguous -- Python gets 61743,
upstream reports 61743, and this port reports **61730**. There is no tie-break
that explains a wrong answer when there is only one candidate, so this port is
scoring a *different set of chunks*. Its `segment_regions` is grouping the 915
captured chunks with `max_gap = extsize = 200` rather than the region structure
upstream's own chunk loop produced. This is the dominant residual and it is a
grouping problem, not a summit problem.

**2. The lower-median index is off by one for odd tie counts.** Peaks 3 and 5
have 9 and 5 tied maxima, and `(len + 1) // 2 - 1` lands one short of upstream in
both. For an odd count that index is `len // 2`, i.e. the true middle -- so either
upstream's `tsummit` list has one *extra* entry that Python's filter is not
reproducing, or the ordering of tied entries differs. Given peak 7 (tied=1) is
correct, the ordering hypothesis is more likely than a plain off-by-one, and it is
cheap to test: print `mids` for peaks 3 and 5 in full and see whether upstream's
answer is `mids[4]`/`mids[2]` or the element after it.

Neither cause can be chased through the port until (1) is fixed -- with the wrong
chunk set, every comparison downstream is uninformative. So the ordering is:
reproduce upstream's **region grouping** first (using `min_length` and the
`max_gap` that `__chrom_call_treat_ctrl_a_chromosome` actually used, both of which
the capture already records in columns 1 and 2), then re-run this same comparison.

## F65 -- upstream forms 13 regions; F64a's contiguity inference was wrong

F64a inferred 22 regions from contiguity breaks in the capture. That was wrong,
and the way it was wrong is instructive. The instrumentation now writes a
per-region header line (a module-level counter, because a Cython extension type
cannot carry new attributes -- `self._macs3rs_region` silently never incremented
on the first attempt, which is itself a small trap worth noting).

Ground truth for `se_model/realistic`:

```text
regions: 13
chunks/region: [67, 88, 95, 80, 19, 56, 80, 68, 61, 76, 85, 64, 76]
sum: 915
```

**13 regions, 915 chunks**, of which 12 become peaks -- one region is dropped by
the `min_length = 200` length filter. That is the structure this port has to
reproduce, and it is now a hard, checkable target rather than an inference.

**Why the contiguity inference failed.** `__close_peak_wo_subpeaks` is called once
per region, but consecutive entries *within* a region do not satisfy
`tstart[i] == tend[i-1]`, because `ts = pos_array[ti-1]` and `te = pos_array[ti]`
and consecutive `above_cutoff` indices `ti` need not be consecutive array
indices. A break in contiguity therefore does **not** mark a region boundary, and
counting them over-counted. That is the same class of error as F64's 7-of-12: a
structural quantity inferred from a proxy instead of measured. The instrumentation
was already rebuilt twice; writing the region index directly is what it was for,
and the inference should not have been attempted before spending it.

**What this settles and what it does not.** It settles F64 cause 1's question --
this port must produce 13 regions with exactly these chunk counts, and can be
checked in one run. It does not by itself explain peak 7's summit: with upstream's
own region 8 chunks (index 8, 61 chunks -- which is exactly the count my Python
check found for peak 7's range) and upstream's own `ttreat_p`, the answer should
be unambiguous at 61743, and this port reports 61730. So either this port is not
splitting into those 13 regions, or `close_peak_wo_subpeaks` is being handed
something other than the region it thinks.

The next step is therefore narrow and mechanical: have `--peak-content` split the
capture on its `#region` header lines and emit **one peak per region**, bypassing
`segment_regions` entirely. If that yields 12 peaks with all summits exact, the
summit problem is closed and only the chunk *construction* remains; if not, the
fault is inside `close_peak_wo_subpeaks` rather than upstream of it.

Also recorded: `getattr(self, "max_gap", -1)` returns **-1** on
`CallerFromAlignments`, so `max_gap` is not carried on that object at all. It is
supplied by the caller as a plain local in `__chrom_call_treat_ctrl_a_chromosome`,
so it can only be captured by instrumenting *that* method -- another rebuild, and
one that is worth doing in the same pass as any further `peak_content` work.

## F66 -- summit parity achieved: 12/12 peaks, coordinates AND summits

F65's step, done. `--peak-content` now splits the capture on its `#region`
header lines and emits **one peak per region**, bypassing `segment_regions`
entirely, so this port sees exactly the regions upstream did:

```text
regions from capture: 12      (13 captured, 1 dropped by min_length = 200)
rust peaks:     12
upstream peaks: 12
CALLPEAK PEAK MATCH (12 peaks, coordinates + summits)
```

**Every peak matches on `start`, `end` and `abs_summit`.** For the first time
in this project, a `callpeak` fixture is at full parity with upstream including
summit coordinates.

That settles every stage of the narrow-peak path except one. The summit
selection, the lower-median tie-break, the `(tstart + tend + 1) // 2` midpoint,
the `min_length` filter and the region→peak emission are all **proven correct**
against upstream's own captured data. The only thing this port cannot yet
reproduce is the *construction* of the region list and the per-chunk `ttreat_p`
from its own signals -- and that is now a well-posed problem with a binary
pass/fail check:

    oracle/dump_peak_content.py --peak-content ...   ->  CALLPEAK PEAK MATCH

regression suite entry, since it is a fixed oracle artifact; and

* the reference chunk counts per region -- [67, 88, 95, 80, 19, 56, 80, 68, 61,
  76, 85, 64, 76], summing to 915 -- which the port must reproduce exactly before
  any summit is compared again.

**Why this matters beyond one fixture.** F66 converts the G9 summit problem from
"something in the summit code is wrong" into "the chunk list is wrong". Every
earlier hypothesis about the summit stage (F46, F47a, F47b, F50, F51a, F59-F63)
was probing the wrong layer, and the reason is now clear: the summit code was
right from F62 onwards, and only the input list was wrong. The `max_gap` capture
still outstanding (F65) matters only for the construction path, not for the
summit stage itself.

## F67 -- chunk construction: this port makes 522, upstream makes 915; the paired position array is the missing input

With F66 closing the summit stage, the residual is a single measurable gap.
Chunk counts for `se_model/realistic`:

```text
this port (from the coalesced q-track):   522
upstream (peak_content capture):           915   across 13 regions
```

So this port's chunk list is ~57% of upstream's -- it is **too coarse**, which is
the same direction as F59c originally suggested before that was withdrawn, and it
is now measured against the capture rather than inferred.

The cause is known and is a consequence of two things established earlier:

* the q-track **coalesces** (F59f: 4763 p-spans become 732 q-spans), because many
  distinct p-scores map to the same q-score and unmapped ones become 0;
* upstream's `above_cutoff` is `np.nonzero(score_array > cutoff)[0]` over the
  **position array**, not over runs (F45). Every position whose score clears the
  cutoff contributes one entry, so several entries can share a run's `(ts, te)`
  -- which is exactly what the capture shows (`ti=779 -> (7571,7590)`,
  `ti=780 -> (7590,7596)`).

So upstream's chunk count is the number of **above-cutoff positions**, and this
port's is the number of above-cutoff **runs**. Merging positions that share a
score is the entire deficit.

**What is needed to close it.** The chunk is `(pos_array[ti-1], pos_array[ti],
treat_array[ti])` for each above-cutoff `ti`, so the input this port lacks is the
**paired position array** itself -- not the q-track, not the p-track, and not the
`--bdg` (F54). It is built inside `__chrom_pair_treat_ctrl` (F35), which pairs
treatment and control and truncates at the paired extent, and it is `cdef`-private
(F60), so it cannot be read from Python.

That makes it one more instrumented capture: dump `pos_array`, `treat_array` and
`ctrl_array` from `__chrom_pair_treat_ctrl` -- the same one-line pattern as the
three `CallPeakUnit` dumps already in place, gated on the same `MACS3RS_PC`
environment variable. With that array the port can build chunks by index and
should reproduce 915 exactly, which `run_peak_e2e.sh` can then check directly.

This is the last unidentified input in the narrow-peak path. Everything from the
BED through the q-score is verified against upstream (F42, F47, F55, F59g), and
everything from the chunk list through the XLS row is verified (F66).

## F68 -- the paired position array is captured; 10,684 entries, and this is the last input

F67's requested capture is done. A dump was added to
`__chrom_pair_treat_ctrl` (the same one-line, `MACS3RS_PC`-gated pattern as the
three `CallPeakUnit` dumps) and the extension rebuilt:

```text
#paired  10684  10684  10684
paired array rows: 10684, positions 41..93003
first entries: (41, 0.0, 0.18901) (57, 0.0, 0.19094) (91, 0.0, 0.19094) ...
top treatment values (value -> first position):
   11.0 @ 84646   10.0 @ 38601   9.0 @ 7720   8.0 @ 7705   7.0 @ 7664
```

**Why this is the piece that was missing.** The paired array is
`(pos_array, treat_array, ctrl_array)` with one entry per run end, truncated at the
paired extent (F35) -- 10,684 entries ending at **93003**, which is exactly the
paired truncation point this port has been computing by hand all along (F42,
F52a). Everything downstream is now checkable against it directly:

* `peak_content` entries are `(pos_array[ti-1], pos_array[ti], treat_array[ti])`;
* `above_cutoff` selects indices of this array;
* the 915 chunks in 13 regions are a **subset** of 10,684 indices -- this port's
  522 is a coalesced approximation of that same selection.

**Correction to an overclaim in the first version of this finding.** It said the
paired array's treatment reaches higher than anything in the `--bdg`, citing
F54's 8.0. That was wrong: the 8.0 was from `spikes_only`, not this fixture.
Measured on `se_model/realistic`:

```text
max treatment in --bdg:        11
max treatment in paired array: 11.0
```

They **agree**. So for this fixture the `--bdg` treatment and the paired-array
treatment are the same signal, which is reassuring and consistent with F42's
green result -- and it means the `--bdg` differential is a valid check on the
treatment track after all. F54's separate conclusion, that the `--bdg` and the
*score* arrays differ for `spikes_only`, is unaffected and still stands; the two
fixtures were being conflated.

**Consequence for the construction path.** The port no longer has to *derive* the
chunk list from a coalesced track. It can build it by index:

```text
ti  in above_cutoff            (score > cutoff, and treat_array[ti] > 0)
ts = pos_array[ti-1],  te = pos_array[ti],  ttreat_p = treat_array[ti]
region breaks on tl = ts - lastp > max_gap
```

with `max_gap` still the one unrecorded input (F65) -- and `max_gap` can now be
*read off* the capture rather than assumed, by locating the positions where
upstream's own region boundaries fall.

This makes the remaining G9 work a construction exercise with a fixed reference
(10,684 array entries -> 915 chunks -> 13 regions -> 12 peaks, 67/88/95/80/19/56/
80/68/61/76/85/64/76 per region) and a binary pass/fail check already built
(`CALLPEAK PEAK MATCH`). Every stage either side of that boundary is verified.

## F69 -- the last structural difference: this port's RLE coalesces, upstream's pileup does not

Comparing the two arrays that the chunk construction runs on:

```text
upstream paired position array (F68):  10,684 entries, positions 41..93003
this port's retadd walk (F59f):          4,763 spans
```

A 2.2x deficit, and the cause is a single line in `macs-rle`:

```rust
pub fn push_run(&mut self, run: Run<T>) {
    ...
    if let Some(last) = self.runs.last_mut() {
        if last.value == run.value { last.end = run.end; return; }   // <-- coalesce
    }
```

Upstream's `se_all_in_one_pileup` emits a position at **every** read start and
read end, whether or not the depth changes, so its pileup arrays retain every
endpoint. The capture shows this directly -- consecutive entries carry identical
depths:

```text
(41, treat=0.0, ctrl=0.18901)
(57, treat=0.0, ctrl=0.19094)
(91, treat=0.0, ctrl=0.19094)      <- ctrl identical to the previous entry
(117, treat=0.0, ctrl=0.19287)
```

Repeated `ctrl` values across consecutive entries are only possible if
redundant endpoints are kept. This port drops them, so its paired array has fewer
indices, so its `above_cutoff` selection has fewer candidates, so it produces
**522 chunks where upstream produces 915** (F67).

**So the construction path needs a non-coalescing pileup.** Not a change to how
values are computed -- those are verified (F42, F68) -- but a change to how
breakpoints are recorded: the pileup that feeds the score and chunk stages must
retain every endpoint, while the RLE `SignalTrack` stays coalesced everywhere
else it is used (bedGraph output, local-lambda merge, peak reporting). Those are
genuinely different requirements and conflating them is what produced the 522-vs-
915 gap.

That is the last structural difference in the narrow-peak path. The full chain
now reads:

| stage | status |
|---|---|
| BED parse, dedup | verified (1492 files, F19/F38) |
| treatment pileup values | verified (F42, and agrees with the paired array, F68) |
| control local lambda | verified (F42) |
| p-score, p->q walk | verified (F47, F52, F55, F59g) |
| **paired array breakpoints** | **this is the gap (F69)** |
| region split, `min_length` | verified (F65, F66) |
| summit argmax, tie-break, midpoint | verified (F63, F66) |
| XLS emission | verified (F66, byte-identical on all 12) |

### F69a -- correction: the pileup does NOT coalesce away endpoints; F69's diagnosis was wrong

F69 claimed `SignalTrack::push_run`'s coalescing was dropping redundant endpoints
and shrinking the paired index space. It was tested rather than assumed, and it
is **not** the cause:

* `push_exact` was added to both `SignalTrack` and `TrackBuilder` (kept -- it is a
  legitimate, documented alternative, just not the bug);
* the pileup sweep in `macs-pileup` was switched to it;
* result: **`p-spans` stayed at 4763, unchanged.**

The reason is straightforward in hindsight: every event in the endpoint sweep
changes the depth, so consecutive pushes never carry equal values and the
coalescing branch never fires. There was nothing to coalesce away.

So the 4763-vs-10684 gap has another cause, and the arithmetic in F68/F69 does not
account for it. The counts do not fit the "union of two pileups' run ends" model
at all: `se_model/realistic` retains **177 treatment reads and 1989 control
reads** after deduplication, so even counting every read start *and* end as a
distinct position gives on the order of 2 x (177 + 1989) = 4332 -- not 10,684.
Something upstream does emits roughly 2.5x more positions than the deduped reads
can account for, which means the paired arrays are built from a **larger, earlier**
signal than the deduplicated track, or they retain non-deduped structure.

That is the honest state of it: a measured, reproducible 2.2x discrepancy whose
cause is not identified, with the specific arithmetic that the obvious
explanation fails to satisfy. The `push_exact` additions stay (they are correct
and documented), but they are **not** credited with fixing anything.

Worth stating as a rule alongside F59c, F64a and F68: three of the last four
diagnoses were wrong on first formulation. The pattern is consistent -- reach for
the nearest structural explanation, then measure it before writing it up. The
measurement here took one build; the claim took longer than it should have.

## F70 -- the 10,684 positions are the union of three control *scales*, not read endpoints

F69a could not explain 10,684 positions. The arithmetic now pins it.

`se_model/realistic` has **192 treatment reads** and **2000 control reads**
(`total tags in treatment: 192`, `total tags in control: 2000`, both unchanged by
dedup). Every read contributes at most two endpoints, so the *read-endpoint* model
can produce at most `2 * (192 + 2000) = 4384` positions. **10,684 is well above
that**, which is why F69a's reasoning failed.

But the control lambda is not one pileup. `as_pairs` builds **three** scales --
`d = 200`, `slocal = 1000`, `llocal = 10000` (F39) -- each a separate pileup of
the same 2000 reads at a different extension length, hence with breakpoints at
*different* positions. `pileup_a_chromosome_c` folds them with
`over_two_pv_array(..., func="max")`, whose walk advances only the lagging
pointer and therefore emits close to the **union** of the three scales' run ends
(F47a/F51a established that walk is faithful).

That gives a workable estimate: one scale over 2000 reads yields on the order of
3.5k distinct event positions, and three scales unioned give roughly
`3 x 3500 ~ 10,500` -- the same order as the observed 10,684.

**So the paired array's control side is the merged three-scale lambda, and its
breakpoint count is set by the merge, not by the reads.** That is consistent with
everything else: the treatment side is tiny (192 reads), the treatment values in
the capture are small integers, and the `--bdg` control file -- which *is* the
merged lambda -- matched this port's merge to 5e-6 (F42) while the *breakpoint
count* was never checked, because the bedGraph coalesces at 1e-5 (F42) and hides
it. F42's caveat, written when the summit work started, applies here.

The next measurement is therefore specific and cheap: count the runs in the
**merged control track** (not the treatment, not the p-score track) and compare
with 10,684. If it is short, the merge walk is dropping breakpoints the way
`over_two_pv_array` drops tails -- and that is an upstream behaviour already
confirmed present in both implementations (F47a), so the fix would be in how the
scales are folded, not in the walk.

## F71 -- F69 was right about the mechanism and wrong about the call site; reconciled

Measured on `se_model/realistic`:

```text
treat_runs = 376   merged_ctrl_runs = 4755   sum = 5131
upstream paired array = 10684
```

The merged control dominates (4755 of 5131), confirming F70's direction. And the
missing factor is ~2.2x.

**Checking whether the merge itself coalesces:** it does not.
`over_two_pv_array_track` builds via `track_from_pv`, which pushes `Run::new`
directly into a `Vec` -- no `push`, no coalescing. So every position the walk
emits is kept, exactly as upstream's raw-array `over_two_pv_array` does (F51a).
The merge is exonerated.

**Which puts the coalescing upstream of it -- in the per-scale pileups.** And
that is exactly where F69 pointed, but F69a tested the wrong function:

* `push_exact` was applied to `macs-pileup`'s endpoint sweep and changed nothing
  (F69a). That sweep builds the **directional** treatment pileup, where every
  event changes the depth by +/-1, so consecutive values never match and the
  coalescing branch never fires. The test was valid and proved only that the
  treatment path was never the problem.
* the **control** pileup is `bidirectional` (`SingleEndParams::bidirectional`).
  There, a plus read contributes events at `p-d` and `p`, and a minus read at `m`
  and `m+d`; with overlapping reads a start and an end can land on the same
  position or net out, so **consecutive depths can be equal** and coalescing
  does fire. That is why the treatment keeps 376 runs while the merged control
  collapses to 4755 instead of the ~10,300 the three scales imply.

**So the correction is narrow and well-targeted:** build the *control* (local
lambda) pileup with `push_exact`, and leave the treatment directional pileup
coalescing as it is. The two have genuinely different event structure and were
wrongly sharing one code path's default. `push_exact` is already added to both
`SignalTrack` and `TrackBuilder` (F69a) and is unused -- this is its purpose.

The arithmetic to check against: three control scales over 2000 reads, unioned,
should reach the 10,684 upstream reports; this port currently reaches 4755 for
the merged control against ~10,300 expected. If `push_exact` on the bidirectional
path closes that, F67's 522-vs-915 chunk deficit should close with it.

**Note on the three findings this reconciles.** F69 named the right mechanism and
pointed at the wrong function; F69a refuted it by testing only the function where
the mechanism cannot act; F70 correctly located the deficit in the control side by
arithmetic. None of the three was wrong on its own evidence -- the error was in
generalising from one measurement to a whole code path. That is a sharper failure
than "guessed too fast", and the fix is the same discipline: measure each call
site before generalising from it.

### F71a -- the non-coalescing control pileup changes nothing either; F71 is unconfirmed

F71's proposed fix was implemented rather than argued: `SingleEndParams` gained a
`coalesce` flag, set `false` by [`SingleEndParams::bidirectional`] and `true` by
[`SingleEndParams::directional`], and `sweep_to_track` uses it to choose between
`TrackBuilder::push` and the `push_exact` added in F69a. The plumbing is correct
and the intent is documented, but the measurement is null:

```text
merged control runs, before: 4755
merged control runs, after:  4755
p-spans, before and after:   4763
```

So the bidirectional control pileup also has nothing to coalesce away, and F71's
mechanism -- like F69a's -- does not account for the 10,684-vs-4755 gap. The
`coalesce` flag and `push_exact` are kept because the distinction between the two
pileup kinds is real and now explicit in the API, but neither is credited with
fixing anything.

**Where that leaves the deficit.** Three separate coalescing explanations have now
been tested and refuted (F69 on the treatment sweep, F71a on the control sweep,
F47a on the merge itself). The merged control has 4755 runs and upstream's paired
array implies ~10,300, and **nothing measured so far explains the factor of ~2.2**.

The one hypothesis not yet excluded is structural rather than numerical: that the
paired array is not built from the merged lambda's run ends at all, but from the
**concatenation** of the per-scale arrays (`d`, `slocal`, `llocal`) rather than
their pointwise maximum -- which would give roughly `3 x 3500 ~ 10,500` directly.
`as_pairs` returns a list of per-scale tracks, and `__chrom_pair_treat_ctrl`
receives that list; if the pairing walks *all* of them rather than the merged one,
the arithmetic fits exactly and no coalescing question arises. That is a question
about which array is passed, and it is answerable from the capture alone by
checking whether the captured `ctrl` values ever exceed the per-scale maximum
implied by the merge.

Until that is checked, this remains the one open item in the narrow-peak path,
and it should be recorded as open rather than explained.

## F72 -- the paired array's control side is NOT the merged lambda; only 58.8% of positions agree

F71a's structural hypothesis, tested directly. This port's merged control track
was dumped (`CALLPEAK_CTRLDUMP`) and compared, position by position, against the
control values in upstream's captured paired array (F68):

```text
paired positions compared: 10684 of 10684
control value matches: 6285 (58.83%), worst abs diff 1.929e-01
```

**41% of positions disagree, by up to 0.19 in absolute lambda.** So the local
lambda that feeds the p-score is materially different from the merged
three-scale lambda this port computes -- and that is the cause of everything still
unexplained:

* different lambda -> different `get_pscore` -> different `above_cutoff` set ->
  **522 chunks where upstream has 915** (F67);
* and it makes sense of the earlier arithmetic puzzle (F70/F71). The 10,684
  breakpoints are not the merged lambda's ~4,755 runs at all; they belong to a
  *different* track with ~2.2x more of them. No coalescing explanation was ever
  going to produce it, which is why F69, F69a and F71a were all refuted.

**And this exposes a gap in a gate that has been green all along.** F42 compares
the `--bdg` control lambda and reports a match to 5e-6 across all 97 corpus
fixtures. That is still true -- but it is a *different* control array from the one
the p-score is computed on. This is the same divergence F54 found on the treatment
side, now confirmed on the control side: the `--bdg` and the paired arrays are
distinct representations, and only the former was ever compared.

So there are two facts that have to be held together, and only one was checked:

| array | status |
|---|---|
| `--bdg` control lambda | verified, 97/97, <=5e-6 (F42) |
| **paired-array control** | **58.8% match, worst 0.19 (this finding)** |

**What is needed next** is concrete and checkable from data already captured: the
per-scale arrays (`d`, `slocal`, `llocal`) as upstream actually builds them, to
determine whether the paired control is a concatenation of the scales rather than
their pointwise maximum. That was F71a's structural guess; this finding makes it
the leading explanation rather than one of several, because a 0.19 discrepancy and
a ~2.2x breakpoint ratio are both consistent with a different array being passed
to `__chrom_pair_treat_ctrl`.

Worth recording as a process note: F42's green gate was never wrong about what it
measured. It measured the `--bdg`. The error was reading that as coverage of the
scored arrays, which is the same inference error as F52's pseudocount and F64a's
region count -- three times now, a single pattern: **a green gate on one
representation is not evidence about another.** The durable fix is to capture
every array a stage actually consumes, not the ones that are convenient to write.

### F72a -- the discrepancy is symmetric, so it is not a baseline or a scale offset

The sign of the 4,399 mismatches at F72:

```text
mine - upstream:  min -0.19287   max +0.19287   mean +0.00005
all mine > upstream: False        all mine < upstream: False
```

A constant offset is ruled out immediately, which matters because the most
tempting explanation was `lambda_bg`: this port applies
`SingleEndParams::with_baseline(lambda_bg)` to every control scale, and if
upstream applied it to only some, or with a different `lambda_bg`, the difference
would be one-signed. It is symmetric about zero instead.

The symmetric spread also bounds it: `|diff| <= 0.19287` exactly at both ends,
which is the maximum local control value on this fixture -- i.e. at the mismatching
positions the two tracks are reading *different runs* of a piecewise-constant
function of similar range, not one being a constant shift of the other.

**Two readings remain, and they are distinguishable with data already captured:**

1. **Breakpoint misalignment.** This port's merged control has 4,755 runs and
   upstream's paired array 10,684. If the underlying function is the same but the
   breakpoint sets differ, then at a coordinate that is a breakpoint on one side
   and interior to a run on the other, `value_at` can legitimately return
   different values -- and symmetric disagreement is exactly what that produces.
   *Test:* compare the two as **functions**, refining both on the union of their
   breakpoints, the same way `run_e2e.sh` compares bedGraphs (F42). If they agree
   under that comparison, the values are fine and only the breakpoint sets differ.
2. **Genuinely different functions.** The per-scale arrays differ, so the
   lambda that feeds the p-score is not this port's pointwise max.

Reading 1 is cheaper and would explain why 58.8% agree (the coordinates where
both tracks happen to be interior to the same-valued run). It also predicts the
~2.2x breakpoint ratio independently of the value question, which is a strong
joint fit.

The function-level comparison is the right next move and reuses machinery that
already exists and is already trusted (`run_e2e.sh`'s interval refinement).

## F73 -- the merged control lambda's VALUES are correct; only the breakpoint set differs

F72a posed two readings and named the cheaper test. Run: build a step function
from the paired array using the end-indexed convention (F5, `ctrl[i]` is in
effect over `(pos[i-1], pos[i]]`) and compare it against this port's merged
control on the union of both breakpoint sets -- the same interval refinement
`run_e2e.sh` already uses for bedGraphs (F42).

```text
function-level: 9441 intervals compared, 0 differ, worst 4.932e-10
-> values AGREE (breakpoints only differ)
```

**The merged three-scale local lambda is correct to f32 precision.** After 89
findings and a long series of refuted coalescing hypotheses, the lambda was never
the problem.

What remains is exactly one thing, and it is now stated without ambiguity:

> this port's merged control has **4,755 runs**; upstream's paired array has
> **10,684** breakpoints. The *function* is identical to 5e-10; only the
> breakpoint set is smaller.

That also retro-explains every refuted hypothesis at once. F69, F69a and F71a all
attacked *values* or *coalescing* and all measured no change, because no value was
ever wrong. The 522-vs-915 chunk deficit (F67) is entirely a consequence of the
breakpoint count: fewer breakpoints means a smaller index space for
`above_cutoff`, which selects *positions*, which is why upstream's
`above_cutoff` yields 915 chunks from 10,684 indices and this port's yields 522
from 4,763 (F59f).

**So the remaining G9 task is a single, mechanical one:** make the merged control
retain the same breakpoints upstream does. The function is already right, so this
is not a numerical problem -- it is finding which intermediate representation
carries the extra ~5,900 breakpoints. Since the merged control is the pointwise
max of three per-scale tracks and `over_two_pv_array` does not coalesce (F71),
the candidates are now narrow and testable:

* the per-scale tracks are built with fewer breakpoints than upstream's
  (`pileup_a_chromosome_c` writes raw arrays with one entry per event), so the
  fix may be that the per-scale pileups must keep every endpoint -- which is
  `push_exact`, already implemented and *already applied* to the bidirectional
  path (F71a) but measured as making no difference to the merged count. That null
  result is now worth re-examining: the function being identical means the
  per-scale breakpoints *are* present but being merged away somewhere that
  preserves values, which is precisely what `over_two_pv_array` + `track_from_pv`
  do not do -- so the loss must be in the per-scale pileups, and F71a's null
  measurement should be re-examined rather than trusted.

Recorded as open, with the discriminator: compare this port's *per-scale* control
track run counts against upstream's, which needs one more capture of
`pileup_a_chromosome_c`'s individual `prev_pileup` before the merge.

## F74 -- per-scale control tracks measured: the loss is in the per-scale pileup, not the merge

`pileup_a_chromosome_c` now dumps each scale's track *before* folding
(F73's next step). On `se_model/realistic`:

```text
upstream per-scale, d      = 3878 entries
upstream merged after d    = 3878
upstream per-scale, slocal = 3879
upstream merged after both = 7622
upstream per-scale, llocal = 3804
upstream merged after all  = 11085      <- vs this port's 4755
this port  merged control  = 4755
```

**This localises the deficit precisely and exonerates two more suspects.**

*The merge is a union, and this port's merge is faithful.* Upstream's fold goes
3,878 -> 7,622 -> 11,085, i.e. close to the **sum** of the inputs, which is the
lagging-pointer union walk of `over_two_pv_array` (F47a/F51a) -- and this port's
`over_two_pv_array` uses the same `Ordering::Less` / `Greater` / `Equal` advance
(verified by reading it alongside the capture). So the merge is not losing
breakpoints.

*The loss is upstream of the merge, in the per-scale pileups.* For this port's
merged 4,755 to come out of a three-way union, its per-scale tracks must be
roughly 1,600 entries each -- against upstream's ~3,878. That is the ~2.4x, and
it is the same factor as the 10,684-vs-4,763 deficit seen at F59f/F67, which
means both observations have the **same single cause** in the per-scale
bidirectional pileup.

**And this makes F71a's null result the thing to re-examine, not trust.** F71a
applied `coalesce = false` to the bidirectional path and measured no change in
the *merged* count. But the merged count is now known to be a poor proxy: the
merge is a union, so extra breakpoints in one scale propagate. A null there is
consistent with the flag simply not reaching the code that builds the per-scale
tracks. The direct measurement now available -- compare this port's per-scale run
counts against 3878 / 3879 / 3804 -- tests the same hypothesis without going
through the merge at all, and is the correct thing to instrument next.

That is the tightest the gap has been stated anywhere in this investigation: a
named function, three target numbers, and a measurement that does not depend on
any downstream stage.

## F75 -- ROOT CAUSE: `SignalTrack::from_runs` coalesces, and the merge goes through it

F74's measurement first, because it exonerates the pileups outright. This port's
per-scale control tracks, measured against the capture:

```text
this port:  d=200 runs=3878   slocal=1000 runs=3879   llocal=10000 runs=3804
upstream:   d       3878               3879                      3804      <- exact match
```

**The per-scale pileups are bit-for-bit correct.** That was never the problem.

The merge is where it goes wrong:

```text
upstream fold:  3878 -> 7622 -> 11085      (union: close to the sum of inputs)
this port fold: 3878 -> 4566 ->  4755      (drops ~2,900 positions)
```

And the cause is in `macs-rle`:

```rust
// SignalTrack::from_runs
if let Some(last) = t.runs.last_mut() {
    if last.value == run.value {
        last.end = run_end;      // <-- coalesces
        prev_end = run_end;
        continue;
    }
}
```

`track_from_pv` builds its `Vec<Run>` with `runs.push(...)` -- which looked like
no coalescing -- but it hands that `Vec` to `SignalTrack::from_runs`, and
**`from_runs` coalesces equal-valued neighbours**. Upstream's
`over_two_pv_array` writes into raw numpy arrays with no such step, so it keeps
every one of its ~7,622 emitted positions.

**This is the whole deficit.** It explains simultaneously:

* the merged control at 4,755 instead of ~11,085 (F74);
* the paired array at 4,763 instead of 10,684 (F59f/F67);
* the chunk count at 522 instead of 915 (F67) -- fewer breakpoints, smaller
  `above_cutoff` index space;
* every summit that came out low, because the selected chunk is chosen from a
  different breakpoint set.

**And it corrects F71**, which exonerated the merge on the strength of reading
`track_from_pv` in isolation. The coalescing is one call deeper, in the
constructor. That is the same mistake as F52 (verifying a convention against a
function the pipeline never calls) and F64a (inferring a structural quantity from
a proxy): **reading one link of a chain and concluding about the whole.** It has
now happened three times in this investigation, which is enough to treat "I read
the function" as insufficient and require running the pipeline instead.

**The fix is small and needs no new oracle capture:** give `from_runs` a
non-coalescing sibling (the `push_exact` added in F69a is the same idea at the
`push` level) and have `track_from_pv` use it, leaving the coalescing form for
bedGraph output and the local-lambda *reporting* path where it is wanted. The
acceptance check already exists and is exact: per-scale counts stay
3878/3879/3804, the fold should reach ~11,085, the paired array ~10,684, and
`oracle/run_peak_e2e.sh` should then show summits matching on the normal path --
which no configuration has achieved yet.

## F76 -- FIXED: `from_runs_exact`; the local-lambda fold now matches upstream's count exactly

F75's root cause is fixed. `SignalTrack::from_runs_exact` is the coalescing
`from_runs` with value-merge omitted (clamping and zero-width/backwards-run
dropping still apply, so malformed input still cannot produce a corrupt value
function). `track_from_pv` now uses it; the coalescing form is untouched for the
bedGraph and local-lambda reporting paths, whose upstream counterparts do
coalesce.

Every number predicted by F74-F75 now holds:

```text
per-scale control tracks   3878 / 3879 / 3804    unchanged, still exact
two-scale fold             3878 -> 7622          was 4566   <- now exact
three-scale fold           3878 -> 7622 -> 11085 was 4755   <- now exact
```

The three-scale figure is upstream's own captured value from F74
(`#merged 3 11085`), so the merged local lambda is now breakpoint-for-breakpoint
the same size as upstream's. The 10,684-entry paired array is that same track
truncated to the paired treatment/control extent, so the two figures are
consistent rather than competing.

**Peak boundaries were already exact and remain exact** on
`se_model/realistic` (7572-8039, 15159-15802, 22782-23497, 30653-31141,
38462-38772, 45945-46526, ...), which is consistent: boundaries are chosen from
the p-value crossing of a lambda that was already numerically right, so only the
breakpoint set mattered for them.

**Summits improved sharply but are not yet exact** -- now within 3-18 bp
(7726 vs 7732, 15495 vs 15506, 23038 vs 23056, 30980 vs 30992, 38607 vs 38610,
46181 vs 46195) where they had been far worse. That is the expected shape of the
result: summit selection picks *which* chunk in a peak is the summit, so it is
the stage most sensitive to the position set, and it is reading the now-correct
one but still forming chunks differently.

**What remains is the chunk construction, and it is now a well-posed problem.**
F66 already established that feeding upstream's chunk boundaries in via
`--peak-content` reproduces summits exactly, so the lambda track is no longer a
suspect at all. The gap is entirely how `peak_content` pairs top positions into
chunks, applies `min_length`, and emits peak edges from the paired array. The
next step needs no new oracle capture: the F66 capture already contains the
upstream chunk boundaries and the merged local lambda, and this port now produces
the same merged lambda.

## F77 -- the harness merges above-cutoff positions into one chunk; upstream emits one per position

With the lambda fold fixed (F76), peak boundaries are exact and summits are
within 3-18 bp. Reading the two remaining candidate stages against upstream rules
both out as *correct*:

* `segment_regions` splits on `c.start - last_end <= max_gap`, which is exactly
  upstream's inline `pre_p - peak_content[-1][1] <= max_gap` in
  `naive_call_peaks` (`PileupV2.py:1119`) -- same operands, same order.
* The summit loop resets the candidate list on a strictly-greater pileup, appends
  on equality, and takes the lower median `(len + 1) // 2 - 1` -- identical to
  `__close_peak` (`PileupV2.py:1067-1074`).

So the divergence is in how the chunk *list* is built, and it is in the harness,
not the library. `macs_peaks::caller::call_peak` builds chunks from the
**position list**:

```rust
let chunks: Vec<Chunk> = above.iter().map(|i| Chunk {
    start: if *i == 0 { 0 } else { paired.pos[i - 1] },
    end:   paired.pos[*i],
    ...
```

-- one chunk per above-cutoff position, exactly like upstream's
`peak_content.append((pre_p, p, v))`.

`bin/callpeak_e2e.rs` instead iterates RLE **runs** and pushes one chunk per
above-cutoff run. Its own comment (F57) states the requirement correctly --

> do NOT merge contiguous above-cutoff runs ... a new peak region is started only
> by `ts - lastp > max_gap`, not by run adjacency

-- but the code does merge, because it emits one chunk per `r` in the run scan.
With F76's fix the *lambda* track now has one breakpoint per paired position, but
the scan in question walks the **score** track, whose runs can still span several
above-cutoff positions. Every position inside such a run that upstream would have
made its own chunk is silently collapsed into one.

**That is exactly the shape of the residual error.** Region boundaries depend only
on the *first* and *last* chunk, so they stay exact -- which is why F76's
boundaries match. The summit depends on the *whole* chunk list, since it takes a
median over the midpoints of the chunks tied at the maximum pileup. Collapsing k
interior chunks into one replaces k distinct midpoints with a single wide-chunk
midpoint, shifting the median by a few bp. That is the 3-18 bp residual, and it
explains why it is one-sided and why it survived a lambda fix that was otherwise
exact.

**The fix is to build chunks from the position list**, as `caller.rs` already
does. The better fix, and the one to prefer, is to **delete the harness's inline
peak-caller and call `macs_peaks::caller::call_peak`**: the duplicate is the only
reason the harness could disagree with the library, and F77 is the third bug in
this investigation that existed solely because a second implementation of
already-correct logic was maintained elsewhere. Verification is exact and
available: chunk count must reach 915, region count 13, and
`oracle/run_peak_e2e.sh` must show exact summits on the **normal** path with no
`--peak-content` -- the first configuration ever to do so.

## F78 -- the refactor target is named, and the paired position array is already available

F77's fix, specified concretely against the current tree.

**The position array is already in hand.** The merged local lambda's runs *are*
the paired array's positions: F76 measured the fold at 11,085 and F73's capture
put upstream's paired array at 10,684 with positions `41..93003` -- the same
breakpoints, truncated to the paired extent. `PairedSignal`
(`crates/macs-peaks/src/driver.rs:124`) is exactly
`{ pos, treat, ctrl }` in the end-indexed `[p, v]` form, with `pos[i]` covering
`[pos[i-1], pos[i])` (F5), so it can be built directly from the merged lambda
plus the treatment track with no new arithmetic.

**Where the two code paths diverge, precisely:**

* `regions.rs::call_peaks_chromosome` (line 376) takes `chunks` as an argument
  and is faithful -- it segments, closes peaks, and applies `min_length`,
  `pmax`, and the summits path.
* `caller.rs::call_peak` builds `chunks` from `PairedSignal` via the **position
  list**, one entry per above-cutoff position, which is upstream's semantics.
* `bin/callpeak_e2e.rs` builds its own `chunks` by scanning RLE **runs** and
  taking the run maximum, which is not upstream's semantics.

So the divergence is one function -- the harness's chunk loop -- and the fix is
to build `PairedSignal` from the merged lambda and the treatment track, then
construct chunks from `above_cutoff` indices exactly as `caller.rs` does, and
pass them to the same `call_peaks_chromosome` the harness's `--peak-content` path
already reaches.

**Two details that will bite if they are missed:**

1. The per-criterion score arrays passed to `call_peaks_chromosome` are indexed
   by `Chunk::score_index`, which under the new construction is a **position
   index**, not an RLE index. The harness currently passes RLE-backed slices for
   the q-score track; it must pass per-position arrays, or the summit's double
   cutoff re-check reads the wrong entries.
2. Upstream's summit value for a position is `treat_array[ti]` -- the treatment
   **at that position**. Taking the run maximum (what the harness does now, and
   what F47c's comment warns about for a different reason) is a different number
   whenever a run spans more than one position, and it changes which chunks tie
   at the maximum.

**Verification, unchanged and exact:** 915 chunks, 13 regions, 12 peaks, and
`oracle/run_peak_e2e.sh` reporting summits matching on the **normal** path with no
`--peak-content` supplied.

## F79 -- position-list chunks: summits close dramatically, boundaries regress; the two conventions differ

F78's change is implemented. The harness now builds `qpos`/`tpos`/`cpos` over the
merged lambda's runs (F76: 11085 positions, upstream's own fold size), emits one
chunk per above-cutoff position with `start = pos[i-1]`, `end = pos[i]`,
`treat = treat[pos[i]]`, and passes a **per-position** score array to
`close_peak_wo_subpeaks` so `Chunk::score_index` indexes positions (F78 trap 1).

**Summits improved substantially**, which confirms the diagnosis:

```text
summit error, peak 2:  18 bp -> 3 bp      (23053 vs 23056)
summit error, peak 3:  12 bp -> 1 bp      (30991 vs 30992)
```

**But peak boundaries regressed by a few bp**, from exact to:

```text
peak 0: 7566-8035 vs 7572-8039    was exact before this change
peak 1: 15156-15791 vs 15159-15802
```

So the two conventions are genuinely different and I had them entangled. The
run-based loop used `start = prev + 1` -- the first base of the containing
interval -- and got boundaries exact. The position-based loop uses
`start = pos[i-1]`, which is what upstream's `naive_call_peaks` writes into
`peak_content`, and moves boundaries 3-6 bp early.

Reading upstream's `naive_call_peaks` more carefully explains why that is not
contradictory: there, `pre_p` is the previous **array position**, and because the
loop assigns `pre_p = p` on *every* iteration including below-cutoff ones,
`pre_p` is the last array position before the first above-cutoff one -- which for
a below-cutoff position `p` is *also* the end of the interval that begins at
`pos[i-1]`. Upstream's `peak_content` therefore records `ts` values that are
**lambda breakpoints**, whereas this port's run-based `lo = prev + 1` is the
first base *inside* the interval. The two differ by however many positions were
below cutoff immediately before the peak, which is exactly the few-bp scale
observed.

**The conclusion is that one convention cannot serve both stages**, which is
worth stating because it is the sort of thing that invites a "just use one" fix:

* a peak's **edges** come from the lambda breakpoints bounding the first/last
  above-cutoff position;
* a peak's **summit** comes from `(tend + tstart + 1) // 2` per position, where
  `tstart`/`tend` are those same breakpoints.

So both should use the same breakpoints -- which means the regression is not a
convention problem but evidence that `pos[i-1]` is being computed against the
**lambda** breakpoints while the correct reference for `ts`/`te` may be the
**paired** breakpoints, and those differ by the 11085-vs-10684 truncation (F76).
The discriminating measurement: compare this port's first above-cutoff position
against upstream's captured `peak_content` first-chunk `tstart` for each of the
13 regions; if they agree, the offset is in `pos[i-1]`'s indexing, and if they
differ by the truncation, the lambda track is being used where the paired array
is required.

**No configuration yet has exact boundaries and near-exact summits together**, so
this is progress but not a fix. The prior state (boundaries exact, summits
3-18 bp) and this state (boundaries 3-6 bp early, summits 1-3 bp) trade off along
one axis, and resolving it needs the discriminating measurement above rather than
another attempt at a single convention.

## F80 -- ROOT CAUSE: the paired position array is the union of treatment AND control breakpoints

F79's discriminating measurement, run against the F73/F66 captures.

Upstream's first chunk of region 1 is `ts=7571 te=7590 ti=779 tp=3.0`, then
`ts=7590 te=7596 ti=780 tp=4.0` -- two chunks, so **7590 is a breakpoint**.
This port produced `ts=7566 te=7596 ti=776 tp=4` then `ts=7596 te=7604 ti=777` --
one 30 bp chunk where upstream has 19 bp and 6 bp, so this port is **missing a
breakpoint at 7590**.

The paired array settles what it should be:

```text
paired positions near 7560-7610: [7566, 7571, 7590, 7596, 7604, 7608]
```

Every one of upstream's chunk boundaries is in the paired array, **including 7571
and 7590 which this port's chunks skip.** So the paired array is a strict
superset of the local-lambda breakpoints in this region, and `s.ctrl.runs()` --
which is what F79 used as the position list -- is not it.

**Why: the lambda is the union of the *control* per-scale tracks only.** F74
measured those at 3878/3879/3804 and F76 measured their fold at 11085. But
`__chrom_pair_treat_ctrl` pairs the *treatment* track against the control lambda,
so its positions are the union of **both** tracks' breakpoints. Positions like
7571 and 7590 are treatment-pileup breakpoints that carry no control breakpoint,
which is exactly why the lambda is missing them and exactly why they still appear
in `peak_content`.

**This also retracts the significance of F76's count match.** This port's fold is
11085 and upstream's captured fold is 11085, but that agreement is a
**coincidence of totals**: this port is short by the treatment-only breakpoints
in some regions and correspondingly long elsewhere. The value function still
agrees to 5e-10 (F73) because coalescing a breakpoint changes the *breakpoint
set*, not the function -- which is the same distinction F73 drew, and the reason
a count match was never going to be sufficient evidence. Counting equal totals
across differently-placed breakpoints is the weakest kind of check, and it is the
one F76 leaned on.

**The fix:** build the position list as the union of the treatment track's
breakpoints and the merged lambda's breakpoints, sorted, rather than using the
lambda alone. `treat.value_at` / `ctrl.value_at` then sample both tracks at each
union position, which is precisely what `__chrom_pair_treat_ctrl` produces. That
also subsumes F78's `score_index` concern: positions and scores become the same
index space again.

**Verification:** my first chunk must become `ts=7571 te=7590 ti=779 tp=3.0`,
matching upstream's capture line for line, and the region count must reach 13.

## F81 -- F80's union works, and exposes a second, independent off-by-one

F80's fix is implemented: the position list is now the sorted union of the
treatment track's and the merged lambda's breakpoints, sampled at each.

**It reproduces upstream's chunk line for line.** Upstream's captured first two
chunks of region 1 versus this port's:

```text
upstream  ts=7571 te=7590 ti=779 tp=3.0   |  ts=7590 te=7596 ti=780 tp=4.0
this port ts=7571 te=7590 ti=779 tp=3     |  ts=7590 te=7596 ti=780 tp=3
```

`ts`, `te`, `ti` and the first `tp` all match exactly. Before F80 the same
comparison produced `ts=7566 te=7596 ti=776` -- wrong on every field.

**Getting `tp` right required a second fix**, independent of F80. The paired
array is end-indexed (F5): `treat_array[ti]` is the value over
`[pos[ti-1], pos[ti])` -- the run whose end **is** `pos[ti]`. This port sampled
with a point lookup at `p`, which reads the interval to the *right* of `p` and
shifts every treatment value by one run. Sampling instead by "the run with the
largest end <= p" -- which coincides with the run ending at `p` when `p` is a
breakpoint, and with the covering run to its left when `p` is a
treatment-only position -- reproduces `tp=3.0` exactly.

The fallback matters: returning `0.0` for treatment-only positions (where the
control has no run ending) makes the p-score lambda non-positive and panics
`pscore.rs` with `NonPositiveLambda(0.0)`. The control value over
`[prev, p)` is non-zero, which is precisely what `__chrom_pair_treat_ctrl`
writes.

**Two residual mismatches remain, and both are now small and specific:**

1. This port emits a leading chunk `ts=7566 te=7571 ti=778 tp=2` that upstream
   does not -- upstream's `peak_content` opens at 7571. Position 7566 is a paired
   position with `treat=2.0`, so the chunk is being selected by the `tpos > 0`
   filter; the q-score at 7566 must be crossing this port's threshold earlier
   than upstream's. Since `qpos` is sampled by point lookup at `p` while the
   paired array is end-indexed, **the same off-by-one that F81 fixed for
   treatment almost certainly applies to the q-score track**, and that is the
   next thing to check.
2. `tp` for `ti=780` is 3 here versus 4.0 upstream, i.e. this port's treatment
   track has no breakpoint at 7596 while upstream's does. Consistent with (1):
   a single missing treatment breakpoint explains both.

**Peak boundaries are 6-9 bp early and summits 3-23 bp off**, so this is progress
(F79 had boundaries 3-6 bp and summits 1-3 bp; summits got worse, boundaries
similar) but **not** a fix, and I am not claiming otherwise. The chunk indices
now line up, which is what makes the remaining two items diagnosable rather than
diffuse: with `ti` matching upstream, the treatment track and the q-score track
can each be compared against the captured `peak_content` position by position.

## F82 -- union composition reconciles exactly; the residual is a genuine treatment-pileup divergence

F81's residual 2 needs a measurement rather than more reading: how do this port's
merged-lambda breakpoints relate to upstream's paired positions, set-wise?

```text
upstream paired positions                     10684  (all distinct)
paired positions NOT breakpoints of my lambda    335
my lambda breakpoints NOT in paired             736
```

**The arithmetic closes exactly.** `736 - 335 = 401`, which is precisely the
lambda-versus-paired truncation already identified at F76 (`11085 - 10684 =
401`). So the two sets reconcile:

* 10,349 paired positions are lambda breakpoints in both;
* 335 paired positions come from the **treatment** track alone -- the
  treatment-only positions F80 predicted, including the `7571` class;
* 736 lambda breakpoints lie beyond the paired extent, which is the truncation.

**This confirms F80's union is structurally right and complete.** It is not a
coincidence of totals this time -- it is set-wise agreement with an independent
third number (401) that was not used to construct it.

**It also settles the `tp=3` versus `tp=4.0` at `ti=780` against reading A.**
`7596` is a *lambda* breakpoint in this port, not a treatment one, so
`end_val(&s.treat, 7596)` finds no treatment run ending there and falls back to
the run ending at `7590`, value `3`. Upstream's `treat_array[780] = 4.0` at
`pos = 7596` means the value over `[7590, 7596)` is `4.0` -- so **upstream's
treatment pileup has a breakpoint at 7596 that this port's does not.** That
confirms the end-indexed reading (F81) is right and that the earlier point-lookup
reading was wrong for a different reason than F81 assumed.

**So the remaining G9 gap is a treatment-pileup parity defect, not a peak-calling
one.** That is a meaningful relocation: G5's 1,776 bit-exact vectors did not cover
this corpus, and the acceptance criteria require treatment-pileup parity across
the whole corpus anyway. The discriminating measurement is direct and needs no new
oracle capture: compare this port's treatment track breakpoints against the
positions that are in `paired` but not in `lambda` (the 335), and for each,
compare `treat_array[i]` from the F73 capture against this port's treatment value
over the same interval.

**Expected shape if this reading is right:** disagreements confined to the 335
treatment-only positions, with this port's treatment *constant* across intervals
where upstream's steps. If instead the disagreements are spread across all
10,684 positions, the treatment track is wrong more broadly and the per-scale
control numbers from F74 -- which matched exactly -- are the anomaly worth
re-examining.

## F83-F86 -- G9 narrow peaks reach exact parity on the normal path

Four fixes, in dependency order. The end result is the first configuration in
this investigation with exact peak coordinates **and** summits produced from the
production path, with no `--peak-content` input.

**F82's conclusion was wrong, and its own measurement said so.** F82 inferred a
missing treatment breakpoint from `tp=3` versus upstream's `4.0`. Comparing the
two treatment signals as step functions over the union of their intervals -- the
F73 method -- gave **0 of 10,965 shared intervals differing**. The treatment
pileup was bit-exact all along; the bug was in the lookup F81 had introduced.
`end_val` took the run with the largest `end <= p`, which for `p = 7596` skips the
run `(7590, X, 4.0)` because `X > 7596`, and returns the previous run's `3.0`.
The value over `[lo, p)` is the run **containing `lo`** (`a <= lo < b`, runs are
half-open) -- F83's `containing`.

**F84 applied the same rule to the q-score track,** which F81 had left as a point
lookup. That removed a spurious leading chunk (`ts=7566 te=7571 tp=2`) that
upstream's `peak_content` does not contain, because the q value to the *right* of
a position was being compared instead of the value to its *left*. After F84 the
first chunks match the capture exactly:

```text
upstream  ts=7571 te=7590 ti=779 tp=3.0  |  ts=7590 te=7596 ti=780 tp=4.0  |  ts=7596 te=7604 ti=781
this port ts=7571 te=7590 ti=779 tp=3    |  ts=7590 te=7596 ti=780 tp=4    |  ts=7596 te=7604 ti=781
```

**With the chunks now identical, what remained was a uniform 1 bp offset on every
peak start and on most summits** -- the signature of a single convention rather
than a chunk-set error.

**F85** is the XLS start convention: upstream reports `tstart + 1`. `tend` needs no
adjustment, being already the inclusive XLS end. **F86** is the summit: it is
`(tend + tstart + 2) // 2`, i.e. the midpoint of the *reported* (already-shifted)
edges, not of the raw chunk edges. F62 had recorded `(tend + tstart + 1) // 2`,
which is right when `tstart` is the raw chunk start but is 1 low once the start
is reported one base higher.

**Result on `se_model/realistic`, production path, no captured inputs:**

```text
CALLPEAK PEAK MATCH (12 peaks, coordinates + summits)
```

All 12 peak starts, ends and summits are byte-identical to pinned MACS3 3.0.5.

**The lesson worth keeping from F82.** I inferred a data defect from one mismatched
scalar, and stated it as a relocated root cause with a predicted signature. The
step-function comparison that would have settled it in one command had already
been written for F73 and was sitting in the findings file. Two measurements --
total counts, then a single value -- both pointed at the data; only the interval
comparison pointed at the code. Cheap, well-scoped differentials keep getting
deferred in favour of plausible one-scalar stories, and the cost here was a full
turn spent asserting something the next turn's first measurement refuted.

## F87-F88 -- G9 narrow peaks: corpus-wide parity, and two harness bugs that were hiding it

Extending F83-F86 from one fixture to the whole corpus required fixing the gate
harness itself. Both bugs were under-reporting success, which is the dangerous
direction: they made a green result look like a red one and would equally have
made a red result look green.

**F87: the comparator read full agreement as failure.** On parity the harness
prints a single `CALLPEAK PEAK MATCH` line and no per-peak lines, so the
per-fixture regex found nothing and reported `0 0 0` -- counted as zero peaks
agreeing rather than all of them. `se_model/realistic` showed `12 / 12 / 0 / 0 /
0`. Now `PEAK MATCH` counts as full parity for the fixture. A differential whose
success path emits a *different shape* than its failure path will silently
misreport; the harness must treat a whole-fixture match as a first-class result.

**F88: the harness replayed flags that were not the fixture's.** It hardcoded
`--gsize 2000000 --extsize 200 --slocal 1000 --llocal 10000` for every fixture,
but each golden was generated with that fixture's own command line. The effective
genome size feeds the p-value -> q-value map, so the six `gonechrom_*` fixtures
generated with `-g 12000` were being scored against a q map for a 2 Mb genome,
which moved the cutoff crossing: each peak came out **1-8 bp wider on the left**
while its summit and right edge stayed exact. That signature -- one edge wrong,
the other exact -- is what identified it, and it is why the earlier runs showed
`summits` matching while `bnd_ok` did not. The harness now reads each fixture's
real flags from `command.json`.

**Corpus result, `poisson` mechanism (upstream's SE default), single-end:**

```text
TOTAL (poisson)   rust 83  upstream 83  bnd_ok 83  summits 83  worst_bnd 0
  fixtures compared         : 152  (skipped upstream-reject: 0)
  fixtures with full parity : 97
```

Every one of the 83 peaks the corpus produces matches pinned MACS3 3.0.5 exactly
in start, end and summit. The remaining 55 fixtures call no peaks on either side.
The `subtract` mechanism is reported alongside precisely so it cannot be
mistaken for the default: it overcalls (121 vs 83) and is not at parity, which
is consistent with F48.

**The signal gate is unaffected**, which matters because F76 changed the lambda
track's construction:

```text
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

So `from_runs_exact` fixed the peak pipeline without disturbing the control
bedGraph, which keeps its coalescing form.

## F89-F90 -- paired-end: harness support, and one argument-order bug

The corpus has **136 paired-end fixtures** (each with 13 golden variants) that
had never been driven: `run_peak_e2e.sh` skipped anything with a `treat.bedpe`
because the harness only knew how to read `treat.bed`. Those fixtures would have
reported a vacuous 0-vs-0, which is the same failure mode as F87 -- a skip that
looks like agreement.

**F89 -- one implementation, two modes.** Rather than add a second peak-caller,
the single-end and paired-end paths now build their per-chromosome signals
differently (position pileup vs fragment pileup) and converge on one shared
`finish` function holding the p-score, q histogram, chunk construction, peak
closing, and XLS comparison. F77 was precisely the cost of maintaining a second
implementation of already-correct logic, and `ChromSignals` now carries an
explicit `chrom` id so the shared pipeline needs no genome handle. The harness
auto-detects `treat.bedpe`, or takes `--format BEDPE`; the gate replays `-f` from
each fixture's `command.json` alongside the other flags.

`d` is the mean template length of the retained fragments, which is what upstream
records as `# d = ...` under `--nomodel` -- verified on
`gtiny_mpe_d1200_w600_noc_310`, whose 509 retained fragments are all length 180
and whose XLS header reads `# d = 180`.

**F90 -- `rlength` is not `d`, and passing `d` silently deleted two thirds of the
genome.** `pileup_from_fragments(chrom, fragments, rlength, scale_factor, base)`
takes the track *length* third and clamps every fragment endpoint to it. Passing
the fragment length `d` (180) truncated all fragments at 180 and left a 76-entry
position list covering `0..180` of a genome whose fragments run to 1077, so
paired-end called **no peaks at all**. The symptom -- a position list that stops
exactly where `d` stops -- is what identified it; the argument is also named
`rlength` rather than `d`, which is a hint worth heeding before the next call.

**Paired-end is now driven and partially green.** Corpus, `poisson`:

```text
TOTAL (poisson)   rust 135  upstream 173  bnd_ok 99  summits 100  worst_bnd 44
  fixtures compared         : 288
  fixtures with full parity : 110
```

Decomposing against the single-end result (83 peaks, all exact): paired-end
contributes 52 of upstream's 90 peaks, with 16 boundary and 17 summit matches.
So **38 peaks are still missed** and boundary disagreements reach 44 bp -- real
remaining work, now visible rather than hidden behind a skip. Single-end parity
is unaffected at 83/83.

**What I would check next**, in order: (1) whether the missing peaks are rejected
by `min_length`, which in paired-end mode should be `d` rather than `--extsize`,
and the harness still passes `a.extsize` into `finish`; (2) whether
`max_gap` should likewise be `d`; (3) the 44 bp boundary outliers, which are
large enough to suggest the treatment track's scale factor rather than a
convention.

## F91-F93 -- paired-end: `d` is the fragment length, truncated; ENV.lock restored

Resuming after the auto-continue limit, the first concrete step from the handoff.

**F91: `min_length` and `max_gap` are `d`, not `--extsize`.** Upstream says so in
its own banner ("The minimum length of peaks is assigned as the predicted
fragment length d") and sets both from the estimated `d`. In single-end `d` *is*
`--extsize`, so passing `a.extsize` was invisible there -- which is exactly why it
survived. Paired-end has `d = 180` against `--extsize = 200`, and regions 180-200
bp long were being dropped that upstream keeps. `finish` now takes `d` as a
parameter.

**F92: upstream truncates `d`, it does not round.** `gmini_mpe_d1200_w60_ctrl_012`
has mean template length 149.999167 and upstream records `# d = 149`; rounding
gives 150. Confirmed against `gtiny_mpe_d1200_w600_noc_310`, whose mean is
exactly 180.0 and where the two agree, so a single fixture could not have
distinguished them -- both were checked.

**F92b was tried and reverted.** I hypothesised that `d` comes from the track
*before* duplicate filtering (upstream prints it before filtering, while the
pileup uses the retained set), changed `load_bedpe` to capture the pre-filter
mean, and it cost 4 peaks. The measurement that should have come first showed it
could not be the explanation: that fixture's peak is 160 bp long, which clears
both 148 and 149, so `min_length` cannot be what rejects it either way. Reverted
rather than left in on a plausible-sounding story -- the same trap as F82, caught
this time by looking before editing.

**F93: `oracle/ENV.lock` recorded `MACS3_VERSION=unknown`**, which was failing
`oracle_lock_records_the_pinned_version`. It is read from
`MACS3/Utilities/Constants.py` (`MACS_VERSION = "3.0.5"`) at pinned commit
`c5443190e3edfeb301cc94acf450e2b2c026a223`, and the lock now records that commit
so the pin is checkable instead of asserted. The file also documents the standing
hermeticity caveat: `CYTHON=3.3.0` is present because the oracle is instrumented
in place, so the tree is **not** byte-identical to `c544319` and must be restored
from `/tmp/CallPeakUnit.py.orig` and `/tmp/oracle_so_backup/` before any release
claim.

**Where the corpus stands after F91-F92:**

```text
TOTAL (poisson)   rust 141  upstream 155  bnd_ok 99  summits 100  worst_bnd 48
  fixtures compared         : 288
  fixtures with full parity : 109
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

Upstream's 155 peaks split 83 single-end (all exact) and 72 paired-end; this port
now calls 58 of the 72. **14 paired-end peaks are still missed.**

**The next lead is the control scale factor, and the evidence is quantitative.**
For `gmini_mpe_d1200_w60_ctrl_012`, upstream's `fold_enrichment` of 1.23882 with
`pileup` 87 implies its local lambda at the summit is **70.228**. This port's raw
control depth there is **153**, so upstream applied a scale of
`70.228 / 153 = 0.4590`. None of the candidate rules reproduces it: `n_C/n_T`
gives 0.4353, the mean-length-weighted sum ratio gives 0.6511, and the
extsize-weighted ratio gives 0.3238. The window geometry matters here -- this
fixture's control fragments average 99.45 bp against a `d` of 149, so a
window-averaged lambda is not the raw depth at all, and the arithmetic above is
the wrong model rather than merely the wrong constant.

The discriminating measurement, needing no new capture: dump this port's
per-scale control tracks for a PE fixture (the same instrumentation added to
`pileup_a_chromosome_c` at F74) and compare each scale's value at the summit
against upstream's implied 70.228. That isolates whether the error is in the
per-scale windowing, in the fold, or only in the `1/ratio` rescaling applied
afterwards.

## F94-F95 -- paired-end scaling read from source; G9 narrow peaks green in both modes

The handoff's lead was the paired-end control scale factor. Rather than keep
reverse-engineering it from `fold_enrichment`, this reads `PeakDetect.py:155-181`
and `callpeak_cmd.py` directly. Three things were wrong, and all three were
invisible in single-end -- which is why they survived the entire G9 campaign.

**F94: the paired-end control is counted at *both* fragment ends.** Upstream's own
comment: "entire fragment is counted as 1 in treatment whereas both ends of
fragment are counted in control/input", and `callpeak_cmd.py:216` doubles `c1` for
the same reason. So `control_sum = control.total * 2 * avg_template_length`, while
the treatment keeps one count per whole fragment via `treat.length` (the sum of
fragment *lengths*, not `total * d`). The control must therefore be piled up from
single-end tags at both endpoints, not as fragments.

Passing reversed `(r, l)` pairs to `pileup_from_fragments` does **not** work: it
assumes `start <= end`, so the sweep misplaces them and the control collapses to
nothing (PE dropped to 0 peaks). The control is piled up with
`pileup_from_positions(bidirectional(d))` over the fragment starts and ends, which
is the same shape as the single-end control.

**F92b, which I reverted last turn, was right after all -- and the reason I gave
for reverting it was sound but drew the wrong conclusion.** I had verified the
160 bp peak cleared both candidate thresholds, so `min_length` could not explain
that fixture's rejection. True -- but that only showed `min_length` was not the
cause *there*. `callpeak_cmd.py:165` sets `options.d = options.tsize` under
`--nomodel` in PE mode, and `options.tsize = tp.d` (`callpeak_cmd.py:362`) is the
mean template length of the track **as read**, before duplicate filtering. Three
fixtures confirm it independently:

| fixture | all 1200 frags, truncated | retained, truncated | upstream `# d =` |
|---|---|---|---|
| `gmini_mpe_d1200_w180_ctrl_039` | **145** | 122 | **145** |
| `gmini_mpe_d1200_w60_ctrl_012` | **149** | 148 | **149** |
| `gtiny_mpe_d1200_w600_noc_310` | 180 | 180 | 180 |

**F95: upstream uses two *different* fragment lengths in PE mode.** `ctrl_d_s` is
seeded with `self.d` (= `opt.d` = `tsize`, the pre-filter truncated mean), which is
the first scale's extension length and `min_length`/`max_gap`. But the
`slocal`/`llocal` scale factors use a *local* `d` bound at `PeakDetect.py:156` to
`self.treat.average_template_length` -- the **post-filter** mean, as a float. Here
that is 145 versus 122, a factor of 1.19 on every wide-window lambda. Collapsing
the two is what left four `w180` fixtures overcalling and five `w600` fixtures 9-12
bp out.

**Result. `TOTAL (poisson)` over 288 fixtures:**

```text
rust 155   upstream 155   bnd_ok 152   summits 153   worst_bnd 1
  fixtures with full parity : 174
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

Peak **counts now match exactly** across the whole corpus, and single-end remains
83/83 exact. Paired-end went from 52 of 72 peaks (F89) to 70 of 72, and worst
boundary disagreement from 48 bp to **1 bp**. Five fixtures remain, all off by
exactly one base at a region edge -- three where this port's first above-cutoff
position is one later than upstream's (start 74 vs 73) and two where the summit
chunk differs by one position (153 vs 152). Those are single-position differences
in where the q cutoff is crossed, not conventions.

**The lesson from the revert.** I reverted F92b on evidence that was correct but
narrow: it proved `min_length` was not the cause of *one* fixture's rejection, and
I generalised that to "pre-filter `d` is wrong". The right conclusion was "this
fixture cannot discriminate, find one that can" -- and `gmini_mpe_d1200_w180_ctrl_039`,
with a 23 bp gap between its pre- and post-filter means, was sitting in the corpus
the whole time. Locating a discriminating fixture is cheap; discarding a hypothesis
because the nearest test was insensitive is not.

## F96 -- the fragment pileup was coalescing; correct, but not the 1 bp cause

`pileup_from_fragments` passed `coalesce = true` to `sweep_to_track`. Upstream's
fragment sweep is `_pileup_sorted_unit_as_list` (`PileupV2.py:344`), which appends
one entry per sweep step into raw arrays and never merges neighbours -- the same
property that forced `from_runs_exact` in F76. Set to `false`.

**It changed nothing measurable** (`bnd_ok` and `summits` identical before and
after), because these fixtures have no equal-valued adjacent fragment runs. The
change is kept because it is what upstream does and because leaving a
coalescing default on a path whose contract is explicitly non-coalescing is the
kind of latent bug that only fires on other data. Recorded as a null result rather
than as a fix.

**The 1 bp residuals are a lambda *breakpoint* difference at region edges, and
the dump localises it.** For `gmini_mpe_d1200_w600_ctrl_066` the port's first
above-cutoff position is `i=72, pos=74`, so its first chunk starts at 73. Upstream
reports peak start 73, i.e. chunk `tstart = 72`. The two position arrays therefore
differ:

```text
this port:  ... 72:73  73:74  74:75 ...
upstream :  ... 71:72  72:74  74:75 ...   (implied: has 72, lacks 73)
```

Positions `72` and `73` are both well inside the peak, and the port's control is
flat at `54.76142` across the whole neighbourhood while treatment climbs
`70, 71, 72, ...`. So the port's merged lambda is missing a breakpoint at 72 that
upstream's has -- the *same* class as F80/F82, where the position list is the union
of treatment and control breakpoints and the control side contributes one the
treatment side does not.

**The discriminating measurement** is the per-scale control capture added at F74,
run on a paired-end fixture: compare each scale's breakpoints in the 65-80 window
against the port's `ends` union. That separates "the per-scale fragment pileup has
extra breakpoints upstream lacks" from "the fold drops them" -- and, unlike
inferring it from a 1 bp peak shift, it does not depend on the peak caller at all.

**Five fixtures remain**, all off by exactly one base: three where the port's first
above-cutoff position is one later than upstream's, two where the summit chunk
differs by one position (153 vs 152) with boundaries exact. Corpus stands at:

```text
rust 155   upstream 155   bnd_ok 152   summits 153   worst_bnd 1
  fixtures with full parity : 174
```

## F97 -- the paired-end lambda is 0.16 % high; the scales now mirror upstream exactly

F95's implementation assembled the extents and factors with an extra rescale
whose algebra was wrong (`d_extent/avg_tl` where upstream means `avg_tl/scale`).
The scales are now built the way `PeakDetect.py:192-220` builds them: extents
`[opt.d, sregion, lregion]`, first factor `ratio`, wide-window factors
`d/sregion * ratio` on the **post-filter** `d`. The two are deliberately different
bases and are now separate expressions rather than one number reused.

**The fix changed nothing measurable** (`bnd_ok` 152, `summits` 153, `worst_bnd` 1
before and after). That is informative: the wide-window scales are not the
argmax anywhere in this corpus, so their factor is irrelevant to the result and my
algebra error was latent. Recorded as a null result; the code is kept because it
now matches the source rather than coincidentally agreeing.

**The residual is a 0.16 % lambda offset, and the numbers pin it precisely.** For
`gmini_mpe_d1200_w600_ctrl_066` upstream's row reports `pileup 169`,
`fold_enrichment 1.35317`, implying a local lambda at the summit of

```text
169 / 1.35317 = 124.897      (upstream)
this port's ctrl at pos 189   = 125.101524
difference                   =   0.2045   = 0.16 %
```

The treatment matches exactly (both 169), and the summit coordinate matches, so
the offset is entirely in the control. It is **not** a depth difference: 0.2045
divided by this port's `ratio` (0.47208122) is 0.433, which is not an integer, so
the control depth at the summit is right and the **scale factor** is 0.16 % too
large. The needed ratio is 0.471306 against 0.47208122 computed.

Since `ratio = treat.length / (control.total * 2 * avg_template_length)` and the
treatment pileup is already exact, the candidates are narrow: a 0.16 % difference
in `treat.length` (the sum of fragment lengths), or in the post-filter mean used
as `avg_template_length`, or `lambda_bg = treat.length / gsize` landing 0.2 low
against upstream's. Each is a single number that can be printed and compared
against upstream's own, and the `fold_enrichment` column gives an independent
check at every summit, so the corpus can be swept for the pattern rather than one
fixture examined at a time.

A 0.16 % lambda shift is enough to move the q cutoff across a position: at
`pos 73` this port's q sits just under `1.30103`, and a slightly lower lambda puts
it over. That accounts for the three leading-edge fixtures exactly.

**Corpus is unchanged and near-closed:**

```text
rust 155   upstream 155   bnd_ok 152   summits 153   worst_bnd 1
  fixtures with full parity : 174
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

**Next measurement, and it is a sweep rather than a single fixture:** print
`treat.length`, `control.total`, the post-filter mean and the implied ratio for
every paired-end fixture, and compare each implied lambda against this port's at
the summit. Whether the 0.16 % is constant or varies with `control.total` decides
whether it is a rounding artefact (`0.16 %` is roughly `f32` epsilon accumulated
over a windowed average) or a genuine arithmetic difference.

## F98 -- sweeping every paired-end summit shows the 0.16 % is not a constant

F97 predicted the offset might be constant (f32 epsilon) or arithmetic. It is
neither: it affects a **subset**. The harness now emits `LAM\tchrom\tsummit\tpileup\tfold`
and `oracle/sweep_pe_lambda.py` compares every paired-end summit against the
`fold_enrichment` column of the golden XLS, which is an independent read of the
same lambda. Over **70 compared summits**:

```text
sweep/gmini_mpe_d400_w600_ctrl_057    up=1.27193  mine=1.26675  -0.4074%
sweep/gmini_mpe_d400_w600_ctrl_219    up=1.30591  mine=1.30081  -0.3903%
sweep/gmini_mpe_d400_w180_ctrl_030    up=1.31890  mine=1.31376  -0.3896%
sweep/gmini_mpe_d1200_w600_ctrl_066    up=1.35317  mine=1.34812  -0.3732%
sweep/gmini_mpe_d4000_w600_ctrl_075   up=1.34569  mine=1.34074  -0.3678%
sweep/gmini_mpe_d4000_w600_ctrl_399   up=1.34569  mine=1.34074  -0.3678%
sweep/gmini_mpe_d4000_w600_ctrl_390   up=1.30886  mine=1.30419  -0.3567%
sweep/gmini_mpe_d1200_w600_ctrl_228    up=1.30707  mine=1.30244  -0.3541%
...
sweep/gmini_mpe_d4000_w600_ctrl_237   up=1.33576  mine=1.33579  +0.0025%
sweep/gmini_mpe_d400_w180_ctrl_354    up=1.33521  mine=1.33523  +0.0018%
sweep/gmini_mpe_d4000_w60_ctrl_183    up=1.23144  mine=1.23146  +0.0018%
sweep/gmini_mpe_d4000_w60_ctrl_327    up=1.25093  mine=1.25094  +0.0011%
```

**The majority are exact to within 0.003 %**, which is f32 noise. Twelve are low
by 0.35-0.41 % -- three times the noise floor and a factor of two narrower than
the 0.16 % F97 extrapolated from a single fixture, so that estimate was itself an
artefact of one sample.

**A constant formula error is ruled out.** If `ratio`, `lambda_bg` or the mean were
wrong by a fixed proportion, every summit would be off by that proportion. The
split is per-fixture, and `gmini_mpe_d4000_w600_ctrl_237` is *exact* while its
siblings `075` and `399` from the same generator are off by 0.37 % -- same
parameters, same code path, different answer.

**The shape of the error points at max-selection, not arithmetic.** The affected
group is uniformly about -0.38 %, which is what you get when the argmax over the
three lambda scales picks a *different scale* than upstream's does, or picks the
same scale from a track that differs in one narrow window. An arithmetic error
would move all summits; a scale-selection error moves only those summits whose
lambda is decided near a tie between two scales, which is a minority by
construction. That also explains why it is concentrated in the high-coverage
fixtures: their lambdas sit closer together across scales, so ties are more
likely.

**The discriminating measurement is now cheap and per-fixture:** for one affected
fixture (`gmini_mpe_d400_w600_ctrl_057`), print all three per-scale control values
at the summit alongside the merged max, and compare which scale wins against
upstream's own `fold_enrichment`. If upstream's value equals this port's *second*
largest scale, it is max-selection; if it is between scales, it is a per-scale
track difference. Both are one dump away, and the sweep script makes the affected
set enumerable rather than anecdotal.

## F99 -- the control is sampled at `pos[i]`, the treatment at `pos[i-1]`

F98's per-scale probe paid off immediately, and the answer was not the one the
max-selection hypothesis predicted.

**The merge was correct all along.** Probing an affected fixture at four adjacent
positions shows the merged lambda tracking the `d`-scale exactly:

```text
at=177  d=119: 90.947372  slocal: 13.268001  llocal: 3.790857  merged: 90.947372
at=178  d=119: 90.573097  slocal: 13.268001  llocal: 3.790857  merged: 90.573097
at=179  d=119: 90.198830  slocal: 13.268001  llocal: 3.790857  merged: 90.198830
at=180  d=119: 89.824562                                merged: 89.824562
```

No max-selection error: the merged value equals the largest scale at every
position. F98's hypothesis is refuted, and it is recorded as refuted rather than
quietly dropped.

**The real defect was a one-position asymmetry in sampling.** Upstream's
`fold_enrichment` at the summit implies a lambda of **90.2000**, which is the
control value **at 179**. This port was reporting **90.5731**, the value at
**178** -- the summit chunk's *start*. So the control was being sampled over
`[pos[i-1], pos[i])` when upstream samples at `pos[i]`.

The fix is one line, and it is deliberately **not** applied to the treatment:

```rust
cpos.push(ctl.map_or(0.0, |c| containing(c, p)));   // was containing(c, lo)
```

**The treatment keeps `containing(treat, lo)`** because it is already right --
`pileup 115` matched upstream exactly on every fixture checked, and F83
established `treat_array[ti]` as the value over `[pos[ti-1], pos[ti])`. The two
tracks genuinely differ, and that is not a contradiction: `d_pileup_d` and
`ctrl_d_pileup_d` are separate arrays with separate breakpoint sets, and the
control's array is indexed by the *paired* index while being evaluated in the
control's own interval space. Treating them as symmetric is exactly the
assumption F99 removes.

**Verified on the fixture that motivated it:**

```text
before:  fold 1.266748   (upstream 1.27193, -0.41 %)
after:   fold 1.271946   (upstream 1.27193, +0.0013 %)
```

**And the F98 group is gone.** Re-running the 70-summit sweep, the entire
0.35-0.41 % cluster has disappeared; the worst remaining disagreement is a
*different* fixture at +2.4 %, with a spread of roughly -1.0 % to +2.4 % and no
sign bias. So F98's "subset, not constant" observation was right about the
*shape* and wrong about the *mechanism*: the subset was exactly the set of
fixtures whose summit lambda changed materially under a one-position shift.

Note what did **not** move: corpus `bnd_ok` 152, `summits` 153, `worst_bnd` 1,
174 fixtures at full parity. The lambda is materially more correct and the five
1 bp residuals are unchanged, so those residuals are caused by something the
lambda fix does not touch. The sweep is now the instrument that would find it:
it compares 70 independent lambda readings rather than 5 peak coordinates, and
it is far more sensitive to a one-position error than the peak gate is.

**Verified:** 333 passed, 0 failed, 0 clippy, fmt clean, signal gate 97/97.

## F100 -- the paired-end control depth is off by exactly +-1 count at the summit

F99 removed the percentage-scale error, and the sweep that remained is small enough
to be measured in **units of the control count** rather than as a percentage. For
each of the 70 paired-end summits, divide `(lambda_this_port - lambda_upstream)` by
this port's `ratio` -- which converts a lambda difference into a **number of
control tags**:

```text
sweep/gmini_mpe_d4000_w600_ctrl_237   depth_units = -1.0068
sweep/gmini_mpe_d4000_w60_ctrl_183    depth_units = -1.0058
sweep/gmini_mpe_d400_w180_ctrl_354    depth_units = -1.0046
sweep/gmini_mpe_d4000_w60_ctrl_345    depth_units = -1.0037
sweep/gmini_mpe_d1200_w60_ctrl_336    depth_units = -1.0032
sweep/gtiny_mpe_d1200_w600_ctrl_067   depth_units = +1.0002
sweep/gtiny_mpe_d400_w60_ctrl_166     depth_units = +1.0000
sweep/gtiny_mpe_d400_w180_ctrl_193    depth_units = -1.0000
sweep/gtiny_mpe_d1200_w600_ctrl_229   depth_units = -1.0000
```

**Every one of the 70 summits is within a thousandth of a whole tag.** The
remaining paired-end error is not a scale factor, not a rounding artefact, and not
a windowing difference: the control pileup at the summit carries exactly one tag
more (mostly) or one tag fewer (occasionally) than upstream's.

That also reinterprets the whole of F97-F99. The "0.16 %" was one tag on a small
lambda; the "0.35-0.41 %" cluster was one tag on a different-sized lambda; and the
sign flips between fixtures because the error is genuinely `+-1` rather than a
bias. F99's one-position sampling fix was real and necessary -- it moved the
control from the wrong interval to the right one -- but it addressed the *interval*,
not the *count*.

**The cause is an endpoint convention in the paired-end control pileup.** Upstream
builds it with `pileup_from_PN_shifted(P_array, N_array, five_shift, three_shift,
...)` from the fragment starts and ends, and this port uses
`SingleEndParams::bidirectional(d, 0, rlength, factor)`. Both place a plus event
and a minus event per fragment end, so they agree in shape; they disagree by one
count on positions where an event lands exactly on a window boundary, which is why
the sign is not constant -- `d/2` rounds differently for odd `d` in the two
formulations.

**The discriminating measurement is now a single fragment.** Pick one summit with
`depth_units = +1` and one with `-1`, find a control fragment whose start or end
lands within one base of `d/2` modulo the summit position, and compare the
`pileup_from_PN_shifted` event positions against `bidirectional`'s for that one
fragment. Both are pure functions of `(position, d, shifts)` with no global state,
so they can be compared directly rather than inferred from peaks.

**Why this matters beyond paired-end.** A one-tag discrepancy in the local lambda
is exactly the scale of error that moves a q cutoff across a position, which is
what the five remaining 1 bp boundary residuals look like. It is also the kind of
defect that a percentage-based gate hides: at a lambda of 125 a single tag is
0.09 %, comfortably inside any tolerance that looks reasonable, but it is a
systematic endpoint error rather than noise. The unit-based sweep is a strictly
better instrument and is kept in `oracle/sweep_pe_lambda.py`.

**Corpus unchanged** -- `rust 155 / upstream 155, bnd_ok 152, summits 153,
worst_bnd 1, 174 fixtures at full parity` -- because a one-tag lambda error moves
summits only when it changes which chunks tie at the maximum, which is rare. The
lambda is nevertheless now measurably wrong in a *characterised* way rather than an
unexplained one, which is the difference between a bug and a mystery.

## F101 -- the paired-end control is symmetric, and G9 narrow peaks close corpus-wide

F100 said the cause was an endpoint convention in the paired-end control pileup
and pointed at `pileup_from_PN_shifted`. Reading it confirmed the geometry and
found the actual difference, which is *not* in `pileup_from_PN_shifted` itself but
in what the two callers pass it.

`FixWidthTrack.pileup_a_chromosome_c` (single-end control) passes:

```python
five_shift  = d//2 - end_shift
three_shift = end_shift + d - d//2
```

which for `end_shift = 0` is the **asymmetric** split `d//2` / `d - d//2` -- exactly
what [`SingleEndParams::bidirectional`] produces, which is why single-end has been
exact throughout.

`PairedEndTrack.pileup_a_chromosome_c` (paired-end control) passes instead:

```python
five_shift  = d//2
three_shift = d//2
```

**Both floored**, giving a span of `d` bases for even `d` and `d - 1` for odd `d`.
The port was using the single-end split for the paired-end control, so every odd
`d` produced a control track one base too long on the minus strand -- F100's
one-tag discrepancy.

The fix is a new constructor, `SingleEndParams::symmetric(d, rlength, factor)`,
used only by the paired-end control. `bidirectional` is untouched, because
single-end genuinely does want the asymmetric split.

**G9 narrow peaks are now closed across the entire corpus, both modes:**

```text
TOTAL (poisson)   rust 155   upstream 155   bnd_ok 155   summits 155   worst_bnd 0
  fixtures compared         : 288
  fixtures with full parity : 179
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

Every peak the corpus produces matches pinned MACS3 3.0.5 exactly in start, end
**and** summit. The five 1 bp boundary residuals and the two summit residuals that
survived F99 are gone; the gate went from `bnd_ok 152 / summits 153 / worst_bnd 1`
to `155 / 155 / 0`.

**A residual remains, and it is honest to say so:** the unit sweep still reports
about one tag of disagreement at the summit on a subset of paired-end fixtures, so
the `fold_enrichment` column is not yet exact even though every coordinate is. The
peak gate does not compare that column, so this is invisible to the gate and would
be a silent gap in the XLS byte-comparison acceptance criterion.

It is the same class of defect and the same next step as before: the surviving
differences are `-1.00x` (a small minority `+1.00x`), which after F101 can only be
the `five_shift`/`three_shift` *clamping* at contig edges or the `d - 1` span for
odd `d` interacting with `fix_coordinates`, since the interior geometry is now
provably identical. The sweep is retained as the instrument; it is strictly more
sensitive than the peak gate, which is exactly why it found this after the peak
gate had already gone green.

**Verified:** 333 passed, 0 failed, 0 clippy, fmt clean.

## F102 -- the merged paired-end control lambda is now breakpoint-identical to upstream

F101 left a residual that the peak gate could not see, so it needed a real
comparison rather than inference from `fold_enrichment`. `PairedEndTrack.pileup_a_chromosome_c`
is now instrumented to dump the merged control's values under `MACS3RS_PECTRL`, and
the harness dumps its own under `CALLPEAK_PECTRL`, so the two can be compared
position by position.

One methodological note worth recording: the first attempt produced nonsense
(upstream appearing to have *no* control coverage below position 5000). The cause
was a capture harness bug -- the upstream invocation omitted `-c`, so it ran
without a control track and the dump was the auto-generated global lambda. The
tell was that the fixture's control fragments span `0..200`, which no coordinate
shift could explain. A capture that disagrees with a fact already known from the
fixture is a broken capture, not a finding.

**With `-c` supplied, the comparison is:**

```text
upstream rows 191   this port's rows 191
breakpoint sets equal: True
worst value difference: at 144, upstream 156.503937, mine 156.500000, diff 0.0039
```

**The breakpoint sets are identical, element for element.** That is the structural
result F80-F102 were working toward: the paired-end control lambda now has exactly
upstream's breakpoints, not merely the same count and a similar function. It is
the same property single-end reached at F76/F82, now confirmed by direct capture
rather than by reconciling counts.

**The residual is 2.5e-5 in relative terms** -- `156.503937 / 156.5 = 1.0000252`,
and the same factor holds at every position (`29.000729 / 29`, `30.000755 / 30`).
A constant *relative* factor, not a constant additive one, which rules out the
baseline and `lambda_bg` (`5.682857` here, nowhere near the `0.000729` gap).

Since `ratio = treat.length / (2 * control.total * avg)` and
`avg = treat.length / treat.total`, the ratio reduces *identically* to
`n_T / (2 n_C) = 200 / 400 = 0.5` regardless of the lengths involved, so the
factor cannot come from the ratio arithmetic at all. The remaining candidates are
therefore narrow: either the depth in one of the two wide-window scales is off by
a fraction of a tag, or a wide-window scale is marginally above the `d`-scale here
and upstream's max picks it while this port's picks the `d`-scale. Both point at
the same place -- the merged max across scales -- and both are worth one dump of
the per-scale values at a single position, which the `CALLPEAK_ATSCALE` probe
already provides.

**G9 narrow peaks remain closed:**

```text
TOTAL (poisson)   rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
  fixtures with full parity : 179
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

## F103 -- end-indexed sampling needs the run *ending* at `p`, not the run containing it

F102's probe produced an apparent contradiction: the merged dump said the value at
position 144 was `156.5`, while `value_at(144)` and the `CALLPEAK_ATSCALE` probe
both said `156.0`. Both were correct, and the discrepancy is the F5 end-indexed
convention meeting the wrong helper.

`containing(t, p)` returns the run with `start <= p < end` -- in end-indexed
coordinates that is the *next* interval, because a run ending at `p` covers
`(prev, p]` and its successor covers `(p, next]`. Sampling "the value at `p`" in
upstream's coordinates therefore wants the run whose end is the largest `<= p`,
which is a different function.

F99 switched the control from `containing(c, lo)` to `containing(c, p)`. That was
right in effect -- it moved the sample from the previous interval to `p` -- but it
used a helper that is only accidentally equivalent when the track happens to have a
breakpoint at `p`, and it leaves the port unable to express the rule. It is now
written directly:

```rust
// end-indexed sample at `p` -- the run whose end is the largest <= p
let at_end = |t: &SignalTrack<f32>, p: Coord| -> f32 { ... };
cpos.push(ctl.map_or(0.0, |c| at_end(c, p)));
```

**Effect on the sweep.** F103 removes the entire 0.35-0.41 % cluster that F101
left behind: `gmini_mpe_d4000_w600_ctrl_075` and `_399` move from `-0.3678 %` to
`+0.0026 %`, i.e. into f32 noise. Of the 14 worst summits the sweep reports, **12
are now within 0.01 %** and 2 remain outside (`-1.16 %` and `+0.37 %`, both
`gonechrom_mpe_*`). The lambda column is materially closer, and what remains is two
fixtures rather than a dozen.

**Two of them are `gonechrom` fixtures, which is a hint rather than a conclusion.**
"gonechrom" fixtures have a chromosome removed relative to their generator
siblings, so their contig set differs; the summit coordinate still matches and the
peak gate is unaffected, which points at the lambda rather than at peak selection.

**Gates unchanged and still green:**

```text
TOTAL (poisson)   rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
  fixtures with full parity : 179
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

**Method note.** Two of the three remaining helper-semantics bugs in this
investigation (F83's `end_val`, F99's `containing`, F103's `at_end`) have the same
root shape: two functions that differ by one interval on a half-open/end-indexed
boundary, both reading as "obviously the position lookup". The
`value_at`/dump/probe disagreement is what exposed it, and that is an argument for
printing the same quantity through two paths whenever a numeric result is hard to
explain -- a single probe would have let F99 stand indefinitely.

## F104 -- the sweep could not tell "almost all good" from "a dozen bad, rest unknown"

F103 reported "12 of the 14 worst summits within 0.01 %", and flagged that the
sweep only prints the worst 14 of 72 -- so the other 58 were never measured. That
is an instrument that cannot distinguish two very different worlds, and it is the
same failure mode as F76 (a count that agreed while the breakpoints did not) and
F87 (a comparator that read success as failure). The sweep now reports the whole
distribution:

```text
compared summits: 72
  <= 0.005%      70  (97.2 % of total)
  <= 0.01%       70  (97.2 %)
  <= 0.05%       70  (97.2 %)
  <= 0.25%       70  (97.2 %)
  <= 1%          71  (98.6 %)
  > 1%           72  (100 %)
median |dev| 0.000093 %    p90 0.001593 %
outside 0.01 %: 2
```

**The measurement is much better than F103's tail suggested and much better than
"12 bad fixtures" would have implied.** 70 of 72 paired-end summits agree to within
0.005 %, a median deviation of 9.3e-5 % and a p90 of 1.6e-3 % -- f32 noise
throughout. F103's "12 of 14 within 0.01 %" was reading the tail of a distribution
whose body it had not measured; the tail happened to contain both the good and the
bad.

**The two real outliers are characterised, not just located:**

```text
sweep/gonechrom_mpe_d400_w600_ctrl_383   -1.1578 %
sweep/gonechrom_mpe_d1200_w180_ctrl_041   +0.3731 %
```

Both are `gonechrom_*`, and both carry **exactly one chromosome** (`chr1` only)
with `n_T == n_C` (400/400 and 1200/1200). That is the distinguishing property:
with treatment and control the same size, `to_control` is false, `ratio` reduces to
`n_T / (2 n_C) = 0.5`, and the whole scaling rests on the control-doubling of F94
being exactly right. A `+/- 1 %` error in that one factor shows up only when the
two arms are balanced, which is why no fixture with unequal arms ever exposed it.

That is a testable prediction rather than a story: a synthetic paired-end fixture
with `n_T == n_C` and several `d` values should reproduce the deviation, and
correcting the doubling in that regime should move both fixtures to f32 noise. It is
also the last known gap to byte-exact `fold_enrichment`.

**G9 narrow peaks remain closed, unchanged by F104:**

```text
TOTAL (poisson)   rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
  fixtures with full parity : 179
```

**Verified:** 333 passed, 0 failed, 0 clippy, fmt clean.

## F105 -- F104's "balanced arms" prediction is REFUTED; it used pre-dedup counts

F104 predicted that the two remaining `fold_enrichment` outliers were isolated to
the `n_T == n_C` regime, on the basis that both fixtures have a single chromosome
and equal raw fragment counts. **The raw counts are equal; the retained counts are
not**, and it is the retained counts that enter the arithmetic:

| fixture | raw n_T / n_C | retained n_T / n_C | d | lambda_bg |
|---|---|---|---|---|
| `gonechrom_mpe_d400_w600_ctrl_383` | 400 / 400 | **286 / 396** | 180 | 4.290 |
| `gonechrom_mpe_d1200_w180_ctrl_041` | 1200 / 1200 | **209 / 1128** | 180 | 3.135 |

So neither fixture is in the balanced regime; the second is in fact strongly
*control-heavy* (209 vs 1128). F104's prediction is withdrawn.

**This is the third time in this investigation that a pre-filter count has been
mistaken for a post-filter one** (F92b's revert, F95's split of `d`, and now
F104). The reason is mechanical: `wc -l` on the fixture is the cheap, obvious
measurement, and the retained count requires running the deduplicator. Every one
of those three errors came from reaching for the cheap one under budget pressure,
and each cost a turn.

**What the arithmetic does say.** Both fixtures have a very small lambda in which
`lambda_bg` dominates:

```text
383: upstream lambda 4.549989, lambda_bg 4.290  ->  depth contributes 0.259989 (5.7 %)
041: upstream lambda 3.485185, lambda_bg 3.135  ->  depth contributes 0.350185 (10.0 %)
```

With the depth contributing under 10 % of the total, a *one-tag* difference in the
control depth -- F100's residual, still not closed -- moves the lambda by a large
**relative** amount while moving it by a small absolute one. That is consistent
with the observation and does not require the balanced-arm story: these two
fixtures are simply the ones where `fold_enrichment` is most sensitive to the
unresolved tag.

So the open item is unchanged from F100 and is **not** narrowed by F104: a single
control tag at the summit, still unaccounted for, now visible only where lambda is
close to its baseline. The discriminating test is unchanged too -- a synthetic
paired-end fixture with a deliberately small `lambda_bg` (large `--gsize`) and
balanced *retained* counts, which would amplify any one-tag error into plain view.
The point to carry forward is that the amplifying condition is **low signal-to-
baseline**, not arm balance.

**Gates unchanged:** `rust 155 / upstream 155, bnd_ok 155, summits 155,
worst_bnd 0, 179 fixtures at full parity`; signal gate 97/97.

## F106 -- the merged paired-end control is now BIT-IDENTICAL; the fold's lambda is not it

F104 said the two `fold_enrichment` outliers amplify a one-tag error because their
`lambda_bg` dominates. That implies the error is in the control, so the decisive
test is a direct capture on one of those amplifying fixtures --
`gonechrom_mpe_d400_w600_ctrl_383`, where a 1 % relative error is only 0.26 in
absolute terms.

```text
upstream merged control rows : 2820
this port's merged rows      : 2820
breakpoints equal            : True
positions differing          : 0 of 2820
```

**The merged paired-end control lambda is bit-identical to upstream's, at every
one of 2820 breakpoints.** F101 (symmetric extension) and F103 (end-indexed
sampling) together closed the control track completely. F100's one-tag discrepancy
is gone at the track level, not merely reduced.

**And yet this fixture's fold still differs by 1.16 %:**

```text
LAM  chr1 3081 112 20.124667          (this port)
upstream                                 20.3604
implied lambda  = 113/fold - 1 :  mine 4.61187   upstream 4.55001
```

Both lambdas exceed `lambda_bg = 51480/12000 = 4.29`, so both have a depth
contribution (`0.32187` and `0.26001`), and neither is `ratio = 286/(2*396) =
0.361111` times a whole number of tags. Since the merged control is bit-identical,
**the lambda entering `fold_enrichment` is not the merged control value at the
summit.** That is a genuine, specific contradiction rather than a tolerance
question, and it relocates the last lambda defect out of the control track and
into whatever builds `paired.ctrl` -- the per-position resampling of the merged
track that the peak caller reads, or the `lambda_bg` added into each scale before
the max rather than after it.

It also means F104's and F105's framing ("a one-tag control error amplified by low
signal-to-baseline") was the wrong explanation for a residual that no longer
exists at the track level. What survives is a resampling or baseline-placement
difference, and it is confined to fixtures where the depth term is a small share
of the total -- which is why 70 of 72 summits are unaffected and these two are not.

**The discriminating measurement** is a single line: print `paired.ctrl[i]` for the
summit index `i` alongside `merged_ctrl.value_at(pos[i])` and
`merged_ctrl.value_at(pos[i-1])`. If `paired.ctrl[i]` equals neither, the
resampling is dropping or adding the baseline; if it equals the `pos[i-1]` value,
F103's `at_end` did not reach the code path the peak caller uses.

**Gates unchanged:**

```text
TOTAL (poisson)   rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
  fixtures with full parity : 179
lambda sweep: 70 of 72 summits within 0.005 %; 2 outside
```

## F107 -- F103 is correct on the corpus; the two outliers fold a *different* array

F106 said the lambda entering `fold_enrichment` is not the merged control value at
the summit, and named a one-line test. For `gonechrom_mpe_d400_w600_ctrl_383` the
merged control (bit-identical to upstream at all 2820 breakpoints) reads:

```text
upstream implied lambda at summit 3081 = 4.549989
  at_end(3080)  = 4.615000    containing(3080) = 4.550000
  at_end(3081)  = 4.615000    containing(3081) = 4.550000
```

`containing(3081)` matches upstream to six decimals; `at_end(3081)` -- what F103
installed -- does not. On this evidence F103 looks like a regression.

**It is not, and testing it was worth the turn anyway.** Reverting F103 to
`containing` and re-running the full sweep:

```text
at_end (F103, current)      70 of 72 summits within 0.005 %
containing (F103 reverted)  37 of 72 summits within 0.005 %
```

`at_end` is right for 33 fixtures and wrong for at most 2; `containing` is wrong
for 35. **F103 stands.** The single-fixture comparison that appeared to refute it
was a fixture-specific coincidence -- the same trap F76 and F104 describe: one
fixture is not a distribution.

**What the coincidence means.** `containing` happens to reproduce this fixture's
*fold* while `at_end` reproduces the merged *track's* value. Those are different
quantities, and upstream's `fold_enrichment` reads the **paired** control array
`ctrl_d_pileup_d`, which `__chrom_pair_treat_ctrl` builds by pairing treatment
against control -- not the merged local lambda. That is the array-identity
distinction F35 and F80 established, resurfacing at the last stage: where the two
arrays agree both helpers agree; where they differ, only the helper that samples
the paired array reproduces the fold.

The port has no `ctrl_d_pileup_d` counterpart -- `cpos[i]` samples the merged track
-- correct for 70 of 72 summits and wrong for the 2 whose paired array upstream
builds differently from its merged lambda. **The remaining defect is the absence
of that array**, not a wrong sampling helper. Building it explicitly -- pairing
the treatment and control arrays as `__chrom_pair_treat_ctrl` does -- is the fix,
and it is the same shape as the F80 union fix that closed the position list.

**Gates unchanged:** `rust 155 / upstream 155, bnd_ok 155, summits 155,
worst_bnd 0, 179 at full parity`; lambda sweep 70 of 72 within 0.005 %.

## F108 -- the paired arrays need upstream's pointer walk; first attempt regressed and was reverted

F107 concluded the port lacks a `ctrl_d_pileup_d` counterpart. Reading
`__chrom_pair_treat_ctrl` (`CallPeakUnit.py:646-752`) gives the algorithm exactly,
and it is not what the port approximates. Upstream walks the two input arrays
**together**, emitting one row per step:

```python
while it < lt and ic < lc:
    if   t_p[it] < c_p[ic]:  emit (t_p[it], t_v[it], c_v[ic]); it += 1
    elif t_p[it] > c_p[ic]:  emit (c_p[ic], t_v[it], c_v[ic]); ic += 1
    else:                    emit (t_p[it], t_v[it], c_v[ic]); it += 1; ic += 1
```

The emitted position is the **smaller** of the two current positions, and the
control value is the one at the control pointer -- which may be the run whose end
is the *next* control breakpoint after the emitted position. The answer therefore
depends on whether the arrays are in step, so no per-position lookup reproduces it.
That is why F107 measured `containing` right on one fixture and `at_end` right on
seventy: neither is wrong, each is a snapshot of the walk at a different alignment.

**The walk was implemented, and it regressed badly. It has been reverted.**

```text
before (union + lookups)   rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
walk attempt               rust 155  upstream 264  bnd_ok 0    summits 152  worst_bnd 33
```

`upstream 264` against `rust 155` means the port started calling peaks upstream
does not, so the walk as first written is wrong in a way the gate catches
immediately. Two candidate causes, both identifiable from the algorithm rather than
guessed:

1. **Truncation.** The loop ends when either array is exhausted, so the position
   list is cut short wherever the treatment and control arrays differ in length.
   That changes which positions exist at all, and therefore which clear the cutoff.
2. **Repeated treatment values.** A step that advances only the control pointer
   re-emits `t_v[it]`. Sampling the q-score track independently at those repeated
   positions is not the same thing as upstream computing the score from the emitted
   triples, so the cutoff test sees a different array than upstream's does.

Cause 2 is the deeper one: the score track must be a function of the emitted
`(pos, treat, ctrl)` triples, not of a separately-sampled q track. Reinstating the
walk means deriving the score from the triples, which is a larger change than a
substitution in the chunk loop and is not something to attempt without room to
verify it.

**Reverted, and the state is restored exactly:** `rust 155 / upstream 155`,
`bnd_ok 155`, `summits 155`, `worst_bnd 0`, 179 fixtures at full parity; lambda
sweep 70 of 72 within 0.005 %; 333 tests passing, 0 clippy, fmt clean.

**The lesson is F107's, applied to myself.** The walk looked like a mechanical
substitution -- twenty lines replacing twenty lines -- and it was wrong in a way
that inverted a perfect gate. Getting F107's *diagnosis* right did not make F108's
*implementation* right, because the diagnosis described the control array while the
implementation also changed the position list and the score array's inputs. When a
finding names one stage, check what else that stage feeds before editing it.

## F109 -- the gate now drives `--call-summits`, which shows a large summit gap

The corpus has 425 golden variants in each of `default`, `call_summits` and
`broad`, but `run_peak_e2e.sh` only ever drove `default` -- so an entire
acceptance mode was untested while the gate reported green. `run_peak_e2e.sh
--variant <name>` now selects the variant, reading both the XLS and `command.json`
from it, so the existing machinery covers all three.

**`--call-summits`, `poisson`, 288 fixtures:**

```text
TOTAL (poisson/call_summits)   rust 166  upstream 166  bnd_ok 160  summits 19  worst_bnd 0
  fixtures compared : 288
  fixtures with full parity : 65
```

**Peak counts match exactly and the boundaries are exact** -- `worst_bnd 0` means
not one reported peak boundary is displaced. The `bnd_ok` of 160 against 166 is
count bookkeeping across fixtures (six peaks exist on one side only), not
displacement.

**Summits are the gap: 19 of 166.** That is the expected consequence of F66's
finding, and of the honest caveat that has been carried since F85: the
`--call-summits` path was *aligned* to the F85/F86 conventions but never
*verified*, and alignment is not correctness. The summit mode does not take the
midpoint of the max-pileup chunk at all -- it writes the region's pileup into a
dense array, finds Savitzky-Golay smoothed maxima (`maxima`), and filters them
through `enforce_peakyness`, with three separate fallbacks to the non-summit path.
None of that has been exercised against the oracle until now.

So the real state of G9/G10 is narrower than "narrow peaks closed": **narrow
single-end and paired-end are closed on every coordinate including summits;
`--call-summits` has exact counts and boundaries but a summit implementation that
is not yet correct; `--broad` has not been driven at all.**

The next step for `--call-summits` is a stage-by-stage comparison of the summit
path alone -- the dense array, the smoothed maxima, and the `enforce_peakyness`
filters -- using the same capture-then-compare method that closed the narrow path,
with the narrow path held fixed as the control.

## F110 -- `--call-summits` summit index is +1; 160 of 166 now exact

F109 found `--call-summits` reporting exact counts and boundaries but only 19 of
166 summits, every one of them exactly **1 low**. The summit in that mode is not a
chunk midpoint: the region's pileup is written into a dense array over
`[peak_start - 10, peak_end + 10]`, Savitzky-Golay smoothed maxima are taken, and
the offsets are mapped back as

```rust
summit: start + *offset as Coord,
```

so the coordinate is `array_index + array_origin` with no offset correction --
whereas the narrow path's summit carries F86's `+1` from the XLS start convention.

Adding that single `+ 1` moves `--call-summits` from 19 to **160 of 166**, with
narrow mode unchanged (`bnd_ok 155 / summits 155 / worst_bnd 0`).

```text
TOTAL (poisson/call_summits)   rust 166  upstream 166  bnd_ok 166  summits 160  worst_bnd 0
TOTAL (poisson/default)        rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
```

**The six remaining fixtures overshoot in the opposite direction** -- one
summit 1 *high* rather than 1 low -- and all six are paired-end with a control
(`w600_ctrl_*`). That asymmetry is the useful part: a uniform index shift would
move every summit the same way, so the `+1` is not a blanket correction but a
convention that applies at the boundary cases differently. The likely cause is
plateau handling in `sg::maxima` -- when the smoothed maximum sits on a flat run,
which end upstream reports depends on whether the plateau is clipped by the
region edge, and a constant `+1` gets that wrong in the opposite direction.

**This is now a small, well-localised gap rather than a broken mode**, and the next
measurement is a single fixture's smoothed-maximum offsets on both sides:
instrument `maxima` in the oracle for one of the six, compare the offset arrays,
and settle whether the discrepancy is the plateau rule or the padding clip.

## F111-F114 -- `--broad` implemented; all three acceptance modes now driven

F109's caveat that `--broad` was "a gate-configuration item rather than an
algorithmic one" was **wrong**, and worth correcting: the library had no broad
implementation at all. `grep broad crates/macs-peaks/src` found nothing outside a
comment. Broad mode is a real feature and it is now built.

**What broad mode is** (`CallPeakUnit.py:1888-1966`, `PeakDetect.py:255-268`):
the region set is called **twice** -- once at the q-value cutoff into `lvl1`
(strong peaks) and once at `--broad-cutoff` into `lvl2` (broad regions) -- with
`lvl1_max_gap = maxgap` and `lvl2_max_gap = maxgap * 4`, and the two are then
combined, each lvl2 region collecting the lvl1 peaks inside it to form the
hierarchical gappedPeak blocks. The XLS columns all come from the lvl2 region.

**Four things were wrong, each found by measurement rather than reading:**

* **F111 — the XLS parser.** The broad XLS has *no `abs_summit` column*:
  `chr start end length pileup -log10(pvalue) fold_enrichment -log10(qvalue) name`.
  The fixed column index read `pileup` as the summit, so every broad comparison
  was against garbage. `read_xls` now detects the header.
* **F111 — the broad close.** `close_peak_for_broad_region` reports the start as
  `tstart + 1`, the **same** XLS convention as narrow (F112), every column is a
  **length-weighted mean** via `mean_from_value_length(value, tend - tstart)`, and
  there is **no summit** (upstream stores `0`). Peak carries a `broad` flag so the
  two paths cannot be confused.
* **F113 — `lvl2_max_gap` is `maxgap * 4`.** Using `maxgap` made every broad
  region too narrow.
* **F114 — the broad level needs its own above-cutoff chunk list.** Re-filtering
  the narrow list cannot recover positions that never entered it, so lvl2 regions
  could never extend past the narrow cutoff -- the opposite of what
  `--broad-cutoff` means. `chunks_at(cutoff)` now builds a list per cutoff.

**Result, `--broad`, `poisson`, 288 fixtures:**

```text
TOTAL (poisson/broad)   rust 156  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
  fixtures with full parity : 178
```

155 of 156 broad peaks match exactly in start and end, with **worst boundary gap
0**. One peak is over-called (`sweep/gmini_mpe_d400_w180_ctrl_192`).

**All three acceptance modes are now driven, and none regressed:**

```text
TOTAL (poisson/default)        rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
TOTAL (poisson/call_summits)   rust 166  upstream 166  bnd_ok 166  summits 160  worst_bnd 0
TOTAL (poisson/broad)          rust 156  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

**A tooling bug worth recording:** the first broad `--all` run reported
`bnd_ok 24 / worst_bnd 576` and a `NameError` traceback per fixture -- the flag
extractor contained a dict literal with single quotes inside a shell single-quoted
string, which terminated it and turned the keys into bare names. The `24/576`
number was therefore garbage from a broken comparator, not a real regression, and
it took the traceback in the output to notice. A differential that errors should
fail loudly rather than emit numbers; this one printed a plausible table anyway.

## F115 -- upstream's literal summit formula does NOT reproduce the corpus, and that is informative

F110's `+1` on the summit was empirical: summits went from 19 of 166 to 160.
Reading upstream says it should not be there. `__close_peak_with_subpeaks`
(`CallPeakUnit.py:1475-1570`) is unambiguous:

```python
peak_start = peak_content[0][0]          # no +1
start      = max(peak_start - 10, 0)
...
summit=start + summit_offset             # no +1
```

and `maxima` (`SignalProcessing.py:29-44`) is reproduced faithfully in
`crates/macs-peaks/src/sg.rs` -- `window_size//2*2+1`, `.round(16)`, `np.sign`,
`np.diff(sign) <= -1` yielding the *left* index of the falling edge -- so the
offset arrays are the same function of the same input.

**Implementing the formula literally is much worse:**

```text
peak_start = region[0].start,  summit = start + offset   (upstream, verbatim)
   -> bnd_ok  17 of 166, summits 17 of 166
current (peak_start = region[0].start + 1, summit = start + offset + 1)
   -> bnd_ok 166 of 166, summits 160 of 166
```

**Both adjustments are load-bearing, and they cannot both be right if the chunk
list is upstream's.** Two `+1`s that cancel in the summit but not in the reported
start means one of them is compensating for a systematic one-base difference
*before* the summit arithmetic -- namely that this port's chunk
`start = pos[i-1]` is **one lower than upstream's
`above_cutoff_startpos`**. F80 established that upstream's
`above_cutoff_startpos = pos_array[above_cutoff-1]`, which is the same expression,
so either `pos_array` here is not upstream's `pos_array`, or the two ports
disagree about which index `above_cutoff` names.

The likelihood is that upstream's `above_cutoff` indexes the **paired** array
`pos_array` (F35/F80/F107's array-identity distinction, now biting for a fourth
time), while this port's `i` indexes the **union** position list built in the
harness. Where the two arrays agree the offsets cancel; where they differ by one
position -- 6 of 166 summits -- they do not, which is exactly the observed
pattern: 160 correct, 6 overshooting in the same direction, all on
paired-end-with-control fixtures.

**So F110's `+1` is a compensating correction, not the rule**, and it should be
replaced once the chunk index is taken from the paired array rather than the union
list -- the same fix F108 identified for `fold_enrichment`, and now shown to
affect the summit path too. F108's attempt to build the paired arrays by upstream's
pointer walk regressed for the two reasons recorded there (early truncation, and a
q-score track that is not a function of the emitted triples); this finding is
independent evidence that the same root cause is behind three separate residuals.

**Current state across all three modes, unchanged and verified:**

```text
TOTAL (poisson/default)       rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
TOTAL (poisson/call_summits)  rust 166  upstream 166  bnd_ok 166  summits 160  worst_bnd 0
TOTAL (poisson/broad)         rust 156  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
```

## F116 -- the CLI has a real, derived argument surface

`macs-cli` was a list of 14 names with an empty `main.rs`, and
`macs-bedgraph`, `macs-model`, `macs-hmmratac`, `macs-callvar` were one-line
placeholders. `oracle/flag_matrix.tsv` -- the auto-derived argparse surface, 315
rows across all 14 subcommands -- was already the authoritative list, so the
parser could be **derived from it** rather than hand-maintained.

`crates/macs-cli/src/flags.rs` embeds that matrix with `include_str!` (so the
shipped binary keeps no runtime data dependency) and parses against it with
argparse's observable semantics:

```text
$ macs3-rs --version
macs3-rs 3.0.5
$ macs3-rs badcmd
macs3-rs: error: argument command: invalid choice: badcmd
(choose from callpeak, bdgpeakcall, bdgbroadcall, bdgcmp, bdgopt, cmbreps, bdgdiff, filterdup, predictd, pileup, randsample, refinepeak, hmmratac, callvar)
$ echo $?
2
$ macs3-rs callpeak
usage: macs3-rs callpeak [-h] ...
macs3-rs callpeak: error: the following arguments are required: -t
$ echo $?
2
```

Covered: `--flag=value` and `--flag value`, short forms, defaults applied before
parsing, `store_true`/`store_false`/`count`, `choices` rejection, missing required
flags, unknown flags, `--help`, and the exit-code convention (2 for usage, 1 for
everything else). 17 tests in `macs-cli` exercise these, including one that
asserts every one of the 14 subcommands has rows in the matrix -- so a subcommand
cannot be added to the list without the surface following.

**What this does not do, stated plainly.** No subcommand writes output yet,
because the verified peak-calling pipeline lives in the *differential harness*
(`crates/macs-peaks/src/bin/callpeak_e2e.rs`), not in the library. That is the
real structural blocker for all 14 commands, and it is worth naming as such: the
oracle-verified implementation is currently reachable only through a test binary.
Promoting it -- moving the track construction, the score/q stage and the peak
caller into `macs-peaks` behind a `run_callpeak(options)` entry point, with the
harness reduced to a thin wrapper -- is the next structural step, and every
command after `callpeak` needs the same shape.

Two implementation notes: `FlagSpec` originally borrowed `&'static str` fields
from the embedded matrix via an interning helper, which forced a `Mutex<Vec>`
into a data-only table for no benefit, so the fields own `String`s instead; and
an early version of the flag extractor in `run_peak_e2e.sh` put a dict literal
with single quotes inside a shell single-quoted string, which silently turned the
keys into bare names and made the broad comparator report `bnd_ok 24 / worst 576`
(F111). Both were caught by reading the emitted output rather than the code.

## F117 -- the XLS/narrowPeak writers, with `%.5g` checked against CPython

F116 established the CLI's flag surface but no command could write output,
because no writer existed at all -- `grep` for an XLS writer across `crates/`
returned nothing. `crates/macs-io/src/peakout.rs` now provides them, transcribing
`MACS3/Commands/callpeak_cmd.py`:

* `xls_row` / `xls_body` -- the XLS, with `length = end - start + 1` (inclusive)
  and the broad header that **omits `abs_summit`** (F111);
* `narrowpeak_row` -- 0-based half-open, `chromStart = xls start - 1`,
  `chromEnd = xls end`, scores as *signed* `-log10(p)` with five decimals;
* `summit_row` -- the `*_summits.bed` line;
* `format_g` -- the `%.5g` the score columns use.

**`format_g` is the part that needed care, and it is exactly the kind of defect a
numeric comparison cannot see.** The acceptance criterion is *byte-identical*
XLS, and the score columns are formatted rather than printed, so a plausible but
wrong transcription produces numerically identical output and a differing file.
The transcription is Python's rule -- fixed notation while
`10^-4 <= |v| < 10^sig`, exponent form outside, trailing zeros stripped, and the
exponent written `e+05` with a two-digit minimum -- and it is checked by
`crates/macs-io/tests/fmt_golden.rs` against 13 strings produced by CPython's own
`"%.5g" %`, plus the three boundary cases (exponent below and above the range,
and the last fixed-notation value `99999`).

That test would have caught the two most likely transcription errors: Rust's
`{:e}` writes `1.2346e-5` where Python writes `1.2346e-05`, and Rust does not
strip trailing zeros the way `%g` does.

**Nothing downstream regressed:**

```text
TOTAL (poisson/default)       rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
TOTAL (poisson/call_summits)  rust 166  upstream 166  bnd_ok 166  summits 160  worst_bnd 0
TOTAL (poisson/broad)         rust 156  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
```

The remaining structural step for a working `macs3-rs callpeak` is now just the
assembly: load the inputs, build the tracks, run the library's existing
`pair_treat_ctrl` + `call_chromosome_paired`, and write through these writers.
Every piece of that exists; only the glue is missing, which is the F116 finding
confirmed from the other end.

## F118 -- attempt to lift the pipeline into the library, reverted

F117 left the assembly glue as the only thing between the verified pipeline and
a working CLI. The move was attempted: extract the ~298-line per-chromosome body
of `finish` from `bin/callpeak_e2e.rs` into `macs-peaks/src/pipeline.rs` behind
`call_peak_chromosomes(sig, table, qtracks, options, d)`, with a `ChromInput`
and a `CallpeakOptions`, then let both the CLI and the harness call it.

**It was reverted.** The extraction is mechanical in principle and was not
mechanical in practice: the body references harness-only items -- the
`--peak-content` oracle override, `ChromSignals.chrom`, the `CALLPEAK_*` debug
probes, `broad_level`/`combine_broad`, and the `macs_peaks::` self-import that
only works from inside the crate's own bin -- and stripping each one in sequence
left the file inconsistent three attempts in. By the fourth the file had lost two
closing braces to the scripted surgery itself, and continuing was no longer
defensible with the remaining budget.

**Reverted and the tree is verified back to green:**

```text
TOTAL (poisson/default)       rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
TOTAL (poisson/call_summits)  rust 166  upstream 166  bnd_ok 166  summits 160  worst_bnd 0
TOTAL (poisson/broad)         rust 156  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
cargo test: 347 passed, 0 failed; 0 clippy; fmt clean
```

**What the attempt established, which is worth keeping.** The glue is *not* a
thin wrapper. The pipeline body depends on harness-only affordances, and three of
them are load-bearing for the gates rather than incidental:

* `--peak-content` substitutes upstream's captured `peak_content` for the chunk
  list. That is an **oracle hook** and must not ship, but the gates' `--peak-content`
  mode is how F60-F66 isolated the summit stage, so the hook has to stay in the
  harness even after the production path is in the library.
* The `CALLPEAK_*` probes (`CHUNKDUMP`, `NEAR`, `SCALEVAL`, `ATSCALE`, `PECTRL`,
  `TREATDUMP`) are **how every finding from F72 to F107 was localised**. Deleting
  them with the move would remove the instrument that found the bugs.
* `s.chrom` exists because the harness interns chromosomes from whichever input
  track came first; the library needs that resolved before the move, not during.

So the right shape is: the *production* path moves into `macs-peaks`, and the
oracle probes stay in the harness as a separate, thinner layer that can still
substitute captured upstream inputs. That is a deliberate split rather than a
verbatim move, and it needs room to do properly.

## F119 -- byte-comparing the XLS, which immediately found a q-score coverage gap

F117 added the writers. F119 uses them: the harness renders its peaks through
`macs_io::peakout` under `CALLPEAK_WRITE_XLS`, and the result is diffed against
the golden XLS. This is the `byte-identical *.xls` acceptance criterion, and it is
the first test that can see column order, the inclusive `length`, the number
format, the header and the peak-name prefix.

**Two things had to be corrected in the writer first.**

* Upstream writes the score columns with **`%.6g`**, not `%.5g`
  (`MACS3/IO/PeakIO.py:815-840`). My `%.5g` rendered a fold enrichment of
  `7.20635` as `7.2064`. F117's golden test used the *right* values at the *wrong*
  precision -- it validated the transcription against CPython but chose the wrong
  precision to validate. The test now carries both: 13 cases at `%.5g` and 8 real
  golden values at `%.6g`.
* `pileup` is `round(pileup, 2)` *before* formatting.

The peak-name prefix is `<name>_peak_<n>`, from `name_prefix % name` with
`name_prefix = "%s_peak_"` -- so it depends on `-n`, which the gate does not
currently replay.

**The result is close, and it found a real defect the coordinate diff structurally
cannot see.** Five fixtures, data rows only (headers are upstream's, which this
port does not yet reproduce):

```text
IDENTICAL  se_basic/gauss_two_peaks
IDENTICAL  se_shapes/bias
differs    se_model/realistic       1 of 12 rows
differs    se_model/spikes_only    3 rows
differs    tiny/two_contigs        1 row
```

`realistic` is **11 of 12 rows byte-identical**, and the single differing row is
peak 8 -- whose summit chunk is one of the two with the known lambda residual
(F107/F108), so the coordinates match and only the score columns differ. That is
the byte-level confirmation that F115's diagnosis was right: the coordinates are
exact, the lambda behind one summit is not.

**`spikes_only` found something new:**

```text
mine     chr20 15254 15792 539 15610 8 10.5485  6.99203 0        spikes_only_default_peak_2
upstream chr20 15254 15792 539 15610 8 10.685   7.04753 5.71769  spikes_only_default_peak_2
```

`-log10(qvalue)` is **0** where upstream has `5.71769`. A q-score of zero is the
documented "no entry in the table" fallback, so this port's q-score histogram is
**missing p-score keys** that upstream's table contains, at those summits. The
p-score column differs too (`10.5485` vs `10.685`), so the summit chunk's
lambda differs -- but the *zero* is the sharper signal: whatever p-score this port
computed, its canonical key was not in the table, which means the histogram was
built from a different position set than the one the peak caller reads. That is
the same array-identity family as F107/F108/F115, now visible in a fourth place
and the first one that affects a *written column*.

The next step is therefore concrete: dump the p-score keys at that summit from
both sides and find which are absent, rather than reasoning about the lambda
again.

## F120 -- the q-table and the peak caller read *different* p-score arrays

F119 found `-log10(qvalue)` written as `0` on `se_model/spikes_only`. The
`CALLPEAK_PQ` dump resolves it completely.

```text
se_model/spikes_only, peak 2
  mine      chr20 15254 15792 539 15610 8 10.5485  6.99203 0
  upstream  chr20 15254 15792 539 15610 8 10.685   7.04753 5.71769
```

The dump of this port's own q-score track settles it:

```text
max pscore in q-track dump: 9.406459808   with qscore 5.717694759
```

**Two facts fall out.**

1. **`5.71769` is already present in this port's table** -- it is the q-score of
   the span whose p-score is `9.406`. So the q *value* upstream writes is one this
   port already knows; it is the *lookup* that fails, not the table's content.
2. **The summit's p-score is `10.5485`, which exceeds every p-score in the
   table.** `qscore_or_zero` therefore returns the documented zero fallback, and
   a `0` is written into a result column.

The p-score the peak caller uses and the p-scores the table was built from are
**not the same array**. The table comes from the `paired_pscore` track -- the
coincident-boundary walk over the treatment and control tracks (F59) -- while
`close_peak_wo_subpeaks` recomputes `pscore(treat[chunk], ctrl[chunk])` from the
chunk's own values, which are sampled from the union position list (F80). Where
those two agree the lookup succeeds and everything matches; where they differ the
summit's p-score can fall **outside the whole table**, and the failure is not a
tolerance question but a missing key.

That is the F107/F108/F115 array-identity root cause with the sharpest possible
symptom: it is the only failure mode found so far that puts a visibly wrong
number in an output file rather than shifting a coordinate by a base. It also
explains why the residual looked like "a small lambda difference" -- the lambda
difference is real (F107), but the *consequence* in this fixture is a p-score
outside the histogram.

**The fix follows directly and is not optional bookkeeping:** the chunk's treatment
and control values, the score track, and the table must all be derived from the
same paired arrays. Concretely, `close_peak_wo_subpeaks` should read the summit's
score out of the table by index rather than recomputing `pscore` from independently
sampled values -- so that a peak can never report a score the histogram never saw.
Upstream has no such failure mode because it has exactly one array:
`d_pileup_d` / `ctrl_d_pileup_d`, from which the score track, the histogram and the
peak columns are all derived.

**Verification, once made:** `se_model/spikes_only` and `tiny/two_contigs` should
join `gauss_two_peaks` and `bias` as fully byte-identical, and `realistic` should go
from 11 of 12 rows to 12 of 12.

## F121 -- the summit-score-by-index path exists but is wrong as wired; the real defect is the score track's range

F120 proposed taking the summit's p/q from the score track by index, so a peak can
never report a score the histogram never saw.
`close_peak_wo_subpeaks_with_p` and `call_peaks_chromosome_with_p` implement that,
and the harness builds the per-position p array. **Wiring it in regressed the gate
from 155 to 264 peaks**, so the harness keeps the proven recomputation and the new
entry points sit unused. Recorded as a partial implementation with a known reason,
not a fix.

**The measurement that explains why, and it is bigger than anything found so far.**

```text
CALLPEAK_HIST dump of this port's p-score track, se_model/spikes_only
  runs 4783   max p 5.9993   runs >= 9.5: 0
upstream peak 2 of that fixture: -log10(p) 10.685, -log10(q) 5.71769
```

**This port's p-score track tops out at 6.0 on a fixture where upstream's reaches
10.685.** So the summit's p-score of 10.55 is not merely absent from the
histogram -- it is *more significant than any value the score track contains
anywhere*. That is not a table-coverage problem; it means the score track is built
over a position set whose p-values are systematically lower than the chunks'.

The cutoff still works (155 of 155 peaks, `worst_bnd 0`) because the cutoff is
`-log10(0.05) = 1.301`, and everything above that agrees. The divergence only
becomes visible where the histogram is asked for a *high* p-score, which is exactly
what the summit column does. That is why F107's lambda residual looked like 0.16 %
for so long: it was measured on fixtures whose peaks sit low in the score range,
where a 0.16 % lambda error moves the p-score by a thousandth.

**So F120's diagnosis needs correcting in its mechanism.** The "two arrays" framing
is right -- the chunks and the score track are built from different position sets --
but the failure is not that the lookup misses a key that is nearly present. It is
that the two disagree by a factor that grows with signal strength, and the summit
column is the first place that samples the top of the range.

That also explains why `tiny/two_contigs` writes `-log10(q) 0` with
`-log10(p) 19.3145` (upstream `18.8895` / `16.685`): a p-score of 19 is far outside
any plausible score-track range here.

**The next measurement is therefore about the score track, not the table:** dump
this port's p-score track and upstream's `pscore_track` for one fixture and compare
their maxima and their value at a peak summit. If upstream's track genuinely reaches
10.685 where this port's reaches 6.0, the paired walk that builds the score track
(F59) is dropping or under-sampling the highest-signal positions -- and that would
also explain F107, F115 and the `--call-summits` residuals as one defect.

**Gates unchanged and verified:**

```text
TOTAL (poisson/default)       rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
TOTAL (poisson/call_summits)  rust 166  upstream 166  bnd_ok 166  summits 160  worst_bnd 0
TOTAL (poisson/broad)         rust 156  upstream 155  bnd_ok 155  summits 155  worst_bnd 0
e2e pileup: 97 passed, 0 failed, 3 skipped (upstream rejects) (of 100)
cargo test: 348 passed, 0 failed; 0 clippy; fmt clean
```

## F122 -- ROOT CAUSE: upstream has no separate score track; both the histogram and the peak caller read the paired arrays

F121 measured that this port's p-score track tops out at 6.0 where upstream's
reaches 10.685, and guessed the paired walk was dropping high-signal positions.
Reading the source shows something simpler and more fundamental: **there is no
separate score track in upstream at all.**

`__pre_computes` builds the histogram from the paired arrays:

```python
# MACS3/Signal/CallPeakUnit.py:896-935
[pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl
score_array = self.__cal_pscore(treat_array, ctrl_array)
for n in range(len(tmplist)):
    ...
    above_cutoff = np.nonzero(score_array > cutoff)[0]
```

and the peak caller computes the *same* expression over the *same* arrays:

```python
# CallPeakUnit.py:1191-1215
self.pileup_treat_ctrl_a_chromosome(chrom)
[pos_array, treat_array, ctrl_array] = self.chr_pos_treat_ctrl
for i in range(len(scoring_function_s)):
    score_array_s.append(self.__cal_pscore(treat_array, ctrl_array))
    score_array_s.append(self.__cal_qscore(treat_array, ctrl_array))
...
above_cutoff = np.nonzero(apply_multiple_cutoffs(score_array_s, score_cutoff_s))[0]
```

`__cal_qscore` then maps each position's own p-score through the table
(`CallPeakUnit.py:1629`), so a peak's q-score is the table entry for the p-score
of *its own* paired position.

**So upstream has exactly one array, and it feeds everything.** The histogram, the
cutoff, the chunks (`treat_array[ti]`, `ctrl_array[ti]`), the summit's p-score and
the summit's q-score are all derived from `d_pileup_d` / `ctrl_d_pileup_d`. Because
the table is built from the very score array the caller indexes, **a peak can never
report a score the histogram never saw** -- the failure F119 found is structurally
impossible there.

**This port has two geometries.** The histogram is fed from `paired_pscore`, a
walk over the treatment and control *tracks* emitting only **coincident** run
boundaries (F59), while the chunks sample the **union** position list (F80). At the
summit those two disagree, the chunk's p-score can exceed every value the
coincident-boundary track contains, and `qscore_or_zero` returns its zero fallback
-- which is the `0` in `se_model/spikes_only`'s XLS and the six 1 bp summit
residuals, all one defect.

**The fix is now fully specified and is an ordering problem, not an arithmetic
one.** The union position list must be built *before* the histogram, and the p-score
track must be built over that list rather than over coincident boundaries:

```text
current:  tracks -> paired_pscore -> sink/histogram      (line 724)
          finish -> union pos/tpos/cpos -> chunks        (line 1010)
required: union pos/tpos/cpos -> pscore track over it -> sink/histogram
          -> chunks, summit scores, q lookup
```

Today the union is assembled 300 lines downstream of where the histogram is fed,
which is exactly why the two could drift apart for so long without any gate
noticing -- the coordinate comparison only reads the chunks, and the byte
comparison only reads the summit columns.

`__cal_pscore` is already bit-exact (117 vectors) and `PqTable` already mirrors
upstream's histogram, so nothing new is needed arithmetically: the score track
simply has to be computed from the same arrays the chunks use.

**Verified unchanged:** 348 tests, 0 clippy, fmt clean; all three peak gates and
the signal gate still at `worst_bnd 0`.

## F123 -- MEASURED: naively unifying the geometry REGRESSES the gate; the real blocker is an exact `__chrom_pair_treat_ctrl`

F122 says upstream feeds one geometry to the histogram, the cutoff and the
summit. The obvious way to realise that here is to build the histogram from the
same union arrays the chunks already use (`containing`/`at_end` sampling) and
look the q-score up per position. I implemented exactly that, behind the
existing sink, and measured it:

```
before (retadd_pscore histogram + containing/at_end chunks)
  rust 155  upstream 155  bnd_ok 155  summits 155  worst_bnd 0    (baseline)

after  (histogram built from the containing/at_end union, qpos = table[g.p[i]])
  rust 155  upstream 155  bnd_ok  23  summits 136  worst_bnd 32
```

**Peak count is unchanged at 155, but boundary agreement collapses 155 -> 23
with peaks off by up to 32 bp.** So the above-cutoff position set moved. That is
diagnostic, not cosmetic: the chunk cutoff is `qpos[i] > -log10(qvalue)`, and
`qpos` is now the table lookup of `pscore(containing(treat), at_end(ctrl))`
instead of the q-score sampled from the `retadd_pscore` track. The two p-score
geometries differ enough to move which positions clear `-log10(0.05) = 1.301`.

**This tells us which of the two approximations is closer to upstream.**
`retadd_pscore` -- the `min(p1,p2)` pointer walk over the two tracks -- reproduces
upstream's above-cutoff positions exactly (155/155). The
`containing`/`at_end` union does not. So `retadd_pscore` is the better model of
`__chrom_pair_treat_ctrl`'s `treat_array`/`ctrl_array` (and therefore of
`score_array`), even though F83/F84 showed `containing`/`at_end` gives the right
chunk *pileup values*.

That is the crux: this port currently approximates the single upstream array
**two different ways**, and each is right where the other is wrong:

```
                    position set        chunk treat/ctrl value
retadd_pscore       upstream-exact      (not used for chunks)
containing/at_end   diverges            bit-exact pileup (F83/F84)
```

The summit's p/q columns are stuck between them. `close_peak_wo_subpeaks`
recomputes `pscore(tpos[i], cpos[i])` from the *chunk* geometry, while the
histogram holds *retadd*-geometry p-scores, so the two can disagree by more than
the histogram's range -- exactly the F119 `0` and F121 `6.0 vs 10.685` symptom.
Patching the histogram to the chunk geometry makes the summit's q well-defined but
moves the cutoff, so it cannot be done in isolation.

**The correct fix is not histogram plumbing: it is implementing
`__chrom_pair_treat_ctrl` exactly**, so one array serves the cutoff selection,
the chunk values and the summit scores, and every consumer indexes that array.
F108's reverted pointer walk was the right idea but terminated when either array
ran out and repeated a treatment value across control-only steps, changing the
above-cutoff set. Reinstating it is now worth doing, because unlike F108 there is
a concrete target: it must reproduce `retadd_pscore`'s positions exactly (that
is the measured upstream-exact position set) **and** F83's `containing` values for
the chunk piles, rather than being accepted or rejected on a peak count.

Until that exists, the callpeak narrow coordinate gate is at its practical floor
(155/155, worst_bnd 0) and the residual defects are: the summit p/q columns in
the XLS, six `--call-summits` summit coordinates, and one broad overcall -- all
downstream of the same unimplemented paired-array walk.

**Verified after revert:** 155/155, `bnd_ok 155`, `summits 155`, `worst_bnd 0`.

## F124 -- G6 ported: `bedGraphTrackI.overlie` is bit-exact against upstream for every reachable function

G6 (`macs-bedgraph`) is now a real port of `MACS3/Signal/BedGraph.py`
(`bedGraphTrackI`) rather than a one-line placeholder, and `overlie` is verified
against the pinned oracle rather than only against hand-written expectations.

A bedGraph is a per-chromosome list of `(end, value)` pairs -- `value[i]` covers
`[end[i-1], end[i])` -- which is exactly `SignalTrack`'s run representation. So
the upstream pointer walks port directly onto the run slices, and `overlie`'s
union-breakpoint loop is:

```text
lowest = min(current end of each input)
emit [pre, lowest) with op(current values)      # merging equal neighbours
advance every input sitting at `lowest`
stop as soon as any input is exhausted            # upstream StopIteration
```

**Verification.** For this input pair (dumped from upstream with `%.6g`):

```text
a.bdg: chr1 0 100 0 / 100 200 3 / 200 300 4
b.bdg: chr1 0 150 1 / 150 250 2 / 250 300 4
```

upstream produces, and this port now reproduces exactly:

```text
max      1,3,4                     (three runs, not five)
sum      1,4,5,6,8
subtract 1,-2,-1,-2,0
product  0,3,6,8,16
mean     0.5,2,2.5,3,4
fisher   0.481146,2.99096,3.90264,4.82928,6.71174
```

Two details the differential nailed, both easy to get wrong:

* **`max` yields three runs, not five.** Upstream's docstring lists the five
  intervals `1,3,3,4,4`, but `overlie` writes through `add_loc`, which
  coalesces an equal-valued predecessor -- so the two `3`s (and the two `4`s)
  merge. A run-length representation must merge too, or the bedGraph written out
  has spurious breakpoints.
* **`divide` is not reachable with two tracks.** `divide_func` is `x[1]/x[2]`
  (`BedGraph.py:74-77`), but `overlie` passes a value tuple of length 2 when
  there is one overlaid track, so upstream raises `IndexError: list index out of
  range`. Reproducing a crash is not useful, and the "zero panics" criterion
  forbids it, so this port returns `0.0` for the under-length case and matches
  upstream on every function that actually runs. `subtract` is likewise
  `x[1] - x[0]` and is **not** "corrected" to `v[0] - v[1]`.

**A real cross-track bug surfaced while porting.** `ChromId` is an index into a
per-input `Genome`, so comparing ids across two different `BedGraph`s is
meaningless: chromosome `x` from file A and chromosome `y` from file B both get
id `0`. The first implementation intersected the id sets and therefore emitted a
merged `chr`-less track whenever the two files shared *no* chromosome names but
both used id 0. `overlie` now intersects by chromosome **name** and resolves each
name to an id per track. The same hazard applies anywhere a `ChromId` crosses a
`Genome` boundary, so it is worth remembering for the remaining commands.

**Still open for a full G6:** `ScoreTrack.TwoConditionScores` (bdgdiff),
`make_ScoreTrackII_for_macs` (bdgcmp), and `call_broadpeaks`/`refine_peaks`
over bedgraphs.

**Verified:** 358 workspace tests pass (10 new in `macs-bedgraph`), 0 clippy,
`cargo fmt` clean; the callpeak default gate is unchanged at
`155/155, bnd_ok 155, summits 155, worst_bnd 0`.

## F125 -- `callpeak` now runs end-to-end and writes byte-identical files; three narrowPeak/summits writer bugs fixed

Two things landed, and the second one was found only because the first forced a
real byte comparison.

### The production path exists and the harness shares it

`macs-peaks/src/callpeak.rs` now owns the per-chromosome core
(`call_chromosome`) and the full single-end pipeline (`run_callpeak_se`), and the
differential harness's normal path *calls that library function* rather than its
own inline copy. The oracle-only `peak_content` injection path stays in the
harness, since it substitutes an upstream capture and has no production
equivalent. Consequence: the harness and a shipping binary cannot drift, and the
golden gates exercise the shipped code directly. All three modes still measure
identically after the lift (`155/155`, `166/166`, `156/155`, `worst_bnd 0`).

`crates/macs-peaks/src/bin/callpeak.rs` is the first command that produces files:

```text
$ macs-callpeak -t treat.bed -c ctrl.bed -f BED -g 2000000 --nomodel \
      --extsize 200 --outdir out -n realistic_default
callpeak: 12 peaks written to out
$ diff realistic_default_peaks.narrowPeak out/realistic_default_peaks.narrowPeak
```

**One row differs** (peak 8), and only in the summit `-log10(p)`/`-log10(q)`
values -- coordinates, summit, fold enrichment and every other peak are
byte-identical. That row is the known F123 `__chrom_pair_treat_ctrl` residual.

### Three real bugs in the narrowPeak / summits writers

F119 byte-compared the **XLS only**, so these were never exercised. Reading the
real golden files against the writers exposed all three:

1. **score column.** The writer emitted `pileup`. Upstream writes
   `int(10 * peak[score_column])`, and `score_column` is `pscore` when
   `--log_pvalue` is given, else `qscore` (`callpeak_cmd.py:293-298`). The golden
   shows `... 79 . 7.20635 ...` where `int(10 * 7.9977) == 79` -- the qscore, not
   the pileup.
2. **summit column.** The writer emitted an absolute summit. Upstream writes the
   **offset** `peak['summit'] - peak['start']` (`PeakIO.py:760`). The golden's
   `160` for a peak spanning `7571..8039` with summit `7732` is exactly that
   offset, not a coordinate.
3. **signedness.** `narrowPeak`/`summits` p/q columns were written `-X.XXXXX`.
   Upstream writes plain `%.6g` (`PeakIO.py:761-770`); the golden has `10.828`
   and `7.9977`, unsigned.

`summits.bed` was also missing its score column entirely (4 columns instead of
5). All three files are now pinned byte-for-byte to the real golden in unit
tests (`narrowpeak_row_is_byte_identical_to_golden`,
`summit_row_is_byte_identical_to_golden`).

Two more upstream details came out of the same exercise:

* **XLS peak numbering also groups by end.** `write_to_xls` uses the same
  `groupby(end)` + `subpeak_letters` scheme as narrowPeak/summits, so a
  `--call-summits` peak's sub-peaks are `name_peak_1a`/`1b` in the XLS too. The
  writer numbered per-row, which is identical for the default mode (every peak
  has a unique end) but wrong for call-summits. All three writers now share one
  grouping.
* **`subpeak_letters(26)` is `"ba"`, not `"aa"`.** The upstream docstring
  advertises `a..z, aa, ab, ...`, but the implementation is
  `subpeak_letters(i // 26) + chr(97 + i % 26)`, which yields `ba`. The code is
  the contract; the test pins the real behaviour.

### A default-flag bug caught by the byte compare

The first production run disagreed with the harness on peak *end* (8068 vs
8039). The cause was mine: I had defaulted `--slocal 0 --llocal 1000`, while
upstream's argparse defaults are **1000 and 10000** (confirmed against the
auto-derived `oracle/flag_matrix.tsv`). A different local-lambda window changes
the merged lambda, hence the q-track, hence which positions clear the cutoff.
Fixed. Worth remembering: an unstated default is a numerical parameter here, not
a cosmetic one.

**Verified:** 363 tests pass (5 new), 0 clippy, fmt clean; default/call_summits/
broad gates unchanged; the signal gate is 97/97.

## F126 -- `callpeak` is wired into the unified `macs3-rs` CLI; a real argparse bug fixed

`callpeak` moved from a standalone `macs-callpeak` binary into
`macs-cli/src/commands/callpeak.rs`, dispatched from the `macs3-rs` entry point.
A standalone binary is not drop-in CLI compatibility -- the goal is `macs3
callpeak ...`, which now works and is covered by
`crates/macs-cli/tests/callpeak_cli.rs` (a real `parse_flags` call feeding the
shared pipeline and asserting all three files appear with the right column
counts).

**The derived-parser choice bug.** The first `macs3-rs callpeak -f BED` failed:

```text
macs3-rs callpeak: error: argument -f: invalid choice: BED
    (choose from AUTO,BAM,SAM,BED,...,FRAG)
```

`BED` is *in* the list. The matrix stores choices **comma**-separated
(`oracle/flag_matrix.tsv`, field `choices`) but the parser split on `|`:

```rust
let allowed: Vec<&str> = spec.choices.split('|').collect();  // one bogus element
```

so every legitimate value was rejected. This was latent because nothing had
exercised a `-f/--format` choice end-to-end through the CLI before. Fixed to
`split(',')`. Worth a dedicated test now exists
(`a_valid_choice_flag_is_accepted_by_the_derived_parser`), which also pins that a
genuinely bad choice is still a usage error.

**On the p-score zero-lambda panic (F127):** the CLI end-to-end test's first
synthetic fixture had *no* control coverage anywhere near the treatment
clusters, so `pscore(treat>0, ctrl==0)` hit the degenerate Poisson and aborted.
That is not a realistic MACS input (control is sequenced genome-wide), and the
golden corpus never triggers it -- every fixture has `lambda_bg > 0`. See F127 for
the precise characterization and why the obvious fix is not safe yet.

**Verified:** 366 tests pass (3 new CLI integration), 0 clippy, fmt clean;
`macs3-rs callpeak` reproduces the golden `se_model/realistic` output byte-for-byte
except the single known F123 summit-p/q row.

## F127 -- RESOLVED: the zero-lambda abort was the *symptom*; no-control runs take a different upstream code path

`compute_uncached(k, 0.0)` aborted, which was reachable on **every** no-control
fixture: `se_edge/no_control` and all 216 `sweep/*_noc_*` fixtures. The earlier
write-up of this finding concluded the degenerate lambda was "load-bearing for
parity" and deferred it, on the evidence that a substitute value moved the narrow
gate from 155/155 to 264. That inference was wrong, and the reason it was wrong
is the real content of this finding.

### Instrumenting the abort found the exact input

Replacing the `expect` with a logger showed the calls are always
`expectation == 0.0` exactly, reached from `close_peak_wo_subpeaks_with_p` -- the
**summit** p-score recomputation, not the p-score track. Mapping every fixture
that hit it gave a clean partition: 216 of 288, every one a no-control fixture.

### Root cause: without `-c`, upstream does not use the with-control method at all

`PeakDetect.call_peaks` (`PeakDetect.py:109`) dispatches:

```python
if self.control:  self.peaks = self.__call_peaks_w_control()
else:             self.peaks = self.__call_peaks_wo_control()
```

`__call_peaks_wo_control` (`PeakDetect.py:1104`) differs in three ways, all of
which the port was getting wrong by running the with-control formulas against an
absent control:

| | with control | without control |
|---|---|---|
| control scales | `d`/`slocal`/`llocal` ladder | **one** `lregion` window |
| SE scale factor | `d/lregion` | same, but no ladder |
| PE scale factor | per-scale | `treat_length / (lregion * treat_total * 2)` |
| `treat_scale` | `1/ratio` | **hard-coded `1.0`** |
| `lambda_bg` | `control_sum/gsize` or `treat_sum/gsize` | SE `d*total/gsize`, PE `treat_length/gsize` |

The port computed `to_control = t_total > c_total`, which with `c_total == 0` is
always true, so it took `lambda_bg = control_sum / gsize = 0` -- a zero lambda at
every base -- and `treat_scale = 1/ratio = 1/0`. Hence `get_pscore(k, 0.0)`, which
upstream treats as a fatal `AssertionError`.

Fixing the branch (not the Poisson) resolves it. `se_edge/no_control` now matches
(`CALLPEAK PEAK MATCH`), and `--nolambda` is wired to the empty-scale-list case,
which is upstream's `ctrl_pv = [treat_pv[0][-1:], np.array([lambda_bg], dtype="f4")]`
one-element array (`CallPeakUnit.py:621`).

A second latent bug surfaced on the way: `at_end`, which reads the control value
at a position, returned `0.0` below the control track's first run end. Upstream's
`__chrom_pair_treat_ctrl` emits a triple carrying the *current* value of both
pointers, so a position before the control's first end still sees that track's
first value. Now falls back to the first run.

### Why the earlier "obvious fix" looked load-bearing

Returning a finite stand-in for the degenerate lambda changed 155 to 264 peaks --
which reads like "the degenerate branch is reached by the real corpus and its
faithful value matters". It was not: the same 264 is what the *correct* run now
produces, because upstream really does call 264 peaks on this corpus. The
substitute was silently inventing peaks on the no-control fixtures, and the
155/264 comparison was reading that invention as parity.

### Gate movement

```text
                     before          after
poisson/default       155/155        264/264    (288/288 fixtures full parity)
poisson/call_summits  166/166        275/275    (278 fixtures full parity)
poisson/broad         156/155        265/264    (287 fixtures full parity)
```

`--call-summits` and `--broad` still miss their upstream counts by a small
number of fixtures; `poisson/default` is now exact on every fixture, with
`worst_bnd 0`.

The `subtract` mechanism remains at 0 fixtures full parity (`worst_bnd 11475`) --
it is not upstream's mechanism (F46) and is retained only as a diagnostic
alternative, so this is a known non-goal rather than a regression.

**Verified:** 408 tests passed, 0 clippy, fmt clean.

## F128 -- `bdgopt` and `cmbreps` land, byte-identical; the bedGraph writer used the wrong float format

Two more commands are live in `macs3-rs`, both built on the `macs-bedgraph`
library from F124:

```text
bdgopt   multiply, add, max, min, p2q   -> BYTE-IDENTICAL to upstream (all 5)
cmbreps  max, mean, fisher              -> BYTE-IDENTICAL to upstream (all 3)
```

Verified by running upstream's own `bedGraphTrackI` on the same inputs and
diffing the files, not just by unit expectations.

### The bedGraph writer had the wrong float format

The bedGraph writer is **`%.5f`** (fixed five decimals), with a specific `track`
line (`BedGraphIO.py:97-143`):

```text
track type=bedGraph name="<NAME>" description="<DESC>" visibility=2 alwaysZero=on
```

My writer emitted a `%.6g` (the *peak* writers' format, which is genuinely
different) and a placeholder `track type="bedGraph"` line with no name or
description. Every one of the eight comparisons above would have failed on that
alone. Fixed, with name/description escaped as upstream does.

### `p2q`: an upstream quirk that is deliberately preserved

`bedGraphTrackI.p2q` (`BedGraph.py:904-984`) is a base-weighted AFDR walk, but
**the rank `k` is never incremented** inside the loop -- `nhcal` is, and is
unused. So `log10(k)` is always `0` and the mapping collapses to

```text
q = max(0, min(pre_q, p - log10(N)))   over distinct p, descending
```

This is *not* the same walk as the callpeak `pqtable` (which advances the rank by
length). Writing the textbook BH version would diverge from the oracle, so the
quirk is reproduced and the test pins upstream's bytes. `p2q` output for the test
input is `0.00000` then `0.19897`, which is what upstream writes.

### Two more derived-parser bugs, both found by using the commands

1. **`nargs *` / `nargs +` were silently dropping values.** The parser stored one
   string per destination, so `cmbreps -i a.bdg b.bdg c.bdg` kept only `c.bdg` --
   and `cmbreps` then failed with "needs at least two -i". Options now carry a
   `lists` map alongside `values`, with `get_all(dest)`, and the parser consumes
   the variable-length run. Pinned by `a_multi_value_flag_keeps_every_value`.
2. **A choice bug for these commands too.** `bdgopt -m p2q` and `cmbreps -m
   fisher` would have been rejected for the same comma/`|` reason F126 fixed; the
   shared `check_choice` now covers the multi-value path as well.

Neither was reachable before because no CLI command took a choice or a
multi-value flag end to end.

**Verified:** 371 tests pass (5 new), 0 clippy, fmt clean; the three peak gates
are unchanged (`155/155`, `166/166`, `156/155`, `worst_bnd 0`).

## F129 -- `--help` did not short-circuit required-flag validation

`macs3-rs callpeak --help` printed usage **and then** still failed with
`error: the following arguments are required: -t`, exiting non-zero. argparse
resolves `-h`/`--help` before any other validation, so `--help` must always
succeed and exit `0`. The derived parser validated required destinations after
the whole argument loop with no `help` escape. Fixed by returning early once
`--help` is seen, pinned by `help_short_circuits_required_validation`.

This one only surfaced because F128 wired real commands into the dispatch: until
then no command reached the `if o.help` branch with required destinations unset,
so nothing had exercised the combination. The accepted behaviour is unchanged
without `--help` (`macs3-rs callpeak -g 1000` still reports the missing `-t`).

Worth noting as a pattern: **three of the four bugs found this session
(F126 choice separator, F128 multi-value, F129 help) were all latent defects in
the derived argparse parser**, invisible until a command actually drove it
end to end. The parser was built from the matrix and unit-tested in isolation,
but nothing had ever run `macs3 <cmd> --flag value` through the whole chain.
Adding the three CLI end-to-end tests (`callpeak_cli.rs`, `bedgraph_cmds.rs`) is
what closed that gap; more commands will likely pay it off again.

## F130 -- `filterdup` lands, byte-identical on real data in all three `--keep-dup` modes

`macs3-rs filterdup` (`crates/macs-cli/src/commands/filterdup.rs`) reads BED,
limits tags per (position, strand) and writes upstream's `print_to_bed` layout.
Byte-identical to the pinned oracle on the `se_model/realistic` fixture across
every mode:

```text
--keep-dup 1    BYTE-IDENTICAL (192 records)
--keep-dup 2    BYTE-IDENTICAL (192 records)
--keep-dup auto BYTE-IDENTICAL (192 records)
```

and on a mixed-strand unit fixture that exercises the minus path.

### The strand convention is the subtle part

`FixWidthTrack.print_to_bed` (`FixWidthTrack.py:490-531`) writes

```text
+ reads:  chrom  pos      pos + fw  . . +
- reads:  chrom  pos - fw  pos      . . -
```

where `pos` is the *stored* coordinate -- and BED parsing keeps the **start** for
`+` but the **end** for `-` (`parse_bed_line`, mirroring `fw_parse_line`). So a
minus read at `[200, 210)` is stored at 210 and written back as `[200, 210)`,
reproducing its input interval, while a plus read at `[100, 110)` is stored at
100 and written as `[100, 110)` via `pos + fw`. Getting this backwards is
invisible on plus-only input and silently shifts every minus record by one
`tsize`, so the tests mix both strands deliberately.

Two further details reproduced:

* **Chromosome order.** Upstream iterates `get_chr_names()`, a Python `set`, so
  its multi-chromosome order is arbitrary (hash order, version-dependent). This
  port writes name order -- deterministic, and identical for the single- and
  small-chromosome cases the corpus uses. This is a deliberate, documented
  divergence from an upstream non-determinism rather than an oversight.
* **The `<= 1` bypass quirk.** With `--keep-dup auto` on 6 tags and
  `-g 2000000`, the binomial inverse `binomial_cdf_inv(1 - 1e-5, 6, 1/2e6)`
  yields `max_dup = 0`, which would drop everything -- yet both upstream and this
  port keep the lone `-` tag, because a strand holding at most one position
  bypasses the filter (and the total update). `macs-track`'s `filter_dup` already
  documents and reproduces this; the test pins it.

`cal_max_dup_tags` reuses `macs-stats`' `binomial_cdf_inv` directly, which is
already validated against upstream's binomial, so no new numerics were needed.

**Verified:** 376 tests pass (4 new), 0 clippy, fmt clean; all three peak gates
and the signal gate unchanged.

## F131 -- `randsample` lands byte-identical, confirming the NumPy RNG stream was already exact

`macs3-rs randsample` is live and byte-identical to the pinned oracle. The hard
part was never the command logic but the RNG: upstream seeds NumPy's **global**
stream and calls `np.random.shuffle`, so the surviving tags are reproducible only
if the shuffle matches exactly. `macs-stats`' `NumpyRng` and `macs-track`'s
`sample_percent` already reproduced it (F-era G4 work), and this command is what
finally *proves* it end to end -- previously nothing had driven that code path.

On the `se_model/realistic` fixture:

```text
-p 50 --seed 42   BYTE-IDENTICAL (95 records)
-p 30 --seed 42   BYTE-IDENTICAL (56 records)
-p 50 --seed 7    BYTE-IDENTICAL (95 records)
-p 30 --seed 7    BYTE-IDENTICAL (56 records)
-n 50 --seed 3    BYTE-IDENTICAL (49 records)
```

and on a 20-tag mixed-strand fixture the output is pinned to upstream's exact
bytes (`a_fixed_seed_reproduces_upstream_exactly`), including *which* tags the
shuffle kept -- a count match would not have caught a stream divergence, only
the byte comparison does.

Port notes (`randsample_cmd.py:28-78`):

* `-n/--number` overrides `-p/--percentage` (upstream recomputes the percentage
  from it) and `number > total` is rejected **before** any output file is
  created -- the acceptance criterion that errors precede output.
* `--seed >= 0` seeds NumPy globally; a negative value uses the ambient state.
  With no seed the result is intentionally non-reproducible, matching upstream.
* The `print_to_bed` layout is shared with `filterdup` (F130): `[pos, pos+fw)`
  for `+`, `[pos-fw, pos)` for `-`, with BED parsing keeping start for `+` and
  **end** for `-`.

The one thing deliberately *not* reproduced is upstream's per-chromosome `set`
iteration order for the write (see F130): this port emits chromosomes in name
order, which is deterministic and identical for the single-chromosome corpus.

**Verified:** 381 tests pass (5 new), 0 clippy, fmt clean.

## F132 -- `pileup` lands byte-identical; and a reminder that the parser seeds argparse defaults

`macs3-rs pileup` is live and byte-identical to the oracle across both extension
modes and two fragment sizes on the `se_model/realistic` fixture:

```text
directional  --extsize 100  BYTE-IDENTICAL (378 rows)
directional  --extsize 200  BYTE-IDENTICAL (376 rows)
-B           --extsize 100  BYTE-IDENTICAL (378 rows)
-B           --extsize 200  BYTE-IDENTICAL (380 rows)
```

No new numerics were needed: `SingleEndParams::directional` (`five=0`,
`three=d`) and `::bidirectional` (`five=d/2`, `three=d-d/2`) are precisely the
shift pairs `PileupV2.pileup_and_write_se` selects when `halfextension=False`,
which is the branch this command takes (`pileup_cmd.py:76-82`). `--both-direction`
doubles `--extsize` and switches to the symmetric pair, matching upstream. Values
are floored at the baseline *after* scaling and written `%.5f`.

### A test that failed for an instructive reason

The first version of the `pileup` test asserted that omitting `--extsize` is an
error. It isn't -- and the reason matters: **the derived parser seeds each
flag's argparse `default` before parsing**, so `--extsize` arrives as `200`
(upstream's documented default) even when the user never passes it. The test was
asserting the wrong thing about a correct implementation.

This is the mirror image of F126/F128/F129. Those were flags the parser *failed*
to honour; this is a flag it honours by supplying a default, which is correct
argparse behaviour but easy to forget when hand-writing a command. The test now
pins that the default (200) is used and agrees with an explicit `--extsize 200`.
The general lesson for the remaining commands: **a command that reads a flag must
ask what the matrix default is**, because silently inheriting 200 where upstream
would too is correct, and inheriting a wrong one is a silent numerical bug.

**Verified:** 384 tests pass (3 new), 0 clippy, fmt clean.

## F133 -- `bdgpeakcall` and `bdgbroadcall` land byte-identical; and a `store_false` parser bug

Two more commands are live, both ports of `bedGraphTrackI.call_peaks` /
`call_broadpeaks` (`BedGraph.py:405-668`), byte-identical to the oracle across
cutoffs, min-lengths, gaps, two chromosomes and `--no-trackline`:

```text
bdgpeakcall   -c 5 -l 100 -g 30              BYTE-IDENTICAL
bdgpeakcall   -c 3 -l  50 -g 10              BYTE-IDENTICAL
bdgpeakcall   -c 8 -l 200 -g 50              BYTE-IDENTICAL
bdgpeakcall   ... --no-trackline              BYTE-IDENTICAL
bdgbroadcall  -c 5 -C 2 -l 100 -g 30 -G 800  BYTE-IDENTICAL
bdgbroadcall  -c 4 -C 1 -l  50 -g 10 -G 100  BYTE-IDENTICAL
```

### A real parser bug: `store_false` was treated as `store_true`

The derived parser handled `store_true | store_false | count` identically,
inserting `flags[dest] = true` when the flag was seen. But `--no-trackline` is
`store_false` with default `true`: seeing it must set the destination to
**false**. Treating it as store_true made a negated flag indistinguishable from
the default, so `bdgpeakcall` emitted a track line even under `--no-trackline`.

This is the fourth distinct parser defect (after F126 choice separator, F128
multi-value nargs, F129 help short-circuit). All four share a root cause worth
naming: **the parser was built and unit-tested against the matrix in isolation,
but nothing ever ran a complete `macs3 <cmd> --flag` invocation through the whole
chain** until commands were actually wired in. Each command has now paid for
that gap once. The end-to-end tests (`callpeak_cli.rs`, `bedgraph_cmds.rs`,
`filterdup_cli.rs`, `randsample_cli.rs`, `pileup_cli.rs`, `bdgpeakcall_cli.rs`)
are the durable fix.

### Three output details the byte comparison forced

1. **Peak naming is not uniform across commands.** `callpeak` passes
   `name_prefix = "%s_peak_"` so its peaks are `NAME_peak_1`, but `bdgpeakcall`
   passes `oprefix + "_narrowPeak"` (no format template), so `peakprefix % name`
   raises and upstream falls back to the literal prefix -- giving `t_narrowPeak1`,
   *not* `t_narrowPeak_peak_1`. The writer now appends a bare counter and each
   caller supplies its own prefix.
2. **`name` and `name_prefix` are different strings.** The narrowPeak track line
   uses `name` (the oprefix, `t`) while peak rows use `name_prefix`
   (`t_narrowPeak`); the gappedPeak track line uses neither and defaults to
   `name="peak" description="peak"`. Conflating them produces a one-token
   difference on exactly one line.
3. **gappedPeak emits literal `0 0 0`, not a computed thickStart.** The 15-field
   row is `chrom start end name score . 0 0 0 blockNum blockSizes blockStarts
   fc pscore qscore` -- there is no thickStart field, so emitting the computed
   one added a spurious column. `blockNum`/`blockSizes`/`blockStarts` follow
   `__add_broadpeak` exactly, including the 1bp left/right complement blocks and
   the empty-set case (two 1bp blocks spanning the region).

**Verified:** 388 tests passed (4 new), 0 clippy, fmt clean.

## F134 -- `refinepeak` lands byte-identical (Watson/Crick tag depth), 9 of 14

`macs3-rs refinepeak` is live, byte-identical to the oracle across window sizes
and both cutoff regimes on the `se_model/realistic` fixture:

```text
-w 200 -c 5    BYTE-IDENTICAL     -w 50  -c 5    BYTE-IDENTICAL
-w 500 -c 5    BYTE-IDENTICAL     -w 200 -c 100  BYTE-IDENTICAL   (all _F)
-w 200 -c 1    BYTE-IDENTICAL
```

### It is not the `refine_peaks` I expected

`refinepeak_cmd` never calls `bedGraphTrackI.refine_peaks`. It reads tags, and for
each peak collects the tags in `[start - w, end + w]`, then evaluates a running
**Watson/Crick tag-depth** statistic (`find_summit`, `refinepeak_cmd.py:67-102`):

```text
wtd(j) = 2 * sqrt(watson_left(j) * crick_right(j)) - watson_right(j) - crick_left(j)
```

where `watson_*`/`crick_*` are the plus/minus tag counts in sliding half-windows
of size `w` around `j`. The argmax over `j` in the window is the refined summit,
suffixed `_R` when it beats `--cutoff` and `_F` otherwise. Reproduced exactly,
including the incremental slide (`left_forward` / `right_forward`) rather than a
per-position recount.

The names are `peakname_R` / `peakname_F` -- the peak's own name from the
`--bedfile` with the suffix appended, and the output is `\n`-**joined** with no
trailing newline (`b"\n".join([...])`), which the byte compare caught: my first
version emitted a trailing newline and `p1R` instead of `p1_R`.

The value can legitimately be **negative** (a 10-tag toy fixture yields `-1.00`),
so the `_R`/`_F` threshold is meaningful only against the actual statistic -- a
test that assumes "low cutoff always refines" is wrong, and mine was until I
checked the real value. The implementation was correct; the assumption was not.

**Verified:** 391 tests passed (3 new), 0 clippy, fmt clean.

## F135 -- `--help` is NOT a valid "is this command implemented?" probe

A direct consequence of the F129 fix, and worth stating plainly because it
quietly invalidated a check I had been using.

Before F129, `macs3-rs <cmd> --help` failed for any command with required
destinations, so its exit code happened to correlate with whether the command
was wired. After F129 `--help` short-circuits before validation -- which is the
correct argparse behaviour and is required for drop-in compatibility -- so it now
**exits 0 for all 14 subcommands, including the 5 that are not implemented**.

Running the inventory that way reported all 14 as live. The real check is either
to inspect the dispatch table:

```text
$ grep -oE '"[a-z]+" => macs_cli' crates/macs-cli/src/main.rs
callpeak bdgopt cmbreps filterdup randsample pileup bdgpeakcall bdgbroadcall refinepeak
```

or to make a real invocation and look for output. Current state:

| command | status |
|---|---|
| callpeak, bdgopt, cmbreps, filterdup, randsample, pileup, bdgpeakcall, bdgbroadcall, refinepeak | **live, byte-identical to oracle** |
| bdgcmp, bdgdiff, predictd, callvar, hmmratac | not implemented (`not yet implemented`, exit 1) |

The lesson generalises past this repo: **a smoke test that only checks "does the
CLI accept this invocation" can be silently voided by a correctness fix.** The
F129 change made `--help` more correct and simultaneously made my inventory
check wrong. Assertions about implemented surface should query the thing that
actually implements it.

**Verified:** 391 tests passed, 0 clippy, fmt clean; all three peak gates and the
signal gate unchanged; all 9 live commands produce real output files.

## F136 -- `predictd` / G8: PeakModel ported; d and alternative_d exact, `*_model.r` matches except a 1e-12 float tail

`macs-model` is no longer a placeholder. It ports `PeakModel.build` end to end
(`PeakModel.py:98-513`): naive pileup -> `naive_call_peaks` summits ->
`find_paired_peaks`/`find_pair_center` -> `model_add_line`/`count` strand
profiles -> normalise -> cross-correlate -> smooth -> local maxima -> `d`. And
`macs3-rs predictd` is live and dispatch-confirmed, closing **gate G8**.

On a synthetic corpus with a known separation the model recovers `d` and
`alternative_d` **exactly**, and the emitted R script matches upstream
line-for-line except the `ycorr` array:

```text
p (plus profile)     IDENTICAL
m (minus profile)    IDENTICAL
xcorr (lag axis)     IDENTICAL
altd (alternative d) IDENTICAL
d                    210  (exact; upstream also prints 210)
ycorr                max relative diff ~1e-12, same length, same argmax
```

### Three float-order bugs that byte-comparison caught

`model.r` is an acceptance criterion (`*_model.r`), so the last digits matter,
and each of these was invisible until the file was diffed:

1. **`np.correlate` is not `sum_j a[j+l]*b[j]` over `i+j=k`.** My first version
   paired indices by `i + j == k`, which mirrors the lag axis; the correlation
   peak landed at lag ~0 instead of the true lag and `d` came out 546 instead of
   210. The correct form is `out[k] = sum_j a[j + lag] * b[j]` with
   `lag = k - (len(b) - 1)`. Verified against `np.correlate([1,2,3],[10,20])`
   == `[20,50,80,30]`.
2. **`np.linspace` is not the integers `-(n/2)..=n/2`.** Upstream builds the lag
   axis as `np.linspace(len//2*-1, len//2, num=len)` -- exactly `len` *fractional*
   points. Using `(-half..=half)` produced one element too many and shifted every
   candidate `d`.
3. **`np.mean`/`np.std` use NumPy's pairwise summation, not a sequential `sum`.**
   Float addition is not associative, so the two disagree in the last bits, and
   that survives normalisation -> correlation -> smoothing. `numpy_pairwise_sum`
   (8-way unrolled, `PW_BLOCKSIZE = 128`) reproduces it exactly. Likewise
   `np.convolve` accumulates **forward** (`s[k+j]`), not reversed -- the flat
   window is symmetric so the sum is mathematically identical either way, but the
   digit sequence is not.

With those fixed, `p`, `m`, `xcorr` and `altd` are byte-exact and only `ycorr`
retains a ~1e-12 relative difference -- NumPy's correlate/convolve dispatch for
this particular length does not reduce in an order I could reproduce. It is
functional noise: same length, same argmax index (810), same `d`.

Also fixed: `find_pair_center`'s pairing rule is not a tag-count-weighted
average. Upstream requires `0.5 < pn/mn < 2` **and** `pp < mp`, and then takes
the plain integer midpoint `(pp + mp) / 2`. Getting that wrong was worth 336 bp of
`d` (546 vs 210).

**Verified:** 394 tests passed (3 new), 0 clippy, fmt clean; all three peak gates
and the signal gate unchanged.

## F137 -- `bdgcmp` / `bdgdiff`: `ScoreTrackII` and `TwoConditionScores` ported; 12 of 14 commands live

`bdgcmp` and `bdgdiff` are now implemented and byte-identical to upstream:
`macs3-rs bdgcmp` across all 8 methods x 2 pseudocounts x 12 random
2-chromosome fixtures (192 comparisons), and `macs3-rs bdgdiff` across 3 depth
scalings x 12 random fixtures (36 comparisons, all three BED files each). Both
are dispatch-confirmed, bringing live commands to **12 of 14**; only `callvar`
(fermi-lite assembler) and `hmmratac` remain.

`macs_score::ScoreTrack2` ports `ScoreTrackII` (`ScoreTrack.py:166-770`) and
`macs_score::TwoScores` ports `TwoConditionScores` (`ScoreTrack.py:1378-1960`).

### The float-width contract, pinned empirically

`bdgcmp` turned out to be a much sharper float probe than the peak caller,
because bedGraph pileups are *fractional* rather than integral, so every
pseudocount add and every division is observable. Four distinct width rules, each
pinned against upstream rather than reasoned about:

1. **`logLR_asym` / `logLR_sym` evaluate entirely in `f64`.** The parameters are
   `cython.float`, but in a C expression every `float` operand is promoted to
   `double`, so both logarithms, the products, the sum and the `log10(e)` scale
   all run in `f64` and the result is narrowed once on assignment. Computing in
   `f32` fails 24/76 rows; computing the logs in `f32` and the rest in `f64`
   still fails. With `treat=12.3, ctrl=1.9, pseudocount=0.1`, upstream scores
   `5.308994770050049` -- the `f32` path gives `5.308995723724365`, which
   prints `5.30900` instead of `5.30899`.
2. **`LOG10_E` is an `f64`, not the `f32` its annotation claims.**
   `LOG10_E: cython.float = 0.43429448190325176` sits at *module scope*, where
   Cython stores it as a Python attribute and does not enforce the annotation.
   Using the `f32` narrowing mismatches 21/76 `logLR` and 20/76 `slogLR` rows.
3. **The pseudocount is added in `f64` before the treatment count is truncated.**
   `19.899999618530273 + 0.10000000149011612` is `20.0` in `f32` but
   `19.999999620020389` in `f64`; upstream counts 19, not 20. This is fixed in
   `macs_score::pseudocounted_inputs`, shared with the peak caller. It changes
   no peak-caller result (those counts are integral) but fixes 3 of 59 rows on
   a fractional fixture.
4. **`FE` widens the division, `logFE` does not.** `compute_foldenrichment`
   evaluates `(p[i] + pseudocount)/(c[i] + pseudocount)` in `f64` and narrows
   once; `compute_logFE` divides in `f32` -- because the ratio is the argument to
   `get_logFE`, whose parameter is `cython.float` -- but takes the `log10` in
   `f64`. Both roundings are visible in the fifth decimal.

### Two upstream quirks reproduced deliberately

- **`make_pq_table` always zeroes the last score.** The back-fill loop starts at
  `i`, the loop counter, and `i` is left at `len - 1` when the loop completes
  without breaking. So the final distinct p-score is overwritten with 0 whether
  or not the `q <= 0` break fired. On the 76-row fixture, `p = 0.00639` keeps
  `q = 0.00156` while the immediately following `p = 0.00537` is forced to 0
  even though its computed `q` is a healthy `0.0053 > 0`.
- **`logLR_sym` is not symmetric.** Its `y > x` branch reuses `log(x) - log(y)`
  rather than swapping the arguments, so `sym(2, 10) = -3.5153444` while
  `asym(2, 10) = -2.0764158`. `asym` is not antisymmetric either: the two
  branches are different expressions, so swapping the arguments does not negate
  the result. Both pinned as `f32` bit patterns.

### Argparse mutually-exclusive groups (F138)

The flag matrix records each flag's attributes but not its *group*, so
`crates/macs-cli/src/flags.rs` now carries `MUTEX_GROUPS`: seven groups
transcribed from `bin/macs3` at the pinned commit. This fixed a real
accept/reject divergence -- `macs3-rs randsample -i x -f BED` previously
**accepted** a run upstream rejects with exit 2. Also pinned: argparse names the
later flag first in the conflict message (`-p 1 -q 1` reports `-q/--qvalue`), and
`bdgcmp`'s group lists `--o-prefix` before `-o/--ofile` while `add_output_group`
lists them the other way round.

`group_postprocessing` in callpeak looks like a mutual-exclusion group but holds
only `--call-summits`, so it constrains nothing; upstream accepts
`callpeak --broad --call-summits` together. Verified rather than assumed.

### Known divergence: `bdgcmp -m ppois -p 0` crashes

With the default `-p 0.0` and a control track containing explicit zeros,
upstream raises `AssertionError: Lambda must > 0` and Rust panics at the same
point. Same behaviour, but Rust's exit code is 101 rather than upstream's
traceback, so this is recorded rather than claimed as parity. It shares its root
cause with the open F127 peak-caller issue: a floor on the lambda side upstream
of `pscore`.

**Verified:** 407 tests passed (5 new bdgcmp/bdgdiff, 1 new parser test), 0
clippy, fmt clean; all three peak gates and the signal gate unchanged.

## F139 -- G0 hermeticity closed: the pinned oracle is no longer instrumented

The porting plan's stage 1 requires a *hermetic* oracle. An earlier phase
satisfied that only by option 1a -- instrumenting MACS3 in place
(`CallPeakUnit.py`, `PairedEndTrack.py`, `FixWidthTrack.py`) to dump the paired
position/treatment/control arrays, which are otherwise unreachable because
`__chrom_pair_treat_ctrl` is `cdef`-private (F60/F64a). That is a real deviation:
the working tree was not byte-identical to c5443190.

It is no longer necessary. `oracle/dump_stages.py` drives the *real* pipeline and
records what each stage computed -- the pipeline, the arithmetic and the ordering
are upstream's, only the observation is ours -- so it reaches everything the
instrumentation exposed without touching the oracle.

Restored the three sources from `/tmp/CallPeakUnit.py.orig` and
`/tmp/oracle_so_backup/CallPeakUnit...so`, then re-verified **all four gates are
unchanged**:

```text
default       155/155  bnd 155  summits 155  worst 0
call_summits  166/166  bnd 166  summits 160  worst 0
broad         156/155  bnd 155  summits 155  worst 0
signal         97 passed, 0 failed, 3 skipped
```

`CallPeakUnit.so`, `Pileup.so` and `PileupV2.so` are byte-identical to the
pre-instrumentation builds, so nothing was silently left instrumented.

New guard: `oracle/verify_oracle_clean.sh` checks the commit, `git status`, the
absence of `MACS3RS_` hooks in the Python sources, and the `.so` checksums.
`macs-stats`'s `oracle_tree_is_unmodified` test runs it, so CI fails if the
oracle is re-dirtied. Verified it actually detects a dirty tree (appended a line
to `BedGraph.py`, saw exit 1, restored, saw exit 0).

The three stale Cython `.c` artifacts still mention the hooks. They are
gitignored, never imported, and are reported as a note rather than a failure;
the loaded modules are the `.so` files, which are checksummed.

**Verified:** 408 tests passed, 0 clippy, fmt clean; four gates unchanged with
the oracle now hermetic.

## F140/F141/F142 -- `--call-summits`: the padded window must not inherit the `+1` XLS convention, and the SG derivative is decided by float cancellation

Three defects in the subpeak path, found by diffing the 288-fixture
`--call-summits` gate. Together they took it from 166/166 to **275/275
boundaries** (280 of 288 fixtures at full parity, `worst_bnd 0`).

### F140: the `+ 1` belongs to the reported start, not the padded window

`close_peak_wo_subpeaks` reports `region[0].start + 1` (F85), because the XLS
`start` column is written one base above `tstart`. `close_peak_with_subpeaks`
was reusing that value to *place its padded window*:

```python
peak_start = peak_content[0][0]      # CallPeakUnit.py:1438 -- no + 1
```

Adding one shifts `start_boundary` by one, so every smoothed maximum lands one
base late. On `sweep/gmini_mpe_d400_w600_ctrl_219` upstream's summit is 105 and
the `+1` form produced 106. The window now uses the raw first-chunk start and the
reported peak start still gets the `+ 1`.

### F141: `np.linalg.pinv` leaves a 2e-20 residue at the SG centre coefficient

`SignalProcessing.py:271-272` builds the derivative weights with
`np.linalg.pinv(b)[1]`, not the closed form `3k / (h*(h+1)*(2h+1))`. The two agree
to ~1e-19 relative but are not bit-identical: at `window_size = 180`, **153 of the
179 coefficients differ in the last bit**, and `pinv` leaves `2.18e-20` at the
centre where the closed form is exactly `0.0`.

That residue is load-bearing, because `maxima` calls `.round(16)` and the
derivative over a *flat-topped* pileup region is itself ~1e-16. Whether the flat
stretch reads as sign-0 is what decides one maximum versus two. The `pinv`
coefficients are therefore pinned for `window_size = 179`, the size
`--call-summits` actually uses (`smoothlen = min_length = d`,
`CallPeakUnit.py:1477`), and the closed form is kept for other sizes.

Also fixed here: the right-edge pad. Upstream is
`signal[-1] + abs(signal[-half-1:-1][::-1] - signal[-1])`, and that slice reversed
runs from `n-half-1` up to `n-2`. The port was walking `signal[n-2]` first, i.e.
`signal[::-1]`, which agrees only where the tail is flat.

### F142: the residual summit gap is float *cancellation order*, and the obvious fix is a trap

8 `--call-summits` fixtures still differ, all with the peak present and boundaries
exact but the summit off -- e.g. `sweep/gonechrom_mpe_d4000_w600_ctrl_077`, where
upstream reports summit 3075 and this port reported 3119.

Root cause, measured rather than guessed. Replaying upstream's own `maxima` on this
port's chunk list finds two candidate maxima, offsets 442 and 527; `enforce_peakyness`
discards **both** (their clipped regions have one unique value, failing the
`np.unique(s) >= 6` test), so upstream falls through to
`__close_peak_wo_subpeaks` and picks the median max-pileup chunk midpoint, 3075.
This port kept offset 527, so the two paths disagree about which summit exists.

The derivative over the flat stretch between 442 and 527 is a near-total
cancellation -- every term the same magnitude, result ~1e-16 -- so the *order of
the additions* decides the sign of the residue:

```text
sequential sum, closed-form coefficients -> +1.318390e-16  (does not round to 0)
upstream                                    -> -4.857226e-17 (rounds to 0)
8-way partial-sum accumulation             -> -3.469447e-17 (rounds to 0)
```

The 8-way partial sums fix this fixture and **regress the whole gate**, from 280
fixtures at full parity to 58 (`summits 12`). So numpy's convolution kernel is
neither sequential nor this partial-sum form; matching it means reproducing its
exact reduction, which is not derivable from the source. Reverted; the sequential
sum is kept and this gap is recorded rather than papered over.

Worth stating plainly: "the flat region reads as zero" is load-bearing for
`--call-summits` in a way that is invisible at `%.5f` output precision and shows
up only as an occasional one-maximum-versus-two summit decision.

**Verified:** 403 tests passed, 0 clippy, fmt clean. `poisson/default` 264/264
(288/288 fixtures), `poisson/broad` 265/264 (287 fixtures),
`poisson/call_summits` 275/275 boundaries with 267/275 summits.

## F143 -- `macs-hmmratac`: the HMMRATAC library is ported; hmmlearn's two passes use *opposite* matrix orientations

`crates/macs-hmmratac` is no longer a placeholder. It ports the pieces of
`hmmratac` that are deterministic:

| component | upstream |
|---|---|
| `pnorm2` | `Prob.py:62` |
| `generate_weight_mapping` | `HMMR_Signal_Processing.py:53` |
| `pileup_from_LR_hmmratac` | `PileupV2.py:680` |
| `make_bdg_of_bins_from_regions` | `HMMR_Signal_Processing.py:232` |
| `extract_value_hmmr` | `BedGraph.py:1074` |
| `extract_signals_from_regions` | `HMMR_Signal_Processing.py:146` |
| `generate_states_path` | `hmmratac_cmd.py:575` |
| `save_accessible_regions` | `hmmratac_cmd.py:621` |
| `hmm_model_save` / `hmm_model_init` | `HMMR_HMM.py:82`/`:118` |

Self-training is the one declared deviation (it is a `hmmlearn` Baum-Welch fit
seeded from `numpy.random.MT19937`, i.e. an SVD, a QR and a specific draw
sequence). **Inference is exact**, and `tests/hmm_inference.rs` pins it against
`hmmlearn.GaussianHMM.predict_proba` on a real 3-state, 4-feature model.

### hmmlearn 0.3.3's forward and backward disagree about the matrix

This is the finding. `predict_proba` calls `_score_log`, which runs
`_hmmc.forward_log` and `_hmmc.backward_log`. Both index `transmat`, but **not
the same way**:

```text
forward_log:  alpha[t][i] = logsumexp_j ( alpha[t-1][j] + log A[j][i] ) + logB[t][i]
backward_log: beta [t][j] = logsumexp_i ( log A[j][i]   + logB[t+1][i] + beta[t+1][i] )
```

The forward is the textbook form. The backward is its **transpose**: the state
leaving `j` is written as `A[j][i]`, not `A[i][j]`. Both were measured against
hmmlearn on the fixture; the transposed backward is what reproduces it:

```text
beta at t=3, hmmlearn   : [0.37969795, 2.8e-15, 0.62030205]
  textbook  A[i][j]     : [0.03015013, 3.0e-32, 0.96984987]
  transposed A[j][i]    : [0.37969795, 2.8e-15, 0.62030205]   <-- exact
```

Using the textbook form in both passes leaves the argmax right (so every
hmmratac output would still look correct) but moves the non-dominant states by
about 2%, which is enough to flip a tie in `generate_states_path` and shift a
region boundary.

The scaled variant (`_score_scaling`) is *not* what `predict_proba` uses, and it
raises `ValueError: forward pass failed with underflow` on these models -- a
trained ATAC HMM drives whole states to ~1e-88, which is well past what the
probability-space recursion survives. Worth knowing before reaching for it.

### Result

Worst total-variation distance over the 10-frame reference sequence:
**1.4e-6**, with the argmax identical on every frame. `generate_states_path`
labels each bin by its most probable state, so the argmax is the only thing that
reaches the output; the tolerance is set at 1e-5 with the argmax checked exactly.

### Two upstream details worth recording

- `save_accessible_regions` collects the whole `nuc`-`open`-`nuc` triple and then
  filters the `nuc` entries out, so a reported region is the **middle `open` run
  only** -- not the triple's span. The length test, however, is on the triple's
  span (`states_path[i+2][2] - states_path[i][1] > openregion_minlen`).
- The contiguity test is `states_path[i][2] == states_path[i+1][1]`, i.e. the
  *end* of one run must equal the *start* of the next.

**Verified:** 416 tests passed (12 new in `macs-hmmratac`), 0 clippy, fmt clean;
all four gates unchanged.

## F144 -- `macs3 hmmratac`: four independent bugs the cutoff-analysis report exposed

`hmmratac`'s `--cutoff-analysis-only` path exercises the fragment pileup, the
digest, the fragment-length EM and the fold-change normalisation in one shot, so
it is the cheapest end-to-end probe of the whole command. Getting it
byte-identical needed four separate fixes, each of which is load-bearing.

### 1. The FRAG count column was being dropped on the floor

`FragParser.build_petrack` always builds a `PETrackII`, whose rows carry a
multiplicity in a `u2` column. That column reaches everything downstream:

* `pileup_from_LRC_as_list` passes it as the **endpoint weight**, so depth is
  `sum(count)`, not the row count;
* `pileup_bdg_hmmr` explodes it with `np.repeat(arange(n), counts)` *before*
  applying the weight mapping -- so a row with count 7 adds the class weight
  seven times, which in `f32` is not the same as adding `7 * w` once;
* `PETrackII.finalize` sets `total = sum(counts)`, and `length` is
  count-weighted too, so `average_template_length` -- which hmmratac uses as
  `min_length` everywhere -- is `sum((end-start)*count) / sum(count)`.

Reading counts as 1 put `min_length` at 572 instead of 143 and made the whole
fold-change track ~2% high. The digest, the EM and the bin extraction were all
correct; only the multiplicity was wrong.

`--max-count` is a **cap**, not a collapse: `if max_count: count =
min(count, max_count)`. `PETrackII.add_loc` appends, and the per-barcode merge
lives in the `--barcodes` whitelist, which is applied at read time.

### 2. `bedGraphTrackI.cutoff_analysis` reports nothing when the first run is above the cutoff

```python
above_cutoff_startpos = pos_array[above_cutoff - 1]
```

When the leading run clears the cutoff, `above_cutoff[0]` is `0`, so the `-1`
wraps and `pos_array[-1]` is the **last** end in the track. The first chunk then
reads `(last_end, first_end)`, every subsequent chunk is appended (the gap test
`ts - lastp <= max_gap` never trips), and the whole track collapses into one
chunk whose measured length is negative. `peak_length >= min_length` therefore
never holds and *every* threshold reports zero peaks.

This is why the `--pileup-short` report is header-only: after normalising, the
short signal's smallest value is `5.9e-9`, `round(5.9e-9, 3)` is `0.0`, and the
leading run is above `0.0`. The default report is not affected because its
minimum is `0.17385` and `round(0.17385, 3)` is `0.174` -- just above the
leading run, so the leading run is excluded and `above_cutoff[0] > 0`.

Reproduced: `crates/macs-bedgraph/src/peakcall.rs`.

### 3. The ladder step is `f64`, but the span is `f32`

`minv`, `maxv` and `maxv - minv` are `cython.float`; only then is the difference
widened to a Python float for `/ steps`, and only `np.arange` widens again. Doing
the division in `f32` shifts the whole ladder by 0.01 -- which, on the
`mfrag_d400_w180` fixture, moves eleven of 48 reported cutoffs and makes two
`lpeaks` columns disagree.

### 4. `summary`'s standard deviation is computed over a broken `pre_p`

`BedGraph.py:374-399` resets `pre_p = 0` inside the first (per-chromosome) loop
but not in the second, so the variance pass measures every run against the *last
chromosome's* final `pre_p`. The mean is computed over correct run lengths and
the variance is not. hmmratac only reads `mean_v`, so nothing downstream changes,
but it is reproduced rather than corrected.

**Result:** `*_cutoff_analysis.tsv` is byte-identical to the oracle for both the
default and `--pileup-short` paths on `gmini_mfrag_d400_w180_ctrl_033`.

## F145 -- NumPy's `SeedSequence` is reproducible; its MT19937 array seeding is not (yet)

`PETrackII.sample_percent_copy` seeds with
`np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))`, and
`PETrackI.sample_percent_copy` with `np.random.default_rng(seed)` -- two
different generators from the same seed, so the two track classes cannot share a
sampler.

`SeedSequence` itself is now ported (`crates/macs-stats/src/seed_sequence.rs`) and
verified against the installed NumPy for the mixed pool *and* for
`generate_state(4)`. Two details had to be measured rather than read, and both
contradict the obvious transcription:

* **All four pool slots are hashed**, including the ones past the end of the
  entropy. Leaving them at zero gives a completely different pool. Checked for
  seeds 0, 1, 42 and 10151.
* **`generate_state` returns `hash(pool_word)`**, with no extra `data_val ^=`
  folded in. `hash` already XORs `hash_const` into the value; adding the pool
  word again gives a different state.

`MT19937(SeedSequence(...))` then feeds the generated key to `init_by_array`, and
**that half does not match NumPy**. The textbook `init_by_array` (key expansion
with `1812433253`, then the `1664525` and `1565833941` folds, then the top-bit
mask) reproduces `np.random.seed(1)` -- already covered by
`NumpyRng::seeded` -- but reproduces `np.random.seed([1])` and
`np.random.seed([1,2,3])` for no seed tried, under any of: the modular key
expansion, a zero-padded one, or either length of the second fold. NumPy's array
seeding must hash the key first. The previous `init_by_array` was worse than
wrong: it indexed `state[(i - 1).wrapping_sub(623)]` and panicked.

**Consequence.** The EM's 10% down-sample is a different permutation, so the
fitted `mono`/`di`/`tri` means differ by ~3% (155.7 vs upstream's 149.7 on the
fixture above) and therefore so does the weight mapping and the digest. This is
separate from, and additional to, the declared self-training deviation. It is
invisible with `--no-fragem`, where `means`/`stddevs` come from the user and are
used verbatim -- that path is exact.

**Verified:** 470 tests passed, 0 clippy, fmt clean. callpeak gates unchanged
(default 264/264 and 288/288 fixtures; `--call-summits` 275/275 boundaries,
267/275 summits, 280/288 fixtures; broad 265/264, 287/288 fixtures); signal gate
97 passed / 0 failed; oracle hermetic; hmmratac `*_cutoff_analysis.tsv`
byte-identical.

## F146 -- BAM: the BGZF block-size field is the member size *minus one*

Reading a BAM means reading BGZF, which is gzip with one extra gzip subfield.
`MACS3.IO.BAM` reads it as

```python
self.bamfile.seek(10, 1)                                    # past ID1 ID2 CM FLG MTIME XFL OS
xlen   = unpack("H", self.bamfile.read(2))[0]               # XLEN
extra  = self.bamfile.read(xlen)
bsize  = unpack("H", extra[extra.index(b'BC\x02\x00') + 4:])[0]
cdata  = self.bamfile.read(bsize - xlen - 19)
block  = decompress(cdata, -MAX_WBITS)                      # raw DEFLATE
self.bamfile.seek(8, 1)                                     # CRC32 + ISIZE
```

The arithmetic looks like it consumes the wrong number of bytes -- `bsize -
xlen - 19` rather than `bsize - xlen - 20` -- but it is exact, because the `BC`
field holds the **total member size less one**. So after the `seek(8, 1)` the
cursor lands on `bsize + 1`, which is the next block's start. Writing the field
as the member size instead puts every subsequent block one byte off, and the
reader's "bad gzip magic" check then fails on the second block.

`cdata` is *raw* DEFLATE (`-MAX_WBITS`), not zlib-wrapped, so `flate2`'s
`DeflateDecoder` is the right primitive and the header must be stripped first.

Two further traps, both of which produced a "valid" BAM that no reader could
use:

* **Records must not straddle a BGZF block.** Upstream's own comment notices
  this (`"I am not sure if i_bytes+entrylength can be outside of
  dcdata_length..."`) and its slice-then-unpack raises. samtools flushes at
  record boundaries for exactly this reason.
* **No alignment may start in block 0.** A BAI chunk with compressed offset 0 is
  indistinguishable from "no chunks at all": `get_coffset_by_region` returns 0
  for both and the caller returns an empty list immediately. The header therefore
  gets its own block, padded to a 4-byte multiple, as samtools writes it.

And `BAIFile.__load_bins` does `self.bins[i].pop(37450)` unconditionally, so a
BAI **without** the pseudo-bin raises `KeyError` on open -- even though the
pseudo-bin carries no offsets. `oracle/make_bam_fixtures.py` emits it.

## F147 -- BAM record filtering and `rightmost`

`__fw_binary_parse` (`BAM.py:757`) rejects in this order:

1. `flag & (4 | 512 | 256 | 2048)` -- unmapped, QC-fail, secondary, supplementary;
2. if `flag & 1` (paired), also require `flag & 2` (proper pair) and `!(flag & 8)`
   (mate mapped);
3. `MAPQ < 1 or MAPQ == 255`;
4. then a **hard failure** if the `MD` tag is absent.

`rightmost` sums only the reference-consuming CIGAR operations -- `M`, `D`, `N`,
`=`, `X` -- so `I`, `S`, `H` and `P` do not advance it. A soft-clipped read
therefore has a *shorter* span here than in a pileup: in the fixture,
`5S 20M 15S` at 1400 reports `rpos = 1420`, not 1440, while `10M 5N 10M`
reports 1425. Getting this wrong silently shifts every variant-calling window.

`__get_SEQ_QUAL` trims a trailing `=` -- the padding nibble for an odd `l_seq` --
and then asserts `len(seq) == len(qual)`, so the decoded sequence must be exactly
`l_seq` bases.

`get_n_edits` counts CIGAR `I` and `S` lengths plus every `MD` byte in the raw
range `'A'..='Z'`. That range includes `N`, so a `^NN` deletion contributes **two**
edits. `RACollection.remove_outliers` filters on it, so the count is
load-bearing.

`get_reads_in_region` overlaps with `rpos > left` and stops at the first record
with `lpos > right`, and suppresses duplicates by comparing
`(lpos, rpos, strand, cigar)` on **consecutive** records, keeping the first
`maxDuplicate`.

**Verified:** `crates/macs-io/src/bam.rs` is byte-identical to
`BAMaccessor.get_reads_in_region` across 8 region queries covering both
contigs, whole-contig scans and `--keep-dup 1`/`2`/`4` (1218 records compared) --
`oracle/run_bam_e2e.sh`, plus 11 unit/golden tests in `crates/macs-io/tests/bam_parity.rs`.
492 workspace tests, 0 clippy, fmt clean.

## F149 -- the callpeak gate measures coordinates only, hiding a paired-end p-score residual

Wiring the `callpeak` CLI to paired-end input (`-f BEDPE`/`BAMPE`/`FRAG`) made it
possible to compare the **whole** `*_peaks.xls` against the oracle rather than
just the peak coordinates. Doing so over all 133 `*_mpe*` sweep fixtures:

| result | count |
|---|---|
| `*_peaks.xls` byte-identical, including every score | 110 |
| coordinates and summits identical, scores differ | 23 |

Worst observed deviations:

```text
max |d -log10(pvalue)|   = 6.690e-01   gonechrom_mpe_d400_w600_ctrl_383 chr1:2707
max |d fold_enrichment|  = 2.357e-01   same
max |d -log10(qvalue)|   = 4.920e-01   same
```

so the residual is far outside the `<= 1e-9` acceptance criterion. The other 22
are much smaller (4th-5th decimal, e.g. `1.67134` vs `1.6713`).

**This is a pre-existing library gap, not a CLI regression.** Verified by running
the `callpeak-e2e` harness on the same fixture:

```text
CALLPEAK_WRITE_XLS=... macs-callpeak-e2e gonechrom_mpe_d400_w600_ctrl_383 ...
  pscore = 111.283      # identical to the CLI
upstream               = 111.952
```

`oracle/run_peak_e2e.sh` reports `264/264` for this fixture because its
`read_xls` only parses `chrom start end length abs_summit` and compares those.
The score columns were never checked for *any* mode, so this survived.

### Neither p-score construction reproduces it

Two candidates, both implemented:

* `retadd_pscore` -- upstream's `retadd` union walk over `min(p1, p2)`, with
  the treatment truncated to an integer and **no** pseudocount (F51/F55). This
  finds the peak but gives `111.283`.
* `paired_pscore` (F59) -- boundaries only where *both* tracks change, which is
  upstream's `pos_array[above_cutoff]` over paired arrays. This **loses the peak
  entirely** on `gonechrom_mpe_d400_w600_ctrl_383` (and on
  `gmini_mpe_d1200_w600_ctrl_066` when driven through the CLI).

`paired_boundaries` therefore defaults to `false` in `PeResult`, matching the
gate's calibration, and F59's variant stays reachable behind the harness's
`CALLPEAK_BOUNDARY=paired` switch. Since neither matches, the truth is likely a
mix -- union boundaries (so no peak is lost) carrying paired *values* -- which is
the next thing to measure. Recorded as an open gap, not closed.

### Also fixed while wiring the CLI

* **`--extsize` is ignored in PE mode.** `callpeak_cmd.py:165` sets
  `options.d = options.tsize` under `--nomodel` for `PE_MODE`, with no
  `info("#2 Use %d as fragment length")` line. `-s/--tsize` is the only override.
  Using `--extsize` here silently changed `d`, `min_length` and `max_gap`.
* **`-g` needs the species shortcut table.** `EFFECTIVEGS` is `hs`/`mm`/`ce`/`dm`
  only -- and `hs` is GRCh38 (2 913 022 398), *not* the `hg38` row of
  `genomesize::TABLE`, so the two lists must not be merged.
* **Python's `"%.2e"` pads the exponent**: `4.82e+04`, not Rust's `4.82e4`.
* **`# control file = None`**, not `[]`: `-c` is `nargs="*"` with no default, so
  an absent control interpolates as the object `None`.
* **Zero peaks is not an error.** Upstream exits 0 and writes a header-only xls.
* The whole `# ARGUMENTS LIST:` block is conditional per flag, and its lines are
  byte-compared; `mfold` interpolates as a Python list repr (`[5, 50]`).

**Verified:** `*_peaks.xls` now byte-identical (including all scores) on 110 of
133 PE fixtures; `*_peaks.narrowPeak` and `*_summits.bed` byte-identical on the
ones checked. All three callpeak gates unchanged (default 264/264 + 288/288
fixtures, `--call-summits` 275/275 + 280/288, broad 265/264 + 287/288).

## F150 — the control lambda is written from the *paired union walk*, not from the pileups

`callpeak -B` writes `<name>_treat_pileup.bdg` and `<name>_control_lambda.bdg`.
Both are **two coalesced views of the same `chr_pos_treat_ctrl` rows**
(`__write_bedGraph_for_a_chromosome`, `CallPeakUnit.py:1717`), not
serialisations of the two pileups. Four byte-visible details:

* **no `track` line** — bare `fprintf`s to a file handle, unlike
  `bedGraphIO.write_bedGraph` which emits one;
* the run start is a running `pre` initialised to `0`, not the track's own start,
  so the first row spans `[0, pos[0])`;
* a new row only when the value changed by **more than `1e-5`** ("precision is 5
  digits") — which is what turns ~2800 union rows into ~336;
* values are `"%.5f"` with trailing zeros **kept** (`2.2555000782` →
  `2.25550`), and the final row is emitted unconditionally.

`--SPMR` divides by the deeper sample's million-normalised total; otherwise the
denominator is `1.0`.

Serialising the raw pileups instead produced 2821 runs against upstream's 336 for
`gonechrom_mpe_d400_w600_ctrl_383`, with identical values but different
breakpoints. The 4-column form (`chrom`, `start`, `end`, `value`) is required;
a 3-column bedGraph is not interchangeable.

### The PE control projection is symmetric, and the source reads otherwise

`pileup_treat_ctrl_a_chromosome` (`CallPeakUnit.py:606`) calls
`pileup_a_chromosome_c(chrom, ctrl_d_s, ctrl_scaling_factor_s,
baseline_value=lambda_bg)` for `PE_MODE` and **omits** `directional`, whose
declared default on `FixWidthTrack` is `True` — which would give
`five_shift = 0`, `three_shift = d`.

It is nevertheless symmetric, because `PE_MODE` dispatches to
`PETrackI.pileup_a_chromosome_c` (`PairedEndTrack.py:829`), which has **no
`directional` parameter at all** and hardcodes `five_shift = three_shift = d // 2`.

Switching the PE control to `directional` makes
`gonechrom_mpe_d400_w600_ctrl_383` worse (`-log10(pvalue)` 100.756 against
upstream's 111.952, versus 111.283 for the symmetric projection), so the
symmetric projection stands — but the F101 conclusion should cite
`PairedEndTrack.py:829`, not `FixWidthTrack`, because the two classes disagree.

## F151 — `--format FRAG` is a `PETrackII`, and it projects differently again

`--format FRAG` gives a counted `PETrackII`, whose two pileup entry points are
**not** `PETrackI`'s:

* treatment: `pileup_from_LRC_as_list` — `[l, r)` coverage weighted by the count;
* control: `pileup_a_chromosome_c` → `pileup_from_LRC_centers_as_list`
  (`PairedEndTrack.py:1458`) — one **`d`-wide window centred on each fragment
  *end***, `[x - d//2, x - d//2 + d)`, carrying the count.

Neither is the unit-depth `d//2`-shift projection used for BEDPE/BAMPE. Three
consequences, all of which had to be fixed together:

* the treatment needed `pileup_from_weighted_fragments` — the unit sweep made
  every FRAG fixture's depth ~`mean(count)` times too small and produced **no
  peaks at all**;
* the control needed `pileup_from_weighted_positions`;
* **pooling dropped the counts.** `pool_pe` re-pushed every fragment through the
  unit `push`, so a single `-t` — which goes through the pooling path — lost all
  multiplicities. That is why the summit pileup read 947 where upstream reports
  3783.

### FRAG also forces duplicate filtering off

`OptValidator.validate_options` (`MACS3/Utilities/OptValidator.py:110-112`):

```python
if options.format == 'FRAG' and options.keepduplicates != "all":
    info("Since the FRAG format specified, all duplicates in reads will be kept.")
    options.keepduplicates = "all"
```

FRAG rows already carry a multiplicity, so deduplicating them discards the
weights — and, because the xls count block is gated on
`keepduplicates != "all"`, it also emits four "after filtering" lines upstream
never writes for FRAG. `# total fragments in treatment:` is the **sum of the
count column** (4794 on `gmini_mfrag_d1200_w180_ctrl_042`, not its 1200 rows).

### FRAG peak starts can be negative

`pileup_from_LRC_centers_as_list` does **not** clip: with `llocal = 10000`,
`x - llocal/2` is `-5000` for a fragment starting at 0, and the union walk
carries those positions. `gmini_mfrag_d1200_w600_ctrl_069` therefore has a peak
whose xls `start` is `-13` (`length 194`, narrowPeak start `-14`). This is still
open — see Next.

## F152 — no control **file** is not `--nolambda`

Two different code paths were being collapsed onto one:

* `--nolambda` empties `ctrl_d_s`/`ctrl_scaling_factor_s`, which sets
  `no_lambda_flag`, and the control pileup becomes the **one-element**
  `[treat_pv[0][-1:], [lambda_bg]]` array;
* with no control *file*, `PeakDetect.call_peaks` (`PeakDetect.py:109`) dispatches
  to `__call_peaks_wo_control` (`PeakDetect.py:1104`), which builds a
  **single-scale control from the treatment itself**:

  ```python
  ctrl_scale_s      = [float(self.d) / self.lregion,]
  ctrl_d_s          = [self.lregion,]
  ```

  — one `llocal`-wide window scaled by `d/llocal`.

Collapsing both onto the one-element array left a single paired index, so the
q-value could never clear `-log10(0.05)` and **no peak was ever called on a
no-control fixture**: all 100+ `sweep/*_noc_*` plus `se_edge/no_control`.

Note the second case also depends on the single scale landing in
`LambdaScales::d`: `scales == vec![llocal]` puts `llocal` in the `d` slot, so
`lscales.llocal == 0` is *not* a valid "no lambda" test.

## F153 — `# tag size is determined as N bps` is the mean of the **first ten** tags

`GenericParser.get_tag_size` (`Parser.py:275-293`) does not average the file:

```python
while n < 10 and m < 10000:
    m += 1
    thisline = self.fhd.readline()
    this_taglength = self.tlen_parse_line(thisline)
    if this_taglength > 0:
        s += this_taglength
        n += 1
self.fhd.seek(0)
self.skip_first_commentlines()
if n != 0:
    self.tag_size = cython.cast(cython.int, (s/n))
```

so `options.tsize` — and the xls header line, `%d`, so truncated — is the
truncated mean of the **first ten** tags with a positive length.
`gmini_mse_d1200_w60_ctrl_009` reports 190 while its whole-file mean is 151.

Two further points, both of which had to be right simultaneously:

* this is the **measured** size. `--extsize` still overrides `d` for the pileup
  extension afterwards, which is why a 190 bp header pairs with a 200 bp
  extension. Reporting `--extsize` in the header is wrong; so is using the
  measured mean as `d`.
* `tlen_parse_line` is `atoi(f[2]) - atoi(f[1])`, computed **before** the
  strand swap, so a minus-strand line's length is `start - end` — the absolute
  difference of the parsed 5' coordinate and the end column.

## F154 — the summit p/q come from the chunk's paired values, and the control is indexed by *cursor*

Two independent bugs in `callpeak`'s single-end path, both invisible to a
coordinate diff and both fixed by the golden byte gate.

### F154a — the summit p/q are recomputed, not read off the score track

`__close_peak_wo_subpeaks` (`CallPeakUnit.py`):

```python
summit_treat = peak_content[summit_index][2]
summit_ctrl  = peak_content[summit_index][3]
summit_p_score = get_pscore(cython.cast(cython.int, summit_treat), summit_ctrl)
summit_q_score = self.pqtable[summit_p_score]
```

The score **track** is a different array with a different indexing (it is built
over score spans; `peak_content[i][4]` is `ti`, an index into the paired
arrays). F120 had the summit read `pscore_track[ti]`. Where the two indexings
disagree the recomputed score can land **outside the whole histogram**, so
`qscore_or_zero` returns its documented zero fallback and a literal `0` is
written into the XLS, the narrowPeak `qscore` column and the summits `score`
column.

This was the single largest defect in the port: `qscore 0` against upstream's
`31.4677` on `onechrom/single_contig` and `16.685` on `tiny/two_contigs`, i.e.
most of the corpus. `close_peak_wo_subpeaks` now recomputes from
`(int(summit_treat), summit_ctrl)`, consulting the track only when no cache is
available (the library's own unit tests pass `cache: None`).

### F154b — `chr_pos_treat_ctrl` pairs by **cursor**, not by "last end at or before"

The single-end path had built its paired arrays as a sorted, deduplicated union
of run ends, sampling the control as *the last run ending at or before `p`*.
`__chrom_pair_treat_ctrl` (`CallPeakUnit.py:519`) instead walks head-to-head and
emits `(min(t_p[i], c_p[j]), t_v[i], c_v[j])` -- the values each track's cursor
currently carries, so a position contributed by the treatment still shows the
control run *ahead* of it, not the one behind it -- and stops as soon as either
input is exhausted (`while it < lt and ic < lc`), dropping the tail.

On `tiny/two_contigs` the summit reported ctrl `7.4` (fold `5.35269`) where the
union value at that position is `7.6` (fold `5.22821`, upstream's answer). Two
sub-rules had to be reproduced together, since either alone still misaligns the
paired index space: the cursor semantics above, and the shared
`above_cutoff_startpos = pos_array[above_cutoff - 1]` with the
`if above_cutoff[0] == 0` fix-up.

`__cal_qscore(treat_array, ctrl_array)` also means the q-value is looked up
**per paired index** from the shared table, not resampled from a separate
q-score track.

## F155 — the weighted pileup coalesces, and `end = start + d`

`_pileup_sorted_weighted_as_list` (`PileupV2.py`) ends each step with

```python
if scaled_z == pre_z:
    ret_p_ptr[c-1] = p          # extend the previous run's end
else:
    ret_p_ptr[c] = p; ret_v_ptr[c] = scaled_z; c += 1; pre_z = scaled_z
```

so an equal-valued neighbour **overwrites the previous run's end** rather than
adding a breakpoint. Unlike the unweighted sweep (F71, where the rule is
deliberately *off*) this is unconditional: a counted event's `+w`/`-w` pair can
land on one position and leave the depth unchanged. `--format FRAG`'s control
lambda therefore has far fewer breakpoints than its raw event list, and
`SignalTrack::push`'s coalescing is the right behaviour, not a bug.

`pileup_from_LRC_centers_as_list` also builds `end = start + d` from
`start = x - d//2`, where the unweighted positional path builds `x + d//2`.
Those agree only when `d` is even.

## F156 — two different template lengths in `PeakDetect.__init__`

There are two in scope:

```python
if self.PE_MODE:
    d = self.treat.average_template_length      # the post-filter mean
    control_total = self.control.total * 2
    treat_sum   = self.treat.length
    control_sum = control_total * d             # <-- the local mean
    self.ratio_treat2control = float(treat_sum)/control_sum
...
ctrl_d_s = [self.d]                            # <-- options.tsize
ctrl_scale_s = [ratio]                                  if not tocontrol
ctrl_scale_s += [float(self.d)/sregion * ratio]         if sregion
ctrl_scale_s += [float(self.d)/lregion * ratio]         if lregion
```

Only `control_sum` uses the local mean; every scale factor uses `self.d`, the
truncated as-read `tsize`. F95 had both using the post-filter mean, which is the
same number only when duplicate filtering removes nothing.

## F157 — `--mfold` is `nargs=2`, so `--mfold 3 20` is two argv tokens

The flag matrix records `nargs` as the string `"2"`, and the parser only handled
`*`/`+` (variable) and the implicit single value. `-m 3 20` consumed `3` as the
value and then rejected `20` as an unrecognised argument, so **every case in the
`mfold_3_20` matrix variant exited 2** instead of 0.

A fixed `nargs=N` is now consumed as `N` tokens, with the whole list kept so
`# model fold = [3, 20]` can be rendered. The subtlety that made the obvious
version silently wrong: the index must be advanced **past** the last value, or
the outer loop re-reads it as a flag.

## F158 — the golden byte gate, and what it found that coordinates hid

`oracle/run_peak_e2e.sh` compares peak **coordinates and summits** against
upstream's XLS. That is a real gate but a weak one for this command: a score
that moves in the fourth decimal, a header count that is wrong, or a `-B`
bedGraph that has the right values at the wrong breakpoints all pass it.

`oracle/run_golden.py` (+ `run_golden.sh`) instead replays **every** recorded
argv in `tests/golden/<group>/<fixture>/<variant>/command.json` with `macs3-rs`
and byte-compares every file the recording lists. Three header lines are
normalised away, because they quote the invocation and cannot match a replay:

* `# Command line:` (different interpreter, and `--outdir` moves);
* `# ChIP-seq file = [...]` / `# control file = [...]` in the `# ARGUMENTS LIST`
  block, which echo whatever spelling of the path the caller used.

Every other line -- all counts, `d`, scale factors, and every peak row -- is
compared verbatim, so nothing can hide behind the normalisation.

Measured over the whole corpus: `default` **267/425** cases byte-identical (181
before this round), all recorded variants **2906/7685**. The failures cluster by
matrix variant, which is exactly the shape the acceptance criteria ask about:

| variant | failing cases |
|---|---|
| `B`, `broad`, `mfold_3_20` | 425 each |
| `nomodel_shift`, `shift_only` | 274 each |
| `scale_to_large` | 272 |
| `keepdup_all`, `keepdup_auto` | 262, 259 |
| `nolambda` | 256 |
| `call_summits`, `q05`, `bw300`, `default`, `gsize_numeric`, `nomodel_extsize`, `spmr` | 158-169 each |
| `q001` | 147 |
| `slocal_500_llocal_2000` | 106 |
| `keepdup1` | 103 |

So the remaining work is dominated by **flag-matrix** coverage rather than by the
`default` numerics: broad mode, `-B`'s leading row, `--shift`, `--scaleto large`,
`--keep-dup`, and `--nolambda` each fail on most of the corpus, and the FRAG
family (132 of the 158 `default` failures) still differs in the weighted control
lambda's 5th significant figure.

## F159 — broad mode: the two missing writers, and an empty p→q table

`--broad` never wrote `<name>_peaks.broadPeak` or `<name>_peaks.gappedPeak` at
all, so **every case in the `broad` matrix variant failed** (425/425) however
correct the peak rows were.

### The writers

`PeakIO.write_to_broadPeak` (`PeakIO.py:1470`) -- nine columns, **no track line**
(`trackline=options.trackline`, off by default):

```c
"%s\t%d\t%d\t%s%d\t%d\t.\t%.6g\t%.6g\t%.6g\n"
```

The **column names are misleading**: field 6 is `fc`, field 7 `pscore`, field 8
`qscore`, and `score` is `int(10 * qscore)`. So the golden
`905 . 21.2389 91.1085 90.5768` has `score = int(10 * 90.5768)`. Rows are
numbered per **group of equal `end`**, taking the group leader -- unlike the
narrowPeak/XLS writers, which number every row. Both files use the 0-based
half-open convention, so the XLS `start` (already 1-based, F85) shifts down by
one.

`PeakIO.write_to_gappedPeak` (`PeakIO.py:1343`) -- BED12 + 3, again no track line,
gated on `peak["thickStart"] != b"."` so a peak with no block structure is
**skipped entirely**. It hardcodes `0 0 0` for thickStart/thickEnd/itemRgb; only
the block columns carry information.

The block structure is `__add_broadpeak` (`BedGraph.py:600`): one block per
strong (lvl1) sub-peak, positioned relative to the broad region, plus a **1 bp**
block on each side when the strong set does not reach the broad region's own
ends. With no strong sub-peaks at all it still emits **two** blocks (`1,1` at `0`
and `end-start-1`). `gmini_mse_d1200_w180_noc_117` gives one lvl1 peak of length
345 at start 4 inside `[3,353)`, so both flanks are added and the result is
`3` blocks, sizes `1,345,1`, starts `0,1,349`.

To get the lvl1 set to the writers, `call_chromosome` now returns
`macs_peaks::callpeak::Called { peak, lvl1 }` instead of a bare `Peak`.

### The empty table

`broad_level` was calling `close_peak_for_broad_region(r, params, scores,
&PqTable::empty(), None)`. That closer derives the peak's q-score with
`__cal_qscore(tarray_pileup, tarray_control)` -- a **per-position lookup into
`self.pqtable`** (`CallPeakUnit.py`, `__close_peak_for_broad_region`) -- so an
empty table made every lookup take its zero fallback and every broad peak was
written with `-log10(qvalue) = 0` and `score = int(10 * 0)` on the whole corpus.

The real table and the p-score cache are now threaded through. The rest of the
closer was already right: `peak_score`, `pileup`, `pscore` and `fc` are all
length-weighted means over the region's chunks, via `mean_from_value_length`
(`sum(v_i * l_i) / sum(l_i)`, accumulated in `f64` and cast to `f32` at the end),
and the summit is `0`.

## F160 — `--shift` is single-end only, and the single-end pileup takes it

`pileup_treat_ctrl_a_chromosome` (`CallPeakUnit.py:588-620`) passes `end_shift`
to `self.treat.pileup_a_chromosome(..., directional=True, end_shift=self.end_shift)`
in the **non-PE** branch, and to **neither** pileup in the PE branch --
`PETrackI.pileup_a_chromosome` and `PETrackI.pileup_a_chromosome_c` take no shift
argument at all. So `--shift N` moves every 5' end by `N` before the
single-end extension and is a **no-op** in paired-end mode, whatever the help
text's "for BAM/BEDPE input" says.

`--shift -100` on `gmini_mse_d1200_w180_noc_117` moves the reported peak from
`5-349` to `1-253`; with the flag ignored it stayed at `5-349`.

## F161 — `--keep-dup auto` is `binomial_cdf_inv(1e-5 complement, N, 1/gsize)`

`cal_max_dup_tags` (`callpeak_cmd.py:331`) is

```python
return binomial_cdf_inv(1-p, tags_number, 1.0/genome_size)
```

i.e. the smallest `x` whose binomial CDF reaches `1 - p` for
`Binomial(tags_number, 1/gsize)`. It needs the **retained treatment count**, so
it can only be resolved after loading; the control is then filtered with the
*same* number (`callpeak_cmd.py:141` uses `treatment_max_dup_tags`, not a
control-specific one). `macs-stats::binomial` already had `binomial_pdf` and
`binomial_cdf_inv`; only the wiring was missing, and the `keepduplicates` flag was
returning a usage error -- so every case in the `keepdup_auto` variant exited 1.

## F162 — `--nolambda` is decided by the **scale lists**, and applies with a control file present

`PeakDetect.call_peaks` empties `ctrl_d_s`/`ctrl_scale_s` when `nolambda`, and
`CallerFromAlignments.__init__` turns that into `no_lambda_flag`. The flag is then
tested **before** the control is looked at:

```python
if not self.no_lambda_flag:
    ctrl_pv = self.ctrl.pileup_a_chromosome_c(chrom, ctrl_d_s, ctrl_scaling_factor_s, baseline_value=self.lambda_bg)
else:
    # a: set global lambda
    ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg,], dtype="f4")]
```

Two things follow, and both were wrong in the port:

* it is keyed on the **empty scale lists**, not on the absence of a control file,
  so it applies with `-c` present too. Checking `--nolambda` only in the
  no-control arm meant it was silently ignored whenever `-c` was supplied and
  the dynamic lambda was used instead;
* the substitute is a **one-entry** array, so the union walk's cursor rule pins
  `ic` at 0 until the treatment is exhausted and every paired row carries
  `lambda_bg`. Building it as a run-length track with
  `SignalTrack::empty(chrom, t.end(), t.end())` + `push(t.end(), lambda_bg)`
  produces an **empty** run -- `push` drops a run that does not advance the
  cursor -- so the control came out with no runs and nothing to pair against.
  The faithful equivalent is a single run `[0, t.end())` carrying `lambda_bg`,
  which is what it now builds.

Verified against `frag_basic/barcode_fragments --nolambda`, whose golden control
bedGraph is one row `[0, 200) = 16.32253` (i.e. exactly `lambda_bg`) over a
94-row treatment track.

## F163 — `# total %ss in control` is not gated on `--keep-dup`

`callpeak_cmd.py:126-155`:

```python
if control is not None:
    tagsinfo += "# total %ss in control: %d\n" % (tag, c0)
    if options.keepduplicates != "all":
        tagsinfo += "# %ss after filtering in control: %d\n" ...
```

Only the three lines *after* the total are gated. Tying the whole control block to
"was a duplicate limit computed" dropped the total line for every `--keep-dup all`
run -- which is exactly where `--format FRAG` forces `all` (F151), so **every FRAG
fixture** lost the line.

## F164 — `--scaleto large`: the source reading does not reproduce the oracle

`callpeak_cmd.py:242-259` reads

```python
if options.scaleto == "large":
    options.tocontrol = (t1 <= c1)      # post-filter, c1 already doubled in PE
else:
    options.tocontrol = (t1 > c1)
```

That makes `tocontrol` true on `frag_basic/barcode_fragments` (`t1 == c1 == 7995`)
and `treat_scale = 1/ratio = 2`, doubling every reported pileup against the
oracle's unscaled value.

**Superseded -- the rule as written is correct and the earlier reading was not.**
`t1`/`c1` are the totals *after* `filter_dup` (`callpeak_cmd.py:123` and `:146`),
or `t0`/`c0` when `keepduplicates == "all"`, and `c1` is then doubled in PE mode at
`:215`. For `frag_basic/barcode_fragments`, `--format FRAG` forces
`keepduplicates = "all"` (`OptValidator.py:110`), so `t1 == c0 == 7995` and
`c1 == 15990`: `--scale-to large` takes the `t1 <= c1` branch, `to_control` is
true and `treat_scale = 1/ratio = 2` -- which is exactly what the oracle reports
(summit pileup 15978 instead of 7989). Implemented per source, both the default
and `--scale-to large` variants are byte-identical and the whole family moved with
them. The error was in reading `t1 > c1` as "scale control to treatment" with
*undoubled, pre-filter* counts.

## F165 — event weights sort with their positions, not independently

`pileup_from_LRC_centers_as_list` sorts `start_poss` with `np.argsort` and
permutes `start_weights` with the *same* argsort. A port that sorts the positions
and leaves the weights in input order attaches each count to the wrong position.
The totals stay plausible -- most counts are similar -- so only the summit pileup
and the `-log10(p)` column move.

## F166 — `tsize_exact` is the unweighted row mean

`PeakDetect.__init__` takes `d=options.d`, while `__call_peaks_w_control` reaches
`float(self.d)` in the wide-window control factors. `options.tsize` is
`tp.d`, the parser's **pooled** mean (count-weighted for `--format FRAG`), while
`self.d` is the **unweighted** mean over rows. The two differ as soon as a row
carries a multiplicity, so the wide-window control scale needs both numbers.

## F167 — a counted control's projection is centred and weighted

`PETrackII.pileup_a_chromosome_c` calls
`pileup_from_LRC_centers_as_list`, which makes one `d`-wide window per fragment
*end* at `x - d//2` carrying the count. `PETrackI` instead uses
`pileup_from_PN_shifted`: unit-depth ends shifted `d//2` both ways. The two classes
need different projections, not one shared one.

## F168 — `over_two_pv_array` drops the longer array's tail

Measured against the compiled extension: the walk stops when either input is
exhausted, so the extra tail of the longer array is discarded rather than
flushed. An implementation that flushes the tail produces a longer merged lambda
and extra union rows.

## F169 — AFDR histogram lengths are signed

`__cal_pvalue_qvalue_table` measures each run as `pos_array[j] - pre_p` with
`pre_p = 0`. A centred local-lambda window reaches below zero, so the first run's
length is negative and `N = sum(pscore_stat.values())` can be much smaller than
the chromosome. Reproducing this needs `i64` lengths.

## F170/F171 — negative FRAG coordinates need a shift

A `--format FRAG` control projects unclipped `d`-windows centred on each fragment
end (`x - d//2`), which is negative for contig-local fragments. A `u64` position
cannot hold it, so a counted run shifts every coordinate by `max(d, slocal,
llocal) / 2` before building signals and the caller subtracts it when writing.
`XlsRow`'s peak coordinates become `i64`.

## F172 — a zero local lambda is a runtime error, never a panic

`PeakDetect.__call_peaks_w_control` divides by the local lambda
(`se_edge/contig_edges` is the recorded example). Upstream raises
`ZeroDivisionError: float division`, so the observable behaviour is exit status
1 with **no output file** -- not a Rust panic (exit 101) and not a silent
continuation. `compute_uncached` now answers `+inf` and records the event on
`PScoreCache::hit_bad_lambda`; `callpeak` turns that into a runtime error before
any writer is opened.

## F175 — no control file: the "control" is the treatment, projected per class

`CallerFromAlignments.__init__` sets `self.ctrl = treat` when there is no control
file, so the local lambda is a projection **of the treatment track itself** --
`ctrl_d_s = [lregion]`, `ctrl_scale_s = [treat_length / (lregion * treat_total *
2)]`. Which projection depends on the track class: `PETrackI` gets
`pileup_from_PN_shifted`, `PETrackII` (`--format FRAG`) gets
`pileup_from_LRC_centers_as_list`. Using the ends projection for a counted track
collapsed the whole local lambda onto a single run, so every no-control FRAG
fixture reported the wrong `-log10(p)` (`sweep/gmini_mfrag_d1200_w180_noc_123`
reported 9634.71 against upstream's 6775.87).

`__call_peaks_wo_control` also sets the score calculator's `d = 0` and
`lambda_bg = treat_length / gsize` in PE mode -- not `d * treat_total / gsize`.

## F176/F177 — the peak start is `pos_array[above_cutoff[0] - 1]`, and the
## degenerate first region is dropped by `min_length`

`__chrom_call_peak_using_certain_criteria` fixes the first chunk's start to `0`
only when `above_cutoff[0] == 0`; otherwise the chunk starts at the **previous**
paired position. When the q-score of the very first row is finite but every row
after it is `NaN`, the first chunk becomes a one-row region spanning to a
negative coordinate; `__close_peak_wo_subpeaks` rejects it because
`peak_length < min_length`, and the reported peak starts at the *next* above-cutoff
row's predecessor. Observed on `sweep/gmini_mfrag_d1200_w180_noc_123`: upstream
starts the peak at 16, a naive implementation starts it at 0.

## F178/F179 — the p->q walk produces `NaN` q-values when the rank goes negative

`__cal_pvalue_qvalue_table` is

```
N = sum(pscore_stat.values()); k = 1; f = -log10(N); pre_q = 2147483647
for v in sorted(unique_values, reverse=True):
    q = v + (log10(k) + f)
    if q > pre_q: q = pre_q
    if q <= 0: q = 0; break
    pqtable[v] = q; pre_q = q; k += l
```

`log10` is `cython.cimports.libc.math.log10`, and `k` is a plain integer that goes
**negative** as soon as the negative first span (F169) is added. `log10(negative)`
is `NaN`, so `q` is `NaN`, `q > pre_q` is false, `q <= 0` is false, and the `NaN`
is stored. `NaN > cutoff` is false, so those positions are excluded from the
above-cutoff set. A textbook monotone-FDR implementation never produces `NaN`
here, and therefore places the peak start 16 bp early on every FRAG fixture whose
control window reaches below zero.

## F180 — the p-score track is the **union** walk, not the coincident boundaries

`__cal_pvalue_qvalue_table`, `__cal_pscore`, `__cal_qscore` and
`pos_array[above_cutoff]` all index `self.chr_pos_treat_ctrl`, which
`__chrom_pair_treat_ctrl` fills with one row per **union** position carrying
whatever each cursor holds. Building the histogram on only the positions where both
tracks change both loses rows (175 against 306 on
`sweep/gmini_mfrag_d1200_w180_noc_123`) and removes the negative first span, so
`k` never goes negative and F179's `NaN` window never appears. The peak-calling
path already used the union walk, so this was an internal inconsistency: the same
`chr_pos_treat_ctrl` was being reconstructed two different ways.

## F181 — `pre_p = 0` suppresses the run at coordinate zero

`_pileup_sorted_weighted_as_list` keeps `pre_p = 0` and only advances it **inside**
`if p != pre_p`, so deltas at coordinate 0 are applied but no run is emitted
there. A sweep that starts at its first boundary emits a leading `depth == 0` run
instead: one extra run, one extra union row, and a different signed first span in
the histogram. The suppression has to land on upstream's raw zero, so a
coordinate-shifted caller passes its shift as the sweep origin.

## F182 — `apply_multiple_cutoffs` returns a count, not a mask

`ret = arrays[0] > cutoffs[0]; ret += arrays[i] > cutoffs[i]`. With a single
criterion the result is the boolean array itself, so `np.nonzero` behaves as
expected; with several criteria (`['p','q']`) the result is an integer count and
`np.nonzero` selects every position that satisfies *any* cutoff. Recorded because
the multi-criterion path is reachable through the flag matrix.

## F183 — negative result: `average_template_length` is not the source of the
## residual PE control-lambda shift

The remaining paired-end failures are a uniform ~1.0e-5 relative shift in the
merged control lambda, visible in the `-B` bedGraph and worth a fourth decimal in
`-log10(p)`/`-log10(q)`. Upstream's `PETrackI/II.average_template_length` and
`Parser.d` are both `cython.float`, so rounding the mean to f32 is the obvious
suspect: `control_sum = control.total * 2 * average_template_length` would inherit
it. Measured, it is not: applying the f32 rounding **loses** 20 golden cases
(6358 -> 6338 of 7685) and leaves `sweep/gmini_mpe_d1200_w600_ctrl_066` unchanged,
so `ratio_treat2control = 0.4720812240955015` and `lambda_bg = 0.3850829875518672`
agree with upstream to the last printed digit while the in-run merged lambda does
not. A probe that calls the compiled `pileup_a_chromosome_c` with those exact
arguments reproduces **this** port's value, not the corpus's, so the divergence is
inside the real run and not in the factor inputs.

## F184 -- `control_sum` is truncated to an integer before the division

This was the ~1e-5 relative shift in every paired-end control lambda, and it is
**not** in the `.py`. `PeakDetect.__init__` declares

```
  long __pyx_v_control_total;
  long __pyx_v_treat_sum;
  long __pyx_v_control_sum;
```

so the compiled module assigns `control_sum` through `__Pyx_PyLong_As_long`,
which goes via `__Pyx_PyNumber_IntOrLong` and **truncates the product**:

```
    control_sum = int(control_total * average_template_length)
    ratio_treat2control = float(treat_sum) / control_sum
```

`average_template_length` is a `cython.float`, so in PE mode the truncation bites.
On `sweep/gmini_mpe_d1200_w600_ctrl_066` that is `int(394 * 99.79032135009766)`
= 39317 rather than 39317.386..., moving `ratio_treat2control` from
0.4720812240955015 to **0.47208586616476333**. Every control lambda scale derives
from it, and `lambda_bg` divides the same truncated value.

The single-end branch truncates too (`treat_sum = long(treat.total * self.d)`,
`control_sum = long(control.total * self.d)`), but `self.d` is an integer there,
so it is a no-op.

## F185 -- `--keep-dup auto` writes the "after filtering" header block

The block is gated on `options.keepduplicates != "all"`
(`callpeak_cmd.py:96`), not on the value being a literal, so `auto` resolves to a
number and prints all four lines -- `# <tag>s after filtering in treatment: N`,
`# maximum duplicate <tag>s in treatment = <resolved>`, and
`# Redundant rate in treatment: R`. Suppressing them for `auto` left every
`keepdup_auto` fixture's `*.xls` three lines short.

## F186 -- the `--nolambda` one-entry control still yields one row per treatment position

`CallerFromAlignments.pileup_treat_ctrl_a_chromosome` substitutes

```
    ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg], dtype="f4")]
```

which reads like "one row". It is not: the single entry's **position** is the
treatment's *last* position, so the ordinary union walk emits one row per
treatment position -- the treatment cursor always wins `min(t_p[i], c_p[0])` and
`lambda_bg` rides along -- until the treatment is exhausted. Collapsing it to one
row left exactly one scorable position and `--nolambda` called **no peaks at all**
on every fixture whose treatment has more than one run.

Separately, `--nolambda` sets `no_lambda_flag` from the **scale lists being
empty**, not from the absence of a control file, so the guard has to be on the
flag alone; matching it on `None` as well let a `--nolambda` invocation *with* `-c`
keep the real local lambda.

## F187 -- the counted control's window is one base wider than the uncounted one

The two paired-end track classes project the control differently and must not
share a parameter set:

| class | call | window for an end at `x` | span |
|---|---|---|---|
| `PETrackI` | `pileup_from_PN_shifted(l, r, d//2, d//2, ...)` | `[x - d/2, x + d/2)` | `d` even, `d - 1` odd |
| `PETrackII` (`--format FRAG`) | `pileup_from_LRC_centers_as_list` | `[x - d/2, x - d/2 + d)` | `d` always |

The second is `start_poss = concat(l - d//2, r - d//2)` with
`end_poss = start_poss + d`, i.e. `five_shift = d//2` and
`three_shift = d - d//2` -- the same split `FWTrack.pileup_a_chromosome` uses in
its non-directional branch. Using the narrow split for FRAG left every
right-hand window one base short: on `sweep/gmini_mfrag_d1200_w600_ctrl_231`
(`d = 119`) 137 of the 480 `-B` control-lambda rows were wrong and
`-log10(pvalue)` came out 1.18269 against upstream's 1.17793.

## F188 -- `maxgap` is `opt.tsize`, not `d`

```
    if opt.maxgap:
        self.maxgap = opt.maxgap
    else:
        self.maxgap = opt.tsize
    if opt.minlen:
        self.minlen = opt.minlen
    else:
        self.minlen = self.d
```

callpeak defines no `--max-gap`, so the merge gap is the measured tag size while
`min_length` is the predicted fragment length. They coincide in paired-end mode
and whenever `--extsize` happens to equal the mean read length, which is why the
distinction is invisible on most fixtures. `max_length = d` is still correct.

## F189 -- `-p/--pvalue` selects the p-value scoring function

`OptValidator` (`MACS3/Utilities/OptValidator.py:127`) resolves the two cutoffs
mutually exclusively, with `-p` winning:

```
    if options.pvalue:
        options.log_qvalue = None
        options.log_pvalue = log(options.pvalue, 10) * -1
    else:
        options.log_qvalue = log(options.qvalue, 10) * -1
        options.log_pvalue = None
```

and `PeakDetect.call_peaks` dispatches on which is set:

```
    if self.log_pvalue is not None:
        peaks = scorecalculator.call_peaks(['p',], [self.log_pvalue], ...)
    elif self.log_qvalue is not None:
        peaks = scorecalculator.call_peaks(['q',], [self.log_qvalue], ...)
```

Note `-p` takes the **natural** p value and the cutoff applied is `-log10(p)`;
`-p 1e6` therefore has a cutoff of -6, which every position passes, and upstream
calls peaks for it. The XLS header echoes the *flag* value, not the log.

`__cal_pvalue_qvalue_table` still runs either way -- `call_peaks` builds it before
dispatching -- so the XLS q column is populated under `-p` too; only the *criterion*
and the BED score column change. `callpeak_cmd.py:298` picks the score column the
same way:

```
    if options.log_pvalue is not None:  score_column = "pscore"
    elif options.log_qvalue is not None: score_column = "qscore"
```

and every BED-family writer emits `int(10 * peak[score_column])`, including
`write_to_broadPeak` and `write_to_gappedPeak`.

This port always scored against the q-value, so every `-p` invocation silently
used `-log10(0.05)`. The golden corpus contains no `-p` variant, so the gate could
not see it. Verified against the oracle afterwards on
`se_basic/gauss_two_peaks`: `-p 2/3/5/10` and `-p 1e6`/`-p 1e-6`/`-p 1e-40` all
match byte-for-byte, including the narrowPeak score column. `-p 1` still differs
in the peak boundary, which is the same region-edge defect as the remaining
corpus failures.

## F190 -- the counted projection has **no plus/minus mirror**

`pileup_from_LRC_centers_as_list` (the `--format FRAG` control projection) builds

```
    start_poss = np.concatenate((LRC['l'] - half_d, LRC['r'] - half_d))
    end_poss   = start_poss + d
```

so a fragment's left **and** right end each get the window
`[x - d//2, x - d//2 + d)`. `pileup_from_PN_shifted` -- the single-end and
`PETrackI` projection -- treats the two lists as plus and minus *strands* and
mirrors the second one (`N - three_shift`, `N + five_shift`). Sharing one
endpoint builder between them moved every counted window one base further left,
adding depth on the left of each right end and losing it on the right.

On `sweep/gmini_mfrag_d1200_w600_ctrl_231` that left 21 of the 480 `-B`
control-lambda rows wrong and `-log10(pvalue)` at 1.18269 instead of 1.17793.
`SingleEndParams::centred` now selects the unmirrored variant, and
`endpoints_from_positions` and `pileup_from_weighted_positions` both honour it.

## F207 -- `pileup_PV` did not merge runs, so every digested interval was one line

`pileup_PV` (`PileupV2.py:653-672`) does not emit one record per endpoint. It emits a
record only when the depth *changes*, and when the incoming value equals the previously
emitted one it extends the previous record's end instead:

```python
for i in range(PV_array.shape[0]):
    e = PV_array[i]['p']; v = PV_array[i]['v']
    if e != s:
        if z == pre_z:
            pileup_PV[c-1]['p'] = e        # extend the run
        else:
            pileup_PV[c] = (e, z)          # start a new one
            c += 1
            pre_z = z
    z += v
    s = e
```

This port walked the endpoint list and pushed a record per bin. On the ATAC fixture that
produced **23085** records for `a_digested_short.bdg` where upstream has 8720 -- which is
why upstream's file contains `chr1 510 1680 0.00000`, a single line covering 170 bins.
The two files were never comparable line-for-line, and `--save-likelihoods` was
structurally wrong for the same reason (78246 records against upstream's 3830).

**Fix.** `pileup_from_lr_hmmratac` now returns merged `(end, value)` runs exactly as
`pileup_PV` stores them, and `pileup_bdg_hmmratac` builds its tracks from those.

Two details that a "just merge equal neighbours" rewrite gets wrong, both caught by the
new tests:

* **Zero-depth runs are emitted, not skipped.** The depth is 0 between clusters, so a
  record is produced for the gap. Dropping it would splice two unrelated runs into one
  interval.
* **The merge is on the accumulated value, not adjacency.** Two *abutting* fragments of
  equal weight give depth 1.0 across the join and merge into a single run; two fragments
  separated by a gap give a zero run between them and do not.

**Result**: record counts now track upstream closely (8677 vs 8720, 18419 vs 18429,
20812 vs 20813, 13481 vs 13468), and the residual per-value deviation is <=1e-5.

## F212 -- `shrink_to_fit` on the signal tracks makes peak RSS *worse*

F210 shrunk the reads' `Vec`s. The obvious next step was the same for the tracks, which
are the larger object. Measured on the benchmark:

| build | after signals | peak RSS |
|---|---|---|
| reads shrunk only | 162.3 MB | **223.8 MB** |
| reads + tracks shrunk | 147.0 MB | **235.1 MB** |

Steady-state went *down* by 15 MB and peak RSS went *up* by 11 MB. Reverted.

**Why.** `shrink_to_fit` on a multi-hundred-megabyte allocation is not free: glibc either
`mremap`s in place or allocates a fresh block, copies, and frees the old one. Either way
the peak briefly holds **both** buffers, and the freed one is not returned to the OS, so
RSS stays at the high-water mark. That is the opposite of what "release the slack" assumes.

So F210's win is real but small (232 -> 224) because the reads' buffers are per-chromosome
and right-sized by the allocator. The tracks are not a slack problem at all -- they are
structurally 2x too wide. Fixing them means making `Run<f32>` 8 or 12 bytes, not asking
`Vec` to be tidier. Reverting; see F211 for the layout that actually has to change.

## F213 -- `callvar`'s likelihood model: two bugs a golden caught in ten minutes

`MACS3/Signal/VariantStat.py` is 549 lines of Cython and had no Rust equivalent at
all, so `callvar` could only ever write a VCF header. Ported as
`crates/macs-callvar/src/variant_stat.rs`, with `oracle/check_variant_stat.py`
generating a 410-case golden from the compiled `.so`.

Two defects, both in the port, both invisible to any tolerance-based test:

**1. `log1p(-e)` is not `ln(1 - e)`.** Upstream writes
`log1p(-exp(-q * LN10_tenth))`. The natural transcription is
`(1.0 - e).ln()`, and that is what I wrote first. They agree in exact arithmetic and
differ in the last `f64` bit. For a Q37 base `e ~ 1e-4`, so `1 - e` retains only ~12
significant digits and the error is *in the fourth significant digit of the answer*.
The golden failed immediately on the single-base cases. `ln_1p` it is now.

**2. Upstream *raises* where Rust returns `-inf`.** These functions run in CPython,
where `math.log1p(-1.0)` and `math.log(0.0)` raise `ValueError: math domain error`.
A Phred quality of `0` gives `e == 1.0`, so `log1p(-e)` is `log1p(-1.0)` -- upstream
**aborts the call**. A straight port returns `-inf`, writes it into a VCF's GQ
computation, and carries on producing garbage. Every logarithm is now range-checked
first and the failure is propagated as `Err`, so a `q=0` base (what several aligners
emit for a masked base) refuses the call exactly as upstream does.

**Also worth recording, because "fixing" either would be wrong:**

* `max_allowed_ar` is a Cython `cython.float`, i.e. **f32**. `0.99` therefore arrives
  as `0.9900000095367432`, and in `calculate_ln` the expression `1 - max_allowed_ar` is
  evaluated in C `float` arithmetic before `log` widens it. That subtraction decides
  which side of the clamp `r` falls on, so it is done in `f32` here too.
* `calculate_GQ` returns **343** for a far-from-diploid likelihood and
  `calculate_GQ_heterASsig` returns **386**. The `(int)` cast does not clamp to the
  VCF's 0-255 range, and neither does this port. Clamping would be an improvement and
  it would be a compatibility break, so it is asserted against in a test rather than
  tidied away.
* Upstream's `GreedyMaxFunctionAS` ends both sweep directions with an
  `if btemp: ... else: ...` whose two arms are character-for-character identical. That
  is very likely an upstream bug, and it has no effect on the result, so it is
  collapsed here with a comment rather than "fixed" -- changing it would cost a round
  trip to prove that nothing moved.

**Status: 410/410 cases bit-identical** (`to_bits()` equality on 6 floats each, plus
8 integer fields per case), including 23 cases where both sides must error.

## F214 -- `callvar`: an end-to-end oracle anchor, and what the model actually prefers

`PosReadsInfo`, `PeakVariants` and `RACollection` are all `cython.cclass`. Their
methods are not reachable from Python, so unlike `VariantStat` (F213) there is **no
way to test them one function at a time** -- any differential has to drive a whole
`callvar` run. Before writing 2000 lines against no feedback, that was worth fixing.

`oracle/check_callvar.sh` now does it, against upstream's own test data
(`CTCF_PE_ChIP_chr22_50k.bam`, `CTCF_PE_CTRL_chr22_50k.bam`,
`callvar_testing.narrowPeak`, all read out of the pinned oracle tree so no 10 MB BAM
has to be committed). The oracle run is deterministic -- verified by running it twice
and diffing -- and yields **16 variant records**, committed as
`crates/macs-callvar/tests/data/callvar_variants.golden`. That is a complete target:
every field (POS/REF/ALT, QUAL, `M`, `MT`, per-allele depth, the `SB` quadruple, all
four BICs, `AR`, and `GT:DP:GQ:PL`) is checked.

**`PosReadsInfo` is now ported** (`crates/macs-callvar/src/pos_reads_info.rs`), with
37 unit tests. Two findings worth writing down:

**1. Dict insertion order decides the ALT column.** Upstream picks the top two alleles
with `sorted(self.n_reads, key=self.n_reads_T.get, reverse=True)[:2]`. `sorted` is
stable and `reverse=True` does *not* reverse equal elements, so ties fall back to the
insertion order of `n_reads` -- `{ref_allele, A, C, G, T, N, *}` with the reference
first. With `ref = T` and an exact 10-vs-10 tie between `A` and `AGGG`, the ALT column
comes out `A,AGGG`, not `AGGG,A`, and `MT` is correspondingly `SNV,Insertion` rather
than `Insertion,SNV`. The port keeps a `Vec` of entries rather than a `HashMap` for
exactly this reason.

A second consequence: `sorted` runs over *all* seven seeded alleles, so `top2` is very
often an allele with **zero** treatment reads. That is why upstream tests
`n_reads[top2] == 0` in three separate places, and why a reference-only position ends
up `homo_ref` and filtered.

**2. The 2-nat BIC margin is hard to clear, so most calls are `heter_unsure`.** The
oracle's 16 records split 10 `heter_unsure` / 6 `heter_noAS` and **not one** `homo` or
`1/2`. Measuring the model explains it: with 10 reads of one allele and 10 of another
and no control, `BIC_homo_major` = 138.2 against `BIC_heter_noAS` = 34.2 -- a wide
win. But with 11-vs-1, `BIC_homo_major` = 13.8 and `BIC_heter_AS` = 13.4, a 0.4 gap,
so the `+ 2` guard rejects every model and the position becomes `unsure`/`unsure` and
is filtered. My first unit tests "expected" `1/1` and `1/2` calls from lopsided input;
they were wrong about upstream, not the port. All 16 golden records also share the
signature `PL = <x>,0,<y>` -- the middle entry is zeroed for every `heter_*` call
because the winning model is the reference -- which is now a unit-tested invariant
rather than something noticed by eye.

**One real bug the tests caught**: `self.type.startswith("heter")` is evaluated
*after* the BIC ladder, once a type exists. Reading it before the ladder (as a
"capture the previous value" optimisation) sees the empty initialiser and silently
drops every heterozygous call -- `gt` stays `unsure` and `ALT` stays `.`.

**Still open**: `PeakVariants` (402 lines) and `RACollection` (906 lines), then the
fermi-lite assembler. `is_kernel_implemented()` remains `false` and the command still
refuses with exit 1 *before* creating any file, which is the correct behaviour for an
incomplete kernel.

## F215 -- `PeakVariants`: the `f32` narrowing is observable, and one stale-read bug

`PeakVariants` is the last piece before `RACollection`. It is a plain data class plus
the per-peak indel repair, and both halves carry traps.

**1. The six `%.2f` fields are `f32`, and it changes the output.** `Variant.__init__`
takes `deltaBIC`, the four BICs and `AR` as `cython.float`, but `PosReadsInfo.toVariant`
hands them `double`s. `toVCF` then prints them with `%.2f`. `f32` has ~7 significant
digits, so the narrowing lands *inside* the second decimal for values of that
magnitude, and the rendered string changes:

| value | as `f64` | as `f32` |
|---|---|---|
| 23.215 | 23.21 | **23.22** |
| 13.535 | 13.54 | **13.53** |
| 138.185 | 138.19 | **138.18** |
| 103.985 | 103.98 | **103.99** |

Every one of those is a `DBIC=`/`BIC*=.*=` value that could appear in a VCF, so a port
that keeps them as `f64` produces a visibly different INFO column. `Variant` stores
them as `f32` and widens only inside `format!`, mirroring Cython passing a C `float`
to `printf("%f")`. The test asserts both renderings, so a regression to `f64` fails
rather than passing unnoticed.

**2. `is_indel` is a substring test; `is_only_del` is equality.** `find("Deletion")`
against `== "Deletion"`. So `Deletion,SNV` is an indel but is neither "only a
deletion" nor "only an insertion" -- and all three `fix_indels` passes key off the
`is_only_*` predicates. A multiallelic indel is therefore deliberately left alone,
which is only obvious once you notice the two functions disagree.

**3. Pass 2 reads the *original* variant, not the copy it just rewrote.** Upstream:

```python
self.d_Variants[p-1] = copy(self.d_Variants[p])
self.d_Variants[p-1]["ref_allele"] = rs[p - self.start] + ...
self.d_Variants[p-1]["alt_allele"]  = rs[p - self.start]
if self.d_Variants[p].top1isreference:      # <-- the OLD one, at p
```

My first port tested the rewritten copy. That is wrong almost always: prepending the
anchor base to REF makes `ref != top1`, so the branch never fired and `top1`/`top2`
were left holding the stale `*`. The port now captures `top1isreference`/`top2isreference`
before the mutation.

**4. Two coordinate conventions in one module.** `toVCF` emits `str(p + 1)` while
`fix_indels` indexes `rs[p - start]`. Both are 1-based vs 0-based in the same file.

**5. `fix_indels` can raise, and must not panic here.** Pass 1 evaluates
`d_Variants[p0]` with `p0 == -1` and `p == p1 + 1 == 0` -- a `KeyError` upstream, for a
pure deletion at absolute position 0. `fix_indels` returns `Result` and surfaces it,
which keeps the "zero panics" contract without inventing output upstream never
produces. A short `refseq` is handled the same way.

**Status: 51 unit tests across the three ported modules.**

## F216 -- every deletion is a "tip", and `callvar` writes duplicate records

Porting the read-level query (`ReadAlignment.get_variant_bq_by_ref_pos`, now
`ReadAlignment::variant_bq_by_ref_pos` in `macs-io`) surfaced two upstream
behaviours that a "sensible" port would have silently corrected.

**1. `pos` is never assigned on the deletion path, so `pos == 0` and `tip` is always
true.** The walk assigns `pos` only inside the `M`/`=`/`X` branch -- and that branch
`break`s. So the `D`/`N` branch is reachable only when no match was ever recorded, and
`pos` reads back as `0`. Since `tip` is `pos == 0 or pos == len - 1`, **every deletion
is reported as a read tip**, at any depth in the alignment.

This is not cosmetic. `tip` feeds `update_top_alleles`'s
`n_t[allele] - n_tips[allele]` alt-allele count, so a `*` allele is discounted as if
every read supporting it were a read end -- which makes deletions far harder to call
than the code appears to intend. Replicated exactly, with the reasoning in the source,
because the alternative is a plausible-looking port that silently diverges.

**2. Two byte-identical input peaks produce two identical VCF records.** The oracle's
16 records cover only **13 distinct positions**. `callvar_testing.narrowPeak` contains
`run_callpeak_bampe_narrow_peak_7a` and `run_callpeak_bampe_narrow_peak_7b` at the same
coordinates; each is processed independently and writes its own lines, and there is no
cross-peak deduplication. Three positions consequently appear twice.

That fell out of the 16-record golden, and it immediately exposed a hole in my own
harness: `oracle/check_callvar.sh` compared with `grep -qxF`, i.e. set membership, so a
port that emitted each record exactly once would have *passed*. It now compares as a
multiset (counts must match) and reports records we emit that the oracle does not.

**Also: `res < op_l - 1`, not `res < op_l`.** The M branch uses an off-by-one-looking
split (`res < op_l - 1` takes a base mid-block, `res == op_l - 1` takes the last base
and then looks for a following `I`). Writing it as `res < op_l` looks equivalent and
shifts every base by one; 14 differential tests now cover the six operator branches and
their boundary conditions.

## F217 -- `RACollection`: the consensus is built by a slice assignment that resizes

The last unported piece is `RACollection`, and it turned out to hide the single most
surprising thing in `callvar`.

**`callvar` has no genome.** `__get_peak_REFSEQ` builds the peak's reference sequence
**from the reads themselves**: every gap in the treatment pileup is back-filled with the
*previous* read's reference bases, then the control pileup back-fills over the top.
There is no FASTA and no `--genome`; `peak_refseq[i - left]` indexes that consensus.

**The fill writes the previous read's bases into the current read's span.**

```python
read = ralist[i]
if read["lpos"] > prev_r:
    read = ralist[i - 1]                  # previous read
    read_refseq = read.get_REFSEQ()       # its bases
    ind   = read["lpos"] - start          # current read's lpos
    ind_r = ind + read["rpos"] - read["lpos"]
    seq[ind: ind_r] = read_refseq
```

The span is sized by one read and filled with another's, so the two lengths routinely
disagree -- and `seq` is a `bytearray`, so Python's slice assignment **changes its
length**, shifting every later coordinate. My first port rejected the mismatch, and it
failed on the very first peak of upstream's own `callvar_testing` BAM: a 100-byte slice
receiving 101 bytes. The port now reproduces the resize (`py_splice`), including
Python's clamping of both slice ends and its treatment of a reversed slice as empty.

**Also in this module:**

* `remove_outliers` computes its percentile over **both samples combined**, so noisy
  control reads can raise the threshold and evict clean treatment reads. The test that
  pins this discriminates combined-vs-per-sample rather than just checking a count.
* It runs *after* the collection is built, so the consensus still comes from the
  unfiltered reads. Re-deriving the consensus afterwards would change the reference
  bases the calls are made against.
* `get_REFSEQ` walks `MD` with raw byte-range tests, so a stray space is an error.
  Upstream raises; `refseq()` returns `Err`, keeping the no-panic contract.

**State: ported and unit-tested (17 tests), wired into the CLI, but not yet exact.**
`--fermi off` runs end to end and produces a VCF; against the oracle it currently emits
2890 records where upstream emits 22, so the reproduced consensus still differs
somewhere. `callvar` therefore continues to **refuse, before writing any file** -- a
VCF full of wrong variants is worse than no VCF, because it looks like an answer. The
two remaining `-F off` targets are committed
(`crates/macs-callvar/tests/data/callvar_variants_nofermi.golden`, 22 records) so the
next attempt is measured against a known answer rather than eyeballed.

**Also measured: `--fermi` changes the answer.** Same inputs, three settings:

| `-F` | records |
|---|---|
| `off` | 22 |
| `auto` (default) | 16 |
| `on` | 15 |

`auto` *replaces* the no-assembly calls for any peak carrying an indel or a
reference-biased het, so it is not a superset of `off`.

## F218 -- `__fill_refseq` was not mixing the two reads; I was

Last turn this port emitted 2890 variant records where the oracle emits 22, and the
consensus looked wrong everywhere. The build-up — `VariantStat` bit-identical on 410
cases, `PosReadsInfo`, `PeakVariants`, `variant_bq_by_ref_pos`, `get_REFSEQ` — all
verified, yet the answer wrong *everywhere*, not at one peak. That shape is the
tell: a systematic offset, not a statistics bug.

So before touching anything else I dumped the intermediate the VCF diff could not see.
`RACollection.__getitem__` exposes `peak_refseq` and `peak_refseq_ext`, and the class
is a `cdef`, so those are otherwise unreachable from Python. `oracle/dump_callvar_refseq.py`
builds a `RACollection` per peak inside the pinned oracle and writes both buffers;
`macs-callvar`'s `dump_refseq` example emits the same fields. Comparing them:

```
ext len 447/447 first diff @0 | ora 'ATTATTTAGCACTG' ours 'NNNNNNNNNNNNNN'
```

Ours was `N` at index 0 on **every** peak. `start = min(RAs_left, left)`, and the
first gap-filling write should land at `ralist[0].lpos - start == 0`.

And here is the bug, which was mine. I had read this:

```python
read = ralist[i]
if read["lpos"] > prev_r:
    read = ralist[i - 1]        # rebind
    read_refseq = read.get_REFSEQ()
    ind   = read["lpos"] - start
    ind_r = ind + read["rpos"] - read["lpos"]
```

as "the previous read's bases go into the *current* read's span" — a genuine upstream
oddity worth a careful note. It is not that. The rebinding of `read` is immediately
followed by *every* remaining use of `read`, including `lpos` and `rpos`, so the write
places read `i-1`'s own bases at read `i-1`'s own coordinates. The rebinding is
redundant, not a mix-up. I ported the mixing version, which is off by the distance
between consecutive reads — hence `N` at the start and a wrong reference base everywhere
downstream.

The same applies to `prev_r = read["rpos"]`: it is read `i-1`'s `rpos`, and it is
deliberately *not* advanced for overlapping reads, so once a gap is found the condition
stays true for a while. That is reproduced too, with a comment.

**After the fix:**

```
10/10 peaks byte-identical in peak_refseq_ext
10/10 peaks byte-identical in peak_refseq
22/22 variant records identical (--fermi off), including the 3 emitted twice
```

`--fermi off` is now **byte-identical** to the pinned oracle, all ten peaks, checked by
`oracle/check_callvar.sh` and pinned by `crates/macs-callvar/tests/peak_refseq_diff.rs`.

**What this says about the last turn's note.** I wrote that the fill "writes the
previous read's bases into the current read's span" and built an elaborate mechanism
(`py_splice`) around it. The mechanism is still needed — Python's `bytearray` slice
assignment really does resize, and on three of the ten peaks our `ext_len` differed —
but the premise was wrong, and had I compared the intermediate before writing any of
it I would have spent the effort on the right bug. A VCF that is wrong *everywhere* is
evidence of a systematic offset; check the stage that feeds it before theorising about
the stages below.

**Also pinned, so a future divergence is not 2890 records of noise:**
`crates/macs-callvar/tests/data/callvar_peak_refseq.golden` holds, per peak, both
buffers in hex. And `--fermi auto` is now asserted to be *refused with no output file*
— upstream re-calls indel and reference-biased peaks from the assembly rather than
adding to them (`-F off` 22, `auto` 16, `on` 15), so answering `auto` with the
no-assembly result would be a different wrong answer, not a partial one.

## F219 -- bridging fermi-lite: two ownership rules, both of which corrupt memory if missed

`--fermi on|auto` needs a de novo assembler, and `PORTING_PLAN.md` says exactly how to
get one: *"callvar may initially bridge to upstream fermi-lite via thin C FFI; a pure-Rust
assembler is the last item, after every other gate is green."*

So fermi-lite is now vendored verbatim at the revision the oracle pins
(`MACS3/fermi-lite`, MIT, `FML_VERSION "r53"`, recorded in `vendor/fermi-lite/PINNED.md`)
and compiled by `crates/macs-callvar/build.rs`. The FFI surface is five functions and
three structs: `fml_opt_init`, `fml_assemble`, `fml_utg_destroy`, `fml_opt_t`,
`bseq1_t`, `fml_utg_t`. Nothing else is reachable.

**Two ownership rules, both learned the hard way.**

**1. Qualities go in as Phred+33.** `fml_assemble` expects ASCII and upstream does
`cqual[j] = tmpq[j] + 33` on raw BAM qualities. Passing raw values silently changes
error correction and therefore every unitig -- no crash, just a different answer.

**2. The `bseq1_t` array must come from libc `malloc`, not a Rust `Vec`.**
`fml_assemble` frees the array it is handed (upstream comments it twice: *"we rely on
fermi-lite to free this mem"*), so a `Vec`'s buffer gets freed by libc and again by the
`Vec`'s own drop. My first version did exactly that and the test binary aborted with
`free(): double free detected in tcache 2`. The same applies to each `seq`/`qual` pair.

**Two of my own assertions were wrong, and the library corrected them.**

* I wrote `assert_eq!(size_of::<Bseq1T>(), 4 + 2 * size_of::<*mut c_void>())` -- 20. C
  pads that `int32_t` out to pointer alignment, so it is **24**, and `repr(C)`
  reproduces the padding rather than assuming it away. Guessing the layout instead of
  copying the header would have produced silent misreads.
* I passed `min_asm_ovlp = 10_000` to test "an impossible overlap". That value is used
  **directly as a k-mer length** in `bfc_ch_init`, which asserts `k <= 63`, so the
  process aborted. The reachable range is 1..=63 and the test now walks it.

**What I did not pin, and why.** My first fixture -- 40 reads of 100 bp tiling a
template every 15 bp -- produced **zero** unitigs, while a denser set produced one.
Rather than reverse-engineer fermi-lite's cleaning thresholds (which derive `min_elen`
from total length and read count, and interact with `trim_len`/`max_bdist`), the tests
now pin the **FFI contract**: struct layout, no double free over many calls, the
reachable overlap range, and that assembled bases come from the input (which catches a
buffer read after being freed). Counting unitigs against synthetic data would have been
pinning an accident. The behavioural contract is the `-F auto` VCF, compared end to end.

**State: the bridge builds, links and runs clean. The orchestration around it is not
ported** -- `build_unitig_collection`, `align_unitig_to_REFSEQ` and `UnitigCollection`
remain, so `is_assembly_implemented()` is still `false` and `-F auto` still refuses.
That is the honest boundary: the assembler is available, the code that feeds and reads
it is not.

## F220 -- `align_unitig_to_REFSEQ`: the second vendored C file, and a `main` that had to be renamed

`--fermi auto` needs a unitig placed back on the peak consensus. Upstream does that with
`MACS3/Signal/swalign.c`'s `smith_waterman` -- local alignment, **without** an affine gap
penalty -- which is a second vendored C file alongside fermi-lite. Vendored verbatim
(`diff` against the oracle is empty) and built by the same script.

**The scoring scheme is not Smith-Waterman's usual one.** Per `swalign.c`: match **+2**,
mismatch **-3**, gap **-5**, second gap in a run **-2**. `verify_alns`'s docstring gives
the operational reading: a score of 150 over 100 bp is "10 mismatches within 100bps".
Re-deriving that matrix in Rust would be a different algorithm with different
tie-breaking, and the placement feeds straight into the calls, so the C stays.

**The vendored file ships a CLI, and its `main` collides with the binary's.** The obvious
fixes are editing the file or moving it, and both are worse than they look: the point of
vendoring is that `vendor/` stays byte-identical to the oracle so a diff shows nothing.
So the symbol is renamed at build time instead:

```rust
.flag("-Dmain=swalign_cli_main_unused")
```

A build flag, not an edit. Its own `print_alignment` becomes unreferenced, which `-w`
hides. The vendored bytes are untouched, and that is asserted by a `diff` against the
oracle tree.

**The forward/reverse-complement branch mutates the caller's list.** Each unitig is
aligned twice and the better score wins; on a tie the **forward** alignment is kept
(strict `>`). When the reverse complement wins, upstream overwrites `unitig_list[i]`
with the revcomp, because `remap_RAs_w_unitigs` later re-derives read positions against
the possibly-flipped unitig. A port that aligned both and kept only the scores would be
subtly wrong in the next step rather than visibly wrong here.

**One of my own assertions was wrong and the library corrected me.** I asserted that a
disjoint pair would score *below* zero, reasoning that mismatches must dominate. But
Smith-Waterman is **local**, so its score is never negative -- the worst case is a single
matching base, which is exactly the `+2` observed. The test now asserts `>= 0` and that
it is far below a real alignment, which is the property `verify_alns` actually relies on.

**State.** fermi-lite assembly and unitig-to-consensus alignment are both ported and
tested (8 alignment tests, including 200 repeated calls to exercise the free path).
Still missing: `verify_alns`, `filter_unitig_with_bad_aln`, `remap_RAs_w_unitigs`,
`add_to_unitig_list`, `build_unitig_collection` and `UnitigCollection` (324 lines).
`is_assembly_implemented()` stays `false` and `-F auto` still refuses -- unchanged and
correct for the state of the port.

## F221 -- the assembly path now runs, and is gated on *measured* exactness

The remaining fermi orchestration is ported: `verify_alns`, `remap_RAs_w_unitigs`,
`add_to_unitig_list`, `build_unitig_collection`, `UnitigRAs`, `UnitigCollection`
(`crates/macs-callvar/src/unitig.rs`, 13 tests). `--fermi auto` now runs end to end.

**Two rules decide everything downstream, and both are easy to get subtly wrong.**

`verify_alns` normalises by **alignment columns, not unitig bases**:

```python
if aln_scores[i] * 100 / len(markup_alns[i]) < min_score_100:   # 150
```

Gap columns are in the denominator. A test now pins exactly this: a 10-column
alignment of which 5 are gaps and which scores 12 gives `12*100/10 = 120` and is
**dropped**, whereas normalising by the 5 non-gap bases would have given 240 and kept
it. Getting the denominator wrong changes which unitigs survive, so it changes the
calls.

A read is assigned to the **first** unitig containing it, on an exact substring test
(`tmp_ra_seq in unitig`) -- no alignment, no partial match -- and unmapped reads go
round again: first at the full overlap, then *unconditionally* at half. New unitigs are
placed **before** existing ones, which reorders the first-match-wins assignment.

`UnitigRAs.get_variant_bq_by_ref_pos` walks residues counting only **non-gap reference**
columns, so `index_aln` can land on a `-` in the unitig -- which is exactly the deletion
case, reported as `*` at quality 93 for every mapped read. An insertion is a gap in the
*reference*, not the unitig; I had that backwards first and silently got the matched base
alone.

**Where it stands, measured.** `oracle/check_callvar.sh --fermi auto` matches **7 of 16**
records. The three peaks needing assembly are called at the **right coordinates** with
well-formed 10-column records, but their treatment depths differ from upstream's -- so
our unitig-to-read assignment diverges somewhere in the chain.

**So `--fermi on/auto` still refuses, and `is_assembly_implemented()` still returns
`false`.** I made that a one-line flip, documented as such, with the 7/16 measurement
recorded next to it.

**A mistake worth recording.** When I wired the path in, I *replaced* the refusal gate
with `let use_assembly = ... && is_assembly_implemented()` and deleted the early return
above it. With the flag false, `-F auto` then fell through to the `-F off` code and wrote
a VCF -- reporting **22** records under `-F auto`, where the oracle gives 16. That is the
precise failure mode this whole command has been guarded against: not an error, not a
crash, but a confident wrong answer. The checker caught it ("`--fermi auto` was accepted
but needs fermi-lite assembly"). The gate is now a separate early return that cannot be
bypassed by the feature flag below it.

## F222 -- one missing branch, and a divergence narrowed to a single function

Last turn's assembly path was at 7/16. Two things moved it to **10/16**, and the
remaining failure is now narrow enough to name a function.

**1. A branch I had not implemented at all.** `--fermi auto` has a second use of the
assembly that is easy to miss, because the log line never mentions it:

```python
if (fermi == "auto" and (not peak_variants.has_indel()) and peak_variants.has_refer_biased_01()):
    for i in peak_variants.get_refer_biased_01s():
        PRI = unitig_collection.get_PosReadsInfo_ref_pos(i, ref_nt, Q=minQ)
        PRI.update_top_alleles(...); PRI.call_GT(...); PRI.apply_GQ_cutoff(...)
```

When the assembly ran but the peak has **no** indel, upstream re-calls every
reference-biased `0/1` from the unitig consensus instead of re-calling the peak. That is
the whole point of the stage: a reference-biased het is exactly the artefact local
assembly corrects, so the no-assembly answer is known to be biased toward the reference.
Two details: the order is `update_top_alleles` -> `call_GT` -> `apply_GQ_cutoff` with
**no** `raw_read_depth == 0` skip (that belongs to `call_variants_at_range`), and
`has_refer_biased_01()` takes no argument upstream, so `is_refer_biased_01`'s own default
of `0.85` is what applies.

**2. The remaining divergence is in `remap_RAs_w_unitigs`, and the evidence says so.**

All six still-differing records are the same three positions counted twice -- the input's
two byte-identical peaks. And at each of them:

| position | | oracle | ours |
|---|---|---|---|
| 17652406 | DPT / DPC | 15 / **3** | 12 / **3** |
| 17652461 | DPT / DPC | 19 / **2** | 8 / **2** |
| 17653030 | DPT / DPC | 17 / **1** | 11 / **1** |

**The control depths match exactly.** `DPC` matching to the base across three positions
means the unitig list, the alignments, the `lpos`/`rpos` arithmetic, the indel repair and
the *control* assignment are all correct. Only treatment reads are missing, so the loss is
in the treatment half of `remap_reads_with_unitigs` -- reads that should have matched a
unitig by exact substring and did not.

Ruled out already: `SEQ` construction. Upstream's `__get_SEQ_QUAL` appends both nibbles
then trims a single trailing `=` if present, which is exactly what `ReadAlignment::sequence`
does, so a length or content mismatch there is not the cause.

**Still gated off.** `is_assembly_implemented()` returns `false`, `-F on|auto` refuses,
and the refusal text now carries the 10/16 figure and the DPC/DPT localisation so the
next attempt starts from a diagnosis rather than a hunch.

## F223 -- the reads exist and are assigned; the question is *which* unitig

Continuing F222 rather than starting something new, because the remaining `--fermi auto`
divergence is one function and the VCF cannot show why.

**What the diagnostic shows** (`cargo run -p macs-callvar --example unitig_stats`, added
for this). On upstream's own peak 7:

```
reads after outlier removal: T=147 C=17
unitigs: 7
  lpos=17652311 rpos=17652552 span=241 readsT=4  readsC=3
  lpos=17652312 rpos=17652553 span=241 readsT=6  readsC=0
  lpos=17652391 rpos=17652553 span=162 readsT=7  readsC=0
  lpos=17652416 rpos=17653027 span=611 readsT=87 readsC=7
  ...
reads assigned to some unitig: T=147 C=17
```

**All 147 treatment reads are assigned to some unitig**, and **every one of the 17
control reads too** -- so nothing is lost in `build_unitig_collection`'s read handling,
and `remap_reads_with_unitigs` is not dropping reads outright.

Three unitigs overlap position 17652406 and hold 4 + 6 + 7 = **17 candidate treatment
reads** between them. Upstream reports DPT 15 there; we report 12. The shortfall is
therefore produced by the *per-position* filter inside
`UnitigRAs::get_variant_bq_by_ref_pos`:

```python
ra_pos = index_unitig - self.seq.find(ra_seq) - 1
if ra_pos < l_read and ra_pos >= 0:
```

**The hypothesis this leaves.** `remap` assigns a read to the **first** unitig that
contains its sequence, and three of these unitigs overlap with near-identical spans
(17652311, 17652312, 17652391). A read assigned to unitig A rather than B gets a
different `ra_pos` -- `find` is relative to `A.seq`, and `index_unitig` is relative to
A's alignment -- so it may or may not satisfy `0 <= ra_pos < l_read` at this position.
The three depths differ by exactly 3, and moving three reads between overlapping unitigs
is the size of change that produces.

That points at the unitig **list order**, which is set by `add_to_unitig_list` placing
second-round unitigs at the *front* (`new_unitig_list.extend(unitig_list)`), and at
whether our second-round unitigs are the same set as upstream's -- which in turn depends
on the unmapped-read list from the first round. Neither is verifiable from our side
alone, since upstream's unitig list is not exposed.

**State unchanged and gated:** `is_assembly_implemented()` is `false`, `-F on/auto`
refuses, `-F off` is 22/22. The diagnostic example is kept so the next attempt starts
from these numbers rather than from the VCF diff.

## F224 -- `callvar` is done: the bug was branch structure, not arithmetic

`--fermi auto` went from 10/16 to **16/16** by fixing a control-flow error, not a
numeric one.

**What upstream actually does when assembly is attempted.** The no-assembly variants are
computed and then held **unwritten**, because the assembly check comes *before* any write:

```python
if (fermi == "auto" and (peak_variants.has_indel() or peak_variants.has_refer_biased_01())) or fermi == "on":
    unitig_collection = ra_collection.build_unitig_collection(fermiMinOverlap)
    if unitig_collection == -1:  continue                 # peak dropped entirely
    elif unitig_collection == 0: ... fall back to previous results
    # auto, no indel, refer-biased -> revisit ONLY the reference-biased positions
    # otherwise                    -> reset and re-call the whole peak from unitigs
else:
    # no assembly attempted -> write now
```

So there are **four** outcomes, and in the "revisit" case the no-assembly variants survive
with only the reference-biased positions replaced. `fix_indels` is applied at write time,
once, in whichever branch is taken -- never before.

**My version replaced the whole peak with the assembly-derived one** in both `auto`
sub-cases. Two consequences, and the second is what cost the most time:

1. Correct no-assembly calls were discarded whenever the assembly ran.
2. `fix_indels` was applied eagerly, so the variants the revisit was supposed to edit had
   already been merged and had positions removed.

**The symptom pointed somewhere else entirely.** The wrong branch showed up as DPT 12
where upstream reports 15 -- a *read-assignment* symptom. I chased it through
`remap_reads_with_unitigs` and then through the per-position `ra_pos` range filter inside
`get_variant_bq_by_ref_pos`, and wrote two findings (F222, F223) doing so. Both were
reasoning about the right layer from the wrong layer's data: the reads were never lost,
because those reads belonged to a no-assembly call that should never have been thrown
away in the first place. F223's diagnostic -- "all 147 reads are assigned, three unitigs
hold 17 candidates, we surface 12" -- was accurate, and still described a peak whose
variant set came from the wrong source.

**The lesson, and it is not about care.** The evidence said "depth shortfall at three
positions, control depths exact". The correct next question was not "which unitig does
this read belong to" but "**which code path produced these records at all**". A depth
that is too low is equally consistent with "the wrong variant set" as with "the right
variant set, missing reads", and I resolved that ambiguity by assuming the latter because
it was the hypothesis I had already invested in. Comparing the `-F off` and `-F auto`
goldens directly -- 22 records versus 16 -- would have shown that peak 7's records are
*supposed* to come from the no-assembly path, and that was visible before F222.

**State: `callvar` is byte-identical in both modes**, on upstream's own fixtures:

```
-F off    22/22 records identical
-F auto   16/16 records identical
```

`is_assembly_implemented()` is now `true`. The early-return gate in `run()` is kept
anyway, as an invariant rather than a dead branch: it is the thing that stops a future
regression from silently emitting a wrong VCF, which is exactly what happened in F221
when the check was folded into a feature flag.

## F226 -- a `BedGraph` with the wrong `Genome` panicked on valid input

`hmmratac --save-digested` on a BEDPE run aborted with
`index out of bounds: the len is 0 but the index is 1` in `Genome::name`.

`BedGraph::new(baseline)` starts with an **empty** `Genome`, and the digest writer built
one that way and then inserted tracks keyed by `petrack`'s `ChromId`s -- a different
dictionary. `chroms_sorted` then indexes the empty `names` vector. Ids are positional, so
the two dictionaries are interchangeable *only if they are the same one*.

`BedGraph::with_genome` now exists for that. Two of the three sites in `hmmratac` had the
same defect; `apply_fold_change` was the third and would have panicked next.

The backtrace was misleading: the frame nearest the panic was attributed to
`pileup_bdg_hmmratac`'s sort closure, and `genome.is_empty()` reported `false` with 17
names present. The real frame was `BedGraph::chroms_sorted`, several levels up, reached
after the digest tracks were already built. **A panic blamed on the wrong function is
worth re-deriving from the full trace rather than the first frame.**

## F227 -- `_accessible_regions.narrowPeak` had the right summit offset and the wrong file

Columns 1-5 matched upstream exactly and column 10 (the summit offset) matched too. Columns
6-9 did not:

```
oracle  chrI 20660 21150 MACS_peak_1 33 . 0        0        0        341
ours    chrI 20660 21150 MACS1_atr    33 . 3.386310 20660    21150     341
```

`PeakIO.write_to_narrowPeak` is ENCODE narrowPeak (BED6+4):

```python
"%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n"
  chrom start end  name int(10*score) .    fc     pscore qscore summit-start
```

Columns 7-9 are three **`%.6g` values**, not coordinates. `save_accessible_regions` adds
its peaks with only `(chromosome, start, end)` and never sets `fc`, `pscore` or `qscore`,
so all three are the default **0** -- the fold change appears only in column 5 as
`int(10 * score)`.

Two smaller things in the same function:

* **The peak name ignores `--name`.** `write_to_narrowPeak(fhd)` is called with its
  defaults, so the prefix is `b"MACS_peak_"` and the dataset label is the default
  `b"MACS"`. Ours interpolated the user's name (`MACS1_atr`).
* `%.6g` is **six significant digits**, not six decimals -- `{:.6}` renders `3.38631` as
  `3.386310`. The existing `macs_io::peakout::format_g` already implements it exactly.

Columns 1-5 are now byte-identical to the oracle.

## F228 -- G13 is further from green than the status claimed

`oracle/check_hmmratac.py` now runs both sides on upstream's own `yeast500k` BEDPE data
with an exported model (the inference path the criteria require to be exact) and reports
every component, so the number is reproducible rather than measured by hand.

```
regions:          oracle 9286  ours 9278  shared 7371
Jaccard (region): 0.6585   requirement >= 0.98
Jaccard (base):   0.6488
unmatched regions: oracle-only 1915, ours-only 1907
one-edge-matched regions: 1782, far edge off by median 10 / max 1800 bp
summit shift:     max 0 bp over 7385 identically placed regions
```

**The recorded "0.973" was measured on a 946-region corpus; on yeast500k it is 0.66.**
Both numbers were real; the smaller one was not representative, and nothing recorded
that. That is the argument for the harness rather than for a different number.

What the shape says: **summit offsets are exactly right** (0 bp over 7385 identically
placed regions) and 1782 regions match on one edge with the other edge off by a **median
of 10 bp** -- one bin. So this is a states-path bin-edge problem, not an HMM or statistics
problem. The inference path is close; the region boundaries are not yet.

## F229 -- where G13 actually is: EM matches exactly, the digest does not

Localising the 0.66 properly, rather than leaving it as a number.

**What is now exactly right.**

* The EM fragment-length fit, on the same input, both sides:

  ```
  oracle  means:  50  166.7  331.2  460.3   stddevs: 20  48.1  34.1  34.9
  ours    means:  50.0000  166.7000  331.2000  460.3000
          stddevs: 20.0000  48.1000  34.1000  34.9000
  ```

  And upstream uses those EM means for `generate_weight_mapping` **even when
  `--model` is given** -- the model file only replaces the HMM, not the digest weights.
  So running EM unconditionally here is correct, not a bug.
* `_make_bdg_of_bins_from_regions`, including the global-grid floor `s//binsize*binsize`
  that a "bins relative to the region start" implementation would get wrong.
* The states-path argmax tie-break: Python's `max` over `((open,0),(nuc,1),(bg,2))` keeps
  the **first** maximal, which strict `>` reproduces.
* Columns 1-5 of `_accessible_regions.narrowPeak` (F227).

**Where it is not.** The digested signals diverge at *scale*, not in the last bits:

| track | oracle records | ours | shared starts | max abs dev | shared with different end |
|---|---|---|---|---|---|
| short | 336509 | 336453 | 336182 | **2.14** | 536 |
| mono | 451893 | 451887 | 451886 | 1.0 | 8 |
| di | 266801 | 266830 | 266437 | 0.139 | 555 |
| tri | 113800 | 113721 | 112826 | 0.861 | 1401 |

An earlier note recorded `max_abs_dev 1e-5` here. That was measured at shared bin
*midpoints* on a different corpus; on shared record starts the deviation is two orders of
magnitude larger. **Both measurements were real and they are not comparable**, which is
the second time this project has had two honest numbers for the same quantity and no
record of which corpus produced which. Hence the harness.

Since the EM inputs agree, the divergence is inside `generate_weight_mapping` /
`pileup_from_lr_hmmratac` / the digest writer's record segmentation -- and because the
HMM consumes the digest, a depth error of ~2 moves state boundaries, which is exactly the
observed 10 bp median boundary offset and the ~1900 regions present on one side only.

So the order of work is: **digest values first, boundaries second.** Fixing the states
path without fixing the digest would be rearranging the deck.

## F230 -- pnorm2 is a float32 subtraction, and four tests were pinning a dropped run

Continuing F229's "digest values first". The digest is a weighted endpoint sweep,
so the weights are `pnorm2` divided by the sum of four `pnorm2` values. Two real
defects came out of that.

### 1. `pre_z` starts at -10000, so the leading zero-depth run is real

`PileupV2.py:651` seeds `pre_z = -10000`. That is the whole point: at the first
endpoint `z` is still 0, so the comparison `z == pre_z` **fails** and the first
record is emitted. With `pre_z = 0` it succeeds, and the `c > 0` guard upstream
turns that into a silent no-op -- the leading zero-depth run disappears.

This was invisible on yeast500k only because its first fragment starts at
position 0, so `z` had already left 0 by the time `e != s` was first true. On a
fragment `(10, 110)` weighted 1.0 the two disagree outright:

```
pre_z = -10000 (upstream, correct)   [(10, 0.0), (110, 1.0)]
pre_z = 0     (what we did)           [(110, 1.0)]
```

**Four existing tests asserted the buggy answer** --
`the_hmm_pileup_sweeps_weighted_fragments`, `..._merges_adjacent_equal_runs`,
`..._keeps_zero_depth_gaps` and `..._drops_zero_depth_runs`. They were not
"tests that pass"; they were a bug reproduced four times and named. Fixed, and the
zero-weight case re-framed: a missing weight gives `+0` then `-0`, every
`z == pre_z` succeeds and one record survives with its end extended.

### 2. `pnorm2` subtracts in float32

`Prob.pnorm2` is annotated `float32_t`, and it is compiled (`Prob.cpython-312*.so`).
So `x - u` is a C float subtraction: **rounded to f32 before it is squared**.
Widening first and subtracting in double is algebraically identical and is not
bit-identical.

Calling the compiled function directly is the only reliable way to settle this,
because the same expression read as *source* is ambiguous. Measured agreement of
raw f32 bit patterns over the 2208 `(fraglen, class)` pairs of yeast500k:

| how `x - u` is done | agreement |
|---|---|
| f32 subtraction, rest f64 | **2208/2208** in the numpy-scalar harness, **10/12** on the sampled table pinned in `PNORM2_ORACLE` |
| f64 subtraction | 2081/2208 |

Those two numbers disagree, and the disagreement is itself the finding: the same
expression evaluated through numpy `float64` scalars and through Rust's own
arithmetic does not agree with the compiled oracle on the same inputs. `pow` vs
`mul`, `exp` vs `expf`, `sqrt` vs `sqrtf`, and float-vs-double intermediates were
each tried; none of them, alone, accounts for it. **pnorm2 is therefore improved
but not yet bit-exact, and F230 stays open.**

`pnorm2_open_gap_is_recorded` asserts the current count exactly, so closing this
cannot happen silently. `pnorm2_agrees_with_compiled_where_it_matters` guards the
regression that actually matters -- someone "simplifying" `x - u` back to double
costs ~25% of the agreement.

### What the weight rounding still is

Upstream stores `p_s / s` in a plain `dict`, so it is nominally a double, yet all
2208 weights come back exactly f32-representable. No model reproduces all of them:

| `s` summation | quotient | weights reproduced (of 2140 non-zero) |
|---|---|---|
| f32, step-rounded | f32 | **1085** |
| f32, one rounding at the end | f32 | 1024 |
| f64 | f32 | 858 |
| any | f64 | 205 -- and *only* class 0 |

The last row is the informative one: an f64 quotient is exactly right for 205
class-0 weights and for **zero** of the 1078 in classes 1-3, so the quotient is not
uniformly f64. Solving for the implied `s` per class spreads by 1.5e-7, more than
f32 rounding of the weight can explain, so `s` is not a single shared value in the
way the source reads. The residual concentrates where the component probabilities
are **denormal** f32 (`p_tri ~ 4e-40`), i.e. exactly where double rounding and
denormal handling diverge. Until that is pinned, the best-measured model is kept
rather than a tidier one.

### Effect

`mono` record counts now match the oracle **exactly** (451893), first records align
on all four signals, and Jaccard is unchanged at 0.6585 -- so this moved the digest
but not yet the regions. Order of work is still: digest values first, boundaries
second.

## F231 -- the square is `powf`, and the digested tracks are now byte-identical

F230 got `pnorm2` nearly right and then stalled: 2208/2208 in one harness, 10/12
in another, same expression. The resolution was to stop reading the `.py` and read
the generated C, which is committed next to it.

### What the C actually does

`Prob.py:69` reads

```python
ret = 1.0/sqrt(6.283185307179586 * v) * exp(-(x-u)**2 / (2.0 * v))
```

and `Prob.c` (line 14430) says:

```c
__pyx_t_7  = PyFloat_FromDouble((6.283185307179586 * __pyx_v_v));   /* sqrt arg, double */
__pyx_t_9  = (-powf((__pyx_v_x - __pyx_v_u), 2.0));                 /* <-- powf, NOT pow */
__pyx_t_10 = (2.0 * __pyx_v_v);                                     /* double */
__pyx_t_4  = PyFloat_FromDouble(exp(((double)__pyx_t_9) / __pyx_t_10));
__pyx_t_7  = __Pyx_PyNumber_Multiply_object_float(__pyx_t_6, __pyx_t_4);
__pyx_v_ret = __Pyx_PyFloat_AsFloat(__pyx_t_7);                     /* single rounding */
```

So `(x-u)**2` is **`powf` in single precision** -- both the subtraction and the power
round to f32 before the widening division. Squaring in f64 is algebraically the
same expression and is not the same function. That was the entire F230 gap.

### The three float widths, all different

`generate_weight_mapping` is the same story, and every part had to be read off the C:

| site | C | width |
|---|---|---|
| `p_s/p_m/p_d/p_t` | `pnorm2(...)` return assigned to `float` | f32 |
| `s = ((p_s+p_m)+p_d)+p_t` | `float __pyx_v_s` | f32, **rounded at every step** |
| `ret_mapping[i][fl] = p_i / s` | `PyFloat_FromDouble((__pyx_v_p_s / __pyx_v_s))` | f32 **division** (`divss`) |

The last row looked impossible: `float / float` in C is a double, so the stored
weight should be a general double -- yet all 2208 weights come back exactly
f32-representable. That is only consistent with a single-precision divide, which is
what GCC emits for `float` operands on SSE. This also explains why an f64 quotient
was right for 205 class-0 weights and **0** of 1078 others (F230): the tiny classes
lose the ulp to the rounding, the dominant one does not.

### Effect: the digest is exact

With those three widths and libm `powf`/`exp`/`sqrt`, agreement is **2208/2208 on
`pnorm2` and 2208/2208 on the weights**, and the whole digest is byte-identical:

```
short records oracle=336509 ours=336509 IDENTICAL=True
mono  records oracle=451893 ours=451893 IDENTICAL=True
di    records oracle=266801 ours=266801 IDENTICAL=True
tri   records oracle=113800 ours=113800 IDENTICAL=True

identical records: 1169003/1169003
```

`oracle/check_hmmratac.py` now treats digest identity as a hard gate rather than a
printed line, because a digest that differs moves every HMM state boundary and would
make the Jaccard number stop localising.

### A correction worth recording

For most of F229/F230 I was chasing a "target" of `[5893, 8115, 4789, 2030]` record
counts. Those numbers came from **my own Python re-implementation**, not from the
oracle's output; the oracle's real chrI `short` count is **5899**. So a long stretch
of that work was optimising toward a target I had manufactured. The digest is only
now demonstrably right, and the honest marker for that is the byte comparison, not a
count. Lesson recorded because it cost several turns: a self-written reimplementation
must never be used as the reference for another reimplementation.

### What is left

Jaccard is still 0.6587. With the digest exact, the gap is downstream: the states
path differs (`chrXVI 330-600` is `open` upstream and stays `open` here through 680),
so the divergence is in `hmm_predict` -- hmmlearn's forward-backward and the
`argmax` tie-break over it, not in the signals feeding it.

## F232 -- the bin reads the run *before* it, and the walk stops when the signal does

F231 left one thing open and it turned out to be two defects, both in
`extract_value_hmmr`'s two-pointer walk.

### The value is off by one run

`get_data_by_chr` returns **run end** positions (`add_chrom_data_PV` is fed the
`(end, value)` pairs `pileup_PV` produces), and the walk advances only `while p1 < p2`.
So the value reported for a bin at `p2` is the one in force **before** the bin, not
the one covering it. Against chrXVI's actual runs:

```
[0,8)=0.83172   [8,10)=0.85724   [10,11)=1.75845   [18,20)=11.15583   [20,21)=12.32470

bin at 10:  upstream 0.85724  (run [8,10))   containing-run lookup gives 1.75845
bin at 20:  upstream 11.15583 (run [18,20))  containing-run lookup gives 12.32470
```

We were reading the covering run. Fixing the cursor to stop on the first run whose
`end >= bp2` is what took the bin count from 607603 to **605941**, exactly upstream's.

### `StopIteration` is load-bearing

Upstream's loop is a bare `while True` with `except StopIteration: pass`. When a
signal track runs dry the **entire walk ends** -- that bin and every later one on the
chromosome are silently dropped. We clamped the cursor to the last run instead.

That truncation is not incidental, it is why upstream has ~99 fewer bins per
chromosome than a straightforward implementation:

```
chrom        oracle     ours  diff      (before the fix)
chrIV        73758    73853     95
chrXV        52135    52234     99
chrXVI       43881    43980     99
TOTAL       605941   607603   1662      (17 chromosomes x ~1 per region)
```

All four signals end at the same position (they are built from the same fragments),
which is why upstream's `assert nn == len(extracted_data[i])` still holds.

## F233 -- a states run must never span two chromosomes

With the bins fixed, `*_accessible_regions.narrowPeak` came out byte-identical but
`*_states.bed` still differed at the very last line of chrXVI:

```
oracle  chrXVI 947660 948020 open | chrXV 0 710 open
ours    chrXVI 947660   710 open | (chrXV 0-710 open lost)
```

`generate_states_path` branches on `chrname != prev_chrom_name`; the chromed variant
had no such branch, so its "same label, just extend" path ran across chromosomes.
Because the emission order is descending by name and every chromosome starts at bin 0,
the first bin of the next chromosome usually wins the *same* label as the last bin of
the previous one -- so the runs merged and the previous chromosome's end was
overwritten with a coordinate from the next.

`narrowPeak` survived it only by luck: the merged run still happened to span a valid
open region. That is exactly the kind of thing that passes a Jaccard check and fails
a byte compare.

### Effect: G13 closed

```
regions:          oracle 9286  ours 9286  shared 9286
Jaccard (region): 1.0000   requirement >= 0.98
Jaccard (base):   1.0000
unmatched regions: oracle-only 0, ours-only 0
summit shift:     max 0 bp over 9310 identically placed regions

digested short/mono/di/tri:  IDENTICAL
*_accessible_regions.narrowPeak:  BYTE-IDENTICAL
*_states.bed:                    BYTE-IDENTICAL

PASS
```

The CI step for G13 is now strict; it carried `continue-on-error` for exactly as long
as it was red.

### How the two were found

Not by reading the source -- F230 had already established that reading `BedGraph.py`
was not enough. The sequence that worked:

1. Compare the bin *positions* against the oracle's `cr_bins`. They diverged only at
   the end of each chromosome, ~99 bins apiece. Systematic per-region, not per-base.
2. Dump our candidate regions and the oracle's `Regions`. **Identical, 1523/1523** --
   which eliminated the whole peak-calling half of the command.
3. That left the walk. Dump our bins *and their values*, and read the oracle's `cr_data`.
   Same positions, different values, on the first two bins.

Step 3 is the one worth keeping: matching *positions* while *values* differ points
immediately at the value lookup, and no amount of staring at the source had suggested
an off-by-one-run.

## F234 -- F191 is not the Poisson function: 1220 points agree exactly

F191 is the last large class of golden failures (~208 cases, plus 14
`_control_lambda.bdg` files). Every one of them is the **5th decimal** of a
printed value:

```
_peaks.narrowPeak  8.92562  vs  8.92563
_control_lambda.bdg  15.74799  vs  15.74800
```

An earlier note recorded that a 1410-point sweep showed 125 one-ulp disagreements
and concluded they were "the rounding boundary, not a different result". That
conclusion was drawn without a harness, so it was worth re-testing properly.

### The Poisson path is exact

`pscore = -1 * poisson_cdf(observed, lambda, False, True)`
(`ScoreTrack.py:83`), and `poisson_cdf(..., log10=True)` ends in
`round((residue-lbd)/log(10), 5)` -- so the p-score is quantised to 1e-5 *before*
anything compares it. A disagreement of one unit in the 5th decimal is a byte
difference in `*_peaks.xls` even when the pre-rounding values agree to 1e-12.

`crates/macs-stats/tests/poisson_vs_oracle.rs` is a new differential harness
(ignored by default) over a capture from the compiled oracle. Over **1220 points**
-- `k` from 1 to 5000, `lambda` from 0.3 to 4000 geometrically -- agreement is
**1220/1220**, bit-for-bit on the returned f64. So `log10_poisson_cdf_q_large_lambda`,
`logspace_add` and `py_round` are all correct, and F191 is **not** in the Poisson
function. It is in `lambda`.

The same harness also settles a transcription question: a faithful Python
re-implementation of upstream's loop agrees with the compiled extension on 165/165
points, so the source really is readable here -- unlike `pnorm2` (F230/F231), where
it is not.

### One real defect fixed along the way

`coalesce_into` (`callpeak.rs`) divided the f32 array value by the denominator **in
f32**, and `denominator` itself was truncated to f32. Upstream's generated C:

```c
float __pyx_v_pre_v_c;
__pyx_v_pre_v_c = ((__pyx_v_ctrl_array_ptr[0]) / __pyx_v_denominator);
```

`ctrl_array_ptr` is a `cython.float*` but `denominator` is a **Python float**, so
the `/` is the Python operator: a **double** division, narrowed to f32 exactly once
on the assignment. Dividing in f32 double-rounds, and narrowing the denominator
first loses a second bit. Both are fixed.

This did **not** move the golden totals, because the failing cases are non-SPMR
where `denominator` is exactly `1.0` and the division is a no-op. It is recorded
because it is provably correct from the C, and it would have bitten every `--SPMR`
bedGraph.

### Where it actually is

With the Poisson path exonerated, the residual has to be in the lambda track. The
`_control_lambda.bdg` mismatches are all on `chrM` at **negative** coordinates
(`-4910`), i.e. the `F170` coordinate-shift path, and the lambda value entering
`%.5f` is an f32 sitting within one ulp of a `%.5f` boundary.

Next step, and it is the same technique that closed F232: capture upstream's
`ctrl_array` at full precision (it is only ever written at 5 decimals) and diff it
against ours run for run.

## F235 -- the weighted sweep applies starts before ends; and F191 is one f32 ulp in the lambda

Continuing F234 on the `_control_lambda.bdg` mismatches. Two results: one real defect
fixed, and the residual pinned to a single bit.

### Fixed: starts before ends at every position

`_pileup_sorted_weighted_as_list` (`PileupV2.py:278`) keeps `start_poss` and
`end_poss` in **separate sorted arrays** and folds them in with two pointers:

```python
while i_s < ls and start_ptr[i_s] == p:  z += start_w_ptr[i_s]
while i_e < le and end_ptr[i_e]   == p:  z -= end_w_ptr[i_e]
```

Our `sweep_weighted_to_track` merged the two into one list sorted by `(pos, value)`
and swept that. Sorting by value puts the **negative** end weights first at a tie --
the reverse of upstream. `z` is an f32 accumulator, so the two orders are not
interchangeable, and because the sweep coalesces on `scaled_z == pre_z` (an exact f32
equality) a last-bit difference moves the **breakpoints**, not just the printed value.
Now a two-pointer walk, starts first.

This did not move the golden totals (see below), but it is unambiguously correct from
the source and it removes a real order-of-operations divergence.

### The residual is exactly one ulp, and it is in the lambda

New dump: `MACS3_RS_DUMP_CTRL_LAMBDA` writes the exact f32 values entering `%.5f`,
which are otherwise unobservable -- the file only ever records five decimals.

For `sweep/gmini_mfrag_d1200_w180_ctrl_042_B`, the row `[chrM, -4910, -4909)`:

```
ours     15.74799538   bits 0x1.f7ef9401d5175p+3   -> %.5f 15.74800
golden                                             -> %.5f 15.74799
```

f32 ulp near 15.748 is `2^-20 = 9.54e-7`, and the `%.5f` boundary sits at
`15.747995`. So upstream's value is the adjacent f32, `15.74799443`. The difference
is **one bit of one multiply**, `scaled_z = z * scale_factor`, in the control's
windowed maximum over the `d`/`slocal`/`llocal` ladder.

What is *not* the cause, now excluded by measurement rather than by argument:

* `poisson_cdf` -- 1220/1220 exact (F234).
* `_control_lambda.bdg` is the **only** intermediate that differs; `_treat_pileup.bdg`
  is byte-identical for the same fixtures. So the treatment-side pileup, the model
  build and `lambda_bg` are all exact, and the divergence is confined to the control's
  scale ladder.
* `denominator` (F234) is a non-issue here: these runs are non-SPMR, so it is `1.0`.
* Accumulation order (above) is now correct but was not the cause.

Two loose ends worth writing down rather than rediscovering:

1. `ratio_treat2control` compares **weighted** totals, so it is not `1.0` even though
   both `#1 total fragments` log lines read `4794`. Our lambda value
   `15.74799538` is not `z * {1.0, 0.2, 0.02}` for any integer `z`, which is the
   evidence that the ratio is fractional.
2. `lambda_bg = treat_sum / gsize` reproduces the golden's first row
   (`14.28046`) exactly, so the floor is right; the ramp above it is the windowed
   max, as the log's "Range for calculating regional lambda is: 1000 bps and 10000
   bps" indicates.

Next step is unchanged from F234 and is now cheap: dump `ctrl_scale_s` and the
per-scale `z` for one fixture and compare the three products against the oracle's
control array at full precision.

## F236 -- F191's ladder: four of the five inputs are now confirmed, one is not

Continuing F235. The control lambda is `max over i of (window max_i) * factor_i`
for `i` in the `d`/`slocal`/`llocal` ladder. Instrumenting
`sweep/gmini_mfrag_d1200_w180_ctrl_042_B` (`MACS3_RS_DUMP_CTRL_LAMBDA` now reports the
whole ladder, not just the winner) gives:

```
d=143 slocal=1000 llocal=10000 ratio=0.5000003632044805
factors=[0.50000036, 0.07177755, 0.0071777552]
tsize_exact=143.555 control_sum=1376635 treat_sum=688318
```

### Confirmed correct, and worth writing down because each was a plausible suspect

* **`average_template_length` is f32 upstream.** `PairedEndTrack.py:244` is
  `self.average_template_length = cython.cast(cython.float, self.length) / self.total`
  -- an explicit cast to f32. So `ratio = 688318 / (9588 * f32(143.579...)) =
  0.5000003632044805`, **not** the clean `0.5`. Ours reproduces that exactly. Anyone
  "simplifying" this to f64 would break the control ladder.
* **`lambda_bg = treat_sum / gsize = 688318 / 48200 = 14.28046`**, which is the
  golden's first control-lambda row to five decimals. The floor is right.
* **`self.d` is 143, not 200.** `--extsize 200` is *ignored* in PE mode
  (`callpeak_cmd.py:165`, `options.d = options.tsize`), and `--nomodel` makes tsize
  the predicted `d`. So `ctrl_d_s = [self.d]` is `[143]` and our window matches.
* **`control_total = control.total * 2`** -- the log's "total fragments in treatment:
  4794" against 1200 file lines is the *count* sum, and the doubling is why.

### The one open input

`tsize_exact = 143.555` is the **per-line** mean fragment length; upstream's `tp.d`,
which becomes `options.tsize` (`callpeak_cmd.py:362`) and therefore the `float(self.d)`
in the wide-window factors, is **143.579** -- the weighted mean `length/total`. Both
print as `143.6`, which is exactly why this survived so long: the log cannot
distinguish them.

Changing `tsize_exact` to the weighted mean was tried and the golden suite did not
move (7079/7685, 22245/22456 unchanged), so the change is **not** justified by
evidence and was reverted rather than shipped. What is needed is upstream's `tp.d`
definition read directly, not a guess that happens to be neutral on 7685 cases.

### Why the residual is small

The failing row is `chrM -4910 -4909`, value `15.74799538` (bits `0x417bf7ca`)
against upstream's `15.74799`, i.e. `0x417bf7c9` or lower -- **3 ulps**, not 1.
Since none of the three factors divides our value to an integer depth
(`/0.50000036 = 31.49599`, `/0.07177755 = 219.388`, `/0.0071777552 = 2193.88`),
the reported value is not a single scale's `z * factor` for an integer `z`; the
window maxima themselves are involved. That is consistent with `over_two_pv_array`
merging the paired control's two projections and with the ladder's max being taken
on the scaled values.

Next step: read `tp.d` in upstream's FRAG parser and pin it against ours. That is a
single definition, and it is the last unconfirmed input on this path.

## F237 -- `tp.d` and `average_template_length` are different numbers, and both are f32

F236 left exactly one unconfirmed input on the control-lambda ladder: whether
`tsize_exact` should be the per-line mean (143.555, what we use) or the weighted mean
(143.579). It is the per-line mean. `Parser.py:1496`, in `FragParser`:

```python
self.d = cython.cast(cython.float, m) / i
```

`m` accumulates `right_pos - left_pos` once per **line** and `i` counts **lines**.
Counts in column 5 are not multiplied in. For
`sweep/gmini_mfrag_d1200_w180_ctrl_042` that is `f32(172266) / 1200 = 143.555`, which
is our `tsize_exact` exactly.

So upstream carries two different means and they are not interchangeable:

| quantity | definition | value on this fixture | used for |
|---|---|---|---|
| `options.tsize` = `tp.d` | `f32(unweighted length sum) / line count` | 143.555 | `self.d`, `ctrl_d_s[0]`, `options.d` |
| `treat.average_template_length` | `f32(weighted length) / total count` (`PairedEndTrack.py:244`) | 143.579 | `control_sum`, `ratio_treat2control` |

Both are printed as `143.6` by the oracle's own log, which is exactly why the
distinction survived: **no log line in the entire upstream output can tell them
apart.** A wrong choice here is invisible until a factor's last bits move.

This closes the last unconfirmed input from F236, and confirms that reverting the
weighted-mean experiment was right -- it was neutral on 7685 golden cases, which is
what a wrong-but-harmless change looks like, not what a fix looks like.

### Also re-verified while narrowing

* `over_two_pv_array` (`PileupV2.py:994`) takes the max **before** comparing
  positions, so when `a1_pos < a2_pos` it emits `a1_pos` carrying `max(v1, v2)` --
  a value belonging to a different coordinate. Upstream's merge is genuinely
  misaligned; `crates/macs-peaks/src/merge.rs` reproduces it (`val.push` precedes the
  position comparison) and drops the longer array's tail per F168. Nothing to change,
  but it is the kind of thing that looks like a bug to be "fixed" and must not be.
* `ctrl_total = control.total * 2` -- the log's "total fragments in treatment: 4794"
  against 1200 file lines is the count sum, and the doubling is real.

### Where the 3 ulps now stand

Every input to `max_i(window_max_i * factor_i)` is confirmed, and the merge that
combines the three scales is confirmed. The residual is therefore in the *values*
each scale's pileup produces -- but it is not any single scale's `z * factor` for an
integer `z`:

```
15.74799538 / 0.50000036  = 31.49599
15.74799538 / 0.07177755  = 219.388
15.74799538 / 0.0071777552 = 2193.88
```

A counted control's depth is an integer sum of counts, so one of those quotients
should be an integer. It is not, which means the surviving discrepancy is upstream of
the multiply -- in the accumulated depth itself, before scaling. That is the next
place to look, and it is a much narrower target than the five inputs just closed.

## F238 -- `tp.d` is f32, and that one ulp was 12 golden cases

F237 established that `options.tsize` is `tp.d`, the **unweighted** mean. It did not
settle the *width*. It is `f32`, and that is the whole residual.

`Parser.py:1496`:

```python
self.d = cython.cast(cython.float, m) / i
```

`cython.cast(cython.float, m)` makes the numerator a C float, and C promotes the `int`
denominator to float -- so this is `divss`, an **f32** division. `options.tsize` is
therefore a float32, and `PeakDetect.py:209`'s

```python
tmp_v = float(self.d)/self.lregion*self.ratio_treat2control
```

widens it without loss. We were keeping the mean in f64.

At the failing row on `sweep/gmini_mfrag_d1200_w180_ctrl_042_B` the control depth is
exactly **2194** -- a counted control's depth is an integer sum of counts, so it is
exact under either rounding and cannot be the source. The factor is:

```
f64 mean 143.555                   -> 0.007177755236625671  2194 * f = 0x417bf7ca -> 15.74800  (ours)
f32 mean 143.55499267578125       -> 0.007177754770964384  2194 * f = 0x417bf7c9 -> 15.74799  (upstream)
```

One ulp in one factor, found by solving for `(z, factor)` instead of guessing at it:
brute-forcing the integer depth against each candidate factor identified `z = 2194` on
the `llocal` scale and nothing else, which is what made the factor the only suspect.

### Effect

```
cases byte-identical : 7079 -> 7091   (+12)
output files matched : 22245 -> 22257 (+12)
_control_lambda.bdg mismatches: 14 -> 2
```

The two survivors are a **different** problem: they are coordinate differences, not
value ones (`chrM -300 -15` against `-300 -14`, and `-300 -9` against `-300 -8`), so
they are a run-boundary question in the control track, not arithmetic.

`crates/macs-peaks/tests/f238_control_factor_f32.rs` pins the invariant directly,
including the two bit patterns, so this cannot regress silently.

### How this was actually found

Four earlier hypotheses on this path were all plausible and all wrong: the Poisson
function (exonerated, 1220/1220), accumulation order (fixed anyway, F235), the
weighted mean (reverted, F237), and the SPMR denominator (irrelevant here). What
finally worked was refusing to guess and instead **solving the equation** -- enumerating
integer depths against each candidate factor until exactly one pair produced the
observed bits. Everything before that was narrowing by intuition; that was arithmetic.

The corollary is the useful habit: when a float disagrees by one ulp, ask what integer
times what reproduces the bits. It converts a search over causes into a search over
inputs.

## F239 -- where the remaining golden failures now are, measured rather than estimated

With F238 in, the golden corpus reads **7091/7685** cases and **22257/22456** files.
This is a census of what is left, so the next pass starts from data.

### 147 mismatching files across ~50 fixture/variant combinations

Every one of them is still the same *class*: the fifth significant digit of
`-log10(pvalue)`, in both `*_peaks.xls` and `*_peaks.narrowPeak` (and the
`*_summits.bed` that inherits it). No coordinate, ordering, count or exit-status
mismatch remains anywhere.

```
36.2831 vs 36.2832     5.97002 vs 5.97001     8.83045 vs 8.83046    21.3881 vs 21.388
```

Two exceptions, and they are a different problem:

* `gmini_mfrag_d400_w600_ctrl_384_B_control_lambda.bdg` -- **coordinate** only
  (`chrM -300 -15` against `-300 -14`): a run boundary in the control track, not a value.
* `gmini_mfrag_d400_w600_ctrl_384_*_summits.bed` -- follow-on from the above.

### The clusters split cleanly by code path

| fixture | control | mode | variants failing |
|---|---|---|---|
| `realistic` | yes | BED (SE) | 7 |
| `spikes_only` | yes | BED (SE) | 3 |
| `disjoint_chromosomes` | yes | BED (SE) | 1 |
| `gmini_mfrag_d400_w600_ctrl_384` | yes | FRAG (PE) | 14 |
| `gonechrom_mpe_d4000_w180_noc_131` | no | PE | 12 |
| `gtiny_mfrag_d400_w600_noc_142` | no | FRAG (PE) | 11 |

So there are **two** remaining causes, not one: a with-control SE path, a with-control
PE path, and a no-control PE path. F238 fixed the with-control **FRAG/PE** path; the
`gmini_mfrag_d400_w600_ctrl_384` survivors are the run-boundary case above, and the SE
and no-control PE paths are untouched by it.

### The no-control PE path is a width question again

`PeakDetect.__call_peaks_wo_control`:

```python
if self.PE_MODE:
    lambda_bg = treat_length / self.gsize
else:
    lambda_bg = float(d) * treat_total / self.gsize
```

and `lambda_bg: cython.float`, so the f64 expression is narrowed once. Ours computes
`(treat_sum as f64 / gsize) as f32`, which matches -- but `treat_sum` must be exactly
`treat_length`, and with `--shift` the extended lengths are what `treat.length()`
accumulates. Worth a direct check rather than an assumption.

### The SE with-control path is the same ulp class

For `realistic_call_summits_peak_4` the summit depth is 489 and the printed p-score is
5.97002 (upstream) against 5.97001 (ours). Solving the oracle's own `poisson_cdf` for
the lambda that produces each printed value:

```
lambda ~= 392.2   d(pscore)/d(lambda) ~= 0.111 per unit
1e-5 of p-score    <=>  ~9e-5 of lambda  <=>  ~2-3 f32 ulps at 392
```

Same signature as F238: a float32 lambda landing on the wrong side of a `%.5f`
boundary. On this path `self.d = options.extsize` is an argparse **int**, so there is
no `tp.d` width question -- the candidate is whichever intermediate is being kept in
f64 where upstream keeps it in f32.

The method that found F238 applies directly here: instrument the effective local
lambda at one failing summit, then solve for `(integer depth, factor)` against the
observed bits rather than reasoning about which line is at fault.

## F240 -- G0: duplicates are reachable after all, and the scale ladder provably is not

The porting plan's step 1 requires intermediates for "reads, duplicates, `d`, model
arrays, treatment pileup, control pileup, per-scale lambda, merged lambda, p-score
track, q-score table, candidate/merged intervals, summits". Two of those were
documented as unreachable. One of the two turned out to be reachable, and the way to
prove the other unreachable is worth keeping.

### Duplicates: the retained counts were already there, the rates were not

`dump_stages.py` said duplicates could not be captured because `filter_dup`'s rate is
recomputed inside compiled `FWTrack.filter_dup`. That is true of the *structure* and
irrelevant to the *value*: upstream **reports** it.

```
#1  total tags in treatment: 2500
#1  tags after filtering in treatment: 32
#1  Redundant rate of treatment: 0.99
```

The retained counts were already captured as the pre/post-filter totals; only the
rates were missing. Capturing them needs the log **in process** -- `record_stages.py`
only sees the output after `dump_stages.py` has exited, so the first attempt read a
file that did not exist yet. A `logging.FileHandler` on the root logger, attached
before `callpeak_cmd.run` and removed in `finally`, fixes it.

`duplicates` is now a recorded stage: 26 numeric leaves across the five stage fixtures,
wired into `macs-compare`'s `STAGE_ORDER` between `reads_post_filter` and `totals`.

A wrinkle worth recording: for `-f FRAG` only the totals appear, and that is correct --
upstream forces `--keep-dup all` for FRAG and runs no duplicate filtering at all. An
absent field that means "no filtering happened" must not be read as a missing capture.

### The control scale ladder: a recorded negative result

`ctrl_d_s` / `ctrl_scaling_factor_s` / `lambda_bg` are C attributes of
`CallerFromAlignments` and raise on attribute access. The obvious second route is that
they are *passed into* something, so wrapping that call would observe them.

They are not. `callpeak_cmd.py:265` is `peakdetect.call_peaks()` -- **no arguments** --
because `PeakDetect.__init__` builds the ladder as locals and keeps it on `self`. A
delegating proxy around the `PeakDetect` object records the call and finds nothing.

So the stage is recorded with an **empty** ladder, deliberately: a negative result in
the committed corpus is worth more than a silently absent one, and it is a regression
guard -- if a future upstream ever passes the ladder as an argument, the stage fills in
and the diff notices. The ladder remains diffable a different way: it is what
`stages_control_lambda.bdg` is written from, and that is captured.

### Model arrays

Still not capturable as an intermediate (`PeakModel` is opaque), but `*_model.r` -- the
file MACS3 writes from those same arrays -- is byte-compared by the golden layer. That
covers the content as an output rather than as an intermediate, which is weaker evidence
but is real evidence, and it is now stated that way instead of being called a gap.

## F241 -- the remaining clusters are not all ulp boundaries, and that matters

F239 split the 594 remaining cases into clusters. This checks whether they share a
cause, and the answer is no -- which changes what to do next.

### The no-control PE cluster is a systematic difference, not a rounding edge

For `sweep/gonechrom_mpe_d4000_w180_noc_131/default`, a single row differs:

```
golden  chr1 2877 3309 433 3090 180 213.755 32.3792 210.233
ours    chr1 2877 3309 433 3090 180 213.755 32.3793 210.233
```

Solving the oracle's `poisson_cdf` for the lambda behind each printed value, at
depth 180:

```
golden  pscore 32.3792  ->  lambda 63.83760
ours    pscore 32.3793  ->  lambda 63.83748
relative difference 1.9e-6   (~30 f32 ulps at 63.8, not 1)
```

F238 was one ulp. This is **thirty**. So this cluster is not the F238 class and cannot
be closed by hunting a single rounding site -- something is computed differently, by
about two parts in a million.

The obvious inputs do not explain it either, which is itself informative:

```
fragments 4000, treat_length 720000 (mean 180)
lambda_bg                       = 720000 / 12000 = 60.0      (exact)
PE ctrl_scale_s[0]              = 720000 / (1000*4000*2) = 0.09
```

The observed 63.84 is above `lambda_bg` and above `0.09 * anything integral`
(`63.8376 / 0.09 = 709.3`, not an integer), so the value is not a single scale's
`z * factor` here either -- the paired control's `over_two_pv_array` merge and the
windowed maximum are both in play, and neither is dumped for this path.

### What this means for prioritising

| cluster | signature | class |
|---|---|---|
| with-control FRAG/PE (F238) | 1 f32 ulp | fixed |
| with-control SE (`realistic`, `spikes_only`, `disjoint_chromosomes`) | ~2-3 ulps | same class as F238, different site |
| no-control PE (`gonechrom_mpe_*`, `gtiny_mfrag_*_noc_*`) | ~30 ulps | **different cause** |
| `gmini_mfrag_d400_w600_ctrl_384` | coordinate only | run boundary, not arithmetic |

So the with-control SE cluster is worth attacking next with the F238 technique, and
the no-control PE cluster needs the effective local lambda dumped for that path before
any numeric hypothesis is worth forming. Treating all 594 cases as one "last decimal"
problem would have sent the next pass after rounding sites that cannot be at fault.

## F242 -- PE without a control builds an all-zero ladder here, and "fixing" it made things worse

Continuing F241 on the no-control PE cluster. There is a real structural asymmetry
between our code and upstream's, and it is *not* the cause of the 30-ulp difference.

### The asymmetry

`PeakDetect.__call_peaks_wo_control` (`PeakDetect.py:340-346`) never builds the
`d`/`slocal`/`llocal` ladder. It builds a **single-scale** control, and in PE mode the
two sub-branches are not even the same shape:

```python
if not self.nolambda:
    if self.PE_MODE:
        ctrl_scale_s = [float(treat_length) / (self.lregion*treat_total*2),]
    else:
        ctrl_scale_s = [float(self.d) / self.lregion,]
    ctrl_d_s = [self.lregion,]
```

Our PE path builds the with-control ladder unconditionally. With no control,
`ratio = treat_sum / 0 = 0`, so every factor comes out **0.0**:

```
PE lambda scales: d=180 slocal=1000 llocal=10000 ratio=0.0 factors=[0.0, 0.0, 0.0]
```

An all-zero control track is structurally wrong, and it is surprising that the golden
suite barely notices -- which says the zeroed control simply never wins the local
maximum on these fixtures (the treatment's own d-scale dominates).

### Implementing the upstream form made the output worse

Replacing it with upstream's `[treat_length/(lregion*treat_total*2)]` and
`ctrl_d_s = [lregion]` moved the peak boundaries and the p-score *away* from the
oracle:

```
golden  chr1 2877 3309 433 3090 180 213.755 32.3792 210.233
after   chr1 2873 3311 439 3090 180 227.756 37.513 224.234
```

So the change was **reverted**. Two reasons, and the second is the important one:

1. My factor indexing was wrong -- I left the single factor at index 0 while
   `LambdaScales` reads `factors.get(2)`, so all three slots still received `0.0` and
   only the *window* changed (180/1000/10000 became 10000 alone). The comparison above
   is therefore not even a clean test of the hypothesis.
2. More fundamentally: a change that moves output away from the oracle is not evidence
   of anything until it is understood. The observation stands on its own -- our PE
   no-control ladder is not upstream's shape -- but the fix needs the window/factor
   placement pinned before it is worth attempting again, and that needs the effective
   local lambda dumped for this path.

Kept as a documented asymmetry rather than a half-applied change. The golden corpus is
back at 7091/7685 with the original single-row difference.

## F243 -- the SE ladder is instrumented, and its factors are now the suspect

F239 ranked the with-control SE cluster (`realistic`, `spikes_only`,
`disjoint_chromosomes`) as the same class as F238 and therefore the right next target.
F242 then demonstrated the prerequisite: guessing at a fix without the ladder printed
is guessing. So the F236 diagnostic was extended to the SE path.

### The SE ladder, printed

For `se_model/realistic/call_summits` (`-f BED --nomodel --extsize 200`, with a
control):

```
SE lambda scales: d=200 slocal=1000 llocal=10000
  factors=[0.09643395, 0.01928679, 0.0019286791]
  lambda_bg=0.0192  treat_scale=1.0
```

### Applying the F238 technique there

The failing row is summit depth 489 with p-score 5.97002 (upstream) against 5.97001
(ours), which the oracle's `poisson_cdf` places at lambda ~ 392.2. F238's method is to
enumerate integer depths against each candidate factor until exactly one pair
reproduces the observed bits.

**No pair reproduces either value.** Against all three factors:

```
392.2 / 0.09643395 = 4067.03      392.2 / 0.01928679 = 20335.16
392.2 / 0.0019286791 = 203351.61
```

and an exhaustive sweep of `f32(depth) * f32(factor)` over depths 1..400000 for all
three factors produces **neither** 5.97002 nor 5.97001.

Since the depth of a BED control track is an integer count, that is decisive: **our
factors are wrong**, not our depths and not the Poisson function. This is the F238
signature relocated -- there, `tsize_exact` needed to be f32; here, something in the SE
`ratio`/factor chain is being kept at a different width than upstream keeps it.

Note what this rules out cheaply. Had the lambda been far from any single scale's
`z * factor`, the natural reading would have been "the windowed maximum is wrong". The
sweep distinguishes those: a wrong maximum would still be some integer times some
factor, and none is. So the maximum is fine and the factor is not.

The next step is exact and small: dump `ratio_treat2control` for this fixture and
compare its f32 value against the printed `factors[0] = 0.09643395`, then do for
`d/sregion` and `d/lregion` what F238 did for `tsize_exact`.

## F244 -- retracting F243: the factors are right, and the differing column is not the p-score

F243 concluded "our factors are wrong" from an exhaustive sweep finding no
`(integer depth, factor)` pair that reproduces the printed value. **That conclusion was
wrong**, and it was wrong for two independent reasons.

### The factors are exactly right

Upstream's own log for `se_model/realistic/call_summits`:

```
#1  total tags in treatment: 192
#1  tags after filtering in treatment: 192
#1  total tags in control: 2000
#1  tags after filtering in control: 1991
```

Nine control duplicates are removed, so `ratio = 192/1991 = 0.09643395278754395`, and:

```
f32(ratio)                    = 0.09643395   <- our factors[0]
f32((200/1000)  * 192/1991)   = 0.01928679   <- our factors[1]
f32((200/10000) * 192/1991)   = 0.001928679  <- our factors[2]
```

All three match the printed ladder exactly. The sweep missed because it assumed
`ctrl.total` was the file's line count (2000) rather than the retained count (1991) --
the very thing the `duplicates` stage added in F240 now records.

### The differing column is fold_enrichment, not -log10(p)

The failing row is

```
chr20 30653 31141 489 30978 9 8.5323 5.97002 6.45005  realistic_call_summits_peak_4
```

and the `*_peaks.xls` header is `chr start end length abs_summit pileup
-log10(pvalue) fold_enrichment -log10(qvalue) name`. Mapping those columns:

| value | column |
|---|---|
| 489 | peak **length** |
| 30978 | abs_summit |
| 9 | pileup |
| 8.5323 | `-log10(pvalue)` -- identical on both sides |
| **5.97002** | **fold_enrichment** -- the column that differs |
| 6.45005 | `-log10(qvalue)` -- identical |

So the p-score and q-score were never wrong on this fixture; the whole difference is in
`fold_enrichment`.

The **column identity** is certain -- `PeakIO.py:797` writes the header
`... pileup, -log10(pvalue), fold_enrichment, -log10(qvalue) ...`, which is the file I
diffed. The **formula** is an inference and is flagged as one: `fold_enrichment`
appears in `PeakDetect.py` only inside two docstrings, and is computed in compiled
Cython (`CallPeakUnit`), so it cannot be read from source the way `pnorm2` could not be
in F230. Treating it as `pileup / local_lambda` gives

```
golden  9 / 5.97002 = 1.5075303
ours    9 / 5.97001 = 1.5075348     ~3e-6 relative
```

which is consistent with a small systematic lambda difference -- but confirming that
needs the peak-tuple construction in `CallPeakUnit.c`, not an assumption. It is recorded
here as a lead, not as a result.

I had read `489` as the pileup and `5.97002` as `-log10(p)` because the narrowPeak
column order differs from the xls one, and then built a numeric argument on a depth
489 that was never a depth. Every sweep that followed was answering the wrong question.

### What this does and does not change

It does not explain the residual. But it removes a false lead and corrects the record: F238's fix was
found by the right method (solve for `(depth, factor)` against observed bits), and F243
applied that method to a misidentified column, so its negative result said nothing about
the factors.

The habit worth keeping is narrower than I claimed: the technique is only as good as
the column mapping underneath it. `fold_enrichment` is derived from the scored peak rather than being a score in its
own right, so a fixture failing only there is not automatically evidence about the
Poisson path -- the two have to be separated before any sweep is run. Which is what
went wrong.

## F245 -- naming the column changes the census: 86 files differ in fold_enrichment *only*

F243 was retracted in F244 for reading a differing line positionally. That error was
made twice before it was caught, so the harness now names the column instead of the
line: `_named_diff` in `oracle/run_golden.py`, plus `oracle/check_peak_columns.py` as a
standalone checker. Re-running the corpus with that changes the picture materially.

```
realistic_call_summits_peaks.xls:34 columns=fold_enrichment
  golden=chr20 30653 31141 489 30978 9 8.5323 5.97002 6.45005 realistic_call_summits_peak_4
  got=   chr20 30653 31141 489 30978 9 8.5323 5.97001 6.45005 realistic_call_summits_peak_4
```

### The corrected census of all 147 mismatching files

| differing column(s) | files | what it means |
|---|---|---|
| `fold_enrichment` / `signalValue` only | **86** | p- and q-scores **identical** |
| `start` (plus knock-on `length`, `-log10(p)`, `qValue`, `peak`) | 34 | a real peak-boundary shift |
| `name` | 2 | peak index differs |
| unnamed (`.bdg` / `.bed`) | 25 | lambda bedGraphs, summits |

F239 reported 594 failing *cases* "all still the 5th significant digit of
`-log10(pvalue)`". That was wrong in the same way F243 was: it read line positions, not
columns. **On 86 of them the p-score and q-score match exactly** -- only
`fold_enrichment = f(treatment, local_lambda)` differs.

### Why this matters for the acceptance criteria

The release criteria separate *numeric* fidelity from *byte* fidelity:

* `p-score max abs error <= 1e-9` and `q-score <= 1e-6` -- **satisfied on those 86**;
  the printed values are identical to the last digit.
* byte-identical `*.xls` / `*_peaks.narrowPeak` -- **not** satisfied, because
  `fold_enrichment` is a byte.

So the remaining work splits cleanly:

1. `fold_enrichment` (86 files) is one derived quantity, `treatment / local_lambda`,
   computed in compiled Cython and therefore not readable from source. Closing it is one
   number, not 86.
2. The 34 `start` shifts are a *different and more serious* class than F239 suggested --
   F239 implied these were all arithmetic; they are coordinate differences, and a
   coordinate difference in `start` legitimately changes `length`, the summit and the
   q-score. That class deserves attention first on severity grounds, not on count.

### The process fix

The guard is deliberately cheap and total: it cannot be skipped, because the column
names are validated against upstream's own literals (`PeakIO.py:797` for narrowPeak,
`:968` for the xls reader) and any file whose header does not match is reported as a
failure rather than indexed by position. `check_peak_columns.py` reports the same
information for an ad-hoc pair of directories, and cross-references narrowPeak
`signalValue` to the xls `fold_enrichment` name so the two files are reported once.

The lesson generalises past this project: **a diff that reports positions should be
converted to names before any reasoning happens on top of it.** Two findings and one
retraction came from not doing that.

## F246 -- `fold_enrichment` is now a verified formula, and the residual is `summit_ctrl`

F244 flagged `fold_enrichment = pileup / local_lambda` as an **inference**, because
the value is assembled in compiled Cython. It was close but wrong in a detail that
matters. `CallPeakUnit.py:1385`:

```python
fold_change=(summit_treat + self.pseudocount) / (summit_ctrl + self.pseudocount),
```

with `pseudocount: cython.double = 1`. So there is a **pseudocount on both sides** --
`9/lambda` would be wrong, and for small values the difference is large.
`crates/macs-peaks/src/regions.rs:414` already matches this exactly.

### The residual is `summit_ctrl`, and it is quantified

For `se_model/realistic/call_summits` peak 4, `summit_treat = 9` on both sides and the
p-score and q-score are byte-identical, so the whole difference is the control term:

```
golden  (9+1)/(c+1) = 5.97002  ->  c = 0.6750362645351273
ours    (9+1)/(c+1) = 5.97001  ->  c = 0.6750390702863143
relative difference 4.2e-6
```

`summit_ctrl` is the control value at the summit -- the quantity the p-score uses as
its lambda. A 4.2e-6 relative difference is ~35 f32 ulps at 0.675: far too large for a
rounding site and far too small to be a wrong window.

### What is now pinned, and what is not

Ruled out by measurement, on this cluster:

* the ladder factors -- `f32(192/1991)` and its `lregion` siblings match the printed
  values bit-for-bit (F244);
* the Poisson path -- `poisson_cdf` is 1220/1220 exact (F234), and the printed p-score
  is identical anyway;
* the pseudocount formula -- verified above and already implemented identically;
* the peak boundary -- `start`, `end`, `length`, `abs_summit` and `pileup` all match;
  the 34 files that *do* differ in `start` are a separate cluster (F245).

Left: the value of the control track at the summit differs by ~4e-6 relative. Since the
factors are exact and a counted depth is an integer, the candidates are the **windowed
maximum selection** (which scale wins, and over which half-window) or the width of one
intermediate inside the window computation -- not the factor, not the depth, not the
output formatting.

That is a narrower target than any previous statement of this residual, and it is
derived from the file's own bytes plus two source lines rather than from a guess.

## F247 -- retracting F246's magnitude: this is a print boundary, not a 4e-6 error

F246 reported "a 4.2e-6 relative difference" in `summit_ctrl`. That number is wrong,
and it is wrong in the same way as F243 and for the same reason: **it inverts a rounded
print as though it were exact.**

`fold_enrichment` is written with `%.6g`, so `5.97002` denotes an *interval* of values,
not a point. Inverting it gives a bucket, and F246 took the bucket's centre as if it
were the value:

```
fold 5.97002  ->  c+1 in [1.675034862, 1.675037667]
fold 5.97001  ->  c+1 in [1.675037667, 1.675040473]
```

Doing the comparison the other way -- from our own computed value -- gives the real
picture:

```
ours  c+1 = f32(7 * f32(192/1991)) = 1.675037682056427   fold 5.970014948 -> %.6g 5.97001
f64   c+1 = 7 * (192/1991)          = 1.6750376695128075   fold 5.970014993 -> %.6g 5.97001
boundary (fold exactly 5.970015)    = 1.675037667409546
```

The two plausible computations of `summit_ctrl` differ from each other by **1.25e-8**
(1.9e-8 relative), and **both sit within 1.5e-8 of the `%.6g` rounding boundary**.
Upstream's value is on the other side of that boundary.

### So what the 86 `fold_enrichment` files actually are

`c+1 = 1.675037682` versus upstream's `<= 1.675037667`: the two sides agree to
roughly **1e-8**, and the printed digit flips only because the value sits that close to
a rounding boundary. That is:

* **far inside** any numeric tolerance the release criteria state;
* **still a byte difference**, so the `*.xls` / `*.narrowPeak` byte criterion is not met;
* and, importantly, **not** a systematic error of the kind F238 and F244 chased.

The remaining question is a genuine last-ulp one -- which of the two candidate roundings
upstream performs -- and it is worth at most the 86 files' single printed digit.

### The pattern, stated plainly

Three times now this project has built a numeric argument on a value read out of a
formatted output: F243 (read a line position as a column), F246 (inverted a `%.6g`
bucket as a point), and the F229 "target" counts that were this machine's own
reimplementation rather than the oracle's. Each produced a confident, wrong, specific
number.

The rule that follows, and which is worth more than any of the three findings:

> **Never reason numerically on a formatted value. Compute the value, or widen the
> formatting and read it back.**

The harness now names *columns* (F245), which removed the first failure mode. This one
needs the second guard: when inverting an output to recover an internal quantity, take
the interval and report the interval, not its midpoint.

## F248 -- the `fold_enrichment` cluster is one f32 ulp in the SE `ratio`

F247 localised the 86 `fold_enrichment` files to a value sitting within 1.5e-8 of a
`%.6g` boundary. That is enough to finish the job, because the candidate set is now
small enough to enumerate.

### The boundary, and what satisfies it

`fold = (treat + 1)/(ctrl + 1)` with `treat = 9`, so `5.97002` requires

```
ctrl + 1 <= 10 / 5.970015 = 1.675037667409546   ->   ctrl <= 0.675037667409546
```

Our `ctrl` is an **f32** (`regions.rs:61`), matching upstream's `cython.double` holding
an f32 product, so the comparison is at f32 granularity. Two candidate roundings:

```
ours   factor = f32(192/1991)          = 0.09643395245075226
       7 * factor = 0.675037682056427  -> ABOVE the boundary  -> 5.97001  (what we print)

       factor = 1 ulp below f32(192/1991) = 0.09643394500017166
       7 * factor = 0.6750376224517822  -> below the boundary -> 5.97002  (what upstream prints)
```

So the whole cluster reduces to: **upstream's SE control factor is exactly one f32 ulp
below `f32(192/1991)`, and ours is `f32(192/1991)`.** One ulp, in one number, explaining
86 files.

### Why this is not F238 again

F238 was `tsize_exact` held in f64 where upstream holds f32. Here the *inputs* are
already f32 on both sides -- `192` and `1991` are integers, the division is correctly
rounded, and `f32` of it is a single well-defined value. There is no width to
re-discover. What differs is the **quotient**: upstream's `ratio` is
`float(treat_sum)/control_sum` where `treat_sum = total * d = 38400` and
`control_sum = control.total * d = 398200`, i.e. `38400/398200`, not `192/1991`.

Those are the same real number, and a correctly-rounded f64 division returns the same
double for both -- but only if both operands are exact. They are. So `f32` of them is
one specific value, and it is the one we have. **Therefore upstream's factor is not
produced by narrowing that quotient at all**, which means the ratio reaching
`ctrl_scale_s` is not `38400/398200`.

The candidates worth testing next, in order of likelihood:

1. `ratio` computed from the **truncated** `control_sum`. F183/F184 already established
   that `control_sum` truncates on this path in some configurations; a truncation one
   unit high would push the ratio down by ~2.5e-6, which is ~20 f32 ulps -- too much --
   so this is likely wrong, but it is cheap to check.
2. `ctrl_scale_s[0]` being narrowed **twice**, or `ratio_treat2control` being carried as
   a `cython.float` from an earlier assignment.
3. The retained control count differing from `1991` in the *control* track used for the
   ladder versus the one used for the ratio -- the F240 `duplicates` stage records both
   and is the right instrument.

That is a short, checkable list rather than a guess, and each entry names the stage that
now holds the evidence.

## F249 -- checking F248's own premise: the ratio is identical, and all three scales give the same f32

F248 concluded that "upstream's `ratio` reaching `ctrl_scale_s` is not `38400/398200`".
That was an inference from the fact that our factor reproduces the printed value while
the observed fold needs a value one ulp lower. The premise is wrong, and it is worth
recording that the check was one line.

```
38400/398200 == 192/1991 (f64 bitwise): True     0.09643395278754395 both
f32 of each identical:                    True     0.09643395245075226 both
```

So `treat_sum/control_sum` and `treat.total/control.total` are the same double, and
narrowing gives the same f32. Our factor is upstream's factor.

### And the three scales are indistinguishable here

```
scale0 = ratio                  -> 0.09643395245075226
scale1 = (d/sregion) * ratio    -> 0.01928679086267948
scale2 = (d/lregion) * ratio    -> 0.0019286790629848838

boundary ctrl+1 = 1.675037667409546
  scale0: z=   7  ctrl+1 = 0.675037682056427   ABOVE
  scale1: z=  35  ctrl+1 = 0.675037682056427   ABOVE
  scale2: z= 350  ctrl+1 = 0.675037682056427   ABOVE
```

All three land on the **same** f32 and all sit a quarter of one f32 ulp above the
boundary. So which scale wins is irrelevant on this fixture, which removes F248's
candidate 3 as well.

### Where that leaves the 86 files

* the ratio is right (verified two ways);
* the factors are right (F244, and this note);
* the product `z * scale_factor` is the only remaining term, and it is **exactly one f32
  ulp** from what upstream must compute -- with both operands already f32 and exactly
  specified.

Two f32 values with exactly specified operands and one correctly-rounded multiply have
one answer. So upstream's `z` or its `scale_factor` is **not** the value we are using,
even though ours is the correct f32 of the correct ratio. The remaining possibility is
that `summit_ctrl` does not reach `fold_enrichment` as `z * scale_factor` at all -- that
`peak_content[summit_index][3]` is some other intermediate (a `max` against `d` or
`lambda_bg`, or a scaled value) whose f32 happens to differ by one ulp.

That is the same conclusion F248 reached by a longer route, with one candidate now
eliminated and the premise corrected. It needs the construction of `peak_content`
inside `CallPeakUnit` -- compiled, and not readable from source, so it must be captured
behaviourally the way `over_two_pv_array` was in F168 rather than read.

## F211 -- where the peak RSS actually goes, measured rather than argued

The criterion is `<=50%` of upstream, and the port was at ~1.6-2.1x upstream's RSS.
The previous note attributed it to "read count, not genome size" and prescribed "a
smaller read representation and a smaller `Run<f32>`". That was a hypothesis. Here is the
measurement, via an `MACS3_RS_RSS_TRACE` probe reading `/proc/self/statm` at three
points, on the documented benchmark (24 chromosomes x 200 kbp, 3 M treatment + 1.5 M
control, `--nomodel --extsize 200`, 1 thread):

```
rss[entry]             3.0 MB
rss[signals built]   162.3 MB     <- loading reads + building the two pileups
rss[after qtable]    188.8 MB     <- lambda, p-value and q-value tracks
peak                 223.8 MB
```

So of 224 MB:

| component | size | note |
|---|---|---|
| reads | ~36 MB | 4.5 M `Coord` (u64) positions, two `Vec<Coord>` per chromosome |
| treatment pileup track | ~80 MB | 3 M reads -> up to 6 M runs at 16 bytes each |
| control pileup track | ~40 MB | 1.5 M reads -> up to 3 M runs |
| lambda / p / q tracks | ~27 MB | already trimmed in narrow mode by `build_with_tracks` |

**The signal tracks are two thirds of the footprint, not the reads.** `Run<f32>` is
`{ end: u64, value: f32 }`, which `repr(Rust)` pads to **16 bytes**: 4 bytes of every run
are padding. Halving that to 8 by making the in-track position a `u32` (a chromosome is
never 4 Gbp, and upstream's own PV arrays already store positions as `i4`) would take the
tracks from ~120 MB to ~60 MB and the peak to roughly 130 MB -- parity with upstream,
still short of the 66 MB the criterion asks for.

**Also landed (F210):** `SortedPositions::finalize` now calls `shrink_to_fit`. `Vec`
doubles its allocation and neither `sort_unstable` nor duplicate filtering gives the
slack back, and the reads are alive for the whole pileup stage. Worth 8 MB
(232 -> 224), which is real but not where the problem is.

**Not attempted this turn:** the `Run<f32>` repack is a change to `macs-rle`'s core type
and every consumer, and the two-pass per-chromosome restructuring below it is larger
still. Both need their own turn with the full golden gate behind them; doing either in
hurry near a gate is how a 3 MB win becomes a 700-case regression.

## F209 -- the stage harness claimed stages it never captured

`dump_stages.py`'s own table listed `treatment_pileup`, `control_pileup`,
`lambda_merged`, `pvalue_track`, `qvalue_track` and `pvalue_stat`. None of them were in
`stages.json`: the harness had 7 stage keys, and a claim in a docstring is not evidence.

**Why they were missing.** They live in cdef attributes of `CallerFromAlignments`, an
immutable Cython extension type: `compute_pvalue` and `compute_qvalue` are not even
readable from Python, and `ScoreTrackII` cannot be rebound either (it is not a subclass
and has no `__dict__`). Intercepting them the way the harness intercepts `PeakDetect` does
not work.

**Fix.** The tracks are obtained the only other way there is -- by asking the engine to
persist them. `callpeak -B` writes `<name>_treat_pileup.bdg` and
`<name>_control_lambda.bdg`, which are the treatment pileup and the merged local lambda,
as the same bytes a user receives. `dump_stages.py` now passes `-B` and parses both.

**And the Rust side was extended to match**, so the comparison is real rather than
one-sided: `stagedump::record_tracks` emits `treat_pileup` in the same
`{chrom: [[start, end, value], ...]}` shape at `%.5f`, which is what
`bedGraphIO.write_bedGraph` writes.

Intervals are kept rather than just values, deliberately: a differential comparing only
values would miss a segmentation difference, which is exactly the F207 class of bug.

**Result**: `treat_pileup` is byte-identical to the oracle on all five recorded fixtures
(SE, three PE, and FRAG) -- 1086 numeric leaves, zero deviation:

```
  treat_pileup             5    1086       0      0.000e0      0.000e0  ok
```

**A bug this surfaced.** The FRAG fixture first differed by exactly 5000 on every
coordinate. That is not a library defect: a counted (`--format FRAG`) run is computed in
a shifted frame, `coord_shift = max(d, slocal, llocal) / 2` (F181), and *every writer
subtracts it again* before emitting a coordinate. The new dump was the one writer that
forgot, so it reported the whole FRAG pileup displaced. Fixed by subtracting `coord_shift`
in `record_tracks`, which is why the SE branch passes a literal `0` with a comment saying
so rather than leaving it looking like a placeholder.

**Still genuinely unreachable**, and now documented as such rather than claimed: the
control pileup (upstream's `-B` writes the control *lambda*, not the control pileup), the
p-score track (it exists only between two conversions inside an immutable type), and the
per-*scale* lambda arrays (locals in `call_peaks`; `xls_header` carries the window sizes
instead). `ctrl_pileup` is therefore not emitted by the Rust side either -- a stage with
no oracle counterpart would report "present only in ours" on every fixture forever, which
is the kind of permanent noise that trains people to ignore the comparator.

## F208 -- the whole fuzz layer had never been built, and it found a panic immediately

The L6 layer is wired into CI as six `cargo-fuzz` targets. `oracle/fuzz` was not:

* `fuzz/Cargo.toml` declared a `parse_frag` target whose source file **did not exist**;
* `fuzz_targets/pileup_pipeline.rs` passed `Vec<u32>` positions where
  `pileup_from_weighted_positions` takes `&[Coord]` (= `u64`).

Both are compile errors, so `cargo fuzz build` -- step one of the nightly job -- has been
failing for as long as the layer existed. No fuzz target had ever been executed.

**Fix, and what it found.** After building all six, the first run of the new
`parse_frag` target crashed within seconds on
`attempt to multiply with overflow` at `crates/macs-track/src/frag.rs:83`:

```rust
self.length += frag.end.saturating_sub(frag.start) * u64::from(count);
```

`PETrackII.add_loc` accumulates that product into a C `long`, so **upstream wraps**. Rust's
`*` panics under the overflow checks the debug and fuzz builds enable, and a FRAG record
can supply any coordinate and count. That is a direct violation of the "zero panics on
malformed input" acceptance criterion. Now computed with `wrapping_mul`/`wrapping_add`,
with `a_huge_count_on_a_long_fragment_does_not_overflow` pinning it.

### Four assertions that were wrong, not four bugs

Worth recording, because the temptation was to "fix" the library instead. Every one of
these fired almost immediately and **the library was right**:

| assertion | reality |
|---|---|
| `parse_frag`: `frag.end > frag.start` | `BEDPEParser` stores `atoi(f[1])`/`atoi(f[2])` verbatim and `PETrackI.add_loc` does not normalise, so an inverted fragment is reachable -- and upstream carries it into `fraglengths()` too. |
| `parse_bed`: `length >= 0` | Same class: the parser may accept `end < start` (upstream does), so the minus-strand length `pos - end` can be negative. |
| `parse_bedgraph`: `end >= start` | `bedGraphIO.add_loc` stores the pair verbatim. |
| `parse_bedgraph`: `value.is_finite()` | `atof("inf")` is infinity, and C's `atof` -- what upstream calls -- does the same. |

`pqtable_pipeline` needed more than a relaxation. It asserted that q-score is monotone in
p-score, and fired three separate ways, each for a different bad reason:

1. over a **sweep of off-grid p-scores**, where `qscore_or_zero` is an *exact-key* lookup:
   a p-score that never occurred in the histogram answers 0, so interleaving real keys
   with absent buckets looks non-monotone when it is only sparse;
2. with **negative p-scores**, which no caller can produce: `pvalue_stat` is `-log10(p)`
   and therefore `>= 0` by construction;
3. over the table's `entries()`, which are **not sorted** -- they are keyed on f32 *bits*,
   whose order is not the order of the values.

Monotonicity is a property of a *valid* AFDR table, and it is already checked where that
is meaningful: against oracle-derived tables in the crate's own tests (23 bit-exact
p->q tables). Asserting it against random bytes is not a stronger test, it is a wrong one,
so the target now restricts its input to the real domain and asserts what the acceptance
criteria actually require: nothing panics, every stored key and value is finite, and
lookups terminate across the whole range the criteria cover (including `p <= 1e-9`,
`q <= 1e-6`).

**Verified**: all six targets build under `cargo +nightly fuzz build` and complete 40000
executions each with no crash.

## F206 -- the HMMRATAC fragment down-sample used the wrong NumPy generator

**Root cause of the remaining G13 divergence.** Found by adding upstream's own
diagnostic log lines (`# Downsampled N fragments ...`, `# means: ...`, `# stddevs: ...`)
to our stderr so the EM could be compared against the oracle's log instead of inferred:

```
upstream  Downsampled 1040 fragments ...  means: 50 186.3 388.9 583  stddevs: 20 52.8 50.2 63.8
macs3-rs  Downsampled 1041 fragments ...  means: 50 182.6 381.4 580.8  stddevs: 20 53.7 46.6 68.0
```

One extra fragment in the down-sample moves the fitted nucleosome means by several bp.
Those means feed `generate_weight_mapping`, which builds the four digested signals, which
are the HMM's input features -- so the decoded state path and every accessible region
downstream is displaced. Nothing after the EM was wrong.

`PETrackI.sample_percent_copy` (`PairedEndTrack.py:653-659`) seeds a **modern** NumPy
generator:

```python
if seed >= 0:
    rs = np.random.default_rng(seed)
else:
    rs = np.random.default_rng()
rs_shuffle = rs.shuffle
```

`default_rng` is **PCG64**. This port's `macs_stats::NumpyRng` is the legacy
`RandomState` MT19937 stream -- correct, and bit-exact, for the `randsample`/`filterdup`
shuffles, but a completely unrelated algorithm from the one used here. Note the two APIs
are not even both "seeded NumPy": `FWTrack.sample_percent` really does use
`np.random.seed(seed)`, which is why only the *fragment* sampler had to change.

`macs_stats::pcg64` now implements PCG64 and is verified against the installed NumPy.
Four details were needed, and three of them are not guessable:

1. **The PCG64 seeding pair is big-endian.** `pcg128_t` takes `w[0]` as the *high*
   half: `initstate = w[0]<<64 | w[1]`. Little-endian pairing yields a plausible but
   different 128-bit constant and a wholly different stream.
2. **The output permutation is applied after the LCG step.** `pcg_output_xsl_rr_128_64`
   permutes the *new* state; permuting the pre-step state shifts the whole stream by one.
3. **`random_interval` draws `next_uint32`, not `next_uint64`, and the bound is
   inclusive.** `distributions.c`:

   ```c
   uint64_t random_interval(bitgen_t *bitgen_state, uint64_t max) {
     if (max == 0) return 0;
     mask = max;  mask |= mask >> 1; ... mask |= mask >> 32;
     if (max <= 0xffffffffUL)
       while ((value = (next_uint32(bitgen_state) & mask)) > max) ;
     else
       while ((value = (next_uint64(bitgen_state) & mask)) > max) ;
     return value;
   }
   ```

   A shuffle's bounds are element indices, so the 32-bit branch is *always* taken.
   Using the 64-bit draw, or Lemire's `bounded_lemire_uint64` from
   `random_bounded_uint64`, silently produces a different permutation.
4. **`next_uint32` is a 32-bit buffer over the 64-bit stream, low half first.**
   Recoverable exactly from Python because
   `rng.integers(0, 2**32, dtype=uint32)` hits `random_bounded_uint32`'s
   `rng == 0xFFFFFFFF` fast path and returns `next_uint32` unfiltered. Over seeds 0, 1,
   7, 999, 10151, 12345 the sequence is `step; lo32; hi32; step; lo32; hi32; ...`.

   The obvious candidates are all wrong: a second LCG step through PCG's
   `pcg_output_xsl_rr_128_32`, the low half alone, the high half alone, and one step per
   32-bit value. `random(dtype=float32)` is *not* a usable probe -- it exposes only the
   top 24 bits, and several wrong variants agree there.

**And a fifth bug, in the truncation, which turned out to be the whole residual.**
`sample_percent_copy` computes

```python
num = cython.cast(cython.uint, round(loc.shape[0] * percent, 5))
```

`round(x, 5)` returns a **float** rounded to five decimal places, and the Cython cast to
`uint` then *truncates*. So `11859 * 0.1` becomes `1185.9` and then **1185** -- not the
1186 that integer rounding gives. This port used `py_round(n * percent, 0)`.

**Verified**: the EM now agrees exactly with the oracle on every mean and stddev
(50/186.3/388.9/583 and 20/52.8/50.2/63.8, downsampled 1040), and the shuffle
permutation matches NumPy at n = 5, 12, 1186 and 11859 -- the last being the size
HMMRATAC actually uses, checked because `random_interval`'s rejection rate depends on the
bound.

**Effect on the gate**: base-level Jaccard on accessible regions rose from 0.958 to
**0.975** (940 of 946 regions). Still short of the 0.98 requirement.

### Where the residual G13 gap actually is

With the EM now identical, the digested signals were measured against upstream directly
(shared bin midpoints, since record segmentation is a separate issue now fixed by F207):

```
a_digested_short.bdg   shared  8560   max_abs_dev 1.0e-05
a_digested_mono.bdg    shared 18391   max_abs_dev 1.0e-05
a_digested_di.bdg      shared 20810   max_abs_dev 1.0e-05
a_digested_tri.bdg     shared 13441   max_abs_dev 1.0e-05
```

An independent NumPy re-implementation of `generate_weight_mapping` +
`pileup_from_LR_hmmratac` reproduces upstream to the same <=1e-5, so the weight mapping,
`pnorm2`'s float32 semantics and the pileup formula are all correct. The floor is float32
**accumulation order**: `PV.sort(order='p')` is numpy's introsort, which is *not* stable,
so endpoints sharing a position are summed in an order that depends on the sort's
internal pivoting, while `z += v` is sequential float32. This port uses a stable
`sort_by_key`, giving a different (equally valid) order and a different last bit.

Matching it would mean reproducing numpy's introsort for a structured key. That is
possible but disproportionate: the deviation is 1e-5 on a signal whose downstream effect
is to flip a handful of near-tied Viterbi decisions. Recorded as the open G13 residual
rather than papered over.

## F205 -- `hmmratac` decodes the whole genome; upstream decodes only candidate regions

Still open, and now the sole remaining cause of the G13 region divergence.

`hmmratac_cmd.py:434-470` builds `candidate_regions` by applying `--prescan-cutoff` to
the digested signals, then loops

```python
while candidate_regions.total != 0:
    cr = candidate_regions.pop(options.decoding_steps)
    [cr_bins, cr_data, cr_data_lengths] = extract_signals_from_regions(
        digested_atac_signals, cr, binsize=options.hmm_binsize, ...)
    prob_data = hmm_predict(cr_data, cr_data_lengths, hmm_model)
```

so only candidate regions are ever decoded. This port bins the whole genome instead. On
the ATAC fixture that is 78246 decoded bins against upstream's 3830, and the extra bins
are *not* harmless padding: they change the sequence the forward-backward pass runs over,
so the posteriors move (shared bins agree to ~1e-4 relative, e.g. upstream `0.00383`
against ours `0.00376`) and with them the accessible regions.

**Correction.** This diagnosis was wrong, and the measurements that disproved it are
worth recording. Our candidate-region construction already matches the oracle exactly
(330 candidate peaks -> expand by 1000 -> merge -> **17 regions**, same on both sides),
and our decoded bin count matches upstream's. The 78246-vs-3830 gap in
`a_bg.bdg` is not decode scope at all: `bedGraphIO` accumulates each bin through
`add_loc`, which merges runs of equal values into a single record, and `write_bedGraph`
then emits one line per *record*. Upstream's `chr1 510 1680 0.00000` is 170 merged bins
printed as one line. F204 fixed the columns and the precision but not this merging, so
`--save-likelihoods` still differs structurally from upstream's.

The real cause of the region divergence is F206.

## F199 -- `randsample` rejected every invocation without an explicit `-s`

Found by extending the base-invocation sweep in `oracle/audit_accepted_flags.py` to
more subcommands: getting a valid base for `randsample` failed on our side, which is
itself the bug.

`-s/--tsize` is **optional**. `randsample_cmd.py:87-89` only falls back to the
parser's own estimate when the user did not supply it:

```python
if not options.tsize:           # override tsize if user specified --tsize
    ttsize = tp.tsize()
    options.tsize = ttsize
```

so omitting `-s` is the normal path, not an error. This port required `tsize > 0` and
rejected the invocation:

```
$ macs3    randsample -i treat.bed -f BED -p 50 -o a.bed   # rc=0
$ macs3-rs randsample -i treat.bed -f BED -p 50 -o a.bed   # rc=1 --tsize must be > 0
```

**Fix.** `macs_io::detect_tsize` transcribes `Parser.tsize()` (`Parser.py:267-293`
for text, `Parser.py:1024-1070` for BAM). Two details are load-bearing:

* The text path accumulates into a `cython.int`, so `s/n` is **C integer division**. A
  mean read length of 149.1 reports `149`, not `149.1`.
* The BAM path accumulates into a `cython.double`, so its `s/n` is **float** division
  truncated on the cast. The same file therefore truncates in two different directions
  depending on the container format.

The scan is also bounded twice -- at most 10 *successful* lines and at most 10000 lines
read -- and seeks back to 0 afterwards, so detection never consumes a record.

**Verified**: `randsample` is now byte-identical for the auto-detected size, an explicit
`-s`, and `-n` (3 cases in `oracle/check_predictd_randsample.sh`).

## F200 -- `predictd` exited 1 where upstream exits 0

Two independent accept/reject divergences, both found while looking for a valid
`predictd` base.

**Not enough paired peaks is a warning, not an error.** `predictd_cmd.py:61-83` wraps
the whole model fit:

```python
try:
    peakmodel = PeakModel(treatment=treat, max_pairnum=MAX_PAIRNUM, opt=options)
    peakmodel.build()
    ...
except NotEnoughPairsException:
    warn("# Can't find enough pairs of symmetric peaks to build model!")
```

`macs3 predictd` therefore **exits 0 and writes no file**. This port surfaced the
condition as `MacsError::ModelBuildFailed` and exited 1, so every input upstream treats
as a successful empty run was a hard failure here.

**Paired-end mode is not "unsupported".** `predictd_cmd.py:54-60` short-circuits before
the model fit entirely:

```python
if options.PE_MODE:
    d = treat.average_template_length
    info("# Average insertion length of all pairs is %d bps" % d)
    return
```

The port rejected `--format BAMPE`/`BEDPE` outright. Both paths now warn/print and
return 0 with no output file. The insertion length is printed through `%d`, so a mean
of 98.307 is reported as `98`.

**Fix.** Added `MacsError::NotEnoughPairs { found, needed }` so callers can tell this
apart from a genuine model-build failure, and handled both short-circuits.

**Verified**: `predictd` byte-identical on a no-model dataset (exit 0, no files), a
paired-end dataset (exit 0, no files), and a dataset that does fit (see F201).

## F201 -- `*_model.r` was not byte-identical: `smooth` scaled the wrong side

`predictd` fitting a model exposes `*_model.r`, which the acceptance criteria list as a
byte-identical file. It was not: `p`, `m` and `xcorr` matched exactly, but `ycorr`
differed in 807 of 1200 elements by up to 1.6e-15 relative -- last-bit noise from a
different summation order.

`PeakModel.py:505-506` is:

```python
w = np.ones(window_len, 'd')
y = np.convolve(w/w.sum(), s, mode='valid')
```

The window is divided **first**, so every product is individually scaled and rounded,
and then summed:

```python
w = np.ones(window_len)
y = np.convolve(w / w.sum(), s, mode='valid')   # sum of (1/n) * s[k+j]
```

This port summed the window and divided once:

```rust
let mut acc = 0.0;
for j in 0..window_len { acc += s[k + j]; }
*yk = acc / window_len as f64;                   // sum of s[k+j], then / n
```

Those are different operations, not two spellings of one: on a 50-element window they
disagree in 31 of 50 outputs. The existing comment claimed forward `j` order was
required "to match `np.convolve`" -- the *order* was right, the *placement of the
division* was wrong, which is why a correct-looking loop still diverged.

Fixed by scaling each term before accumulating. The index set is unchanged (the flat
kernel is symmetric), and forward `j` order is confirmed correct against NumPy.

**Verified**: `predictd_model.R` byte-identical, and a regression test
(`smooth_flat_scales_each_term_before_summing`) asserts both the exact expectation and
that the old ordering is distinguishable, so the test cannot silently stop testing
anything.

## F197 -- `bdgpeakcall --cutoff-analysis` was accepted and wrote nothing

**Found by `oracle/audit_accepted_flags.py`**, the generalisation of how F195 was found.
That script sweeps every `store_true` flag in the auto-derived matrix, runs the
subcommand with and without it, and reports flags whose output is byte-identical. The
golden corpus cannot find this class of bug: `run_oracle.py` records only invocations
someone thought to try, so "this flag does nothing" is never checked.

`bdgpeakcall_cmd.py:44-56` is an **`if/else`**, and the peak call is in the `else`:

```python
if options.cutoff_analysis:
    cutoff_analysis_result = btrack.cutoff_analysis(...)
    if options.ofile:
        fhd = open(os.path.join(options.outdir, options.ofile), 'w')
    else:
        fhd = open(os.path.join(options.outdir,
                    "%s_l%d_g%d_cutoff_analysis.txt" % (options.oprefix,
                                                        options.minlen,
                                                        options.maxgap)), 'w')
    fhd.write(cutoff_analysis_result)
else:
    peaks = btrack.call_peaks(...)
```

So the port had to get three things right, and got none of them:

1. `--cutoff-analysis` **replaces** peak calling. It does not accompany it. We called
   peaks *and* would have written a report.
2. The report goes to `-o`'s path when `-o` was given -- **overwriting** the file the
   peaks would have used. Only the no-`-o` case uses the
   `<oprefix>_l<minlen>_g<maxgap>_cutoff_analysis.txt` name.
3. The sweep is data-derived: `minv = max(min_score, track.minvalue)`,
   `maxv = min(track.maxvalue, --cutoff-analysis-max)`, `s = (maxv-minv)/steps`, and the
   ladder is `np.arange(minv, maxv, s)` with `round(v, 3)`.

**Fix.** `macs_peaks::callpeak::bedgraph_cutoff_analysis` transcribes
`bedGraphTrackI.cutoff_analysis` (`BedGraph.py:1262-1377`), reusing the same chunk walk
as the callpeak version, and `bdgpeakcall` branches on the flag.

**Verified** against the oracle on a bedGraph with a real score spread:

```
$ macs3  bdgpeakcall -i vary.bdg -c 5.0 -l 200 -g 30 --cutoff-analysis -o rep.txt
$ macs3-rs bdgpeakcall ... (same arguments)
rows: upstream 100 ours 100
IDENTICAL
```

Both output paths (`-o` and `--o-prefix`) are byte-identical, as are the plain and
`--call-summits` narrowPeak runs.

## F198 -- `bdgpeakcall -o NAME` wrote the wrong track and peak names

Same audit, next line of the same file (`bdgpeakcall_cmd.py:60`):

```python
if options.ofile:
    options.oprefix = options.ofile
```

`-o` **overwrites the prefix**, and that prefix is what `write_to_narrowPeak` uses for
both the `track` line and every peak name. So `-o peaks.narrowPeak` yields
`name="peaks.narrowPeak"` and `peaks.narrowPeak_narrowPeak1`.

This port kept `oprefix = "bdgpeakcall"`, producing a track line called `bdgpeakcall`
inside `peaks.narrowPeak`. The peak coordinates, scores and summit offsets were all
correct -- the entire diff was the two name fields, which is exactly the kind of
divergence a coordinate-only differential misses.

**Verified**: identical for the plain, `--call-summits` and `--o-prefix` forms.

## F195 -- `--cutoff-analysis` was accepted and silently did nothing

**Found by auditing the stage corpus**, not by a corpus diff: `oracle/record_stages.py`
writes `NAME_cutoff_analysis.txt`, and the port produced no such file at all.

The flag was in `oracle/flag_matrix.tsv` and parsed, so `macs3-rs callpeak
--cutoff-analysis` exited 0 -- while `macs3` wrote a 5-column report. That is an
**accept/accept** divergence, which is the same class as F172: the invocation is
accepted by both, but the observable result differs. Nothing in the golden corpus
caught it because `run_oracle.py` never passed the flag.

**Fix.** Transcribed `CallerFromAlignments.__cal_pvalue_qvalue_table`
(`CallPeakUnit.py:875-1023`):

* the ladder is `np.arange(0.3, 10.0, 0.3)` with `round(x, 5)`, descending. `round(x, 5)`
  rounds to **five decimals**, not to an integer. Rounding to integers collapses the
  ladder onto `{0, 1, ..., 10}` -- still 33 values, still descending, so it looks right
  and compares the wrong cutoffs. That was the first version's bug, and
  `cutoff_ladder_is_0_3_to_9_9_in_steps_of_0_3` now pins it.
* the per-cutoff scan runs on the **same p-score track** the pipeline already builds,
  folded in during the q-table pass, so no stage is computed twice.
* the cutoffs are then seeded into the AFDR histogram as **zero-length buckets**
  (`pscore_stat[cutoff] = 0`, line 978) *before* the q walk. This is the load-bearing
  part: `N = sum(values)` is unchanged, but `unique_values` gains a key, so the q-table
  acquires an entry and every q-score in the run shifts. A run with the flag therefore
  produces a legitimately different q-table from the same input -- which is why this
  cannot be implemented as a post-hoc report.
* `pos_array[above_cutoff - 1]` **wraps to the last element** in NumPy when the first
  above-cutoff index is `0`. Reproduced; otherwise the first chunk of such a peak starts
  at the wrong coordinate and `lpeaks` is wrong.

**Verified.** `oracle/check_cutoff_analysis.sh` compares our output against upstream's
recorded files: **5/5 byte-identical**, including the 34-row counted-FRAG case and the
single-header-line SE case where nothing clears even the 0.3 cutoff.

## F194 -- `Interval::overlaps` disagreed with `Interval::overlap_len`

Found by a `proptest` invariant (`overlap_len_is_bounded_by_both_inputs`) asserting
`overlaps(o) == (overlap_len(o) > 0)`, which the two methods did not satisfy.

```rust
pub fn overlaps(&self, other: &Interval) -> bool {
    self.start < other.end && other.start < self.end              // says true
}
pub fn overlap_len(&self, other: &Interval) -> Len {
    min(self.end, other.end).saturating_sub(max(self.start, other.start))  // says 0
}
```

For an **empty** interval the strict inequalities still hold -- `[781, 781)` against
`[0, 782)` gives `781 < 782 && 0 < 781` -- so `overlaps` returned `true` while
`overlap_len` returned `0`. The doc comment says "True when the two intervals share at
least one base", and an empty interval shares none, so the implementation contradicted
its own contract.

**Fix.** Guard on `!is_empty()` for both operands, so the two methods agree. Note the
distinction from `gap`, which is deliberately *looser*: `gap == 0` means "overlapping
**or touching**", and that is the test the peak walk uses to decide whether to merge.
`overlaps` means genuinely shared bases.

**Verified**: 563 tests pass; the full corpus is unchanged at 7079/7685 with 0 exit
mismatches, confirming no call site depended on the lenient behaviour.

## F193 -- the Savitzky-Golay kernel was a closed form, and only the SVD is right

**The dominant remaining `--call-summits` residual**, found by cross-checking
`sg::maxima` against NumPy on 300 random pileup-shaped signals rather than by a
corpus diff.

`SignalProcessing.savitzky_golay_order2_deriv1` does **not** use the closed-form
Savitzky-Golay weights:

```python
b = np.array([[1, k, k**2] for k in range(-half_window, half_window+1)], dtype='i8')
m = np.linalg.pinv(b)[1]
```

Row 1 of the **pseudo-inverse of an integer design matrix**, computed through an SVD.
Comparing the resulting `maxima` index sets against upstream:

| coefficient source | mismatches / 300 |
|---|---|
| `np.linalg.pinv(b)[1]` | **0** |
| closed form `3k / (h(h+1)(2h+1))` | 113 |
| `b.T @ inv(b @ b.T)` (mathematically exact) | 299 |

The *mathematically exact* alternative is worse than the closed form, which is the
whole point: `maxima` rounds the smoothed derivative to 16 decimals and then takes
its sign, so what is being decided is the last bit. The derivative of a flat-topped
pileup stretch is ~1e-16; `pinv` puts it at -4.9e-17 (which `round(16)` sends to
`0.0`, i.e. sign 0, i.e. **one broad summit**) and the closed form at +1.3e-16
(which survives as `1e-16`, i.e. sign +1, i.e. the flat stretch **splits into two
summits**). No mathematically-equivalent formula reproduces that, because the answer
is determined by the SVD's rounding.

**F141 had already found the closed form was wrong and "fixed" it by hardcoding
`pinv` for window size 179 only**, keeping the closed form for every other size --
on the reasoning that `smoothlen = min_length = d` takes only a few values and 180
is the most common. That is wrong for every other `d`: the corpus's `--extsize`/`d`
values map to windows all over the range, and those runs kept the broken closed
form. A special case for the value that happens to be most frequent is not a fix.

**Fix.** `crates/macs-peaks/src/sg_coeffs.tsv` now carries `np.linalg.pinv(b)[1]`
verbatim for every odd window size 3..=513, generated by
`oracle/gen_sg_coeffs.py` (committed, so the table is reproducible and auditable
rather than a wall of unexplained literals) and embedded with `include_str!`.
`deriv1_coefficients` looks the window up and the closed form is only a fallback for
a window the table does not cover -- which `callpeak` cannot request, since
`maxima` maps `d` to `d//2*2+1 <= 513` for any `d <= 512`.

**Verified.** The numpy cross-check drops from 117/300 mismatches to **4/300**, and
the corpus from 7071 to **7079** byte-identical cases (22221 -> 22245 files) with
**0 exit-status mismatches** -- no regression.

The remaining 4/300 are the same class: `np.convolve` accumulates through BLAS
`ddot`, whose summation order depends on the build's vector width and unrolling, so
the last bit of a near-zero derivative is not reproducible in portable Rust. Where
that flips a result the consequence is *bounded* rather than wrong: `maxima`
returning empty sends the caller down upstream's own failsafe
(`__close_peak_wo_subpeaks`), which is a legitimate upstream path.

## F192 -- `d` counted reversed BED records (`col3 < col2`)

**Found by** writing the L2 `proptest` invariants, not by the corpus. Upstream's
`GenericParser.d()` ([`Parser.py:279-293`]) samples tag lengths like this:

```python
while n < 10 and m < 10000:
    m += 1
    this_taglength = self.tlen_parse_line(thisline)
    if this_taglength > 0:
        s += this_taglength
        n += 1
```

and `BEDParser.tlen_parse_line` ([`Parser.py:432-439`]) is

```python
thisfields = thisline.split(b'\t')
return atoi(thisfields[2]) - atoi(thisfields[1])
```

The difference is **signed**, and it is taken between the two coordinate columns
**exactly as written** -- there is no strand handling at all in `tlen_parse_line`.
So a record with `col3 < col2` yields a negative length and is skipped by the
`> 0` gate.

This port instead used `(rec.end - rec.pos).unsigned_abs()` ("a minus-strand line
stores `end` as its 5' coordinate, so the length is `|end - pos|`"). That is wrong
in two ways at once, and they happen to cancel on ordinary input, which is why no
fixture caught it:

* for a plus record `col3 > col2` both forms agree;
* for a minus record `rec.pos` is already column 3 and `rec.end` is column 2, so
  `rec.end - rec.pos` is `-(col3 - col2)` and `unsigned_abs` restores the
  magnitude;
* the divergence needs a record with `col3 < col2` **and** a large magnitude.

**Repro.** Four records on `chr1`: three of `col3 - col2 = 100` and one reversed
`col3 - col2 = -800`.

```
$ macs3 callpeak ... -t t.bed -f BED --nomodel --extsize 200
# tag size is determined as 100 bps        <- upstream

$ macs3-rs callpeak ... (before the fix)
# tag size is determined as 275 bps        <- 1100 / 4 instead of 300 / 3
```

`d` is not just a header line: `options.tsize` becomes `maxgap` (F188) and the
minimum peak length, so a whole run shifts.

**Fix.** Take the signed difference of the columns as written and keep the `> 0`
gate. Since `rec.pos` is column 3 for a minus-strand record and column 2 otherwise,
`col3 - col2` is `pos - end` on minus and `end - pos` on plus:

```rust
let tlen = match rec.strand {
    Strand::Minus => rec.pos - rec.end,
    _ => rec.end - rec.pos,
};
if tlen > 0 { sum += tlen as f64; n += 1; }
```

**Verified**: upstream and this port now both report `100`; a second fixture whose
first ten tags are `100..109` confirms the cap at ten and the `cython.cast(cython.int,
s/n)` truncation (upstream `104`, this port `104`). Five unit regressions in
`crates/macs-cli/src/commands/callpeak.rs` (`tsize_tests`) cover the reversed
record, the minus strand, the ten-tag cap with truncation, and the `-1` sentinel.
The full corpus is unchanged at 7071/7685, i.e. the fix removes a real divergence
without disturbing any matching case.

**Lesson recorded**: the recorded corpus only covers well-formed input, so a
differential gate cannot be the only defence. This is the second bug found by an
invariant rather than by a diff (the first was F190's missing mirror), and both
were in code paths where the wrong answer was *plausible* rather than absurd.

## F191 -- OPEN: the compiled `poisson_cdf` disagrees in the 5th decimal

**Status: unresolved, and the largest remaining item.** 208 of the 614 remaining
corpus failures are this: `-log10(pvalue)` differs from upstream by exactly
`1e-5` while every coordinate, summit, pileup, length and fold change matches.
Because the p-score is quantised to `1e-5`, a disagreement there is one unit in the
last printed decimal.

Measured against the pinned oracle with a 1410-point sweep covering
`-log10(p)` from 0 to 1000:

| transcription | mismatches / 1410 |
|---|---|
| this port (Rust) | **125** |
| a literal Python transcription of `Prob.py` | 179 |
| `m` starting at `k` instead of `k+1` | 1187 |
| `m` incremented after computing `logy` | 1397 |

Ruled out by measurement:

* the convergence threshold of `log10_poisson_cdf_Q_large_lambda` -- `1e-5` is
  best; `2e-5`, `5e-5`, `1e-4`, `2e-4`, `5e-4` are monotonically worse (78 -> 534
  -> 1277 -> 2486 -> 3890 on the `k = 1` sweep);
* the loop structure, `m`'s initial value, and the `sum_ln_m` range (the C and the
  `.py` agree with each other and with this port on all three);
* `logspace_add`: `math.log1p` vs `math.log(1 + exp .)`, and the `d < -745`
  underflow guard;
* `Prob.c` and `Prob.py` are *mutually consistent* (`c_constants = {0.5, 1.0,
  1e-5, 1e-10}`, `fabs(pre_residue - residue) < 1e-5`, `return round((residue-lbd)
  / log(10), 5)`), yet no faithful transcription reproduces the compiled
  extension's rounding. The unrounded values differ by ~1e-7 in the last place --
  far more than `O(sqrt(k) * eps)` accumulation error over `k` terms, so the
  compiled module was built from a `Prob.py` that differs from both artefacts in
  the tree.

Reproduction: `MACS3/Signal/Prob.c`'s `log10_poisson_cdf_Q_large_lambda` and
`logspace_add` next to the `.so`'s behaviour for a fixed `(k, lambda)`.

## Current measured state

`oracle/run_golden.sh`, release build, whole recorded corpus:

| scope | start of this round | now |
|---|---|---|
| `default` variant, cases byte-identical | 344 / 425 | **416 / 425** |
| all 7685 recorded cases byte-identical | 5660 / 7685 | **7071 / 7685** |
| output files matched | 18169 / 22456 | **22221 / 22456** |
| exit-status mismatches | 56 | **0** |
| workspace tests | 496 | **579** |

Every one of the 7685 recorded invocations reproduces upstream's exit status, and
every rejection path writes **no output file**.

Remaining failures by fixture family (across all variants):

| family | failing cases |
|---|---|
| `*_mfrag` (FRAG, counted) | 96 |
| `*_mpe` / `*_mse` / others | 78 |
| p-score 5th-decimal class (F191) | 208 of the 614 |

The paired-end families are closed. What is left is 208 cases of one float
(F191) plus 96 FRAG cases and 8 others spread across the run-boundary and
subpeak-naming edges.

### F250 (negative result): the `summit_ctrl` behavioural capture that F249 asks for is NOT reachable.
    F249 concluded that `peak_content[summit_index][3]` is a cdf local inside the compiled
    `__chrom_call_peak_using_certain_criteria`, so it must be captured rather than read. It cannot be:
      1. `CallPeakUnit.call_peaks` is a cdef method, so `self` inside it is the *real* `PeakDetect`.
         An object proxy around the factory is therefore never consulted (F240's `lambda_ladder`
         empty result rests on this, so that negative result is also non-authoritative).
      2. `peaks.add(chrom, start, end, summit=, peak_score=, pileup=, pscore=, fold_change=, qscore=)`
         IS a plain `def` and does receive upstream's `fold_change` as a full f64, with `pileup` as
         the summit treatment value -- so `summit_ctrl = (pileup + 1) / fold_change - 1` would recover
         it exactly. But `MACS3.IO.PeakIO.PeakIO` is an *immutable cdef class*:
             TypeError: cannot set 'add' attribute of immutable type 'MACS3.IO.PeakIO.PeakIO'
         Class-level rebinding is refused, so `add` cannot be intercepted either.
    Conclusion: the 86-file `fold_enrichment` residual (one f32 ulp, ours a quarter ulp above the
    `%.6g` boundary) has no observational handle inside the oracle. It stays open, and any future
    claim about it must not rest on inverting the formatted column. Our ratio is confirmed correct
    (`38400/398200 == 192/1991` bitwise in f64) and all three ladder scales give the same f32.

### F251: the `peak_content` capture is unreachable, and for a more specific reason than F250.
    `oracle/dump_peak_content.py` existed and claimed to do it. It never worked: it wrote a
    header-only TSV and reported success. Two separate defects, both now fixed/explained:
    - F251a: it rebound `callpeak_cmd.CallerFromAlignments`, a name `callpeak_cmd` does not
      mention, so the `Spy` subclass was never constructed. The real sites are
      `MACS3/Signal/PeakDetect.py:234` (with control) and `:350` (without).
    - F251b: `MACS3.Signal.PeakDetect` is a compiled `.so`, and Cython caches a module-level
      `from X import Y` in a *C* global at module init. Rebinding `PeakDetect.CallerFromAlignments`
      after the fact reads back as `Spy` yet the compiled code still builds the base class --
      verified directly. Rebinding `CallPeakUnit.CallerFromAlignments` *before* `PeakDetect` is
      first imported does work: `Spy.__init__` now provably fires once per region.
    With the subclass constructed, `_dump` still never runs, because:
        @cython.cclass
        class CallerFromAlignments:            # CallPeakUnit.py:351-352
        def __close_peak_wo_subpeaks(...)      # CallPeakUnit.py:1313
        getattr(C, '_CallerFromAlignments__close_peak_wo_subpeaks', None)  ->  None
    Cython gives *private* methods of a cdef class C-level names, so the mangled attribute does
    not exist in the type dict at all: it is neither readable nor overridable, and Python
    subclass dispatch never reaches it.
    So F250's conclusion holds, now with a mechanism and two independent confirmations
    (`PeakIO.add` immutable, `__close_peak_wo_subpeaks` C-level-private). The 86-file
    `fold_enrichment` residual has no handle inside the oracle. The harness now exits
    non-zero on an empty capture instead of reporting success -- a silently-empty
    intermediate capture is worse than none, because it reads as evidence.

### F252: `summit_ctrl` is unobservable inside the oracle -- third and final confirmation.
    F250 (immutable `PeakIO.add`) and F251 (C-level-private `__close_peak_wo_subpeaks`)
    both closed off the *input* side. This closes the *output* side, which is the one I
    actually wanted, and it does so by finding the values genuinely reachable:
      - `PeakDetect.call_peaks` is a public `def`, so wrapping `callpeak_cmd.PeakDetect`
        (an ordinary Python module -- the `dump_stages.py` trick) really does intercept
        the real call. Verified: 13 `PeakContent` objects captured on `se_model/realistic`.
      - `PeakIO.peaks` is `cython.declare(dict, visibility="public")` (`PeakIO.py:236`), and
        its live value is a dict of `{chrom: [PeakContent]}` -- a dict of LISTS, which is
        itself a trap: assuming dict-of-dicts raises `'list' object has no attribute 'items'`.
      - Subclassing `CallerFromAlignments` does NOT work for the public `call_peaks`
        either (F251c): Cython calls `self.call_peaks(...)` through a static vtable slot,
        so a Python override never dispatches -- `Spy.__init__` fires, zero rows.
    But `PeakContent` is `@cython.cclass` with no public surface whatsoever:
        dir(PeakContent)   ->  []
        pk.pileup          ->  AttributeError: 'PeakContent' object has no attribute 'pileup'
    So the tuple values are unreachable, and the derived
        summit_ctrl = (pileup + 1) / fc - 1
    cannot be recovered. `oracle/dump_peak_scores.py` keeps the working interception and
    records what remains observable -- the peak count -- and reports the opacity explicitly
    instead of writing a silently-empty file.
    Consequence for the plan: the 86-file `fold_enrichment` residual (one f32 ulp, ours a
    quarter ulp above the `%.6g` boundary, ratio and all three ladder scales verified
    identical) is NOT fixable from oracle evidence. It is recorded as an accepted
    unresolvable-without-upstream-patch difference, not as an open bug, and no further
    effort should be spent on it. The 34-file `start` shift class is directly observable
    from `_peaks.narrowPeak` and remains fully actionable.

### F253: a real off-by-one in the control-lambda track -- the first *observable* peak-adjacent bug
    found since the golden reporting fix. Reproducible on 2 files.

    Fixture `sweep/gmini_mfrag_d400_w180_ctrl_357` variant `B`
    (`callpeak -g 48200 -f FRAG --nomodel --extsize 200 -B`, treatment and control both 1597
    fragments, mean frag 144.2 bp, so no scaling is applied). `treat_pileup.bdg`,
    `_peaks.narrowPeak` and `_summits.bed` are byte-identical; only the control lambda and the
    (empty-peak) xls differ. Golden vs ours:

        golden  chrM  -300  -15  230.21603      ours  chrM  -300  -14  230.21603
        golden  chrM   -15  -14  234.50050      ours  chrM   -14  -13  234.50050
        golden  chrM   -14  -13  238.00050      ours  chrM   -13  -12  238.00050
        golden  chrM   -13  -11  238.50052      ours  chrM   -12  -10  238.50052

    So every boundary from -15 onward is displaced by exactly +1 in ours, while the track is
    identical before -300. F255 below shows this "shift" framing is misleading -- the values
    themselves differ by a 1-3% under-count at 129/5117 positions, and the apparent shift is an
    artifact of the bedGraph's run boundaries. The second affected file is `sweep/gmini_mfrag_d400_w600_ctrl_384`
    variant `B`, same signature (`-300 -9 179.14348` vs `-300 -8 179.14348`).

    Ruled out by inspection, not assumed:
    - Not a fragment-end inclusive/exclusive convention: no `ctrl.frag` fragment ends anywhere
      in `[-17, -13]`, and the lag starts 285 bases into the run.
    - Not a global window-offset convention: that would displace the whole track, not leave
      `[-303, -15)` exact.
    - Not the treatment path: `treat_pileup.bdg` is byte-identical, so the control-only lambda
      ladder or its window is implicated, not the pileup.

    Two readings fit the output equally so far -- an event at -15 is *lost*, or the track is
    *shifted* by 1 for x >= -15 -- because the value changes at nearly every boundary here, so
    a shift and a dropped run are indistinguishable in the bedGraph alone. Next experiment: dump
    our control ladder before coalescing (`MACS3_RS_DUMP_CTRL_LAMBDA`) for this fixture and
    check whether a run boundary exists at -15 and is lost during coalescing, or never formed.
    That distinguishes "coalescing drops a run" from "lambda is built with a 1-base lag".
    CAUTION: do not compare that dump against upstream's bedGraph position-by-position as-is.
    F254 shows the dump's `pos` column labels each run by its END, so `dump(p)` holds the value
    upstream writes for the interval `[p-1, p)`. Aligning the columns without that offset
    manufactures 448 spurious mismatches out of 535 positions -- the same class of phantom as
    F243 and F246.

    Note for the `fold_enrichment` class: the 86-file residual is *bidirectional* (ours
    `5.97001` < golden `5.97002`, but ours `8.83046` > golden `8.83045`), which contradicts the
    "we are exactly 1 ulp high" framing of F247-F249. A control track that disagrees at some
    positions -- as here -- explains that sign flip far better than any final-multiply rounding
    story does, so the residual should be re-scoped against control-track parity before more
    arithmetic is attempted.

### F254: `MACS3_RS_DUMP_CTRL_LAMBDA` labels each run by its END position, not its start

Ours-only, measured on `sweep/gmini_mfrag_d400_w180_ctrl_357` variant `B` by comparing the
diagnostic dump against our own written bedGraph:

    dump(p) == our_bedGraph(p-1)   535 / 535 positions
    dump(p) == our_bedGraph(p)      95 / 535 positions

So the dump holds the correct final control track values -- it is NOT a different intermediate, as
I first assumed when it disagreed with upstream far from the real bug -- but its `pos` column
carries `end` where the bedGraph interval `[start, end)` carries `start`. In other words `dump(p)`
is the value upstream would write for the interval `[p-1, p)`.

Consequence: F253's proposed experiment must align the columns with that offset. Comparing them
naively manufactures 448 spurious mismatches out of 535 positions, which is exactly the kind of
phantom that F243 and F246 were. Left in place deliberately rather than "fixed", because the
writer's indexing is not obviously wrong on its own terms -- it writes `pos[i]` for run `i`, and
`pos` is an array of boundaries -- so the offset is a property of how the columns line up, not
necessarily a defect. Anyone using this hook must offset by one.

This also retires a wrong intermediate conclusion: I had recorded the dump as a pre-coalescing
stage distinct from the written track. It is not.

### F255: F253's control-lambda defect is an under-count, not a coordinate shift -- and the
    obvious parameter tweak is a no-op

Sharpening F253 with measurement rather than inference.

**The lag was a coincidence of the region.** Matching our value at `p` against upstream's at
`p-k` for `k` in 0..3 over all 5117 shared `chrM` positions gives an offset that alternates 0/1
sporadically from -15 to 199, with positions matching nothing at all. Comparing the *values*
directly is decisive:

| pos | upstream | ours | diff |
|---|---|---|---|
| 20 | 367.00079 | 357.50076 | -9.50003 |
| 24 | 381.00082 | 371.50079 | -9.50003 |
| 25 | 392.00085 | 381.00082 | -11.00003 |
| 14 | 350.00076 | 350.00076 | 0 |

So we are systematically **lower**, by a variable 1-3%, with exact agreement at scattered
positions. F253's "one-base lag" and "event lost at -15" readings are both wrong as descriptions of
the track; the values themselves are wrong. Only `129/5117` positions differ, which is why the
bedGraph still looks like a boundary shift.

**Isolating the window.** Comparing against freshly generated oracle baselines for the *same* flags
(an earlier attempt compared against the default-flag golden and was therefore meaningless):

| flags | mismatches | span | ours - upstream at first |
|---|---|---|---|
| default | 129/5117 | -15 .. 128 | -4.28447 |
| `--slocal 200 --llocal 200` | 82/298 | 0 .. 128 | -1.50497 |

Two conclusions. The mismatch is confined to a ~129-position window and does **not** reach the
chromosome's negative extension once the windows are small enough to stay inside `chrM`; and it
survives with every window reduced to 200 bp, so it is not caused by wide-window clipping at
negative coordinates.

**The obvious fix does nothing.** Widening `SingleEndParams::symmetric` by one base
(`three_shift = h + 1`) changed the output by **zero bytes** -- 129 mismatches and
`sum|diff| 1366.788` both before and after. The reason: `ctrl.frag` carries five columns, so
`cc.has_counts()` is true and this fixture takes `pe_ctrl_params_counted` -> `SingleEndParams::centred`,
not `pe_ctrl_params` -> `symmetric`. Two lessons: FRAG fixtures with a count column must be tested
against the counted path, and a parameter tweak that changes nothing is evidence the hypothesis was
never exercised, not that the bug is insensitive to it. The change was reverted and
`crates/macs-pileup/src/lib.rs` is byte-identical to its pre-experiment state.

CLOSED BY F257: the under-count's cause was the `d`-scale window **width**, which was one integer
short because the ladder used the local post-filter mean instead of `self.d`. Next step at the time
of writing -- isolate which scale dominates -- turned out to be unnecessary once the coverage
arithmetic was done directly: coverage 734 (`half_d` 72) vs 715 (`half_d` 71) identified it.

### F256 RETRACTED (superseded by F257): the counted control lambda is NOT the LRC two-window
    projection -- the F151 comment misleads here

F251/F255 pointed at `PileupV2.pileup_from_LRC_centers_as_list` as the obvious suspect, because its
docstring and our own F151 note both say the counted paired-end control builds one window per
fragment *end*:

    start_poss = np.concatenate((LRC_array['l'] - half_d, LRC_array['r'] - half_d))   # PileupV2.py:918
    end_poss   = start_poss + d                                                        # :919
    weights    = np.concatenate((LRC_array['c'], LRC_array['c']))                     # :920

We were passing the fragment's own span `[l, r)` once, with weight `c`. So I implemented upstream's
construction verbatim -- two `d`-wide windows, one anchored on each end, each carrying the full
count, with `end_poss` derived from the *unclipped* start and both clipped independently
afterwards. It compiled and it was wrong:

| control-lambda variant | chrM mismatch | ours - upstream at first mismatch |
|---|---|---|
| single window `[l, r)`, weight `c` (current) | 178/269 | **-9.00002** (under-count) |
| two windows, weight `c` each (upstream's literal code) | 270/270 | **+275.50060** (over-count) |

Two `d`-wide windows per fragment against one fragment-length window, at the same `scale_factor`,
over-counts by ~2.8x -- which is the ratio of `2*d` to the mean fragment length. Upstream's lambda is
within 1-3% of our single-window value, so upstream cannot be multiplying the coverage by two at
this scale. The literal reading of `pileup_from_LRC_centers_as_list` is therefore **not** what
produces the control lambda at these scales, whatever F151's note says.

Reverted: `crates/macs-peaks/src/callpeak.rs` is back to the single-window form and reproduces the
original 178/269 exactly. This is a *negative* result, and a useful one -- it removes the most
plausible-sounding hypothesis and, more importantly, it means the `l`/`r`/`c` semantics of the LRC
array and the `scale_factor` normalisation are entangled in a way that has to be read from
`PETrackII`/`PeakDetect` before another window experiment is worth running. F255's under-count
remains open and unexplained; do not re-attempt the two-window form without first resolving where
upstream halves it.

### F256 RETRACTED: the two-window experiment was invalid, and the real cause was the window width

F256 reported that upstream's counted control is *not* the LRC two-window projection, on the
evidence that implementing it over-counted by ~2.8x. **That test was wrong and the conclusion is
withdrawn.** `macs_pileup::pileup_from_weighted_positions` takes 5' *positions* and applies the
params' shifts to build each window -- it does not accept explicit `[start, end)` bounds. So the
replacement code passed pre-shifted boundaries as positions and the window was applied *twice*, once
from my explicit bounds and once from `centred(d)`'s shifts. The 2.8x was my bug, not upstream's
behaviour.

We were already implementing the LRC construction correctly: `cstarts` receives each fragment's `l`
and `cends_v` its `r`, and `SingleEndParams::centred` applies the *same* window to minus positions as
to plus ones (`macs-pileup/src/lib.rs:401`, "no plus/minus mirror"), which is exactly upstream's
"both ends, same window" (`PileupV2.py:918-920`).

The actual fault was the window **width**. See F257.

### F257: the paired-end lambda ladder's width must be `self.d`, not the local post-filter mean

**The fix.** `crates/macs-peaks/src/callpeak.rs`, `build_signals_pe`:

```rust
let scales: Vec<i64> = std::iter::once(cfg.tsize_exact as i64)   // was: d as i64
```

**Why.** `PeakDetect.py:192` builds `ctrl_d_s = [self.d]`, where `self.d` is `options.tsize`, the
truncated *as-read* mean fragment length. But `PeakDetect.py:155` introduces a **different** `d` in
the same function:

```python
if self.PE_MODE:
    d = self.treat.average_template_length   # post-filter mean, a local
```

`control_sum` and `ratio` use that local `d`; `ctrl_d_s` uses `self.d`. F95 had already fixed this
for the *factors* (they are built from `cfg.tsize_exact`) but left the ladder *width* on the local
mean, so the two disagree whenever duplicate filtering removes anything.

**Why one integer of `d` moves so much.** The counted paired-end control is two `d`-wide windows per
fragment, one anchored on each end, each spanning `[anchor - d//2, anchor + d - d//2)`. So one
integer of `d` moves `d//2` by half a base and narrows **both** windows.

**The arithmetic that pins it.** On `sweep/gmini_mfrag_d400_w180_ctrl_357` variant `B`:

| quantity | value |
|---|---|
| local `d` (post-filter mean), truncated | **143** (`half_d` 71) |
| `cfg.tsize_exact` = `options.tsize`, truncated | **144** (`half_d` 72) |
| `ratio_treat2control` | 0.5000010873654657 |
| upstream lambda at chrM:20 | 367.00079 |
| ours | 357.50076 |
| coverage implied by upstream | **734** |
| coverage implied by ours | **715** |

Evaluating the real `ctrl.frag` under each model settles it: coverage at chrM:20 is **734** with
`half_d = 72` and **715** with `half_d = 71`, and nothing else in the neighbourhood produces either
number (738 for `half_d` 73, 196 for the fragment span alone, 38 for a single centred window). That
is why this is the fault and not a coincidence.

**Verified effect.** `chrM` control lambda now matches upstream on **5117/5117** positions for the
default flags and **270/270** with `--slocal 0 --llocal 0`. Corpus-wide (`oracle/run_golden.sh
--jobs 14`, 7685 cases / 22456 files):

| metric | before | after |
|---|---|---|
| output files matched | 22257/22456 | **22291/22456** (+34) |
| cases byte-identical | 7091/7685 | 7091/7685 |
| cases with differences | 147 (est.) | **113** |
| `start`-column differences | 34 | **2** |

The `start`-column class collapsing from 34 files to 2 confirms the diagnosis: F253's "mysterious
one-base shift" was this under-count all along, and the peak-boundary class was a downstream symptom
of the control lambda being slightly low, not a separate coordinate bug.

**Residual.** 113 cases still differ: 43 in `fold_enrichment`/`signalValue`, 2 in `start`, 2 in peak
`name`, 1 in `score`, 1 in `abs_summit`. The `fold_enrichment` class is the one F250-F252 established
is unobservable from inside the oracle, and F255 noted its sign is *bidirectional* -- which this fix
plausibly explains, since a lambda that is low at some positions and exact at others moves
`fold_enrichment` in both directions.

**Regression test.** `crates/macs-peaks/tests/f257_pe_ladder_width.rs` (3 tests) recomputes the
coverage from the real fixture and pins 734 vs 715, so a regression to the local mean fails with a
number rather than a bedGraph diff.

### F258: the reported paired-end `d` is the **unweighted** row mean, not the count-weighted one

Second bug from the same root confusion as F257, and it was corrupting the xls *header* of every
counted paired-end fixture.

`Parser.py:1496` computes upstream's `self.d` unweighted:

```python
m += right_pos - left_pos      # once per ROW, no multiplicity
...
self.d = cython.cast(cython.float, m) / i
```

and `callpeak_cmd.py:382` prints `options.tsize` from it. We were passing the **count-weighted** mean
(`t.length() / t.total()`), which is a different number as soon as any row carries a count -- which is
the situation in every `--format FRAG` fixture. Measured on
`sweep/gmini_mfrag_d400_w180_ctrl_357` variant `broad`:

```text
golden  # fragment size is determined as 144 bps     # d = 144
ours    # fragment size is determined as 143 bps     # d = 143
count-weighted mean -> 143.xx
unweighted row mean  -> 144.155 -> 144
```

`mean_row` is the quantity F166 already routes to `cfg.tsize_exact`, so after the fix the header and
the scale factors agree on one `d` instead of two.

**Why the harness hid this.** `oracle/run_golden.py::_first_diff` skips lines where *both* sides start
with `#`, and falls through to `"bytes differ"` when the only differences are comment lines. So these
cases were reported as an uninformative `bytes differ` while the real fault sat in the
`# fragment size ...` and `# d = ...` lines. Reading the reported column census alone would have
sent this class nowhere; it only surfaced by reproducing one case by hand and diffing it.

**Verified effect.** Corpus-wide (`oracle/run_golden.sh --jobs 14`, 7685 cases / 22456 files):

| metric | before F258 | after F258 |
|---|---|---|
| cases byte-identical | 7091/7685 | **7156/7685** (+65) |
| output files matched | 22291/22456 | **22356/22456** (+65) |
| cases with differences | 113 | **48** |

Combined with F257 the session total is +65 cases and +99 files.

**Residual after both fixes: 48 cases**, of which 43 are last-digit-only `fold_enrichment` /
`signalValue`, plus 2 `start`, 2 `name`, 1 `score`, 1 `abs_summit`. The 43 are now concentrated in
`gonechrom_mpe_d4000_w180_noc_131` (13), `gtiny_mfrag_d400_w600_noc_142` (12) and
`se_model/realistic` (6), i.e. the no-control paired-end path plus single-end.

**The residual is signed per path**, which is worth recording because it rules out a single
systematic rounding error:

| fixture | golden | ours | sign |
|---|---|---|---|
| `gtiny_mfrag_d400_w600_noc_142` (PE, no control) | 8.92562 | 8.92563 | ours high |
| `gonechrom_mpe_d4000_w180_noc_131` (PE, no control) | 32.3792 | 32.3793 | ours high |
| `se_model/realistic` (SE, with control) | 5.97002 | 5.97001 | ours low |

The two no-control paired-end cases share a sign and a path; single-end has the opposite sign. That
is consistent with F252's conclusion that `summit_ctrl` cannot be read from inside the oracle, but it
also means the no-control PE path deserves its own look: its factor already matches upstream's
expression (`treat.length / (lregion * treat.total * 2)`, `PeakDetect.py:342`), so whatever remains
is in the window projection rather than the factor.

**Not claimed.** F252's "unobservable" verdict was reached while the control lambda was wrong (F257),
so it should be re-examined now that the lambda is right -- but nothing in this turn reopens it, and
the bidirectional sign split means it is not a one-ulp story.

### F259: `fold_enrichment` is stored narrowed to **f32** -- this was the whole `fold_enrichment` class

The last-digit residual that F246-F249, F252 and F255 all chased, and which F252 had concluded was
unobservable from inside the oracle. It was neither unobservable nor a wrong operand: both operands
were right and the **width of the stored result** was wrong.

`CallPeakUnit.py:1385`:

```python
fold_change = (summit_treat + self.pseudocount) / (summit_ctrl + self.pseudocount)
```

`summit_treat` and `summit_ctrl` hold C floats and `pseudocount` is a `cython.double`, so C's usual
arithmetic conversions promote the whole expression to **double** -- an f64 divide. But the value
lands in `PeakContent.fc`, declared `cython.float` (`PeakIO.py`, alongside `pileup`, `pscore` and
`qscore`), so it is **narrowed to f32** before any writer sees it. We kept f64 all the way to
`%.6g`.

Measured on `sweep/gonechrom_mpe_d4000_w180_noc_131`, whose control lambda, treatment pileup and
summit position are *already byte-identical to upstream*, leaving this as the only remaining step:

| step | value | `%.6g` |
|---|---|---|
| f64 divide, kept as f64 | 32.379250536484484 (f64 `0x4040308b4815987d`) | 32.3793 -- **wrong** |
| f32 arithmetic | 32.379249572753906 (f32 `0x4201845a`) | 32.3792 -- upstream |
| f64 divide, then f32 | 32.379249572753906 (f32 `0x4201845a`) | 32.3792 -- upstream |

The gap is ~3e-8 relative -- a single f32 ulp -- but `%.6g` puts the two candidates on opposite
sides of a rounding boundary. That is the entire mechanism behind a residual that looked like a
data problem for five findings.

Note that "f32 divide of f32 sums" and "f64 divide narrowed to f32" agree here. They are not
equivalent in general, so the implementation does the f64 divide and narrows the stored value, which
is what the C semantics actually describe.

**Why five findings missed it.** Every previous attempt tried to explain the residual as a *wrong
`summit_ctrl`*: F246 inverted a rounded print as if exact (4.2e-6, wrong by orders of magnitude),
F247 concluded "exactly 1 ulp high" from a single case, F248 built a theory about the `ratio` and
then disproved it, F249 proved `summit_ctrl` unobservable and F252 proved the peak tuple
unobservable too. All of that was sound work on an unanswerable question. The decisive move was to
stop inverting the printed value and instead ask what the printed value is a function *of*: with the
lambda and pileup already byte-identical, `fold = (pileup + 1) / (ctrl + 1)` has only two unknowns,
the **width** and the operands -- and the operands were pinned by the byte-identical tracks. F247's
own warning ("do not read a precision out of a `%.6g` field") is what made the bucket-width
realisation available; it just needed to be applied to the *output* rather than the *input*.

**Verified effect.** Corpus-wide (`oracle/run_golden.sh --jobs 14`, 7685 cases / 22456 files):

| metric | before F259 | after F259 |
|---|---|---|
| cases byte-identical | 7156/7685 | **7200/7685** (+44) |
| output files matched | 22356/22456 | **22445/22456** (+89) |
| cases with differences | 48 | **4** |

All 43 `fold_enrichment` / `signalValue` cases are resolved. Session total across F257, F258 and
F259: **+109 cases, +188 files**, and differing cases 147 -> 4.

**Regression tests.** Two in `crates/macs-peaks/src/regions.rs` pin the f32 narrowing on the real
operands (including the f64 bits that must *not* be stored) and check both the summit and broad
paths round-trip as f32.

**The 4 remaining cases are four distinct edge-case bugs**, none of them numeric drift:

| fixture | symptom |
|---|---|
| `gmini_mfrag_d4000_w600_ctrl_078/call_summits` | we emit two sub-summits (`peak_1a`, `peak_1b`), upstream one (`peak_1`) |
| `gmini_mfrag_d400_w600_ctrl_222/nolambda` | peak starts at -5000 instead of 0; length 5200 vs 200 |
| `gmini_mfrag_d400_w600_ctrl_384/call_summits` | summit one base late (45 vs 46), pileup 734 vs 731 |
| `gmini_mpe_d400_w180_ctrl_192/broad` | upstream calls one broad peak, we call none |

The sub-peak case has a concrete lead: `CallPeakUnit.py:1440` computes `start = max(peak_start - 10,
0)`, and our `peak_start.saturating_sub(10)` does **not** clamp at zero (the comment above it claims
it does). On a peak starting at -12 upstream pads from 0 and then hits numpy's negative-index slice
semantics in `peakdata[m:n]`, which our `as usize` cast skips instead of wrapping. Not yet fixed --
`Coord` is `u64`, so the negative coordinates these fixtures exercise come from a shift applied
somewhere upstream of this code, and that needs to be understood before changing the indexing.

### F260: the sub-peak padding clamp and `start_boundary` are signed -- fixed a summit off-by-one

`CallPeakUnit.py:1440` and `:1444`:

```python
start = max(peak_start - 10, 0)
start_boundary = peak_start - start
```

Both depend on upstream's **signed** coordinates. A peak starting before the contig start makes the
clamp bite, and `start_boundary` then goes **negative** and is used as a signed bound by the
`np.searchsorted` pair that trims maxima out of the padding.

We hold coordinates as `u64`. For counted paired-end tracks, upstream's negative excursion is
emulated by adding `coord_shift = max(d, slocal, lregion) / 2` on the way in and subtracting it again
when writing, so upstream's `0` is `coord_shift` here. Two consequences, both fixed:

1. `peak_start.saturating_sub(10)` did **not** clamp at zero -- `saturating_sub` only protects
   against underflow, so a peak 12 bases left of the contig got the full 10 bp of padding and
   reached further left than upstream, putting extra maxima inside `maxima()`. The clamp floor is
   now `CallParams::clamp_floor`, set from `ChromCall::clamp_floor` (the shifted `coord_shift` on the
   paired-end path, `0` everywhere else). Note the comment above that line already *claimed* the
   clamp existed; the code did not do it.

2. `start_boundary = peak_start - start` was computed in `u64`, so for a peak at -12 it wrapped to
   ~2^64 and the padding filter kept a maximum from the **padding** instead of the peak body. That
   is how this class first presented: a summit reported at 180 where upstream reports 64. It is now
   `i64` throughout, with the `partition_point` bounds compared as signed values.

**Verified effect.** `sweep/gmini_mfrag_d400_w600_ctrl_384/call_summits` is **fixed outright** -- its
summit moved from 46 to upstream's 45, and the pileup/q-score followed. Corpus-wide:

| metric | before F260 | after F260 |
|---|---|---|
| cases byte-identical | 7200/7685 | **7201/7685** |
| output files matched | 22445/22456 | **22448/22456** |
| cases with differences | 4 | **3** |

**Not fixed: `gmini_mfrag_d4000_w600_ctrl_078/call_summits`.** Its *first* sub-summit is now
byte-identical to upstream in every field (start, end, length, abs_summit, pileup, p/q scores,
fold_enrichment) and only the name differs, because we still emit a second, spurious sub-summit at
abs_summit 178 with `peak_1a`/`peak_1b` where upstream emits a single `peak_1`.

The remaining cause is numpy's **negative-index slice** semantics, which the fix above exposed rather
than resolved. With `start = 0` and a chunk starting at -12, upstream evaluates
`peakdata[-12:180] = tscore` on an array of length 190, which numpy resolves to `peakdata[178:180]`
-- the chunk's value lands at the *end* of the array. Our fill loop casts to `usize`, so that chunk
is skipped instead (`if m >= width { continue }`), and `peakindices` there stays `-1` for the wrapped
chunk. Replicating this faithfully means implementing negative-slice resolution for `m` and `n`
independently and in numpy's order, and it interacts with which chunk writes last.

**Closed by F261 below**, by transcribing numpy's slice rule rather than fitting the fixture.

### F261: numpy's negative-index slice, faithfully transcribed -- closes the sub-peak case

F260 fixed the first sub-summit of `sweep/gmini_mfrag_d4000_w600_ctrl_078/call_summits` but left a
spurious second one at `abs_summit 178`. The cause is numpy slice semantics.

`CallPeakUnit.py:1446-1450` fills the padded arrays with plain numpy slicing:

```python
m = tstart - start
n = tend - start
peakdata[m:n] = tscore
peakindices[m:n] = i
```

With `start = max(peak_start - 10, 0)` clamped to 0 and a chunk starting at -12, `m` is -12. numpy
resolves a negative slice bound by **adding the length and then clipping** (`slice.indices`), so on
an array of length 190, `peakdata[-12:180]` addresses `peakdata[178:180]` -- the chunk's value lands
at the far end of the array. We were casting to `usize`, which sent such chunks past the
`m >= width` guard and dropped them, leaving `peakindices` at -1 there.

The fix transcribes the rule exactly:

```
start = bound + len if bound < 0 else bound,   then clamp to [0, len]
```

applied independently to `m` and `n`. That is numpy's documented behaviour rather than a fit to the
fixture, which is what F260 said was needed before touching the loop.

**Verified.** `sweep/gmini_mfrag_d4000_w600_ctrl_078/call_summits` now matches on all three files
(`_peaks.xls`, `_peaks.narrowPeak`, `_summits.bed`). Corpus-wide:

| metric | before F261 | after F261 |
|---|---|---|
| cases byte-identical | 7201/7685 | **7202/7685** |
| output files matched | 22448/22456 | **22451/22456** |
| cases with differences | 3 | **2** |

### F262 (attempted, reverted): `--nolambda`'s substitute control run is not the source

`sweep/gmini_mfrag_d400_w600_ctrl_222/nolambda` reports the peak at -5000 where upstream reports 0,
uniformly offset by exactly `coord_shift` = 5000 (start -4999 vs 1, length 5200 vs 200, summit offset
5178 vs 178). The hypothesis was that the `--nolambda` substitute control run
(`SignalTrack::from_runs_exact(chrom, 0, t.end(), ..)`) begins at `0` while the treatment track begins
at `coord_shift`, so `__chrom_pair_treat_ctrl` -- which pairs from `min(t_p[0], c_p[0])` -- pulls the
whole paired array a full shift early, and the `- coord_shift` applied on write then lands 5000 low.

Changing the run's start to `coord_shift` compiled and changed **nothing**: the peak still reported
-4999. So the region start is not set by the control track on this path -- see **F263** below, which
pins it on the q-track instead. Reverted -- an inert change
that cannot be justified by evidence is a liability, and this session's standing rule has been to
leave the tree in a verified state.

The uniform 5000 offset means the shifted frame is correct in magnitude but something on this path
never enters it. Next step: dump the treatment track's own start (`t.start()`) and the paired array's
first position for this fixture and compare against `coord_shift`, rather than reasoning about which
track *should* set it.

**Two cases remain**, both now characterised:

| fixture | symptom |
|---|---|
| `sweep/gmini_mfrag_d400_w600_ctrl_222/nolambda` | uniformly `coord_shift` low; summits.bed already byte-identical |
| `sweep/gmini_mpe_d400_w180_ctrl_192/broad` | upstream calls one broad peak, we call none |

### F263 (REVERTED): a "obviously correct" pileup change that cost 474 cases

Kept because the failure is instructive, not because the change was good.

Measurement on `sweep/gmini_mfrag_d400_w600_ctrl_222 --nolambda` (F262) showed
`coord_shift = 5000`, fragments correctly shifted to `(5000, 5180)`, but the resulting treatment
track reporting `start() == 0`. The sweep in `pileup_from_weighted_fragments_from` plainly begins at
`origin`, so the track *claiming* to start at 0 while its first breakpoint sits at 5000 looked like a
straightforward defect:

```rust
let mut b = TrackBuilder::with_capacity(chrom, origin, rlength, deltas.len());  // was: 0
```

A temporary `MACS3_RS_PROBE_FRAME` probe confirmed it took effect: `t.start()` became 5000. And it
changed **nothing** for the fixture it was aimed at -- the peak still reported -4999.

**Then the corpus measurement:**

| metric | with F263 | reverted |
|---|---|---|
| cases byte-identical | 6728/7685 | **7202/7685** |
| output files matched | 21257/22456 | **22451/22456** |
| cases with differences | 476 | **2** |

474 cases broken by a change that looked self-evidently correct and passed a targeted fixture test.
Reverted; the tree is back to 7202/22451 with 2 differing cases.

What this says about the earlier work: `TrackBuilder::with_capacity`'s `start` argument is evidently
*not* "the coordinate the track begins at" in the sense the name suggests -- most likely it is an
internal cursor floor that the sweep relies on, and anchoring it at `origin` suppresses the leading
zero-depth run that upstream's arrays do carry. That is a real invariant of the builder, not
something inferable from one call site.

**Method note.** This is the second time this session a hypothesis was validated only against the
fixture that motivated it. F262's control-run change was caught the same way (inert), but this one
slipped through because it *did* change the probed value -- which felt like progress. A targeted
fixture test cannot distinguish "fixed the case" from "moved a number the probe reports"; only the
corpus can distinguish those. Every candidate fix from here runs the full golden before being kept.

### Where the `nolambda` case actually points

The probes did narrow it, though. With `t.start() == 5000` **and** the substitute control run also
starting at `coord_shift`, the peak still reported internal 0. The measurement that explains it:

```text
PROBE3 chrom=chrM treat.start=5000 ctrl.start=Some(5000) qtrack.start=0 qtrack.len=0
```

`qtrack.start == 0` with `qtrack.len == 0` -- that is the `EMPTY_TRACK` fallback reached through
`qtracks.get(k).unwrap_or(&EMPTY_TRACK)`, not a built q-track. So under `--nolambda` no q-track is
produced for this chromosome, and whatever anchors the chunk positions is anchored at 0, which is why
the `- coord_shift` on write lands 5000 low. The next step is therefore to find out **why the q-track
list is short enough for the fallback to trigger**, which is a `build_qtable` question, not a
coordinate-frame question at all.

### F265 (negative, reverted): the broad-peak case is a one-base region boundary, not a cutoff

`sweep/gmini_mpe_d400_w180_ctrl_192 --broad` is the one remaining case where we emit output and
upstream emits none. (Earlier notes had this backwards; the golden `broadPeak` and the golden `xls`
are both empty of data rows, and *we* produce
`chrM 55 200 146 103 5.37852 1.57585 4.76168`.)

Upstream calls no peak in **any** of this fixture's 20 variants, and our narrow path agrees -- so the
divergence is specific to the broad path.

**Hypothesis 1, tested and rejected: a missing early return.** `__chrom_call_broadpeak_using_certain_criteria`
does compute lvl1 `above_cutoff` first and `return` when it is empty:

```python
above_cutoff = np.nonzero(apply_multiple_cutoffs(score_array_s, lvl1_cutoff_s))[0]
if above_cutoff.size == 0:
    return
```

Adding the equivalent guard (`if chunks.is_empty() { return }`) changed nothing, and correctly so: the
guard is about *positions clearing the cutoff*, and this fixture has some -- our narrow path drops
them later on `min_length`, not before. The guard was removed again rather than left inert.

**Hypothesis 2, and where it actually lands.** The broad caller's only real filter is in
`__close_peak_for_broad_region` (`CallPeakUnit.py:2193`):

```python
peak_length = peak_content[-1][1] - peak_content[0][0]
if peak_length >= min_length:      # too small -> reject
```

For this fixture `d = 145` (and `--extsize 200` is ignored in PE mode), so `min_length = 145`. Our
region is `55 .. 200`, i.e. `peak_length = 200 - 55 = 145` -- **exactly** `min_length`, so `>=`
accepts it. For upstream to reject it, upstream's region must be at least one base shorter
(`144 < 145`), which puts the difference squarely on the **broad region's start boundary**, not on the
cutoff, the score, or the length rule.

That is consistent with `--broad-cutoff 0.1` being *less* stringent than the 0.05 q-value: at lvl2 a
146-base region exists for us and evidently not for upstream, so either upstream's lvl2 chunk starts
one base later, or upstream's `lvl1_max_gap` merge differs (`broad_level` uses `params.max_gap * 4`
whereas the caller passes `lvl1_max_gap // 2`).

Not fixed. The next step is to dump the lvl2 above-cutoff positions for this fixture and compare the
first chunk's start against `pos_array[above_cutoff-1]`, taking upstream's
`if above_cutoff[0] == 0: above_cutoff_startpos[0] = 0` special case into account. No code change
should be made until that first position is known, for the reason F263 recorded: a change that looks
obviously correct and passes only the fixture that motivated it cost 474 cases.

### F266: the broad-path level sets are the missing intermediate, and the divergence is in lvl1

`MACS3_RS_DUMP_BROAD_LEVELS` now dumps the per-chromosome lvl1 and lvl2 above-cutoff chunk sets on
the `--broad` path (`crates/macs-peaks/src/callpeak.rs`). Added because the porting plan lists
"candidate/merged intervals" as a captured stage, the broad level sets were not among them, and
`sweep/gmini_mpe_d400_w180_ctrl_192 --broad` is undiagnosable without them: when the two sides
disagree the golden `*_peaks.xls` is simply empty, so nothing downstream carries the information.

**Measured on that fixture** (`min_length = 145`, `max_gap = 145`, `lvl1_cut = 1.301030`,
`lvl2_cut = 1.000000`):

| level | chunks | gap | merged region | `peak_length` | verdict |
|---|---|---|---|---|---|
| lvl1 | `57..128`, `147..200` (n=116) | 19 | `57..200` | 143 | **reject** (`< 145`) |
| lvl2 | `54..128`, `142..200` (n=124) | 14 | `54..200` | 146 | **accept** (`>= 145`) |

Both gaps are far below `max_gap`, so each level merges to a single region. That is exactly our bug:
**lvl1 is empty, lvl2 has one region, so we emit one broad peak and upstream emits none.**

**Why this cannot be resolved from the outputs.** `CallPeakUnit.call_broadpeaks` consumes the two
level sets like this:

```python
try:
    lvl1 = lvl1peakschrom_next()
    for i in range(len(lvl2peakschrom)):
        lvl2 = lvl2peakschrom[i]
        ...
except StopIteration:
    self.__add_broadpeak(broadpeaks, chrom, lvl2, tmppeakset)   # `lvl2` assigned only inside the loop
```

`lvl2` is bound only inside the `for`, so an **empty lvl1 set raises `StopIteration` before `lvl2`
is ever assigned** and the handler dereferences an unbound local. Upstream exits 0 on this fixture,
which means its lvl1 set is *not* empty -- so upstream's lvl1 region is at least 145 long, i.e. it
starts at 55 or earlier, against our 57. Yet upstream also emits no broad peak, which its
`__add_broadpeak` path (`CallPeakUnit.py:2247`, which emits even for an empty lvl1 set) says should
not happen.

So at least one of upstream's two level sets differs from ours, and the argument above shows we
cannot tell which from the outside: every combination we can construct predicts either output or a
crash, and upstream produced neither.

**Blocked on oracle-side capture.** Upstream's lvl1/lvl2 peak sets live in `PeakIO` objects built by
`__close_peak_for_broad_region` and assembled in `call_broadpeaks`. `PeakIO` is an immutable cdef
class (F250) and `call_broadpeaks` is a method on a cdef class reached through a static vtable slot
(F251c), so neither can be wrapped -- the same wall that made `summit_ctrl` unobservable. Closing
this case needs either a patched oracle build or a level-set dump added to upstream's source, which
breaks the "pristine oracle tree" rule this project keeps.

Practical next step that does *not* need the oracle: dump our q-score for `chrM:50..60` and compare
against `lvl1_cut = 1.301030`. Our lvl1 starts at 57 and lvl2 at 54, so the q-track in `54..56` sits
close to the boundary; if it is within an ulp of 1.301 there is a single-ulp q-track defect here that
would explain a 2-base lvl1 region difference and is fixable with evidence.

### F267: the `--broad` case is a q-table difference, and the "single ulp" hypothesis is dead

Continuing F266. The F266 next step was to check how close our q-track sits to `lvl1_cut` in
`54..56`, on the theory that a near-boundary single-ulp q-track defect would explain a 2-base lvl1
region difference. **It does not, and the theory is refuted by the numbers.**

`MACS3_RS_DUMP_BROAD_LEVELS` now also prints the per-position `(pos, treat, ctrl, p, q)` around the
lvl1/lvl2 divergence. For `sweep/gmini_mpe_d400_w180_ctrl_192 --broad`
(`lvl1_cut = 1.301030`, `lvl2_cut = 1.000000`):

| pos | treat | ctrl | p | q | |
|---|---|---|---|---|---|
| 54 | 50 | 40.769798279 | 1.169990 | 0.992812 | |
| 55 | 51 | 41.154418945 | 1.240570 | 1.053483 | above lvl2 only |
| 56 | 52 | 41.154418945 | 1.369760 | 1.175940 | |
| 57 | 53 | 41.539039612 | 1.445320 | 1.244661 | |
| 58 | 54 | 41.923660278 | 1.522420 | **1.314812** | first above lvl1 |

The lvl1 boundary is nowhere near a rounding edge: pos 57 is `0.0564` **below** the cutoff and pos 58
is `0.0138` **above** it. Moving the boundary would need a ~4% q-score error, so this is not an
arithmetic-width issue of any kind.

**What the numbers actually force.** Our q-track reaches **2.340059** on this chromosome and has
**13** positions above `lvl1_cut` in the sampled window. Work back through F266's control flow:

- If upstream had any position above `lvl1_cut`, `above_cutoff.size != 0`, the function proceeds,
  and since lvl2's chunk set is a superset of lvl1's, a lvl1 region that clears `min_length`
  forces the lvl2 region to clear it too -- so `__add_broadpeak` would emit something.
- If upstream had none, the `if above_cutoff.size == 0: return` early-out fires, `call_broadpeaks`
  is never reached, and upstream emits nothing with rc 0.

Upstream emits nothing with rc 0. Therefore **upstream has zero positions above `lvl1_cut` on
`chrM`**, against our 13 (and our maximum of 2.34). That is not a rounding difference; it is a
different q-value table.

**Why it is almost certainly the AFDR histogram, not the score.** The p-scores are plausible and
smooth (1.17 -> 1.52 across the boundary, monotone in treatment depth), and the control lambda is
exact everywhere else in the corpus. `q = self.pqtable[pscore]` is a *table lookup*, so the only way
to move every q-score downward by this much is for the table itself to differ. This fixture is a
degenerate case for that construction: `chrM` is 200 bp inside a `gsize` of 48200, so the p-score
histogram has very few occupied bins, and `__cal_pvalue_qvalue_table` normalises by total span.

**Next step, with evidence:** dump `self.pqtable` (the p->q map) for this fixture and compare its
entries against our `PqTable`. Upstream's is reachable without touching the oracle -- `--cutoff-analysis`
writes the p/q cutoff ladder that the table is seeded from (`callpeak_cmd` +
`__cal_pvalue_qvalue_table`, line 978), and `oracle/check_cutoff_analysis.sh` already exists to
compare it. If the ladders agree, the divergence is in the histogram accumulation, not the seeding.

### F268: the q-table hypothesis is refuted, and the remaining contradiction is in my own transcription

F267 proposed that `sweep/gmini_mpe_d400_w180_ctrl_192 --broad` diverges because upstream has zero
positions above `lvl1_cut` while we have 13 (max 2.34) -- i.e. a different p->q table. **That is
wrong, and the refutation is cheap.**

`--cutoff-analysis` writes the p/q cutoff ladder the table is seeded from
(`__cal_pvalue_qvalue_table`, line 978), and `oracle/check_cutoff_analysis.sh` already exists to
compare it. Running both sides on this fixture:

```text
pscore  qscore  npeaks  lpeaks  avelpeak
1.20    1.02    1       146     146.00
0.90    0.76    1       149     149.00
0.60    0.49    1       154     154.00
0.30    0.21    1       161     161.00
```

**Byte-identical.** And the ladder is consistent with our per-position dump: our `q` at `p = 1.16999`
is `0.992812`, so `p = 1.20 -> q ~ 1.02` is exactly what our table predicts. The p->q map agrees, the
p-scores agree, the control lambda agrees, and therefore **our q-track agrees with upstream's**, along
with the lvl1/lvl2 chunk sets derived from it. F267's inference -- "upstream has no position above
`lvl1_cut`" -- is therefore false, and with it the claim that the residual is a q-table defect.

**What survives is a contradiction, and it points at my own reading.** Our measured level sets:

| level | chunks | merged region | `peak_length` |
|---|---|---|---|
| lvl1 | `57..128`, `147..200` | `57..200` | 143 |
| lvl2 | `54..128`, `142..200` | `54..200` | 146 |

lvl2's chunk set is a strict **superset** of lvl1's (lower cutoff, `1.0` vs `1.301`), so
`lvl2_region ⊇ lvl1_region` and the lvl2 region is always at least as long. Then:

| `min_length` | lvl1 | lvl2 |
|---|---|---|
| 145 or 146 | reject | **accept** |
| >= 147 | reject | reject |

So no value of `min_length` yields "lvl1 non-empty, lvl2 empty" -- and that combination is the only
one consistent with upstream's observed behaviour. Reading `call_broadpeaks`
(`CallPeakUnit.py:1888-1922`) as I transcribed it:

- lvl1 empty, lvl2 non-empty -> `lvl1 = lvl1peakschrom_next()` raises `StopIteration` before `lvl2`
  is bound, and the `except` handler dereferences it;
- both empty -> the same unbound `lvl2`;
- lvl1 non-empty, lvl2 empty -> the `for` never runs, no exception, **no output, rc 0** (matches);
- both non-empty -> `__add_broadpeak` emits (does not match).

Only the third matches, and the arithmetic above says it is unreachable. Since the numerics are now
verified against the oracle, the remaining suspect is **my transcription of the broad assembly** -- most
plausibly `lvl2` is a Cython-typed local that is initialised rather than left unbound, or the
`except StopIteration` body differs from what I read. That is a claim about upstream's *source*, and
it is checkable directly; it is not a claim about our code.

**Next step:** re-read `CallPeakUnit.py:1875-1925` with the Cython declaration block in view (the
`lvl1`/`lvl2`/`tmppeakset` type declarations sit above the `try`), and establish what the
`except StopIteration` path actually does when both level sets are empty. Do not change our code
until that is settled -- the numbers say our level sets are right, so any change to
`chunks_at`, `broad_level` or `min_length` would be fitting to a symptom whose cause is upstream-side.

Related note on method: `--cutoff-analysis` is `store_true`, not a value. Passing
`--cutoff-analysis 0.0,0.5,...` makes argparse fail, and both sides produced *no output at all* --
which, compared naively, looks like agreement. Two "identical" results in a row from a command that
never ran.

### F269: the broad chromosome loop is driven by **lvl1**, so an empty lvl1 means silence

Closes `sweep/gmini_mpe_d400_w180_ctrl_192 --broad`, and it took reading F268's contradiction
literally rather than trying another numeric tweak.

`CallPeakUnit.call_broadpeaks` assembles the broad output like this:

```python
chrs = lvl1peaks.get_chr_names()          # <-- the LEVEL 1 peaks choose the chromosomes
for chrom in sorted(chrs):
    lvl1peakschrom = lvl1peaks.get_data_from_chrom(chrom)
    lvl2peakschrom = lvl2peaks.get_data_from_chrom(chrom)
    ...
```

The chromosome loop is driven by **lvl1peaks**, not by the input chromosome list. A chromosome with
no lvl1 peak is therefore skipped outright, and its lvl2 regions are never looked at.

That resolves everything F266-F268 could not:

- Our lvl1 region is `57..200` = **143** bases against `min_length` **145**, so lvl1 is empty.
- Upstream's chromosome never enters the loop, so the `except StopIteration` branch that would
  dereference an unbound `lvl2` is **never reached** -- upstream is silent at rc 0 rather than
  crashing. F268 had deduced that "lvl1 non-empty, lvl2 empty" was the only rc-0 outcome reachable,
  and found it arithmetically impossible; the missing piece was that the loop is keyed on lvl1, so
  "lvl1 empty" does not lead to that branch at all.
- Our `combine_broad(&l2, &l1)` had no such condition and emitted the lvl2 region (`54..200` = 146
  >= 145) regardless.

**The fix** is one guard in `crates/macs-peaks/src/callpeak.rs`, placed *after* `l1` is computed:

```rust
if l1.is_empty() {
    return Vec::new();
}
```

**Why this is the general fix and not a special case.** lvl2's above-cutoff set is a strict superset
of lvl1's (cutoff `1.0` versus `1.301`), so "lvl1 empty while lvl2 is not" is the *normal* outcome
whenever the strong level yields nothing below `min_length`. Any fixture where `--broad-cutoff` is
looser than the q-value cutoff and no strong peak clears `min_length` was affected.

**Verified.** `_peaks.broadPeak` and `_peaks.gappedPeak` are byte-identical, and the `_peaks.xls`
difference was only the `-n` name from a manual run. Corpus-wide, on the full replay:

| metric | before | after |
|---|---|---|
| cases byte-identical | 7202/7685 | **7203/7685** |
| output files matched | 22451/22456 | **22454/22456** |
| cases with differences | 2 | **1** |

**What the three findings in this chain got wrong, for the record.** F266 measured the level sets
correctly and inferred the wrong cause. F265's first guess (an early `above_cutoff.size == 0` return)
was real upstream code but fired on the chunk set rather than the peak set, so it was inert. F267
inferred a q-table defect from "upstream must have zero positions above `lvl1_cut`" -- true as logic,
false as fact, and refuted in one command by `--cutoff-analysis`. F268 then proved the numerics were
right and localised the fault to my own reading of the assembly. The actual cause was one line of
upstream source (`chrs = lvl1peaks.get_chr_names()`) that was in the same function I had been reading
for three turns, four lines above the `try` I kept quoting.

### F270 (reverted): the last case is a shifted-frame mapping, and two attempts each cost ~450 cases

`sweep/gmini_mfrag_d400_w600_ctrl_222 --nolambda` is the only remaining golden difference. Upstream
reports `chrM 0 200` with summit offset 178; we report `chrM -5000 200` with offset 5178 --
uniformly `coord_shift` (`max(d, slocal, lregion) / 2` = 5000) low in both columns, while
`summits.bed` already matches byte-for-byte.

**The fix that works on the fixture.** Upstream's broad/narrow chunk builder has an explicit
first-position special case:

```python
above_cutoff_startpos = pos_array[above_cutoff-1]
if above_cutoff[0] == 0:
    above_cutoff_startpos[0] = 0
```

That assigns the literal **coordinate 0**, which in our shifted frame is `coord_shift`, not our array
origin (`cc.treat.start()` == 0). Threading a `zero_coord` field through `ChromCall` and setting it to
`coord_shift` on the paired-end path makes `_peaks.narrowPeak` and `_summits.bed` **byte-identical**
for this fixture.

**Why it was reverted anyway.**

| metric | with F270 | reverted |
|---|---|---|
| cases byte-identical | 6759/7685 | **7203/7685** |
| output files matched | 21121/22456 | **22454/22456** |
| cases with differences | 445 | **1** |

442 cases broken by a change that is correct for the fixture it was written for and is grounded in an
explicit upstream source line. Together with F263 (474 cases) that is two independent attempts at the
same frame mapping, each costing about the same amount, and both wrong.

**What that pair of failures establishes.** Upstream's literal `0` is not a single constant in our
frame. `coord_shift` is added on the way in to make `u64` stand in for signed coordinates and
subtracted on the way out, so *which* physical coordinate a literal `0` denotes depends on the
coordinate that the surrounding arithmetic is relative to -- the track origin for the pileup builder
(F263), the array origin for the chunk builder (F270), and neither of those is the other. Our
`u64`-plus-shift emulation is faithful for values that are genuinely negative upstream (it produces
`-12` correctly on `sweep/gmini_mfrag_d4000_w600_ctrl_078/call_summits`) and wrong at the frame
boundary, where upstream's value is neither clearly shifted nor clearly unshifted.

**Do not attempt a third variant by picking a different constant.** The next move has to be a
measurement of the mapping, not another choice of number:

1. Capture upstream's `pos_array` for this fixture (its first and last elements). `--cutoff-analysis`
   already proves the p/q side agrees, and `MACS3_RS_DUMP_BROAD_LEVELS` gives ours, so the pair pins
   the array origin exactly.
2. With the origin known, decide whether `chunks_at`'s `i == 0` branch wants that origin or the
   literal coordinate, and re-run the **full** corpus before keeping the change -- never the single
   fixture, which both of these attempts passed.

This is the last golden difference; the corpus is otherwise 22454/22456 files with 0 exit-status
mismatches.

### F271: `chunks_at`'s `i == 0` branch is right for 442 cases and wrong for one -- and the
    discriminator is `--nolambda`

F270 said: stop picking constants, measure. Measured.

`chunks_at` mirrors upstream's

```python
above_cutoff_startpos = pos_array[above_cutoff-1]
if above_cutoff[0] == 0:
    above_cutoff_startpos[0] = 0
```

so the question is how often that branch is actually taken. Instrumenting it and replaying the
corpus:

```text
counted-PE cases in corpus : 4860
sampled                    : 473
total i==0 firings         : 292
cases firing >= 1          : 282
cases firing 0             : 191
```

For `sweep/gmini_mfrag_d400_w600_ctrl_222 --nolambda` it fires **exactly once**.

**This accounts for F270 exactly.** 282 of 473 sampled cases is ~60%, which scales to the ~442 cases
F270 broke. So the literal `0` currently in `chunks_at` is **correct for every one of them** -- those
cases passed golden before F270 and stopped passing when the branch was remapped to `coord_shift`.
F270 was not "almost right"; it was wrong for the overwhelming majority of the branch's occurrences,
and the single fixture it fixed was the minority case.

**So the discriminator is not the frame, it is the code path.** The failing fixture is the only
remaining `--nolambda` case, and `--nolambda` is exactly the path where upstream's control side is not
a real pileup:

```python
ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg], dtype="f4")]
```

(`CallPeakUnit.py`, the `no_lambda_flag` branch). The paired array is then walked with a
single-element control, so `pos_array` is derived from the treatment side alone rather than from the
`min(t_p[i], c_p[j])` cursor that normally produces it. That is a plausible reason for the first
above-cutoff index, and therefore for the meaning of the `= 0` special case, to differ on this path
only -- but "plausible" is not "measured", and the corpus has only this one case to say.

**What is now known, and what is not.**

- Known: the branch is hot (60% of counted-PE cases), our literal `0` is correct for all of them, and
  `--nolambda` is the only path where it is not.
- Not known: whether the right fix is a `nolambda`-conditional origin, or whether upstream's
  `pos_array` genuinely differs on that path for a deeper reason that would also show up in a
  `--nolambda` fixture that currently passes by luck.

**Do not implement the conditional yet.** A `cfg.nolambda`-keyed origin is a one-line change that
would be justified by exactly one fixture, and the last two one-line changes aimed at this case each
cost ~450 cases. The measurement that settles it is to capture upstream's `pos_array` (first and last
elements, and the index of the first above-cutoff position) for one *passing* `--nolambda` fixture and
for this one. If the passing fixture's first above-cutoff index is also 0, the conditional is wrong
and something upstream-side differs; if it is greater than 0, the conditional is right and the
remaining question is only why our index disagrees.

The probe (`I0_FIRES`) has been removed; the tree is back to 7203/7685 cases and 22454/22456 files
with 1 differing case and 0 exit-status mismatches.

### F272: `--nolambda` is the one path where upstream's literal `0` is `coord_shift` -- golden parity

Closes the last golden difference, and completes the corpus.

`chunks_at` mirrors upstream's first-above-cutoff special case:

```python
above_cutoff_startpos = pos_array[above_cutoff-1]
if above_cutoff[0] == 0:
    above_cutoff_startpos[0] = 0
```

F271 measured that this branch fires in **282 of 473** sampled counted-PE cases and that our literal
`0` is correct for every one of them; F270 had already established that remapping it globally to
`coord_shift` breaks ~442 cases. The remaining failure was `sweep/gmini_mfrag_d400_w600_ctrl_222
--nolambda`, where upstream reports `chrM 0 200` (summit offset 178) and we reported
`chrM -5000 200` (offset 5178) -- uniformly `coord_shift` low.

The discriminator is `--nolambda`, and it is structural rather than incidental: that is the one path
where upstream has no real control pileup, substituting a one-element array
(`ctrl_pv = [treat_pv[0][-1:], [lambda_bg]]`). The paired array is then walked without a
`min(t_p[i], c_p[j])` control cursor, so the first above-cutoff index -- and therefore the meaning of
the literal `0` -- differs there and nowhere else.

`ChromCall` gained `zero_coord`, set to `coord_shift` on the paired-end path **only when
`cfg.nolambda`**, and `0` everywhere else.

**Validated at corpus scale, not on the target fixture.** `--nolambda` appears throughout the sweep
fixtures, so a change keyed on it is immediately falsifiable by the rest of the corpus:

```text
cases with output    : 7204 (byte-identical: 7204/7204)
output files matched : 22456/22456
invalid-invocation   : 481 cases, exit-status mismatch 0
```

**Every recorded output file in the corpus is byte-identical**, every exit status matches, and nothing
regressed.

**Caveat, recorded rather than glossed.** This is justified by one failing fixture plus the fact that
the other several hundred `--nolambda` cases do not break. It is *consistent* with the structural
argument above, and it is validated far more broadly than a single-fixture fix, but the deeper
question F271 raised -- whether upstream's `pos_array` differs on the `--nolambda` path for a reason
that this conditional merely happens to model -- is not settled. The honest description is "the
smallest change that is correct on all available evidence", not "derived from first principles".

### F273: the golden summary reported 93.7% for a corpus that was at 100%

`oracle/run_golden.py` printed `cases byte-identical : 7204/7685` while the corpus was in fact at
parity. The counter `full` requires `total > 0`, and the denominator was `len(cases)` -- which
includes the **481 invalid-invocation cases that legitimately produce no output files** (425 argparse
rejections at rc 2, 56 runtime rejections at rc 1). Those cases can never be byte-identical because
they have nothing to compare, so they silently depressed the ratio.

```text
cases with recorded files   : 7204
cases with NO recorded files: 481   (425 at rc 2, 56 at rc 1)
```

So the true figure was 7204/**7204**, and the exit statuses of all 481 rejection cases were already
being checked (and matched) by `tot_rc` -- it was simply not visible in the same place.

The summary now reports the three quantities separately:

```text
cases with output    : 7204 (byte-identical: 7204/7204)
output files matched : 22456/22456
invalid-invocation   : 481 cases, exit-status mismatch 0
```

This is the same failure mode as F245 (a positional column read that named the wrong field) and F254
(a diagnostic whose columns did not line up): a summary that is technically accurate and
practically misleading. Worth noting that this one survived a long time precisely because the corpus
*was* short of 100% for most of it -- the number only became obviously wrong once it was right.

---

## F274. `predictd -m 100 200` read only its first value, and silently discarded every peak

`bin/macs3`, `argparser_predictd.add_argument("-m", "--mfold", type=int, default=[5, 50], nargs=2)`.

This is not an upstream behaviour but a divergence found while answering "how would the corpus be
widened past `callpeak`". Recording it here because the mechanism is the interesting part.

`nargs=2` means the two bounds arrive as two separate argv entries. Our `Options` keeps the **first**
value in `values` and the whole list in `lists` (`flags.rs`, `take_action`), and `predictd.rs` read it
with `o.get("mfold")` -- which returns `"100"` for `-m 100 200` -- then split that on `,`:

```rust
let mut it = s.split(',');
let lo = it.next()...;                    // 100
let hi = it.next()...unwrap_or(lo);       // no comma -> falls back to lo
```

So `-m 100 200` became `(100, 100)`. In `PeakModel.build` that makes

```text
min_tags = total * lmfold * peaksize / gsize / 2
max_tags = total * umfold * peaksize / gsize / 2
```

equal, and `naive_call_peaks` keeps a summit only when `v > min_v` **and** (in `__close_peak`)
`summit_value < max_v`. With the bounds equal, no value satisfies both, so every peak was discarded.
On the real CTCF fixture:

```text
MACS3   : Total number of paired peaks: 3803, d = 99, writes predictd_model.R
macs3-rs: "can only find 0 paired peaks", writes nothing
```

`oracle/check_predictd_randsample.sh` passed throughout, because it drives `-m` through the same
reader and therefore agreed with the bug. Confirming the arithmetic was *not* the fault: a probe of
`find_paired_peaks` on the same file reproduced upstream's per-chromosome counts exactly
(chr2: 682 plus / 0 minus; chr8: 425 / 403) and its 3,803 paired centres, so only the values
*reaching* it were wrong.

Fixed by reading `o.get_all("mfold")` and pattern-matching both elements, as `callpeak.rs` already
did. `predictd_model.R` is now byte-identical on the 5 M-read fixture.

Two lessons, both already stated elsewhere and both re-confirmed here:

- **A differential gate that drives the parser under test agrees with itself.** It has to compare
  against upstream *and* pin the parsed value, or it cannot see a value that never arrives.
- **The recorded corpus covered `callpeak` only, so `predictd` had no gate at corpus scale.** The
  per-command gate used one synthetic fixture whose bounds happen not to distinguish `(100, 100)`
  from `(100, 200)`. The first non-callpeak invocation on real data found a command-breaking bug in
  about ten minutes.
