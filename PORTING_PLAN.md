# macs3-rs — Porting Plan

**Objective:** a Rust reimplementation of MACS3 3.0.5 that is a drop-in replacement for
the `macs3` Python package (CLI + library semantics), produces bit-comparable scientific
results, and delivers a large performance and memory win. The fermi-lite assembler is a
deliberate C FFI exception, retained beyond v1.0; a pure-Rust replacement is post-v1.0
scope.

**Status:** this is the original plan of record; implementation is active and substantial
portions are complete. See [the live status](docs/status.md) for measured coverage,
including the remaining compatibility and performance work. Historical milestones below
remain useful as targets, not as a claim that no code exists.

---

## 0. Executive summary

MACS3 is not a slow Python script. Its numeric core is Cython + vectorized NumPy, and
3.0.5 already ships an optimized `PileupV2` path. A naive "translate the .py files to
.rs" effort would produce a *slower* program. Therefore this port is specified as a
**behaviour-preserving redesign**:

| Layer | Upstream | Port strategy |
|---|---|---|
| Orchestration | Python | removed entirely; becomes typed Rust |
| Representation | per-base NumPy/`array` chromosome vectors | run-length-encoded signal tracks |
| Pileup | dense array fill | endpoint event sweep → RLE |
| Local lambda | dense max-reduction over arrays | RLE element-wise max |
| p/q scores | dense per-base tables | RLE + (p,q) key caches + histogram-driven p→q map |
| Parallelism | effectively single-threaded | chromosome-level Rayon parallelism |
| I/O | Python file objects | buffered binary readers/writers |

Correctness is not a byproduct of performance work. It is the **gate structure** of this
plan: each gate is a hard, machine-checkable pass/fail. Performance work is only permitted
after the corresponding correctness gate is green, and every optimization is required to
re-pass all previously-green gates.

The end state is a publishable crate set, a compat matrix, a benchmark suite, and a
documented numerical contract — i.e. a software-engineering result, not a Rust exercise.

---

## 1. Definition of "drop-in replacement"

We commit to a precise, testable claim. For every subcommand in §2, on every fixture in
the regression corpus, for a matrix of option combinations:

1. **CLI identity** — same program name, same subcommand names, same long/short flags,
   same aliases, same defaults, same value parsers, same mutual-exclusion rules, same
   error-exit conditions, same generated output filenames.
2. **Exit codes** — `0` on success, non-zero on validation failure, no panics ever.
3. **Content identity** — output records compared by a genomic-aware comparator
   (§9.6), not by `diff`.
4. **Scientific identity** — peak coordinates, summits, ordering, and p/q/FE values equal
   within the tolerance budget of §4.3.

The one intentional deviation, declared up front and documented in the README:

> `macs3-rs` will not reproduce upstream *human-readable error text* verbatim. It will
> reproduce upstream *error conditions* (same input ⇒ same failure, before any output
> file is created).

### 1.1 Success metrics (release v1.0 definition)

| # | Metric | Target |
|---|---|---|
| M1 | Peak coordinate identity vs MACS3 (narrowPeak/broadPeak/bedGraph) | 100% on all golden + differential fixtures |
| M2 | Summit identity | 100%; max allowed shift 0 bp on all fixtures except the explicitly-declared HMMRATAC set |
| M3 | qscore max abs error | ≤ 1e-6 |
| M4 | pscore max abs error | ≤ 1e-9 |
| M5 | Retained read counts / duplicate counts / d | exact integers |
| M6 | Output byte identity (formatting) | 100% for `*.xls`, `*_summits.bed`, `*_peaks.narrowPeak`, `*.bdg` |
| M7 | CLI acceptance parity | 100% on the shared flag matrix; MACS3-valid invocations must succeed, MACS3-invalid must fail |
| M8 | End-to-end speedup, `callpeak` SE/PE/broad, 1–100M reads | ≥ 3× wall-clock (median, same host) |
| M9 | Peak RSS | ≤ 50% of MACS3 on all end-to-end benchmarks |
| M10 | No command slower than MACS3 | 0 regressions, or documented justification |
| M11 | Test suite | zero ignored/failing tests on supported platforms |
| M12 | Reproducibility | 1 thread ≡ N threads, byte-identical output |

M12 is absolute, not aspirational: parallelism may change runtime, never output.

---

## 2. Scope — the 14 subcommands

Compatibility target is the full MACS3 3.0.5 CLI surface:

`callpeak`, `bdgpeakcall`, `bdgbroadcall`, `bdgcmp`, `bdgopt`, `cmbreps`, `bdgdiff`,
`filterdup`, `predictd`, `pileup`, `randsample`, `refinepeak`, `callvar`, `hmmratac`.

### 2.1 Compatibility level per subcommand

| Subcommand | Level | Notes |
|---|---|---|
| `bdgopt`, `bdgcmp`, `cmbreps`, `bdgdiff` | **Exact** | pure run-length arithmetic |
| `bdgpeakcall`, `bdgbroadcall` | **Exact** | + documented `max_gap` boundary semantics |
| `filterdup`, `randsample` | **Exact** | `randsample` requires MACS3-compatible RNG |
| `pileup` | **Numeric** | SPMR/BDG float formatting |
| `predictd` | **Numeric** | cross-correlation tolerance |
| `callpeak` (narrow) | **Numeric** | + Exact coordinates |
| `callpeak` (broad) | **Numeric** | |
| `refinepeak` | **Numeric** | |
| `callpeak -f FRAG` | **Numeric** | |
| `hmmratac` (inference) | **Exact** | fixed exported model |
| `hmmratac` (training) | **Oracle-checked numeric** | Gaussian and Poisson self-training are implemented; see G13 evidence in `docs/status.md` |
| `callvar` | **Exact on gated fixtures** | Both explicit and `--fermi` paths are tested; existing fermi-lite FFI is retained beyond v1.0 scope |

### 2.2 Non-goals (v1.0)

Python/R bindings, AnnData/H5AD public API, bigWig output, GPU/SIMD hand-tuning,
distributed execution, new compression formats, GUI. AnnData support is *permitted* as an
optional feature but must never block CLI parity.

---

## 3. Phase 0 — Freeze the Oracle (M0)

**Nothing is written in Rust before this phase is complete.** A port without a
byte-level oracle is a rewrite, and a rewrite cannot be validated.

### 3.1 Pin the reference

Create `oracle/ENV.lock` recording: `MACS3_VERSION`, upstream git commit SHA, Python
minor version, exact `numpy`/`scipy`/`pandas`/`scikit-learn`/`hmmlearn` versions, OS, and
`ldd`-level detail for the compiled Cython extensions. Build MACS3 from that commit in a
hermetic venv (`uv`/conda lockfile committed). The oracle image is reproducible in CI.

### 3.2 Build the fixture corpus

`tests/fixtures/` is generated, not hand-written, plus a curated set:

- **Real**: MACS3's own test corpora (e.g. the 200K CTCF ChIP/control pair) plus
  additional public ChIP, ATAC/BAMPE, broad histone, FRAG/scATAC, and HMMRATAC inputs.
- **Synthetic (differential)**: a seeded generator producing miniature genomes
  (1–5 contigs, 1e4–1e6 bp), Gaussian/narrow/broad/bimodal peak shapes, Poisson read
  counts, controlled duplicate rates, controlled local-bias spikes, barcode multiplicity,
  and degenerate edge cases. Every synthetic fixture is fully described by a seed, so the
  corpus is regenerable and growable.
- **Edge-case**: single read on a contig; reads with no common chromosomes between
  treatment and control; zero control reads; all-duplicates; `d` at the lower limit;
  control lambda exactly equal to treatment; contiguous plateau p-scores; equal-height
  summits; summit at interval boundary; fragment length 0/1; 2bp fragment; fragment
  spanning a contig end; contig with 1 base; 4,294,967,296 offsets (u32 overflow probe).

Target: ≥ 400 fixtures at M0 exit, growing monotonically.

### 3.3 The oracle harness

`scripts/oracle.py` runs upstream MACS3 and captures **both final outputs and every
intermediate representation** the reference implementation materialises:

parsed reads · filtered reads · duplicate counts · predicted `d` · model arrays
(`plus_line`, `minus_line`, `xcorr`, `ycorr`, `alternative_d`) · treatment pileup ·
control pileup · per-scale local lambda tracks · merged local lambda · p-score track ·
q-score track and the p→q table · candidate intervals · merged intervals · summits ·
final records.

Everything is dumped in a **stable, documented text or fixed binary format**
(`docs/fixtures.md`) with checksums, so upstream can be upgraded and the corpus re-derived
deterministically.

> The purpose is to convert *"I think my algorithm is MACS3"* into
> *"my stage 6 diverges from MACS3 at coordinate 18,353,022."*

### 3.4 CLI flag matrix

`tests/fixtures/cli_matrix.json` enumerates every flag, its aliases, its default, and a
set of invocations per subcommand — including invalid ones that must fail identically.
This matrix is generated by introspecting upstream's argparse configuration and then
hand-audited.

### Gate G0 — Oracle frozen

**Pass criteria (all mechanical):**
- `oracle` image builds from a committed lockfile on a clean machine.
- `macs3 --version` matches `MACS3_VERSION`; commit SHA recorded.
- Corpus ≥ 400 fixtures; every fixture has a recorded SHA-256 and ≥ 1 oracle artefact.
- 100% of the corpus runs to completion on the oracle (any upstream crash becomes a
  documented known-failure with a reproduction, not a silent skip).
- `cli_matrix.json` covers 100% of upstream flags (auto-diffed against argparse at runtime).
- `docs/numerical-contract.md` (§4) written and reviewed.
- `docs/compatibility.md` (§2.1) committed.

**Exit rule:** no Rust code merges until G0 is green.

---

## 4. Numerical contract

Written *before* the first algorithm, because compatibility is a floating-point problem
long before it is a data-structure problem.

### 4.1 Type policy

```rust
type Coord  = u64;   // 0-based genomic offset
type Count  = u64;
type Len    = u64;
type Score  = f64;
```

Chosen deliberately:

* `u64` for coordinates so contigs > 4.29 Gbp and > 2^32 offsets are impossible to
  overflow.
* Counting in `u64` avoids the 32-bit overflow class of bug that upstream had to patch.
* `f64` internally, but **upstream `f32` behaviour is emulated wherever it is observable**
  (a documented list of call sites). Silent "improvement" to `f64` is a compatibility
  regression, because a threshold crossing changes a genomic interval and therefore an
  entire peak.

### 4.2 Required content

`docs/numerical-contract.md` documents, per stage: input types, accumulation type,
promotion points, final output type, rounding, formatting, NaN/Inf handling, sort
stability and tie-breaking, and equality semantics (exact vs epsilon).

### 4.3 Tolerance budget

| Quantity | Tolerance | Rationale |
|---|---|---|
| Coordinates, summits, lengths, counts | exact (0) | integers; any difference is a bug |
| `d`, fragment sizes | exact (0) | integers |
| p-score | ≤ 1e-9 abs | single probability evaluation |
| q-score | ≤ 1e-6 abs | cumulative-min over a monotone p→q map |
| fold enrichment / SPMR | ≤ 1e-6 rel | division + formatting |
| correlation arrays (`predictd`) | ≤ 1e-6 abs | direct O(n²) accumulation order |
| HMMRATAC posteriors | ≤ 1e-4 | EM-conditioned; biological level |

Tolerances are **budgets, not targets**. The implementation goal is exact agreement;
`macs-compare` reports the observed max error so regressions are visible.

### 4.4 Deliberate format rules

* BED-family coordinates: 0-based, half-open `[start, end)`.
* XLS output: MACS3's own 1-based-ish convention preserved exactly (this is a known
  divergence between XLS and BED in upstream and must be replicated, not fixed).
* Float formatting: replicate upstream's `repr`/format spec, including trailing-zero
  behaviour, because M6 requires byte identity.
* Sorting: reproduce upstream's stable-sort and tie-break behaviour explicitly; Rust's
  `sort_unstable_by` is **forbidden** anywhere output order is observable.

### Gate G1 — Numerical contract accepted

**Pass:** document exists, covers all 9 pipeline stages, has a machine-checkable
`tolerance_budget` table consumed by `macs-compare`, and has been validated by
recomputing 3 hand-derived examples from upstream source (Poisson tail, binomial
inverse-CDF, AFDR step function) in exact arithmetic.

---

## 5. Workspace architecture

```
macs3-rs/
├── Cargo.toml                    # workspace
├── crates/
│   ├── macs-core/                # coords, intervals, errors, chromosomes, interning
│   ├── macs-io/                  # readers: BED/BEDPE/FRAG/SAM/BAM/bedGraph; writers
│   ├── macs-track/               # SingleEndTrack, PairedEndTrack, FragTrack
│   ├── macs-rle/                 # SignalTrack<T> run-length primitive + algebra
│   ├── macs-pileup/              # endpoint sweep, single-end extension, PE fragments
│   ├── macs-stats/               # Poisson, binomial, inverse binomial, normal
│   ├── macs-score/               # p-score, q-score, FE, LLR, p→q histogram
│   ├── macs-model/               # PeakModel: strand profiles + cross-correlation
│   ├── macs-peaks/               # narrow caller, broad caller, summits, deconvolution
│   ├── macs-bedgraph/            # bdgopt/cmp/cmbreps/diff/peakcall/broadcall
│   ├── macs-hmmratac/            # fragment EM, 4-channel signal, 3-state HMM
│   ├── macs-callvar/             # BAM query, candidate, likelihood, assembly, VCF
│   ├── macs-cli/                 # clap surface, 14 subcommands, orchestration
│   └── macs-compare/             # genomic-aware differential comparator (dev tool)
├── tests/
│   ├── fixtures/                 # inputs (git-LFS or generated)
│   ├── golden/                   # oracle artefacts + manifests
│   ├── regression/               # one tiny fixture per historical/upstream bug
│   ├── differential/             # generated + property tests
│   └── fuzz/
├── benches/                      # criterion micro + e2e
├── oracle/                       # pinned upstream env, harness, flag matrix
└── docs/                         # numerical-contract, compatibility, architecture, fixtures
```

### 5.1 Crate rules

1. `macs-core`, `macs-rle`, `macs-stats` have **zero I/O** and zero heavy dependencies.
   They are the audit surface; everything else may be audited through them.
2. Dependencies are additive-only. No transitive dependency may influence numeric results
   that we also compute ourselves (this is why statistics are hand-written, §7).
3. Every crate exposes a `#![forbid(unsafe_code)]` unless a documented `unsafe` island
   exists (BAM decoding via `noodles` is the only expected exception).
4. `macs-io` depends on a `trait AlignmentReader`; `noodles` is the default backend and
   `rust-htslib` is an optional backend used only if exact BAM/BAI random-access
   semantics for `callvar` require it. Core algorithms never see a backend type.
5. Public API is documented, semver-stable, and re-exported from `macs` (facade crate).
6. Every algorithmic crate has a `par` module boundary: the deterministic single-threaded
   implementation is the reference; the parallel implementation is a thin wrapper that
   must produce identical output (property-tested, §9.3).

### 5.2 Core representations

```rust
struct Interval { start: Coord, end: Coord }          // [start, end)
struct ChromId(u32);
struct Genome { names: Vec<Box<[u8]>>, lookup: HashMap<Vec<u8>, ChromId>, lens: Vec<Len> }
struct Peak { start, end, summit: Coord, pileup, pscore, qscore, fe: Score }
struct Run<T>  { end: Coord, value: T }              // run is [prev_end, end)
struct SignalTrack<T> { chrom: ChromId, runs: Vec<Run<T>> }
```

Chromosome names are **interned**; no `String` is attached to a genomic record. Reads are
stored as sorted primitive vectors per chromosome and strand; fragments as
`(start, end)` pairs; barcodes as `BarcodeId(u32)` into a dictionary.

---

## 6. Pipeline architecture (the whole design in one diagram)

```
   input files
        │
        ├── parse ─────────► interned chromosomes, sorted compact tracks
        │                        │
        │                   filterdup / randsample / exclusion / max-count
        │                        │
        │                   predictd ─────────► d  (SE only)
        │                        │
        │                   pileup (endpoint sweep, RLE)
        │                        │
        │    control ──► per-scale lambda tracks (d, slocal, llocal)  RLE
        │                        │
        │                   element-wise max  ──────────►  λ_local   RLE
        │                        │
        │                   treatment vs λ_local  ─────►  p-score   RLE
        │                        │
        │                   p-score histogram ────────►  p→q map
        │                        │
        │                   q-score   RLE
        │                        │
        │                   threshold ─ merge ─ min-len ─ summit ─ deconvolve
        │                        │
        └────────────────────────► writers (xls, narrowPeak, summits, bdg, model)
```

**Invariants, all property-tested:**

* `∫ pileup dpos == Σ fragment_length × weight` exactly.
* Adjacent equal-valued runs are always coalesced ⇒ canonical RLE form.
* Runs are strictly increasing in `end`; no empty runs.
* No stage ever materialises a per-base genome array unless an algorithm provably
  requires it (documented exception list, reviewed).
* `RLE(merge_adjacent(parse(x))) == parse(x)` valuewise for every track type.
* 1 thread output `==` N threads output, byte for byte.

---

## 7. Milestones and gates

Each milestone: **deliverable → gate → hard pass criteria → what becomes permanently
immutable afterwards.** A gate is only marked green with the stated command output
attached to the PR. Regressions against a green gate block merge (CI-enforced).

---

### M1 — Core + RLE algebra  → **Gate G2**

**Deliverable:** `macs-core`, `macs-rle`, `macs-stats` (Poisson pdf/cdf/sf/`-log10`,
large-λ path, binomial pdf/cdf, inverse binomial CDF, normal cdf/pdf), RLE algebra
(`zip`, `map`, `max`, `sub`, threshold, coalesce, interval-query, iterators).

**G2 pass criteria:**
- `macs-stats` validated against an independent reference (`scipy` in the oracle env) on
  ≥ 10,000 parameter points incl. λ ∈ {0, 1e-8, 1, 10, 100, 1e3, 1e5, 1e7} and extreme
  tails, max abs error ≤ 1e-12 vs. scipy, ≤ 1e-9 vs. upstream `MACS3.Signal.Prob`.
- Inverse binomial CDF: monotone, round-trips (`icdf(ppf(x)) == x` for all tested x).
- RLE: 10,000 proptest cases; canonical-form and zip-commutativity properties pass.
- `cargo fuzz run` (≥ 60 s) on `bed` and `rle` targets: no panics, no OOM.
- 100% `cargo llvm-cov` line coverage for `macs-stats` and `macs-rle`.

**Becomes immutable:** coordinate/type policy, RLE canonical form, statistical functions.

---

### M2 — Parsers  → **Gate G3**

**Deliverable:** BED, BEDPE, FRAG, bedGraph, narrowPeak, peak-BED, barcode list, SAM,
BAM, BAMPE. Legacy: BOWTIE, ELAND, ELANDMULTI, ELANDEXPORT, AUTO detection.

**G3 pass criteria:**
- Record-level identity: `MACS3 parsed record == macs3-rs parsed record` (chrom, start,
  end, strand, name/order, multiplicity, flags) on 100% of the ≥ 400 corpus fixtures,
  plus the dedicated edge-case set: reverse strand, soft/hard clip, I/D, secondary,
  supplementary, unmapped, mate-on-other-chrom, duplicate, QC-fail, unmapped-mate,
  zero-length, MAPQ 0, unmapped-with-mapped-mate, >2^32 offsets, >2^31 file offsets,
  truncated gzip, corrupt records, empty file, header-only file, CRLF, no trailing
  newline, extra columns, missing columns.
- Every malformed input produces a typed error, **never a panic**; `catch_unwind` fuzz
  harness confirms 0 panics over ≥ 1e6 mutated inputs.
- Round-trip property: for every track type, `serialize(parse(s))` then `parse` is
  value-identical.
- BAM backend parity: `noodles` and `htslib` backends produce identical records on all
  BAM fixtures.
- Throughput baseline recorded (not yet a gate): BED ≥ 200 MB/s/core, BAM ≥ 60 MB/s/core.

**Becomes immutable:** the record abstraction (`Alignment`, `FragmentRecord`,
`BedGraphRecord`) used by every downstream stage.

---

### M3 — Tracks  → **Gate G4**

**Deliverable:** `SingleEndTrack`, `PairedEndTrack`, `FragTrack` (separate types, not one
generic), with append, finalize/sort, totals, genomic coverage, mean fragment length,
duplicate detection, random sampling, down-sampling, exclusion, barcode filter, max-count.

**G4 pass criteria:**
- Totals (`total`, `total_span`, `avg_fragsize`, `max_fragsize`, per-chromosome stats)
  exact integer match on 100% of the corpus.
- Duplicate counts and `max_dup` (auto mode via inverse binomial CDF) exact match.
- Sampling: `--keep-dup auto` output identical on all fixtures; seeded
  `randsample --down-sampled-size`/percentage output identical for seeds {1, 42, 12345,
  2^31-1} and sizes {0, 1, N/2, N-1, N, N+1}; percentage {0, 0.5, 100}.
  (Where NumPy RNG order cannot be reproduced, we implement an explicit
  MACS3-compatible RNG/shuffle rather than accepting divergence.)
- `FragsTrack` handles multiplicity and barcode identity; property: sum of counts equals
  the multiplicity-weighted total.
- Memory property: no per-read heap allocation after finalize (measured via allocation
  counting in a test binary).

**Becomes immutable:** track types and the duplicate-definition semantics (explicitly
*MACS3's* definition, not BAM flags).

---

### M4 — Pileup engine  → **Gate G5**

**Deliverable:** endpoint-sweep pileup producing RLE; single-end extension (forward,
reverse, bidirectional, `shift`, `extsize`), paired-end true fragments, scaling
(`--scale-factor`, SPMR), `--max-count`, barcode selection, `--B`.

**G5 pass criteria:**
- Integral conservation (property, 10,000 proptest cases, exact integer equality).
- Pileup values and run boundaries identical to the oracle on 100% of the corpus, in
  both 1-thread and N-thread modes.
- SE-vs-PE: paired-end pileup uses real fragment intervals and ignores estimated `d`
  (verified by a fixture where `d` is deliberately absurd).
- No per-base allocation; peak RSS for a 100M-read pileup < 1.5× input record size.
- Benchmark: ≥ 1.5 GB/min of BAM input per core (recorded as baseline, becomes a
  release gate at G13).

**Becomes immutable:** pileup semantics.

---

### M5 — bedGraph suite  → **Gate G6**

**Deliverable:** `bdgopt`, `bdgcmp`, `cmbreps`, `bdgdiff`, `bdgpeakcall`, `bdgbroadcall`.

**G6 pass criteria:**
- Byte-identical output on 100% of the corpus for all six commands across the full flag
  matrix (all `--method`/`-m` modes for `bdgcmp`: `ppois`, `qpois`, `subtract`, `logFE`,
  `FE`, `logLR`, `slogLR`, `max`; all `--operator` modes for `cmbreps`; all
  `--transformation` modes for `bdgdiff`).
- `bdgpeakcall` gap boundary triple tested explicitly for every fixture: gap =
  `max_gap-1` (merged), `max_gap` (merged), `max_gap+1` (not merged) — the exact
  upstream comparison operator is nailed down and documented.
- `bdgbroadcall`: level2 max gap = 4 × level1 max gap, verified; `broadPeak` and
  `gappedPeak` emitted from one internal object; invalid cutoff ordering rejected
  identically.
- File merging semantics (`--cat`, multiple inputs, differing score fields, overlapping
  records, zero-score records, negative coordinates rejected) all byte-identical.
- Micro-benchmark: track zip ≥ 1e9 run-elements/s/core.

**Strategic value:** this milestone delivers user value early, exercises the entire RLE
stack, and de-risks the hardest part (q-values) in a small surface.

**Becomes immutable:** bedGraph I/O formatting and the peak-merging comparators.

---

### M6 — p/q score engine  → **Gate G7** *(highest-risk gate)*

**Deliverable:** p-score RLE, p-score histogram, p→q map, q-score conversion, score caches.

This is the single most compatibility-sensitive component. Upstream does **not** apply a
generic BH correction; it accumulates the number of base pairs at each p-score and
constructs the p→q mapping from those genomic lengths. A 1e-12 error here shifts
thousands of peaks.

**G7 pass criteria:**
- The p→q table is golden-tested independently of `callpeak` (dump table, diff).
- Table is monotone non-decreasing in q for descending p; pscore ≥ qscore pointwise
  (property test, 10,000 cases).
- q-score RLE identical to oracle on 100% of the corpus, max abs error ≤ 1e-6, **and**
  identical peak sets when thresholded.
- p-score cache is a pure memoisation: with cache disabled vs enabled, output is
  byte-identical (test both ways on every fixture).
- Cache-key canonicalisation is total: no two semantically different inputs share a key
  (property test using random (λ, k) pairs).
- Edge cases: zero local lambda everywhere; λ = 0 at a position; k = 0 with λ > 0;
  plateaus; single-base genome; control larger than treatment.
- Caches never introduce quantization that upstream does not perform (explicitly
  tested: quantization off by default; any future quantization requires a gate).

**Becomes immutable:** the p→q map algorithm.

---

### M7 — predictd / PeakModel  → **Gate G8**

**Deliverable:** strand profiles, naive ± peak finding, +/− pairing, normalisation,
**direct O(n²) cross-correlation**, smoothing, local-max selection, `d` choice,
model fitting, `--shift`/`--extsize` bypass, `MODEL.r` output.

Port the upstream sequence literally, in order: naive + peaks → naive − peaks → pair
compatible + / − peaks (accept only `0.5 < plus_count/minus_count < 2`) → aggregate
strand profiles → normalise → cross-correlate → smooth → local maxima → drop `d ≤ d_min`
→ choose strongest surviving `d`. Reproduce `peaksize = 2 × bandwidth`, tag expansion
of 10, and the MFOLD thresholds exactly.

**G8 pass criteria:**
- `plus_line[]`, `minus_line[]`, `xcorr[]`, `ycorr[]`, `alternative_d[]` arrays match the
  oracle within §4.3 **element-by-element** on the real CTCF corpus and 100% of
  synthetic fixtures with ≥ 100 paired peaks.
- Predicted `d` exact integer match on 100% of the corpus.
- The "< 100 paired peaks ⇒ refuse to model" condition is reproduced exactly, including
  the error condition.
- `MODEL.r` byte-identical.
- Deliberate *not* done: FFT cross-correlation. Vectors are small; direct O(n²) buys
  numerical parity far more cheaply. FFT is a post-v1.0 optimization that must
  reproduce the same arrays within tolerance.

**Becomes immutable:** the model algorithm.

---

### M8 — callpeak (narrow)  → **Gate G9**

**Deliverable:** full orchestration: filterdup, exclusion, model, scaling, local lambda
(`λ = max(λ_global, λ_d, λ_slocal, λ_llocal)`), p/q, peak calling, `--call-summits`,
`--nolambda`, `--nomodel`, `--extsize`, `--shift`, `--bw`, `--mfold`, `--slocal`,
`--llocal`, `--mfold` validation, `--scale-to small|large|random`, `--SPMR`, `--BAMPE`,
`-f BED|BAM|BEDPE|BAMPE|FRAG`, all output writers.

**G9 pass criteria:**
- Identical retained read counts, `d`, treatment pileup, λ_local segmentation, peak
  coordinates, summit coordinates, peak ordering (M1–M5, M12) on **100% of the corpus**,
  for both `--call-summits` off and on.
- Byte-identical `*.xls`, `*_peaks.narrowPeak`, `*_summits.bed`, `*_model.r`,
  `*_treat_pileup.bdg`, `*_control_lambda.bdg`.
- Validation parity: identical accept/reject for the full invalid-invocation matrix
  (e.g. `--slocal < d`, `--mfold` length ≠ 3, mismatched chromosome naming, `PE format
  with SE options`, `genome size` not in table and not numeric, missing BAI, `--keep-dup`
  out of range, ...). Every failure occurs **before any output file is created**.
- Normalization semantics: paired-end treatment coverage vs. control background handled
  exactly as upstream, verified by dedicated fixtures where the two differ.
- `--nolambda` and `--nomodel` paths each fully covered.
- Threading: output byte-identical for `--threads 1/2/4/8/16/32`.

**Becomes immutable:** the callpeak pipeline.

**Gate command:** `oracle/run_golden.sh`. It replays every recorded
`tests/golden/<group>/<fixture>/<variant>/command.json` argv with `macs3-rs` and
byte-compares every output file the recording lists, which is what "byte-identical
`*.xls`, `*_peaks.narrowPeak`, `*_summits.bed`, `*_model.r`,
`*_treat_pileup.bdg`, `*_control_lambda.bdg`" means. Three header lines are
normalised because they quote the invocation and cannot match a replay
(`# Command line:`, and the `# ChIP-seq file = [...]` / `# control file = [...]`
lines of the `# ARGUMENTS LIST` block); everything else, including all counts,
`d`, scale factors and every peak row, is compared verbatim.

`oracle/run_peak_e2e.sh` remains the fast loop, and it compares **coordinates and
summits only**. That is a real gate but a strictly weaker one: a score that moves
in the fourth decimal, a wrong header count, or a `-B` bedGraph with the right
values at the wrong breakpoints all pass it. Use `run_golden.sh` for the G9
pass criteria.

---

### M9 — callpeak (broad) + `--call-summits` deconvolution  → **Gate G10**

**Deliverable:** strong/weak two-level caller, `broadPeak`/`gappedPeak`, plus the
standalone summit-deconvolution module and the pathological summit fixture set.

**G10 pass criteria:**
- Broad peak calls identical to the oracle on 100% of the corpus (coordinates and
  boundaries), for the full matrix of broad cutoff ratios.
- `broadPeak` and `gappedPeak` byte-identical.
- Summit deconvolution identical on a dedicated pathological set: single flat summit;
  two equal maxima; two nearby peaks; three nearby peaks; maximum at boundary; short
  supporting signal; sub-cutoff gap; plateau; maximum with an asymmetric shoulder; a
  maximum whose supporting region is clipped by a peak boundary. Each case has an
  individually-reviewed expected result derived from upstream source, not from observed
  output (this is the only place where "whatever upstream does" is not an acceptable
  spec — we record *why* upstream does it, in a code comment citing the upstream line).
- q-score recomputation after summit changes matches oracle.

**Becomes immutable:** summit algorithm.

---

### M10 — pileup / randsample / filterdup CLIs + refinepeak  → **Gate G11**

**Deliverable:** these subcommands' full flag surfaces, plus `refinepeak` reusing the
track + interval-query primitives (no duplicated pileup infrastructure).

**G11 pass criteria:**
- Byte-identical output on 100% of the corpus for `pileup`, `randsample`, `filterdup`.
- All `pileup` options exercised: `--extsize`, `--shift`, `--direction`, `--bidir`,
  `--scale-factor`, `--SPMR`, `--barcode`, `--max-count`, `--zero-based`, `--keep-dup`.
- `refinepeak`: identical refined peak regions and `refined_summits.bed` on all fixtures;
  verified also on inputs that were *not* produced by MACS3 (candidate regions from an
  external caller), since this is a common real-world use.

**Becomes immutable:** the auxiliary subcommand surfaces.

---

### M11 — FRAG / single-cell  → **Gate G12**

**Deliverable:** FRAG parsing with multiplicity, barcode interning + dictionary,
`--max-count`, barcode-aware single-cell calling, optional AnnData/H5AD I/O (feature
gated, non-blocking).

**G12 pass criteria:**
- FRAG corpus parity: byte-identical pileups and peak calls, including multiplicity
  weighting and barcode filtering.
- Scale: 500k fragments × 50k barcodes completes in < 60 s and < 2 GB RSS on one core,
  with chromosome-level parallelism scaling linearly (measured 1/4/8/16 threads).
- `--max-count` boundary behaviour exact (count == N included, count == N+1 excluded;
  verified against upstream, whichever it is, with the upstream line cited).

**Becomes immutable:** FRAG semantics.

---

### M12 — HMMRATAC  → **Gate G13** (split into two gates)

**Deliverable:** fragment-size EM (short / mono-nucleosomal / di / tri), class
assignment, 4-channel genomic signal, 3-state HMM, disk-backed decode store, output
writers. **No generic HMM crate** — 3 states × 4 features is small enough to implement
forward, backward, posteriors, Viterbi, and Baum-Welch directly, which gives complete
control over numerics, initialisation, randomness, and convergence.

**G13a — inference parity (Exact):**
- Export a trained upstream model (start probs, transition matrix, means, covariances /
  Poisson λ, bin size) to JSON; the Rust binary loads exactly that file.
- Posterior probabilities match upstream within 1e-4; state assignments and output
  regions identical. This isolates the *inference* implementation from *training*.

**G13b — training (oracle-checked):**
- Gaussian and Poisson self-training are implemented and checked against a fresh run of
  the pinned oracle by `oracle/check_hmm_training.py`. On the yeast 500k fixture, training
  regions and feature rows match (feature values within 2.53e-12); fitted model values
  are within 2.52e-10 (Gaussian) and 3.87e-12 (Poisson), and accessible-region intervals
  match exactly for both types. These are measured fixture results, not a claim of all-
  flags or all-corpus parity.
- Decode store uses buffered binary blocks in a temp file, not JSON: bounded RSS on a
  1M-candidate run (< 500 MB), verified.

---

### M13 — callvar  → **Gate G14**

**Deliverable:** split into `bam_query`, `pileup`, `candidate`, `likelihood`, `assembly`,
`realign`, `vcf`.

**Strategy — two stages, deliberately:**

1. **Bridge (v0.x):** thin C FFI to upstream `fermi-lite`, so callvar is functionally
   complete early. A de novo assembler must never be on the critical path of the peak
   caller.
2. **Pure Rust (post-v1.0 scope):** read correction, overlap discovery, unitig graph, unitig
   extraction, local realignment — reimplemented, each step gated by a unitig-level
   oracle fixture set (graph shape + contig sequences), which is a *much* easier
   equivalence test than end-to-end VCF.

**G14 pass criteria:**
- SNV records: identical REF/ALT, QUAL within 1e-3 relative, genotype concordant.
- INDEL records: identical REF/ALT and left-normalisation on gated fixtures; the
  fermi-lite assembler bridge is retained, with a pure-Rust assembler remaining
  post-v1.0 scope.
- VCF header (##fileformat, ##contig, ##INFO, ##FILTER, ##FORMAT, column order) byte-
  identical; sample-column formatting identical.
- Known-issue list (`docs/callvar-known-issues.md`) is empty of unreproduced upstream
  behaviour.

---

### M14 — Performance programme  → **Gate G15**

Performance is a **gated deliverable**, not an afterthought. It begins only after G15
depends on all prior gates being green.

**Work list, in order (profiling-driven, not intuition-driven):**
1. Parsing (memchr-based tokenisation, zero-copy fields, `unsafe` islands where measured)
2. Track finalize/sort (radix sort by coordinate; `sort_unstable` allowed *only* for
   internal state, never for output order)
3. Pileup endpoint generation + merge (arena reuse, branchless sweep, PGO)
4. Control-λ generation (windowed RLE max, single pass)
5. p-score cache (open-addressing; typed key packs to 16 bytes; per-chromosome)
6. q-value histogram (sort-free bucketing where it does not change the map)
7. Chromosome-level Rayon parallelism everywhere it is deterministic
8. Output buffering (large `BufWriter`, hand-rolled float formatting to match upstream)
9. HMM inner loops (`f64`→`f32` only where the tolerance budget permits and it is
   *verified*, never assumed)
10. Assembly in callvar

**G15 pass criteria (release gate):**
- M8, M9, M10, M12 from §1.1 met on the benchmark matrix.
- Zero commands slower than MACS3 (§1.1 M10) on any benchmark.
- `cargo flamegraph` profile archived per optimization PR; no optimization lands without
  a before/after benchmark and an unchanged golden-test result.
- Criterion microbenchmarks with recorded regression thresholds (5% on tracked benches).
- LTO + `codegen-units=1` + PGO in the release profile; release build reproducible.

---

## 8. Test strategy (the six layers)

| Layer | What | Tool | Gate |
|---|---|---|---|
| L1 Unit | statistics, interval algebra, RLE, parsers, scoring, merging | `#[test]`, `approx_eq` | every PR |
| L2 Property | invariants of §6, cache totality, round-trips, merge-equivalence | `proptest` | every PR |
| L3 Differential | macs3-rs vs. pinned MACS3 on the whole corpus, stage by stage | `macs-compare` + oracle harness | every PR (nightly: full matrix) |
| L4 Golden | real biological datasets, byte-comparison | CI artifacts | every PR |
| L5 Regression | one tiny fixture per bug ever found, incl. mined upstream fixes | `tests/regression/` | every PR |
| L6 Fuzz | parsers, RLE, HMM decoder, VCF writer, CLI arg handling | `cargo-fuzz` (nightly) | every PR on changed targets |

### 8.1 L3 is the load-bearing layer

Stage-by-stage differential testing is mandatory. It is not enough that the final peaks
match — we assert equality at *each* pipeline stage against the oracle artefact, so a
divergence localises to exactly one stage. CI reports, per fixture per stage:
`equal` / `mismatch@coordinate` / `count-only-mismatch` / `missing`.

### 8.2 Regression discipline

Every bug found becomes:
1. a minimal 2-line fixture,
2. a test that fails without the fix,
3. a comment on the fix explaining the *semantic* reason, not just the symptom.

Regression fixtures are additionally mined from upstream's own changelog fixes (summit
scoring, peak padding, HMM refinement, FRAG pileups, integer overflow) — those are known
divergence-prone areas and each gets a dedicated regression case.

### 8.3 Determinism testing

A dedicated CI job runs the *entire* corpus at `--threads 1` and `--threads 8` and
requires byte-identical outputs. Failure = hard stop.

### 8.4 Fuzzing budget

Nightly ≥ 4 CPU-hours: 1h parsers, 1h RLE/score algebra, 1h HMM decode store, 1h CLI
fuzzing (random flag combinations must error cleanly, never panic), plus corpus
mutation fuzzing (bit-flips, truncations, coordinate corruption) to prove graceful
degradation.

---

## 9. Tooling built for the port

### 9.1 `oracle` harness (in-repo, Python)
Runs pinned upstream, dumps intermediates + final outputs + checksums.

### 9.2 `macs-compare` (Rust dev tool)
Genomic-aware comparator. Reports, per fixture:

```
Fixture: ctcf_200k
  retained reads     MACS3 412,331   Rust 412,331        IDENTICAL
  d                  MACS3 178       Rust 178            IDENTICAL
  pileup             8,402 runs      8,402 runs          IDENTICAL
  lambda_local       3,117 runs      3,117 runs          IDENTICAL
  p-score            max|Δ| 3.1e-15                    PASS (≤1e-9)
  q-score            max|Δ| 2.7e-12                    PASS (≤1e-6)
  peaks              4,927 / 4,927    coords 4,927/4,927  PASS
  summits            4,925 / 4,927    max shift 1 bp      PASS
  output bytes                                                    IDENTICAL
```
For HMMRATAC: Jaccard, precision, recall, boundary/summit displacement distributions.
This tool is the project's primary development instrument and is itself released.

### 9.3 Threading-equivalence harness
Property test + CI job asserting 1-thread ≡ N-thread output on a rotating fixture subset.

### 9.4 Benchmark suite
Criterion micro-benchmarks (parser throughput, track finalize, 1M-fragment pileup, RLE
zip, Poisson p-score, p→q conversion, peak merge) + end-to-end benchmarks (callpeak SE /
PE / broad / FRAG, pileup, all bdg commands, hmmratac, callvar) recording wall, CPU,
peak RSS, bytes in/out, records/s. Environment pinned; results published per release.

### 9.5 `docs/compatibility.md`
Machine-readable compat matrix: command × option × fixture-set × status × max observed
deviation. Generated by CI, not maintained by hand.

### 9.6 Comparator semantics
Comparison is on parsed records, not text, so that floating formatting differences do not
mask coordinate differences. Interval comparison is exact; the comparator reports
*first divergence* with genomic context.

---

## 10. CI matrix

| Job | Trigger | Content |
|---|---|---|
| `fast` | every push | fmt, clippy `-D warnings`, `m4` (L1+L2, < 5 min) |
| `golden` | every push | L4 golden byte-compare on ~15 small fixtures |
| `corpus` | merge to main / nightly | L3 full corpus (400+ fixtures, all stages, oracle-backed) |
| `threads` | nightly | 1 vs 8 thread determinism over the whole corpus |
| `fuzz` | nightly | cargo-fuzz budget §8.4 |
| `perf` | nightly + release tags | benchmark suite, regression thresholds, RSS limits |
| `bench` | release tags | full end-to-end matrix vs. pinned MACS3, published |
| `msrv` | weekly | minimum supported Rust version |
| `platform` | weekly | Linux x86_64/aarch64, macOS (aarch64), optional Windows |

Pinned MACS3 runs in its own hermetic container; the Rust side never depends on Python.

---

## 11. Release plan

| Version | Scope | Gate |
|---|---|---|
| `0.1` | `macs-core`, `macs-rle`, `macs-stats` (published to crates.io, no CLI) | G2 |
| `0.2` | `macs-io` + `bdgopt`, `bdgcmp`, `cmbreps`, `bdgdiff` | G3, G6-partial |
| `0.3` | `pileup`, `filterdup`, `randsample`, `bdgpeakcall`, `bdgbroadcall` | G4, G5, G6, G11 |
| `0.4` | `predictd` | G8 |
| `0.5` | `callpeak` narrow | G9 |
| `0.6` | `callpeak` broad, `--call-summits`, `refinepeak` | G10, G11 |
| `0.7` | FRAG / single-cell | G12 |
| `0.8` | `hmmratac` (inference and Gaussian/Poisson training oracle-checked) | G13 |
| `0.9` | `callvar` via fermi-lite bridge | G14-bridge |
| **`1.0`** | **drop-in replacement**: 14/14 commands, M1–M12 met | G15 |
| `1.1+` | optional pure-Rust assembler, FFT correlation, library ergonomics, bindings | post-v1.0 |

**v1.0 announcement is gated on:** 14/14 commands, 100% corpus parity at declared
tolerances, ≥ 3× end-to-end speedup on the headline matrix, ≤ 50% peak RSS, byte-identical
golden outputs, and published benchmark methodology.

**Deprecation policy for Python:** MACS3 is not modified. Drop-in-ness is proven by the
comparator over the same corpus, and the README ships a migration table
(`macs3 …` → `macs3-rs …`) plus a `macs3-rs-compat` environment shim that re-exports the
CLI names so existing shell scripts and Nextflow/Snakemake pipelines run unchanged.

---

## 12. Risk register

| # | Risk | Likelihood | Impact | Mitigation | Early-warning signal |
|---|---|---|---|---|---|
| R1 | p/q numerical parity | High | Critical | G7 as a standalone gate before any caller; golden p→q table; histogram-based map mirrors upstream construction, not textbook BH | q-score divergence on any fixture |
| R2 | Summit/boundary exactness | High | Critical | G10 pathological fixture set; upstream-line citations; summit deconvolution isolated as its own module | any 1 bp summit shift |
| R3 | λ_local semantics (d/slocal/llocal, PE vs SE normalization) | High | High | Oracle intermediates for every stage; dedicated PE normalization fixtures | segmentation mismatch in G8 |
| R4 | `predictd` exact `d` | Medium | High | Literal port; direct O(n²) correlation; array-level comparison | `xcorr` element mismatch |
| R5 | BAM/BAMPE edge cases | High | High | Dual-backend cross-check; heavy fuzzing; borrow `noodles` semantics | backend divergence |
| R6 | HMMRATAC training divergence | High | Medium | Split inference (exact) from training (biological); declare deviation in writing | Jaccard < 0.98 |
| R7 | callvar + fermi-lite | High | Medium | FFI bridge first; pure Rust last; never on the critical path | unitig graph mismatch |
| R8 | Upstream changes during the port | Medium | Medium | Pinned commit; re-derivation script; compat matrix is per-commit | oracle rebuild diff |
| R9 | "Compatible" claims unsupported by evidence | Medium | Critical (reputational) | Every claim is a CI artifact; `docs/compatibility.md` is generated, not asserted | any hand-edited matrix |
| R10 | Rust underperforms Cython+NumPy | Medium | High | G15 is a hard release gate; profiling-driven work list; abandon the "no per-base arrays" rule only with documented evidence | first e2e benchmark < 1× |
| R11 | Scope creep into 14 commands too early | Medium | Medium | Milestone ordering; v1.0 is a hard gate, not a date | 3+ gates open simultaneously |
| R12 | Numerical drift from third-party crates | Medium | High | Statistics hand-written; no default-`f64`-assuming distribution crates; numerics audit in code review | value differs from scipy/upstream |

---

## 13. Effort and sequencing

Rough person-months, single engineer, with the assumption that the oracle harness and
comparator are built first and used religiously:

| Phase | Months | Cumulative gates |
|---|---|---|
| M0 oracle + contract | 1.5 | G0, G1 |
| M1 core/RLE/stats | 1.5 | G2 |
| M2 parsers | 2.0 | G3 |
| M3 tracks | 1.5 | G4 |
| M4 pileup | 1.0 | G5 |
| M5 bedGraph suite | 1.5 | G6 |
| M6 p/q engine | 1.5 | **G7** |
| M7 predictd | 1.0 | G8 |
| M8 callpeak narrow | 2.0 | **G9** |
| M9 broad + summits | 1.0 | G10 |
| M10 aux CLIs + refinepeak | 1.0 | G11 |
| M11 FRAG | 1.0 | G12 |
| M12 HMMRATAC | 2.0 | G13 |
| M13 callvar (bridge) | 2.0 | G14 |
| M14 performance programme | 2.5 | G15 |
| **Total to v1.0** | **~21.5** | |

Critical path: M0 → M1 → M2 → M3 → M4 → M6 → M8 → M14. Everything else can be
parallelised or deferred. Pure-Rust fermi-lite and FFT correlation are explicitly
post-v1.0.

---

## 14. What "amazing performance" concretely means

Claimed only with evidence, published per release, measured on a pinned host with
`taskset` pinning and identical inputs:

1. **Throughput** — BAM/BED parsing and pileup within 2–4× of raw I/O bandwidth;
   peak calling is then memory-bandwidth-bound, not compute-bound.
2. **Algorithmic** — RLE pipeline avoids the O(genome) work that dominates upstream;
   expected several-fold wall-clock reduction on large genomes, plus dramatically lower
   memory (RSS bounded by the largest single chromosome, not the whole experiment).
3. **Parallelism** — chromosome-level parallelism that upstream does not have at all;
   near-linear scaling until memory bandwidth saturates.
4. **Elision** — 1-thread and N-thread paths share code; no duplicated implementations
   to drift.
5. **Measured, not modelled** — every number in the announcement comes from the
   `bench` CI artifact, and the benchmark methodology (host, pinning, input hashes,
   thread counts, warm/cold) is published so others can reproduce it.

---

## 15. Definition of Done — v1.0

The project is done when, for every fixture in the ≥ 400-fixture corpus, across the
documented flag matrix:

```bash
macs3    callpeak -t chip.bam -c input.bam -f BAM -g hs -n test -B -q 0.01 --threads 8
macs3-rs callpeak -t chip.bam -c input.bam -f BAM -g hs -n test -B -q 0.01 --threads 8
```

yield, per `macs-compare`:

* identical retained read counts, duplicate counts, `d`;
* identical treatment pileup and identical λ_local segmentation;
* identical p-scores (≤ 1e-9) and q-scores (≤ 1e-6);
* identical peak coordinates, summit coordinates, and peak ordering;
* byte-identical `*.xls`, `*_peaks.narrowPeak`, `*_summits.bed`, `*_model.r`;
* identical output at `--threads 1` and `--threads 32`;

and `macs3-rs` does it with **≥ 3× lower wall-clock and ≤ 50% peak RSS** across the
release benchmark matrix. HMMRATAC training is checked against the pinned oracle as
documented under G13; full release claims still depend on the complete matrix.

All 14 subcommands present. Zero panics. Zero failing or ignored tests. Published
benchmarks. Published compatibility matrix generated by CI. Drop-in for existing
`macs3` command lines.
