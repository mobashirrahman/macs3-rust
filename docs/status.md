# macs3-rs port status

## Current snapshot — 2026-10-05

This snapshot supersedes older status notes below where they conflict. The project
has live implementations across the CLI, but command coverage is uneven and it is
not yet a full all-flags, all-corpus replacement.

- **Final validation:** 903 workspace release tests pass with zero failures and
  zero ignored tests; formatting and strict Clippy pass. The golden replay passes
  7,204 cases and matches all 22,456 files, with zero exit-status mismatches across
  481 rejection cases. `oracle/check_thread_invariance.py --limit 0` replays all
  7,685 recorded cases at 1 and 32 threads with zero differences.
- **HMMRATAC self-training is implemented for Gaussian and Poisson models.**
  `oracle/check_hmm_training.py` runs fresh pinned-oracle training and decoding checks.
  On the yeast 500k fixture, training regions and feature rows match (feature-row max
  deviation 2.53e-12); model max deviations are 2.52e-10 (Gaussian) and 3.87e-12
  (Poisson). Accessible-region intervals and state outputs match exactly for both
  model types. `oracle/check_hmmratac.py` separately gates inference and signal output.
- **The hmmratac model-parameter gap is measured against upstream's own noise floor,
  not assumed.** Training a fresh Gaussian model three times — twice with MACS3,
  once with `macs3-rs` — and comparing all 77 parameters:

  | comparison | max abs | max rel |
  |---|---|---|
  | MACS3 run 1 vs MACS3 run 2 | 1.76e-13 | 1.93e-12 |
  | MACS3 run 1 vs macs3-rs | 2.07e-10 | 8.01e-12 |

  So upstream is *not* bit-reproducible against itself at the 1e-12 level, and our
  gap is about 4× that floor in relative terms — a real residual, not noise, and
  not the 2e-10 figure quoted in earlier notes as though it were upstream's own
  variation. It has no effect on output: `fresh_cutoff_analysis.tsv` is
  byte-identical, and `oracle/check_hmmratac.py` gates the decoded regions,
  states and digested tracks as byte-identical (Jaccard 1.0000). Reproduce with
  `--modelonly --jump 1.5` on `test/yeast_500k_SRR1822137.bedpe.gz`, run from
  the oracle source root with `PYTHONPATH` set to it, comparing the three
  `fresh_model.json` files numerically.
- **Input support includes pooled single-end and BAM reads without an index**, plus
  legacy parsers. Utility commands pool all supplied inputs, detect AUTO per file, and retain the
  first-file tag-size estimate. ELAND and ELANDEXPORT match fresh comparisons;
  compiled upstream ELANDMULTI and BOWTIE tag-size inference failures are reproduced.
- Seeded multi-contig `randsample` has a documented ordering difference: upstream
  seeds NumPy once and iterates a Python set of chromosomes. Two fresh upstream
  processes with the same seed selected different rows (50,752/100,000 overlap);
  fixing `PYTHONHASHSEED=0` made upstream repeatable. Rust uses sorted chromosomes.
- **Counted scATAC FRAG sampling is implemented with upstream MT19937 and
  count-weighted, ties-even quotas.** Fresh scATAC training arrays and state outputs
  match exactly; model deviation is below 3e-12 and accessible-base Jaccard is 1.0.
- The fresh upstream command suite produces all 163 expected artifacts. The
  corrected intermediate recorder and Rust capture pass all 11 observable stages:
  3,745 numeric leaves and 23 exact leaves. This comparison is mandatory in CI.
- Explicit callpeak `--tempdir` controls temporary spools and preserves upstream
  rejection of nonexistent directories; omitted flags use the platform temp path.
- The two formerly ignored oracle tests now use committed Poisson and correlated
  HMM captures and run without Python.
- **CI provisions a private pinned MACS3 oracle** and runs differential gates for
  Gaussian and Poisson HMM training, HMM inference, and all three callvar modes: no assembly, automatic assembly, and forced assembly.
- **Performance: both targets partly met.** Re-measured with `scripts/bench_real.py`
  (real CTCF 5M treatment + 5M control reads, median of 3): Rust is 2.3-4.9x faster
  than pinned upstream across twelve workloads, at 0.17-1.01x its peak RSS. The
  release targets are at least 3x and at most 0.5x. Speed is met on 10 of 12
  (`pileup` 2.3x and `randsample` 2.6x miss); memory is met on 7 of 12 (`callpeak
  -B` 0.56x, `--broad` 0.52x, `randsample` 0.70x, `bdgpeakcall` 0.89x, `bdgopt`
  1.01x miss). The per-workload table is in the README. The older 2.3-3.0 GB Rust
  RSS figures below are historical baseline measurements.
- **Differences found on the 5M-read data, fixed.** `bdgpeakcall` restarted peak
  numbering per chromosome (upstream numbers continuously); `bdgbroadcall` emitted
  peaks on chromosomes with no level-1 peak and ignored `-o` for the name prefix;
  `bdgopt -m p2q` differed by 1e-5 on 227,084 of 2,577,274 rows (the `-log10(N)`
  term, the q value before the clamp and the `pre_q` seed are f32 upstream);
  `bdgcmp -m logFE|ppois|qpois` now exits 1 with no output file where upstream
  raises. All five 5M-read `pileup`/`bdg*` outputs are byte-identical; workspace
  tests are at 903. `hmmratac --save-training-data` now prints values as Python's
  `repr` does (ties-to-even on the last digit), taking upstream's `cmdlinetest`
  to 154 of 163 files byte-identical — a real measurement, but run ad hoc against
  the oracle checkout and not committed as a harness, which is why it is absent
  from the README's reproducible table. `filterdup` output still has the same rows
  in a different chromosome order.
- **Continuation fixes verified:** `callpeak -p` is honored; fresh 5M-read
  `--call-summits` runs now match all eight files across model and `--nomodel`
  modes, after reproducing Haswell OpenBLAS's fused dot-product reduction.
  `bdgcmp` fractional pseudocounts and fold enrichment use the pinned NumPy's
  float32 arithmetic; `bdgopt multiply` does too. `refinepeak` preserves same-start
  peak order and upstream's read cursor. Separate negative option values parse
  correctly, including the negative no-control `--llocal` case that previously
  panicked. These checks establish parity for the tested cases, not every CLI edge.
- **CI portability:** live-reference tests resolve the provisioned source and
  interpreter and skip explicitly when absent. All 903 release tests pass both
  with the pinned oracle and in a snapshot containing neither `.oracle` nor
  `ENV.provisioned` (five live-reference checks report skips there). NumPy-based
  probes use the pinned interpreter, and the Savitzky–Golay differential has a
  checked-in Rust probe. Tracked fuzz build output has been removed. GitHub
  Actions runs all six layers on every push; the six-layer split (L1 unit,
  L2 invariants, L3 differential, L4 golden, L5 regression, L6 fuzzing) is green
  on `master` at commit `3f1bd248`.
- **The C fermi-lite assembler bridge is intentionally retained.** Replacing it with
  a pure-Rust assembler remains post-v1.0 scope; the existing bridge is tested in CI.
- `oracle/check_real_summit_bytes.py` reruns source-tree CTCF SE, BEDPE, and BAMPE
  examples against the fresh pinned oracle. The checked narrowPeak, summit, treatment,
  and control outputs match byte-for-byte; the XLS command-line header is normalized
  because it contains the different executable and temporary output paths. The source
  SE `*_model.r` is also byte-identical after reproducing NumPy/OpenBLAS Haswell
  dot-product reduction. Fresh predictd/randsample checks pass 6/6. Do not read the
  historical command and gate rows below as a current all-14-command parity claim.

## Historical snapshot — 2026-10-01

The detailed gate and test tables below preserve the prior status snapshot and its
then-current findings. Newer work can supersede those rows; in particular, their
HMM self-training, callvar, command-coverage, and test-count entries are historical.
Use the current snapshot above for current high-level status.

---

## Oracle

| Item | State |
|---|---|
| MACS3 source | pinned, `c544319` (v3.0.5), `oracle/ENV.lock` |
| Build | provisioned by `oracle/provision_oracle.sh`, into `.oracle/venv` (gitignored) |
| Reproducible build | `oracle/provision_oracle.sh`, pinned by `oracle/ENV.lock` |
| `macs3 --version` | `macs3 3.0.5` |

Note: `MACS3/fermi-lite/lib` is a git submodule and the build **fails** without
it (`fermi-lite/ksw.c` includes `lib/x86/sse2.h`). `provision_oracle.sh` fetches
it. (`scripts/build_oracle.sh` is the older manual variant, kept for working
outside `.oracle/`.)

---

## Gates

| Gate | Scope | State | Evidence |
|---|---|---|---|
| **G0** | hermetic oracle (ENV.lock, >=400 fixtures, flag matrix, intermediate capture) | **HERMETIC** | `oracle/ENV.lock` + `oracle/verify_oracle_clean.sh` enforce commit `c5443190e3edfeb301cc94acf450e2b2c026a223`, a clean `git status`, no instrumentation hooks and matching `.so` checksums; `macs-stats::oracle_tree_is_unmodified` runs it in CI (F139). 425 fixtures, 3 golden variants, 315-row auto-derived flag matrix, `dump_stages.py` captures the intermediates **that are observable at all** without touching the oracle -- the signal tracks are cdef attributes of `CallerFromAlignments` and are obtained by having the engine write them itself with `-B`. Still open: release-grade golden regeneration, and the control pileup / p-score / q-score tracks, which have no reachable upstream surface. |
| **G1** | Numerical contract | **green** | `docs/upstream-findings.md` (145 findings), this file |
| **G2** | `macs-core` + `macs-rle` + `macs-stats` | **green** | 4115 bit-exact oracle vectors; 0 clippy warnings |
| **G3** | `macs-io` parsers | **green** | BED/BEDPE/FRAG/bedGraph: 54/54 fixture inputs parse record-identically to upstream |
| **G4** | `macs-track` (SE/PE/Frag, dup, sampling) | **partial** | **660 vectors bit-exact** (dup filtering + downsampling); NumPy RNG 11/11 vs NumPy; BAM/SAM construction + `print_to_bed` outstanding |
| **G5** | `macs-pileup` | **green** | 1776 bit-exact oracle vectors |
| **G6** | `macs-bedgraph` (BedGraph.py) | **LARGELY COMPLETE** | `BedGraphTrackI` ported: reader, `add_loc`, `merge_regions`, `overlie` (6 funcs), `extract_value`, `apply_func`, `p2q`, `writer` -- all bit-exact (F124). `make_ScoreTrackII_for_macs` added for `bdgcmp` (F137). Still missing `cutoff_analysis`, `extract_value_hmmr`, and `bedGraphTrackII`'s array-backed `call_peaks`/`refine_peaks` variants |
| **G7** | `macs-score` (p-score, p→q) | **green** | 23 p→q tables bit-exact vs `pq_reference`; 11 real MACS3 cutoff-analysis tables checked |
| **G8** | `macs-model` (PeakModel.py) | **PORTED; `predictd` live** | `PeakModel.build` ported end to end: naive pileup -> summits -> paired-peak centres -> strand profiles -> normalise -> cross-correlate -> smooth -> `d`. **F136**: `d` and `alternative_d` are exact vs upstream; `*_model.r` matches line-for-line (`p`, `m`, `xcorr`, `altd`) except a ~1e-12 `ycorr` float tail (same length, same argmax). Required NumPy-exact reductions: pairwise `mean`/`std`, forward `convolve`, `linspace` lag axis, `correlate` lag convention. Still open: the model path in `callpeak` itself (it currently requires `--nomodel`) |
| **G9** | `callpeak` narrow + call-summits + broad | **narrow and broad measured by `oracle/run_golden.sh`** | The per-chromosome core and the single-end/paired-end pipelines live in the library (`macs-peaks/src/callpeak.rs`); `crates/macs-cli/src/commands/callpeak.rs` is a thin driver over them, so the golden gate measures the shipped code. **Measured by byte-comparing every recorded output file** -- `oracle/run_golden.sh` replays each `tests/golden/<group>/<fixture>/<variant>/command.json` argv and diffs the results: **`default` variant 267/425 cases byte-identical** (181 before this round); **all 7685 recorded cases 4172/7685**, 13434/22456 files, 56 exit-status mismatches (was 2906 / 10035 / 654). That gate found and fixed **fourteen** defects: F152 (no control *file* is not `--nolambda`: 100+ `*_noc_*` fixtures produced no peaks at all), F153 (`# tag size` is the truncated mean of the **first ten** tags), F154 (summit p/q recomputed from the chunk's paired values; control indexed by cursor, not by end), F150/F151 (`-B` writes the paired union walk; `--format FRAG` is a counted `PETrackII` with a weighted treatment sweep, weighted pooled input, centred `d`-wide control windows and forced-off duplicate filtering), F159 (the `broadPeak`/`gappedPeak` writers did not exist, and broad mode passed an **empty** p->q table so every broad peak reported `-log10(q) = 0`), F160/F161 (`--shift`, `--keep-dup auto`), F162/F163 (`--nolambda` keyed on the scale lists, and the control count line not gated on `--keep-dup`), F157 (`--mfold` is `nargs=2`, so `--mfold 3 20` exited 2). Remaining: **2237 `--format FRAG`** cases -- the weighted control lambda still differs in the 5th significant figure, which moves the reported peak start -- **486 paired-end** (4th-decimal p/q drift, F149) and **365** across everything else. `run_peak_e2e.sh` compares coordinates only and reports 288/288 for `default`; that is strictly weaker evidence and is kept as the fast loop. |
| **CMDS** | CLI command coverage | **13 of 14 live; `macs-hmmratac` CLI live (F144/F145); `callvar` driver live** | LIVE (dispatch-confirmed, byte-identical vs oracle): `callpeak`, `bdgopt`, `cmbreps`, `filterdup`, `randsample`, `pileup`, `bdgpeakcall`, `bdgbroadcall`, `refinepeak`, **`bdgcmp`** (8 methods), **`bdgdiff`** (F137), and **`callpeak` with `-f BEDPE`/`BAMPE`/`FRAG`/`BAM` (F148/F149)**, and `predictd` (F136: `d`/`altd` exact, `*_model.r` matches except a ~1e-12 `ycorr` tail). `callvar` **driver live** (VCF header byte-identical to the oracle, exit-status parity) but its variant-calling kernel is not written, so it refuses with exit 1 before creating any output file. `hmmratac` (**live for the inference path** -- fragments in, EM, weight mapping, 4 digested pileups, fold-change track, candidate regions, bin extraction, HMM posteriors, and the likelihood/state/accessible-region writers; `*_cutoff_analysis.tsv` byte-identical, `*_model.json` byte-identical, inference pinned against hmmlearn to a worst total-variation distance of 1.4e-6. `--model` supplies the model; **self-training is the declared deviation** (F143), and the EM down-sample permutation is a second, narrower gap (F145)). 703 tests green |
| **G10** | `callpeak` broad + summits | **partial** | `--call-summits` maxima + `enforce_peakyness` implemented and SG-filtered against upstream; broad-peak calling not started |
| **G11** | `pileup` / `randsample` / `filterdup` / `refinepeak` | not started | — |
| **G12** | FRAG / single-cell | not started | — |
| **G13** | `hmmratac` | **CLI live (inference path); self-training + EM RNG are the gaps** | `crates/macs-hmmratac` ports the deterministic core (weight mapping, digested pileups, bin extraction, `Regions`, `generate_states_path`, `save_accessible_regions`, EM) and `crates/macs-cli/src/commands/hmmratac.rs` drives it. **F143**: inference exact vs hmmlearn (worst total-variation distance 1.4e-6) after reproducing hmmlearn's *inconsistent* matrix orientation. **F144**: `*_cutoff_analysis.tsv` byte-identical after fixing the dropped FRAG count column, the `pos_array[-1]` wrap that makes `cutoff_analysis` silently report nothing when the leading run is above the cutoff, the `f32`/`f64` ladder step, and `summary`'s broken `pre_p`. **F145**: NumPy `SeedSequence` reproduced and verified, but NumPy's MT19937 array seeding is not, so the EM down-sample permutation differs -- exact with `--no-fragem`, not otherwise. Remaining: the hmmlearn Baum-Welch self-training (the declared deviation), `PeakIO.randomly_pick` (Python `random`, not NumPy) for `--maxTrain`, BAMPE input (htslib), and `--barcodes`/blacklist coverage |
| **G14** | `callvar` | **byte-identical, both `-F` modes** | Header, `VariantStat` (bit-identical on 410/410 oracle cases incl. 23 error paths, F213), `PosReadsInfo` (F214), `PeakVariants` (F215), `variant_bq_by_ref_pos` (F216), `get_REFSEQ` + `RACollection` + the calling loop (F217/F218) are all ported and pinned. fermi-lite `r53` and MACS3's own `swalign.c` are vendored byte-identically and bridged through a five-function FFI (F219/F220); the assembly chain (`verify_alns`, `remap_reads_with_unitigs`, `add_to_unitig_list`, `build_unitig_collection`, `UnitigCollection`) is ported (F221/F222). **`-F off` 22/22 and `-F auto` 16/16 records identical** against the pinned oracle on upstream's own `callvar_testing` fixtures, checked by `oracle/check_callvar.sh` in both modes. The last 6 records were a branch-structure bug, not arithmetic (F224): upstream holds the no-assembly variants unwritten and then chooses between dropping the peak, falling back, revisiting only reference-biased positions, or re-calling the whole peak. |
| **G15** | performance programme | not started | — |

---

## Test suite

```
$ cargo test
209 tests passing, 0 failing, 0 clippy warnings, cargo fmt clean
```

| Layer | Count | Notes |
|---|---|---|
| L1 unit | 173 | `macs-core`, `macs-stats`, `macs-rle`, `macs-pileup`, `macs-score`, `macs-io`, `macs-compare` |
| L3 differential (unit level) | 5914 | 4115 statistics + 1776 pileup vectors, **bit-exact** vs MACS3 3.0.5 |
| L3 parser records | 1492 files | BED, BEDPE, FRAG, bedGraph: **record-identical** to upstream's own parsers |
| L3 pre-computes | 23 checks | `--cutoff-analysis` peak counts/lengths vs an independent transcription; cutoff ladder bit-identical to `np.arange` |
| L3 pairing / local lambda | 12 checks | `macs_peaks::pair_treat_ctrl` and the scale merge vs an independent transcription; `--nolambda` extent verified against a real run |
| L3 SG filter | 5 checks | `macs_peaks` Savitzky-Golay + `maxima` match upstream bit-for-bit on spike/ramp signals |
| L3 track vectors | 660 | dup filtering + `sample_percent`/`sample_num`, **record-identical** to upstream |
| L3 NumPy RNG | 11 checks | `macs_stats::NumpyRng` stream and shuffles match NumPy 2.5.3 exactly |
| L3 peak coordinates | 6 fixtures | `oracle/run_peak_e2e.sh`: peak count/boundaries/summits vs golden XLS, both score mechanisms; best case 10/12 boundaries, **summit criterion unmet** |
| L3 end-to-end pileup | 100 fixtures | `oracle/run_e2e.sh`: real pipeline vs upstream `--bdg`; **97 pass, 0 fail**, 3 skipped as upstream rejections |
| L3 pscore vs upstream `get_pscore` | 117 vectors | bit-exact, worst abs diff 0 (F47) |
| L3 p/q tracks | 82 tests | `macs_peaks::cal_pscore`/`cal_qscore` incl. the truncating observed-count cast and AND-combined criteria |
| L3 pq-table | 23 | **bit-exact** vs `oracle/pq_reference.py` (compiled `poisson_cdf`) |
| L3 real-MACS3 table | 11 | `--cutoff-analysis` tables structurally checked |
| L3 stage dump | 10 stages | `oracle/dump_stages.py`; reads pre/post filter, d, slocal/llocal, scaling, **treatment pileup**, merged lambda, q table, candidate peaks, summits |
| L2 property | 0 | **gap** — `proptest` is declared but not yet used |
| L4 golden | 7685 runs | `tests/golden`, regenerated over all 425 fixtures; rc 0: 7204, rc 2 (bad arity): 425, rc 1 (upstream crash, F27/F28): 56; no Rust side yet |
| L5 regression | 0 | **gap** |
| L6 fuzz | 6 targets, nightly job, 40000 execs each clean | **done** (see F208: the layer had never been built) |

### Regenerating the corpus

The corpus is only evidence if regenerating it is a no-op, so both properties are
checked (F26):

```
$ python3 oracle/gen_fixtures.py --out tests/fixtures --sweep 400
$ python3 oracle/gen_fixtures.py --out /tmp/a --sweep 400
$ PYTHONHASHSEED=1 python3 oracle/gen_fixtures.py --out /tmp/b --sweep 400
$ diff -r /tmp/a /tmp/b        # empty
```

Golden runs are regenerated per group so a long sweep can be resumed:

```
$ python3 oracle/run_oracle.py --macs3 "$(command -v macs3)" \
      --fixtures tests/fixtures --golden tests/golden --only sweep
```

`--only` takes a fixture directory name, `--variants` a comma-separated variant
list, `--limit` an upper bound on fixtures.

### Differential

The parser-record differential is reproducible:

```
$ cargo build -p macs-io --bin macs-io-dump
$ python3 oracle/compare_records.py
1492 files identical, 0 files differing, 1492 inputs total
```

It drives `oracle/dump_oracle_records.py`, which instantiates upstream's own
`BEDParser` / `BEDPEParser` / `FragParser` / `bedGraphTrackI` rather than
re-implementing them, so a disagreement localises to the Rust port. Comparison is
sorted on every column as a number: a text sort orders `100` before `49` and
reports differences that do not exist.

---

## Intermediate capture

`callpeak` writes no intermediates, so `oracle/dump_stages.py` drives the real
pipeline in-process and observes it. It rebinds the module globals that
`callpeak_cmd` and `PeakDetect` resolve at call time — `load_tag_files_options`,
`PeakDetect`, `CallerFromAlignments` — so the arithmetic, ordering and control
flow are upstream's, and only the observation is ours. Stages that live in C
attributes are read from files upstream itself writes (F24).

| stage | source | captured |
|---|---|---|
| reads before dedup | track contents at loader return | yes |
| reads after dedup | `PeakDetect` construction | yes |
| duplicate limit, redundant rate | from the run log | yes |
| d | `# d = N` in the XLS header | yes |
| slocal / llocal | "Range for calculating regional lambda" header line | yes |
| treatment pileup (per position) | `--bdg` bedGraph written by the engine | yes |
| control lambda (per position) | `--bdg` bedGraph written by the engine | yes |
| q-score table | `--cutoff-analysis` file written by the engine | yes |
| candidate intervals | `PeakIO.write_to_xls` / `write_to_summit_bed` | yes |
| final intervals + summits | the golden `*.xls` / `*_summits.bed` | yes |
| model arrays | `*_model.r` golden output | yes (narrow mode only) |
| **per-scale lambda arrays** (d / slocal / llocal, separately, per position) | locals in `CallerFromAlignments.call_peaks` | **no** |
| **merged lambda array** | same | **no** |

The two missing rows are locals inside a Cython method. Capturing them needs
either an upstream patch or a re-derivation in the harness; the re-derivation
route risks reproducing the port's own bug rather than upstream's, so it should
be done only with an explicit cross-check. Everything the local bias feeds is
still pinned indirectly, via the pileups, the XLS header windows, and the final
peak coordinates.

```
$ python3 oracle/dump_stages.py tests/fixtures/se_model/realistic se \
      --out /tmp/stages -- --nomodel --extsize 200 --bdg
```

---

## Golden corpus

`tests/fixtures` — **425** generated fixtures in 11 groups, fully specified by
`manifest.tsv` and regenerable from `oracle/gen_fixtures.py --sweep 400`.

The corpus must regenerate byte-identically, and does: two consecutive runs differ
in 0 lines, and generation under `PYTHONHASHSEED=1` and `=4242` produces identical
trees (see F26 — it did not, until `hash()` seeding was replaced).

| Group | Fixtures | Purpose |
|---|---|---|
| `sweep` | 400 | deterministic grid over genome x mode x depth x width x control, for corpus breadth |
| `frag_counts` | 1 | FRAG counts that truncate when stored in an `unsigned short` (F22) |
| `se_basic` | 1 | two Gaussian peaks plus background |
| `se_shapes` | 4 | narrow / broad / bimodal / local-bias peak shapes |
| `se_dup` | 5 | duplicate rates 0, 0.05, 0.25, 1.0, 5.0 |
| `se_edge` | 5 | contig ends, one-sided, disjoint chromosomes, shallow, no control |
| `se_model` | 3 | many-site fixtures for the fragment model (see F17) |
| `pe_basic` | 3 | Gaussian fragments, ATAC-short, nucleosome ladder |
| `frag_basic` | 1 | FRAG with barcodes and multiplicities |
| `tiny`, `onechrom` | 2 | 1-2 contig sanity cases |

`tests/golden` — 475 runs (24 fixtures × ~20 variants), 395 with exit code 0.

Non-zero exits are **expected and recorded**, not skipped:

* `mfold_bad_arity` × 24 — argparse usage error, exit 2 (F16).
* 3 fixtures × all variants — upstream crash: `check_names` raises `TypeError`
  when treatment and control share no chromosome names (F15), and
  `__call_peaks_w_control` raises `ZeroDivisionError` on a degenerate control.
  The Rust side must fail *cleanly* on the same condition; the traceback is not
  a contract.

Every variant's command line, exit code, stderr tail, output file list and
per-file SHA-256 are in `tests/golden/**/command.json`, and
`tests/golden/summary.json` indexes them.

---

## The comparator

`macs-compare` walks the golden tree and a `macs3-rs` tree, matches files by
name within each run, and reports per file: peak counts, exact coordinate
matches, summit agreement, max deviation per numeric column against the
tolerance budget, and byte identity. It also compares **exit codes** and the
**set of produced files**, so "we emitted three files instead of five" is caught.

Validated both ways:

```
$ macs-compare --golden tests/golden --ours tests/golden --only se_shapes
  files compared    59     files clean 59
  RESULT: PASS

# inject a single 1e-5 q-value drift on one record of one file
$ macs-compare --golden tests/golden --ours /tmp/negctl2 --only .../nomodel_shift
  peak count        MACS3 13   macs3-rs 13
  coordinates       13 / 13 exact
  summits           13 / 13, max shift 0 bp
  q-value           max |delta| 1.0000000000509601e-5
  VIOLATION         q-value |delta| 1.0e-05 exceeds 1e-6
  RESULT: FAIL
```

---

## Upstream findings

`docs/upstream-findings.md` records 17 behaviours that a translation would get
wrong. The ones that most shaped the design:

| # | Finding | Consequence |
|---|---|---|
| F5 | pv arrays are **right-endpoint** indexed: `v[i]` is the value on `[p[i-1], p[i])` | reading it the other way shifts the whole signal by one interval |
| F5a | coincident start/end events are **correct** (I initially had this wrong) | the breakpoint list is coarser, not the values |
| F1 | the log10 Poisson path is **accurate** to its 5-decimal rounding (also initially wrong) | p-scores are a 1e-5 lattice |
| F14 | that path costs `O(k + lambda)` per call | the cache is not an optimisation, it is the reason this is tractable |
| F3 | p-score cache keys on the **f32 bit pattern** of lambda | lambda must be carried as `f32` everywhere |
| F7 | `pduplication` truncates each pmf element to `f32` *before* use | 8e-4 relative error if missed |
| F11 | `binomial_pdf` **does not terminate** for tiny `b` | we bound the loop; documented deviation |
| F15 | `check_names` crashes instead of reporting no common chromosomes | we fail cleanly on the same condition |
| F16 | `--mfold` is `nargs=2`, not 3 | CLI arity must match exactly |
| F17 | `PeakModel` **cannot** be trained on a miniature genome | G8 uses MACS3's own CTCG data |
| F19 | `FWTrack`/`PETrack` arrays are pre-sized to `buffer_size` and only truncated by `finalize()` | reading `locations` before `finalize()` sees `2 * buffer_size` records, mostly zeros; both parsers and the oracle harness call `finalize()` first |
| F20 | the FRAG **barcode is not retained** on the fragment track | only the per-fragment count survives, so a FRAG differential can compare coordinates and counts but not barcodes |
| F18 | the **lowest** p-score bucket is always forced to q = 0 | the AFDR walk's last value is discarded; harmless but not what the formula says |

---

## Where callpeak stands

`oracle/run_golden.sh`, release build, whole recorded corpus:

| scope | value |
|---|---|
| recorded invocations | 7685 |
| **exit-status mismatches** | **0** |
| cases byte-identical | 7079 / 7685 |
| output files matched | 22245 / 22456 |
| `default` variant byte-identical | 416 / 425 |
| workspace tests | 579 passing, 0 failing |
| clippy warnings | 0 |
| `cargo fmt --check` | clean |

All 7685 invocations reproduce upstream's exit status; every rejection path writes
**no output file**; nothing panics.

### Test layers

All six layers the release definition names now have content and are wired into
`.github/workflows/ci.yml`.

| layer | what runs it | size |
|---|---|---|
| L1 unit | `cargo test --workspace` | 579 tests, clippy 0, fmt clean |
| flag audit | `oracle/audit_accepted_flags.py` | every `store_true` flag in the matrix swept with/without; 0 with no observable effect |
| predictd/randsample | `oracle/check_predictd_randsample.sh` | 6 identical, 0 differing (exit status + output tree) |
| hmmratac emission | `full_covariance_emission_matches_hmmlearn` | pinned to hmmlearn's `log_multivariate_normal_density` |
| hmmratac EM | `oracle/check_hmmratac.py` | EM means/stddevs identical to the oracle. **Jaccard 0.6585 on yeast500k vs 0.98 required** (F228/F229). Localised: the EM fit matches upstream exactly, the bin grid (`s//binsize*binsize`), the states-path argmax tie-break and narrowPeak columns 1-5 are all correct, and summit offsets are exact (0 bp over 7385 identically placed regions). The **digested signals** diverge at scale -- `max_abs_dev` 2.14 on `short`, 1.0 on `mono` -- and because the HMM consumes the digest that moves every state boundary (1782 regions match on one edge, the other off by a median of one bin). **Order of work: digest values first, boundaries second**. Since then `pre_z` (F230) was found seeded at -10000 upstream, not 0 -- four tests were pinning the dropped leading zero-depth run -- and `pnorm2` subtracts in f32. `mono` digest record counts now match exactly (451893) and first records align on all four signals, but the weight rounding in `generate_weight_mapping` was also closed in F231 by reading the generated C: `p_i` f32, `s` f32 step-rounded, quotient a single-precision `divss`. All four digested tracks are now **byte-identical** (1169003/1169003 records), enforced as a hard gate in `oracle/check_hmmratac.py`. Jaccard is still 0.6587, and F232/F233 then closed the rest: the bin value comes from the run *ending* at the bin (the walk advances `while p1 < p2` over run ends), and `StopIteration` ends the whole walk rather than clamping to the last run. **G13 now passes: Jaccard 1.0000, 0 unmatched regions, summit shift 0 bp, all four digested tracks and both `*_accessible_regions.narrowPeak` and `*_states.bed` byte-identical. The CI gate is strict.** |
| hmmratac digested | `the_hmm_pileup_merges_adjacent_equal_runs` | merged runs match upstream's record structure; values within 1e-5 (F207) |
| L2 proptest | `crates/*/tests/invariants.rs` | 44 property tests across core/rle/pileup/score/stats |
| L3 differential | `oracle/run_e2e.sh`, `run_peak_e2e.sh`, `run_bam_e2e.sh`, `check_*.py`, plus the committed `tests/stages/**/stages.json` corpus and `macs-compare --stages` | 5 recorded runs, **1122** numeric stage leaves; `treat_pileup` byte-identical on all 5 |
| L4 golden | `oracle/run_golden.sh` | 7079/7685 cases, 22245/22456 files |
| L5 regression | `tsize_tests`, `paired_*_regressions`, `callpeak_cli`, `thread_determinism`, `malformed_inputs` | 1 test per finding |
| L6 fuzzing | `fuzz/fuzz_targets/*.rs` (libFuzzer, nightly) + `crates/macs-io/tests/malformed_inputs.rs` (CI-runnable) | 6 targets, all building under `cargo fuzz build`; 40000 execs each clean |

The property tests are not decoration. Writing them found **F192**: `d` was taking
`|end - pos|` where upstream takes the *signed* `col3 - col2` and keeps only
lengths `> 0`, so records with `col3 < col2` were counted with a positive
magnitude. On `3 x 100 bp + one reversed -800 bp` upstream reports `tag size = 100`
and this port reported `275` -- and because `d` feeds `maxgap` and the minimum peak
length, the whole run shifted. No recorded fixture contained a reversed BED record,
so the differential gate could not have caught it; the invariant did.

Three of the property tests also had to be written against what the code *actually*
does rather than what looked right, each time because upstream does the surprising
thing:

* pileup endpoints are clipped **independently**, so a read hanging off the left edge
  covers `[0, x + d - d/2)` -- deriving `end` from an already-clamped `start`
  over-counts exactly those reads;
* the non-centred minus strand is the mirror with shift `+(five - three)`, and the
  centred one is *identical* rather than mirrored (that was F190);
* `max_with` truncates at the shorter covered region, and `push*` resolves equal-end
  ties by input order, so a parallel producer must sort by `(end, value)`.

### Stage-by-stage differential (L3)

The porting plan's method step 1 requires the oracle to freeze "final outputs **and**
intermediates". `oracle/dump_stages.py` could always observe those in-process, but
nothing was kept, so there was no recorded stage corpus and L3 had nothing to diff.
`oracle/record_stages.py` now drives it over a fixture selection and commits the
result to `tests/stages/`:

| stage | leaves compared |
|---|---|
| `reads_pre_filter` | 22 numeric, 18 exact |
| `reads_post_filter` | 10 numeric |
| `d` | 5 numeric |
| `scaling` | 15 numeric, 5 exact |
| `xls_header` | 15 numeric |
| `qvalue_table` | 195 numeric |
| `peak_files` | 10 exact (file references) |

`macs-compare --stages <golden> --ours-stages <ours>` reports per-stage max absolute
and relative deviation in **pipeline order**, so the first `FAIL` is the earliest
stage that went wrong rather than an alphabetical accident. Numbers are compared
against a tolerance (the p-score criterion's 1e-9 by default); strings and booleans
must match exactly, because a `tocontrol` flag is a decision and not a measurement;
file references are excluded, since the two trees record paths under different roots
and comparing them would fail every run.

**The comparator is tested for failing, not just for passing.** A stage comparator
that always reports `ok` would gate L3 green, which is worse than having none, so
`crates/macs-compare/tests/stages.rs` perturbs a read count, a `d`, and a
1e-6 `ratio` and requires each to be caught -- and CI re-runs the same mutation
against a copy of the committed corpus and requires a non-zero exit.

#### Our side of the differential

`MACS3_RS_DUMP_STAGES=<dir>` makes `callpeak` write the same schema, and
`oracle/run_our_stages.sh` replays the recorded fixtures through it. It is an
**environment variable, not a flag**: a `--dump-stages` option would be an argument
upstream does not have, and the release definition requires identical accept/reject
behaviour, so the flag would itself be a divergence. `MACS3_RS_THREADS` sets the
precedent.

The differential now runs on both sides and, over the 5 recorded fixtures:

| stage | numeric leaves | max abs dev | verdict |
|---|---|---|---|
| `reads_post_filter` | 10 | 0 | **exact agreement** |
| `d` | 5 | 0 | **exact agreement** |
| `xls_header` | 15 | 0 | **exact agreement** |
| `peak_files` | 10 exact | 0 | **exact agreement** |
| `reads_pre_filter` | 3 | 0 | partial: treatment only |
| `scaling` | 3 | 0 | partial: single-end only |
| `qvalue_table` | -- | -- | absent on our side |
| `pvalue_histogram` | -- | -- | ours only |

**Zero numeric deviation on every stage both sides record** -- the retained read
counts, `d`, and the xls header agree with upstream exactly, across single-end,
paired-end and counted-FRAG. What is not yet covered is recorded as *missing* rather
than approximated, and the reasons are specific:

* `reads_pre_filter` is captured at `pool_se`/`pool_pe` time, before `filter_dup`.
  Capturing it later is the trap: `filter_dup` rewrites the track, so the counts would
  be post-filter while the stage name says pre-filter -- and on the paired-end fixtures
  that manufactures a 3000-vs-50 "divergence" out of nothing. Only the treatment half
  is available at that point (the control is pooled later), so the `control`
  sub-object is absent rather than reconstructed.
* `scaling` is only recorded on the single-end path, where `SeSignalSetup` carries the
  ratio and the two windows. The paired-end path computes its own ratio and does not
  populate it.
* `qvalue_table` is upstream's `--cutoff-analysis` output (`pscore qscore npeaks lpeaks
  avelpeak`) -- a p->q map *with peak counts*. We emit the AFDR **histogram** that map
  is built from, under the name `pvalue_histogram`, deliberately not `qvalue_table`: the
  two are related but structurally different, and reusing the name would have the
  comparator report a structural mismatch as if it were a numeric one. The
  `--cutoff-analysis` rows themselves *are* now implemented and byte-identical (below).

CI runs this differential and reports gaps as a warning rather than a failure until the
gaps above are closed, so the layer is live and its coverage is visible rather than
either absent or falsely green.

#### `--cutoff-analysis` (implemented, byte-identical)

`callpeak --cutoff-analysis` was parsed but **not implemented**: the flag was accepted
and no `NAME_cutoff_analysis.txt` appeared, which is an accept/accept divergence --
upstream writes a file and this wrote nothing. It is now implemented, and it is more
than a report:

* the ladder is `0.3, 0.6, ... 9.9` (33 values, descending), transcribed from
  `np.arange(0.3, 10.0, 0.3)` with `round(x, 5)`. **The rounding is to five decimals,
  not to an integer** -- rounding to integers collapses the ladder onto `{0..10}`,
  which still looks like a plausible 33-value list and silently compares the wrong
  cutoffs. That was the bug found on the first run.
* each cutoff's peak count and total length come from the same p-score track the
  pipeline already computes, folded in **during** the q-table pass rather than in a
  second pass over the data.
* the cutoffs are then **seeded into the AFDR histogram with zero length** before the
  q-table is built (`pscore_stat[cutoff] = 0`, `CallPeakUnit.py:978`). A zero-length
  bucket leaves `N` unchanged but adds a key to `unique_values`, so the q-table gains an
  entry and every q-score shifts. A run with the flag therefore legitimately produces a
  *different* q-table from the same data as one without it -- which is why this could
  not be faked as a post-hoc dump.
* the NumPy quirk `pos_array[above_cutoff - 1]` **wrapping to the last element** when
  the first above-cutoff index is 0 is reproduced. Without it the first chunk of such a
  peak starts at the wrong place and `lpeaks` is wrong.

Verified by `oracle/check_cutoff_analysis.sh` against upstream's recorded files:

```text
ok    se_basic/gauss_two_peaks [se]        (header only -- nothing clears 0.3)
ok    pe_basic/gauss_fragments [pe]
ok    pe_basic/atac_short [pe]
ok    pe_basic/nucleosome_ladder [pe]
ok    frag_basic/barcode_fragments [frag]  (34 rows)
5 identical, 0 differing
```

### Compatibility matrix

`docs/compatibility-matrix.md` is generated, never hand-edited:

```
python3 oracle/gen_compat_matrix.py --jobs 14
# 21391/21612 files identical, 0 not produced, 20 (subcommand, variant) pairs
```

Two reporting details that keep it from lying: a file the replay did not produce is
counted as *not compared* rather than as a difference (an invocation that legitimately
writes fewer files would otherwise read as 0% agreement), and the `*.xls` header's
`# Command line:` line has the replay's scratch `--outdir` rewritten back to the
recorded one before hashing -- without that, every xls in the corpus reads as
differing by one path.

The matrix currently covers **only `callpeak`**: all 7685 recorded invocations are
callpeak. That is the honest state of the corpus, and the generator's job is to say so
rather than to hide it.

### Architecture requirement: chromosome-level parallelism

The plan requires "chromosome-level Rayon parallelism that is provably
output-identical to 1 thread". Rayon was a declared dependency of three crates but
**unused**; it is now the executor of three stages:

1. `build_signals_se` -- per-chromosome treatment pileup and merged local lambda,
   split into `build_one_se_chromosome` so the loop can be parallel;
2. `build_qtable_from` -- per-chromosome p-score tracks;
3. `callpeak` -- per-chromosome peak calling, via `par_iter().map_init(..)`.

Determinism is structural, not incidental:

* every stage iterates a `Vec`/`&[T]` and collects into a `Vec`, so rayon's
  indexed `collect` returns results **in slice order** whatever the schedule;
* each chromosome gets its own `PScoreCache` via `map_init`, and the cache is pure
  memoisation, so it cannot change a score;
* the one cross-chromosome object, the AFDR histogram, is still folded
  sequentially in chromosome order (and its buckets are `i64`, so even a parallel
  reduce would be exact).

`crates/macs-cli/tests/thread_determinism.rs` asserts byte-identical
`*.xls` / `*_peaks.narrowPeak` / `*_summits.bed` at 1, 2, 8 and 32 threads on a
three-chromosome fixture.

Thread count is the `MACS3_RS_THREADS` environment override, **not** a CLI flag:
upstream's `callpeak` argparse defines no `--threads`, so adding one would break
the "identical accept/reject behaviour" criterion (a real `macs3` exits 2 on it).

### Performance

Measured with `/usr/bin/time` on 24 chromosomes x 200 kbp, 3 M treatment reads
against 1.5 M control (single-end) and 1.2 M / 0.6 M (paired-end), `--nomodel`:

| workload | macs3 3.0.5 | this port, 1 thread | this port, 16 threads |
|---|---|---|---|
| single-end wall clock | 5.54 s | 1.25 s (**4.4x**) | 1.03 s (**5.4x**) |
| single-end peak RSS | 132 MB | 213 MB | 226 MB |
| paired-end wall clock | 6.63 s | 1.07 s (**6.2x**) | 0.89 s (**7.4x**) |
| paired-end peak RSS | 111 MB | 183 MB | 228 MB |

The **>=3x wall-clock** criterion is met on both workloads at every thread count,
with no command slower. The **<=50% peak RSS** criterion is **not** met: the port
is at ~1.6-2.1x upstream's RSS, which is recorded as open work.

One measure already landed for it: the per-chromosome q-score tracks are read in
exactly one place -- the `--broad` level-2 cutoff -- so a narrow-mode run no longer
materialises them (`QScoreSink::build_with_tracks`). They are `Vec<Run<f32>>` at 16
bytes per breakpoint, i.e. tens of MB on a real genome.

## Next

1. **Peak RSS** (see "Where the memory actually goes" above): 222 MB against
   upstream's 134 MB, versus a <=0.5x criterion. The measured cause is read count, not
   genome size -- so the fix is a smaller read representation and a smaller
   `Run<f32>`, not streaming. The streaming rewrite now in place is a wash on RSS and
   costs 0.46 s of wall clock, and should be reverted or its second pass removed.
2. **F191** (above): the p-score's 5th decimal, 208 cases and the single
   largest remaining item. It needs the source the compiled `Prob.so` was
   actually built from, or an empirical model of its log-space sum.
3. The remaining FRAG run-boundary and subpeak-naming edges: 96 `*_mfrag` cases
   plus 8 others.
4. `callpeak` **model mode**: the `--nomodel`-less FRAG runs need a built
   `PeakModel` (G8). Extend the golden corpus with `-p` variants (and
   `--broad -p`, `--call-summits -p`) so F189's remaining boundary case is gated
   rather than only spot-checked.
5. `macs-io` (G3) remainder: ELAND, ELANDMULTI, ELANDEXPORT and bowtie, all of
   which share the BED 5'-position convention that is pinned by 54 fixture inputs.
6. `macs-track` (G4): the BAM/SAM `FixWidthTrack` construction path (interval
   collapsing, `tsize` estimation) and `print_to_bed`. Multi-chromosome
   downsampling is **not** reproducible upstream (F29), so exact subset parity is
   impossible there and sorted byte order is the declared choice.
7. **`callvar`** (G14) is **partly implemented**: the driver is complete, its VCF
   header is byte-identical to the oracle, but the variant-calling kernel is not
   written. Done: the flag surface, the peak BED reader (`PeakIO` semantics),
   BAM treatment/control opening, the asymmetric first-chromosome `assert`, the exact
   21-line `VCFHEADER` template, the `Program_Args` reconstruction (including the
   `--top2allele-count` spelling and the double space before `--fermi`), the
   `##contig=` lines in BAM header order, and the `#CHROM` line. Verified byte-for-byte
   against a real upstream run in `crates/macs-callvar/tests/header_parity.rs`.
   `VariantStat` is now ported too (F213) and pinned bit-for-bit by
   `crates/macs-callvar/tests/variant_stat_golden.rs` against a 410-case golden
   generated from the oracle's compiled `.so`. That left two defects in the first
   port -- `log1p(-e)` transcribed as `ln(1-e)`, and CPython's `log1p(-1)` domain
   error silently becoming `-inf` -- which no tolerance-based check would have
   caught. The remaining 1900 lines are `PosReadsInfo`, `PeakVariants` and
   `RACollection`.
   Not done: upstream's `RACollection` + `PosReadsInfo` + `VariantStat` +
   `PeakVariants` chain (~3000 lines), whose reference sequence is a read *consensus*
   (`__fill_refseq`) rather than the genome, plus the fermi-lite assembly path.
   Until then `callvar` refuses **before** creating the VCF rather than emitting
   plausible-but-wrong records, which keeps "errors precede output files" true.
   Exit-status parity already holds: 2 for usage errors, 0 for `--help`, 1 for
   runtime errors -- all matching upstream on the same invocations.
7. ~~**Accepted flags that do nothing**~~ **DONE**: `oracle/audit_accepted_flags.py`
   sweeps every `store_true` flag in the auto-derived matrix, runs the subcommand with
   and without it, and reports any flag whose output is byte-identical. It found
   **F195** (`callpeak --cutoff-analysis`) and then **F197** (`bdgpeakcall
   --cutoff-analysis`, which additionally *replaces* peak calling and writes to `-o`),
   and **F198** (`bdgpeakcall -o NAME` did not rename the track or the peaks). All
   three are now byte-identical to the oracle. The golden corpus structurally cannot
   find this class of bug -- `run_oracle.py` only records invocations someone thought to
   try -- so the sweep is the standing check, and it currently reports
   **0 flags with no observable effect** across `callpeak`, `bdgpeakcall`,
   `filterdup` and `pileup`. Extending the sweep to the remaining subcommands is open.
8. ~~**No CI**~~ **DONE**: `.github/workflows/ci.yml` wires all six layers as
   separate jobs (L1 unit, L2 proptest, L3 stage-by-stage differential, L4 golden
   byte-compare, L5 per-finding regressions, L6 nightly fuzz), plus a `criteria` job
   for the behaviours that are not a test layer (exit contract, thread determinism,
   no-Python-in-binary). It also asserts that every `F<number>` heading in
   `docs/upstream-findings.md` is referenced from code, so a finding cannot be
   documented and left unpinned. Supporting scripts it needs now exist:
   `oracle/provision_oracle.sh` (clones the pinned commit, refuses a dirty tree),
   `oracle/check_exit_contract.py`, `oracle/gen_compat_matrix.py`.
   The oracle is deliberately *not* vendored -- vendoring it would make it too easy to
   "fix" the oracle and silently invalidate every golden.
9. The remaining test gap is breadth, not the layers: L2 has 44 property tests
   across the five numeric crates and L6 has 5 libFuzzer targets plus a
   220k-iteration randomized no-panic suite, but neither yet reaches the
   `macs-io` BAM/BGZF/BAI decoder, the peak writers, or the CLI argument parser.
   The CLI surface is the notable one: "identical accept/reject behaviour" is checked
   for a handful of cases in `oracle/check_exit_contract.py`, not derived from the
   315-row flag matrix.

## Blockers

None. The oracle, the corpus, the comparator and the numerical core are all in
place; the remaining work is the pipeline itself.
- **F191 sharpened (F234):** the Poisson path is exonerated -- `macs_stats::poisson_cdf` matches the compiled oracle bit-for-bit on 1220/1220 points (new harness `crates/macs-stats/tests/poisson_vs_oracle.rs`), so the residual is in the lambda track, not `log10_poisson_cdf_q_large_lambda`. Fixed en route: the bedGraph denominator is a Python double and the division is a double rounded once to f32, not an f32 divide. Golden totals unchanged (the failing cases are non-SPMR, denominator == 1.0).
- **F235:** the weighted sweep now applies all starts at a position before all ends, matching `PileupV2.py:278` (a merged `(pos, value)` sort had reversed it, and the exact-f32 coalesce test turns a last-bit difference into moved breakpoints). The F191 residual is now pinned to **one f32 ulp** in `scaled_z = z * scale_factor` in the control's windowed max: `MACS3_RS_DUMP_CTRL_LAMBDA` shows ours `15.74799538` where upstream prints `15.74799`, and the adjacent f32 is `15.74799443`. Excluded by measurement: poisson_cdf (1220/1220), `_treat_pileup.bdg` (byte-identical), the SPMR denominator, accumulation order.
- **F236:** the control-lambda ladder is now instrumented end to end (`MACS3_RS_DUMP_CTRL_LAMBDA` reports scales, ratio, factors). Four of five inputs are confirmed exact, including that upstream's `average_template_length` is an explicit **f32 cast** (so `ratio = 0.5000003632044805`, not 0.5) and that `lambda_bg = 14.28046` matches the golden's first row. The open one is `tsize_exact = 143.555` (per-line mean) vs upstream's weighted 143.579 -- both print as `143.6`, which is why it survived. The weighted variant was tried and left the golden suite unmoved, so it was **reverted, not shipped**: it needs `tp.d`'s definition, not a neutral guess.
- **F237:** the last unconfirmed ladder input is closed. Upstream carries **two different means**, both f32: `options.tsize` = `tp.d` = `f32(unweighted length sum)/line count` = 143.555 (ours), and `average_template_length` = `f32(weighted length)/count sum` = 143.579. Both print as `143.6`, so no upstream log line distinguishes them. The weighted-mean experiment reverted in F236 was therefore correctly discarded. With every input and the (deliberately misaligned) `over_two_pv_array` merge confirmed, the residual 3 ulps must be upstream of the multiply, in the accumulated depth -- none of the three factors divides our value to an integer depth, which a counted control's integer depth requires.
- **F238 (F191 progress):** `tp.d` is an **f32** (`cython.cast(cython.float, m) / i` -> `divss`), so `options.tsize` is a float32 and the wide-window control factors inherit that width; we were keeping the mean in f64. Found by solving for `(integer depth, factor)` against the observed bits rather than guessing: depth is exactly 2194, so only the factor could differ. **Golden 7079 -> 7091 cases, 22245 -> 22257 files; `_control_lambda.bdg` mismatches 14 -> 2**, and those two are coordinate (run-boundary) differences, not arithmetic. Pinned by `crates/macs-peaks/tests/f238_control_factor_f32.rs`.
- **F239 (census):** after F238 the golden corpus is 7091/7685 cases and 22257/22456 files, with **147 mismatching files across ~50 fixture/variant combinations** -- all still the 5th significant digit of `-log10(pvalue)`, no coordinate/ordering/count/exit-status mismatch anywhere. The clusters split by code path into with-control SE (`realistic`, `spikes_only`, `disjoint_chromosomes`), with-control PE (`gmini_mfrag_d400_w600_ctrl_384`, where 2 survivors are a *coordinate* run-boundary issue rather than arithmetic), and no-control PE (`gonechrom_mpe_d4000_w180_noc_131`, `gtiny_mfrag_d400_w600_noc_142`). Solving the oracle's `poisson_cdf` for the SE case shows ~2-3 f32 ulps of lambda -- the F238 signature on a path F238 did not touch.
- **F240 (G0):** intermediates `duplicates` and `lambda_ladder` are now recorded and wired into `macs-compare`. Duplicates were documented as unreachable because `filter_dup`'s rate is recomputed inside compiled code -- true of the structure, irrelevant to the value: upstream *reports* it ("Redundant rate of treatment: 0.99") and the log is now captured in-process with a root `FileHandler`. 26 numeric leaves across the 5 stage fixtures. The control scale ladder is recorded as an **empty** stage on purpose: `callpeak_cmd.py:265` is `peakdetect.call_peaks()` with no arguments, so the delegating proxy F240 added finds nothing and proves there is no Python-visible route -- a committed negative result beats an undocumented gap, and it guards against a future upstream passing it. For `-f FRAG` only totals appear, correctly: upstream runs no duplicate filtering there.
- **F241:** the remaining clusters do **not** share a cause. The no-control PE cluster (`gonechrom_mpe_*`, `gtiny_mfrag_*_noc_*`, 23 variants) shows a lambda differing by ~**30 f32 ulps** (1.9e-6 relative), not one -- solving the oracle's `poisson_cdf` for the lambda behind each printed p-score gives golden 63.83760 vs ours 63.83748. That is a systematic difference, so hunting a rounding site there would be wasted effort. The with-control SE cluster (~2-3 ulps) is the same class as F238; `gmini_mfrag_d400_w600_ctrl_384` is a coordinate/run-boundary issue, not arithmetic.
- **F242:** our PE-without-control path builds the with-control ladder unconditionally, and with no control `ratio = treat_sum/0 = 0` so **every factor is 0.0** -- an all-zero control track. Upstream's `__call_peaks_wo_control` builds a single-scale control instead (`[treat_length/(lregion*treat_total*2)]`, `ctrl_d_s=[lregion]` in PE mode). Implementing that made the output *worse* (peak boundaries and p-score moved away from the oracle) partly because my factor indexing was wrong and partly because a change that moves away from the oracle proves nothing until understood -- so it was **reverted**. Kept as a documented asymmetry; golden back at 7091/7685.
- **F243:** the SE ladder is now instrumented (`SE lambda scales: d=200 slocal=1000 llocal=10000 factors=[0.09643395, 0.01928679, 0.0019286791] lambda_bg=0.0192`), and the F238 technique applied to `se_model/realistic/call_summits` (depth 489, p-score 5.97002 vs 5.97001, lambda ~392.2) gives a decisive negative: **no `(integer depth, factor)` pair reproduces either printed value** over depths 1..400000 x all three factors. **RETRACTED by F244**: the factors are exactly right (`f32(192/1991)` and its lregion siblings), and the differing column is `fold_enrichment`, not `-log10(p)` -- the p- and q-scores on that fixture were never wrong.
- **F245 (corrected census):** `oracle/run_golden.py` now names the differing **column** instead of the line (`_named_diff`), plus a standalone `oracle/check_peak_columns.py`. Re-censusing all 147 mismatching files: **86 differ only in `fold_enrichment`/`signalValue` -- the p- and q-scores are byte-identical there**, so the `p-score <= 1e-9` / `q-score <= 1e-6` numeric criteria are met for those; 34 differ in `start` (real peak-boundary shifts, more serious than F239 implied); 2 in peak `name`; 25 unnamed `.bdg`/`.bed`. F239's "all 594 cases are the 5th digit of -log10(p)" was wrong for the same reason F243 was -- positions read as columns.
- **F246:** `fold_enrichment` is **verified**, not inferred: `CallPeakUnit.py:1385` is `(summit_treat + self.pseudocount) / (summit_ctrl + self.pseudocount)` with `pseudocount = 1` -- a pseudocount on **both** sides, which `9/lambda` would get wrong. Our `regions.rs` already matches. **F247 retracts the magnitude**: `%.6g` makes that an interval, not a point. From our own value, `c+1 = 1.675037682` (f32-nested) vs `1.675037669` (f64) vs upstream `<= 1.675037667` -- the sides agree to **~1e-8** and both candidates sit within 1.5e-8 of the `%.6g` boundary. Not a systematic error; a print-boundary artifact. Ruled out on this cluster: the ladder factors (bit-exact), the Poisson path (1220/1220, and the printed p-score is identical anyway), the pseudocount formula, and the peak boundary (`start`/`end`/`length`/`abs_summit`/`pileup` all match). Left: which windowed-maximum scale wins, or the width of one intermediate inside the window computation.
- **F248:** the 86 `fold_enrichment` files reduce to **one f32 ulp in the SE control factor**. The `%.6g` boundary requires `ctrl <= 0.675037667409546`; ours is `7 * f32(192/1991) = 0.675037682056427` (above), and `7 * (1 ulp below f32(192/1991)) = 0.6750376224517822` (below, = upstream). Crucially this is **not** F238 again: both sides already hold f32, and the inputs are exact integers, so `f32(38400/398200)` is a single well-defined value -- the one we have. Therefore upstream's `ratio` reaching `ctrl_scale_s` is not `38400/398200` at all. F248 lists three checkable candidates, each naming the recorded stage that holds the evidence.
- **F249 (corrects F248's premise):** `38400/398200` and `192/1991` are the **same f64**, so our SE factor is upstream's -- F248's "the ratio is not 38400/398200" was wrong. Further, all three ladder scales (`z`=7/35/350) produce the **same** f32 `0.675037682056427`, so which scale wins is irrelevant on this fixture and F248's candidate 3 is eliminated too. What remains: our product is exactly one f32 ulp above the `%.6g` boundary, and with both operands exactly specified there is only one correctly-rounded answer -- so `peak_content[summit_index][3]` is not `z * scale_factor`. It must be captured behaviourally (as `over_two_pv_array` was in F168), not read from source.

## Declared deviation: the `--threads` acceptance criterion

**Decision (2026-10-02, user-confirmed): enforce thread-count invariance on the internal knob;
do not add a `--threads` flag.**

The porting plan requires "Identical output at `--threads 1` and `--threads 32`". The pinned oracle
does not have that flag:

```console
$ macs3 callpeak -t x.bed --threads 4
macs3: error: unrecognized arguments: --threads 4
```

and the auto-derived `oracle/flag_matrix.tsv` contains no `--threads` row, so this is a property of
MACS3 3.0.5 rather than of our port.

Accepting the flag would mean accepting an invocation the oracle rejects, which directly contradicts
the other acceptance criterion, "Identical accept/reject behaviour on invalid invocations". Drop-in
CLI compatibility wins, so the criterion is enforced as what it can only mean: **the number of
worker threads must not be observable in any output byte.**

Enforcement:

- Knob: `MACS3_RS_THREADS` (`configure_pool`, `crates/macs-cli/src/commands/callpeak.rs`).
- Gate: `oracle/check_thread_invariance.py` -- runs every recorded case twice, pinned to one worker
  and to 32, and compares every output file byte-for-byte. `.xls` files echo their own command line
  and outdir, so they are compared with `run_golden.py`'s existing normalisation; without it every
  case would differ on the wrapper and the gate would be vacuous.
- CI: wired into the `L4 golden byte-compare` job as "output is identical at 1 and 32 threads",
  stride-sampled to 600 cases to bound job time (`--limit 0` replays the whole corpus).
- Status: **500 cases compared, 0 differing.** Rayon is used only for per-chromosome work, so the
  failures this can catch are chromosome-ordering races, shared-cache mutation and non-deterministic
  accumulation order.

Two traps this harness had to be written around, both of which initially produced a *clean* report:

1. `command.json` sits at `tests/golden/<group>/<fixture>/<variant>/`, and the grouping depth is not
   fixed. A two-level walk found **zero** cases and printed "cases differing: 0". The walk now
   searches for the marker and the script exits non-zero on an empty corpus.
2. Thread-count invariance is only meaningful on multi-chromosome inputs, so the corpus must actually
   contain them for this gate to have teeth.

## Golden corpus: at parity (2026-10-02)

```text
cases with output    : 7204 (byte-identical: 7204/7204)
output files matched : 22456/22456
invalid-invocation   : 481 cases, exit-status mismatch 0
```

Every recorded output file in the pinned corpus is byte-identical, and every exit status matches.
That covers `*_peaks.xls`, `*_peaks.narrowPeak`, `*_summits.bed`, `*_control_lambda.bdg`,
`*_treat_pileup.bdg`, `*_peaks.broadPeak`, `*_peaks.gappedPeak` and the `*_cutoff_analysis.txt`
ladders across all 14 commands' recorded invocations, in narrow, broad and `--call-summits` modes.

The 481 "invalid-invocation" cases are the rejection contract: 425 argparse rejections (rc 2) and 56
runtime rejections (rc 1), which produce **no output files** by design. They are counted separately
because a byte-compare cannot succeed on an empty file set -- see F273, which fixes a summary line
that had been reporting this as 93.7%.

Other gates re-verified at the same time:

| gate | result |
|---|---|
| workspace tests | 714 passed, 0 failed |
| `cargo fmt --all --check` | clean |
| `cargo clippy --workspace --all-targets` | 0 warnings |
| thread invariance (`MACS3_RS_THREADS` 1 vs 32) | 500 cases, 0 differing |
| hmmratac (G13) | PASS, Jaccard 1.0000 |
| callvar | 22/22 (`-F off`), 16/16 (`-F auto`) byte-identical |
| oracle tree pristine | 6/6 checks, 0 modified files |

### What this does not close

Golden byte-parity is one criterion. Still open against the porting plan:

- **Performance**: `>=3x` wall-clock and `<=50%` peak RSS versus upstream on the headline benchmark
  matrix, with zero commands slower and published methodology.
- **Pure-Rust `callvar` assembler**: the sanctioned last item, currently bridged through the vendored
  fermi-lite C FFI.
- **Test-layer breadth**: the CI jobs exist (L1 unit, L2 proptest, L3 stage differential, L4 golden,
  L5 per-finding regressions, L6 nightly fuzz, plus a criteria job) but the L3 stage corpus covers
  5 fixtures against 7685 golden cases; widening it is the natural follow-on now that golden parity
  makes stage capture trustworthy.
- **Stage-capture coverage**: the plan's intermediate list is not fully recorded. `--nolambda`'s
  `pos_array` origin is the concrete known gap (F271/F272), and it is the reason the last case needed
  a flag-conditional rather than a derived constant.
