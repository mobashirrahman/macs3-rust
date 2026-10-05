# macs3-rs

A Rust reimplementation of [MACS3](https://github.com/macs3-project/MACS), the
standard peak caller for ChIP-seq and ATAC-seq. Same command line, **byte-identical
output**, 3.4–4.0× faster.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/parity-dark.svg">
  <img src="docs/figures/parity-light.svg" alt="7,204 of 7,204 recorded runs byte-identical; 22,456 of 22,456 output files byte-identical; 0 exit-status mismatches; 3.4 to 4.0 times faster on real ChIP-seq">
</picture>

> **Development project.** All 14 subcommands are implemented and the recorded
> corpus is at byte parity, but it has not been validated on every flag of every
> command. See [what is still open](#status).

## Accuracy

MACS3's own test suite accepts a Jaccard index above 0.99 on peak files. This
project compares every output byte instead, because a score that moves in the
fourth decimal is a different answer.

| check against pinned MACS3 3.0.5 | result |
|---|---|
| `callpeak` flag matrix, 7,204 runs across 425 fixtures and 19 option variants | **22,456 / 22,456 files byte-identical** |
| Invalid invocations, 481 cases | **0 exit-status mismatches** |
| Real CTCF data: single-end, BEDPE, BAMPE | peaks, summits and signal tracks byte-identical |
| `hmmratac` on yeast ATAC-seq | accessible regions identical (Jaccard 1.0) |
| `callvar`, with and without assembly | VCF records identical |
| 1 thread vs 32 threads | byte-identical output |

Getting there meant reproducing upstream's arithmetic exactly: float32 widths,
NumPy's summation order, its Mersenne Twister stream, and a number of upstream
bugs. The [findings log](docs/upstream-findings.md) records each one.

## Performance

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/performance-dark.svg">
  <img src="docs/figures/performance-light.svg" alt="Speedup over MACS3 and peak memory relative to MACS3 for five callpeak workloads">
</picture>

Real CTCF ChIP-seq, 5.0 M treatment + 5.0 M control reads, human `hs` genome
size, release build, single run:

| `callpeak` workload | MACS3 3.0.5 | macs3-rs | speedup | peak memory |
|---|---|---|---|---|
| Single-end, bedGraph out (`-B`) | 44.3 s / 301 MB | 12.9 s / 169 MB | **3.4×** | 56% |
| Single-end, `--SPMR` | 32.8 s / 303 MB | 8.7 s / 142 MB | **3.8×** | **47%** |
| Single-end, broad | 51.7 s / 309 MB | 14.9 s / 168 MB | **3.5×** | 54% |
| Single-end, no control | 23.1 s / 173 MB | 7.2 s / 87 MB | **3.2×** | **50%** |
| Paired-end BAM | 1.41 s / 80 MB | 0.41 s / 31 MB | **3.4×** | **39%** |

The same host, on the rest of the release benchmark matrix:

| workload | MACS3 3.0.5 | macs3-rs | speedup | peak memory |
|---|---|---|---|---|
| `pileup`, single-end | 7.8 s / 138 MB | 3.4 s / 33 MB | 2.3× | **24%** |
| `filterdup` | 5.7 s / 138 MB | 1.6 s / 23 MB | 3.6× | **17%** |
| `randsample` | 4.2 s / 138 MB | 1.5 s / 52 MB | 2.8× | **37%** |
| `bdgpeakcall` | 7.5 s / 107 MB | 1.7 s / 81 MB | 4.4× | 76% |
| `bdgopt -m p2q` | 10.9 s / 106 MB | 2.4 s / 117 MB | 4.6× | 110% |
| `cmbreps -m max` | 44.6 s / 401 MB | 10.1 s / 199 MB | 4.4× | **50%** |
| `bdgcmp -m ppois` | 60.4 s / 626 MB | 19.6 s / 199 MB | 3.1× | **32%** |

Two honest notes. The speed target (≥3×) is **missed** on `pileup` (2.3×) and
`randsample` (2.8×), which are I/O-bound, not compute-bound. The memory target
(≤50%) is met on seven of the eleven workloads; `callpeak -B`/`--broad` sit at
54–56%, and the two still above the line are input-bound: `bdgpeakcall` (76%)
and `bdgopt` (110%) each parse the same bedGraph upstream parses, and that
parse *is* the peak (72 MB of the 107 MB upstream total), so halving it means
streaming the input, not shrinking the output.

Where the memory goes, measured:

* `callpeak` SE: ~43 MB is the two resident position arrays, then one
  chromosome's signal build plus the q-table histogram — the parallel window's
  working set. `CHUNK` chromosomes are built at once and the window is the
  memory/speed dial.
* the bedGraph and spool bodies are streamed to per-chromosome temp files and
  concatenated, so `-B` no longer buffers a chromosome's ~38 MB of text.
* every writer streams. Holding whole files as strings is what made `pileup`
  (313 MB), `filterdup` (158 MB), `cmbreps` (1096 MB) and `bdgcmp` (1383 MB)
  *worse* than upstream; all are fixed.
* `bdgcmp` and `cmbreps` no longer build their whole-genome result. The
  bedGraph merge emits one row per breakpoint of *either* input — 25 M rows on
  the benchmark, more than both inputs — and every scorer but `qpois` is a pure
  function of `(treat, ctrl)`, so the row is scored and written during the walk
  instead. `cmbreps`' combined track is likewise written as it is produced.
  Together those took `bdgcmp` from 598 MB to 199 MB and `cmbreps` from 376 MB
  to 199 MB.
* `MALLOC_ARENA_MAX=1` and `MALLOC_TRIM_THRESHOLD_=0` are applied by re-exec
  (glibc reads them before `main`). The arena cap stops per-thread retention;
  trimming on free returns the chunked passes' dropped chromosomes to the OS
  instead of letting RSS creep up across the pass. A user-set value wins.

## Quick start

```sh
cargo install --path crates/macs-cli
macs3-rs callpeak -t chip.bed.gz -c input.bed.gz -f BED -g hs -n sample
```

Subcommands and flags are spelled as in MACS3. Thread count is set with
`MACS3_RS_THREADS` rather than a flag, because MACS3 rejects `--threads` and so
must a drop-in replacement.

## Status

- **Open:** single-end peak memory, as above.
- **Open:** the byte-parity corpus covers `callpeak`; the other 13 commands are
  checked by smaller command-specific comparisons.
- **Declared deviations:** seeded `randsample` on multiple contigs (MACS3 itself
  is not repeatable there), and `callvar` still calls the original fermi-lite
  assembler through a small C bridge.

Details: [compatibility by command](docs/compatibility.md) ·
[current status](docs/status.md) · [porting plan](PORTING_PLAN.md)

## How it is built

14 crates, from `macs-core` and `macs-stats` up to `macs-cli`, under `crates/`.
`oracle/` holds the differential harness that replays recorded MACS3 runs, and
`tests/` holds the fixtures and recorded outputs. MACS3 itself is pinned, not
vendored, so the reference cannot be edited to make a test pass.

```sh
cargo test --workspace     # 753 tests
oracle/run_golden.sh       # replay every recorded run and compare bytes
```

## Licence

BSD-3-Clause, matching MACS3. The vendored fermi-lite sources keep their MIT
licence.
