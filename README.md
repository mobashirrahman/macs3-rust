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

Real CTCF ChIP-seq, 5 M treatment + 5 M control reads, median of 3 runs:

| `callpeak` workload | MACS3 3.0.5 | macs3-rs | speedup | peak memory |
|---|---|---|---|---|
| Single-end, bedGraph out (`-B`) | 42.7 s / 304 MB | 11.9 s / 186 MB | **3.6×** | 61% |
| Single-end, `--SPMR` | 32.3 s / 299 MB | 8.1 s / 173 MB | **4.0×** | 58% |
| Single-end, model mode | 44.4 s / 296 MB | 11.9 s / 183 MB | **3.7×** | 62% |
| Single-end, no control | 22.9 s / 169 MB | 6.7 s / 88 MB | **3.4×** | 52% |
| Paired-end BAM | 1.41 s / 80 MB | 0.41 s / 31 MB | **3.4×** | 39% |

The speed target (≥3×) is met on every workload. The memory target (≤50% of
upstream) is met for paired-end only; single-end is still 2–12 points above it.

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
