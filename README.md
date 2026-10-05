# macs3-rs

A Rust reimplementation of [MACS3](https://github.com/macs3-project/MACS), the
standard peak caller for ChIP-seq and ATAC-seq. Same command line, **byte-identical
output** on the recorded `callpeak` corpus, 2.3–4.9× faster.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/parity-dark.svg">
  <img src="docs/figures/parity-light.svg" alt="7,204 of 7,204 recorded runs byte-identical; 22,456 of 22,456 output files byte-identical; 0 exit-status mismatches; 2.3 to 4.9 times faster across 12 workloads">
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
| `callpeak` on 5 M + 5 M real CTCF reads, five modes | 17 / 17 files byte-identical |
| `callpeak --call-summits` on 5 M + 5 M reads, model and `--nomodel` | 8 / 8 files byte-identical |
| `pileup`, `bdgcmp`, `cmbreps`, `bdgopt`, `bdgpeakcall` on the same data | 5 / 5 files byte-identical |
| MACS3's own `cmdlinetest` suite, 14 subcommands | passes; 154 / 163 output files byte-identical |
| `hmmratac` on yeast ATAC-seq | 1,605 / 1,605 accessible regions identical |
| `callvar`, assembly off / auto / on | 22 / 22, 16 / 16, 15 / 15 VCF records identical |
| 1 thread vs 32 threads, all 7,685 recorded runs | byte-identical output |

Fresh real-data comparisons normalize the XLS `# Command line:` header, which
contains the different executable and output paths. All other output bytes are
compared unchanged.

Of the other 9, seven are `hmmratac` model files that are numerically
equivalent (parameters agree to 2 × 10⁻¹⁰ and the peak calls are identical),
and two depend on MACS3's unordered chromosome iteration, which differs between
MACS3's own runs.

Getting there meant reproducing upstream's arithmetic exactly: float32 widths,
NumPy's summation order, its Mersenne Twister stream, and a number of upstream
bugs. The [findings log](docs/upstream-findings.md) records each one.

## Performance

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/performance-dark.svg">
  <img src="docs/figures/performance-light.svg" alt="Speedup over MACS3 and peak memory relative to MACS3 for twelve workloads">
</picture>

Real CTCF ChIP-seq, 5 M treatment + 5 M control reads, median of 3 runs on a
Ryzen 7 3700X. Reproduce with [`scripts/bench_real.py`](scripts/bench_real.py).

| workload | MACS3 3.0.5 | macs3-rs | speedup | peak memory |
|---|---|---|---|---|
| `callpeak`, single-end `-B` | 42.9 s / 294 MB | 12.6 s / 165 MB | 3.4× | 56% |
| `callpeak`, single-end `--SPMR` | 32.4 s / 303 MB | 8.4 s / 139 MB | 3.9× | **46%** |
| `callpeak`, single-end `--broad` | 38.0 s / 291 MB | 10.5 s / 152 MB | 3.6× | 52% |
| `callpeak`, no control | 17.3 s / 165 MB | 4.6 s / 76 MB | 3.8× | **46%** |
| `callpeak`, paired-end BAM | 1.22 s / 80 MB | 0.34 s / 26 MB | 3.6× | **32%** |
| `pileup` | 7.5 s / 134 MB | 3.2 s / 32 MB | 2.3× | **24%** |
| `filterdup` | 5.4 s / 134 MB | 1.56 s / 22 MB | 3.5× | **16%** |
| `randsample` | 4.5 s / 134 MB | 1.71 s / 94 MB | 2.6× | 70% |
| `bdgpeakcall` | 20.1 s / 225 MB | 4.5 s / 200 MB | 4.5× | 89% |
| `bdgopt -m p2q` | 29.1 s / 231 MB | 5.9 s / 234 MB | 4.9× | 101% |
| `cmbreps -m max` | 43.8 s / 391 MB | 10.0 s / 194 MB | 4.4× | **50%** |
| `bdgcmp -m ppois` | 56.6 s / 611 MB | 18.7 s / 195 MB | 3.0× | **32%** |

Every workload is faster than MACS3, by 2.3× to 4.9×. Bold marks memory at or
under the project's target of half of MACS3's; seven of twelve meet it.
`callpeak -B` and `--broad` are just over, and `randsample`, `bdgpeakcall` and
`bdgopt` hold their whole parsed input.

Streaming output instead of buffering it accounts for the largest gains:
`bdgcmp` fell from 1,383 MB to 195 MB and `pileup` from 314 MB to 32 MB.
[Where the memory goes](docs/compatibility.md#where-the-memory-goes) has the breakdown.

## Quick start

```sh
cargo install --path crates/macs-cli
macs3-rs callpeak -t chip.bed.gz -c input.bed.gz -f BED -g hs -n sample
```

Subcommands and flags are spelled as in MACS3. Thread count is set with
`MACS3_RS_THREADS` rather than a flag, because MACS3 rejects `--threads` and so
must a drop-in replacement.

## Status

The recorded `callpeak` corpus matches byte for byte. Fresh 5 M-read
`--call-summits` runs also match in model and `--nomodel` modes, including the
model script and cutoff analysis. The other 13 commands are checked by
smaller comparisons; running them on the 5 M-read data found three output
differences in `bdgpeakcall`, `bdgopt` and `bdgcmp`, all now fixed and covered
by tests.

- **Known difference:** on multi-chromosome input `filterdup` writes the same
  rows in a different chromosome order, and seeded `randsample` picks different
  reads. MACS3's own output is not repeatable in either case.
- **Open:** five workloads exceed the 50% memory target.
- **Open:** `callvar` still calls the original fermi-lite assembler through a
  small C bridge.

Details: [compatibility by command](docs/compatibility.md) ·
[current status](docs/status.md) · [porting plan](PORTING_PLAN.md)

## How it is built

14 crates, from `macs-core` and `macs-stats` up to `macs-cli`, under `crates/`.
`oracle/` holds the differential harness that replays recorded MACS3 runs, and
`tests/` holds the fixtures and recorded outputs. MACS3 itself is pinned, not
vendored, so the reference cannot be edited to make a test pass.

```sh
cargo test --workspace     # 897 tests
oracle/run_golden.sh       # replay every recorded run and compare bytes
```

## Licence

BSD-3-Clause, matching MACS3. The vendored fermi-lite sources keep their MIT
licence.
