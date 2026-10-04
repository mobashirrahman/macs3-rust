# macs3-rs

A pure-Rust reimplementation of [MACS3](https://github.com/macs3-project/MACS3)
(Model-based Analysis for ChIP-Seq), targeting **byte-for-byte output parity**
with MACS3 3.0.5.

> **Status: not yet a drop-in replacement.** See [Compatibility](#compatibility)
> for exactly which invocations are byte-identical and which are not. Read it
> before substituting this for `macs3`.

## What it is

MACS3 is a Cython/Python toolkit; this is a from-scratch Rust implementation of
the same numerics, validated differentially against the pinned upstream at
commit `c5443190e3edfeb301cc94acf450e2b2c026a223` (v3.0.5).

The goal is stricter than upstream's own test suite: upstream's `test/cmdlinetest`
accepts a Jaccard index > 0.99 on peak files, whereas this project compares
output files byte-for-byte, because a score that moves in the fourth decimal is
a different answer.

## Layout

| crate | contents |
|---|---|
| `macs-core` | coordinates, intervals, chromosome interning, errors |
| `macs-stats` | Poisson / binomial, log-space paths, NumPy-compatible RNG |
| `macs-rle` | run-length-encoded genomic signal |
| `macs-io` | BED / BEDPE / FRAG / bedGraph / BAM parsers and writers |
| `macs-track` | single-end, paired-end and fragment tracks; dedup; subsampling |
| `macs-pileup` | directional and bidirectional pileup |
| `macs-score` | p-score and p-to-q-score tables |
| `macs-model` | fragment-size model (`PeakModel.py`) |
| `macs-peaks` | peak calling, broad calling, summits, cutoff analysis |
| `macs-bedgraph` | `BedGraph.py` port |
| `macs-hmmratac` | HMM peak calling for ATAC-seq / scATAC-seq |
| `macs-callvar` | variant calling in peaks, incl. a fermi-lite FFI bridge |
| `macs-cli` | the `macs3-rs` binary: all 14 subcommands |
| `macs-compare` | differential comparator used by the test gates |

## Install

```sh
cargo install --path crates/macs-cli
```

## Usage

The binary mirrors upstream's CLI, including subcommand names and flag spellings:

```sh
macs3-rs callpeak -t chip.bed.gz -c input.bed.gz -f BED -g hs -n sample
```

Thread count is controlled by the `MACS3_RS_THREADS` environment variable, not a
flag: upstream's `callpeak` defines no `--threads`, and accepting an invocation
upstream rejects would break drop-in compatibility. Output is byte-identical at
1 and 32 threads.

## Compatibility

Verified differentially against upstream 3.0.5 on two corpora: a generated
7685-invocation `callpeak` flag matrix, and upstream's own `test/cmdlinetest`
run on real CTCF ChIP-seq data (5.0 M + 5.0 M reads, and the yeast ATAC BAM).

| corpus | result |
|---|---|
| generated flag matrix, 7685 invocations | **7204 / 7204 cases byte-identical**, 22456 / 22456 files, 0 exit-status mismatches |
| upstream `test/cmdlinetest`, 163 output files | **103 byte-identical**, 17 differing, 43 not produced |
| CLI accept/reject, 255 malformed-value probes | **1 exit-status mismatch** (an upstream crash on `--tempdir`) |

**Byte-identical:**

- `callpeak` single-end, `--nomodel` (`--extsize`, `--shift`, `--nolambda`, `--SPMR`)
- `callpeak` single-end, **model mode** — the default invocation; `d`, `*_model.r`,
  `# alternative fragment length(s)` and all five outputs match on real data
- `callpeak` paired-end with `-f BEDPE` and `-f BAMPE`, narrow and broad
- `callpeak -f FRAG`, with and without `--barcodes`
- `predictd`, single-end (BED/BAM/SAM) and both paired-end modes
- `refinepeak` (single-end BED/BAM/SAM; paired-end is SE-only upstream)
- `pileup`, single-end (BED/BAM/SAM) and paired-end (BEDPE/BAMPE/FRAG)
- `filterdup`, single-end (BED/BAM/SAM) and paired-end (BEDPE/BAMPE)
- `randsample`, single-end (BED/BAM/SAM) and paired-end (BEDPE/BAMPE single-chromosome byte-identical)
- `callvar`, both `-F off` and `-F auto`
- the `bdgopt` / `bdgcmp` / `cmbreps` / `bdgdiff` bedGraph family
- `bdgpeakcall`, `bdgbroadcall`, gzipped input throughout
- every `<subcommand> --help` (captured from the pinned oracle)

**Known differences (17 files on upstream's own test):**

| area | state |
|---|---|
| summit tie-break | 1 peak in 724–735 differs on SE `--call-summits` and both PE modes (8 bp or a lower-scoring summit) |
| `*_model.r` | 1 line in 17: ~1e-12 relative tail in the `ycorr` vector |
| contig ordering | `filterdup` on 50k contigs: content identical, order differs (upstream hash-randomized, non-deterministic run to run) |
| `callvar` VCF | `##Program_Args` echoes the replay's outdir, otherwise identical |
| `hmmratac` | 43 output files not produced: BAM input now loads and `--model` inference is byte-identical, but hmmlearn self-training is unimplemented so the *default* invocation refuses |
| `hmmratac` | 43 output files not produced: hmmlearn Baum-Welch self-training unimplemented (see below) |


## Input formats

Single-end `BED`, `BAM`, `SAM`, `ELAND`, `BOWTIE` and paired-end `BEDPE`, `BAMPE`,
`FRAG` are supported where upstream supports them. Two notes:

- Upstream's SAM parser crashes on any minus-strand read (`TypeError` in CIGAR
  parsing), making `-f SAM` effectively unusable upstream. This port parses SAM
  correctly and is byte-identical on inputs upstream can handle.
- SE BAM 5' ends use the exclusive rightmost directly, matching
  `bam_fw_binary_parse`; BAMPE fragments use `abs(TLEN)` with leftmost-only
  proper pairs, matching `bampe_pe_binary_parse`. No `.bai` index is required.

## Performance

Measured on real CTCF ChIP-seq data (5.0 M treatment + 5.0 M control reads,
human `hs` genome size), median of 3, release build:

| workload | macs3 3.0.5 | macs3-rs | speedup | peak RSS |
|---|---|---|---|---|
| SE `--nomodel --extsize 200 -B` | 45.3 s / 297 MB | 13.3 s / 2967 MB | **3.4x** | **10.0x** |
| SE `--nomodel --SPMR` | 34.0 s / 294 MB | 4.8 s / 2304 MB | **7.0x** | **7.9x** |
| SE `--nomodel`, no control | 23.9 s / 171 MB | 7.2 s / 1391 MB | **3.3x** | **8.1x** |

Wall clock is 3-7x better. **Peak memory is 8-10x worse than upstream**, which
is a serious regression on real workloads and the main outstanding engineering
problem. Threading scales poorly here (~7% from 1 to 16 threads), so a
sequential stage dominates.

## Testing

```sh
cargo test --workspace          # 716 unit / differential / regression tests
cargo clippy --workspace --all-targets
cargo fmt --check
```

The differential corpora under `tests/` and the oracle harness under `oracle/`
replay recorded upstream invocations and compare every output byte. The oracle
itself is deliberately **not** vendored — vendoring it would make it too easy to
"fix" the oracle and silently invalidate every recorded result.

## Licence

BSD-3-Clause, matching upstream MACS3. `vendor/fermi-lite` is MACS3's own
submodule, vendored verbatim under its own MIT licence (`vendor/fermi-lite/LICENSE.txt`).
