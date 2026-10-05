# macs3-rs

A Rust implementation of [MACS3](https://github.com/macs3-project/MACS)
(Model-based Analysis for ChIP-Seq), with a narrow C FFI bridge to fermi-lite,
targeting **byte-for-byte output parity** with MACS3 3.0.5.

> **Status: not yet a drop-in replacement.** See [Compatibility](#compatibility)
> for exactly which invocations are byte-identical and which are not. Read it
> before substituting this for `macs3`.

## What it is

MACS3 is a Cython/Python toolkit; this is a from-scratch Rust implementation of
the same numerics, validated differentially against the pinned upstream at
commit `c5443190e3edfeb301cc94acf450e2b2c026a223` (v3.0.5). Oracle runs pin
`OPENBLAS_CORETYPE=Haswell` alongside the numerical dependency versions in
`oracle/ENV.lock`; Rust reproduces that reduction order without a runtime BLAS
dependency.

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
| `macs-cli` | the `macs3-rs` binary: dispatch for all 14 upstream subcommands (coverage varies by command) |
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

Verified differentially against upstream 3.0.5 on a generated 7685-invocation
`callpeak` matrix, upstream's `test/cmdlinetest`, and command-specific oracle
checks. This is not yet full parity across every flag of all 14 commands; the
current measured results and open differences are tracked in [status](docs/status.md).

| corpus (recorded baseline) | result |
|---|---|
| generated flag matrix, 7685 invocations | **7204 / 7204 cases byte-identical**, 22456 / 22456 files, 0 exit-status mismatches |
| fresh upstream `test/cmdlinetest`, 20 command groups | **All 163 artifacts produced**; current scATAC, Gaussian, and Poisson checks pass (see status) |
| CLI accept/reject, 255 malformed-value probes | Historical: **1 mismatch** for `--tempdir`; custom-directory handling is now implemented and checked against upstream |

**Byte-identical:**

- `callpeak` single-end, `--nomodel` (`--extsize`, `--shift`, `--nolambda`, `--SPMR`)
- `callpeak` single-end, **model mode** — on checked real data, `d`, alternative
  fragment length(s), peak outputs, and `*_model.r` match
- `callpeak` paired-end with `-f BEDPE` and `-f BAMPE`, narrow and broad
- `callpeak -f FRAG`, with and without `--barcodes`
- `predictd`, single-end (BED/BAM/SAM) and both paired-end modes
- `refinepeak` (single-end BED/BAM/SAM; paired-end is SE-only upstream)
- `pileup`, single-end (BED/BAM/SAM) and paired-end (BEDPE/BAMPE/FRAG)
- `filterdup`, single-end (BED/BAM/SAM) and paired-end (BEDPE/BAMPE)
- `randsample`, single-chromosome SE and PE; multi-contig sampling uses sorted
  chromosomes, while upstream uses hash-randomized set iteration
- `callvar`, `-F off`, `-F auto`, and `-F on`
- fresh real-data `callpeak --call-summits` runs for SE, BEDPE, and BAMPE; the
  summit-related files are byte-identical after normalizing the XLS command-line
  header (see `oracle/check_real_summit_bytes.py`)
- the `bdgopt` / `bdgcmp` / `cmbreps` / `bdgdiff` bedGraph family
- `bdgpeakcall`, `bdgbroadcall`, gzipped input throughout
- every `<subcommand> --help` (captured from the pinned oracle)

**Historical differences and current checks:**

| area | state |
|---|---|
| summit tie-break | The older recorded comparison had one differing peak per tested mode; the fresh real-data SE/BEDPE/BAMPE check now matches summit-related outputs byte-for-byte |
| `*_model.r` | The older correlation rounding gap is resolved; fresh predictd and source-SE model files are byte-identical against the pinned Haswell oracle |
| contig ordering | `filterdup` on 50k contigs: content identical, order differs (upstream hash-randomized, non-deterministic run to run) |
| `callvar` VCF | `##Program_Args` echoes the replay's outdir, otherwise identical on the recorded corpus |
| `hmmratac` | This older corpus predates self-training support. Gaussian and Poisson default training now run; fresh yeast oracle checks cover training and decoded output (see [status](docs/status.md)). |


## Input formats

Single-end `BED`, `BAM`, `SAM`, `ELAND`, `ELANDMULTI`, `ELANDEXPORT`, `BOWTIE`
and paired-end `BEDPE`, `BAMPE`, `FRAG` loaders are available. Pooled single-end
and BAM inputs do not require an index. Format-specific behavior and command
coverage still vary; see [status](docs/status.md). Parser details:

- Upstream's SAM parser crashes on any minus-strand read (`TypeError` in CIGAR
  parsing), making `-f SAM` effectively unusable upstream. This port parses SAM
  correctly and is byte-identical on inputs upstream can handle.
- Legacy parsing also preserves pinned upstream failures: ELANDMULTI rejects its
  bytes-to-integer conversion, and BOWTIE tag-size inference rejects fewer than ten successful
  records. ELAND and ELANDEXPORT remain usable.
- SE BAM 5' ends use the exclusive rightmost directly, matching
  `bam_fw_binary_parse`; BAMPE fragments use `abs(TLEN)` with leftmost-only
  proper pairs, matching `bampe_pe_binary_parse`. No `.bai` index is required.

## Performance

The table below is the original benchmark baseline. On the same real 5M-read
single-end CTCF workload (`--nomodel --extsize 200`, one thread), a newer Rust
run took 11.31 s and used 265,788 kB peak RSS; pinned upstream took 32.0 s and
used 297,628 kB. That's about 2.8x speed and 0.89x upstream RSS, so it does not
meet the release targets of at least 3x speed and at most 0.5x RSS on this workload.
The fresh summit check matches narrowPeak and summits byte-for-byte; the XLS
difference is only the output-directory command line header. Performance work
continues, and this single workload is not the full benchmark matrix. See
[status](docs/status.md) for current progress.

Measured on real CTCF ChIP-seq data (5.0 M treatment + 5.0 M control reads,
human `hs` genome size), median of 3, release build:

| workload | macs3 3.0.5 | macs3-rs | speedup | peak RSS |
|---|---|---|---|---|
| SE `--nomodel --extsize 200 -B` | 45.3 s / 297 MB | 13.3 s / 2967 MB | **3.4x** | **10.0x** |
| SE `--nomodel --SPMR` | 34.0 s / 294 MB | 4.8 s / 2304 MB | **7.0x** | **7.9x** |
| SE `--nomodel`, no control | 23.9 s / 171 MB | 7.2 s / 1391 MB | **3.3x** | **8.1x** |

At the time of this baseline, wall clock was 3–7x better and peak memory was
8–10x higher than upstream. Those RSS figures have since been substantially
reduced and must not be read as current measurements.

## Testing

```sh
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check
```

The differential corpora under `tests/` and the oracle harness under `oracle/`
replay recorded upstream invocations and compare every output byte. The oracle
itself is deliberately **not** vendored — vendoring it would make it too easy to
"fix" the oracle and silently invalidate every recorded result.

## Licence

BSD-3-Clause, matching upstream MACS3. `vendor/fermi-lite` is MACS3's own
submodule, vendored verbatim under its own MIT licence (`vendor/fermi-lite/LICENSE.txt`).
