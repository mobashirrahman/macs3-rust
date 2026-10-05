# Compatibility

Verified differentially against upstream 3.0.5 on a generated 7685-invocation
`callpeak` matrix, upstream's `test/cmdlinetest`, and command-specific oracle
checks. This is not yet full parity across every flag of all 14 commands; the
current measured results and open differences are tracked in [status](status.md).

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
| `hmmratac` | This older corpus predates self-training support. Gaussian and Poisson default training now run; fresh yeast oracle checks cover training and decoded output (see [status](status.md)). |


## Input formats

Single-end `BED`, `BAM`, `SAM`, `ELAND`, `ELANDMULTI`, `ELANDEXPORT`, `BOWTIE`
and paired-end `BEDPE`, `BAMPE`, `FRAG` loaders are available. Pooled single-end
and BAM inputs do not require an index. Format-specific behavior and command
coverage still vary; see [status](status.md). Parser details:

- Upstream's SAM parser crashes on any minus-strand read (`TypeError` in CIGAR
  parsing), making `-f SAM` effectively unusable upstream. This port parses SAM
  correctly and is byte-identical on inputs upstream can handle.
- Legacy parsing also preserves pinned upstream failures: ELANDMULTI rejects its
  bytes-to-integer conversion, and BOWTIE tag-size inference rejects fewer than ten successful
  records. ELAND and ELANDEXPORT remain usable.
- SE BAM 5' ends use the exclusive rightmost directly, matching
  `bam_fw_binary_parse`; BAMPE fragments use `abs(TLEN)` with leftmost-only
  proper pairs, matching `bampe_pe_binary_parse`. No `.bai` index is required.

## Where the single-end memory goes

Measured on the 5 M-read CTCF fixture:

* ~43 MB is the two resident u32 position arrays (10 M reads);
* ~86 MB is one chromosome's signal build plus the q-table histogram, which is
  the parallel window's working set: `CHUNK` chromosomes are built at once and the
  window is the memory/speed dial. `CHUNK = 2` is the shipping setting; `CHUNK = 1`
  reaches 0.54x at the cost of the `-B` workload's speed (2.8x);
* the bedGraph and spool bodies are **streamed to per-chromosome temp files** and
  concatenated at the end, so `-B` no longer buffers a chromosome's ~38 MB of text;
* `MALLOC_ARENA_MAX` is capped to 1 by re-exec, because glibc's per-thread arenas
  retain freed blocks (237 MB against 177 MB at the same wall clock).
