# Compatibility

What has been compared against MACS3 3.0.5 (commit `c5443190`, NumPy 2.5.3),
command by command. "Real data" means MACS3's 5 M-read CTCF ChIP-seq and input
files; "small data" means the chr22 and yeast files MACS3 ships for its tests.

**These are one-off differential comparisons, not gates.** Only `callpeak` has a
recorded corpus that CI replays on every push; every row below it was produced by
running both implementations on the named input and diffing. That is weaker
evidence, and it has already earned its keep — three real defects in
`bdgpeakcall`, `bdgopt` and `bdgcmp` were found this way after the `callpeak`
corpus was already green. Recording these comparisons as a corpus is the top open
item.

| command | compared on | result |
|---|---|---|
| `callpeak` | 7,204 recorded runs; real data in narrow, broad, model, `-p`, `--call-summits`, `--min-length`, `--max-gap`, `--shift`, `--keep-dup`, `-B`, paired-end | byte-identical |
| `pileup` | real data, single-end; small data, paired-end | byte-identical |
| `filterdup` | real and small data | same rows; chromosome order differs on multi-chromosome input |
| `randsample` | small data, single chromosome | byte-identical; multi-chromosome sampling differs |
| `predictd` | real data (5 M-read CTCF, `--mfold 100 200` and the default); small data, paired-end | byte-identical, including `*_model.r` |
| `refinepeak` | real data | byte-identical |
| `bdgcmp` | real data, all eight methods | byte-identical |
| `bdgopt` | real data, `p2q` and `multiply`; `add`, `max`, `min` on a subset | byte-identical |
| `cmbreps` | real data, `max`, `mean`, `fisher` | byte-identical |
| `bdgpeakcall` | real data, including `--cutoff-analysis` | byte-identical |
| `bdgbroadcall` | real data, two settings | byte-identical |
| `bdgdiff` | real data | byte-identical |
| `hmmratac` | yeast ATAC-seq, Gaussian and Poisson, BAM, BEDPE and fragment input | regions, states, signal tracks and cutoff analysis byte-identical; trained model parameters agree to ~8 × 10⁻¹² relative |
| `callvar` | small data, assembly off / auto / on | VCF records identical |

Command-line parsing follows argparse: attached short-option values (`-q0.05`),
unambiguous long-option prefixes, Python's numeric literals, `--`, and the same
accept/reject decision and exit status on invalid values.

## Differences that remain

- **`filterdup` and seeded `randsample` on many chromosomes.** MACS3 iterates
  chromosomes in Python hash order, which changes from run to run; this port
  uses sorted order. Rows are the same for `filterdup`; the sample differs for
  `randsample`.
- **`hmmratac` trained model parameters.** Ours agree with upstream's to about
  8 × 10⁻¹² relative; two runs of upstream agree with *each other* to about
  2 × 10⁻¹², so the gap is roughly 4× upstream's own reproducibility floor rather
  than zero. Measured over all 77 parameters of a freshly trained yeast model
  (see status). Nothing downstream moves: the decoded accessible regions, the
  states and the cutoff analysis are byte-identical.
- **`callvar` VCF header.** `##Program_Args` echoes the output path, so it
  differs whenever the path does.
- **SAM input.** MACS3's SAM parser fails on any minus-strand read; this port
  parses SAM and matches on inputs MACS3 can read.
- **Threads.** `MACS3_RS_THREADS` sets the thread count. There is no
  `--threads` flag because MACS3 rejects one.

## Input formats

Single-end `BED`, `BAM`, `SAM`, `ELAND`, `ELANDMULTI`, `ELANDEXPORT`, `BOWTIE`
and paired-end `BEDPE`, `BAMPE`, `FRAG` loaders are available. Pooled single-end
and BAM inputs do not require an index. Parser details:

- Upstream's SAM parser crashes on any minus-strand read (`TypeError` in CIGAR
  parsing), making `-f SAM` effectively unusable upstream. This port parses SAM
  correctly and is byte-identical on inputs upstream can handle.
- Legacy parsing also preserves pinned upstream failures: ELANDMULTI rejects its
  bytes-to-integer conversion, and BOWTIE tag-size inference rejects fewer than ten successful
  records. ELAND and ELANDEXPORT remain usable.
- SE BAM 5' ends use the exclusive rightmost directly, matching
  `bam_fw_binary_parse`; BAMPE fragments use `abs(TLEN)` with leftmost-only
  proper pairs, matching `bampe_pe_binary_parse`. No `.bai` index is required.

## Where the memory goes

* `callpeak` SE: ~43 MB is the two resident position arrays, then one
  chromosome's signal build plus the q-table histogram — the parallel window's
  working set. `CHUNK` chromosomes are built at once and the window is the
  memory/speed dial.
* the local-lambda merge builds its track directly (`over_two_pv_array_track`),
  sized exactly up front; the old path held five tracks at once per fold.
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
