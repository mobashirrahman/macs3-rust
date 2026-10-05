#!/usr/bin/env python3
"""Dump `BAMaccessor.get_reads_in_region` output as a stable text table.

The Rust port's BAM reader is checked against this file. Everything the reader
has to get right shows up here: which records the flag/MAPQ filter keeps, the
`lpos`/`rpos` pair (`rightmost` sums only reference-consuming CIGAR ops), the
strand, the packed sequence and quality, the CIGAR, and the `MD` tag -- plus the
duplicate-suppression order and count.

Usage:
    dump_bam_region.py BAM CHROM LEFT RIGHT MAXDUP

Output: one tab-separated line per read, in the order returned.
"""

import os
import sys

sys.path.insert(0, os.environ.get("MACS3_SRC", "/scratch/mdra00001/tmp/opencode/macs3-src"))

from MACS3.IO.BAM import BAMaccessor  # noqa: E402


def main():
    bam, chrom, left, right, maxdup = (
        sys.argv[1],
        sys.argv[2].encode(),
        int(sys.argv[3]),
        int(sys.argv[4]),
        int(sys.argv[5]),
    )
    b = BAMaccessor(bam)
    out = sys.stdout
    out.write("### references\t%s\n" % ",".join(c.decode() for c in b.get_chromosomes()))
    out.write("### rlengths\t%s\n" % ",".join(str(v) for v in b.get_rlengths().values()))
    reads = b.get_reads_in_region(chrom, left, right, maxDuplicate=maxdup)
    out.write("### n_reads\t%d\n" % len(reads))
    for r in reads:
        # `binaryseq` is the raw 4-bit-packed bytes, `binaryqual` the raw Phred
        # scores; `SEQ`/`QUAL` are the decoded forms. Both are dumped so a port
        # that decodes correctly but packs incorrectly is still caught.
        out.write(
            "\t".join(
                [
                    r["readname"].decode(),
                    r["chrom"].decode(),
                    str(r["lpos"]),
                    str(r["rpos"]),
                    str(r["strand"]),
                    r["binaryseq"].hex(),
                    r["binaryqual"].hex(),
                    ",".join(str(x) for x in r["cigar"]),
                    r["MD"].decode("latin-1"),
                    r["SEQ"].decode("latin-1"),
                    str(r["n_edits"]),
                    str(r["l"]),
                ]
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
