#!/usr/bin/env python
"""Dump records as parsed by *upstream* MACS3, for record-level differential.

This deliberately goes through the real parser classes rather than
re-implementing them, so a disagreement with the Rust port localises to a
difference in the port. Where a parser's internal ``*_parse_line`` is a cdef
(and therefore not callable from Python), we drive the public entry point
(``build_fwtrack`` / ``build_bedtrack``) and read the resulting track.

Because upstream's track groups by chromosome and its arrays are position-keyed,
this dump is *sorted and grouped*, not in input order. The Rust side is compared
against it after sorting too; a separate in-order check is not meaningful here
because upstream discards input order by construction.

Usage:
    dump_oracle_records.py <bed|bedpe|frag|bedgraph> <file>
"""

import sys

from MACS3.IO.Parser import BEDParser, BEDPEParser, FragParser
from MACS3.IO.BedGraphIO import bedGraphTrackI


def dump_single_end(path):
    """(chrom, pos, strand) for every record upstream's BED reader produced."""
    parser = BEDParser(path, buffer_size=100000)
    track = parser.build_fwtrack()
    # FixWidthTrack stores each strand in a numpy array pre-sized to
    # `buffer_size` and only truncates to the fill count in `finalize()`. Without
    # this call every chromosome reports 2*buffer_size records, all but the real
    # ones being zero padding.
    track.finalize()
    out = []
    for chrom in track.get_chr_names():
        plus, minus = track.get_locations_by_chr(chrom)
        for p in plus:
            out.append((chrom, int(p), b"+"))
        for p in minus:
            out.append((chrom, int(p), b"-"))
    out.sort()
    for chrom, pos, strand in out:
        # `chrom` is bytes; decode so the output matches the Rust side exactly
        name = chrom.decode() if isinstance(chrom, bytes) else chrom
        sign = strand.decode() if isinstance(strand, bytes) else strand
        sys.stdout.write(f"{name}\t{pos}\t{sign}\n")


def dump_fragments(path, parser_cls, with_barcode):
    """(chrom, left, right[, barcode, count]) for every paired-end record.

    Driven through ``build_petrack``, the public entry point. Upstream stores
    fragments in a structured numpy array with fields ``l``/``r`` (and ``c`` for
    the FRAG count) that is pre-sized to ``buffer_size`` and only truncated by
    ``finalize()``, so that call is required before reading.
    """
    parser = parser_cls(path, buffer_size=100000)
    track = parser.build_petrack()
    track.finalize()
    out = []
    for chrom in track.get_chr_names():
        locs = track.get_locations_by_chr(chrom)
        fields = locs.dtype.names or ()
        for rec in locs:
            row = [chrom, int(rec["l"]), int(rec["r"])]
            if with_barcode:
                # barcode is not retained on the track itself; it is only used to
                # key the fragment counts, so it cannot be recovered here.
                row += [int(rec["c"]) if "c" in fields else 1]
            out.append(tuple(row))
    out.sort()
    for row in out:
        sys.stdout.write(
            "\t".join(x.decode() if isinstance(x, bytes) else str(x) for x in row) + "\n"
        )


def dump_bedgraph(path):
    """(chrom, start, end, value) from upstream's bedGraph reader."""
    track = bedGraphTrackI()
    track.read_bedGraph(baseline_value=0)
    out = []
    for chrom in track.chromosome.keys():
        iv = track.value[chrom]
        for start in iv:
            for end, value in iv[start].items():
                out.append((chrom, int(start), int(end), float(value)))
    out.sort()
    for chrom, start, end, value in out:
        name = chrom.decode() if isinstance(chrom, bytes) else chrom
        # repr round-trips exactly for doubles, unlike str() on some values
        sys.stdout.write(f"{name}\t{start}\t{end}\t{value!r}\n")


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: dump_oracle_records.py <bed|bedpe|frag|bedgraph> <file>")
    fmt, path = sys.argv[1], sys.argv[2]
    if fmt == "bed":
        dump_single_end(path)
    elif fmt == "bedpe":
        dump_fragments(path, BEDPEParser, with_barcode=False)
    elif fmt == "frag":
        dump_fragments(path, FragParser, with_barcode=True)
    elif fmt == "bedgraph":
        dump_bedgraph(path)
    else:
        sys.exit(f"unknown format {fmt}")


if __name__ == "__main__":
    main()
