#!/usr/bin/env python3
"""Write small BAM/BAI fixtures for the hmmratac/callvar differential harness.

Uses only the standard library: BGZF is gzip with a ``BC`` extra subfield, and
BAI is a fixed binary layout. Writing them here keeps the harness free of a
samtools dependency, which the oracle image does not have.

The alignments are deliberately *not* trivial -- they cover proper/improper
pairs, secondary and supplementary records, MAPQ 0 and 255, soft-clipped reads
(so ``rightmost`` differs from a naive span), reverse-strand reads, duplicated
alignments, and reads with long CIGARs -- because each of those is a filter rule
in ``MACS3.IO.BAM.__fw_binary_parse`` and a reader that silently mis-handles any
of them will differ from the oracle.

Usage:
    make_bam_fixtures.py OUTDIR
"""

import os
import struct
import sys
import zlib

SEQ_CODE = {c: i for i, c in enumerate("=ACMGRSVTWYHKDBN")}

# BAM CIGAR operation codes, in the order the spec assigns them. Note there is
# no code 6 (`P`, padding) in the reference-consuming set below.
CMATCH, CINS, CDEL, CREF_SKIP, CSOFT, CHARD, CPAD, CMATCH_EQ, CX = range(9)


def bgzf_block(payload: bytes) -> bytes:
    """One BGZF member wrapping `payload`."""
    co = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = co.compress(payload) + co.flush()
    bsize = len(cdata) + 25  # total member size minus one
    if bsize > 0xFFFF:
        raise ValueError("BGZF block too large; lower the block size")
    head = struct.pack(
        "<BBBBIBBHBBHH",
        0x1F, 0x8B, 0x08, 0x04,  # magic, deflate, FEXTRA
        0,                        # mtime
        0, 0xFF,                  # XFL, OS
        6,                        # XLEN
        ord("B"), ord("C"), 2, bsize,
    )
    return head + cdata + struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))


EOF_BGZF = bgzf_block(b"")


def reg2bins(rbeg, rend):
    starts = (0, 1, 9, 73, 585, 4681)
    out = []
    for start, shift in zip(starts, range(29, 13, -3)):
        i = (rbeg >> shift) if rbeg > 0 else 0
        j = (rend >> shift) if rend < (1 << 29) - 1 else ((1 << 29) - 1) >> shift
        out.extend(start + k for k in range(i, j + 1))
    return out


def encode_seq(seq):
    out = bytearray()
    for i in range(0, len(seq), 2):
        hi = SEQ_CODE.get(seq[i], 15)
        lo = SEQ_CODE.get(seq[i + 1], 15) if i + 1 < len(seq) else 0
        out.append((hi << 4) | lo)
    return bytes(out)


def encode_cigar(cigar):
    return b"".join(struct.pack("<I", (n << 4) | op) for n, op in cigar)


def alignment(ref_id, pos, flag, mapq, name, cigar, seq, qual, md, tags=b""):
    """One BAM alignment block body (no `block_size` prefix)."""
    name_b = name.encode() + b"\0"
    body = struct.pack("<ii", ref_id, pos)
    body += struct.pack("<BBH", len(name_b), mapq, 0)          # name len, mapq, bin
    body += struct.pack("<HH", len(cigar), flag)               # n_cigar, flag
    body += struct.pack("<i", len(seq))                        # l_seq
    body += struct.pack("<iii", -1, -1, 0)                     # next_refID/pos/tlen
    body += name_b
    body += encode_cigar(cigar)
    body += encode_seq(seq)
    body += bytes(qual)
    body += b"MDZ" + md.encode() + b"\0"
    body += tags
    return body


def write_bam(path, header_text, refs, records, block_payload=8000):
    """Write a coordinate-sorted BAM. `records` is a list of `(ref_id, pos, body)`.

    The header gets its own BGZF block, padded to a 4-byte multiple, exactly as
    samtools does. That is not cosmetic: a BAI chunk whose compressed offset is 0
    is indistinguishable from "no chunks at all" -- both `get_coffset_by_region`
    and this port return early on it -- so an alignment must never start in block
    0.
    """
    head = bytearray()
    head += b"BAM\1"
    hb = header_text.encode()
    head += struct.pack("<i", len(hb)) + hb
    head += struct.pack("<i", len(refs))
    for name, length in refs:
        nb = name.encode() + b"\0"
        head += struct.pack("<i", len(nb)) + nb + struct.pack("<i", length)
    while len(head) % 4:
        head += b"\0"

    # A BAM record must not straddle a BGZF block boundary -- samtools flushes
    # at record boundaries for exactly this reason, and both this port and
    # upstream's reader would otherwise try to parse a record from a short slice.
    groups = []
    cur = []
    cur_len = 0
    for i, b in enumerate(records):
        need = 4 + len(b[2])
        if cur and cur_len + need > block_payload:
            groups.append(cur)
            cur, cur_len = [], 0
        cur.append(i)
        cur_len += need
    if cur or not groups:
        groups.append(cur)

    compressed = bytearray()
    compressed += bgzf_block(bytes(head))
    virtual_of = {}
    pos = 0
    for g in groups:
        coff = len(compressed)
        chunk = bytearray()
        for i in g:
            virtual_of[i] = (coff << 16) | pos
            chunk += struct.pack("<I", len(records[i][2])) + records[i][2]
            pos += 4 + len(records[i][2])
        compressed += bgzf_block(bytes(chunk))
    compressed += EOF_BGZF
    with open(path, "wb") as fh:
        fh.write(bytes(compressed))
    return virtual_of


def reg2bin(beg, end):
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


def write_bai(path, refs, entries, virtual_of):
    """`entries` is the same `(ref_id, span, body)` list `write_bam` reports, with
    `span` already computed from the reference-consuming CIGAR ops."""
    per_ref = {i: [] for i in range(len(refs))}
    for idx, (ref_id, pos, body) in enumerate(entries):
        per_ref[ref_id].append((pos, virtual_of[idx]))

    out = bytearray()
    out += b"BAI\1"
    out += struct.pack("<I", len(refs))
    for ref_id in range(len(refs)):
        rs = sorted(per_ref[ref_id])
        bins = {}
        for pos, v in rs:
            b = reg2bin(pos, max(pos + 1, pos + 50))
            bins.setdefault(b, []).append(v)
        # Pseudo-bin 37450 must be present: `BAIFile.__load_bins` does
        # `self.bins[i].pop(37450)` unconditionally, so a BAI without it raises
        # KeyError on open. It carries (ref_beg, ref_end, n_mapped, n_unmapped)
        # as two chunk-shaped pairs.
        out += struct.pack("<I", len(bins) + 1)
        for b, vs in sorted(bins.items()):
            out += struct.pack("<I", b)
            # one chunk spanning the whole bin is valid and sufficient here
            out += struct.pack("<I", 1)
            out += struct.pack("<Qq", vs[0], vs[-1])
        # Pseudo-bin 37450 must be present: `BAIFile.__load_bins` does
        # `self.bins[i].pop(37450)` unconditionally, so a BAI without it raises
        # KeyError on open. Its two "chunks" are (ref_beg, ref_end) and
        # (n_mapped, n_unmapped).
        out += struct.pack("<I", 37450)
        out += struct.pack("<I", 2)
        out += struct.pack("<Qq", 0, 0)
        out += struct.pack("<Qq", len(rs), 0)
        n_intv = ((max((p for p, _ in rs), default=0) + 1) >> 14) + 1
        out += struct.pack("<I", n_intv)
        filled = []
        cur = 0
        for w in range(n_intv):
            while cur < len(rs) and (rs[cur][0] >> 14) <= w:
                cur += 1
            if cur < len(rs):
                filled.append(rs[cur][1])
            elif filled:
                filled.append(0xFFFFFFFFFFFFFFFF)
            else:
                filled.append(0)
        for v in filled:
            out += struct.pack("<Q", v)
    with open(path, "wb") as fh:
        fh.write(bytes(out))


HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20000\n@SQ\tSN:chr2\tLN:8000\n"
REFS = [("chr1", 20000), ("chr2", 8000)]


def build_records():
    """A pile of reads over a handful of peak-like intervals on both contigs.

    Coverage per position is varied on purpose so that variant calling has
    something to work with, and so that a reader which silently skips or
    re-orders records changes the answer.
    """
    recs = []

    def add(ref_id, pos, flag, mapq, name, cigar, seq, qual, md):
        recs.append((ref_id, pos, alignment(ref_id, pos, flag, mapq, name, cigar, seq, qual, md)))

    # --- chr1: a well-covered 200 bp "peak" at 1000-1200 -------------------
    base = 1000
    for i in range(200):
        off = base + i
        # varying depth: 6 copies in the middle, 2 at the edges
        depth = 6 if 40 <= i < 160 else 2
        for k in range(depth):
            # alternate a reference base and a single-base substitution every
            # 7th position so a caller has a variant to find
            if i % 7 == 3:
                seq = "A" * 40
                md = "20A20"
            elif i % 11 == 5:
                seq = "A" * 40
                md = "20C20"
            else:
                seq = "A" * 40
                md = "40"
            add(0, off, 0, 60, f"r1_{i}_{k}", [(40, CMATCH)], seq, [35] * 40, md)

    # a few reverse-strand reads
    for k in range(3):
        add(0, 1100 + k, 16, 60, f"rr_{k}", [(40, CMATCH)], "T" * 40, [35] * 40, "40")

    # --- chr1: filter cases the parser must get right ----------------------
    # unmapped (flag 4)
    add(0, 1300, 4, 60, "f_unmapped", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    # secondary (256)
    add(0, 1300, 256, 60, "f_secondary", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    # supplementary (2048)
    add(0, 1300, 2048, 60, "f_supp", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    # QC fail (512)
    add(0, 1300, 512, 60, "f_qcfail", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    # paired, mate unmapped (1|4)
    add(0, 1300, 1 | 4, 60, "f_mateunmapped", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    # paired, not a proper pair (1)
    add(0, 1300, 1, 60, "f_notproper", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    # paired, proper pair with mapped mate (1|2|64)
    add(0, 1300, 1 | 2 | 64, 60, "p_proper", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    # MAPQ 0 and MAPQ 255
    add(0, 1300, 0, 0, "f_mapq0", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    add(0, 1300, 0, 255, "f_mapq255", [(40, CMATCH)], "A" * 40, [35] * 40, "40")
    # soft clip: rightmost must NOT include the clipped bases
    add(0, 1400, 0, 60, "c_soft", [(5, CSOFT), (20, CMATCH), (15, CSOFT)], "A" * 40, [35] * 40, "20")
    # insertion: query-only, rightmost excludes it
    add(0, 1400, 0, 60, "c_ins", [(10, CMATCH), (3, CINS), (10, CMATCH)], "A" * 23, [35] * 23, "20")
    # deletion: reference-consuming, rightmost includes it
    add(0, 1400, 0, 60, "c_del", [(10, CMATCH), (3, CDEL), (10, CMATCH)], "A" * 20, [35] * 20, "10A3A10")
    # skip (N): reference-consuming
    add(0, 1400, 0, 60, "c_skip", [(10, CMATCH), (5, CREF_SKIP), (10, CMATCH)], "A" * 20, [35] * 20, "10A10")
    # hard clip: consumes neither
    add(0, 1400, 0, 60, "c_hard", [(5, CHARD), (20, CMATCH)], "A" * 20, [35] * 20, "20")

    # --- chr1: an exact duplicate run of 4 identical alignments -----------
    for k in range(4):
        add(0, 2000, 0, 60, f"dup_{k}", [(40, CMATCH)], "A" * 40, [35] * 40, "40")

    # --- chr2: a second, sparser peak --------------------------------------
    for i in range(120):
        off = 3000 + i
        depth = 5 if 30 <= i < 90 else 1
        for k in range(depth):
            add(1, off, 0, 60, f"c2_{i}_{k}", [(40, CMATCH)], "G" * 40, [35] * 40, "40")

    return recs


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(outdir, exist_ok=True)
    recs = build_records()
    # BAI needs the alignment *span* for reg2bin; recompute it here the same way
    # the reader does (reference-consuming CIGAR ops only).
    entries = []
    for ref_id, pos, body in recs:
        # recover the span from the body we just built
        n_cigar, _flag = struct.unpack("<HH", body[12:16])
        name_len = body[8]
        i = 32 + name_len
        span = pos
        for k in range(n_cigar):
            op = struct.unpack("<I", body[i + 4 * k: i + 4 * k + 4])[0]
            if op & 15 in (CMATCH, CDEL, CREF_SKIP, CMATCH_EQ, CX):
                span += op >> 4
        entries.append((ref_id, span, body))
    bam = os.path.join(outdir, "reads.bam")
    virtual_of = write_bam(
        bam, HEADER, REFS,
        [(r, max(p, 1), b) for (r, p, b) in entries],
    )
    write_bai(bam + ".bai", REFS, entries, virtual_of)
    n_kept = sum(1 for _, _, b in recs)
    print(f"wrote {bam} ({n_kept} records, {len(REFS)} refs) and {bam}.bai")


if __name__ == "__main__":
    main()