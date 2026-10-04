#!/usr/bin/env python
"""Differential for `macs_track` duplicate filtering against upstream.

Drives upstream's own `FWTrack.filter_dup` and `PETrack.filter_dup` over randomly
generated tracks, and compares the retained positions and the reported total
against a table of expected values. The Rust side is checked by
`crates/macs-track/tests/parity.rs`, which reads the same table, so this script
only has to *produce* the oracle.

Each row records the *input* as well as the retained output: `filter_dup` mutates
in place, so the retained set alone cannot reconstruct the case the oracle ran.

Cases deliberately include the shapes where upstream's implementation has quirks:

* a strand holding exactly one position (bypasses filtering *and* the total);
* runs longer than `maxnum` at the same position;
* `maxnum` of 0;
* `maxnum` of -1 (no-op);
* fragments sharing a start but differing in end (not duplicates);
* chromosomes with a single fragment.

Usage:
    gen_track_vectors.py <out.tsv>
"""

import os
import random
import sys

import numpy as np

from MACS3.Signal.FixWidthTrack import FWTrack
from MACS3.Signal.PairedEndTrack import PETrackI

SEED = 20240611

# Sampling grids. 0.0 and 1.0 are the degenerate ends; the 0.333.. values make
# odd array lengths land on ties that the truncating cast resolves downward.
SAMPLE_PERCENTS = (0.0, 0.1, 0.25, 1.0 / 3.0, 0.5, 0.75, 0.9, 1.0)
SAMPLE_SEEDS = (0, 1, 7, 42, 12345)
SAMPLE_SIZES = (0, 1, 5, 10, 1000)


def make_se_track(rng, chroms):
    """A random single-end track, biased towards duplicates."""
    t = FWTrack(buffer_size=100000)
    for chrom in chroms:
        for _ in range(rng.randint(1, 40)):
            pos = rng.choice([10, 10, 10, 20, 20, 30, 40, 50, rng.randrange(1, 200)])
            strand = rng.choice([0, 1])
            t.add_loc(chrom, pos, strand)
    t.finalize()
    return t


def make_pe_track(rng, chroms):
    t = PETrackI(buffer_size=100000)
    for chrom in chroms:
        for _ in range(rng.randint(1, 30)):
            l = rng.choice([10, 10, 20, 30, rng.randrange(1, 150)])
            r = l + rng.choice([10, 20, 20, 50])
            t.add_loc(chrom, l, r)
    t.finalize()
    return t


def se_strand_counts(t):
    """Per-strand survivor counts, so a consumer never has to infer them."""
    plus = minus = 0
    for chrom in t.get_chr_names():
        p, m = t.get_locations_by_chr(chrom)
        plus += len(p)
        minus += len(m)
    return plus, minus


def read_se(t):
    """Every position currently held, rendered as `chrom:strand:pos`.

    Returns "" when nothing is held, so the TSV field is empty rather than a
    stray separator. A reader cannot distinguish "" from a missing column by
    counting fields alone, so the parser dispatches on the row's `kind`.
    """
    out = []
    for chrom in sorted(t.get_chr_names()):
        plus, minus = t.get_locations_by_chr(chrom)
        for p in plus:
            out.append(f"{chrom.decode()}:+:{int(p)}")
        for p in minus:
            out.append(f"{chrom.decode()}:-:{int(p)}")
    return out


def read_pe(t):
    out = []
    for chrom in sorted(t.get_chr_names()):
        for rec in t.get_locations_by_chr(chrom):
            out.append(f"{chrom.decode()}:{int(rec['l'])}-{int(rec['r'])}")
    return out


def se_row(rng, chroms, maxnum):
    t = make_se_track(rng, chroms)
    before = int(t.total)
    # the *input* must be recorded too: `filter_dup` mutates in place, so the
    # retained set alone cannot be used to reconstruct the case the oracle ran
    orig = read_se(t)
    reported = int(t.filter_dup(maxnum))
    return ["se", str(maxnum), str(before), str(reported),
            ";".join(orig), ";".join(read_se(t))]


def pe_row(rng, chroms, maxnum):
    t = make_pe_track(rng, chroms)
    before = int(t.total)
    orig = read_pe(t)
    reported = t.filter_dup(maxnum)
    reported = int(reported) if reported is not None else int(t.total)
    return ["pe", str(maxnum), str(before), str(reported),
            ";".join(orig), ";".join(read_pe(t))]


def sample_se_row(rng, chroms, percent, seed):
    """`FWTrack.sample_percent`: record input and the surviving positions."""
    t = make_se_track(rng, chroms)
    orig = read_se(t)
    t.sample_percent(percent, seed)
    # `repr` round-trips the exact double. A fixed-width decimal would not:
    # 1/3 written as "0.333333" is a *different* value from 1/3, and
    # `int(round(18 * 0.333333, 5))` is 5 where `int(round(18 * (1/3), 5))` is 6.
    kp, km = se_strand_counts(t)
    return ["se_sample", repr(percent), str(seed), str(len(orig)),
            str(int(t.total)), str(kp), str(km),
            ";".join(orig), ";".join(read_se(t))]


def sample_num_se_row(rng, chroms, samplesize, seed):
    """`FWTrack.sample_num`: a target count, which upstream turns into a fraction."""
    t = make_se_track(rng, chroms)
    orig = read_se(t)
    t.sample_num(samplesize, seed)
    kp, km = se_strand_counts(t)
    return ["se_sample_num", str(samplesize), str(seed), str(len(orig)),
            str(int(t.total)), str(kp), str(km),
            ";".join(orig), ";".join(read_se(t))]


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else "crates/macs-track/tests/track_vectors.tsv"
    rng = random.Random(SEED)
    chrom_sets = [
        [b"chr1"],
        [b"chr1", b"chr2"],
        [b"chrA", b"chrB", b"chrC"],
    ]
    # -1 is a no-op, 0 and 1 hit the quirks, 2..5 exercise ordinary truncation
    maxnums = [-1, 0, 1, 2, 3, 5]

    rows = []
    for i in range(240):
        chroms = chrom_sets[i % len(chrom_sets)]
        maxnum = maxnums[i % len(maxnums)]
        rows.append(se_row(rng, chroms, maxnum))
    for i in range(180):
        chroms = chrom_sets[i % len(chrom_sets)]
        maxnum = maxnums[i % len(maxnums)]
        rows.append(pe_row(rng, chroms, maxnum))

    # Down-sampling. Percentages chosen to hit the awkward cases: 0 and 1 (keep
    # nothing / everything), and values whose product with an odd length lands on
    # a .5 that the truncating cast resolves downward.
    #
    # Single-chromosome inputs only. Upstream's `get_chr_names` returns an
    # unordered `set` (F29), so with two or more chromosomes the order arrays are
    # shuffled in -- and therefore which draws come from the shared RNG stream --
    # depends on CPython's per-process string hash seed. Such a case is not
    # reproducible in upstream itself and cannot be a parity target. The
    # multi-chromosome behaviour is covered by the unit test that asserts the port
    # uses sorted byte order, and the limitation is recorded in F29.
    single = [[b"chr1"]]
    for i in range(150):
        pct = SAMPLE_PERCENTS[i % len(SAMPLE_PERCENTS)]
        seed = SAMPLE_SEEDS[i % len(SAMPLE_SEEDS)]
        rows.append(sample_se_row(rng, single[0], pct, seed))
    for i in range(90):
        size = SAMPLE_SIZES[i % len(SAMPLE_SIZES)]
        seed = SAMPLE_SEEDS[i % len(SAMPLE_SEEDS)]
        rows.append(sample_num_se_row(rng, single[0], size, seed))

    # filter rows:    kind, maxnum, before, reported, orig, kept
    # sampling rows:  kind, arg1, seed, before, reported, orig, kept
    header = "kind\targ1\targ2\tbefore\treported\torig\tkept"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"wrote {len(rows)} track vectors to {out_path}")


if __name__ == "__main__":
    main()
