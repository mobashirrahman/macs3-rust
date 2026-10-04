#!/usr/bin/env python3
"""Generate the synthetic fixture corpus for the differential harness.

Every fixture is a small, fully specified miniature experiment: a handful of
contigs, known peak shapes, controlled read densities, controlled duplicate
rates, controlled local-bias spikes, and deliberate edge cases. They exist
because they are cheap to regenerate, cheap to diff, and can be made to hit one
specific behaviour each.

Layout written under `--out` (default `tests/fixtures`):

    <group>/<name>/
        genome.txt          effective genome size notes (not read by MACS)
        manifest.tsv        the generator parameters, for reproducibility
        treat.bed           single-end treatment
        ctrl.bed            single-end control            (SE fixtures)
        treat.bedpe         paired-end treatment
        ctrl.bedpe          paired-end control            (PE fixtures)
        treat.frag          FRAG treatment                (FRAG fixtures)
        expected/           hand-reviewed expectations    (optional)

The generator is deterministic: a fixture is a pure function of its name.

Usage:
    python3 oracle/gen_fixtures.py --out tests/fixtures
"""

import argparse
import hashlib
import os
import random
import shutil
import sys

# ------------------------------------------------------------------ genome
# Small contigs keep the fixtures fast while still exercising clipping,
# contig-end behaviour and chromosome ordering. Names deliberately include
# numeric, alphabetical and bare contigs so ordering bugs surface.
GENOMES = {
    "mini": [
        ("chr1", 20_000),
        ("chr2", 15_000),
        ("chr10", 8_000),
        ("chrX", 5_000),
        ("chrM", 200),
    ],
    "tiny": [("chrA", 2_000), ("chrB", 1_500)],
    "onechrom": [("chr1", 12_000)],
    # MACS3's `PeakModel` sets
    #     min_tags = round(total * lmfold * peaksize / gsize / 2)
    # so a miniature genome cannot be model-trained: with 100 paired peaks the
    # per-site depth can never exceed min_tags, because min_tags grows with the
    # total. Concretely, with `gz` and S sites of depth D: total = S*D, and the
    # requirement D > total*3000/gz gives gz > S*3000. S = 150 needs gz > 450 kb,
    # and MACS3 also demands S >= 100. `mid` is sized comfortably above that.
    "mid": [(f"chr{i}", 100_000) for i in range(1, 21)],
}


def peak_shape(name, n, rng):
    """Return `n` sorted positions whose density follows a named shape."""
    if name == "uniform":
        return sorted(rng.randrange(0, 1) for _ in range(0))  # unused
    raise ValueError(name)


def sample_around(center, spread, n, rlength, rng):
    out = []
    for _ in range(n):
        p = int(rng.gauss(center, spread))
        if 0 <= p < rlength:
            out.append(p)
    return sorted(out)


def gaussian_peak(center, spread, n, rlength, rng):
    return sample_around(center, spread, n, rlength, rng)


def narrow_peak(center, halfwidth, n, rlength, rng):
    out = []
    for _ in range(n):
        p = center + rng.randrange(-halfwidth, halfwidth + 1)
        if 0 <= p < rlength:
            out.append(p)
    return sorted(out)


def broad_peak(center, halfwidth, n, rlength, rng):
    out = []
    for _ in range(n):
        # a plateau with soft shoulders
        u = rng.random()
        off = int(halfwidth * (u ** 0.35))
        p = center + rng.randrange(-off, off + 1)
        if 0 <= p < rlength:
            out.append(p)
    return sorted(out)


def bimodal_peak(c1, c2, sep, n, rlength, rng):
    half = n // 2
    return sorted(
        gaussian_peak(c1, sep / 4.0, half, rlength, rng)
        + gaussian_peak(c2, sep / 4.0, n - half, rlength, rng)
    )


def strand_peaks(rlength, n_sites, per_site, rng, halfwidth=25, dominant=True):
    """Many sharp, strand-symmetric spike sites.

    MACS3's `PeakModel` finds naive peaks separately in the + and - strand
    profiles and then *pairs* them to estimate the fragment length. A single
    broad Gaussian peak with random strand assignment does not produce enough
    strand peaks: upstream reports "Total number of paired peaks: 0" and exits 1.

    So a model-capable fixture needs many *sharp* sites, each with a dominant
    strand, so that the + profile has >= 3 peaks, the - profile has >= 3 peaks,
    and at least 100 of them pair up.
    """
    out = []
    step = max(2 * halfwidth + 10, rlength // (n_sites + 1))
    for i in range(n_sites):
        centre = (i + 1) * step
        if centre + 3 * halfwidth >= rlength:
            break
        # a sharp Gaussian profile, no background: this is what gives the
        # strand profiles their peaks
        out.extend(sample_around(centre, halfwidth / 2.0, per_site, rlength, rng))
    return sorted(out)


def strand_dominant_peaks(rlength, n_sites, per_site, rng, halfwidth=25):
    """Sites where each site is dominated by one strand.

    MACS3's `PeakModel` finds peaks separately in the + and - strand profiles and
    pairs them, on tags that have already been **de-duplicated** (default
    `--keep-dup 1` keeps one tag per position per strand). Two consequences
    shape this function, and both were found by running the real MACS3:

    * Depth must be spread over many *distinct* positions. A spike of 900 reads
      on 120 distinct positions collapses to 120 after de-duplication, which
      changes `treatment.total` by three orders of magnitude and therefore
      rescales `min_tags` and `max_tags` -- and the model then finds nothing.
    * The per-strand peak depth has to land inside
      `[min_tags, max_tags]`, where `min_tags = total*lmfold*peaksize/gz/2` and
      `max_tags` is the same with `umfold`. So `total` is bounded below by the
      required depth and above by `depth * gz / (umfold * peaksize / 2)`.

    Making alternate sites + and - dominant gives each profile peaks that can be
    paired; spreading the reads over a wide `halfwidth` keeps them distinct.
    """
    out = []
    step = max(2 * halfwidth + 10, rlength // (n_sites + 1))
    for i in range(n_sites):
        centre = (i + 1) * step
        if centre + 3 * halfwidth >= rlength:
            break
        n_plus = per_site if i % 2 == 0 else per_site // 6
        n_minus = per_site if i % 2 == 1 else per_site // 6
        out += sample_around(centre, halfwidth / 2.0, n_plus, rlength, rng)
        out += sample_around(centre, halfwidth / 2.0, n_minus, rlength, rng)
    return sorted(out)


def bias_spike(center, halfwidth, n, rlength, rng):
    """A deliberately non-uniform local bias, to make lambda_local bite."""
    out = []
    for _ in range(n):
        u = rng.random()
        if u < 0.9:
            p = center + rng.randrange(-halfwidth, halfwidth + 1)
        else:
            p = rng.randrange(0, rlength)
        if 0 <= p < rlength:
            out.append(p)
    return sorted(out)


def stable_seed(*parts):
    """A process-independent integer seed from the given parts.

    `hash()` cannot be used: CPython randomises `str`/`bytes` hashing per
    process unless PYTHONHASHSEED is set, so any fixture that seeded from a
    string produced different output on every invocation. A golden corpus that
    regenerates differently is worse than no corpus.
    """
    h = hashlib.sha256("|".join(str(p) for p in parts).encode()).digest()
    return int.from_bytes(h[:8], "big")


def write_bed(path, chrom, positions, rlength, name_prefix):
    with open(path, "w") as fh:
        for i, p in enumerate(positions):
            end = min(p + 200, rlength)
            fh.write(f"{chrom}\t{max(p, 0)}\t{end}\t{name_prefix}_read_{i}\t0\t+\n")


def write_bed_stranded(path, chrom, positions, rlength, name_prefix, rng):
    with open(path, "w") as fh:
        for i, p in enumerate(positions):
            end = min(p + 200, rlength)
            strand = "+" if rng.random() < 0.5 else "-"
            fh.write(f"{chrom}\t{max(p, 0)}\t{end}\t{name_prefix}_read_{i}\t0\t{strand}\n")


def write_bedpe(path, chrom, lefts, rights, name_prefix):
    with open(path, "w") as fh:
        for i, (l, r) in enumerate(zip(lefts, rights)):
            fh.write(f"{chrom}\t{max(l, 0)}\t{r}\t{name_prefix}_frag_{i}\n")


def write_frag(path, chrom, frags, barcodes, name_prefix, counts=None):
    """FRAG, as `MACS3.IO.Parser.FragParser.pe_parse_line` reads it.

    The column order is **chrom, start, end, barcode, count** -- barcode in
    column 4 and count in column 5, which is the opposite of the more obvious
    guess and the first thing that got this wrong. The count is a C `unsigned
    short`, capped at 65535.
    """
    if counts is None:
        counts = [1] * len(frags)
    with open(path, "w") as fh:
        for (l, r), bc, n in zip(frags, barcodes, counts):
            fh.write(f"{chrom}\t{l}\t{r}\t{bc}\t{n}\n")


def dups(seq, rate, rng):
    """Inflate a position list to a chosen duplicate rate."""
    out = list(seq)
    extra = int(len(seq) * rate)
    if extra and seq:
        out += [seq[rng.randrange(len(seq))] for _ in range(extra)]
    return sorted(out)


# ------------------------------------------------------------------ fixtures
def build_fixtures():
    """Yield `(group, name, genome_key, treat, ctrl, options)` descriptors.

    Each of `treat` / `ctrl` is a dict: chrom -> (positions_rle, n_total, dup_rate).
    """
    out = []

    def add(group, name, genome, treat, ctrl, **opts):
        # MACS3 refuses to build a fragment model from fewer than 100 paired
        # +/- strand peaks, and exits 1 when it cannot. Small fixtures are
        # therefore declared model-incapable and the runner must pass
        # `--nomodel --extsize` for them -- the Rust CLI has to reproduce the
        # same accept/reject behaviour, so this is part of the contract.
        opts.setdefault("model_capable", True)
        out.append((group, name, genome, treat, ctrl, opts))

    # ---- group: se_basic ------------------------------------------------
    rng = random.Random(1)
    g = GENOMES["mini"]
    treat = {}
    ctrl = {}
    for chrom, rlen in g:
        treat[chrom] = [gaussian_peak(rlen // 3, rlen / 40.0, 2500, rlen, rng),
                        gaussian_peak(2 * rlen // 3, rlen / 60.0, 1800, rlen, rng),
                        [rng.randrange(0, rlen) for _ in range(600)]]
        ctrl[chrom] = [[rng.randrange(0, rlen) for _ in range(1500)]]
    add("se_basic", "gauss_two_peaks", "mini", treat, ctrl, mode="se", model_capable=False)

    # ---- group: se_model ------------------------------------------------
    # the only fixtures that can build a fragment model
    gmid = GENOMES["mid"]
    rng = random.Random(2)
    t, c = {}, {}
    for chrom, rlen in gmid:
        main = gaussian_peak(rlen // 2, rlen / 50.0, 3000, rlen, rng)
        t[chrom] = [main, strand_dominant_peaks(rlen, 12, 12, rng, halfwidth=400)]
        c[chrom] = [[rng.randrange(0, rlen) for _ in range(2000)]]
    # NOT model-capable. See the note on `mid` and F17: MACS3's PeakModel
    # min/max tag thresholds are scaled by `total / gsize`, and the peaks are
    # found on *de-duplicated* tags, so a synthetic library either has too few
    # distinct positions to clear `min_tags` or too many for `max_tags`. The
    # model gate (G8) therefore uses MACS3's own CTCF test dataset; these
    # fixtures cover everything else.
    add("se_model", "sharp_spikes", "mid", t, c, mode="se", model_capable=False)

    rng = random.Random(3)
    t, c = {}, {}
    for chrom, rlen in gmid:
        t[chrom] = [strand_dominant_peaks(rlen, 12, 12, rng, halfwidth=400)]
        c[chrom] = [[rng.randrange(0, rlen) for _ in range(2000)]]
    add("se_model", "spikes_only", "mid", t, c, mode="se", model_capable=False)

    # a model fixture with realistic strand-symmetric peaks: this is what a real
    # ChIP-seq library looks like to PeakModel
    rng = random.Random(4)
    t, c = {}, {}
    for chrom, rlen in gmid:
        sites = strand_dominant_peaks(rlen, 12, 14, rng, halfwidth=400)
        t[chrom] = [sites, gaussian_peak(rlen // 2, rlen / 30.0, 2000, rlen, rng)]
        c[chrom] = [[rng.randrange(0, rlen) for _ in range(2000)]]
    add("se_model", "realistic", "mid", t, c, mode="se", model_capable=False)

    # ---- group: se_shapes ----------------------------------------------
    for shape, fn in (
        ("narrow", lambda rlen, rng: narrow_peak(rlen // 2, 40, 4000, rlen, rng)),
        ("broad", lambda rlen, rng: broad_peak(rlen // 2, rlen // 4, 6000, rlen, rng)),
        ("bimodal", lambda rlen, rng: bimodal_peak(rlen // 2 - 300, rlen // 2 + 300, 60, 5000, rlen, rng)),
        ("bias", lambda rlen, rng: bias_spike(rlen // 2, rlen // 8, 7000, rlen, rng)),
    ):
        # stable_seed, not hash(): see its docstring. hash() on a str is
        # randomised per process, so these four fixtures changed on every run.
        rng = random.Random(stable_seed("se_shapes", shape))
        t, c = {}, {}
        for chrom, rlen in g:
            t[chrom] = [fn(rlen, rng)]
            c[chrom] = [[rng.randrange(0, rlen) for _ in range(250)]]
        # a single broad peak cannot supply the >= 100 paired strand peaks the
        # model needs, so these run with --nomodel
        add("se_shapes", shape, "mini", t, c, mode="se", model_capable=False)

    # ---- group: se_dup --------------------------------------------------
    for rate in (0.0, 0.05, 0.25, 1.0, 5.0):
        rng = random.Random(1000 + int(rate * 100))
        t, c = {}, {}
        for chrom, rlen in g:
            t[chrom] = [dups(gaussian_peak(rlen // 2, rlen / 50.0, 3000, rlen, rng), rate, rng)]
            c[chrom] = [[rng.randrange(0, rlen) for _ in range(1200)]]
        add("se_dup", f"dup_rate_{rate}", "mini", t, c, mode="se", model_capable=False)

    # ---- group: se_edge -------------------------------------------------
    rng = random.Random(7)
    t = {"chr1": [[0, 1, 2, 19999]], "chrM": [[0, 1, 199]], "chr2": [[14999]]}
    c = {"chr1": [[rng.randrange(0, 20000)]], "chrM": [[rng.randrange(0, 200)]],
         "chr2": [[rng.randrange(0, 15000)]]}
    add("se_edge", "contig_edges", "mini", t, c, mode="se", model_capable=False)

    rng = random.Random(8)
    t = {"chr1": [[], []], "chr2": [[]]}
    c = {"chr1": [[rng.randrange(0, 20000)]], "chr2": [[rng.randrange(0, 15000)]]}
    add("se_edge", "one_sided_no_reads", "mini", t, c, mode="se", model_capable=False)

    # treatment on a chromosome the control lacks
    rng = random.Random(9)
    t = {"chr1": [gaussian_peak(8000, 200, 400, 20000, rng)]}
    c = {"chr1": [[rng.randrange(0, 20000) for _ in range(200)]]}
    add("se_edge", "disjoint_chromosomes", "mini", t, c, mode="se", model_capable=False)

    # a very shallow library
    rng = random.Random(10)
    t = {"chr1": [gaussian_peak(8000, 200, 12, 20000, rng)]}
    c = {"chr1": [[rng.randrange(0, 20000) for _ in range(5)]]}
    add("se_edge", "very_shallow", "mini", t, c, mode="se", model_capable=False)

    # no control at all
    rng = random.Random(11)
    t = {"chr1": [gaussian_peak(8000, 200, 3000, 20000, rng)],
         "chr2": [gaussian_peak(7000, 200, 2500, 15000, rng)]}
    add("se_edge", "no_control", "mini", t, {}, mode="se", no_control=True, model_capable=False)

    # ---- group: pe_basic -----------------------------------------------
    rng = random.Random(21)
    g2 = GENOMES["mini"]
    t, c = {}, {}
    for chrom, rlen in g2:
        def frags(center, n, rlen, rng, width=180):
            out = []
            for _ in range(n):
                l = max(0, int(rng.gauss(center, rlen / 40.0)))
                out.append((l, min(rlen, l + int(rng.gauss(width, 30)))))
            return sorted(out)

        t[chrom] = [frags(rlen // 2, 3000, rlen, rng)]
        c[chrom] = [frags(rlen // 2, 2500, rlen, rng)]
    add("pe_basic", "gauss_fragments", "mini", t, c, mode="pe")

    # ATAC-like: short fragments
    rng = random.Random(22)
    t, c = {}, {}
    for chrom, rlen in g2:
        def atac(center, n, rlen, rng):
            out = []
            for _ in range(n):
                l = max(0, int(rng.gauss(center, rlen / 50.0)))
                w = max(30, int(rng.gauss(80, 25)))
                out.append((l, min(rlen, l + w)))
            return sorted(out)

        t[chrom] = [atac(rlen // 2, 3000, rlen, rng)]
        c[chrom] = [atac(rlen // 2, 3200, rlen, rng)]
    add("pe_basic", "atac_short", "mini", t, c, mode="pe")

    # paired-end where the true fragment length is bimodal (nucleosome ladder)
    rng = random.Random(23)
    t, c = {}, {}
    for chrom, rlen in g2:
        def ladder(center, n, rlen, rng):
            out = []
            for _ in range(n):
                l = max(0, int(rng.gauss(center, rlen / 50.0)))
                w = rng.choice([60, 200, 380, 560])
                out.append((l, min(rlen, l + w)))
            return sorted(out)

        t[chrom] = [ladder(rlen // 2, 3000, rlen, rng)]
        c[chrom] = [ladder(rlen // 2, 3200, rlen, rng)]
    add("pe_basic", "nucleosome_ladder", "mini", t, c, mode="pe")

    # ---- group: frag_basic ---------------------------------------------
    rng = random.Random(31)
    t, c = {}, {}
    for chrom, rlen in g2:
        def scat(center, n, rlen, rng):
            out = []
            for _ in range(n):
                l = max(0, int(rng.gauss(center, rlen / 40.0)))
                out.append((l, min(rlen, l + int(rng.gauss(150, 40)))))
            return sorted(out)

        t[chrom] = [scat(rlen // 2, 2000, rlen, rng)]
        c[chrom] = [scat(rlen // 2, 2000, rlen, rng)]
    add("frag_basic", "barcode_fragments", "mini", t, c, mode="frag")

    # ---- group: frag_counts --------------------------------------------
    # `FragParser` stores the count in an `unsigned short` via a truncating C
    # cast, so 70000 becomes 4464 and -5 becomes 65531. The `except
    # OverflowError` branch that is supposed to cap at 65535 never fires, which
    # makes upstream's own "will be capped" warning dead code. Verified with
    # FragParser.build_petrack; see docs/upstream-findings.md F22.
    #
    # These fragments are deliberately tiny and non-model-shaped: the fixture
    # exists so the record-level differential pins the truncation, not so the
    # fragment model has something to fit.
    rng = random.Random(53)
    t, c = {}, {}
    for chrom, rlen in g2:
        # the `mini` genome's contigs are 200 bp, so the fragments have to fit
        # inside one or `BEDPE`/`FragParser` rejects them outright
        frag = lambda start, n: [(start + i * 7, start + i * 7 + 40) for i in range(n)]
        t[chrom] = [frag(10, 16)]
        c[chrom] = [frag(10, 16)]
    add("frag_counts", "out_of_range", "mini", t, c, mode="frag",
        model_capable=False, frag_counts=[1, 70000, -5, 0, 65535, 65536, 32768, 100])

    # ---- group: tiny ---------------------------------------------------
    rng = random.Random(41)
    g3 = GENOMES["tiny"]
    t, c = {}, {}
    for chrom, rlen in g3:
        t[chrom] = [gaussian_peak(rlen // 2, rlen / 20.0, 120, rlen, rng)]
        c[chrom] = [[rng.randrange(0, rlen) for _ in range(60)]]
    add("tiny", "two_contigs", "tiny", t, c, mode="se", model_capable=False)

    # ---- group: onechrom -----------------------------------------------
    rng = random.Random(51)
    g4 = GENOMES["onechrom"]
    t, c = {}, {}
    for chrom, rlen in g4:
        t[chrom] = [gaussian_peak(rlen // 2, rlen / 25.0, 500, rlen, rng)]
        c[chrom] = [[rng.randrange(0, rlen) for _ in range(250)]]
    add("onechrom", "single_contig", "onechrom", t, c, mode="se", model_capable=False)

    return out


# ------------------------------------------------------------------ writer

# Axes for the sweep group. Each is chosen to move a *different* stage of the
# pipeline rather than to vary noise: depth changes the duplicate threshold,
# width changes the peak width, ctrl_none removes the control branch entirely,
# and the mode selects the BED / BEDPE / FRAG reader.
SWEEP_GENOMES = ("mini", "tiny", "onechrom")
SWEEP_MODES = ("se", "pe", "frag")
SWEEP_DEPTHS = (400, 1200, 4000)
SWEEP_WIDTHS = (60, 180, 600)
SWEEP_CTRL = (True, False)


def sweep_fixtures(n):
    """Deterministic fixture sweep, for corpus breadth rather than for depth.

    The corpus has to be large enough that a systematic porting bug cannot hide
    in a parameter corner. Rather than a product of axes (which would explode),
    index `i` is decomposed into axis values, so any prefix of the sweep is
    itself balanced and the count is exact. Every fixture is a pure function of
    its index, so a failure is reproducible from the fixture name alone.

    Yields the same tuple shape as `build_fixtures`.
    """
    out = []
    n_axes = (len(SWEEP_GENOMES) * len(SWEEP_MODES) * len(SWEEP_DEPTHS) *
              len(SWEEP_WIDTHS) * len(SWEEP_CTRL))
    for i in range(n):
        # mixed-radix decomposition, so consecutive indices differ in the
        # fastest-varying axis and the sweep covers every combination evenly
        j = i
        g = SWEEP_GENOMES[j % len(SWEEP_GENOMES)]; j //= len(SWEEP_GENOMES)
        mode = SWEEP_MODES[j % len(SWEEP_MODES)]; j //= len(SWEEP_MODES)
        depth = SWEEP_DEPTHS[j % len(SWEEP_DEPTHS)]; j //= len(SWEEP_DEPTHS)
        width = SWEEP_WIDTHS[j % len(SWEEP_WIDTHS)]; j //= len(SWEEP_WIDTHS)
        have_ctrl = SWEEP_CTRL[j % len(SWEEP_CTRL)]

        rng = random.Random(10_000 + i)
        genome = GENOMES[g]
        treat, ctrl = {}, {}
        for chrom, rlen in genome:
            # two peaks plus background, so narrow and broad both have something
            # to find, at a depth that varies with the axis
            t = []
            for frac in (0.25, 0.7):
                c = int(rlen * frac)
                w = max(5, width // 4)
                t.append(sorted(max(0, min(rlen - 1, int(rng.gauss(c, w))))
                                for _ in range(depth)))
            t.append(sorted(rng.randrange(0, rlen) for _ in range(depth // 2)))
            treat[chrom] = t
            if have_ctrl:
                ctrl[chrom] = [[rng.randrange(0, rlen) for _ in range(depth)]]

        name = (f"g{g}_m{mode}_d{depth}_w{width}_"
                f"{'ctrl' if have_ctrl else 'noc'}_{i:03d}")
        opts = {"mode": mode, "model_capable": False}
        out.append(("sweep", name, g, treat, ctrl, opts))

    return out[:n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="tests/fixtures")
    ap.add_argument("--clean", action="store_true")
    ap.add_argument("--sweep", type=int, default=0,
                    help="also emit N deterministic sweep fixtures (the G0 "
                         "corpus must reach >=400 fixtures)")
    args = ap.parse_args()

    out = os.path.abspath(args.out)
    if args.clean and os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out, exist_ok=True)

    written = 0
    for group, name, genome_key, treat, ctrl, opts in (build_fixtures() +
                                                        sweep_fixtures(args.sweep)):
        d = os.path.join(out, group, name)
        os.makedirs(d, exist_ok=True)
        genome = GENOMES[genome_key]
        rlengths = dict(genome)
        mode = opts.get("mode", "se")

        with open(os.path.join(d, "genome.txt"), "w") as fh:
            fh.write("# effective genome size for this fixture\n")
            for chrom, rlen in genome:
                fh.write(f"{chrom}\t{rlen}\n")
            fh.write(f"# sum={sum(rl for _, rl in genome)}\n")

        manifest = [
            f"fixture\t{name}",
            f"group\t{group}",
            f"genome\t{genome_key}",
            f"mode\t{mode}",
            f"model_capable\t{int(bool(opts.get('model_capable', True)))}",
        ]

        for label, tracks in (("treat", treat), ("ctrl", ctrl)):
            for chrom, populations in tracks.items():
                rlen = rlengths[chrom]
                for i, positions in enumerate(populations):
                    suffix = "" if i == 0 else f"_{i}"
                    # a population is either a list of 5' positions (converted to
                    # fragments here) or a list of (l, r) fragment tuples
                    if positions and isinstance(positions[0], tuple):
                        frs = sorted((l, min(rlen, r)) for l, r in positions)
                    else:
                        frs = [(p, min(rlen, p + 180)) for p in sorted(positions)]
                    positions = frs if mode in ("pe", "frag") else sorted(positions)
                    if mode == "pe":
                        write_bedpe(
                            os.path.join(d, f"{label}{suffix}.bedpe"), chrom,
                            [l for l, _ in frs], [r for _, r in frs], f"{label}{suffix}")
                    elif mode == "frag":
                        barcodes = [f"BC_{l % 97:03d}" for l, _ in frs]
                        # a spread of multiplicities so --max-count has
                        # something to bite on, unless the fixture pins its own
                        # counts (see the frag_counts group)
                        counts = opts.get("frag_counts")
                        if counts is None:
                            counts = [1 + (i % 7) for i in range(len(frs))]
                        elif len(counts) > len(frs):
                            counts = [counts[i % len(counts)] for i in range(len(frs))]
                        write_frag(
                            os.path.join(d, f"{label}{suffix}.frag"), chrom, frs, barcodes,
                            f"{label}{suffix}", counts)
                    else:
                        # A stable digest, not `hash()`: CPython randomises
                        # string hashing per process (PYTHONHASHSEED), so `hash()`
                        # on a tuple containing strings gave different strand
                        # assignments on every run and the golden corpus was not
                        # reproducible at all.
                        rng = random.Random(
                            stable_seed(name, label, chrom, i))
                        write_bed_stranded(
                            os.path.join(d, f"{label}{suffix}.bed"), chrom, positions, rlen,
                            f"{label}{suffix}", rng)
                    manifest.append(f"{label}{suffix}\t{chrom}\t{len(positions)}")

        with open(os.path.join(d, "manifest.tsv"), "w") as fh:
            fh.write("\n".join(manifest) + "\n")
        written += 1

    sys.stderr.write(f"wrote {written} fixtures under {out}\n")


if __name__ == "__main__":
    main()
