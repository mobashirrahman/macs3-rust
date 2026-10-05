//! Upstream's XLS, narrowPeak and summits writers.
//!
//! These are transcription of `writePeak2bedgraph` / `_write_peaks` in
//! `MACS3/Commands/callpeak_cmd.py`. Two conventions matter and are handled
//! separately, as the numerical contract requires:
//!
//! * **XLS** is 1-based inclusive: `start` is the chunk `tstart + 1` (F85) and
//!   `length` is `end - start + 1`.
//! * **narrowPeak** is 0-based half-open: `chromStart = xls start - 1`,
//!   `chromEnd = xls end`.
//!
//! Upstream formats the score columns with `%.5g`, so the byte output depends on
//! getting the formatting right rather than merely the values.
//!
//! # Peak coordinates are signed
//!
//! Every input coordinate -- read, fragment, pileup, `d` -- is a `u64`, and so is
//! every track position. **Peak** boundaries are `i64` here, because upstream can
//! emit a negative one: `--format FRAG` projects the control as `d`-wide windows
//! centred on each fragment *end*, `x - d//2`, with **no clipping**
//! (`pileup_from_LRC_centers_as_list`). With `--largelocal 10000` a fragment at
//! position 0 contributes a window starting at `-5000`, and a peak can start
//! before the contig -- `gmini_mfrag_d1200_w600_ctrl_069` has an xls `start` of
//! `-13`. Reproducing that is the whole point of the row, so the field cannot be
//! unsigned; reads and fragments stay `u64` throughout.

use macs_core::Result;

/// One peak as written to the XLS.
#[derive(Debug, Clone, PartialEq)]
pub struct XlsRow {
    /// Chromosome name.
    pub chrom: String,
    /// 1-based inclusive start.
    pub start: i64,
    /// 1-based inclusive end.
    pub end: i64,
    /// Absolute summit coordinate.
    pub summit: i64,
    /// Pileup at the summit.
    pub pileup: f32,
    /// `-log10(p)` at the summit.
    pub pscore: f32,
    /// Fold enrichment.
    pub fold_change: f64,
    /// `-log10(q)` at the summit.
    pub qscore: f32,
}

/// The score columns, which upstream writes with **`%.6g`**, not `%.5g`
/// (`MACS3/IO/PeakIO.py:815-840`). Six digits is what distinguishes, for example,
/// a fold enrichment of `7.20635` from `7.2064`; F119 measured that difference on
/// a real golden file.
fn g6(v: f32) -> String {
    format_g(v as f64, 6)
}

/// `pileup` is rounded to two places before formatting.
fn pileup6(v: f32) -> String {
    format_g(round2(v) as f64, 6)
}

/// Round half away from zero to two decimal places, as Python's `round(x, 2)`
/// does for the magnitudes a pileup takes (the binary-exact tie cases that
/// differ from `round-half-to-even` do not arise here, but the helper is named so
/// the intent is explicit).
fn round2(v: f32) -> f32 {
    let s = format!("{v:.2}");
    s.parse().unwrap_or(v)
}

/// narrowPeak/summits `score` is `int(10 * score_column)`; upstream selects the
/// column as `pscore` when `--log_pvalue` is given, else `qscore`
/// (`callpeak_cmd.py:293-298`). For the default qvalue call that is the qscore.
fn score_of(r: &XlsRow, score_col: ScoreCol) -> f32 {
    match score_col {
        ScoreCol::P => r.pscore,
        ScoreCol::Q => r.qscore,
    }
}

/// Which score column upstream writes into narrowPeak/summits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreCol {
    /// `--log_pvalue` was given.
    P,
    /// `--qvalue`/`--log_qvalue` (the default) was used.
    Q,
}

/// The block structure of one broad peak, as the gappedPeak format wants it.
///
/// `__add_broadpeak` (`BedGraph.py:600`): one block per strong (lvl1) sub-peak,
/// positioned relative to the broad region, plus a **1 bp** block on each side
/// when the strong set does not reach the broad region's own ends:
///
/// ```python
/// blockNum    = len(lvl1peakset)
/// blockSizes  = ",".join(x["length"] for x in lvl1peakset)
/// blockStarts = ",".join(x["start"] - start for x in lvl1peakset)
/// if int(thickStart) != start: blockNum += 1; blockSizes = "1," + blockSizes; blockStarts = "0," + blockStarts
/// if int(thickEnd)   != end:   blockNum += 1; blockSizes += ",1";          blockStarts += ",%d" % (end - start - 1)
/// ```
///
/// With no strong sub-peaks at all it still emits **two** blocks (`1,1` at `0`
/// and `end-start-1`), so a gappedPeak always has at least two.
#[derive(Debug, Clone, PartialEq)]
pub struct BroadBlocks {
    pub sizes: Vec<i64>,
    pub starts: Vec<i64>,
    /// `false` when upstream's `thickStart` is the sentinel `b"."`, which makes
    /// `write_to_gappedPeak` skip the peak entirely.
    pub present: bool,
}

/// Derive [`BroadBlocks`] from a broad region's span and its strong sub-peaks.
///
/// `start`/`end` are the broadPeak convention (0-based half-open), and `lvl1`
/// holds `(start, end)` of each strong sub-peak in coordinate order.
pub fn broad_blocks(start: i64, end: i64, lvl1: &[(i64, i64)]) -> BroadBlocks {
    let mut sizes: Vec<i64> = Vec::new();
    let mut starts: Vec<i64> = Vec::new();
    if lvl1.is_empty() {
        sizes = vec![1, 1];
        starts = vec![0, end - start - 1];
        return BroadBlocks {
            sizes,
            starts,
            present: true,
        };
    }
    // `thick_start`/`thick_end` are upstream's intermediate values; only their
    // *equality* with the broad region's own ends matters here (they decide
    // whether a 1 bp flank block is added), and the final values are not written.
    let thick_start0 = lvl1[0].0;
    let thick_end = lvl1[lvl1.len() - 1].1;
    for &(s, e) in lvl1 {
        sizes.push(e - s);
        starts.push(s - start);
    }
    if thick_start0 != start {
        sizes.insert(0, 1);
        starts.insert(0, 0);
    }
    if thick_end != end {
        sizes.push(1);
        starts.push(end - start - 1);
    }
    BroadBlocks {
        sizes,
        starts,
        present: true,
    }
}

/// `PeakIO.write_to_broadPeak` (`PeakIO.py:1470`).
///
/// Nine columns, no track line (`trackline=options.trackline`, off by default):
///
/// ```c
/// "%s\t%d\t%d\t%s%d\t%d\t.\t%.6g\t%.6g\t%.6g\n"
/// ```
///
/// Note the **column names are misleading**: the sixth field is `fc`, the seventh
/// `pscore`, the eighth `qscore`, and `score` is `int(10 * qscore)` -- so the
/// golden `905 . 21.2389 91.1085 90.5768` has `score = int(10 * 90.5768)`.
///
/// Peaks are numbered per **group of equal `end`**, taking the first of each
/// group -- unlike the narrowPeak/XLS writers, which number every row.
pub fn broadpeak_body(rows: &[XlsRow], name: &str, score_col: ScoreCol) -> String {
    let mut out = String::new();
    let prefix = format!("{name}_peak_");
    let mut n_peak = 0usize;
    for g in group_by_end(rows) {
        n_peak += 1;
        let r = &g[0];
        out.push_str(&format!(
            "{}\t{}\t{}\t{}{}\t{}\t.\t{}\t{}\t{}\n",
            r.chrom,
            // the broadPeak/gappedPeak convention is 0-based half-open, so the
            // XLS `start` (already 1-based, F85) shifts down by one
            r.start - 1,
            r.end,
            prefix,
            n_peak,
            (10.0 * f64::from(score_of(r, score_col))) as i64,
            format_g(r.fold_change, 6),
            format_g(r.pscore as f64, 6),
            format_g(r.qscore as f64, 6)
        ));
    }
    out
}

/// `PeakIO.write_to_gappedPeak` (`PeakIO.py:1343`): BED12 + 3 columns, again with
/// no track line.
///
/// The gate is `peak["thickStart"] != b"."`, so a peak with no block structure is
/// **skipped entirely** rather than written with zeros.
pub fn gappedpeak_body(
    rows: &[XlsRow],
    blocks: &[BroadBlocks],
    name: &str,
    score_col: ScoreCol,
) -> String {
    let mut out = String::new();
    let prefix = format!("{name}_peak_");
    for (i, (r, b)) in rows.iter().zip(blocks).enumerate() {
        if !b.present {
            continue;
        }
        out.push_str(&format!(
            "{}\t{}\t{}\t{}{}\t{}\t.\t0\t0\t0\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.chrom,
            r.start - 1,
            r.end,
            prefix,
            i + 1,
            (10.0 * f64::from(score_of(r, score_col))) as i64,
            b.sizes.len(),
            b.sizes
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(","),
            b.starts
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(","),
            format_g(r.fold_change, 6),
            format_g(r.pscore as f64, 6),
            format_g(r.qscore as f64, 6)
        ));
    }
    out
}

/// `itertools.groupby(peaks[chrom], key=itemgetter("end"))`: consecutive rows with
/// the same `end`, in order.
fn group_by_end(rows: &[XlsRow]) -> Vec<&[XlsRow]> {
    let mut out: Vec<&[XlsRow]> = Vec::new();
    let mut i = 0usize;
    while i < rows.len() {
        let mut j = i + 1;
        while j < rows.len() && rows[j].end == rows[i].end {
            j += 1;
        }
        out.push(&rows[i..j]);
        i = j;
    }
    out
}

fn subpeak_letters(i: usize) -> String {
    if i < 26 {
        ((b'a' + i as u8) as char).to_string()
    } else {
        format!(
            "{}{}",
            subpeak_letters(i / 26),
            (b'a' + (i % 26) as u8) as char
        )
    }
}

/// Render one narrowPeak line. 0-based half-open, so `chromStart = start - 1`.
///
/// Upstream (`PeakIO.py:761-770`):
/// `chrom start end name int(10*score) . %.6g(fc) %.6g(pscore) %.6g(qscore) (summit-start)`.
/// The summit column is the **offset** from the peak start, and none of the score
/// columns are signed.
pub fn narrowpeak_row(r: &XlsRow, name: &str, score_col: ScoreCol) -> String {
    let start0 = r.start.saturating_sub(1);
    let s = r.summit - r.start;
    format!(
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t{}",
        r.chrom,
        start0,
        r.end,
        name,
        (10.0 * f64::from(score_of(r, score_col))) as i64,
        format_g(r.fold_change, 6),
        g6(r.pscore),
        g6(r.qscore),
        s
    )
}

/// Render one summits.bed line: `chrom summit summit+1 name %.6g(score)`.
/// Upstream `_to_summits_bed` uses the 0-based summit, so it is `xls summit - 1`.
pub fn summit_row(r: &XlsRow, name: &str, score_col: ScoreCol) -> String {
    let summit0 = r.summit.saturating_sub(1);
    format!(
        "{}\t{}\t{}\t{}\t{}",
        r.chrom,
        summit0,
        r.summit,
        name,
        g6(score_of(r, score_col))
    )
}

/// Format with `%.<sig>g` semantics: `sig` significant digits, trailing zeros
/// removed, exponent form outside a wide range.
pub fn format_g(v: f64, sig: usize) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let exp = v.abs().log10().floor() as i32;
    if exp < -4 || exp >= sig as i32 {
        // Python switches to exponent form outside [1e-4, 1e{sig})
        let s = format!("{:.*e}", sig.saturating_sub(1), v);
        // Rust writes `1.2345e6`; Python writes `1.2345e+06`
        if let Some((m, e)) = s.split_once('e') {
            // `%g` drops trailing zeros from the fraction in exponent form too, so
            // `7.25020e-06` is `7.2502e-06` -- and a mantissa that is nothing but
            // zeros keeps its leading digit: `1.00000e-06` is `1e-06`.
            let m = if m.contains('.') {
                m.trim_end_matches('0').trim_end_matches('.')
            } else {
                m
            };
            let ev: i32 = e.parse().unwrap_or(0);
            return format!("{m}e{}{:02}", if ev < 0 { '-' } else { '+' }, ev.abs());
        }
        return s;
    }
    let decimals = (sig as i32 - 1 - exp).max(0) as usize;
    let s = format!("{v:.decimals$}");
    let s = if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    };
    if s == "-0" {
        "0".to_string()
    } else {
        s
    }
}

/// CPython's `repr(float)` -- what `str()`, `%s` and an f-string's empty `format()`
/// spec all produce for a Python `float`.
///
/// `hmmratac`'s `*_training_data.txt` is written as `f"{p[0]}\t{p[1]}\t{v[0]}..."`
/// (`hmmratac_cmd.py:350`), and `v[k]` there is `numpy.float32`. A `numpy.float32`
/// with an empty format spec widens to `double` and defers to `float.__repr__`, so
/// the column is the shortest round-tripping decimal of the **widened `f64`** --
/// not numpy's shorter, `float32`-round-tripping `str()` (`8.765152`).
///
/// Rust's `Display` is shortest-round-trip too, but not the same digit string: on
/// an exact halfway value it rounds away from zero where CPython's `_Py_dg_dtoa`
/// rounds to even, so `f64::from(8.765152f32)` prints `8.765151977539063` here and
/// `8.765151977539062` there -- 840 of 44757 rows of a real run. Rust's
/// fixed-precision `{:.n}` *is* correctly rounded with ties to even, so the answer
/// is the shortest digit count whose correctly rounded rendering parses back.
///
/// The layout is `format_float_short` (`Python/pystrtod.c`): exponent form when
/// the decimal point sits at `-4` or below or above `16`, and `Py_DTSF_ADD_DOT_0`
/// puts a `.0` on an integral value.
pub fn python_repr(v: f64) -> String {
    if v.is_nan() {
        return "nan".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    // 17 significant digits always round-trip an `f64`, so the loop stops there.
    let mut sig = 1usize;
    while sig < 17 && format!("{v:.*e}", sig - 1).parse::<f64>() != Ok(v) {
        sig += 1;
    }
    let sci = format!("{v:.*e}", sig - 1);
    let (mant, exp) = sci
        .split_once('e')
        .expect("LowerExp always writes an exponent");
    let exp: i32 = exp.parse().expect("LowerExp exponent parses");
    let (sign, mant) = mant.strip_prefix('-').map_or(("", mant), |m| ("-", m));
    // `mant` is `d.ddd`; the digits alone are what the point position indexes.
    let digits = mant.replace('.', "");
    // The decimal point sits one place past the leading digit.
    let point = exp + 1;
    if point <= -4 || point > 16 {
        let rest = &digits[1..];
        return if rest.is_empty() {
            format!(
                "{sign}{digits}e{}{:02}",
                if exp < 0 { '-' } else { '+' },
                exp.abs()
            )
        } else {
            format!(
                "{sign}{}.{}e{}{:02}",
                &digits[..1],
                rest,
                if exp < 0 { '-' } else { '+' },
                exp.abs()
            )
        };
    }
    let mut out = String::from(sign);
    if point <= 0 {
        out.push_str("0.");
        for _ in 0..-point {
            out.push('0');
        }
        out.push_str(&digits);
    } else if (point as usize) >= sig {
        out.push_str(&digits);
        for _ in 0..(point as usize - sig) {
            out.push('0');
        }
        out.push_str(".0");
    } else {
        let point = point as usize;
        out.push_str(&digits[..point]);
        out.push('.');
        out.push_str(&digits[point..]);
    }
    out
}

/// The XLS column header, in upstream's order.
pub const XLS_HEADER: &str =
    "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname";

/// The broad-mode XLS header, which has **no** `abs_summit` column (F111).
pub const XLS_HEADER_BROAD: &str =
    "chr\tstart\tend\tlength\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname";

/// Render one XLS data row.
pub fn xls_row(r: &XlsRow, name: &str, broad: bool) -> String {
    let length = r.end - r.start + 1;
    if broad {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.chrom,
            r.start,
            r.end,
            length,
            pileup6(r.pileup),
            g6(r.pscore),
            format_g(r.fold_change, 6),
            g6(r.qscore),
            name
        )
    } else {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.chrom,
            r.start,
            r.end,
            length,
            r.summit,
            pileup6(r.pileup),
            g6(r.pscore),
            format_g(r.fold_change, 6),
            g6(r.qscore),
            name
        )
    }
}

/// Assemble the whole XLS body (header included) for a set of rows.
///
/// Peak numbering follows upstream's `write_to_xls`: peaks are grouped by end
/// within a chromosome, and a group of several (a `--call-summits` peak) is
/// numbered once with sub-peak letters, exactly like narrowPeak and summits.
pub fn xls_body(rows: &[XlsRow], name_prefix: &str, broad: bool) -> Result<String> {
    let mut out = String::new();
    out.push_str(if broad { XLS_HEADER_BROAD } else { XLS_HEADER });
    out.push('\n');
    let prefix = format!("{name_prefix}_peak_");
    let mut n_peak = 0usize;
    for g in peak_groups(rows) {
        n_peak += 1;
        if g.len() > 1 {
            for (i, r) in g.iter().enumerate() {
                out.push_str(&xls_row(
                    r,
                    &format!("{prefix}{n_peak}{}", subpeak_letters(i)),
                    broad,
                ));
                out.push('\n');
            }
        } else {
            out.push_str(&xls_row(g[0], &format!("{prefix}{n_peak}"), broad));
            out.push('\n');
        }
    }
    Ok(out)
}

/// Group the rows the way upstream's writers do: **by end coordinate**, within a
/// chromosome. `--call-summits` emits several peaks that share an end, and they
/// are numbered as one peak with sub-peak letters. Rows must already be ordered
/// by chromosome then start.
fn peak_groups(rows: &[XlsRow]) -> Vec<Vec<&XlsRow>> {
    let mut groups: Vec<Vec<&XlsRow>> = Vec::new();
    for r in rows {
        let same_chrom = groups
            .last()
            .and_then(|g| g.first())
            .is_some_and(|f| f.chrom == r.chrom);
        if same_chrom && groups.last().is_some_and(|g| g[0].end == r.end) {
            groups.last_mut().expect("non-empty").push(r);
        } else {
            groups.push(vec![r]);
        }
    }
    groups
}

/// Assemble a narrowPeak body, reproducing upstream's peak numbering and
/// sub-peak lettering (`PeakIO.py:724-771`).
pub fn narrowpeak_body(rows: &[XlsRow], name: &str, score_col: ScoreCol) -> String {
    let prefix = format!("{name}_peak_");
    let mut out = String::new();
    let mut n_peak = 0usize;
    for g in peak_groups(rows) {
        n_peak += 1;
        if g.len() > 1 {
            for (i, r) in g.iter().enumerate() {
                out.push_str(&narrowpeak_row(
                    r,
                    &format!("{prefix}{n_peak}{}", subpeak_letters(i)),
                    score_col,
                ));
                out.push('\n');
            }
        } else {
            out.push_str(&narrowpeak_row(
                g[0],
                &format!("{prefix}{n_peak}"),
                score_col,
            ));
            out.push('\n');
        }
    }
    out
}

/// Assemble a summits.bed body (`_to_summits_bed`, `PeakIO.py:577-627`).
pub fn summits_body(rows: &[XlsRow], name: &str, score_col: ScoreCol) -> String {
    let prefix = format!("{name}_peak_");
    let mut out = String::new();
    let mut n_peak = 0usize;
    for g in peak_groups(rows) {
        n_peak += 1;
        if g.len() > 1 {
            for (i, r) in g.iter().enumerate() {
                out.push_str(&summit_row(
                    r,
                    &format!("{prefix}{n_peak}{}", subpeak_letters(i)),
                    score_col,
                ));
                out.push('\n');
            }
        } else {
            out.push_str(&summit_row(g[0], &format!("{prefix}{n_peak}"), score_col));
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row() -> XlsRow {
        XlsRow {
            chrom: "chr1".into(),
            start: 100,
            end: 200,
            summit: 150,
            pileup: 9.0,
            pscore: 10.828,
            fold_change: 7.20635,
            qscore: 7.9977,
        }
    }

    #[test]
    fn xls_length_is_end_minus_start_plus_one() {
        let s = xls_row(&row(), "n", false);
        assert!(s.contains("\t101\t"), "length is inclusive: {s}");
    }

    #[test]
    fn narrowpeak_is_zero_based_half_open() {
        let s = narrowpeak_row(&row(), "n", ScoreCol::Q);
        let f: Vec<&str> = s.split('\t').collect();
        assert_eq!(f[1], "99", "chromStart is start-1");
        assert_eq!(f[2], "200", "chromEnd is the XLS end");
    }

    /// The first three data rows of the real golden files
    /// (`tests/golden/se_model/realistic/default/`), which pin the narrowPeak and
    /// summits conventions that a coordinate diff cannot see.
    fn golden_peak1() -> XlsRow {
        XlsRow {
            chrom: "chr20".into(),
            start: 7572,
            end: 8039,
            summit: 7732,
            pileup: 9.0,
            pscore: 10.828,
            fold_change: 7.20635,
            qscore: 7.9977,
        }
    }

    #[test]
    fn narrowpeak_row_is_byte_identical_to_golden() {
        // golden: chr20 7571 8039 realistic_default_peak_1 79 . 7.20635 10.828 7.9977 160
        let got = narrowpeak_row(&golden_peak1(), "realistic_default_peak_1", ScoreCol::Q);
        let want =
            "chr20\t7571\t8039\trealistic_default_peak_1\t79\t.\t7.20635\t10.828\t7.9977\t160";
        assert_eq!(got, want);
    }

    #[test]
    fn summit_row_is_byte_identical_to_golden() {
        // golden summits: chr20 7731 7732 realistic_default_peak_1 7.9977
        let got = summit_row(&golden_peak1(), "realistic_default_peak_1", ScoreCol::Q);
        let want = "chr20\t7731\t7732\trealistic_default_peak_1\t7.9977";
        assert_eq!(got, want);
    }

    #[test]
    fn narrowpeak_score_column_switches_with_score_col() {
        // pscore column => int(10*10.828)=108, qscore => int(10*7.9977)=79
        let p = narrowpeak_row(&golden_peak1(), "n", ScoreCol::P);
        assert_eq!(p.split('\t').nth(4).unwrap(), "108");
        let q = narrowpeak_row(&golden_peak1(), "n", ScoreCol::Q);
        assert_eq!(q.split('\t').nth(4).unwrap(), "79");
    }

    #[test]
    fn subpeak_letters_match_upstream() {
        // upstream `subpeak_letters` is recursive: i<26 -> chr(97+i), else
        // subpeak_letters(i//26) + chr(97 + i%26). Note this yields "ba" for 26,
        // not the "aa" its docstring implies -- the code is the contract.
        assert_eq!(subpeak_letters(0), "a");
        assert_eq!(subpeak_letters(25), "z");
        assert_eq!(subpeak_letters(26), "ba");
        assert_eq!(subpeak_letters(27), "bb");
        assert_eq!(subpeak_letters(51), "bz");
    }

    #[test]
    fn call_summits_sharing_an_end_are_lettered() {
        // two peaks with the same end (call-summits) become one numbered peak a/b
        let a = XlsRow {
            summit: 150,
            ..row()
        };
        let b = XlsRow {
            summit: 160,
            ..row()
        };
        let body = narrowpeak_body(&[a.clone(), b.clone()], "run", ScoreCol::Q);
        let names: Vec<&str> = body
            .lines()
            .map(|l| l.split('\t').nth(3).unwrap())
            .collect();
        assert_eq!(names, vec!["run_peak_1a", "run_peak_1b"]);
        let sbody = summits_body(&[a, b], "run", ScoreCol::Q);
        assert!(sbody.contains("run_peak_1a"));
        assert!(sbody.contains("run_peak_1b"));
    }

    #[test]
    fn the_broad_header_has_no_summit_column() {
        assert!(!XLS_HEADER_BROAD.contains("abs_summit"));
        assert!(XLS_HEADER.contains("abs_summit"));
    }

    #[test]
    fn significant_digit_formatting_matches_pythons_percent_g() {
        assert_eq!(format_g(0.0, 5), "0");
        assert_eq!(format_g(1.0, 5), "1");
        assert_eq!(format_g(7.20635, 5), "7.2063");
        assert_eq!(format_g(0.0000123456, 5), "1.2346e-05");
        assert_eq!(format_g(123456.0, 5), "1.2346e+05");
    }

    /// Every cell here is a `numpy.float32` from `hmmratac`'s `*_training_data.txt`
    /// that the pinned reference wrote, and each one is a value where Rust's own
    /// `Display` rounds the final digit the other way.
    #[test]
    fn python_repr_rounds_an_exact_halfway_value_to_even() {
        // 8.765151977539062 is `hmmratac_yeast500k_training_data.txt` line 88
        // (yeast_500k_SRR1822137.bam); the rest are the tie cells of
        // `gmini_mfrag_d400_w180_ctrl_033/treat.frag`.
        let cases: [(f32, &str); 4] = [
            (8.765152, "8.765151977539062"),
            (136.04337, "136.04336547851562"),
            (89.1022, "89.10220336914062"),
            (895.99835, "895.9983520507812"),
        ];
        for (f32_value, want) in cases {
            let v = f64::from(f32_value);
            assert_eq!(python_repr(v), want, "for {f32_value}f32");
            // The bug this reproduces: `Display` is shortest-round-trip too, but
            // rounds an exact halfway digit away from zero instead of to even.
            assert_ne!(
                format!("{v}"),
                want,
                "Display should differ for {f32_value}f32"
            );
        }
        // A cell that is not halfway still has to come out unchanged.
        assert_eq!(python_repr(f64::from(0.92852986f32)), "0.9285298585891724");
    }

    /// A cell that *is* integral still needs `.0` (`Py_DTSF_ADD_DOT_0`), and the
    /// exponent window of `format_float_short` is `[-4, 16)` on the position of
    /// the decimal point.
    #[test]
    fn python_repr_layout_matches_cpython() {
        assert_eq!(python_repr(0.0), "0.0");
        assert_eq!(python_repr(-0.0), "-0.0");
        assert_eq!(python_repr(3.0), "3.0");
        assert_eq!(python_repr(-2.5), "-2.5");
        assert_eq!(python_repr(0.0001), "0.0001");
        assert_eq!(python_repr(1e-5), "1e-05");
        assert_eq!(python_repr(1.5e-5), "1.5e-05");
        assert_eq!(python_repr(1e15), "1000000000000000.0");
        assert_eq!(python_repr(1e16), "1e+16");
        assert_eq!(
            python_repr(1.7976931348623157e308),
            "1.7976931348623157e+308"
        );
        assert_eq!(python_repr(5e-324), "5e-324");
        assert_eq!(python_repr(f64::NAN), "nan");
        assert_eq!(python_repr(f64::INFINITY), "inf");
        assert_eq!(python_repr(f64::NEG_INFINITY), "-inf");
    }

    /// Whatever the digit string, the result must parse back to the input.
    #[test]
    fn python_repr_always_round_trips() {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        for _ in 0..20_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let v = f64::from_bits(state);
            if !v.is_finite() {
                continue;
            }
            assert_eq!(python_repr(v).parse::<f64>(), Ok(v), "for {v:?}");
        }
    }
}
