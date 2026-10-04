//! The comparison logic for `macs-compare`.
//!
//! This is the project's primary development instrument. Where a text `diff`
//! reports "line 4 differs", this reports:
//!
//! ```text
//! se_shapes/narrow/default  narrowPeak
//!   peak count          MACS3 4   macs3-rs 4
//!   exact coordinates   4 / 4
//!   summits             4 / 4   max shift 0 bp
//!   p-value             max |delta| 0.0        (limit 1e-9)
//!   q-value             max |delta| 3.2e-06    (limit 1e-6)
//!   bytes               identical
//! ```
//!
//! and, when something does differ, the *first* divergent record with its
//! genomic context, so a failing gate says which stage broke rather than just
//! that something did.

use std::fmt;

/// A parsed genomic record, with every field MACS emits.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// Raw columns, split on tab.
    pub fields: Vec<String>,
    /// Chromosome, column 0.
    pub chrom: String,
    /// Start, column 1.
    pub start: u64,
    /// End, column 2.
    pub end: u64,
    /// Column 3: peak name, or the score in bedGraph.
    pub name: String,
    /// Column 4, when present.
    pub score: f64,
    /// Column 5, when present.
    pub f5: f64,
    /// Column 6, when present.
    pub f6: f64,
    /// Column 7, when present.
    pub f7: f64,
    /// Column 8, when present.
    pub f8: f64,
    /// Column 9, when present.
    pub f9: f64,
}

impl Record {
    /// Parse one line of a MACS output file.
    ///
    /// Non-numeric or missing trailing columns become `f64::NAN` rather than an
    /// error, because a format may legitimately stop early (a 4-column
    /// `summits.bed` has no score).
    pub fn parse(line: &str) -> Option<Record> {
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() || line.starts_with('#') || line.starts_with("track ") {
            return None;
        }
        // bedGraph has a `track` line and possibly a `browser` line
        if line.starts_with("browser ") {
            return None;
        }
        let f: Vec<String> = line.split('\t').map(str::to_string).collect();
        if f.len() < 3 {
            return None;
        }
        let num = |i: usize| {
            f.get(i)
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(f64::NAN)
        };
        Some(Record {
            start: f[1].parse().unwrap_or(0),
            end: f[2].parse().unwrap_or(0),
            chrom: f[0].clone(),
            name: f.get(3).cloned().unwrap_or_default(),
            score: num(4),
            f5: num(5),
            f6: num(6),
            f7: num(7),
            f8: num(8),
            f9: num(9),
            fields: f,
        })
    }

    /// The numeric columns, paired with the label used in reports.
    pub fn numeric_columns(&self) -> Vec<(&'static str, f64)> {
        vec![
            ("col4", self.score),
            ("col5", self.f5),
            ("col6", self.f6),
            ("col7", self.f7),
            ("col8", self.f8),
            ("col9", self.f9),
        ]
    }
}

/// The kind of a MACS output file, which fixes the coordinate and score
/// columns that matter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// `*_peaks.narrowPeak`: chrom start end name score strand signal p q -1
    NarrowPeak,
    /// `*_peaks.broadPeak`: chrom start end name score strand signal p q
    BroadPeak,
    /// `*_summits.bed`: chrom start end name score strand
    Summits,
    /// `*_peaks.xls`: 1-based inclusive start, columns as documented
    Xls,
    /// `*.bdg`: chrom start end value
    BedGraph,
    /// `*_peaks.gappedPeak`: chrom start end name score thickStart thickEnd itemRgb
    GappedPeak,
    /// `*_model.r`: an R script, compared as text
    Text,
}

impl Format {
    /// Classify a file name.
    pub fn of(name: &str) -> Format {
        if name.ends_with(".narrowPeak") {
            Format::NarrowPeak
        } else if name.ends_with(".broadPeak") {
            Format::BroadPeak
        } else if name.ends_with(".gappedPeak") {
            Format::GappedPeak
        } else if name.ends_with(".xls") {
            Format::Xls
        } else if name.ends_with("_summits.bed") {
            Format::Summits
        } else if name.ends_with(".bdg") {
            Format::BedGraph
        } else {
            Format::Text
        }
    }

    /// The human label used in reports.
    pub fn label(&self) -> &'static str {
        match self {
            Format::NarrowPeak => "narrowPeak",
            Format::BroadPeak => "broadPeak",
            Format::GappedPeak => "gappedPeak",
            Format::Summits => "summits.bed",
            Format::Xls => "peaks.xls",
            Format::BedGraph => "bedGraph",
            Format::Text => "text",
        }
    }

    /// The column that carries `-log10 p`, and its tolerance.
    ///
    /// `None` when the format has no p-value.
    pub fn p_column(&self) -> Option<(usize, f64)> {
        match self {
            Format::NarrowPeak | Format::BroadPeak => Some((7, 1e-9)),
            Format::Xls => Some((8, 1e-9)),
            _ => None,
        }
    }

    /// The column that carries `-log10 q`, and its tolerance.
    pub fn q_column(&self) -> Option<(usize, f64)> {
        match self {
            Format::NarrowPeak | Format::BroadPeak => Some((8, 1e-6)),
            Format::Xls => Some((9, 1e-6)),
            _ => None,
        }
    }

    /// The column that carries fold enrichment, and its tolerance.
    ///
    /// narrowPeak column 6 is the pileup above lambda, not a ratio, so the
    /// relative tolerance applies to `xls` column 7.
    pub fn fe_column(&self) -> Option<(usize, f64)> {
        match self {
            Format::Xls => Some((7, 1e-6)),
            _ => None,
        }
    }
}

/// Outcome of comparing one file pair.
#[derive(Debug, Default)]
pub struct FileReport {
    /// Number of records on the MACS3 side.
    pub macs3_count: usize,
    /// Number of records on the macs3-rs side.
    pub ours_count: usize,
    /// Records that matched on `(chrom, start, end)` at the same index.
    pub exact_coordinates: usize,
    /// Records whose summit column matched, among the coordinate matches.
    pub summit_matches: usize,
    /// Largest absolute summit difference, in bases.
    pub max_summit_shift: u64,
    /// Largest absolute difference per named numeric quantity.
    pub max_abs_diff: Vec<(&'static str, f64)>,
    /// Whether the two files are byte-identical.
    pub bytes_identical: bool,
    /// The first divergence, with genomic context.
    pub first_divergence: Option<String>,
    /// Tolerance breaches, as `quantity: got vs limit`.
    pub violations: Vec<String>,
}

impl fmt::Display for FileReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "  peak count        MACS3 {}   macs3-rs {}",
            self.macs3_count, self.ours_count
        )?;
        if self.macs3_count == self.ours_count {
            writeln!(
                f,
                "  coordinates       {} / {} exact",
                self.exact_coordinates, self.macs3_count
            )?;
        }
        if self.summit_matches > 0 {
            writeln!(
                f,
                "  summits           {} / {}, max shift {} bp",
                self.summit_matches, self.macs3_count, self.max_summit_shift
            )?;
        }
        for (name, d) in &self.max_abs_diff {
            if d.is_finite() && *d > 0.0 {
                writeln!(f, "  {name:<18} max |delta| {d:e}")?;
            }
        }
        writeln!(
            f,
            "  bytes             {}",
            if self.bytes_identical {
                "identical"
            } else {
                "differ"
            }
        )?;
        if let Some(d) = &self.first_divergence {
            writeln!(f, "  first divergence  {d}")?;
        }
        for v in &self.violations {
            writeln!(f, "  VIOLATION         {v}")?;
        }
        Ok(())
    }
}

impl FileReport {
    /// True when nothing exceeded a tolerance and the counts agree.
    pub fn is_clean(&self) -> bool {
        self.violations.is_empty() && self.macs3_count == self.ours_count
    }
}

/// Parse a whole output file into records.
pub fn parse_file(text: &str) -> Vec<Record> {
    text.lines().filter_map(Record::parse).collect()
}

/// Compare two output files of the same format.
pub fn compare(macs3_text: &str, ours_text: &str, format: Format) -> FileReport {
    let a = parse_file(macs3_text);
    let b = parse_file(ours_text);
    let mut r = FileReport {
        macs3_count: a.len(),
        ours_count: b.len(),
        bytes_identical: macs3_text == ours_text,
        max_summit_shift: 0,
        ..Default::default()
    };

    // Peak counts are the headline: a count mismatch subsumes every other
    // difference, so report it and stop rather than producing noise.
    if a.len() != b.len() {
        r.first_divergence = Some(format!(
            "record count {} != {}; coordinate comparison skipped",
            a.len(),
            b.len()
        ));
        r.violations
            .push(format!("peak count {} != {}", a.len(), b.len()));
        return r;
    }

    let mut diffs: Vec<(&'static str, f64)> = Vec::new();
    let mut push = |name: &'static str, d: f64| {
        if d.is_finite() && d > 0.0 {
            match diffs.iter_mut().find(|(n, _)| *n == name) {
                Some(slot) => {
                    if d > slot.1 {
                        slot.1 = d;
                    }
                }
                None => diffs.push((name, d)),
            }
        }
    };

    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        let same_coords = x.chrom == y.chrom && x.start == y.start && x.end == y.end;
        if same_coords {
            r.exact_coordinates += 1;
        } else if r.first_divergence.is_none() {
            r.first_divergence = Some(format!(
                "record {i}: MACS3 {}:{}-{} vs macs3-rs {}:{}-{}",
                x.chrom, x.start, x.end, y.chrom, y.start, y.end
            ));
        }
        // name / summit
        if same_coords {
            let dx =
                (x.name.parse::<u64>().unwrap_or(0)).abs_diff(y.name.parse::<u64>().unwrap_or(0));
            r.max_summit_shift = r.max_summit_shift.max(dx);
            if x.name == y.name {
                r.summit_matches += 1;
            }
        }
        for (n, xv) in x.numeric_columns() {
            let yv = y
                .numeric_columns()
                .iter()
                .find(|(m, _)| *m == n)
                .map(|(_, v)| *v)
                .unwrap_or(f64::NAN);
            push(n, (xv - yv).abs());
        }
        if let Some((col, limit)) = format.p_column() {
            let xv = numeric_at(x, col);
            let yv = numeric_at(y, col);
            let d = (xv - yv).abs();
            push("p-value", d);
            if d > limit && d.is_finite() {
                r.violations.push(format!(
                    "p-value |delta| {d:e} at {}:{}-{} exceeds {limit:e}",
                    x.chrom, x.start, x.end
                ));
            }
        }
        if let Some((col, limit)) = format.q_column() {
            let xv = numeric_at(x, col);
            let yv = numeric_at(y, col);
            let d = (xv - yv).abs();
            push("q-value", d);
            if d > limit && d.is_finite() {
                r.violations.push(format!(
                    "q-value |delta| {d:e} at {}:{}-{} exceeds {limit:e}",
                    x.chrom, x.start, x.end
                ));
            }
        }
        if let Some((col, limit)) = format.fe_column() {
            let xv = numeric_at(x, col);
            let yv = numeric_at(y, col);
            let d = (xv - yv).abs() / xv.abs().max(1e-12);
            push("fold enrichment", d);
            if d > limit && d.is_finite() {
                r.violations.push(format!(
                    "fold enrichment relative |delta| {d:e} at {}:{}-{} exceeds {limit:e}",
                    x.chrom, x.start, x.end
                ));
            }
        }
    }

    r.max_abs_diff = diffs;
    r
}

fn numeric_at(r: &Record, col: usize) -> f64 {
    r.fields
        .get(col)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_skips_comments() {
        let t = "# comment\ntrack name=x\nchr1\t10\t20\tpeak_1\t100\t.\t5.0\t3.0\t2.0\t-1\n";
        let recs = parse_file(t);
        assert_eq!(recs.len(), 1);
        let r = &recs[0];
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 10);
        assert_eq!(r.end, 20);
        assert_eq!(r.name, "peak_1");
        assert!((r.f7 - 3.0).abs() < 1e-12, "p column is index 7");
        assert!((r.f8 - 2.0).abs() < 1e-12, "q column is index 8");
    }

    #[test]
    fn missing_columns_become_nan_not_an_error() {
        let recs = parse_file("chr1\t10\t20\tsummit_1\t0\t.\n");
        assert_eq!(recs.len(), 1);
        assert!(recs[0].f7.is_nan());
    }

    #[test]
    fn identical_files_are_clean() {
        let t = "chr1\t10\t20\tp_1\t100\t.\t5.0\t3.0\t2.0\t-1\n";
        let r = compare(t, t, Format::NarrowPeak);
        assert!(r.is_clean(), "{r}");
        assert!(r.bytes_identical);
        assert_eq!(r.exact_coordinates, 1);
    }

    #[test]
    fn a_count_mismatch_is_reported_on_its_own() {
        let a = "chr1\t10\t20\tp_1\t100\t.\t5.0\t3.0\t2.0\t-1\n";
        let b = "chr1\t10\t20\tp_1\t100\t.\t5.0\t3.0\t2.0\t-1\nchr1\t40\t50\tp_2\t100\t.\t5.0\t3.0\t2.0\t-1\n";
        let r = compare(a, b, Format::NarrowPeak);
        assert!(!r.is_clean());
        assert_eq!(r.macs3_count, 1);
        assert_eq!(r.ours_count, 2);
        assert!(r.first_divergence.unwrap().contains("record count"));
    }

    #[test]
    fn a_coordinate_shift_is_located() {
        let a = "chr1\t10\t20\tp_1\t100\t.\t5.0\t3.0\t2.0\t-1\n";
        let b = "chr1\t11\t20\tp_1\t100\t.\t5.0\t3.0\t2.0\t-1\n";
        let r = compare(a, b, Format::NarrowPeak);
        assert_eq!(r.exact_coordinates, 0);
        assert!(r.first_divergence.unwrap().contains("chr1:10-20"));
    }

    #[test]
    fn a_q_drift_inside_the_tolerance_is_not_a_violation() {
        let a = "chr1\t10\t20\tp_1\t100\t.\t5.0\t3.0\t2.0\t-1\n";
        let b = "chr1\t10\t20\tp_1\t100\t.\t5.0\t3.0\t2.0000005\t-1\n";
        let r = compare(a, b, Format::NarrowPeak);
        assert!(r.is_clean(), "{r}");
        assert!(r.violations.is_empty());
    }

    #[test]
    fn a_q_drift_outside_the_tolerance_is_a_violation() {
        let a = "chr1\t10\t20\tp_1\t100\t.\t5.0\t3.0\t2.0\t-1\n";
        let b = "chr1\t10\t20\tp_1\t100\t.\t5.0\t3.0\t2.5\t-1\n";
        let r = compare(a, b, Format::NarrowPeak);
        assert!(!r.is_clean());
        assert!(r.violations.iter().any(|v| v.starts_with("q-value")));
    }

    #[test]
    fn format_classification() {
        assert_eq!(Format::of("x_peaks.narrowPeak"), Format::NarrowPeak);
        assert_eq!(Format::of("x_peaks.xls"), Format::Xls);
        assert_eq!(Format::of("x_summits.bed"), Format::Summits);
        assert_eq!(Format::of("x_treat_pileup.bdg"), Format::BedGraph);
        assert_eq!(Format::of("x_model.r"), Format::Text);
    }

    #[test]
    fn xls_uses_its_own_p_and_q_columns() {
        // chrom start end length absSummit pileup logPileup logLR pscore qscore
        let hdr = "chr1\t1\t100\t100\t50\t20\t1.3\t0.5\t4.0\t3.0\n";
        let r = compare(hdr, hdr, Format::Xls);
        assert!(r.is_clean(), "{r}");
        assert_eq!(Format::Xls.p_column().unwrap().0, 8);
        assert_eq!(Format::Xls.q_column().unwrap().0, 9);
        assert_eq!(Format::Xls.fe_column().unwrap().0, 7);
    }
}
