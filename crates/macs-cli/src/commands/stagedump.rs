//! Emit our own `stages.json`, in the same schema as `oracle/dump_stages.py`.
//!
//! This is what makes the stage-by-stage layer (L3) a real *differential* rather than
//! a golden-against-itself check: `macs-compare --stages` needs both sides, and until
//! now only the upstream side existed.
//!
//! # Why an environment variable and not a flag
//!
//! Adding a `--dump-stages` flag to `callpeak` would be a new argument, and the
//! release definition requires *identical accept/reject behaviour* on invalid
//! invocations. A flag upstream does not have is an accept/reject divergence, however
//! useful it is here. `MACS3_RS_DUMP_STAGES` is invisible to argument parsing, which
//! is the same trick `MACS3_RS_THREADS` uses.
//!
//! # What is recorded
//!
//! Mirrors `dump_stages.py`'s table. Stages upstream cannot reach without patching it
//! (the per-scale local lambda arrays) are not invented here: a stage absent on our
//! side is reported as missing, which is the honest signal.

use std::collections::BTreeMap;
use std::fmt::Write as _;

/// An accumulator for the stages, so the call sites stay one-liners.
#[derive(Debug, Default)]
pub struct StageDump {
    dir: Option<std::path::PathBuf>,
    /// Insertion-ordered so the JSON reads in pipeline order.
    stages: Vec<(String, String)>,
}

impl StageDump {
    /// Read the environment. A no-op when the variable is unset, so the hot path pays
    /// one branch per stage.
    pub fn from_env() -> Self {
        let dir = match std::env::var("MACS3_RS_DUMP_STAGES") {
            Ok(v) if !v.is_empty() => Some(std::path::PathBuf::from(v)),
            _ => None,
        };
        StageDump {
            dir,
            stages: Vec::new(),
        }
    }

    pub fn enabled(&self) -> bool {
        self.dir.is_some()
    }

    /// Record a stage as pre-rendered JSON. The caller owns the escaping, which keeps
    /// this module free of a JSON writer for the few shapes involved.
    pub fn put(&mut self, name: &str, json: &str) {
        if self.dir.is_some() {
            self.stages.push((name.to_string(), json.to_string()));
        }
    }

    /// Write `stages.json`. A write failure is reported on stderr rather than turned
    /// into a non-zero exit: the dump is a debugging aid, and failing the run because
    /// a scratch directory is read-only would be worse than a missing file.
    pub fn finish(&self, extra: &[(&str, String)]) {
        let Some(dir) = &self.dir else { return };
        if let Err(e) = std::fs::create_dir_all(dir) {
            eprintln!("macs3-rs: stage dump: cannot create {}: {e}", dir.display());
            return;
        }
        let mut body = String::from("{\n");
        // Keys must be quoted: `dump_stages.py` writes JSON, and a bare `name:` is not
        // JSON at all -- which the comparator's parser (correctly) rejects.
        for (name, json) in &self.stages {
            let _ = writeln!(body, " \"{}\": {json},", esc(name));
        }
        for (name, json) in extra {
            let _ = writeln!(body, " \"{}\": {json},", esc(name));
        }
        body.push_str(" \"_argv\": {}\n}\n");
        let path = dir.join("stages.json");
        if let Err(e) = std::fs::write(&path, body) {
            eprintln!("macs3-rs: stage dump: cannot write {}: {e}", path.display());
        }
    }
}

/// JSON for a count map: `{"total": N, "chroms": {name: {plus, minus}}}`.
pub fn counts_json(total: u64, per_chrom: &BTreeMap<String, (u64, u64)>) -> String {
    let mut s = format!("{{\"total\": {total}, \"chroms\": {{");
    for (i, (c, (plus, minus))) in per_chrom.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        let _ = write!(
            s,
            "\"{}\": {{\"plus\": {plus}, \"minus\": {minus}}}",
            esc(c)
        );
    }
    s.push_str("}}");
    s
}

/// JSON for a map whose values are already rendered, so a boolean stays a boolean.
///
/// `scaling.tocontrol` is a flag upstream records as `false`, not as `0`; comparing
/// `false` against `0.0` would report a type difference on every run and hide the
/// numbers that matter.
pub fn mixed_json(pairs: &[(&str, String)]) -> String {
    let mut s = String::from("{");
    for (i, (k, v)) in pairs.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        let _ = write!(s, "\"{}\": {v}", esc(k));
    }
    s.push('}');
    s
}

/// JSON for a string->number map.
pub fn num_map_json(pairs: &[(&str, f64)]) -> String {
    let mut s = String::from("{");
    for (i, (k, v)) in pairs.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        let _ = write!(s, "\"{}\": {}", esc(k), num(*v));
    }
    s.push('}');
    s
}

/// JSON for the control scale parameters actually passed to the peak score engine.
pub fn lambda_ladder_json(
    ctrl_d_s: &[i64],
    ctrl_scaling_factor_s: &[f64],
    lambda_bg: f64,
    treat_scale: f64,
) -> String {
    let d_s = ctrl_d_s
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    let factors = ctrl_scaling_factor_s
        .iter()
        .map(|&v| num(v))
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "{{\"ctrl_d_s\": [{d_s}], \"ctrl_scaling_factor_s\": [{factors}], \"lambda_bg\": {}, \"treat_scale\": {}}}",
        num(lambda_bg),
        num(treat_scale)
    )
}

/// `dump_stages.py` records lengths as plain integers, and `qvalue_table` carries the
/// **signed** totals, so they are emitted as-is rather than clamped.
pub fn num(v: f64) -> String {
    if v.is_finite() && v == v.trunc() && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        format!("{v}")
    }
}

/// Escape a string for a JSON string literal.
pub fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

/// Record `reads_pre_filter`, `reads_post_filter`, `d` and `xls_header`.
///
/// Called from inside each input branch of `callpeak`, where the tracks and the
/// post-filter counters are all still in scope -- they are consumed and shadowed by
/// the per-mode branches otherwise, so there is no single later point that can see
/// both.
/// Record `reads_pre_filter`, `reads_post_filter`, `d`, `xls_header` and `peak_files`.
///
/// `pre_filter` is the already-rendered `reads_pre_filter` payload, captured at load
/// time by [`pre_json_se`] / [`pre_json_frag`] -- before `filter_dup` rewrote the
/// track. Rendering it here instead would report post-filter counts under a
/// `pre_filter` name, which manufactures a divergence on every paired-end fixture.
#[allow(clippy::too_many_arguments)]
pub fn record_reads(
    sd: &mut StageDump,
    pre_filter: Option<&str>,
    ctrl_present: bool,
    tsize: f64,
    t1: u64,
    c1: u64,
    slocal: i64,
    llocal: i64,
) {
    if !sd.enabled() {
        return;
    }
    // No snapshot means the stage is simply absent -- recorded as absent rather than
    // guessed, because a missing stage is reported by the comparator and a wrong one is
    // not.
    if let Some(j) = pre_filter {
        sd.put("reads_pre_filter", j);
    }
    if ctrl_present {
        sd.put(
            "reads_post_filter",
            &format!("{{\"treatment\": {t1}, \"control\": {c1}}}"),
        );
    }
    sd.put("d", &num(tsize));
    sd.put(
        "xls_header",
        &num_map_json(&[
            ("d", tsize),
            ("slocal", slocal as f64),
            ("llocal", llocal as f64),
        ]),
    );
    // `peak_files` names the run's own outputs, relative to the dump directory. An
    // absolute path would differ between the golden and "ours" roots by construction.
    sd.put(
        "peak_files",
        &format!(
            "{{\"candidate\": [\"{}\", \"{}\"], \"final\": [\"{}\", \"{}\"]}}",
            esc("candidate_peaks.peaks.xls"),
            esc("candidate_peaks.peaks.narrowPeak"),
            esc("final_peaks.peaks.xls"),
            esc("final_peaks.peaks.narrowPeak")
        ),
    );
}

/// The single-end `reads_pre_filter` payload, rendered at the moment the track is
/// loaded -- i.e. **before** `filter_dup` rewrites it.
/// Only the treatment half is available at `pool_se` time: the control is pooled a
/// few lines later, and reordering the pipeline to capture both would mean touching
/// the parity-critical path for a debugging aid. The `control` sub-object is therefore
/// absent, and `macs-compare` reports those leaves as missing on our side rather than
/// silently comparing nothing.
pub fn pre_json_se(treat: &macs_track::SingleEndTrack) -> String {
    let tm = per_chrom_se(treat);
    format!("{{\"treatment\": {}}}", counts_json(treat.total(), &tm))
}

/// Append the control snapshot captured before its duplicate filter to an earlier
/// treatment snapshot. `callpeak` filters treatment before it loads control, so the
/// two pre-filter snapshots necessarily arrive at separate points in the pipeline.
pub fn add_pre_control_se(
    pre_treatment: &str,
    control: Option<&macs_track::SingleEndTrack>,
) -> String {
    let treatment = pre_treatment
        .strip_prefix("{\"treatment\": ")
        .and_then(|v| v.strip_suffix('}'))
        .expect("pre_json_se has a stable treatment wrapper");
    let control = control.map_or_else(
        || "null".to_string(),
        |track| counts_json(track.total(), &per_chrom_se(track)),
    );
    format!("{{\"treatment\": {treatment}, \"control\": {control}}}")
}

/// `reads_pre_filter` payload for paired fragments, as read by the oracle wrapper.
pub fn pre_json_frag(track: &macs_track::FragmentTrack, role: &str) -> String {
    let mut chroms = BTreeMap::new();
    for chrom in track.chroms() {
        let name = String::from_utf8_lossy(track.genome().name(chrom)).into_owned();
        let fields = if track.has_counts() {
            "[\"l\",\"r\",\"c\"]"
        } else {
            "[\"l\",\"r\"]"
        };
        chroms.insert(
            name,
            format!(
                "{{\"fragments\": {}, \"fields\": {fields}}}",
                track.frags(chrom).len()
            ),
        );
    }
    let mut body = String::from("{");
    for (i, (name, details)) in chroms.iter().enumerate() {
        if i > 0 {
            body.push(',');
        }
        let _ = write!(body, "\"{}\": {details}", esc(name));
    }
    body.push('}');
    format!(
        "{{\"{role}\": {{\"total\": {}, \"chroms\": {body}}}}}",
        track.total()
    )
}

/// Add the second paired-end input's pre-filter snapshot to the first.
pub fn add_pre_control_frag(
    pre_treatment: &str,
    control: Option<&macs_track::FragmentTrack>,
) -> String {
    let treatment = pre_treatment
        .strip_prefix("{\"treatment\": ")
        .and_then(|v| v.strip_suffix('}'))
        .expect("pre_json_frag has a stable treatment wrapper");
    let control = control.map_or_else(
        || "null".to_string(),
        |track| {
            let wrapped = pre_json_frag(track, "control");
            wrapped
                .strip_prefix("{\"control\": ")
                .and_then(|v| v.strip_suffix('}'))
                .expect("pre_json_frag has a stable control wrapper")
                .to_string()
        },
    );
    format!("{{\"treatment\": {treatment}, \"control\": {control}}}")
}

/// Duplicate counts/rates reconstructed from the actual input and filtered totals.
/// FRAG tracks carry counts and upstream skips duplicate filtering, so their capture
/// intentionally contains totals only.
pub fn record_duplicates(
    sd: &mut StageDump,
    treat_total: u64,
    treat_after: u64,
    ctrl_total: u64,
    ctrl_after: u64,
    ctrl_present: bool,
    counted_fragments: bool,
) {
    if !sd.enabled() {
        return;
    }
    let mut fields = vec![("treat_total", treat_total as f64)];
    if counted_fragments {
        if ctrl_present {
            fields.push(("ctrl_total", ctrl_total as f64));
        }
    } else {
        fields.push(("treat_after_filter", treat_after as f64));
        fields.push((
            "treat_redundant_rate",
            redundant_rate(treat_total, treat_after),
        ));
        if ctrl_present {
            fields.push(("ctrl_total", ctrl_total as f64));
            fields.push(("ctrl_after_filter", ctrl_after as f64));
            fields.push((
                "ctrl_redundant_rate",
                redundant_rate(ctrl_total, ctrl_after),
            ));
        }
    }
    sd.put("duplicates", &num_map_json(&fields));
}

fn redundant_rate(total: u64, kept: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        ((total.saturating_sub(kept) as f64 / total as f64) * 100.0).round() / 100.0
    }
}

/// The cutoff-analysis report is the upstream stage named `qvalue_table`; it has
/// five fields per emitted cutoff. Parse the same text that the output writer uses so
/// rounding matches the user's report exactly.
pub fn cutoff_table_json(body: &str) -> String {
    let mut out = String::from("[");
    let mut first = true;
    for line in body.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let (Ok(p), Ok(q), Ok(np), Ok(lp), Ok(avg)) = (
            fields[0].parse::<f64>(),
            fields[1].parse::<f64>(),
            fields[2].parse::<u64>(),
            fields[3].parse::<u64>(),
            fields[4].parse::<f64>(),
        ) else {
            continue;
        };
        if !first {
            out.push(',');
        }
        first = false;
        let _ = write!(out, "[{}, {}, {np}, {lp}, {}]", num(p), num(q), num(avg));
    }
    out.push(']');
    out
}

/// Parse `-B`'s actual control-lambda body into the same interval schema used by the
/// oracle's bedGraph reader.
pub fn bedgraph_body_json(body: &str) -> String {
    let mut chroms: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in body.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 4 {
            continue;
        }
        let (Ok(start), Ok(end), Ok(value)) = (
            fields[1].parse::<i64>(),
            fields[2].parse::<i64>(),
            fields[3].parse::<f64>(),
        ) else {
            continue;
        };
        chroms
            .entry(fields[0].to_string())
            .or_default()
            .push(format!("[{start},{end},{}]", num(value)));
    }
    let mut out = String::from("{");
    for (i, (chrom, rows)) in chroms.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(out, "\"{}\": [{}]", esc(chrom), rows.join(","));
    }
    out.push('}');
    out
}

/// Per-chromosome plus/minus read counts, read back through the track API rather than
/// from the counters, so a bug in either path shows up as a disagreement.
pub fn per_chrom_se(t: &macs_track::SingleEndTrack) -> BTreeMap<String, (u64, u64)> {
    let mut m = BTreeMap::new();
    for cid in t.positions().chroms() {
        let name = String::from_utf8_lossy(t.genome().name(cid)).into_owned();
        let p = t.positions().strand(cid, macs_core::Strand::Plus).len() as u64;
        let q = t.positions().strand(cid, macs_core::Strand::Minus).len() as u64;
        m.insert(name, (p, q));
    }
    m
}

/// Record `treat_pileup` and `ctrl_pileup`, matching the shape `dump_stages.py`
/// produces from the engine's own `-B` bedGraphs: `{chrom: [[start, end, value], ...]}`.
///
/// Intervals are kept, not just values. A differential that compared only the values
/// would miss a segmentation difference -- exactly the F207 class of bug, where the
/// depths were right but the runs were split differently.
///
/// Values are rendered at `%.5f` because that is what `bedGraphIO.write_bedGraph`
/// writes, so the two sides compare the digits a user would actually see rather than
/// differing at the 17th significant figure on a value that gets rounded away.
pub fn record_tracks(
    sd: &mut StageDump,
    treat_pileup: &[(String, &macs_rle::SignalTrack<f32>)],
    control: &[(String, &macs_rle::SignalTrack<f32>)],
    coord_shift: i64,
) {
    if !sd.enabled() {
        return;
    }
    sd.put("treat_pileup", &tracks_json(treat_pileup, coord_shift));
    // `control` is only emitted when the oracle can produce a counterpart. It has no
    // `-B` equivalent -- upstream writes the *control lambda*, not the control pileup --
    // so emitting it here would register as "present only in ours" forever, which is
    // noise that trains people to ignore the comparator.
    if !control.is_empty() {
        sd.put("ctrl_pileup", &tracks_json(control, coord_shift));
    }
}

/// `{name: [[start, end, value], ...]}` over a per-chromosome track list.
///
/// `coord_shift` is subtracted because a counted (`--format FRAG`) run is computed in a
/// shifted frame -- `coord_shift = max(d, slocal, llocal) / 2` (F181) -- and every writer
/// subtracts it again before emitting a coordinate. Dumping the raw shifted track would
/// report the whole FRAG pileup displaced by 5000 bases and look like a real divergence.
fn tracks_json(tracks: &[(String, &macs_rle::SignalTrack<f32>)], coord_shift: i64) -> String {
    let mut out = String::from("{");
    for (i, (name, t)) in tracks.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&format!("\"{}\": [", esc(name)));
        // `Run` stores only its end; the start of run `i` is the previous end, and the
        // first run starts where the track does.
        let mut sep = "";
        let mut start = t.start();
        for r in t.runs() {
            out.push_str(sep);
            sep = ",";
            out.push_str(&format!(
                "[{},{},{:.5}]",
                (start as i64 - coord_shift).max(0),
                (r.end as i64 - coord_shift).max(0),
                r.value
            ));
            start = r.end;
        }
        out.push(']');
    }
    out.push('}');
    out
}
