//! Per-stage comparison of two `dump_stages.py` output trees.
//!
//! `oracle/dump_stages.py` drives upstream's real pipeline in-process and records
//! what each stage computed. `macs-compare`'s file comparator answers "do the final
//! outputs agree?"; this answers the earlier question the porting plan actually
//! cares about -- **which stage first diverged** -- by comparing the intermediates
//! directly.
//!
//! # Why the stages are typed rather than all compared as numbers
//!
//! A stage dump mixes three kinds of leaf, and averaging them into one deviation
//! number would hide exactly the bugs this layer exists to catch:
//!
//! * **counts and coordinates** (`reads_pre_filter`, `d`, `candidate_peaks`) must
//!   match *exactly*. A read-count difference of one is a real defect, not a rounding
//!   difference, so it is reported as a mismatch rather than as a small deviation.
//! * **scores** (`scaling`, and the numeric leaves inside the tracks) are compared
//!   with a tolerance, and the worst one is reported per stage.
//! * **file references** (`peak_files`, `candidate_peaks.peaks.xls`) name a sibling
//!   file; those are handed back to the caller's file comparator rather than being
//!   string-compared.
//!
//! # Order matters
//!
//! Stages are reported in pipeline order, not alphabetically, so the first `FAIL` is
//! the earliest stage that went wrong. `STAGE_ORDER` is that order; anything not
//! listed sorts after it, which keeps a newly added stage visible instead of dropped.
//!
//! # Paths are not stage data
//!
//! A stage dump records the command line it was driven with, and a command line names
//! the checkout. The corpus is committed with that path rewritten to the neutral token
//! `<ROOT>` (`oracle/relocate_golden.py`), so a replay -- which writes the live path --
//! and the recording agree on everything that is compared. Exactly two places in a
//! dump are allowed to name a path, and both are outside the comparison: `_argv`,
//! which is a `_`-prefixed top-level key that `compare_trees` skips, and the leaves
//! under [`FILE_KEYS`], which name a sibling file that is already known to live under a
//! different root on each side. `recorded_leaves_name_no_path` in the test suite
//! asserts that for the committed corpus, because a path leaking into a compared leaf
//! would make L3 fail for a reason that has nothing to do with MACS3.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Pipeline order. `dump_stages.py`'s own table is the reference.
pub const STAGE_ORDER: &[&str] = &[
    "reads_pre_filter",
    "reads_post_filter",
    // retained counts and redundant rates, parsed from upstream's own log
    "duplicates",
    "totals",
    "d",
    "scaling",
    // the d/slocal/llocal control ladder. Captured by wrapping `call_peaks`, which
    // records that the ladder is a C attribute of `CallerFromAlignments` and is
    // **not** passed as an argument -- the wrapper returns an empty ladder for every
    // run, which is itself the evidence that there is no Python-visible route.
    "lambda_ladder",
    "xls_header",
    // `dump_stages.py` writes these two from the `-B` bedGraphs the engine itself
    // produces; the control pileup and the p/q tracks are cdf attributes of
    // `CallerFromAlignments` and are not reachable from Python, so they are listed
    // only so that a future capture shows up in the right place rather than at the end.
    "treat_pileup",
    "lambda_merged",
    "ctrl_pileup",
    "pvalue_track",
    "qvalue_track",
    "qvalue_table",
    "peak_files",
];

/// Keys that are never compared as values -- they name a sibling file.
const FILE_KEYS: &[&str] = &[
    "peak_files",
    "candidate_peaks.peaks.xls",
    "candidate_peaks.summits.bed",
    "candidate_peaks.peaks.narrowPeak",
    "final_peaks.peaks.xls",
    "final_peaks.summits.bed",
];

/// The neutral token the recorded corpus spells this checkout with.
///
/// `oracle/relocate_golden.py` rewrites the recording's absolute paths to it, and
/// `oracle/run_golden.py` expands it back per checkout for the golden corpus. Nothing
/// here expands or compares it -- it is named so a leaf that has leaked one is
/// recognisable, which is what [`compared_strings`] is for.
pub const ROOT_TOKEN: &str = "<ROOT>";

/// One stage's verdict, aggregated across every recorded run.
#[derive(Debug, Clone, PartialEq)]
pub struct StageReport {
    pub name: String,
    /// Runs that contributed a report for this stage.
    pub runs: usize,
    /// Leaves compared (numbers and exact-match scalars).
    pub numeric: usize,
    /// Leaves that are integers or strings and had to match exactly.
    pub exact: usize,
    /// Worst absolute difference over the numeric leaves.
    pub max_abs: f64,
    /// Worst relative difference; `inf` when the reference is zero and the value is
    /// not, which is reported as-is rather than silently normalised.
    pub max_rel: f64,
    /// A strict or tolerant mismatch that the tolerance did not cover.
    pub mismatches: Vec<String>,
    /// True when the stage is absent from one side.
    pub missing: Option<String>,
    /// Cached position in [`STAGE_ORDER`], so aggregation does not re-scan it per run.
    rank_cache: usize,
}

impl StageReport {
    fn new(name: &str) -> Self {
        StageReport {
            name: name.to_string(),
            runs: 0,
            numeric: 0,
            exact: 0,
            max_abs: 0.0,
            max_rel: 0.0,
            mismatches: Vec::new(),
            missing: None,
            rank_cache: STAGE_ORDER
                .iter()
                .position(|s| *s == name)
                .unwrap_or(usize::MAX),
        }
    }

    /// Clean when nothing mismatched and the stage was present on both sides.
    pub fn is_clean(&self) -> bool {
        self.mismatches.is_empty() && self.missing.is_none()
    }

    fn rank(&self) -> usize {
        self.rank_cache
    }
}

/// A minimal JSON reader.
///
/// Deliberately not a dependency: this tool ships inside the port, and the stage dump
/// is a small, machine-written subset of JSON (objects, arrays, strings, numbers,
/// booleans, null). `serde_json` would be simpler, but adding it here would pull a
/// parser into the shipped binary for a dev-facing mode.
#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(BTreeMap<String, Json>),
}

impl Json {
    fn walk(&self, path: &str, out: &mut BTreeMap<String, Leaf>) {
        match self {
            Json::Obj(m) => {
                for (k, v) in m {
                    let p = if path.is_empty() {
                        k.clone()
                    } else {
                        format!("{path}.{k}")
                    };
                    v.walk(&p, out);
                }
            }
            Json::Arr(a) => {
                for (i, v) in a.iter().enumerate() {
                    v.walk(&format!("{path}[{i}]"), out);
                }
            }
            other => {
                out.insert(
                    path.to_string(),
                    match other {
                        Json::Num(v) => Leaf::Num(*v),
                        Json::Bool(b) => Leaf::Exact(if *b { "true" } else { "false" }.into()),
                        Json::Str(s) => Leaf::Exact(s.clone()),
                        Json::Null => Leaf::Exact("null".into()),
                        Json::Arr(_) | Json::Obj(_) => unreachable!(),
                    },
                );
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
enum Leaf {
    Num(f64),
    Exact(String),
}

/// Parse a stage dump.
pub fn parse(text: &str) -> Result<Json, String> {
    let b = text.as_bytes();
    let mut i = 0usize;
    let v = parse_value(b, &mut i)?;
    skip_ws(b, &mut i);
    if i != b.len() {
        return Err(format!("trailing input at byte {i}"));
    }
    Ok(v)
}

fn skip_ws(b: &[u8], i: &mut usize) {
    while *i < b.len() && matches!(b[*i], b' ' | b'\t' | b'\n' | b'\r') {
        *i += 1;
    }
}

fn parse_value(b: &[u8], i: &mut usize) -> Result<Json, String> {
    skip_ws(b, i);
    match b.get(*i) {
        None => Err("unexpected end of input".into()),
        Some(b'{') => {
            *i += 1;
            let mut m = BTreeMap::new();
            skip_ws(b, i);
            if b.get(*i) == Some(&b'}') {
                *i += 1;
                return Ok(Json::Obj(m));
            }
            loop {
                skip_ws(b, i);
                let k = match parse_value(b, i)? {
                    Json::Str(s) => s,
                    other => return Err(format!("object key must be a string, got {other:?}")),
                };
                skip_ws(b, i);
                if b.get(*i) != Some(&b':') {
                    return Err(format!("expected ':' after key {k:?}"));
                }
                *i += 1;
                let v = parse_value(b, i)?;
                m.insert(k, v);
                skip_ws(b, i);
                match b.get(*i) {
                    Some(b',') => *i += 1,
                    Some(b'}') => {
                        *i += 1;
                        return Ok(Json::Obj(m));
                    }
                    _ => return Err(format!("expected ',' or '}}' in object, at byte {i}")),
                }
            }
        }
        Some(b'[') => {
            *i += 1;
            let mut a = Vec::new();
            skip_ws(b, i);
            if b.get(*i) == Some(&b']') {
                *i += 1;
                return Ok(Json::Arr(a));
            }
            loop {
                a.push(parse_value(b, i)?);
                skip_ws(b, i);
                match b.get(*i) {
                    Some(b',') => *i += 1,
                    Some(b']') => {
                        *i += 1;
                        return Ok(Json::Arr(a));
                    }
                    _ => return Err(format!("expected ',' or ']' in array, at byte {i}")),
                }
            }
        }
        Some(b'"') => {
            *i += 1;
            let mut s = String::new();
            loop {
                let c = *b.get(*i).ok_or("unterminated string")?;
                *i += 1;
                match c {
                    b'"' => return Ok(Json::Str(s)),
                    b'\\' => {
                        let e = *b.get(*i).ok_or("unterminated escape")?;
                        *i += 1;
                        match e {
                            b'n' => s.push('\n'),
                            b't' => s.push('\t'),
                            b'r' => s.push('\r'),
                            b'b' => s.push('\u{8}'),
                            b'f' => s.push('\u{c}'),
                            b'u' => {
                                let hex = std::str::from_utf8(
                                    b.get(*i..*i + 4).ok_or("short \\u escape")?,
                                )
                                .map_err(|e| e.to_string())?;
                                *i += 4;
                                let cp = u32::from_str_radix(hex, 16)
                                    .map_err(|e| format!("bad \\u escape: {e}"))?;
                                s.push(char::from_u32(cp).unwrap_or('\u{fffd}'));
                            }
                            other => s.push(other as char),
                        }
                    }
                    _ => {
                        // Collect the raw UTF-8 byte run verbatim.
                        let start = *i - 1;
                        let mut end = start + 1;
                        while end < b.len() && b[end] != b'"' && b[end] != b'\\' {
                            end += 1;
                        }
                        s.push_str(std::str::from_utf8(&b[start..end]).map_err(|e| e.to_string())?);
                        *i = end;
                    }
                }
            }
        }
        Some(b't') if b.get(*i..*i + 4) == Some(b"true") => {
            *i += 4;
            Ok(Json::Bool(true))
        }
        Some(b'f') if b.get(*i..*i + 5) == Some(b"false") => {
            *i += 5;
            Ok(Json::Bool(false))
        }
        Some(b'n') if b.get(*i..*i + 4) == Some(b"null") => {
            *i += 4;
            Ok(Json::Null)
        }
        Some(_) => {
            let start = *i;
            while *i < b.len() && matches!(b[*i], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                *i += 1;
            }
            std::str::from_utf8(&b[start..*i])
                .map_err(|e| e.to_string())?
                .parse::<f64>()
                .map(Json::Num)
                .map_err(|_| format!("bad number at byte {start}"))
        }
    }
}

/// Compare one stage.
pub fn compare_stage(name: &str, golden: &Json, ours: &Json, tol: f64) -> StageReport {
    let mut g = BTreeMap::new();
    let mut o = BTreeMap::new();
    walk_stage(name, golden, &mut g);
    walk_stage(name, ours, &mut o);

    let mut r = StageReport::new(name);
    let mut keys: Vec<&String> = g.keys().chain(o.keys()).collect();
    keys.sort();
    keys.dedup();
    for k in keys {
        if is_file_key(k) {
            continue;
        }
        match (g.get(k), o.get(k)) {
            (None, _) | (_, None) => {
                let side = if g.contains_key(k) { "ours" } else { "golden" };
                r.mismatches.push(format!("{k}: present only in {side}"));
            }
            (Some(Leaf::Num(a)), Some(Leaf::Num(b))) => {
                r.numeric += 1;
                let d = (a - b).abs();
                if d > r.max_abs {
                    r.max_abs = d;
                }
                let rel = if *a == 0.0 {
                    if *b == 0.0 {
                        0.0
                    } else {
                        f64::INFINITY
                    }
                } else {
                    d / a.abs()
                };
                if rel > r.max_rel || r.max_rel.is_infinite() && !rel.is_infinite() {
                    r.max_rel = rel;
                }
                if d > tol {
                    r.mismatches
                        .push(format!("{k}: {a} vs {b} (abs {d:e} > tol {tol:e})"));
                }
            }
            (Some(Leaf::Exact(a)), Some(Leaf::Exact(b))) => {
                r.exact += 1;
                if a != b {
                    r.mismatches
                        .push(format!("{k}: {a:?} vs {b:?} (must match exactly)"));
                }
            }
            (Some(x), Some(y)) => r
                .mismatches
                .push(format!("{k}: type differs ({x:?} vs {y:?})")),
        }
    }
    r
}

/// A stage dump's top level maps stage name -> payload, so leaves are keyed
/// `<stage>.<...>`; `walk_stage` only descends into the stage's own payload.
fn walk_stage(stage: &str, v: &Json, out: &mut BTreeMap<String, Leaf>) {
    let mut tmp = BTreeMap::new();
    v.walk("", &mut tmp);
    for (k, leaf) in tmp {
        // A scalar stage (`d`) has an empty inner path, so the naive join would print
        // `d.` -- which reads like a truncated field name rather than the whole value.
        let key = if k.is_empty() {
            stage.to_string()
        } else {
            format!("{stage}.{k}")
        };
        out.insert(key, leaf);
    }
}

/// Does this leaf sit under a key that names a sibling file?
///
/// Matched by *prefix*, not suffix: a leaf under `peak_files` can be
/// `peak_files.candidate[0]` or `peak_files.final_peaks.peaks.xls`, neither of which
/// ends with the key name. Comparing the recorded path would report a failure on
/// every single run, because the golden tree and the "ours" tree live under different
/// roots and so record different absolute paths.
fn is_file_key(k: &str) -> bool {
    FILE_KEYS.iter().any(|f| {
        k == *f
            || k.strip_prefix(f)
                .is_some_and(|rest| rest.starts_with('.') || rest.starts_with('['))
    })
}

/// Compare two stage trees, in pipeline order.
#[derive(Debug, Default)]
pub struct TreeReport {
    pub stages: Vec<StageReport>,
    /// `stages.json` present in the golden tree but not ours, or vice versa.
    pub run_notes: Vec<String>,
}

impl TreeReport {
    pub fn is_clean(&self) -> bool {
        self.stages.iter().all(StageReport::is_clean) && self.run_notes.is_empty()
    }

    pub fn first_failure(&self) -> Option<&StageReport> {
        let mut it = self.stages.iter().filter(|s| !s.is_clean());
        let first = it.next()?;
        it.for_each(|_| {});
        self.stages
            .iter()
            .find(|s| !s.is_clean() && s.rank() <= first.rank())
    }
}

/// Find every `stages.json` under `root`, as paths relative to it.
pub fn stage_dirs(root: &Path) -> Result<Vec<PathBuf>, String> {
    let mut out = Vec::new();
    if !root.exists() {
        return Err(format!("{}: no such directory", root.display()));
    }
    let mut stack = vec![root.to_path_buf()];
    while let Some(d) = stack.pop() {
        let rd = std::fs::read_dir(&d).map_err(|e| format!("{}: {e}", d.display()))?;
        for e in rd {
            let e = e.map_err(|x| x.to_string())?;
            let p = e.path();
            if p.is_dir() {
                stack.push(p);
            } else if p.file_name().and_then(|n| n.to_str()) == Some("stages.json") {
                out.push(p.parent().unwrap_or(root).to_path_buf());
            }
        }
    }
    out.sort();
    Ok(out)
}

/// Compare two stage trees.
pub fn compare_trees(golden: &Path, ours: &Path, tol: f64) -> Result<TreeReport, String> {
    let mut rep = TreeReport::default();
    // Match runs by path **relative to each root**. The two trees live under different
    // directories by construction, so comparing absolute paths would find no pairs.
    let gd = stage_dirs(golden)?;
    let od = stage_dirs(ours)?;
    let mut gmap: BTreeMap<String, PathBuf> = BTreeMap::new();
    for d in &gd {
        gmap.insert(rel(golden, d), d.clone());
    }
    let mut omap: BTreeMap<String, PathBuf> = BTreeMap::new();
    for d in &od {
        omap.insert(rel(ours, d), d.clone());
    }
    for k in gmap.keys() {
        if !omap.contains_key(k) {
            rep.run_notes
                .push(format!("{k}: stage dump missing from ours"));
        }
    }
    for k in omap.keys() {
        if !gmap.contains_key(k) {
            rep.run_notes
                .push(format!("{k}: stage dump not in the golden tree"));
        }
    }
    for (k, d) in &gmap {
        let Some(od) = omap.get(k) else { continue };
        let gtext = std::fs::read_to_string(d.join("stages.json"))
            .map_err(|e| format!("{k}/stages.json: {e}"))?;
        let otext = std::fs::read_to_string(od.join("stages.json"))
            .map_err(|e| format!("{k}/stages.json: {e}"))?;
        let (g, o) = (parse(&gtext)?, parse(&otext)?);
        let mut names: Vec<String> = match (&g, &o) {
            (Json::Obj(a), Json::Obj(b)) => {
                let mut v: Vec<String> = a.keys().cloned().chain(b.keys().cloned()).collect();
                v.sort();
                v.dedup();
                v
            }
            _ => return Err(format!("{k}/stages.json: top level must be an object")),
        };
        names.sort_by_key(|n| {
            (
                STAGE_ORDER
                    .iter()
                    .position(|s| s == n)
                    .unwrap_or(usize::MAX),
                n.clone(),
            )
        });
        for n in names {
            if n.starts_with('_') {
                continue; // `_argv` and friends
            }
            match (gv(&g, &n), gv(&o, &n)) {
                (Some(a), Some(b)) => rep.stages.push(compare_stage(&n, &a, &b, tol)),
                (_a, _b) => {
                    let mut r = StageReport::new(&n);
                    r.missing = Some(if _a.is_some() { "ours" } else { "golden" }.into());
                    rep.stages.push(r);
                }
            }
        }
    }
    // Aggregate across runs: one row per stage, not one per (stage, run). A stage is a
    // property of the pipeline, and a per-run table on a 400-fixture corpus would be
    // 400x too long to read.
    let mut agg: BTreeMap<String, StageReport> = BTreeMap::new();
    for r in std::mem::take(&mut rep.stages) {
        let e = agg
            .entry(r.name.clone())
            .or_insert_with(|| StageReport::new(&r.name));
        e.runs += 1;
        e.numeric += r.numeric;
        e.exact += r.exact;
        if r.max_abs > e.max_abs {
            e.max_abs = r.max_abs;
        }
        if r.max_rel > e.max_rel || (e.max_rel.is_finite() && r.max_rel.is_infinite()) {
            e.max_rel = r.max_rel;
        }
        for m in r.mismatches {
            if e.mismatches.len() < 64 {
                e.mismatches.push(m);
            }
        }
        if e.missing.is_none() {
            e.missing = r.missing;
        }
    }
    rep.stages = agg.into_values().collect();
    rep.stages.sort_by_key(|s| (s.rank(), s.name.clone()));
    Ok(rep)
}

/// The exact-match leaves this comparator would compare, as `(run:key, value)`.
///
/// A diagnostic, not a verdict: everything [`compare_trees`] reaches and [`is_file_key`]
/// does not skip. It exists because a leaf that names a checkout path is a *corpus*
/// defect, not a pipeline divergence, and "L3 fails on `treat_file`" is otherwise
/// indistinguishable from "L3 found a real difference".
pub fn compared_strings(root: &Path) -> Result<Vec<(String, String)>, String> {
    let mut out = Vec::new();
    for dir in stage_dirs(root)? {
        let text = std::fs::read_to_string(dir.join("stages.json"))
            .map_err(|e| format!("{}: {e}", dir.display()))?;
        let Json::Obj(map) = parse(&text)? else {
            return Err(format!(
                "{}/stages.json: top level must be an object",
                dir.display()
            ));
        };
        for (name, payload) in map {
            if name.starts_with('_') {
                continue; // `_argv` and friends: the command line, not stage data
            }
            let mut leaves = BTreeMap::new();
            walk_stage(&name, &payload, &mut leaves);
            for (key, leaf) in leaves {
                if is_file_key(&key) {
                    continue;
                }
                if let Leaf::Exact(v) = leaf {
                    out.push((format!("{}:{key}", rel(root, &dir)), v));
                }
            }
        }
    }
    Ok(out)
}

fn gv(j: &Json, key: &str) -> Option<Json> {
    match j {
        Json::Obj(m) => m.get(key).cloned(),
        _ => None,
    }
}

fn rel(root: &Path, p: &Path) -> String {
    p.strip_prefix(root).unwrap_or(p).display().to_string()
}
