//! `macs-compare` — the genomic-aware differential comparator.
//!
//! Two modes:
//!
//! **File mode** walks a tree of MACS3 golden outputs and a matching tree of
//! `macs3-rs` outputs, matches runs by their relative path, and reports per file:
//! peak counts, exact coordinate matches, summit agreement, the largest deviation in
//! each numeric column, and byte identity. It also compares the *run records* —
//! the exit code and the set of files produced — so a divergence in which files exist
//! at all is caught, not just their contents.
//!
//! ```text
//! macs-compare --golden tests/golden --ours /tmp/ours
//! macs-compare --golden tests/golden --ours /tmp/ours --verbose
//! macs-compare --golden tests/golden --ours /tmp/ours --only se_shapes
//! macs-compare --golden tests/golden --ours /tmp/ours --format narrowPeak
//! ```
//!
//! **Stage mode** compares the intermediate `stages.json` trees written by
//! `oracle/dump_stages.py`, answering "which stage first diverged?" rather than only
//! "do the final files agree?".
//!
//! ```text
//! macs-compare --stages tests/stages --ours-stages /tmp/our_stages
//! ```
//!
//! Exit codes: `0` everything within tolerance, `1` a tolerance breach or a
//! structural difference, `2` a usage or IO problem.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use macs_compare::compare::{compare, Format};
use macs_compare::stages;

#[derive(Debug, Default)]
struct Args {
    golden: PathBuf,
    ours: PathBuf,
    stages_golden: PathBuf,
    stages_ours: PathBuf,
    only: Option<String>,
    format: Option<String>,
    verbose: bool,
    max_report: usize,
    /// Absolute tolerance for the numeric stage leaves.
    stage_tol: f64,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        max_report: 20,
        // The p-score tolerance from the release definition is 1e-9. The stage
        // comparator reports its own max deviations, so this is only the pass/fail
        // line.
        stage_tol: 1e-9,
        ..Default::default()
    };
    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--golden" => a.golden = PathBuf::from(it.next().ok_or("--golden needs a path")?),
            "--ours" => a.ours = PathBuf::from(it.next().ok_or("--ours needs a path")?),
            "--stages" => {
                a.stages_golden = PathBuf::from(it.next().ok_or("--stages needs a path")?)
            }
            "--ours-stages" => {
                a.stages_ours = PathBuf::from(it.next().ok_or("--ours-stages needs a path")?)
            }
            "--stage-tol" => {
                a.stage_tol = it
                    .next()
                    .ok_or("--stage-tol needs a number")?
                    .parse()
                    .map_err(|_| "--stage-tol needs a number")?
            }
            "--only" => a.only = it.next(),
            "--format" => a.format = it.next(),
            "--verbose" | "-v" => a.verbose = true,
            "--max-report" => {
                a.max_report = it
                    .next()
                    .ok_or("--max-report needs a number")?
                    .parse()
                    .map_err(|_| "--max-report needs a number")?
            }
            "-h" | "--help" => {
                println!(
                    "macs-compare --golden <dir> --ours <dir> [--only F] [--format F] [-v] [--max-report N]\n\
                     \x20         --stages <dir> --ours-stages <dir> [--stage-tol X]"
                );
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument {other:?}")),
        }
    }
    if !a.stages_golden.as_os_str().is_empty() {
        if a.stages_ours.as_os_str().is_empty() {
            return Err("--stages requires --ours-stages too".into());
        }
        return Ok(a);
    }
    if a.golden.as_os_str().is_empty() || a.ours.as_os_str().is_empty() {
        return Err("--golden/--ours or --stages/--ours-stages are required".into());
    }
    Ok(a)
}

/// One run directory: the set of files MACS3 produced, and the exit code.
#[derive(Debug, Default)]
struct Run {
    files: BTreeMap<String, String>,
    returncode: Option<i64>,
}

fn read_run(dir: &Path) -> Result<Run, String> {
    let mut run = Run::default();
    let entries = std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    for e in entries {
        let e = e.map_err(|x| x.to_string())?;
        let p = e.path();
        if !p.is_file() {
            continue;
        }
        let name = match p.file_name().and_then(|n| n.to_str()) {
            Some(n) => n.to_string(),
            None => continue,
        };
        if name == "command.json" {
            // The record of the invocation. Its exit status is what "did this run
            // fail" means, and it is not one of the outputs to diff.
            if let Ok(text) = std::fs::read_to_string(&p) {
                run.returncode = extract_returncode(&text);
            }
            continue;
        }
        let digest = sha256_file(&p)?;
        run.files.insert(name, digest);
    }
    Ok(run)
}

/// Pull `returncode` out of a `command.json` without a JSON dependency.
fn extract_returncode(text: &str) -> Option<i64> {
    let key = "\"returncode\"";
    let i = text.find(key)? + key.len();
    let rest = text[i..].trim_start();
    let rest = rest.strip_prefix(':')?.trim_start();
    let end = rest
        .find(|c: char| !c.is_ascii_digit() && c != '-')
        .unwrap_or(rest.len());
    rest[..end].parse().ok()
}

fn sha256_file(p: &Path) -> Result<String, String> {
    use std::io::Read;
    let mut f = std::fs::File::open(p).map_err(|e| format!("{}: {e}", p.display()))?;
    let mut h = sha256::Sha256::new();
    let mut buf = vec![0u8; 1 << 16];
    loop {
        let n = f
            .read(&mut buf)
            .map_err(|e| format!("{}: {e}", p.display()))?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(h.hex())
}

/// A dependency-free SHA-256, so the comparator adds no crates to the shipped binary.
///
/// It has to be here rather than in `compare`: a byte-identity verdict that trusted a
/// cheaper hash would be the one thing in this tool nobody re-checks.
mod sha256 {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    const H0: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];

    #[derive(Clone)]
    pub struct Sha256 {
        h: [u32; 8],
        /// Bytes not yet compressed; always < 64 after `update`.
        buf: Vec<u8>,
        /// Total message length in bytes, needed for the length field in the padding.
        len: u64,
    }

    impl Sha256 {
        pub fn new() -> Self {
            Sha256 {
                h: H0,
                buf: Vec::with_capacity(64),
                len: 0,
            }
        }

        pub fn update(&mut self, data: &[u8]) {
            self.len += data.len() as u64;
            self.buf.extend_from_slice(data);
            let whole = self.buf.len() / 64 * 64;
            // Split the borrow: `chunks_exact` borrows `self.buf`, `compress` needs
            // `&mut self`, so the blocks are copied out first (64 bytes at a time).
            let mut i = 0;
            while i < whole {
                let mut block = [0u8; 64];
                block.copy_from_slice(&self.buf[i..i + 64]);
                self.compress(&block);
                i += 64;
            }
            self.buf.drain(..whole);
        }

        fn compress(&mut self, block: &[u8; 64]) {
            let mut w = [0u32; 64];
            for (i, word) in w[..16].iter_mut().enumerate() {
                *word = u32::from_be_bytes(block[4 * i..4 * i + 4].try_into().expect("4 bytes"));
            }
            for i in 16..64 {
                let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
                let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
                w[i] = w[i - 16]
                    .wrapping_add(s0)
                    .wrapping_add(w[i - 7])
                    .wrapping_add(s1);
            }
            let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = self.h;
            for i in 0..64 {
                let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
                let ch = (e & f) ^ ((!e) & g);
                let t1 = h
                    .wrapping_add(s1)
                    .wrapping_add(ch)
                    .wrapping_add(K[i])
                    .wrapping_add(w[i]);
                let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
                let maj = (a & b) ^ (a & c) ^ (b & c);
                let t2 = s0.wrapping_add(maj);
                h = g;
                g = f;
                f = e;
                e = d.wrapping_add(t1);
                d = c;
                c = b;
                b = a;
                a = t1.wrapping_add(t2);
            }
            for (slot, v) in self.h.iter_mut().zip([a, b, c, d, e, f, g, h]) {
                *slot = slot.wrapping_add(v);
            }
        }

        pub fn hex(mut self) -> String {
            // Pad: 0x80, then zeros to 56 mod 64, then the 64-bit big-endian bit length.
            let bitlen = self.len * 8;
            let mut tail = vec![0x80u8];
            while (self.buf.len() + tail.len()) % 64 != 56 {
                tail.push(0);
            }
            tail.extend_from_slice(&bitlen.to_be_bytes());
            self.update(&tail);
            debug_assert!(self.buf.is_empty(), "padding must fill the final block");
            self.h.iter().map(|v| format!("{v:08x}")).collect()
        }
    }
}

/// Every run directory under `root`, as paths relative to it.
///
/// A "run" is a directory holding a `command.json`, which is what the oracle records.
fn walk(root: &Path, only: Option<&str>) -> Vec<PathBuf> {
    let mut out = Vec::new();
    if !root.exists() {
        return out;
    }
    let mut stack = vec![root.to_path_buf()];
    while let Some(d) = stack.pop() {
        let Ok(rd) = std::fs::read_dir(&d) else {
            continue;
        };
        for e in rd.flatten() {
            let p = e.path();
            if p.is_dir() {
                stack.push(p);
            }
        }
    }
    out.sort();
    if let Some(f) = only {
        out.retain(|d| d.to_string_lossy().contains(f));
    }
    out.retain(|d| d.join("command.json").exists());
    out
}

/// `--stages` mode: per-stage max deviations between two `stages.json` trees.
fn run_stages(a: &Args) -> Result<ExitCode, String> {
    let rep = stages::compare_trees(&a.stages_golden, &a.stages_ours, a.stage_tol)?;
    if rep.stages.is_empty() {
        eprintln!(
            "macs-compare: no comparable stage dumps under {}",
            a.stages_golden.display()
        );
        return Err("nothing to compare".into());
    }
    println!("macs-compare stage summary");
    println!("  golden tree       {}", a.stages_golden.display());
    println!("  ours tree         {}", a.stages_ours.display());
    println!("  numeric tolerance {:.1e}\n", a.stage_tol);
    let mut bad = 0usize;
    for s in &rep.stages {
        let verdict = if s.is_clean() { "ok" } else { "FAIL" };
        if !s.is_clean() {
            bad += 1;
        }
        println!(
            "  {:<20} {:>5} {:>7} {:>7} {:>12} {:>12}  {}",
            s.name,
            s.runs,
            s.numeric,
            s.exact,
            format!("{:.3e}", s.max_abs),
            if s.max_rel.is_finite() {
                format!("{:.3e}", s.max_rel)
            } else {
                "inf".into()
            },
            verdict
        );
        for m in s.mismatches.iter().take(a.max_report) {
            println!("      - {m}");
        }
    }
    for n in &rep.run_notes {
        println!("  note: {n}");
    }
    let nums: usize = rep.stages.iter().map(|s| s.numeric).sum();
    let exact: usize = rep.stages.iter().map(|s| s.exact).sum();
    println!(
        "\n{} stages, {nums} numeric leaves, {exact} exact-match leaves, {bad} failing",
        rep.stages.len()
    );
    if let Some(f) = rep.first_failure() {
        println!("first divergence: {}", f.name);
    }
    Ok(if bad == 0 && rep.run_notes.is_empty() {
        ExitCode::SUCCESS
    } else {
        ExitCode::from(1)
    })
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("macs-compare: {e}");
            return ExitCode::from(2);
        }
    };

    if !args.stages_golden.as_os_str().is_empty() {
        return match run_stages(&args) {
            Ok(code) => code,
            Err(e) => {
                eprintln!("macs-compare: {e}");
                ExitCode::from(2)
            }
        };
    }

    let golden_runs = walk(&args.golden, args.only.as_deref());
    if golden_runs.is_empty() {
        eprintln!(
            "macs-compare: no golden runs found under {}",
            args.golden.display()
        );
        return ExitCode::from(2);
    }

    let mut compared = 0usize;
    let mut clean = 0usize;
    let mut records = 0usize;
    let mut missing = 0usize;
    let mut problems: Vec<String> = Vec::new();

    println!("macs-compare summary");
    println!("  golden runs       {}", golden_runs.len());
    println!("  ours root         {}", args.ours.display());
    println!();

    for gdir in &golden_runs {
        let rel = gdir
            .strip_prefix(&args.golden)
            .unwrap_or(gdir)
            .to_path_buf();
        let odir = args.ours.join(&rel);
        let Some(grun) = read_run(gdir).ok().filter(|_| !is_empty_run(gdir)) else {
            continue;
        };
        if !odir.join("command.json").exists() {
            missing += 1;
            continue;
        }
        let orun = match read_run(&odir) {
            Ok(r) => r,
            Err(e) => {
                problems.push(format!("{}: {e}", rel.display()));
                continue;
            }
        };
        compared += 1;

        // Structural: did the two runs produce the same files and exit the same way?
        let mut structural = Vec::new();
        for name in grun.files.keys() {
            if !orun.files.contains_key(name) {
                structural.push(format!("missing output {name}"));
            }
        }
        for name in orun.files.keys() {
            if !grun.files.contains_key(name) {
                structural.push(format!("extra output {name}"));
            }
        }
        if grun.returncode != orun.returncode {
            structural.push(format!(
                "exit status {:?} vs {:?}",
                grun.returncode, orun.returncode
            ));
        }
        let structural_ok = structural.is_empty();
        if !structural_ok {
            problems.push(format!("{}: {}", rel.display(), structural.join("; ")));
        }

        // Contents.
        let mut file_ok = structural_ok;
        for name in grun.files.keys() {
            let gp = gdir.join(name);
            let op = odir.join(name);
            let Some(fmt) = format_for(name, args.format.as_deref()) else {
                continue;
            };
            let (gt, ot) = match (std::fs::read_to_string(&gp), std::fs::read_to_string(&op)) {
                (Ok(a), Ok(b)) => (a, b),
                _ => continue,
            };
            let rep = compare(&gt, &ot, fmt);
            records += rep.macs3_count;
            let identical = grun.files.get(name) == orun.files.get(name);
            if !rep.is_clean() || (args.format.is_none() && !identical) {
                file_ok = false;
                if args.verbose || !rep.is_clean() {
                    println!("{} :: {}", rel.display(), name);
                    print!("{rep}");
                }
            }
        }
        if file_ok {
            clean += 1;
        }
    }

    println!("  runs compared     {compared}");
    println!("  runs clean        {clean}");
    println!("  records compared  {records}");
    println!("  missing runs      {missing}");
    if !problems.is_empty() {
        println!("\nstructural problems ({}):", problems.len());
        for p in problems.iter().take(args.max_report) {
            println!("  - {p}");
        }
    }
    let ok = clean == compared && missing == 0 && problems.is_empty();
    println!(
        "\n{}",
        if ok {
            "all compared runs agree with MACS3"
        } else {
            "divergences found"
        }
    );
    if ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::from(1)
    }
}

/// Is this a real run directory rather than a fixture or a stray file?
fn is_empty_run(dir: &Path) -> bool {
    std::fs::read_dir(dir)
        .map(|mut rd| rd.next().is_none())
        .unwrap_or(true)
}

/// Pick a record format for a file name, or `None` for a file that carries no
/// records (`command.json`, logs, and the stage dumps).
fn format_for(name: &str, forced: Option<&str>) -> Option<Format> {
    if let Some(f) = forced {
        return Some(Format::of(f));
    }
    // `Format::of` classifies by suffix; only files that carry records are worth
    // comparing, and it returns `Text` for anything else, which is not.
    match Format::of(name) {
        Format::Text => None,
        f => Some(f),
    }
}
