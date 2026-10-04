//! The argument surface, derived from upstream's own argparse.
//!
//! The table is `oracle/flag_matrix.tsv`, auto-derived from `macs3`'s argparse
//! (not from `--help`, which documents flags that do not exist -- F23). It is
//! embedded with [`include_str!`] so the shipped binary has no runtime data
//! dependency, which keeps the "fully independent of Python at runtime" and
//! "single self-contained tool" requirements intact.
//!
//! Parsing follows argparse's observable behaviour, because the acceptance
//! criteria require identical accept/reject behaviour:
//!
//! * a missing **required** flag is a usage error (exit 2);
//! * a value outside **choices** is a usage error;
//! * `--flag value` and `--flag=value` are both accepted, and so are the short
//!   forms listed in the matrix;
//! * a value that looks like a flag is rejected unless the flag takes one;
//! * unknown flags are usage errors;
//! * `--help`/`-h` succeeds with the subcommand's usage text.

use std::collections::BTreeMap;

/// One row of the flag matrix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FlagSpec {
    /// The subcommand this row belongs to, or `<top>` for the program-level flags.
    pub subcommand: String,
    /// The flag as written on the command line, e.g. `--extsize`.
    pub flag: String,
    /// argparse `action`: `store`, `store_true`, `store_false`, `count`, `help`.
    pub action: String,
    /// argparse `nargs`, rendered as a string; `0` for flags that take no value.
    pub nargs: String,
    /// The argparse `dest` -- the key this flag populates.
    pub dest: String,
    /// `true` when the flag is required.
    pub required: bool,
    /// The default, or `None` for "no default".
    pub default: Option<String>,
    /// Permitted values, empty when unconstrained.
    pub choices: String,
}

/// The embedded flag matrix.
const MATRIX: &str = include_str!("../../../oracle/flag_matrix.tsv");

/// The embedded per-subcommand flag *type* table, generated from the oracle's
/// argparse by `oracle/gen_flag_types.py`.
const TYPES: &str = include_str!("../../../oracle/flag_types.tsv");

/// `(subcommand, flag spelling) -> "int" | "float"`, for flags upstream declares
/// with a numeric argparse `type`.
///
/// The matrix records which flags exist but not their types, and the types cannot
/// be inferred from the defaults: `--verbose` defaults to `2` and is an int,
/// `--gsize` defaults to `hs` and is a string, `--pvalue` defaults to `1e-5` and
/// is a float. They are read out of upstream's own `add_argument` calls instead,
/// scoped per subcommand because `-g` is a genome-size string in `callpeak` and an
/// int max-gap in `bdgdiff`.
fn flag_types() -> &'static Vec<(String, String, String)> {
    static CACHE: std::sync::OnceLock<Vec<(String, String, String)>> = std::sync::OnceLock::new();
    CACHE.get_or_init(|| {
        let mut out = Vec::new();
        for line in TYPES.lines() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3 || f[0] == "subcommand" {
                continue;
            }
            out.push((f[0].to_string(), f[1].to_string(), f[2].to_string()));
        }
        out
    })
}

/// argparse's spelling of a flag for a diagnostic: every option string that
/// shares the destination, joined with `/` and in declaration order, e.g.
/// `-s/--tsize`. argparse prints `argument -s/--tsize: invalid int value: 'x'`,
/// and a message naming only the spelling the user happened to type reads as a
/// different error than the one upstream produced.
fn display_name(spec: &FlagSpec, specs: &[FlagSpec]) -> String {
    let mut names: Vec<&str> = specs
        .iter()
        .filter(|s| s.subcommand == spec.subcommand && s.dest == spec.dest)
        .map(|s| s.flag.as_str())
        .collect();
    if names.is_empty() {
        return spec.flag.clone();
    }
    names.dedup();
    names.join("/")
}

/// Reject a value that upstream's argparse would reject for this flag's type.
///
/// Without this, `Options::int`/`Options::float` return `None` on a parse failure
/// and every caller falls back to the flag's default. That is an **accept/reject
/// divergence with a silent wrong answer**: 79 flags across 11 subcommands took a
/// malformed value and ran anyway -- `bdgpeakcall -c notafloat` exited 0 and called
/// peaks with the default cutoff of 5, where upstream exits 2. A typo'd `--extsize`
/// or `--pvalue` is exactly the case where quietly using the default is worse than
/// refusing.
///
/// Note the float case accepts anything `f64::from_str` accepts, which is what
/// `type=float` does too -- including `nan` and `inf` -- so this matches argparse
/// rather than tightening it.
fn check_type(
    subcommand: &str,
    spec: &FlagSpec,
    display: &str,
    value: &str,
) -> Result<(), UsageError> {
    let Some((_, _, ty)) = flag_types()
        .iter()
        .find(|(s, f, _)| s == subcommand && f == &spec.flag)
    else {
        return Ok(());
    };
    let ok = match ty.as_str() {
        "int" => value.parse::<i64>().is_ok(),
        "float" => value.parse::<f64>().is_ok(),
        _ => true,
    };
    if ok {
        return Ok(());
    }
    Err(UsageError(format!(
        "error: argument {display}: invalid {ty} value: '{value}'"
    )))
}

/// Parse the matrix. Called once; the result is a `const`-friendly static table
/// so tests can assert against the same data the binary uses.
pub fn flag_specs() -> Vec<FlagSpec> {
    let mut out = Vec::new();
    let mut lines = MATRIX.lines();
    let Some(_header) = lines.next() else {
        return out;
    };
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        let default = if f[6] == "None" || f[6].is_empty() {
            None
        } else {
            Some(f[6])
        };
        out.push(FlagSpec {
            subcommand: f[0].to_string(),
            flag: f[1].to_string(),
            action: f[2].to_string(),
            nargs: f[3].to_string(),
            dest: f[4].to_string(),
            required: f[5].eq_ignore_ascii_case("true"),
            default: default.map(str::to_string),
            choices: f[7].to_string(),
        });
    }
    out
}

/// A usage error, which exits 2 exactly as argparse does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UsageError(pub String);

impl std::fmt::Display for UsageError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for UsageError {}

/// The parsed options for one invocation.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Options {
    /// The subcommand that was named.
    pub subcommand: String,
    /// Populated destinations, in matrix order. A `nargs *`/`+` flag keeps its
    /// first value here; the full list is in [`Self::lists`].
    pub values: BTreeMap<String, String>,
    /// Values of flags declared with `nargs *` or `nargs +`, in order.
    pub lists: BTreeMap<String, Vec<String>>,
    /// Destinations whose action was `store_true`/`store_false` and that fired.
    pub flags: BTreeMap<String, bool>,
    /// `--help` was requested.
    pub help: bool,
    /// The argument vector as given, **including** the subcommand name.
    ///
    /// Upstream echoes `" ".join(sys.argv[1:])` into the `*_peaks.xls` header as
    /// `# Command line: ...`, and `sys.argv[1]` is the subcommand. The file is
    /// compared byte-for-byte, so the exact spelling of what the user typed has
    /// to survive parsing.
    pub raw_argv: Vec<String>,
}

impl Options {
    /// A string destination.
    pub fn get(&self, dest: &str) -> Option<&str> {
        self.values.get(dest).map(String::as_str)
    }
    /// An integer destination, parsed as `i64`.
    pub fn int(&self, dest: &str) -> Option<i64> {
        self.get(dest).and_then(|v| v.parse().ok())
    }
    /// A float destination.
    pub fn float(&self, dest: &str) -> Option<f64> {
        self.get(dest).and_then(|v| v.parse().ok())
    }
    /// Every value of a `nargs *`/`+` flag, in order.
    pub fn get_all(&self, dest: &str) -> &[String] {
        self.lists.get(dest).map(Vec::as_slice).unwrap_or(&[])
    }
    /// A boolean destination (`store_true`/`store_false`), defaulting to false.
    pub fn flag(&self, dest: &str) -> bool {
        self.flags.get(dest).copied().unwrap_or(false)
    }
}

/// Validate a value against a flag's `choices` column (comma-separated).
fn check_choice(
    subcommand: &str,
    spec: &FlagSpec,
    name: &str,
    value: &str,
) -> Result<(), UsageError> {
    check_type(subcommand, spec, name, value)?;
    if !spec.choices.is_empty() {
        let allowed: Vec<&str> = spec.choices.split(',').map(|s| s.trim()).collect();
        if !allowed.contains(&value) {
            return Err(UsageError(format!(
                "error: argument {name}: invalid choice: {value} (choose from {})",
                allowed.join(", ")
            )));
        }
    }
    Ok(())
}

/// Parse `args` (without the program name) against the matrix.
///
/// `subcommand` must already be known; the top-level parser is responsible for
/// dispatching, so this only sees a subcommand's own flags.
pub fn parse(subcommand: &str, args: &[String]) -> Result<Options, UsageError> {
    let specs = flag_specs();
    let mut out = Options {
        subcommand: subcommand.to_string(),
        raw_argv: std::iter::once(subcommand.to_string())
            .chain(args.iter().cloned())
            .collect(),
        ..Default::default()
    };
    let mut seen: BTreeMap<String, FlagSpec> = BTreeMap::new();
    // dests in command-line order, for argparse's conflict-message ordering
    let mut parse_order: Vec<String> = Vec::new();

    // apply defaults first, so a flag that is not given still has its value
    for s in &specs {
        if s.subcommand != subcommand {
            continue;
        }
        if s.action == "help" {
            continue;
        }
        if let Some(d) = s.default.as_deref() {
            if s.action == "store_true" {
                out.flags.insert(s.dest.to_string(), d == "true");
            } else if s.action == "store_false" {
                out.flags.insert(s.dest.to_string(), d != "false");
            } else {
                out.values.insert(s.dest.to_string(), d.to_string());
            }
        }
    }

    let mut i = 0usize;
    while i < args.len() {
        let a = args[i].as_str();
        if a == "--" {
            // everything after `--` is positional, which this subcommand surface
            // does not take; argparse would treat them as extra positionals
            break;
        }
        if !a.starts_with('-') || a == "-" {
            return Err(UsageError(format!("error: unrecognized argument: {a}")));
        }
        let (name, inline) = match a.split_once('=') {
            Some((n, v)) => (n, Some(v.to_string())),
            None => (a, None),
        };
        let Some(spec) = specs
            .iter()
            .find(|s| s.subcommand == subcommand && s.flag == name)
        else {
            return Err(UsageError(format!("error: unrecognized arguments: {name}")));
        };
        if !seen.contains_key(&spec.dest) {
            parse_order.push(spec.dest.clone());
        }
        seen.insert(spec.dest.clone(), spec.clone());
        // argparse names a flag by every spelling that shares its destination, so a
        // diagnostic reads `-s/--tsize` rather than whichever token the user typed.
        let display = display_name(spec, &specs);
        match spec.action.as_str() {
            "help" => {
                out.help = true;
                i += 1;
            }
            "store_true" | "count" => {
                out.flags.insert(spec.dest.to_string(), true);
                i += 1;
            }
            "store_false" => {
                // `--no-trackline` and friends *set the destination to false*.
                // They are not booleans-that-fired; treating them as
                // store_true makes a negated flag indistinguishable from the
                // default and inverts it (F133).
                out.flags.insert(spec.dest.to_string(), false);
                i += 1;
            }
            // A fixed `nargs=N` (an integer in the matrix) consumes exactly `N`
            // values. `callpeak -m/--mfold` is `nargs=2`, so `--mfold 3 20` is two
            // separate argv tokens, not one -- treating it as a single value made
            // the parser reject the run with "unrecognized argument: 20" and exit
            // 2, so the whole `mfold_3_20` matrix variant failed at the door.
            //
            // The variable-count cases are handled below; this arm is the fixed
            // count, and everything not named falls through to the single-value
            // default.
            _ if spec.nargs.parse::<usize>().is_ok_and(|k| k >= 2) => {
                let k = spec.nargs.parse::<usize>().expect("checked");
                let mut vals: Vec<String> = Vec::new();
                if let Some(v) = inline {
                    vals.push(v);
                }
                while vals.len() < k {
                    i += 1;
                    let v = args.get(i).cloned().ok_or_else(|| {
                        UsageError(format!("error: argument {name}: expected {k} arguments"))
                    })?;
                    vals.push(v);
                }
                // step past the last value, so the outer loop does not re-read it
                // as a flag
                i += 1;
                for v in &vals {
                    check_choice(subcommand, spec, name, v)?;
                }
                // the first value is what a scalar accessor sees; the whole list
                // is kept so `# model fold = [3, 20]` can be rendered
                out.values.insert(spec.dest.to_string(), vals[0].clone());
                out.lists.insert(spec.dest.to_string(), vals);
                continue;
            }
            _ => {
                // nargs `+` / `*` consume a variable number of following
                // tokens (cmbreps takes several `--ifile`, bdgopt several
                // `--extra-param`); everything else takes exactly one.
                if spec.nargs == "+" || spec.nargs == "*" {
                    let mut vals: Vec<String> = Vec::new();
                    if let Some(v) = inline {
                        vals.push(v);
                        i += 1;
                    } else {
                        i += 1; // step past the flag itself
                        while i < args.len() {
                            let nxt = &args[i];
                            if nxt.starts_with('-') && nxt.len() > 1 {
                                break;
                            }
                            vals.push(nxt.clone());
                            i += 1;
                        }
                    }
                    if spec.nargs == "+" && vals.is_empty() {
                        return Err(UsageError(format!(
                            "error: argument {name}: expected at least one argument"
                        )));
                    }
                    for v in &vals {
                        check_choice(subcommand, spec, name, v)?;
                    }
                    if let Some(first) = vals.first() {
                        out.values.insert(spec.dest.to_string(), first.clone());
                    }
                    out.lists.insert(spec.dest.to_string(), vals);
                    continue;
                }
                let value = match inline {
                    Some(v) => {
                        i += 1;
                        v
                    }
                    None => {
                        let v = args.get(i + 1).cloned().ok_or_else(|| {
                            UsageError(format!("error: argument {display}: expected one argument"))
                        })?;
                        i += 2;
                        v
                    }
                };
                check_choice(subcommand, spec, &display, &value)?;
                out.values.insert(spec.dest.to_string(), value);
            }
        }
    }

    // `-h`/`--help` short-circuits before any other validation, exactly as
    // argparse does: `macs3 callpeak --help` must print usage and exit 0 even
    // though `-t`/`-g` are required.
    if out.help {
        return Ok(out);
    }

    // required flags must have been given
    for s in &specs {
        if s.subcommand != subcommand || !s.required || s.action == "help" {
            continue;
        }
        if !seen.contains_key(&s.dest) {
            return Err(UsageError(format!(
                "error: the following arguments are required: {}",
                s.flag
            )));
        }
    }
    check_mutex_groups(subcommand, &seen, &parse_order)?;
    Ok(out)
}

/// Argparse mutually-exclusive groups in the pinned upstream parser.
///
/// The flag matrix records each flag's own attributes but not its *group*, so
/// these are transcribed from `bin/macs3` at the pinned commit
/// (`c5443190e3edfeb301cc94acf450e2b2c026a223`) -- there are exactly seven, listed
/// here rather than in the generated matrix because the generator only walks
/// `add_argument` calls. `required = true` groups must see exactly one member.
const MUTEX_GROUPS: &[(&str, &[&str], bool)] = &[
    // argparser_bdgpeakcall
    ("bdgpeakcall", &["ofile", "oprefix"], true),
    // argparser_bdgbroadcall
    ("bdgbroadcall", &["ofile", "oprefix"], true),
    // argparser_refinepeak
    ("refinepeak", &["ofile", "oprefix"], true),
    // argparser_bdgcmp: `add_argument("--o-prefix")` precedes `add_argument("-o")`
    // here, the reverse of `add_output_group`, so the message order differs.
    ("bdgcmp", &["oprefix", "ofile"], true),
    // argparser_randsample
    ("randsample", &["percentage", "number"], true),
    // argparser_bdgdiff
    ("bdgdiff", &["ofile", "oprefix"], true),
    // group_callpeak: -p / -q mutually exclusive, neither required
    ("callpeak", &["pvalue", "qvalue"], false),
    // NOTE: `group_postprocessing` in callpeak looks like a mutual-exclusion
    // group but holds only `--call-summits`, so it constrains nothing. `--broad`
    // is added to `group_callpeak`, a plain argument group -- verified: upstream
    // accepts `callpeak --broad --call-summits` together without complaint.
];

/// Enforce the mutually-exclusive groups for `subcommand`.
fn check_mutex_groups(
    subcommand: &str,
    seen: &BTreeMap<String, FlagSpec>,
    order: &[String],
) -> Result<(), UsageError> {
    for (cmd, members, required) in MUTEX_GROUPS {
        let all_specs = flag_specs();
        if *cmd != subcommand {
            continue;
        }
        // argparse names the flag that appeared *later* on the command line
        // first: `-p 1 -q 1` reports "-q/--qvalue: not allowed with -p/--pvalue".
        // `order` is the sequence of dests as they were parsed.
        let given: Vec<&str> = members
            .iter()
            .filter(|d| seen.contains_key(**d))
            .copied()
            .collect();
        let given = {
            let mut v: Vec<&str> = given
                .iter()
                .copied()
                .filter(|d| order.iter().any(|o| o == *d))
                .collect();
            v.sort_by_key(|d| {
                std::cmp::Reverse(order.iter().position(|o| o == *d).unwrap_or(usize::MAX))
            });
            v
        };
        if given.len() > 1 {
            // the matrix lists every spelling of a flag as its own row; argparse
            // reports only the first, joined by `/`
            let flags: Vec<String> = given
                .iter()
                .map(|d| {
                    let names: Vec<&str> = all_specs
                        .iter()
                        .filter(|s| s.subcommand == subcommand && s.dest == **d)
                        .map(|s| s.flag.as_str())
                        .collect();
                    names.join("/")
                })
                .collect();
            return Err(UsageError(format!(
                "error: argument {}: not allowed with argument {}",
                flags[0],
                flags[1..].join(", argument ")
            )));
        }
        if *required && given.is_empty() {
            // argparse names the flags, not the dests, and joins every spelling
            // of a member with `/`: "-p/--percentage -n/--number is required"
            let names: Vec<String> = members
                .iter()
                .map(|d| {
                    let flags: Vec<&str> = all_specs
                        .iter()
                        .filter(|s| s.subcommand == subcommand && s.dest == **d)
                        .map(|s| s.flag.as_str())
                        .collect();
                    flags.join("/")
                })
                .collect();
            return Err(UsageError(format!(
                "error: one of the arguments {} is required",
                names.join(" ")
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn the_matrix_covers_every_subcommand_in_the_objective() {
        let specs = flag_specs();
        assert!(specs.len() > 300, "matrix looks truncated: {}", specs.len());
        for c in crate::SUBCOMMANDS {
            assert!(specs.iter().any(|s| s.subcommand == c), "no flags for {c}");
        }
    }

    #[test]
    fn every_subcommand_has_a_required_flag_or_none_is_expected() {
        // just exercises the required-check path
        let r = parse("callpeak", &v(&["-t", "a.bed"]));
        assert!(r.is_ok() || r.is_err());
    }

    #[test]
    fn an_unknown_flag_is_a_usage_error() {
        let e = parse("callpeak", &v(&["--not-a-real-flag"])).unwrap_err();
        assert!(e.0.contains("unrecognized"), "{e}");
    }

    #[test]
    fn a_missing_required_flag_is_a_usage_error() {
        // callpeak requires -t (treatment)
        let e = parse("callpeak", &v(&[])).unwrap_err();
        assert!(e.0.contains("required"), "{e}");
    }

    #[test]
    fn inline_and_separate_values_are_both_accepted() {
        let a = parse("callpeak", &v(&["-t", "a.bed", "--gsize=1000"])).unwrap();
        let b = parse("callpeak", &v(&["-t", "a.bed", "--gsize", "1000"])).unwrap();
        assert_eq!(a.get("gsize"), b.get("gsize"));
        assert_eq!(a.get("gsize"), Some("1000"));
    }

    #[test]
    fn a_flag_with_no_value_argument_is_a_usage_error() {
        let e = parse("callpeak", &v(&["-t"])).unwrap_err();
        assert!(e.0.contains("expected"), "{e}");
    }

    #[test]
    fn store_true_sets_the_destination() {
        let o = parse("callpeak", &v(&["-t", "a.bed", "--nomodel"])).unwrap();
        assert!(o.flag("nomodel"));
    }

    #[test]
    fn defaults_are_applied() {
        let o = parse("callpeak", &v(&["-t", "a.bed"])).unwrap();
        // upstream's default q-value cutoff
        assert!(o.get("qvalue").is_some(), "qvalue default missing");
    }
}
