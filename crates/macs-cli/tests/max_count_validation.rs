//! `OptValidator.py`'s `--max-count` check, and the two other places a negative
//! cap is refused.
//!
//! Three subcommands read `--max-count`, and each handles a negative value
//! differently because it is refused at a different point:
//!
//! * `callpeak` -- `opt_validate_callpeak` refuses it outright for `-f FRAG`
//!   (`OptValidator.py:93-96`) and ignores the option for every other format.
//! * `pileup` -- `opt_validate_pileup` does the same, in its own FRAG branch
//!   (`OptValidator.py:547-552`).
//! * `hmmratac` -- **no validator check at all**. The negative cap survives to
//!   `FragParser.pe_parse_line`, whose `count` is declared `cython.ushort`
//!   (`Parser.py:1416`), so `count = min(count, max_count)` (`Parser.py:1482-1483`)
//!   makes Cython raise `OverflowError: can't convert negative value to unsigned
//!   short` on the first fragment line.
//!
//! `--extsize` is `pileup`'s other validator check (`OptValidator.py:560-563`,
//! `must be > 0`) and is included here because it is the same shape: a value the
//! writer would otherwise accept.
//!
//! Every exit status and file list was measured against the pinned oracle
//! (`OPENBLAS_CORETYPE=Haswell .oracle/venv/bin/macs3 ...`); nothing here needs
//! Python or the reference at run time.

use std::path::PathBuf;
use std::process::Command;

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

fn frag_fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/frag_basic/barcode_fragments")
}

struct Run {
    code: i32,
    files: Vec<String>,
}

fn run(tag: &str, argv: &[String]) -> Run {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/it-max-count")
        .join(tag);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let mut a = argv.to_vec();
    a.push("--outdir".into());
    a.push(dir.display().to_string());
    let st = Command::new(RUST).args(&a).output().unwrap();
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    files.sort();
    Run {
        code: st.status.code().unwrap_or(-1),
        files,
    }
}

fn s(x: &str) -> String {
    x.to_string()
}

#[test]
fn a_negative_max_count_is_refused_by_pileup_for_frag_only() {
    let frag = frag_fixture();
    let base = vec![
        s("pileup"),
        s("-i"),
        frag.join("treat.frag").display().to_string(),
        s("-f"),
        s("FRAG"),
        s("--extsize"),
        s("50"),
        s("-o"),
        s("p.bdg"),
    ];
    let mut neg = base.clone();
    neg.push("--max-count".into());
    neg.push("-1".into());
    let r = run("pileup-neg", &neg);
    assert_eq!(r.code, 1);
    assert_eq!(r.files, Vec::<String>::new());

    // a positive cap is accepted, and `--max-count 0` keeps every count
    for cap in ["1", "5", "0"] {
        let mut a = base.clone();
        a.push("--max-count".into());
        a.push(s(cap));
        let r = run(&format!("pileup-{cap}"), &a);
        assert_eq!(r.code, 0, "--max-count {cap} exited {}", r.code);
        assert_eq!(r.files, vec!["p.bdg".to_string()]);
    }
}

#[test]
fn hmmratac_refuses_a_negative_max_count_on_the_first_frag_line() {
    let frag = frag_fixture();
    let base = vec![
        s("hmmratac"),
        s("-i"),
        frag.join("treat.frag").display().to_string(),
        s("-f"),
        s("FRAG"),
        s("-n"),
        s("p"),
    ];
    let mut neg = base.clone();
    neg.push("--max-count".into());
    neg.push("-1".into());
    let r = run("hmm-neg", &neg);
    assert_eq!(r.code, 1, "exited {}", r.code);
    assert_eq!(
        r.files,
        Vec::<String>::new(),
        "the OverflowError happens while reading, before any output"
    );
}

#[test]
fn pileup_refuses_a_non_positive_extsize() {
    let frag = frag_fixture();
    for ext in ["0", "-5"] {
        let r = run(
            &format!("extsize{ext}"),
            &[
                s("pileup"),
                s("-i"),
                frag.join("treat.frag").display().to_string(),
                s("-f"),
                s("FRAG"),
                s("--extsize"),
                s(ext),
                s("-o"),
                s("p.bdg"),
            ],
        );
        assert_eq!(r.code, 1, "--extsize {ext} exited {}", r.code);
        assert_eq!(r.files, Vec::<String>::new());
    }
}
