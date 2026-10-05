//! `bdgcmp` without `-m` is rejected, and `bdgdiff` keeps going when
//! `--max-gap >= --min-len`.
//!
//! **`bdgcmp`.** `-m/--method` is declared `nargs="+"` with `default="ppois"`
//! (`bin/macs3:894-895`), so when the flag is absent argparse leaves the *string*
//! `"ppois"` in `options.method` -- not a one-element list. `opt_validate_bdgcmp`
//! then runs `for method in set(options.method)` (`OptValidator.py:583-587`), and
//! `set("ppois")` is the five characters `{'p','o','i','s'}`:
//!
//! ```text
//! $ macs3 bdgcmp -t t.bdg -c c.bdg -o out.bdg
//! ERROR ... Invalid method: p
//! $ echo $?
//! 1
//! ```
//!
//! Which character is reported depends on `set` iteration order, i.e. on the
//! hash seed, so only the exit status is reproducible. This port used to treat an
//! absent `-m` as the default method and run it, which is a different result
//! rather than a different message: it exited 0 and wrote the score bedGraph.
//!
//! **`bdgdiff`.** `bdgdiff_cmd.py:43` calls `error(...)` when
//! `options.maxgap >= options.minlen` and then **falls through** -- there is no
//! `sys.exit()`. The reference prints the complaint and runs anyway, exiting 0
//! with all three BED files written; this port refused the run.
//!
//! Every expectation was measured against the pinned oracle
//! (`OPENBLAS_CORETYPE=Haswell .oracle/venv/bin/macs3 ...`); nothing here needs
//! Python or the reference at run time.

use std::path::{Path, PathBuf};
use std::process::Command;

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

/// Two bedGraphs with a strictly positive control, so the only thing that can
/// reject the run is the option check under test.
fn bedgraphs(dir: &Path) -> (String, String, String) {
    let treat = "chr1\t0\t100\t10\nchr1\t100\t200\t4\nchr1\t400\t500\t7\n";
    let ctrl = "chr1\t0\t100\t5\nchr1\t100\t200\t2\nchr1\t400\t500\t3\n";
    let treat2 = "chr1\t0\t100\t9\nchr1\t100\t200\t30\nchr1\t400\t500\t2\n";
    let t = dir.join("t.bdg");
    let c = dir.join("c.bdg");
    let t2 = dir.join("t2.bdg");
    std::fs::write(&t, treat).unwrap();
    std::fs::write(&c, ctrl).unwrap();
    std::fs::write(&t2, treat2).unwrap();
    (
        t.display().to_string(),
        c.display().to_string(),
        t2.display().to_string(),
    )
}

struct Run {
    code: i32,
    files: Vec<String>,
}

struct Run2 {
    code: i32,
    files: Vec<String>,
    stderr: String,
}

fn bdgcmp(tag: &str, extra: &[&str]) -> Run2 {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/it-bdgcmp-no-method")
        .join(tag);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let (t, c, _) = bedgraphs(&dir);
    let mut argv: Vec<String> = vec![
        "bdgcmp".into(),
        "-t".into(),
        t,
        "-c".into(),
        c,
        "--outdir".into(),
        dir.display().to_string(),
    ];
    argv.extend(extra.iter().map(|s| (*s).to_string()));
    let st = Command::new(RUST).args(&argv).output().unwrap();
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        // the three inputs live in the same directory
        .filter(|n| !matches!(n.as_str(), "t.bdg" | "c.bdg" | "t2.bdg"))
        .collect();
    files.sort();
    Run2 {
        code: st.status.code().unwrap_or(-1),
        files,
        stderr: String::from_utf8_lossy(&st.stderr).into_owned(),
    }
}

#[test]
fn bdgcmp_without_a_method_writes_nothing_and_exits_one() {
    for (tag, extra) in [
        ("ofile", vec!["-o", "out.bdg"]),
        ("oprefix", vec!["--o-prefix", "zz"]),
    ] {
        let r = bdgcmp(tag, &extra);
        assert_eq!(r.code, 1, "bdgcmp {extra:?} exited {}", r.code);
        assert_eq!(
            r.files,
            Vec::<String>::new(),
            "bdgcmp {extra:?} wrote {:?}; the method check runs before the score track",
            r.files
        );
        // `set("ppois")` iteration is hash-randomised, so only the shape of the
        // message is fixed: a one-character method name.
        assert!(
            r.stderr.contains("Invalid method: "),
            "unexpected diagnostic: {:?}",
            r.stderr
        );
    }
}

#[test]
fn bdgcmp_with_an_explicit_method_is_unaffected() {
    let r = bdgcmp("explicit", &["-m", "ppois", "-o", "out.bdg"]);
    assert_eq!(r.code, 0);
    assert_eq!(r.files, vec!["out.bdg".to_string()]);
}

fn bdgdiff(tag: &str, extra: &[&str]) -> Run {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/it-bdgdiff-maxgap")
        .join(tag);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let (t, c, t2) = bedgraphs(&dir);
    let mut argv: Vec<String> = vec![
        "bdgdiff".into(),
        "--t1".into(),
        t,
        "--t2".into(),
        t2,
        "--c1".into(),
        c.clone(),
        "--c2".into(),
        c,
        "-o".into(),
        dir.join("a.bed").display().to_string(),
        dir.join("b.bed").display().to_string(),
        dir.join("cc.bed").display().to_string(),
        "--outdir".into(),
        dir.display().to_string(),
    ];
    argv.extend(extra.iter().map(|s| (*s).to_string()));
    let st = Command::new(RUST).args(&argv).output().unwrap();
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".bed"))
        .collect();
    files.sort();
    Run {
        code: st.status.code().unwrap_or(-1),
        files,
    }
}

#[test]
fn bdgdiff_reports_an_inverted_window_pair_but_still_runs() {
    // `-l 100 -g 200`: MAXGAP 200 >= MINLEN 100, which upstream complains about and
    // then ignores. Exit 0, all three BED files written.
    let r = bdgdiff("inverted", &["-l", "100", "-g", "200"]);
    assert_eq!(r.code, 0);
    assert_eq!(
        r.files,
        vec![
            "a.bed".to_string(),
            "b.bed".to_string(),
            "cc.bed".to_string()
        ]
    );
    // ... and a legal pair behaves the same way
    let ok = bdgdiff("legal", &["-l", "200", "-g", "100"]);
    assert_eq!(ok.code, 0);
    assert_eq!(ok.files, r.files);
}
