//! Build fermi-lite from the vendored C sources.
//!
//! `crates/macs-callvar/vendor/fermi-lite` is MACS3's own fermi-lite submodule,
//! vendored verbatim at the commit the pinned oracle pins (see
//! `crates/macs-callvar/vendor/fermi-lite/PINNED.md`), so the build does not
//! depend on the oracle tree being present and the assembler behaves identically
//! to the one `macs3 callvar` uses.
//!
//! The vendored tree lives *inside* this package rather than at the workspace
//! root because Cargo can only package files beneath the package directory. When
//! it sat at the workspace root, `cargo package` shipped a `build.rs` with no C
//! sources beside it, and the build script's `assert!(path.exists())` below would
//! have panicked for every downstream consumer.
//!
//! `PORTING_PLAN.md` allows exactly this: *"callvar may initially bridge to upstream
//! fermi-lite via thin C FFI; a pure-Rust assembler is the last item, after every other
//! gate is green."*
//!
//! The compiled set is fermi-lite's own `OBJS` line verbatim -- all thirteen files.
//! Trimming it is tempting (`unitig.c` and `rle.c` look sufficient) and does not link:
//! `unitig.c` references `rld_extend` from `rld0.o`.
//!
//! Two defines matter and are set here rather than guessed at the call site:
//!
//! * `-DFML_USE_OPENMP` would let fermi-lite spawn threads. Upstream never enables it,
//!   so unitig output order would stop being deterministic; it is deliberately off.
//! * `-w` suppresses the vendored library's warnings, which are not ours to fix and
//!   would otherwise drown out `cargo clippy`.

use std::path::PathBuf;

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fml = root.join("vendor/fermi-lite");

    println!("cargo:rerun-if-changed={}", fml.display());
    println!("cargo:rerun-if-changed=build.rs");

    let mut build = cc::Build::new();
    build
        .include(&fml)
        .flag_if_supported("-w")
        .flag_if_supported("-fPIC")
        // fermi-lite is C99 and uses POSIX headers under `kseq.h`/`kthread.c`.
        .define("_GNU_SOURCE", None)
        .opt_level(3);

    // fermi-lite's Makefile: OBJS = kthread.o misc.o bseq.o htab.o bfc.o rle.o
    //                      rope.o mrope.o rld0.o unitig.o mag.o bubble.o ksw.o
    for src in [
        // `swalign.c` is not part of fermi-lite; it is MACS3's own Smith-Waterman,
        // which `RACollection.align_unitig_to_REFSEQ` uses to place a unitig back
        // on the peak consensus. Vendored alongside because `callvar` needs both,
        // and both are the pinned oracle's own code.
        "kthread.c",
        "misc.c",
        "bseq.c",
        "htab.c",
        "bfc.c",
        "rle.c",
        "rope.c",
        "mrope.c",
        "rld0.c",
        "unitig.c",
        "mag.c",
        "bubble.c",
        "ksw.c",
    ] {
        let path = fml.join(src);
        assert!(
            path.exists(),
            "vendored fermi-lite source missing: {}",
            path.display()
        );
        build.file(path);
    }

    build.compile("fml");

    // `swalign.c` ships a standalone CLI (`usage: swalign TARGET_SEQ QUERY_SEQ`), whose
    // `main` would collide with the Rust binary's. The vendored bytes stay untouched:
    // the symbol is renamed with `-D`, which is a build flag rather than an edit, so the
    // file in `vendor/` remains byte-identical to the oracle's and `git diff` on it shows
    // nothing. Its own `print_alignment` becomes unreferenced, which `-w` hides.
    let mut sa = cc::Build::new();
    sa.include(&fml)
        .flag_if_supported("-w")
        .flag_if_supported("-fPIC")
        .define("_GNU_SOURCE", None)
        .flag("-Dmain=swalign_cli_main_unused")
        .opt_level(3)
        .file(fml.join("swalign.c"))
        .compile("swalign");
    drop(sa);

    for h in [
        "fml.h",
        "internal.h",
        "kseq.h",
        "kvec.h",
        "ksort.h",
        "kmer.h",
        "kstring.h",
        "mrope.h",
        "rle.h",
        "rope.h",
    ] {
        let p = fml.join(h);
        if p.exists() {
            println!("cargo:rerun-if-changed={}", p.display());
        }
    }
    println!("cargo:include={}", fml.display());
}
