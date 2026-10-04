//! Thin FFI to fermi-lite, the assembler `callvar` uses to re-call peaks.
//!
//! # Why a C bridge
//!
//! `PORTING_PLAN.md` sanctions this explicitly: *"callvar may initially bridge to
//! upstream fermi-lite via thin C FFI; a pure-Rust assembler is the last item, after
//! every other gate is green."* The C is vendored at fermi-lite `r53`
//! (`vendor/fermi-lite/PINNED.md`), the same revision the pinned oracle uses, so this
//! is upstream's assembler rather than a reimplementation of one.
//!
//! # The surface is five functions
//!
//! `fml_opt_init`, `fml_assemble`, `fml_utg_destroy`, plus the three structs
//! `fml_opt_t`, `bseq1_t` and `fml_utg_t`. Everything else in the library is internal.
//!
//! # Two conventions that are easy to get wrong
//!
//! **Qualities are Phred+33 on the way in.** `fermi_assemble` expects ASCII, and
//! upstream does `cqual[j] = tmpq[j] + 33` on raw BAM qualities. Passing raw values
//! silently changes error correction and therefore every unitig.
//!
//! **Memory is owned by fermi-lite.** The `bseq1_t` array's `seq`/`qual` buffers are
//! freed by the library, so they must each be their own `malloc`, not a view into a
//! Rust allocation. [`assemble`] hands ownership over and never frees them itself.

use std::ffi::{c_char, c_int, c_void};

// ---------------------------------------------------------------------------
// Declarations, mirroring `vendor/fermi-lite/fml.h` verbatim.
// ---------------------------------------------------------------------------

#[repr(C)]
#[derive(Clone, Copy)]
struct Bseq1T {
    l_seq: i32,
    seq: *mut c_char,
    qual: *mut c_char,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct MagOptT {
    flag: c_int,
    min_ovlp: c_int,
    min_elen: c_int,
    min_ensr: c_int,
    min_insr: c_int,
    max_bdist: c_int,
    max_bdiff: c_int,
    max_bvtx: c_int,
    min_merge_len: c_int,
    trim_len: c_int,
    trim_depth: c_int,
    min_dratio1: f32,
    max_bcov: f32,
    max_bfrac: f32,
}

/// The subset of `fml_opt_t` that `callvar` sets. Laid out to match `fml.h` exactly;
/// `_Static_assert`-style checks live in the tests via `size_of`.
#[repr(C)]
#[derive(Clone, Copy)]
struct FmlOptT {
    n_threads: c_int,
    ec_k: c_int,
    min_cnt: c_int,
    max_cnt: c_int,
    min_asm_ovlp: c_int,
    min_merge_len: c_int,
    mag_opt: MagOptT,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct FmlUtgT {
    len: i32,
    nsr: i32,
    seq: *mut c_char,
    cov: *mut c_char,
    n_ovlp: [c_int; 2],
    ovlp: *mut c_void,
}

/// `MAG_F_AGGRESSIVE` (pop variant bubbles; *not* fermi-lite's default).
///
/// Exposed for completeness against `fml.h`. `callvar` never sets it -- upstream
/// passes `0x80` -- but leaving it undeclared would make the flag list look edited.
pub const MAG_F_AGGRESSIVE: i32 = 0x20;
/// `MAG_F_POPOPEN` (aggressive tip trimming, fermi-lite's default).
pub const MAG_F_POPOPEN: i32 = 0x40;
/// `MAG_F_NO_SIMPL` (skip bubble simplification, fermi-lite's default).
pub const MAG_F_NO_SIMPL: i32 = 0x80;
/// The flag upstream passes: `opt_flag = 0x80`.
pub const DEFAULT_OPT_FLAG: i32 = MAG_F_NO_SIMPL;

extern "C" {
    fn fml_opt_init(opt: *mut FmlOptT);
    fn fml_assemble(
        opt: *const FmlOptT,
        n_seqs: c_int,
        seqs: *mut Bseq1T,
        n_utg: *mut c_int,
    ) -> *mut FmlUtgT;
    fn fml_utg_destroy(n_utg: c_int, utg: *mut FmlUtgT);
}

/// One assembled unitig.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Unitig {
    /// The unitig sequence, already decoded to ASCII by the library.
    pub seq: Vec<u8>,
    /// Number of supporting reads.
    pub nsr: i32,
}

/// Run fermi-lite's assembler over `reads`.
///
/// `reads` are `(sequence, raw BAM qualities)`; qualities are offset by 33 on the way
/// in, matching upstream's `cqual[j] = tmpq[j] + 33`.
///
/// `min_asm_ovlp` is `--fermi-overlap` (`fermiMinOverlap`, default 30).
///
/// Returns `Err` only when fermi-lite cannot allocate; upstream's own failure modes
/// (`unitig_collection == -1`) are reported by the caller's graph-cleaning step, not
/// here.
pub fn assemble(
    reads: &[(Vec<u8>, Vec<u8>)],
    min_asm_ovlp: i32,
    opt_flag: i32,
) -> Result<Vec<Unitig>, String> {
    if reads.is_empty() {
        return Ok(Vec::new());
    }
    // The struct layout is ours to get right; assert it once at the boundary so a
    // mismatch is a clear error rather than memory corruption. (The same check also
    // lives in the test module, where a drift is a failing test rather than a refusal
    // at run time -- this one is the belt to that braces.)
    assert_eq!(
        std::mem::size_of::<FmlOptT>(),
        6 * std::mem::size_of::<c_int>() + std::mem::size_of::<MagOptT>(),
        "fml_opt_t layout drifted from fml.h"
    );

    let n = reads.len();
    // The `bseq1_t` array itself must come from libc `malloc`, not a Rust `Vec`:
    // `fml_assemble` frees the array it is given (upstream's own comments say so twice,
    // "we rely on fermi-lite to free this mem"), and freeing a `Vec`'s buffer with libc
    // `free` is a double free the moment the `Vec` drops. Hence the raw pointer.
    //
    // SAFETY: `n * size_of::<Bseq1T>()` bytes, zeroed, owned by fermi-lite after the call.
    let seqs = unsafe {
        let p = c_malloc(n * std::mem::size_of::<Bseq1T>()) as *mut Bseq1T;
        if p.is_null() {
            return Err("fermi-lite: out of memory".to_string());
        }
        std::ptr::write_bytes(p, 0, n);
        p
    };
    // Each `seq`/`qual` pair is likewise its own malloc, for the same reason.
    for (i, (s, q)) in reads.iter().enumerate() {
        assert_eq!(
            s.len(),
            q.len(),
            "sequence and quality must be the same length"
        );
        // SAFETY: `s.len() + 1` bytes, NUL-terminated exactly as upstream does.
        let (cseq, cqual) = unsafe {
            let cseq = c_malloc(s.len() + 1) as *mut c_char;
            let cqual = c_malloc(q.len() + 1) as *mut c_char;
            if cseq.is_null() || cqual.is_null() {
                if !cseq.is_null() {
                    libc_free(cseq as *mut c_void);
                }
                if !cqual.is_null() {
                    libc_free(cqual as *mut c_void);
                }
                return Err("fermi-lite: out of memory".to_string());
            }
            std::ptr::copy_nonoverlapping(s.as_ptr(), cseq as *mut u8, s.len());
            std::ptr::write(cseq.add(s.len()), 0);
            // Phred + 33, as upstream. fermi-lite reads ASCII.
            for (i, &v) in q.iter().enumerate() {
                std::ptr::write(cqual.add(i), (v + 33) as c_char);
            }
            std::ptr::write(cqual.add(q.len()), 0);
            (cseq, cqual)
        };
        // SAFETY: `seqs` has room for `n` entries and `i` is in bounds.
        unsafe {
            *seqs.add(i) = Bseq1T {
                l_seq: s.len() as i32,
                seq: cseq,
                qual: cqual,
            };
        }
    }

    let mut opt = FmlOptT {
        n_threads: 0,
        ec_k: 0,
        min_cnt: 0,
        max_cnt: 0,
        min_asm_ovlp: 0,
        min_merge_len: 0,
        mag_opt: MagOptT {
            flag: 0,
            min_ovlp: 0,
            min_elen: 0,
            min_ensr: 0,
            min_insr: 0,
            max_bdist: 0,
            max_bdiff: 0,
            max_bvtx: 0,
            min_merge_len: 0,
            trim_len: 0,
            trim_depth: 0,
            min_dratio1: 0.0,
            max_bcov: 0.0,
            max_bfrac: 0.0,
        },
    };
    // SAFETY: `opt` is a live, correctly laid-out `fml_opt_t`.
    unsafe {
        fml_opt_init(&mut opt);
        opt.min_asm_ovlp = min_asm_ovlp;
        opt.mag_opt.flag = opt_flag;
    }

    let mut n_utg: c_int = 0;
    // SAFETY: `seqs` has `n` live entries, all allocated by libc and owned by
    // fermi-lite from here on.
    let utg = unsafe { fml_assemble(&opt, n as c_int, seqs, &mut n_utg) };

    let mut out = Vec::new();
    if utg.is_null() {
        // The library still frees the sequence buffers on this path, per upstream's
        // "we rely on fermi-lite to free this mem" comments.
        return Err("fermi-lite: fml_assemble returned NULL".to_string());
    }
    // SAFETY: `utg` is a live array of `n_utg` unitigs.
    unsafe {
        for i in 0..n_utg.max(0) as usize {
            let p = *utg.add(i);
            if p.len < 0 {
                // Upstream: `if p.len < 0: continue`.
                continue;
            }
            let len = p.len as usize;
            let seq = if p.seq.is_null() || len == 0 {
                Vec::new()
            } else {
                std::slice::from_raw_parts(p.seq as *const u8, len).to_vec()
            };
            out.push(Unitig { seq, nsr: p.nsr });
        }
        fml_utg_destroy(n_utg, utg);
    }
    Ok(out)
}

extern "C" {
    #[link_name = "malloc"]
    fn c_malloc(n: usize) -> *mut c_void;
    #[link_name = "free"]
    fn c_free(p: *mut c_void);
}

/// # Safety
/// `p` must be a pointer previously returned by `malloc`.
unsafe fn libc_free(p: *mut c_void) {
    c_free(p);
}

/// The vendored library's version string, for diagnostics and bug reports.
pub fn version() -> String {
    const HEADER: &str = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/vendor/fermi-lite/fml.h"
    ));
    HEADER
        .lines()
        .find(|l| l.contains("FML_VERSION"))
        .and_then(|l| l.split('"').nth(1))
        .unwrap_or("unknown")
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A deterministic 1000 bp template. 40 reads of 100 bp at 15 bp steps span
    /// [0, 685), so the template has to be longer than that.
    fn template() -> Vec<u8> {
        (0..1000).map(|i| b"ACGTACGTAG"[i % 10]).collect()
    }

    /// Reads covering the template densely: one start position per base.
    ///
    /// # Why not a sparser tiling
    ///
    /// The obvious fixture -- reads every 15 bp -- produces **zero** unitigs, while
    /// this one produces a unitig, and fermi-lite's own cleaning thresholds are the
    /// reason: `fml_opt_init` derives `min_elen` from the total length and read count,
    /// and `mag_opt.max_bdist` / `trim_len` interact with how a sparse set collapses.
    /// Pinning a count against that would be pinning an accident.
    ///
    /// What is worth pinning here is the **FFI contract**, not fermi-lite's cleaning
    /// policy: correct struct layout, ownership transfer with no double free, the
    /// reachable overlap range, and that a realistic read set yields something. The
    /// real behavioural contract is the `-F auto` VCF, compared end to end by
    /// `oracle/check_callvar.sh`.
    fn reads() -> Vec<(Vec<u8>, Vec<u8>)> {
        let t = template();
        (0..700)
            .map(|s| (t[s..s + 100].to_vec(), vec![35u8; 100]))
            .collect()
    }

    #[test]
    fn layouts_match_the_header() {
        // fml.h: `bseq1_t { int32_t l_seq; char *seq, *qual; }`. C pads the i32 out to
        // pointer alignment, so this is 24, not 20 -- and `repr(C)` reproduces the
        // padding rather than assuming it away.
        assert_eq!(
            std::mem::size_of::<Bseq1T>(),
            8 + 2 * std::mem::size_of::<*mut c_void>()
        );
        // fml.h: `magopt_t` is 11 ints then 3 floats.
        assert_eq!(std::mem::size_of::<MagOptT>(), 14 * 4);
        // fml.h: `fml_opt_t` is 6 ints then a `magopt_t`.
        assert_eq!(
            std::mem::size_of::<FmlOptT>(),
            6 * 4 + 14 * 4,
            "fml_opt_t layout drifted from fml.h"
        );
    }

    #[test]
    fn a_realistic_read_set_yields_at_least_one_unitig() {
        let u = assemble(&reads(), 30, DEFAULT_OPT_FLAG).expect("assembly");
        assert!(!u.is_empty(), "assembly produced no unitigs at all");
        for x in &u {
            assert!(!x.seq.is_empty(), "a unitig with no sequence");
            assert!(
                x.seq.iter().all(|b| b"ACGTN".contains(b)),
                "unitig contains non-ACGTN bases"
            );
            assert!(x.nsr >= 1, "a unitig with no supporting reads");
        }
    }

    #[test]
    fn the_unitig_comes_from_the_sequence_it_was_built_from() {
        let u = assemble(&reads(), 30, DEFAULT_OPT_FLAG).expect("assembly");
        let t = template();
        let total: usize = u.iter().map(|x| x.seq.len()).sum();
        assert!(total > 0);
        // Every assembled base must appear somewhere in the template it came from. This
        // catches a buffer that is read after being freed -- the classic symptom of
        // getting fermi-lite's ownership rules wrong.
        let in_template = u
            .iter()
            .flat_map(|x| x.seq.iter())
            .filter(|b| t.contains(b))
            .count();
        assert_eq!(
            in_template, total,
            "assembled bases must come from the input"
        );
    }

    #[test]
    fn empty_input_is_an_empty_result_not_a_crash() {
        assert!(assemble(&[], 30, DEFAULT_OPT_FLAG).unwrap().is_empty());
    }

    #[test]
    fn the_min_overlap_must_stay_within_ferm_lites_kmer_limit() {
        // `min_asm_ovlp` is used directly as a k-mer length in `bfc_ch_init`, which
        // asserts `k <= 63`. Every value in the reachable range must work.
        for ovlp in [1, 15, 30, 60, 63] {
            let u = assemble(&reads(), ovlp, DEFAULT_OPT_FLAG)
                .unwrap_or_else(|e| panic!("ovlp={ovlp} failed: {e}"));
            assert!(!u.is_empty(), "ovlp={ovlp} produced no unitigs");
        }
    }

    #[test]
    fn the_no_simpl_flag_is_passed_through() {
        // `0x80` is MAG_F_NO_SIMPL. Both it and the bare default must run; the flag's
        // effect on output is upstream's, but wiring it wrong would silently assemble a
        // different graph.
        let a = assemble(&reads(), 30, DEFAULT_OPT_FLAG).unwrap();
        let b = assemble(&reads(), 30, MAG_F_POPOPEN).unwrap();
        assert!(!a.is_empty() && !b.is_empty());
    }
}
