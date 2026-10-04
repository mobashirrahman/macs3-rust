//! Reproduce upstream's per-peak `RACollection` consensus for the differential.
//!
//! Emits the identical record `oracle/dump_callvar_refseq.py` produces from the
//! pinned oracle, so a divergence in either can be attributed to a specific field.

use macs_callvar::ra_collection_dump::{dump_peak_records, DumpOptions};
use std::path::Path;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let opt = DumpOptions {
        peaks: Path::new(&a[1]),
        tbam: Path::new(&a[2]),
        cbam: a.get(3).map(Path::new),
        max_duplicate: 1,
    };
    match dump_peak_records(&opt) {
        Ok(recs) => {
            for r in recs {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    String::from_utf8_lossy(&r.chrom),
                    r.left,
                    r.right,
                    r.read_start,
                    r.read_end,
                    r.count_t,
                    r.count_c,
                    r.ext_len,
                    r.seq_hex,
                    r.ext_hex
                );
            }
        }
        Err(e) => {
            eprintln!("dump_refseq: {e}");
            std::process::exit(1);
        }
    }
}
