//! Dump post-dedup single-end positions, for `oracle/check_postdedup.py`.
use macs_core::Strand;
use macs_track::SingleEndTrackBuilder;
use std::io::{BufRead, BufReader};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let path = a.get(1).expect("bed path");
    let chrom = a.get(2).cloned().unwrap_or_else(|| "chrM".into());
    let maxnum: i64 = a.get(3).map(|s| s.parse().expect("maxnum")).unwrap_or(1);
    let f = std::fs::File::open(path).expect("open");
    let mut r = BufReader::new(f);
    let mut line = Vec::new();
    let mut b = SingleEndTrackBuilder::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line).unwrap() == 0 {
            break;
        }
        let Some(rec) = macs_io::parse_bed_line(&line).unwrap() else {
            continue;
        };
        if rec.pos < 0 || rec.chrom.is_empty() {
            continue;
        }
        b.push(&rec.chrom, rec.pos as u64, rec.strand);
    }
    b.finalize();
    let mut t = b.build();
    eprintln!("total={}", t.filter_dup(maxnum).unwrap());
    let p = t.positions();
    let Some(c) = p.genome().get(chrom.as_bytes()) else {
        return;
    };
    for x in p.strand(c, Strand::Plus) {
        println!("+:{x}");
    }
    for x in p.strand(c, Strand::Minus) {
        println!("-:{x}");
    }
}
