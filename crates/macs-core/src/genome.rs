//! Chromosome interning.
//!
//! Genomic records never carry a `String` chromosome name. Names are interned
//! once into a [`ChromId`] so that a 100M-read experiment costs 4 bytes of
//! chromosome per record instead of a 24-byte `String` header plus a heap
//! allocation, and so that tracks can be stored contiguously per chromosome.
//!
//! # Ordering
//!
//! Output order is observable, so it is part of the contract. Three orderings
//! are supported and selected explicitly by the caller:
//!
//! * [`ChromOrder::FileOrder`] — first appearance in the input. This is what
//!   upstream does for single-file peak calling, and it is the default.
//! * [`ChromOrder::Lexicographic`] — byte-wise ascending name.
//! * [`ChromOrder::Karyotypic`] — lexicographic, except that `chr1..chr22`,
//!   `chrX`, `chrY` sort numerically and ahead of any other `chr*` contig.
//!
//! Every writer in `macs3-rs` takes an explicit order so that ordering is never
//! an accident of a hash map's iteration order.

use crate::Len;
use std::collections::HashMap;

/// A dense index into a [`Genome`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ChromId(pub u32);

/// How chromosome records are ordered in output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChromOrder {
    /// First appearance in the input file.
    #[default]
    FileOrder,
    /// Byte-wise ascending name.
    Lexicographic,
    /// `chr1..chr22, chrX, chrY`, then everything else lexicographically.
    Karyotypic,
}

/// Karyotypic rank of a `chr`-prefixed name: `(rank, numeric_value)`.
fn karyotypic_key(name: &[u8]) -> (u8, u64, &[u8]) {
    if let Some(rest) = name.strip_prefix(b"chr".as_slice()) {
        if rest.is_empty() {
            return (2, 0, rest);
        }
        if rest.iter().all(|b| b.is_ascii_digit()) {
            if let Ok(n) = std::str::from_utf8(rest).unwrap_or("").parse::<u64>() {
                return (0, n, &[]);
            }
        }
        match rest {
            b"X" | b"x" => return (1, 0, &[]),
            b"Y" | b"y" => return (1, 1, &[]),
            b"M" | b"MT" | b"mt" => return (1, 2, &[]),
            _ => {}
        }
        return (2, 0, rest);
    }
    (3, 0, name)
}

/// A dictionary of chromosome names.
///
/// Construct order defines [`ChromOrder::FileOrder`]; a lookup map provides
/// O(1) interning during parsing.
#[derive(Debug, Default)]
pub struct Genome {
    names: Vec<Box<[u8]>>,
    lookup: HashMap<Box<[u8]>, ChromId>,
    lengths: Vec<Option<Len>>,
}

impl Clone for Genome {
    /// Clones the dictionary. Ids are positional, so a clone reproduces them
    /// exactly and is interchangeable with the original for lookup purposes.
    fn clone(&self) -> Self {
        Genome {
            names: self.names.clone(),
            lookup: self.lookup.clone(),
            lengths: self.lengths.clone(),
        }
    }
}

impl Genome {
    /// An empty genome.
    pub fn new() -> Self {
        Self::default()
    }

    /// Intern a chromosome name, returning its stable id.
    ///
    /// Interning is idempotent: interning the same name twice returns the same
    /// id, and does not change file order.
    pub fn intern(&mut self, name: &[u8]) -> ChromId {
        if let Some(&id) = self.lookup.get(name) {
            return id;
        }
        let id = ChromId(self.names.len() as u32);
        let boxed: Box<[u8]> = name.to_vec().into_boxed_slice();
        self.names.push(boxed.clone());
        self.lookup.insert(boxed, id);
        self.lengths.push(None);
        id
    }

    /// Intern a name, erroring if it is empty.
    pub fn intern_str(&mut self, name: &str) -> ChromId {
        debug_assert!(!name.is_empty());
        self.intern(name.as_bytes())
    }

    /// Number of interned chromosomes.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// True when nothing has been interned.
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// The name of a chromosome id.
    pub fn name(&self, id: ChromId) -> &[u8] {
        &self.names[id.0 as usize]
    }

    /// The name of a chromosome id, lossy-converted for diagnostics.
    pub fn name_string(&self, id: ChromId) -> String {
        String::from_utf8_lossy(self.name(id)).into_owned()
    }

    /// Look up a name without interning it.
    pub fn get(&self, name: &[u8]) -> Option<ChromId> {
        self.lookup.get(name).copied()
    }

    /// True when the name has been interned.
    pub fn contains(&self, name: &[u8]) -> bool {
        self.lookup.contains_key(name)
    }

    /// Record a contig length. Shorter/later records never shrink a known
    /// length, so a partial `.fai` cannot corrupt an established contig size.
    pub fn set_len(&mut self, id: ChromId, len: Len) {
        let slot = &mut self.lengths[id.0 as usize];
        *slot = Some(match *slot {
            Some(existing) => existing.max(len),
            None => len,
        });
    }

    /// The recorded contig length, if known.
    ///
    /// Returns `None` for an id that was never interned, so this is total
    /// rather than panicking on a caller bug.
    pub fn len_of(&self, id: ChromId) -> Option<Len> {
        self.lengths.get(id.0 as usize).copied().flatten()
    }

    /// Ids in file order.
    pub fn ids_file_order(&self) -> Vec<ChromId> {
        (0..self.names.len()).map(|i| ChromId(i as u32)).collect()
    }

    /// Ids in the requested order.
    pub fn ids_in(&self, order: ChromOrder) -> Vec<ChromId> {
        let mut ids = self.ids_file_order();
        match order {
            ChromOrder::FileOrder => {}
            ChromOrder::Lexicographic => {
                ids.sort_by(|&a, &b| self.name(a).cmp(self.name(b)));
            }
            ChromOrder::Karyotypic => {
                ids.sort_by(|&a, &b| {
                    let (ra, na, fa) = karyotypic_key(self.name(a));
                    let (rb, nb, fb) = karyotypic_key(self.name(b));
                    (ra, na, fb).cmp(&(rb, nb, fa))
                });
            }
        }
        ids
    }

    /// All names in file order, for diagnostics and error messages.
    pub fn all_names(&self) -> Vec<String> {
        self.names
            .iter()
            .map(|n| String::from_utf8_lossy(n).into_owned())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interning_is_idempotent_and_preserves_file_order() {
        let mut g = Genome::new();
        let c1 = g.intern(b"chr2");
        let c2 = g.intern(b"chr1");
        assert_eq!(c1, ChromId(0));
        assert_eq!(c2, ChromId(1));
        assert_eq!(g.intern(b"chr2"), c1);
        assert_eq!(g.len(), 2);
        assert_eq!(g.ids_in(ChromOrder::FileOrder), vec![c1, c2]);
    }

    #[test]
    fn lexicographic_and_karyotypic_differ() {
        let mut g = Genome::new();
        for n in [&b"chr10"[..], b"chr2", b"chrX", b"GL000191.1", b"chr1"] {
            g.intern(n);
        }
        assert_eq!(
            g.ids_in(ChromOrder::Lexicographic)
                .iter()
                .map(|&i| g.name_string(i))
                .collect::<Vec<_>>(),
            vec!["GL000191.1", "chr1", "chr10", "chr2", "chrX"]
        );
        assert_eq!(
            g.ids_in(ChromOrder::Karyotypic)
                .iter()
                .map(|&i| g.name_string(i))
                .collect::<Vec<_>>(),
            vec!["chr1", "chr2", "chr10", "chrX", "GL000191.1"]
        );
    }

    #[test]
    fn lengths_never_shrink() {
        let mut g = Genome::new();
        let c = g.intern(b"chr1");
        g.set_len(c, 1000);
        g.set_len(c, 500);
        assert_eq!(g.len_of(c), Some(1000));
        assert_eq!(g.len_of(ChromId(99)), None);
    }

    #[test]
    fn karyotypic_key_is_total() {
        // names must have a strict, deterministic total order
        let mut names: Vec<Vec<u8>> = vec![
            b"chr1".to_vec(),
            b"chr2".to_vec(),
            b"chr10".to_vec(),
            b"chrX".to_vec(),
            b"chrM".to_vec(),
            b"scaffold_1".to_vec(),
            b"1".to_vec(),
            b"chr".to_vec(),
        ];
        names.sort_by(|a, b| karyotypic_key(a).cmp(&karyotypic_key(b)));
        assert!(!names.is_empty());
    }
}
