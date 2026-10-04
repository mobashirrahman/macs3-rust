//! Read strand.

/// Strand of a read or fragment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    /// Forward / plus.
    Plus,
    /// Reverse / minus.
    Minus,
    /// Strand was not recorded (BED input without a strand column).
    Unknown,
}

impl Strand {
    /// Parse a strand character. `b'+'`, `b'-'` and `b'.'`/`b'?'` are accepted;
    /// anything else is `None` so the caller can raise a typed parse error.
    #[inline]
    pub fn from_byte(b: u8) -> Option<Strand> {
        match b {
            b'+' => Some(Strand::Plus),
            b'-' => Some(Strand::Minus),
            b'.' | b'?' => Some(Strand::Unknown),
            _ => None,
        }
    }

    /// The canonical character for this strand.
    #[inline]
    pub fn as_byte(&self) -> u8 {
        match self {
            Strand::Plus => b'+',
            Strand::Minus => b'-',
            Strand::Unknown => b'.',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        for s in [Strand::Plus, Strand::Minus, Strand::Unknown] {
            assert_eq!(Strand::from_byte(s.as_byte()), Some(s));
        }
        assert_eq!(Strand::from_byte(b'?'), Some(Strand::Unknown));
        assert_eq!(Strand::from_byte(b'x'), None);
    }
}
