//! Single-end storage: sorted per-strand position arrays.

use super::SortedPositions;

/// Duplicate filtering, matching `FWTrack.filter_dup`.
///
/// Upstream's structure, kept because its quirks are observable:
///
/// ```text
/// for each chromosome:
///     if len(plus) <= 1: new_plus = plus        # no filtering, no counting
///     else: run-length truncate to maxnum, total += i_new
///     ... same for minus
/// ```
///
/// Note the `total` accumulation lives *inside* the `else`, so a chromosome with
/// zero or one position on a strand contributes nothing to the returned total even
/// though its position is retained. That asymmetry is upstream's and is
/// reproduced exactly; it is why `total` after `filter_dup(1)` can be smaller than
/// the number of positions still present.
impl SortedPositions {
    /// Filter duplicates, returning upstream's reported total.
    pub(crate) fn filter_dup(&mut self, maxnum: i64) -> u64 {
        let mut total: u64 = 0;
        let n = self.plus.len();
        for idx in 0..n {
            // `len(plus) <= 1` bypasses the whole body, including `total`
            if self.plus[idx].len() > 1 {
                truncate_run(&mut self.plus[idx], maxnum);
                total += self.plus[idx].len() as u64;
            }
            if self.minus[idx].len() > 1 {
                truncate_run(&mut self.minus[idx], maxnum);
                total += self.minus[idx].len() as u64;
            }
        }
        self.total = total;
        total
    }
}

/// Keep at most `maxnum` entries from each run of equal positions.
///
/// The input must be sorted ascending, which `finalize` guarantees. Because equal
/// values are adjacent after sorting, a run of length `k` is truncated to
/// `min(k, maxnum)`.
fn truncate_run(positions: &mut Vec<macs_core::Coord>, maxnum: i64) {
    if positions.len() <= 1 {
        return;
    }
    // `maxnum` is a signed count; upstream only ever calls this with a
    // non-negative value, and a negative one already returned early.
    let cap = if maxnum < 0 { 0 } else { maxnum as usize };
    if cap == 0 {
        // upstream keeps the first entry (`new_plus[0] = plus[0]`) before the run
        // loop, so a maxnum of 0 leaves exactly one position per chromosome
        positions.truncate(1);
        return;
    }
    let mut write = 0usize;
    let mut run = 1usize;
    positions[write] = positions[0];
    write += 1;
    for read in 1..positions.len() {
        let p = positions[read];
        if p == positions[read - 1] && read > 0 {
            run += 1;
        } else {
            run = 1;
        }
        if run <= cap {
            positions[write] = p;
            write += 1;
        }
    }
    positions.truncate(write);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_run_is_truncated_to_its_prefix() {
        let mut v = vec![7, 7, 7, 9, 9];
        truncate_run(&mut v, 2);
        assert_eq!(v, vec![7, 7, 9, 9]);
    }

    #[test]
    fn maxnum_zero_leaves_one_position() {
        // upstream writes plus[0] before entering the run loop
        let mut v = vec![7, 7, 7];
        truncate_run(&mut v, 0);
        assert_eq!(v, vec![7]);
    }

    #[test]
    fn a_run_shorter_than_maxnum_is_untouched() {
        let mut v = vec![1, 2, 3];
        truncate_run(&mut v, 5);
        assert_eq!(v, vec![1, 2, 3]);
    }

    #[test]
    fn all_positions_distinct_is_a_no_op() {
        let mut v = vec![1, 2, 3, 4];
        truncate_run(&mut v, 1);
        assert_eq!(v, vec![1, 2, 3, 4]);
    }

    #[test]
    fn many_runs_of_varying_length() {
        let mut v = vec![1, 1, 1, 2, 3, 3, 4, 4, 4, 4, 5];
        truncate_run(&mut v, 2);
        assert_eq!(v, vec![1, 1, 2, 3, 3, 4, 4, 5]);
    }
}
