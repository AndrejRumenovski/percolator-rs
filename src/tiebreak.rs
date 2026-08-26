//! Order-invariant resolution of exact score ties in target-decoy competition.
//!
//! # Why a tie needs a rule at all
//!
//! Target-decoy competition converts one observed decoy win into one expected
//! incorrect target.  That conversion is only valid if an incorrect target beats
//! its decoy counterpart with the declared probability `p` — 1/2 for a 1:1
//! concatenated search.  When two candidates of the same precursor carry the
//! *exact* same score the competition has no winner, and whatever rule supplies
//! one decides the label of the survivor.
//!
//! Resolving such a tie by row order (keep the first, keep the last, rely on a
//! stable or unstable sort) makes that decision a property of the input file's
//! layout rather than of the measurement.  A PIN that lists the target candidate
//! before the decoy then yields every target as a winner; the same PIN with the
//! two rows swapped yields none.  The reported q-values are computed correctly
//! from a winner list that is itself biased, so the error is invisible
//! downstream.
//!
//! # The rule
//!
//! Ties are broken by a **fair coin**, which is the resolution used by the
//! primary target-decoy literature (Granholm, Noble & Käll 2011, *On using
//! samples of known protein content to assess the statistical calibration of
//! scores assigned to peptide-spectrum matches in shotgun proteomics*): when a
//! target and a decoy tie, a coin decides, which is exactly what keeps the null
//! win probability at 1/2.
//!
//! Reproducibility is preserved by drawing that coin from a hash of the
//! *competition unit's own identity* and the run seed, never from the row order:
//!
//! ```text
//! winner = canonical_order(tied candidates)[ draw(hash(seed, unit identity), k) ]
//! ```
//!
//! * `hash` is a SplitMix64 mixing of the precursor identity (source, scan,
//!   experimental mass) with the run seed.  It sees no labels and no row indices.
//! * `canonical_order` sorts the tied candidates by their own content
//!   (label, peptide, proteins, spectrum id), falling back to the row index only
//!   between rows that are byte-identical in all of those and therefore
//!   interchangeable.
//! * `draw` is uniform over `0..k`, so a `k`-way tie holding `t` targets is won
//!   by a target with probability `t/k`.  With `k = 2` that is the fair coin;
//!   under a 1:1 database each tied candidate is a target with probability 1/2,
//!   so the marginal target-win probability is 1/2 for any `k`.
//!
//! Permuting the rows of a PIN therefore cannot change any winner, any q-value,
//! or any reported count.  Changing the seed re-flips every coin, which is what
//! makes tie sensitivity measurable across seeds instead of hidden.

/// SplitMix64 finalizer: a bijection with good avalanche, so nearby keys give
/// unrelated draws.
#[inline]
const fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A running hash of a competition unit's identity.
///
/// Only label-free identity is ever fed in: which file, which spectrum, which
/// precursor, or which protein pairing key.  Feeding a label in would make the
/// coin depend on the very thing it is meant to decide.
#[derive(Clone, Copy)]
pub struct Coin(u64);

impl Coin {
    /// Start from the run seed.
    #[inline]
    pub fn new(seed: u64) -> Self {
        Coin(mix64(seed ^ 0x243F_6A88_85A3_08D3))
    }

    #[inline]
    pub fn u64(self, value: u64) -> Self {
        Coin(mix64(self.0 ^ mix64(value)))
    }

    #[inline]
    pub fn i64(self, value: i64) -> Self {
        self.u64(value as u64)
    }

    #[inline]
    pub fn u32(self, value: u32) -> Self {
        self.u64(u64::from(value))
    }

    /// Fold a byte string in, length-prefixed so `("ab","c")` and `("a","bc")`
    /// cannot collide.
    #[inline]
    pub fn bytes(self, value: &[u8]) -> Self {
        let mut state = self.u64(value.len() as u64).0;
        let mut chunks = value.chunks_exact(8);
        for chunk in &mut chunks {
            let word = u64::from_le_bytes(chunk.try_into().unwrap());
            state = mix64(state ^ mix64(word));
        }
        let remainder = chunks.remainder();
        if !remainder.is_empty() {
            let mut buffer = [0u8; 8];
            buffer[..remainder.len()].copy_from_slice(remainder);
            state = mix64(state ^ mix64(u64::from_le_bytes(buffer)));
        }
        Coin(state)
    }

    /// Uniform draw in `0..k` (Lemire's multiply-shift reduction).
    ///
    /// `k` is a tie-group size, so it is far below the point where the
    /// modulo bias of this reduction is representable in `f64`.
    #[inline]
    pub fn draw(self, k: usize) -> usize {
        if k <= 1 {
            return 0;
        }
        ((u128::from(self.0) * k as u128) >> 64) as usize
    }

    /// `true` for half of all keys: the fair coin of a two-way competition.
    #[inline]
    pub fn heads(self) -> bool {
        self.draw(2) == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_draw_never_leaves_its_range() {
        for k in 1..64usize {
            for key in 0..512u64 {
                assert!(Coin::new(1).u64(key).draw(k) < k);
            }
        }
    }

    #[test]
    fn two_way_draws_are_close_to_fair() {
        // 20,000 independent precursor identities, one coin each.
        let heads = (0..20_000u64)
            .filter(|&scan| Coin::new(1).u32(0).i64(scan as i64).u64(0).heads())
            .count();
        // A fair coin has SD 70.7 here; 5 SD is 354.
        assert!(
            (heads as i64 - 10_000).abs() < 354,
            "20000 two-way draws gave {heads} heads"
        );
    }

    #[test]
    fn draws_are_uniform_across_a_wider_tie_group() {
        let mut counts = [0usize; 5];
        for scan in 0..20_000i64 {
            counts[Coin::new(1).i64(scan).draw(5)] += 1;
        }
        // SD of a 1/5 bin over 20000 draws is 56.6; 5 SD is 283.
        for (index, &count) in counts.iter().enumerate() {
            assert!(
                (count as i64 - 4_000).abs() < 283,
                "bin {index} got {count} of 20000 draws"
            );
        }
    }

    #[test]
    fn different_seeds_reflip_the_coins() {
        let agree = (0..2_000i64)
            .filter(|&scan| Coin::new(1).i64(scan).heads() == Coin::new(2).i64(scan).heads())
            .count();
        assert!(
            (agree as i64 - 1_000).abs() < 113,
            "seeds 1 and 2 agreed on {agree} of 2000 coins"
        );
    }

    #[test]
    fn the_same_key_always_gives_the_same_draw() {
        for scan in 0..1_000i64 {
            let first = Coin::new(7).u32(3).i64(scan).u64(42).draw(9);
            for _ in 0..4 {
                assert_eq!(Coin::new(7).u32(3).i64(scan).u64(42).draw(9), first);
            }
        }
    }

    #[test]
    fn byte_folding_is_length_prefixed() {
        assert_ne!(
            Coin::new(1).bytes(b"ab").bytes(b"c").0,
            Coin::new(1).bytes(b"a").bytes(b"bc").0
        );
    }
}
