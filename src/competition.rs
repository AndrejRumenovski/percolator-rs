//! Spectrum-level target-decoy competition policy.

use crate::{pin, tiebreak};

/// Spectrum-level target-decoy competition on the rescored values: keep the
/// best-scoring candidate of each precursor and drop the rest.
///
/// A search that reports its top N matches per spectrum hands over N candidates
/// that are not independent hypotheses. They come from one measurement, they
/// compete with each other, and several of them can be accepted together. The
/// target-decoy estimator assumes each spectrum contributes at most the winner
/// of a competition against the decoy database -- that assumption is what makes
/// one observed decoy stand for one incorrect target. Keeping every candidate
/// breaks it, and the reported q-value stops being an FDR estimate.
///
/// This is the reference's `--post-processing-tdc`
/// (`Scores::weedOutRedundantTDC`), which likewise keeps one row per
/// (scan, mass, charge) and erases the rest, and it runs after training on the
/// same rescored values. It is on by default here because percolator-rs makes a
/// calibration claim about its q-values, and that claim is only available on
/// competed input.
///
/// # Exact ties
///
/// A precursor whose best score is attained by more than one candidate has no
/// winner on the evidence. Choosing the earlier or the later row would make the
/// surviving label a property of the file's layout: a PIN listing the target
/// candidate first would report every target as a winner, and the same PIN with
/// the two rows swapped would report none. Ties are therefore drawn with a fair
/// coin keyed on the precursor's own identity and the run seed
/// ([`tiebreak`]), which is the resolution the target-decoy literature
/// prescribes and the only one that keeps the null win probability at the
/// declared `p`. Permuting the rows of the input cannot move the coin.
pub fn winner_indices(ds: &pin::Dataset, score: &[f64], seed: u64) -> Vec<usize> {
    #[derive(Clone, Copy)]
    struct Best {
        score: f64,
        tied: u32,
        row: usize,
    }

    let mut best: ahash::AHashMap<(u32, i64, u64), Best> =
        ahash::AHashMap::with_capacity(ds.n_psm / 2 + 1);
    for i in 0..ds.n_psm {
        let key = ds.spectrum_key(i);
        match best.entry(key) {
            std::collections::hash_map::Entry::Vacant(slot) => {
                slot.insert(Best {
                    score: score[i],
                    tied: 1,
                    row: i,
                });
            }
            std::collections::hash_map::Entry::Occupied(mut slot) => {
                let current = slot.get_mut();
                if score[i] > current.score {
                    *current = Best {
                        score: score[i],
                        tied: 1,
                        row: i,
                    };
                } else if score[i] == current.score {
                    current.tied += 1;
                }
            }
        }
    }

    // Exact ties are rare, so only the tied rows are materialized and ordered.
    // `canonical` describes a row by its own content, never by its position, so
    // the sort below is a function of the data alone.
    let canonical = |row: usize| -> (i8, &str, &str, &str, usize) {
        (
            ds.labels[row],
            ds.peptide[row].as_str(),
            ds.proteins[row].as_str(),
            ds.spec_id[row].as_str(),
            // Reached only between rows that agree in every identifying field
            // and are therefore interchangeable.
            row,
        )
    };
    let mut contested: Vec<((u32, i64, u64), usize)> = Vec::new();
    for i in 0..ds.n_psm {
        let key = ds.spectrum_key(i);
        let entry = best[&key];
        if entry.tied > 1 && score[i] == entry.score {
            contested.push((key, i));
        }
    }
    contested.sort_unstable_by(|(left_key, left), (right_key, right)| {
        left_key
            .cmp(right_key)
            .then_with(|| canonical(*left).cmp(&canonical(*right)))
    });

    let mut winners: Vec<usize> = Vec::with_capacity(best.len());
    for entry in best.values() {
        if entry.tied == 1 {
            winners.push(entry.row);
        }
    }
    let mut start = 0usize;
    while start < contested.len() {
        let key = contested[start].0;
        let mut end = start;
        while end + 1 < contested.len() && contested[end + 1].0 == key {
            end += 1;
        }
        let group = &contested[start..=end];
        let draw = tiebreak::Coin::new(seed)
            .u32(key.0)
            .i64(key.1)
            .u64(key.2)
            .draw(group.len());
        winners.push(group[draw].1);
        start = end + 1;
    }
    winners.sort_unstable();
    winners
}

#[cfg(test)]
mod competition_tests {
    use super::*;

    fn dataset(rows: &[(u32, i64, f64, i8)]) -> pin::Dataset {
        pin::Dataset {
            feature_names: vec!["f".to_string()],
            n_feat: 1,
            n_psm: rows.len(),
            features: vec![0.0; rows.len()],
            labels: rows.iter().map(|r| r.3).collect(),
            spec_id: (0..rows.len()).map(|i| format!("s{i}")).collect(),
            scan: rows.iter().map(|r| r.1).collect(),
            exp_mass: rows.iter().map(|r| r.2).collect(),
            peptide: (0..rows.len()).map(|i| format!("K.P{i}.R")).collect(),
            proteins: (0..rows.len()).map(|i| format!("P{i}")).collect(),
            source: rows.iter().map(|r| r.0).collect(),
            source_names: vec!["a.pin".to_string(), "b.pin".to_string()],
            ensemble: false,
        }
    }

    /// A dataset whose rows carry explicit peptide/protein/spectrum identity, so
    /// a permutation can be described by content rather than by position.
    struct Rows {
        rows: Vec<(u32, i64, f64, i8, String, f64)>, // source, scan, mass, label, peptide, score
    }

    impl Rows {
        fn build(&self, order: &[usize]) -> (pin::Dataset, Vec<f64>) {
            let n = order.len();
            let mut ds = pin::Dataset {
                feature_names: vec!["f".to_string()],
                n_feat: 1,
                n_psm: n,
                features: vec![0.0; n],
                labels: Vec::with_capacity(n),
                spec_id: Vec::with_capacity(n),
                scan: Vec::with_capacity(n),
                exp_mass: Vec::with_capacity(n),
                peptide: Vec::with_capacity(n),
                proteins: Vec::with_capacity(n),
                source: Vec::with_capacity(n),
                source_names: vec!["a.pin".to_string()],
                ensemble: false,
            };
            let mut score = Vec::with_capacity(n);
            for &index in order {
                let row = &self.rows[index];
                ds.source.push(row.0);
                ds.scan.push(row.1);
                ds.exp_mass.push(row.2);
                ds.labels.push(row.3);
                ds.spec_id.push(format!("scan{}_{}", row.1, row.4));
                ds.peptide.push(format!("K.{}.R", row.4));
                ds.proteins.push(row.4.clone());
                score.push(row.5);
            }
            (ds, score)
        }

        /// The winner set described by content, so it can be compared across
        /// permutations that assign different row indices to the same PSM.
        fn winner_identities(&self, order: &[usize], seed: u64) -> Vec<(i64, String)> {
            let (ds, score) = self.build(order);
            let mut identities: Vec<(i64, String)> = winner_indices(&ds, &score, seed)
                .into_iter()
                .map(|row| (ds.scan[row], ds.peptide[row].clone()))
                .collect();
            identities.sort();
            identities
        }
    }

    /// One target and one decoy candidate per spectrum, scoring exactly the same.
    fn tied_pairs(spectra: i64) -> Rows {
        let mut rows = Vec::new();
        for scan in 1..=spectra {
            rows.push((0u32, scan, 500.0, 1i8, format!("TARGET{scan}"), 7.5));
            rows.push((0u32, scan, 500.0, -1i8, format!("DECOY{scan}"), 7.5));
        }
        Rows { rows }
    }

    fn permutation(n: usize, seed: u64) -> Vec<usize> {
        let mut order: Vec<usize> = (0..n).collect();
        let mut state = seed.max(1);
        for index in (1..n).rev() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            order.swap(index, (state % (index as u64 + 1)) as usize);
        }
        order
    }

    /// **Frozen failure case (independent audit C1).**
    ///
    /// 200 spectra, each with one target and one decoy candidate carrying the
    /// exact same score.  The multiset of scores, labels, features and peptides
    /// is identical in every arm; only the order of the rows in the file
    /// changes.  Any competition rule that reads row order turns that metadata
    /// choice into a different scientific result.
    #[test]
    fn exact_ties_survive_every_equivalent_row_permutation() {
        let fixture = tied_pairs(200);
        let n = fixture.rows.len();

        let target_first: Vec<usize> = (0..n).collect();
        let reference = fixture.winner_identities(&target_first, 1);

        // Decoy row before target row, within every pair.
        let pair_reversed: Vec<usize> = (0..n)
            .map(|i| if i % 2 == 0 { i + 1 } else { i - 1 })
            .collect();
        assert_eq!(
            fixture.winner_identities(&pair_reversed, 1),
            reference,
            "reversing each tied target/decoy pair changed the winners"
        );

        // Whole file reversed.
        let reversed: Vec<usize> = (0..n).rev().collect();
        assert_eq!(
            fixture.winner_identities(&reversed, 1),
            reference,
            "reversing the file changed the winners"
        );

        // All targets first, then all decoys, and the other way round.
        let grouped: Vec<usize> = (0..n).step_by(2).chain((1..n).step_by(2)).collect();
        assert_eq!(
            fixture.winner_identities(&grouped, 1),
            reference,
            "grouping targets before decoys changed the winners"
        );
        let grouped_decoys: Vec<usize> = (1..n).step_by(2).chain((0..n).step_by(2)).collect();
        assert_eq!(
            fixture.winner_identities(&grouped_decoys, 1),
            reference,
            "grouping decoys before targets changed the winners"
        );

        // Deterministic shuffles.
        for shuffle_seed in 1..=8u64 {
            assert_eq!(
                fixture.winner_identities(&permutation(n, shuffle_seed), 1),
                reference,
                "shuffle {shuffle_seed} changed the winners"
            );
        }
    }

    /// **Frozen failure case (independent audit C1).**
    ///
    /// The statistical consequence of the ordering attack: whichever label the
    /// rule prefers, it wins every tie, and 200 spectra of pure noise become
    /// 200 confident discoveries.  A fair rule splits the ties, so neither label
    /// can dominate.
    #[test]
    fn exact_ties_do_not_hand_every_precursor_to_one_label() {
        let fixture = tied_pairs(200);
        let order: Vec<usize> = (0..fixture.rows.len()).collect();
        let (ds, score) = fixture.build(&order);
        let winners = winner_indices(&ds, &score, 1);
        assert_eq!(winners.len(), 200, "one winner per spectrum");
        let targets = winners.iter().filter(|&&row| ds.labels[row] > 0).count();
        // A fair coin over 200 two-way ties has SD 7.07; 5 SD is 35.
        assert!(
            (targets as i64 - 100).abs() < 35,
            "{targets} of 200 tied precursors were won by the target label"
        );
    }

    /// Changing the seed must re-flip the coins, so tie sensitivity shows up as
    /// seed variability instead of hiding inside a fixed row order.
    #[test]
    fn a_different_seed_resolves_ties_differently() {
        let fixture = tied_pairs(200);
        let order: Vec<usize> = (0..fixture.rows.len()).collect();
        let first = fixture.winner_identities(&order, 1);
        let second = fixture.winner_identities(&order, 2);
        assert_ne!(
            first, second,
            "two seeds resolved 200 exact ties identically"
        );
    }

    /// Ties among more than two candidates must also be label-fair: with three
    /// targets and one decoy tied, the decoy wins about a quarter of the time.
    #[test]
    fn wider_ties_are_drawn_uniformly() {
        let mut rows = Vec::new();
        for scan in 1..=2_000i64 {
            for candidate in 0..3 {
                rows.push((
                    0u32,
                    scan,
                    500.0,
                    1i8,
                    format!("TARGET{scan}_{candidate}"),
                    7.5,
                ));
            }
            rows.push((0u32, scan, 500.0, -1i8, format!("DECOY{scan}"), 7.5));
        }
        let fixture = Rows { rows };
        let order: Vec<usize> = (0..fixture.rows.len()).collect();
        let (ds, score) = fixture.build(&order);
        let winners = winner_indices(&ds, &score, 1);
        assert_eq!(winners.len(), 2_000);
        let decoys = winners.iter().filter(|&&row| ds.labels[row] < 0).count();
        // Expectation 500, SD 19.4; 5 SD is 97.
        assert!(
            (decoys as i64 - 500).abs() < 97,
            "{decoys} of 2000 four-way ties were won by the decoy"
        );
    }

    /// A strict winner is never displaced by the tie rule.
    #[test]
    fn a_strictly_better_candidate_always_wins() {
        for seed in 1..=16u64 {
            let ds = dataset(&[(0, 1, 500.0, 1), (0, 1, 500.0, -1), (0, 1, 500.0, 1)]);
            assert_eq!(winner_indices(&ds, &[3.0, 4.0, 3.0], seed), vec![1]);
            assert_eq!(winner_indices(&ds, &[5.0, 4.0, 3.0], seed), vec![0]);
        }
    }

    /// Exact equality, not an epsilon, defines a score tie. A one-ULP score
    /// advantage is evidence and must win independently of the lottery seed.
    #[test]
    fn a_one_ulp_near_tie_is_never_randomized() {
        let ds = dataset(&[(0, 1, 500.0, 1), (0, 1, 500.0, -1)]);
        let low = 7.5f64;
        let high = f64::from_bits(low.to_bits() + 1);
        for seed in 1..=64u64 {
            assert_eq!(winner_indices(&ds, &[low, high], seed), vec![1]);
        }
    }

    #[test]
    fn one_winner_per_precursor() {
        // Two spectra, five candidates each, mixed targets and decoys.
        let rows: Vec<(u32, i64, f64, i8)> = (0..10)
            .map(|i| {
                (
                    0u32,
                    (i / 5) as i64,
                    800.0,
                    if i % 2 == 0 { 1i8 } else { -1 },
                )
            })
            .collect();
        let ds = dataset(&rows);
        let score: Vec<f64> = vec![1.0, 5.0, 2.0, 3.0, 4.0, 9.0, 8.0, 7.0, 6.0, 5.5];
        let winners = winner_indices(&ds, &score, 1);
        assert_eq!(winners, vec![1, 5]);
    }

    /// Distinct precursors from one scan must not be collapsed: their
    /// experimental neutral masses differ.
    #[test]
    fn distinct_precursors_of_one_scan_compete_separately() {
        let ds = dataset(&[
            (0, 42, 800.0, 1),
            (0, 42, 800.0, -1),
            (0, 42, 1600.0, 1),
            (0, 42, 1600.0, -1),
        ]);
        let winners = winner_indices(&ds, &[1.0, 2.0, 5.0, 4.0], 1);
        assert_eq!(winners, vec![1, 2]);
    }

    #[test]
    fn joined_files_reusing_scan_numbers_compete_separately() {
        let ds = dataset(&[(0, 7, 900.0, 1), (1, 7, 900.0, 1)]);
        assert_eq!(winner_indices(&ds, &[1.0, 2.0], 1), vec![0, 1]);
    }

    /// The tie draw must be reproducible: repeated calls with the same seed and
    /// the same content give the same winner.
    #[test]
    fn tie_resolution_is_reproducible() {
        let ds = dataset(&[(0, 1, 500.0, 1), (0, 1, 500.0, -1), (0, 1, 500.0, 1)]);
        let winners = winner_indices(&ds, &[3.0, 3.0, 3.0], 1);
        assert_eq!(winners.len(), 1);
        for _ in 0..16 {
            assert_eq!(winner_indices(&ds, &[3.0, 3.0, 3.0], 1), winners);
        }
    }

    /// Input that already holds one candidate per spectrum must pass through
    /// untouched, so competition cannot silently drop rows from a competed file.
    #[test]
    fn already_competed_input_is_unchanged() {
        let rows: Vec<(u32, i64, f64, i8)> = (0..6)
            .map(|i| {
                (
                    0u32,
                    i as i64,
                    700.0 + i as f64,
                    if i % 2 == 0 { 1i8 } else { -1 },
                )
            })
            .collect();
        let ds = dataset(&rows);
        let score = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(winner_indices(&ds, &score, 1), (0..6).collect::<Vec<_>>());
    }

    /// A PIN without an ExpMass column reports 0.0 for every row, which must
    /// collapse a scan to one precursor rather than split it by float noise.
    #[test]
    fn a_pin_without_expmass_competes_per_scan() {
        let ds = dataset(&[(0, 3, 0.0, 1), (0, 3, -0.0, -1), (0, 4, 0.0, 1)]);
        assert_eq!(winner_indices(&ds, &[1.0, 2.0, 0.5], 1), vec![1, 2]);
    }
}
