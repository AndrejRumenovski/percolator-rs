use percolator_rs::stats::{
    qvalues, target_count_at_fdr_into, target_mask_at_fdr_into, QValueWorkspace, Tdc,
};

/// Selection must agree with independently materialized, reverse-minimum
/// q-values, including strict threshold edges and unusual public TDC settings.
#[test]
fn count_and_mask_match_qvalues_at_exact_threshold_edges() {
    let mut settings = vec![Tdc::default()];
    for probability in [
        -1.0,
        0.0,
        f64::from_bits(1),
        0.25,
        0.5,
        0.9,
        1.0f64.next_down(),
        1.0,
        f64::INFINITY,
        f64::NAN,
    ] {
        settings.push(Tdc::training(probability));
        settings.push(Tdc::reported(probability));
    }
    for pi0 in [
        f64::NEG_INFINITY,
        -1.0,
        -0.0,
        0.0,
        0.25,
        2.0,
        f64::MAX,
        f64::INFINITY,
        f64::NAN,
    ] {
        for tdc in [Tdc::reported(0.25), Tdc::training(0.25)] {
            settings.push(Tdc { pi0, ..tdc });
        }
    }

    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut workspace = QValueWorkspace::default();
    let mut mask = vec![7; 400];
    let mut order = vec![usize::MAX; 400];
    for n in [257, 0, 129, 1, 64, 2, 31, 3, 17] {
        for pattern in 0..3 {
            let mut scores = Vec::with_capacity(n);
            let mut labels = Vec::with_capacity(n);
            for row in 0..n {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                scores.push(match (pattern, row % 11) {
                    (0, _) => 1.0,
                    (_, 0) => -0.0,
                    (_, 1) => 0.0,
                    (_, 2) => f64::MAX,
                    (_, 3) => -f64::MAX,
                    (_, 4) => f64::from_bits(1),
                    (_, 5) => -f64::from_bits(1),
                    (1, _) => ((state >> 32) % 7) as f64 - 3.0,
                    _ => f64::from_bits(state & 0xffef_ffff_ffff_ffff),
                });
                // Exercise the full i8 domain: only positivity identifies a
                // target; zero and every negative value identify a decoy.
                labels.push((state >> 56) as i8);
            }
            let mut expected_order: Vec<_> = (0..n).collect();
            expected_order.sort_unstable_by(|&a, &b| scores[b].total_cmp(&scores[a]));
            for &tdc in &settings {
                let q = qvalues(&scores, &labels, tdc);
                let mut thresholds = vec![
                    f64::NEG_INFINITY,
                    -1.0,
                    -0.0,
                    0.0,
                    f64::MIN_POSITIVE,
                    0.01,
                    0.5,
                    1.0,
                    1.01,
                    f64::INFINITY,
                    f64::NAN,
                ];
                for &value in &q {
                    thresholds.extend([value.next_down(), value, value.next_up()]);
                }
                thresholds.sort_unstable_by(f64::total_cmp);
                thresholds.dedup_by(|a, b| a.to_bits() == b.to_bits());
                for threshold in thresholds {
                    let expected_mask: Vec<_> = q
                        .iter()
                        .zip(&labels)
                        .map(|(&value, &label)| u8::from(label > 0 && value < threshold))
                        .collect();
                    let count =
                        target_count_at_fdr_into(&scores, &labels, tdc, threshold, &mut order);
                    target_mask_at_fdr_into(
                        &scores,
                        &labels,
                        tdc,
                        threshold,
                        &mut workspace,
                        &mut mask,
                    );
                    assert_eq!(
                        count,
                        expected_mask
                            .iter()
                            .map(|&value| value as usize)
                            .sum::<usize>(),
                        "n={n}, pattern={pattern}, tdc={tdc:?}, threshold={threshold:?}"
                    );
                    assert_eq!(
                        mask, expected_mask,
                        "n={n}, pattern={pattern}, tdc={tdc:?}, threshold={threshold:?}"
                    );
                    assert_eq!(order, expected_order);
                }
            }
        }
    }
}

#[test]
fn training_mask_rejects_non_finite_scores_before_restricting_rows() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for label in [1, -1] {
            // These training settings permit cutoff selection. Validate every
            // score before selection could discard a non-finite target or decoy.
            let scores = [3.0, 2.0, 1.0, bad];
            let labels = [1, -1, 1, label];
            assert!(
                std::panic::catch_unwind(|| target_mask_at_fdr_into(
                    &scores,
                    &labels,
                    Tdc::training(0.5),
                    0.01,
                    &mut QValueWorkspace::default(),
                    &mut Vec::new(),
                ))
                .is_err(),
                "bad={bad:?}, label={label}"
            );
        }
    }
}

#[test]
fn training_mask_preserves_large_tied_prefixes_and_decoy_extremes() {
    // An odd affine permutation shuffles all 4096 ranks. Each score group has
    // 16 members; the zero group contains both signs and mixed labels.
    let ranks: Vec<_> = (0..4096).map(|row| (row * 109 + 37) % 4096).collect();
    let scores: Vec<_> = ranks
        .iter()
        .map(|&rank| {
            let score = 128.0 - (rank / 16) as f64;
            if score == 0.0 && rank % 2 == 0 {
                -0.0
            } else {
                score
            }
        })
        .collect();
    let mut workspace = QValueWorkspace::default();
    let mut mask = vec![7; 5000];
    for pattern in 0..3 {
        let labels: Vec<_> = ranks
            .iter()
            .map(|&rank| match pattern {
                0 => {
                    if rank < 256 || rank % 3 != 0 {
                        1
                    } else {
                        -1
                    }
                }
                1 => 1,
                _ => -1,
            })
            .collect();
        let targets = labels.iter().filter(|&&label| label > 0).count();
        let zero_decoys = scores
            .iter()
            .zip(&labels)
            .filter(|(score, label)| **score >= 0.0 && **label <= 0)
            .count();
        for probability in [0.5, 0.25, 1e-9, f64::MIN_POSITIVE, 1.0f64.next_down()] {
            let tdc = Tdc::training(probability);
            let q = qvalues(&scores, &labels, tdc);
            // This exact lower bound places the selected decoy inside the
            // signed-zero group for ordinary probabilities and mixed labels.
            let zero_threshold = (zero_decoys as f64 * (probability / (1.0 - probability))
                / targets.max(1) as f64)
                .min(1.0);
            let mut thresholds = vec![
                0.01,
                1.0,
                zero_threshold.next_down(),
                zero_threshold,
                zero_threshold.next_up(),
            ];
            for (row, (&score, &value)) in scores.iter().zip(&q).enumerate() {
                // The final group's q-value makes K == all decoys at p=0.5.
                if row % 509 == 0 || score == 0.0 || score == -127.0 {
                    thresholds.extend([value.next_down(), value, value.next_up()]);
                }
            }
            thresholds.sort_unstable_by(f64::total_cmp);
            thresholds.dedup_by(|a, b| a.to_bits() == b.to_bits());
            for threshold in thresholds {
                target_mask_at_fdr_into(
                    &scores,
                    &labels,
                    tdc,
                    threshold,
                    &mut workspace,
                    &mut mask,
                );
                let expected: Vec<_> = q
                    .iter()
                    .zip(&labels)
                    .map(|(&value, &label)| u8::from(label > 0 && value < threshold))
                    .collect();
                assert_eq!(
                    mask, expected,
                    "pattern={pattern}, probability={probability:?}, threshold={threshold:?}"
                );
            }
        }
    }
}
