#![allow(clippy::items_after_test_module)]

//! Picked-protein target-decoy inference.
//!
//! # What a protein group is here
//!
//! Two proteins of the same target/decoy class are **indistinguishable** when
//! the observed peptide evidence cannot tell them apart — that is, when the set
//! of identified peptides mapping to one is *exactly* the set mapping to the
//! other (Nesvizhskii & Aebersold 2005).  Only those proteins are collapsed
//! into a single reported group.  A target and a decoy with identical evidence
//! remain separate groups so that they can enter target/decoy competition.
//!
//! Grouping by connected components of the peptide-sharing graph is a different
//! and much coarser operation, and it was the previous behaviour of this module.
//! It merges proteins that the data plainly separates: if peptide `AB` maps to
//! both `A` and `B` while peptide `A1` maps only to `A`, then `A` is supported by
//! evidence `B` does not have, `{A}` and `{B}` are distinguishable, and reporting
//! them as one indistinguishable group misstates what was identified.  Sharing a
//! peptide is not the same as being indistinguishable, and a chain of shared
//! peptides can collapse an entire proteome into one "group".
//!
//! # Shared peptides
//!
//! A peptide that maps to several distinguishable groups counts as evidence for
//! every one of them.  No parsimony (razor-peptide assignment, subset removal,
//! or minimal-protein-set selection) is applied.  Subset proteins — where one
//! protein's peptide set is strictly contained in another's — are therefore
//! reported separately, because they *are* distinguishable.  This is a stated
//! limitation, not an inference: a parsimony step would change which proteins
//! are reported and is not implemented.
//!
//! # Error estimates
//!
//! Groups compete target-against-decoy by the picked-protein rule (Savitski et
//! al. 2015): each target group is paired with the decoy group carrying the
//! same decoy-stripped member names, and only the higher-scoring of the pair
//! enters the q-value scan.  Exact score ties between a target group and its
//! decoy counterpart are drawn with a fair coin keyed on the pairing key and the
//! run seed ([`crate::tiebreak`]), never by preferring the target.
//!
//! Picked-protein inference estimates a **cumulative** error rate over protein
//! groups.  It does not estimate a protein-level posterior error probability,
//! so [`ProtGroup::pep`] is `None` on this path and the output column reads
//! `NA`.  Reporting the best peptide's PEP there — the previous behaviour —
//! labelled a peptide-level quantity as a protein-level one.
//!
//! The separate `protein_bayes` module implements probabilistic noisy-OR
//! inference and does produce a protein-level posterior; it fills `pep` in.

use crate::stats;
use crate::tiebreak::Coin;
use std::collections::HashMap;

pub struct ProtGroup {
    pub proteins: Vec<String>,
    pub score: f64,
    pub qval: f64,
    /// Protein-level posterior error probability, when the selected inference
    /// method estimates one.
    ///
    /// `None` for picked-protein FDR, which produces a cumulative error-rate
    /// estimate and no posterior.  A missing value is reported as missing; it is
    /// never filled in from a peptide-level PEP.
    pub pep: Option<f64>,
    pub n_peptides: usize,
    pub is_decoy: bool,
    pub picked: bool, // won (or unpaired in) its target/decoy competition
}

const DECOY_PREFIXES: [&str; 4] = ["DECOY_", "REV_", "RANDOM_", "RANDOM-"];

pub(crate) fn is_decoy_protein(id: &str) -> bool {
    let u = id.to_ascii_uppercase();
    DECOY_PREFIXES.iter().any(|p| u.starts_with(p))
}

/// Strip a decoy prefix (case-insensitive) to recover the paired target name.
fn strip_decoy(id: &str) -> &str {
    let u = id.to_ascii_uppercase();
    for p in DECOY_PREFIXES {
        if u.starts_with(p) {
            return &id[p.len()..];
        }
    }
    id
}

/// Split the raw proteins field (tab- or space-separated protein ids).
#[allow(clippy::manual_pattern_char_comparison)]
pub fn split_proteins(s: &str) -> Vec<&str> {
    s.split(|c: char| c == '\t' || c == ' ' || c == ';')
        .filter(|p| !p.is_empty())
        .collect()
}

/// `entries`: one per peptide-level identification — (score, pep, raw_proteins_field).
///
/// `seed` is the run seed; it decides exact target/decoy picking ties and
/// nothing else.
pub fn infer(entries: &[(f64, f64, String)], seed: u64) -> Vec<ProtGroup> {
    #[cfg(feature = "profiling")]
    let _inference = crate::profile::Scope::with_elements(
        "protein_inference",
        "picked_protein_inference",
        entries.len(),
    );

    // Index protein ids and record, for each protein, the peptides that map to
    // it.  Entry index is the peptide identity: `entries` already holds one row
    // per distinct peptide form.
    let mut id_of: HashMap<&str, usize> = HashMap::new();
    let mut names: Vec<&str> = Vec::new();
    let mut evidence: Vec<Vec<u32>> = Vec::new();
    for (peptide, (_, _, raw)) in entries.iter().enumerate() {
        for protein in split_proteins(raw) {
            let index = *id_of.entry(protein).or_insert_with(|| {
                names.push(protein);
                evidence.push(Vec::new());
                names.len() - 1
            });
            evidence[index].push(peptide as u32);
        }
    }

    // Indistinguishability: proteins of the same target/decoy class with the
    // *same* observed peptide set.  Keeping the class in the key prevents a
    // target and a decoy from disappearing into a single mixed group before
    // target/decoy competition.
    for peptides in evidence.iter_mut() {
        peptides.sort_unstable();
        peptides.dedup();
    }
    let mut group_of: HashMap<(bool, Vec<u32>), usize> = HashMap::new();
    let mut group_members: Vec<Vec<usize>> = Vec::new();
    let mut group_evidence: Vec<(bool, Vec<u32>)> = Vec::new();
    for protein in 0..names.len() {
        let key = (
            is_decoy_protein(names[protein]),
            std::mem::take(&mut evidence[protein]),
        );
        match group_of.get(&key) {
            Some(&group) => group_members[group].push(protein),
            None => {
                group_of.insert(key.clone(), group_members.len());
                group_members.push(vec![protein]);
                group_evidence.push(key);
            }
        }
    }

    let mut out: Vec<ProtGroup> = group_members
        .iter()
        .zip(&group_evidence)
        .map(|(members, (is_decoy, peptides))| {
            let mut proteins: Vec<String> = members
                .iter()
                .map(|&index| names[index].to_string())
                .collect();
            proteins.sort();
            let score = peptides
                .iter()
                .map(|&peptide| entries[peptide as usize].0)
                .fold(f64::NEG_INFINITY, f64::max);
            ProtGroup {
                // Best-peptide score (Savitski picked FDR): the group's best
                // peptide discriminant, continuous and robust to group size.
                score,
                is_decoy: *is_decoy,
                proteins,
                qval: 1.0,
                // Picked-protein FDR estimates no protein-level posterior.
                pep: None,
                n_peptides: peptides.len(),
                picked: false,
            }
        })
        .collect();

    // Deterministic order before any tie-sensitive step, so nothing downstream
    // can depend on hash-map iteration.
    out.sort_by(|a, b| a.proteins.cmp(&b.proteins));
    picked_fdr(&mut out, seed);
    out.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.proteins.cmp(&b.proteins))
    });
    out
}

/// Classic (non-picked) protein FDR: naive target-decoy over *all* groups.
/// Reported only as a comparison point for the picked estimate; neither is a
/// validated protein-level error rate on real data.
pub fn classic_target_q01(groups: &[ProtGroup]) -> usize {
    let scores: Vec<f64> = groups.iter().map(|g| g.score).collect();
    let labels: Vec<i8> = groups
        .iter()
        .map(|g| if g.is_decoy { -1 } else { 1 })
        .collect();
    // Protein-level competition is between a target protein group and its own
    // reversed decoy group, one for one, independently of how many decoy
    // peptides the PSM-level search used.
    let q = stats::qvalues(&scores, &labels, stats::Tdc::reported(0.5));
    q.iter()
        .zip(labels.iter())
        .filter(|(qi, &l)| l > 0 && **qi < 0.01)
        .count()
}

/// Picked-protein FDR (Savitski et al. 2015): pair each target group with its
/// decoy counterpart (matched by decoy-stripped, canonicalized member names),
/// keep only the higher-scoring of the pair ("pick"), and compute q-values over
/// the picked entries.
fn picked_fdr(groups: &mut [ProtGroup], seed: u64) {
    #[cfg(feature = "profiling")]
    let _picked =
        crate::profile::Scope::with_elements("protein_inference", "picked_fdr", groups.len());
    // Pairing key = sorted set of decoy-stripped member names.  It carries no
    // label, which is what lets it seed a label-symmetric coin.
    let key_of = |g: &ProtGroup| -> String {
        let mut ks: Vec<&str> = g.proteins.iter().map(|p| strip_decoy(p)).collect();
        ks.sort_unstable();
        ks.dedup();
        ks.join("|")
    };

    // Bucket group indices by pairing key -> (best target idx, best decoy idx).
    // Within one slot a tie is broken by member names, which are unique per
    // group, so the choice is a function of content rather than of order.
    let mut buckets: HashMap<String, (Option<usize>, Option<usize>)> = HashMap::new();
    for gi in 0..groups.len() {
        let k = key_of(&groups[gi]);
        let e = buckets.entry(k).or_insert((None, None));
        let slot = if groups[gi].is_decoy {
            &mut e.1
        } else {
            &mut e.0
        };
        let replace = match *slot {
            None => true,
            Some(j) => {
                groups[gi].score > groups[j].score
                    || (groups[gi].score == groups[j].score
                        && groups[gi].proteins < groups[j].proteins)
            }
        };
        if replace {
            *slot = Some(gi);
        }
    }

    // One competition entry per bucket: the higher-scoring of target/decoy, with
    // an exact tie decided by a fair coin on the pairing key.
    let mut keys: Vec<&String> = buckets.keys().collect();
    keys.sort_unstable();
    let mut picks: Vec<(usize, f64, bool)> = Vec::with_capacity(buckets.len());
    for key in keys {
        let (t, d) = buckets[key];
        let pick = match (t, d) {
            (Some(ti), Some(di)) => {
                let target_wins = if groups[ti].score == groups[di].score {
                    Coin::new(seed).bytes(key.as_bytes()).heads()
                } else {
                    groups[ti].score > groups[di].score
                };
                if target_wins {
                    (ti, false)
                } else {
                    (di, true)
                }
            }
            (Some(ti), None) => (ti, false),
            (None, Some(di)) => (di, true),
            (None, None) => continue,
        };
        picks.push((pick.0, groups[pick.0].score, pick.1));
    }

    // q-values over the picked list (pi0 = 1)
    let scores: Vec<f64> = picks.iter().map(|p| p.1).collect();
    let labels: Vec<i8> = picks.iter().map(|p| if p.2 { -1 } else { 1 }).collect();
    let q = stats::qvalues(&scores, &labels, stats::Tdc::reported(0.5));
    for (pk, qi) in picks.iter().zip(q.into_iter()) {
        groups[pk.0].picked = true;
        groups[pk.0].qval = qi;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pin;
    use std::fs;
    use std::path::PathBuf;

    fn entry(score: f64, proteins: &str) -> (f64, f64, String) {
        (score, 1.0 / (1.0 + score.max(0.0)), proteins.to_string())
    }

    fn group_sets(groups: &[ProtGroup]) -> Vec<Vec<String>> {
        let mut sets: Vec<Vec<String>> = groups.iter().map(|g| g.proteins.clone()).collect();
        sets.sort();
        sets
    }

    #[derive(Debug, Eq, Ord, PartialEq, PartialOrd)]
    struct GroupSignature {
        proteins: Vec<String>,
        is_decoy: bool,
        score_bits: u64,
        qval_bits: u64,
        pep_bits: Option<u64>,
        n_peptides: usize,
        picked: bool,
    }

    fn inference_signature(groups: &[ProtGroup]) -> Vec<GroupSignature> {
        let mut signature: Vec<_> = groups
            .iter()
            .map(|group| GroupSignature {
                proteins: group.proteins.clone(),
                is_decoy: group.is_decoy,
                score_bits: group.score.to_bits(),
                qval_bits: group.qval.to_bits(),
                pep_bits: group.pep.map(f64::to_bits),
                n_peptides: group.n_peptides,
                picked: group.picked,
            })
            .collect();
        signature.sort();
        signature
    }

    fn deterministic_shuffle<T>(values: &mut [T], mut state: u64) {
        for index in (1..values.len()).rev() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            values.swap(index, (state % (index as u64 + 1)) as usize);
        }
    }

    /// Check a hand-derived graph under every scientifically irrelevant entry
    /// layout used by the protein audit.
    fn assert_graph(name: &str, entries: Vec<(f64, f64, String)>, expected: Vec<Vec<&str>>) {
        let mut expected: Vec<Vec<String>> = expected
            .into_iter()
            .map(|group| group.into_iter().map(str::to_string).collect())
            .collect();
        expected.sort();

        let mut reversed = entries.clone();
        reversed.reverse();
        let mut shuffled_17 = entries.clone();
        deterministic_shuffle(&mut shuffled_17, 17);
        let mut shuffled_91 = entries.clone();
        deterministic_shuffle(&mut shuffled_91, 91);
        let mut target_first = entries.clone();
        target_first.sort_by_key(|entry| {
            split_proteins(&entry.2)
                .iter()
                .all(|protein| is_decoy_protein(protein))
        });
        let mut decoy_first = target_first.clone();
        decoy_first.reverse();
        let reference = inference_signature(&infer(&entries, 31));

        for (arm, arranged) in [
            ("original", entries),
            ("reversed", reversed),
            ("shuffle-17", shuffled_17),
            ("shuffle-91", shuffled_91),
            ("target-first", target_first),
            ("decoy-first", decoy_first),
        ] {
            let groups = infer(&arranged, 31);
            assert_eq!(
                group_sets(&groups),
                expected,
                "{name}/{arm}: grouping does not match the hand-derived evidence graph"
            );
            assert_eq!(
                inference_signature(&groups),
                reference,
                "{name}/{arm}: score, q-value, PEP availability, peptide count, or pick changed"
            );
        }
    }

    /// Independent hand graphs defining the grouping relation. A group is an
    /// equivalence class of proteins with the same observed peptide set and the
    /// same target/decoy class. Sharing, subset relations, or connectedness are
    /// not equivalence.
    #[test]
    fn hand_derived_grouping_graphs_survive_all_entry_permutations() {
        assert_graph(
            "identical peptide evidence",
            vec![entry(10.0, "A B"), entry(9.0, "A B")],
            vec![vec!["A", "B"]],
        );
        assert_graph(
            "partially shared evidence",
            vec![entry(10.0, "A B"), entry(9.0, "B C")],
            vec![vec!["A"], vec!["B"], vec!["C"]],
        );
        assert_graph(
            "one unique peptide",
            vec![entry(10.0, "A B"), entry(9.0, "A")],
            vec![vec!["A"], vec!["B"]],
        );
        assert_graph(
            "multiple unique peptides",
            vec![
                entry(10.0, "A B"),
                entry(9.0, "A"),
                entry(8.0, "A"),
                entry(7.0, "B"),
                entry(6.0, "B"),
            ],
            vec![vec!["A"], vec!["B"]],
        );
        assert_graph(
            "strict subset",
            vec![entry(10.0, "SUB SUPER"), entry(9.0, "SUPER")],
            vec![vec!["SUB"], vec!["SUPER"]],
        );
        assert_graph(
            "disjoint proteins",
            vec![entry(10.0, "A"), entry(9.0, "B")],
            vec![vec!["A"], vec!["B"]],
        );
        assert_graph(
            "target-target indistinguishability",
            vec![entry(10.0, "TARGET_A TARGET_B")],
            vec![vec!["TARGET_A", "TARGET_B"]],
        );
        assert_graph(
            "decoy-decoy indistinguishability",
            vec![entry(10.0, "DECOY_A DECOY_B")],
            vec![vec!["DECOY_A", "DECOY_B"]],
        );
        assert_graph(
            "target-decoy evidence",
            vec![entry(10.0, "MIXED DECOY_MIXED")],
            vec![vec!["DECOY_MIXED"], vec!["MIXED"]],
        );
        assert_graph(
            "exact target-decoy score tie",
            vec![entry(5.0, "PAIR"), entry(5.0, "DECOY_PAIR")],
            vec![vec!["DECOY_PAIR"], vec!["PAIR"]],
        );
        assert_graph(
            "near target-decoy score tie",
            vec![
                entry(5.0, "NEAR"),
                entry(f64::from_bits(5.0f64.to_bits() - 1), "DECOY_NEAR"),
            ],
            vec![vec!["DECOY_NEAR"], vec!["NEAR"]],
        );
    }

    /// A target and decoy with identical peptide evidence are competitors, not
    /// indistinguishable members of one group. Collapsing them destroys the
    /// target/decoy axis before picked competition can operate.
    #[test]
    fn target_and_decoy_proteins_never_collapse_into_one_group() {
        let groups = infer(&[entry(11.0, "MIXED DECOY_MIXED")], 31);
        assert_eq!(
            group_sets(&groups),
            vec![vec!["DECOY_MIXED".to_string()], vec!["MIXED".to_string()]]
        );
        assert_eq!(groups.iter().filter(|group| group.picked).count(), 1);
        assert!(groups.iter().all(|group| {
            group
                .proteins
                .iter()
                .all(|protein| is_decoy_protein(protein) == group.is_decoy)
        }));
        assert_eq!(
            inference_signature(&infer(&[entry(11.0, "DECOY_MIXED MIXED")], 31)),
            inference_signature(&groups),
            "protein order within a peptide mapping changed the competition"
        );
    }

    /// Proteins whose observed peptide sets are identical cannot be told apart
    /// by the data, so they are one group.
    #[test]
    fn indistinguishable_proteins_are_one_group() {
        let groups = infer(
            &[entry(10.0, "A\tB"), entry(8.0, "A\tB"), entry(5.0, "C")],
            1,
        );
        assert_eq!(
            group_sets(&groups),
            vec![
                vec!["A".to_string(), "B".to_string()],
                vec!["C".to_string()]
            ]
        );
    }

    /// **Frozen failure case (independent audit M2).**
    ///
    /// `A` and `B` share peptide `AB`, but `A` also has peptide `A1` that `B`
    /// does not.  Their evidence differs, so they are distinguishable and must
    /// not be reported as one indistinguishable group.  Connected-component
    /// grouping merges them.
    #[test]
    fn a_unique_peptide_keeps_two_proteins_apart() {
        let groups = infer(&[entry(10.0, "A\tB"), entry(9.0, "A")], 1);
        assert_eq!(
            group_sets(&groups),
            vec![vec!["A".to_string()], vec!["B".to_string()]],
            "proteins with different peptide evidence were merged"
        );
        let a = groups.iter().find(|g| g.proteins == ["A"]).unwrap();
        let b = groups.iter().find(|g| g.proteins == ["B"]).unwrap();
        // The shared peptide is evidence for both; the unique one only for A.
        assert_eq!(a.n_peptides, 2);
        assert_eq!(b.n_peptides, 1);
    }

    /// A subset protein is distinguishable from its superset: the superset has
    /// evidence the subset lacks.
    #[test]
    fn subset_proteins_are_not_absorbed() {
        let groups = infer(
            &[
                entry(10.0, "SUB\tSUPER"),
                entry(9.0, "SUB\tSUPER"),
                entry(8.0, "SUPER"),
            ],
            1,
        );
        assert_eq!(
            group_sets(&groups),
            vec![vec!["SUB".to_string()], vec!["SUPER".to_string()]]
        );
        assert_eq!(
            groups
                .iter()
                .find(|g| g.proteins == ["SUPER"])
                .unwrap()
                .n_peptides,
            3
        );
        assert_eq!(
            groups
                .iter()
                .find(|g| g.proteins == ["SUB"])
                .unwrap()
                .n_peptides,
            2
        );
    }

    /// A chain of shared peptides must not collapse a whole component: only the
    /// pairs with identical evidence group.
    #[test]
    fn a_sharing_chain_does_not_collapse_into_one_group() {
        // A-B share p1, B-C share p2, C-D share p3.  Nothing is indistinguishable.
        let groups = infer(
            &[entry(10.0, "A\tB"), entry(9.0, "B\tC"), entry(8.0, "C\tD")],
            1,
        );
        assert_eq!(
            group_sets(&groups),
            vec![
                vec!["A".to_string()],
                vec!["B".to_string()],
                vec!["C".to_string()],
                vec!["D".to_string()]
            ]
        );
    }

    /// Several independent groups coexist, and a shared peptide is evidence for
    /// each group it maps to.
    #[test]
    fn shared_peptides_support_every_group_they_map_to() {
        let groups = infer(
            &[
                entry(10.0, "A1\tA2"), // A1, A2 indistinguishable so far
                entry(9.0, "A1\tA2"),
                entry(8.0, "A1\tA2\tB"), // shared across the A group and B
                entry(7.0, "B"),
            ],
            1,
        );
        assert_eq!(
            group_sets(&groups),
            vec![
                vec!["A1".to_string(), "A2".to_string()],
                vec!["B".to_string()]
            ]
        );
        let a = groups.iter().find(|g| g.proteins.len() == 2).unwrap();
        let b = groups.iter().find(|g| g.proteins == ["B"]).unwrap();
        assert_eq!(a.n_peptides, 3);
        assert_eq!(b.n_peptides, 2);
    }

    /// Grouping is a function of the evidence, not of the order the peptides
    /// arrive in.
    #[test]
    fn grouping_is_invariant_to_entry_order() {
        let base = vec![
            entry(10.0, "A\tB"),
            entry(9.0, "A"),
            entry(8.0, "C\tD"),
            entry(7.0, "C\tD"),
            entry(6.0, "E"),
        ];
        let reference = group_sets(&infer(&base, 1));
        let mut permuted = base.clone();
        permuted.reverse();
        assert_eq!(group_sets(&infer(&permuted, 1)), reference);
        permuted.rotate_left(2);
        assert_eq!(group_sets(&infer(&permuted, 1)), reference);
    }

    #[test]
    fn decoy_proteins_flagged() {
        let groups = infer(
            &[entry(9.0, "sp|P1|REAL"), entry(7.0, "DECOY_sp|P1|REAL")],
            1,
        );
        assert!(groups.iter().any(|g| g.is_decoy));
        assert!(groups.iter().any(|g| !g.is_decoy));
    }

    #[test]
    fn picked_keeps_the_higher_of_a_target_decoy_pair() {
        let groups = infer(
            &[entry(9.0, "sp|P1|REAL"), entry(7.0, "DECOY_sp|P1|REAL")],
            1,
        );
        assert!(groups.iter().find(|g| !g.is_decoy).unwrap().picked);
        assert!(!groups.iter().find(|g| g.is_decoy).unwrap().picked);

        let other = infer(&[entry(3.0, "sp|P2|X"), entry(8.0, "DECOY_sp|P2|X")], 1);
        assert!(other.iter().find(|g| g.is_decoy).unwrap().picked);
        assert!(!other.iter().find(|g| !g.is_decoy).unwrap().picked);
    }

    fn tied_protein_pairs(pairs: usize) -> Vec<(f64, f64, String)> {
        let mut entries = Vec::new();
        for index in 0..pairs {
            entries.push(entry(50.0, &format!("P{index:04}")));
            entries.push(entry(50.0, &format!("DECOY_P{index:04}")));
        }
        entries
    }

    /// **Frozen failure case (independent audit C2).**
    ///
    /// Exact target/decoy protein-score ties must not all go to the target.
    #[test]
    fn exact_picked_ties_do_not_all_go_to_the_target() {
        let groups = infer(&tied_protein_pairs(400), 1);
        let picked_targets = groups.iter().filter(|g| g.picked && !g.is_decoy).count();
        assert_eq!(groups.iter().filter(|g| g.picked).count(), 400);
        // A fair coin over 400 pairs has SD 10; 5 SD is 50.
        assert!(
            (picked_targets as i64 - 200).abs() < 50,
            "{picked_targets} of 400 exactly tied protein pairs were won by the target"
        );
    }

    /// The same tie fixture with the entries permuted must produce the same
    /// picks: input order carries no scientific information.
    #[test]
    fn picked_ties_are_invariant_to_entry_order() {
        let base = tied_protein_pairs(200);
        let reference: Vec<(Vec<String>, bool)> = {
            let mut rows: Vec<(Vec<String>, bool)> = infer(&base, 1)
                .into_iter()
                .filter(|g| g.picked)
                .map(|g| (g.proteins, g.is_decoy))
                .collect();
            rows.sort();
            rows
        };
        let mut permuted = base.clone();
        permuted.reverse();
        let mut reversed: Vec<(Vec<String>, bool)> = infer(&permuted, 1)
            .into_iter()
            .filter(|g| g.picked)
            .map(|g| (g.proteins, g.is_decoy))
            .collect();
        reversed.sort();
        assert_eq!(
            reversed, reference,
            "reversing the entries changed the picks"
        );

        // Interleave decoys first.
        let mut swapped = base.clone();
        for index in (0..swapped.len()).step_by(2) {
            swapped.swap(index, index + 1);
        }
        let mut swapped_picks: Vec<(Vec<String>, bool)> = infer(&swapped, 1)
            .into_iter()
            .filter(|g| g.picked)
            .map(|g| (g.proteins, g.is_decoy))
            .collect();
        swapped_picks.sort();
        assert_eq!(
            swapped_picks, reference,
            "putting the decoy row first changed the picks"
        );
    }

    #[test]
    fn a_different_seed_reflips_protein_ties() {
        let base = tied_protein_pairs(200);
        let picks = |seed: u64| -> Vec<bool> {
            let mut rows: Vec<(Vec<String>, bool)> = infer(&base, seed)
                .into_iter()
                .filter(|g| g.picked)
                .map(|g| (g.proteins.clone(), g.is_decoy))
                .collect();
            rows.sort();
            rows.into_iter().map(|(_, decoy)| decoy).collect()
        };
        assert_ne!(picks(1), picks(2));
    }

    /// A strict winner is never displaced by the coin.
    #[test]
    fn a_strictly_better_group_always_wins_its_pair() {
        for seed in 1..=16u64 {
            let groups = infer(&[entry(9.0, "P"), entry(8.99, "DECOY_P")], seed);
            assert!(groups.iter().find(|g| !g.is_decoy).unwrap().picked);
            let groups = infer(&[entry(8.99, "P"), entry(9.0, "DECOY_P")], seed);
            assert!(groups.iter().find(|g| g.is_decoy).unwrap().picked);
        }
    }

    /// **Frozen failure case (independent audit M3).**
    ///
    /// Picked-protein inference estimates no protein-level posterior, so it must
    /// not emit one — least of all the best peptide's PEP under a protein-level
    /// column name.
    #[test]
    fn picked_inference_reports_no_protein_posterior() {
        let entries = vec![(12.0, 0.123456, "sp|P9|ONLY".to_string())];
        let groups = infer(&entries, 1);
        assert_eq!(groups.len(), 1);
        assert_eq!(
            groups[0].pep, None,
            "picked inference emitted a protein PEP it does not estimate"
        );
    }

    /// Picking removes the loser of every paired competition, so the picked list
    /// holds one entry per pairing key.  This is a structural property of the
    /// method, not a claim that it finds more or better proteins.
    #[test]
    fn picking_leaves_one_entry_per_pairing_key() {
        let entries = vec![
            entry(20.0, "A"),
            entry(5.0, "DECOY_A"),
            entry(18.0, "B"),
            entry(16.0, "C"),
            entry(4.0, "DECOY_X"),
            entry(3.5, "DECOY_Y"),
        ];
        let groups = infer(&entries, 1);
        assert_eq!(groups.iter().filter(|g| g.picked).count(), 5);
        assert!(
            !groups
                .iter()
                .find(|g| g.proteins == ["DECOY_A"])
                .unwrap()
                .picked
        );
    }

    #[test]
    fn synthetic_pin_fixture_groups_and_picks_deterministically() {
        let path = write_synthetic_pin_fixture();
        let ds = pin::parse(path.to_str().unwrap()).expect("synthetic PIN should parse");
        fs::remove_file(&path).ok();

        assert_eq!(ds.n_feat, 1, "fixture should expose one score feature");
        assert_eq!(ds.n_psm, 607, "fixture row count drifted");

        let entries: Vec<(f64, f64, String)> = (0..ds.n_psm)
            .map(|i| {
                let score = ds.row(i)[0];
                (score, 1.0 / (1.0 + score.max(0.0)), ds.proteins[i].clone())
            })
            .collect();
        let groups = infer(&entries, 1);

        // SHARED_A and SHARED_B carry exactly the same three peptides, so they
        // are genuinely indistinguishable and stay one group.
        let shared = groups
            .iter()
            .find(|g| g.proteins.iter().any(|p| p == "SHARED_A"))
            .expect("shared target group should exist");
        assert_eq!(shared.proteins, vec!["SHARED_A", "SHARED_B"]);
        assert_eq!(shared.n_peptides, 3);
        assert!(shared.picked);

        // Every group carries the counts the graph implies, and none reports a
        // protein posterior.
        assert!(groups.iter().all(|g| g.pep.is_none()));
        let picked = groups
            .iter()
            .filter(|g| !g.is_decoy && g.picked && g.qval < 0.01)
            .count();
        let classic = classic_target_q01(&groups);
        // Recorded, not required to move in either direction.
        assert_eq!(classic, 0, "classic q<0.01 count drifted");
        assert_eq!(picked, 121, "picked q<0.01 count drifted");
    }

    fn write_synthetic_pin_fixture() -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "percolator-rs-protein-fixture-{}.pin",
            std::process::id()
        ));
        let mut pin = String::from("SpecId\tLabel\tScanNr\tscore\tPeptide\tProteins\n");
        let mut scan = 1usize;

        for i in 0..120 {
            let target = format!("T{:03}", i);
            let target_score = 1000.0 - i as f64;
            append_group(
                &mut pin,
                &mut scan,
                1,
                target_score,
                &target,
                &[target.as_str()],
                3,
            );

            let decoy = format!("DECOY_{target}");
            let decoy_score = if i < 40 {
                920.5 - i as f64
            } else {
                100.0 - (i - 40) as f64
            };
            append_group(
                &mut pin,
                &mut scan,
                -1,
                decoy_score,
                &target,
                &[decoy.as_str()],
                2,
            );
        }

        append_group(
            &mut pin,
            &mut scan,
            1,
            930.25,
            "SHARED",
            &["SHARED_A", "SHARED_B"],
            3,
        );
        append_group(
            &mut pin,
            &mut scan,
            -1,
            910.25,
            "SHARED",
            &["DECOY_SHARED_A", "DECOY_SHARED_B"],
            2,
        );

        append_group(
            &mut pin,
            &mut scan,
            -1,
            40.0,
            "LONE_DEC",
            &["DECOY_LONE_A"],
            1,
        );
        append_group(
            &mut pin,
            &mut scan,
            -1,
            39.0,
            "LONE_DEC",
            &["DECOY_LONE_B"],
            1,
        );

        fs::write(&path, pin).expect("synthetic PIN fixture should write");
        path
    }

    fn append_group(
        out: &mut String,
        scan: &mut usize,
        label: i8,
        best_score: f64,
        tag: &str,
        proteins: &[&str],
        n_peptides: usize,
    ) {
        for pep_idx in 0..n_peptides {
            let peptide = format!("K.{}_{:03}.R", tag, pep_idx);
            let spec_id = format!("{tag}_{label}_{pep_idx}");
            let score = best_score - pep_idx as f64 * 0.01;
            out.push_str(&spec_id);
            out.push('\t');
            out.push_str(if label > 0 { "1" } else { "-1" });
            out.push('\t');
            out.push_str(&scan.to_string());
            out.push('\t');
            out.push_str(&format!("{score:.2}"));
            out.push('\t');
            out.push_str(&peptide);
            for protein in proteins {
                out.push('\t');
                out.push_str(protein);
            }
            out.push('\n');
            *scan += 1;
        }
    }
}
