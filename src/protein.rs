#![allow(clippy::items_after_test_module)]

//! Picked-protein target-decoy inference.
//!
//! We build the bipartite peptide<->protein graph, collapse proteins that share peptides
//! into indistinguishable groups (union-find = connected components), score each group by
//! its best peptide, and compute protein-level q-values / PEPs by target-decoy competition.
//! Decoy proteins are identified by the `DECOY_` prefix.
//!
//! The separate `protein_bayes` module implements probabilistic noisy-OR inference.

use crate::stats;
use std::collections::HashMap;

pub struct ProtGroup {
    pub proteins: Vec<String>,
    pub score: f64,
    pub qval: f64,
    pub pep: f64,
    pub n_peptides: usize,
    pub is_decoy: bool,
    pub picked: bool, // won (or unpaired in) its target/decoy competition
}

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        // path compression
        let mut c = x;
        while self.parent[c] != r {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pin;
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn groups_proteins_sharing_a_peptide() {
        // peptide1 -> {A, B} (shared) ; peptide2 -> {B} ; peptide3 -> {C}
        let entries = vec![
            (10.0, 0.001, "A\tB".to_string()),
            (8.0, 0.01, "B".to_string()),
            (5.0, 0.2, "C".to_string()),
        ];
        let groups = infer(&entries);
        // A and B collapse into one group; C separate -> 2 groups
        assert_eq!(groups.len(), 2);
        let ab = groups
            .iter()
            .find(|g| g.proteins.contains(&"A".to_string()))
            .unwrap();
        assert!(ab.proteins.contains(&"B".to_string()));
    }

    #[test]
    fn decoy_proteins_flagged() {
        let entries = vec![
            (9.0, 0.001, "sp|P1|REAL".to_string()),
            (7.0, 0.01, "DECOY_sp|P1|REAL".to_string()),
        ];
        let groups = infer(&entries);
        assert!(groups.iter().any(|g| g.is_decoy));
        assert!(groups.iter().any(|g| !g.is_decoy));
    }

    #[test]
    fn picked_keeps_higher_of_target_decoy_pair() {
        // same protein: target scores 9, its decoy scores 7 -> target is picked, decoy dropped
        let entries = vec![
            (9.0, 0.001, "sp|P1|REAL".to_string()),
            (7.0, 0.01, "DECOY_sp|P1|REAL".to_string()),
        ];
        let groups = infer(&entries);
        let t = groups.iter().find(|g| !g.is_decoy).unwrap();
        let d = groups.iter().find(|g| g.is_decoy).unwrap();
        assert!(t.picked, "higher-scoring target should be picked");
        assert!(
            !d.picked,
            "lower-scoring decoy counterpart should be dropped"
        );
        // and if the decoy outscores the target, the decoy is picked instead
        let entries2 = vec![
            (3.0, 0.5, "sp|P2|X".to_string()),
            (8.0, 0.01, "DECOY_sp|P2|X".to_string()),
        ];
        let g2 = infer(&entries2);
        assert!(g2.iter().find(|g| g.is_decoy).unwrap().picked);
        assert!(!g2.iter().find(|g| !g.is_decoy).unwrap().picked);
    }

    #[test]
    fn picked_never_less_sensitive_than_classic() {
        // strong targets, a weaker decoy counterpart (removed by picking), and lone decoys
        let entries = vec![
            (20.0, 1e-8, "A".to_string()),
            (5.0, 0.30, "DECOY_A".to_string()), // paired, target wins -> decoy dropped
            (18.0, 1e-6, "B".to_string()),
            (16.0, 1e-5, "C".to_string()),
            (4.0, 0.40, "DECOY_X".to_string()), // lone decoy
            (3.5, 0.50, "DECOY_Y".to_string()), // lone decoy
        ];
        let g = infer(&entries);
        let picked = g
            .iter()
            .filter(|x| !x.is_decoy && x.picked && x.qval < 0.01)
            .count();
        let classic = classic_target_q01(&g);
        assert!(
            picked >= classic,
            "picked FDR must be >= classic (got {picked} vs {classic})"
        );
    }

    #[test]
    fn synthetic_pin_fixture_shows_picked_gain_on_realistic_groups() {
        let path = write_synthetic_pin_fixture();
        let ds = pin::parse(path.to_str().unwrap()).expect("synthetic PIN should parse");
        fs::remove_file(&path).ok();

        assert_eq!(ds.n_feat, 1, "fixture should expose one score feature");
        assert_eq!(ds.n_psm, 607, "fixture row count drifted");

        let entries: Vec<(f64, f64, String)> = (0..ds.n_psm)
            .map(|i| {
                let score = ds.row(i)[0];
                let pep = 1.0 / (1.0 + score.max(0.0));
                (score, pep, ds.proteins[i].clone())
            })
            .collect();
        let groups = infer(&entries);

        let picked = groups
            .iter()
            .filter(|g| !g.is_decoy && g.picked && g.qval < 0.01)
            .count();
        let classic = classic_target_q01(&groups);
        // Counts under the corrected target-decoy estimator.  Classic FDR used
        // to report 81 here only because the pre-repair estimator let a leading
        // decoy-free run reach q = 0; with the finite-sample safeguard the best
        // achievable estimate on this fixture is 1/81 = 0.0123, above the
        // threshold.  Picking removes the paired decoys, so the picked list still
        // clears 1% at 1/121 = 0.0083.
        assert_eq!(classic, 0, "classic q<0.01 count drifted");
        assert_eq!(picked, 121, "picked q<0.01 count drifted");
        assert!(
            picked > classic,
            "synthetic grouped fixture should favor picked FDR ({picked} vs {classic})"
        );

        let shared = groups
            .iter()
            .find(|g| {
                g.proteins.iter().any(|p| p == "SHARED_A")
                    && g.proteins.iter().any(|p| p == "SHARED_B")
            })
            .expect("shared target group should exist");
        assert_eq!(
            shared.n_peptides, 3,
            "shared target group should aggregate its peptides"
        );
        assert!(
            shared.picked,
            "shared target group should win its picked competition"
        );
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
pub(crate) fn split_proteins(s: &str) -> Vec<&str> {
    s.split(|c: char| c == '\t' || c == ' ' || c == ';')
        .filter(|p| !p.is_empty())
        .collect()
}

/// `entries`: one per peptide-level identification — (score, pep, raw_proteins_field).
pub fn infer(entries: &[(f64, f64, String)]) -> Vec<ProtGroup> {
    #[cfg(feature = "profiling")]
    let _inference = crate::profile::Scope::with_elements(
        "protein_inference",
        "picked_protein_inference",
        entries.len(),
    );
    #[cfg(feature = "profiling")]
    let mut split_vector_calls = 0u64;
    #[cfg(feature = "profiling")]
    let mut split_vector_bytes = 0u64;
    // index protein ids
    let mut id_of: HashMap<&str, usize> = HashMap::new();
    let mut names: Vec<&str> = Vec::new();
    for (_, _, raw) in entries {
        let proteins = split_proteins(raw);
        #[cfg(feature = "profiling")]
        {
            split_vector_calls += 1;
            split_vector_bytes += (proteins.capacity() * std::mem::size_of::<&str>()) as u64;
        }
        for p in proteins {
            if !id_of.contains_key(p) {
                id_of.insert(p, names.len());
                names.push(p);
            }
        }
    }
    let n_prot = names.len();
    let mut uf = UnionFind::new(n_prot);

    // pass 1: union all proteins that co-occur in a peptide (shared-peptide grouping)
    for (_, _, raw) in entries {
        let prots = split_proteins(raw);
        #[cfg(feature = "profiling")]
        {
            split_vector_calls += 1;
            split_vector_bytes += (prots.capacity() * std::mem::size_of::<&str>()) as u64;
        }
        if prots.len() > 1 {
            let first = id_of[prots[0]];
            for p in &prots[1..] {
                uf.union(first, id_of[*p]);
            }
        }
    }

    // pass 2: accumulate per-group evidence. Protein score = Σ −ln(PEP) over its distinct
    // peptides (Fisher combination): real proteins with several good peptides score high;
    // decoys with sporadic single hits stay near zero. Each peptide belongs to exactly one
    // group (all its proteins are unioned), so it contributes once.
    struct G {
        proteins: std::collections::HashSet<usize>,
        best_score: f64, // best peptide SVM score (continuous — drives the protein score)
        best_pep: f64,   // PEP of that best peptide (reported)
        n_peptides: usize,
    }
    let mut groups: HashMap<usize, G> = HashMap::new();
    for (score, pep, raw) in entries {
        let prots = split_proteins(raw);
        #[cfg(feature = "profiling")]
        {
            split_vector_calls += 1;
            split_vector_bytes += (prots.capacity() * std::mem::size_of::<&str>()) as u64;
        }
        if prots.is_empty() {
            continue;
        }
        let root = uf.find(id_of[prots[0]]);
        let g = groups.entry(root).or_insert_with(|| G {
            proteins: std::collections::HashSet::new(),
            best_score: f64::NEG_INFINITY,
            best_pep: 1.0,
            n_peptides: 0,
        });
        for p in &prots {
            g.proteins.insert(id_of[*p]);
        }
        // track the most-identifying peptide (lowest PEP), remembering its SVM score
        let clamped = pep.clamp(1e-12, 1.0);
        if clamped < g.best_pep || (clamped == g.best_pep && *score > g.best_score) {
            g.best_pep = clamped;
            g.best_score = *score;
        }
        g.n_peptides += 1;
    }
    #[cfg(feature = "profiling")]
    crate::profile::allocation_site(
        "protein::infer split protein vectors",
        split_vector_calls,
        split_vector_bytes,
    );

    // materialize
    #[cfg(feature = "profiling")]
    let mut member_sort_duration = std::time::Duration::ZERO;
    #[cfg(feature = "profiling")]
    let mut member_sort_elements = 0u64;
    let mut out: Vec<ProtGroup> = groups
        .into_values()
        .map(|g| {
            let is_decoy = g.proteins.iter().any(|&pi| is_decoy_protein(names[pi]));
            let mut prot_names: Vec<String> =
                g.proteins.iter().map(|&pi| names[pi].to_string()).collect();
            #[cfg(feature = "profiling")]
            let member_sort_start = std::time::Instant::now();
            prot_names.sort();
            #[cfg(feature = "profiling")]
            {
                member_sort_duration += member_sort_start.elapsed();
                member_sort_elements += prot_names.len() as u64;
            }
            ProtGroup {
                // best-peptide score (Savitski picked FDR): the group's best peptide SVM
                // discriminant — continuous (no −ln(PEP=0) saturation ties that scramble the
                // FDR ranking) and robust to group size (unlike a Σ−ln(PEP) combination).
                score: g.best_score,
                proteins: prot_names,
                qval: 1.0,
                pep: g.best_pep,
                n_peptides: g.n_peptides,
                is_decoy,
                picked: false,
            }
        })
        .collect();
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "sort",
        "protein_group_member_order",
        member_sort_duration,
        Some(member_sort_elements),
        None,
    );

    picked_fdr(&mut out);
    #[cfg(feature = "profiling")]
    let output_sort_start = std::time::Instant::now();
    out.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "sort",
        "protein_group_score_order",
        output_sort_start.elapsed(),
        Some(out.len() as u64),
        None,
    );
    out
}

/// Classic (non-picked) protein FDR: naive target-decoy over *all* groups. Returned only
/// to quantify the sensitivity gain from picked FDR (picked >= classic at the same q).
pub fn classic_target_q01(groups: &[ProtGroup]) -> usize {
    let scores: Vec<f64> = groups.iter().map(|g| g.score).collect();
    let labels: Vec<i8> = groups
        .iter()
        .map(|g| if g.is_decoy { -1 } else { 1 })
        .collect();
    let q = stats::qvalues(&scores, &labels, stats::Tdc::reported(0.5));
    q.iter()
        .zip(labels.iter())
        .filter(|(qi, &l)| l > 0 && **qi < 0.01)
        .count()
}

/// Picked-protein FDR (Savitski et al. 2015): pair each target group with its decoy
/// counterpart (matched by decoy-stripped, canonicalized member names), keep only the
/// higher-scoring of the pair ("pick"), and compute q-values over the picked entries.
/// This halves double-counting and is more sensitive & better-calibrated than naive TDA.
#[allow(clippy::useless_conversion)]
fn picked_fdr(groups: &mut [ProtGroup]) {
    #[cfg(feature = "profiling")]
    let _picked =
        crate::profile::Scope::with_elements("protein_inference", "picked_fdr", groups.len());
    #[cfg(feature = "profiling")]
    let mut key_sort_duration = std::time::Duration::ZERO;
    #[cfg(feature = "profiling")]
    let mut key_sort_elements = 0u64;
    #[cfg(feature = "profiling")]
    let mut key_vector_bytes = 0u64;
    // pairing key = sorted set of decoy-stripped member names
    #[allow(unused_mut)]
    let mut key_of = |g: &ProtGroup| -> String {
        let mut ks: Vec<&str> = g.proteins.iter().map(|p| strip_decoy(p)).collect();
        #[cfg(feature = "profiling")]
        {
            key_vector_bytes += (ks.capacity() * std::mem::size_of::<&str>()) as u64;
            key_sort_elements += ks.len() as u64;
        }
        #[cfg(feature = "profiling")]
        let key_sort_start = std::time::Instant::now();
        ks.sort_unstable();
        #[cfg(feature = "profiling")]
        {
            key_sort_duration += key_sort_start.elapsed();
        }
        ks.dedup();
        ks.join("|")
    };

    // bucket group indices by pairing key -> (best target idx, best decoy idx)
    let mut buckets: HashMap<String, (Option<usize>, Option<usize>)> = HashMap::new();
    for gi in 0..groups.len() {
        let k = key_of(&groups[gi]);
        let e = buckets.entry(k).or_insert((None, None));
        let slot = if groups[gi].is_decoy {
            &mut e.1
        } else {
            &mut e.0
        };
        match *slot {
            Some(j) if groups[j].score >= groups[gi].score => {}
            _ => *slot = Some(gi),
        }
    }
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "sort",
            "protein_picked_key_member_order",
            key_sort_duration,
            Some(key_sort_elements),
            None,
        );
        crate::profile::allocation_site(
            "protein::picked_fdr key vectors",
            groups.len() as u64,
            key_vector_bytes,
        );
    }

    // one competition entry per bucket: the higher-scoring of target/decoy
    let mut picks: Vec<(usize, f64, bool)> = Vec::with_capacity(buckets.len()); // (group idx, score, is_decoy)
    for (t, d) in buckets.values() {
        let pick = match (t, d) {
            (Some(ti), Some(di)) => {
                if groups[*ti].score >= groups[*di].score {
                    (*ti, false)
                } else {
                    (*di, true)
                }
            }
            (Some(ti), None) => (*ti, false),
            (None, Some(di)) => (*di, true),
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
