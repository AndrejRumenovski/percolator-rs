//! Fido-style Bayesian protein inference.
//!
//! Proteins have independent presence prior `gamma`. A present protein emits each
//! adjacent peptide with probability `alpha`, while `beta` models an observation
//! arising from noise. Peptide PEPs are converted to likelihood ratios against a
//! configurable peptide prior. Sum-product belief propagation marginalizes protein
//! presence. It is exact for tree-structured components and a deterministic, damped
//! loopy-BP approximation for cyclic components.

use crate::protein::{is_decoy_protein, split_proteins, ProtGroup};
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Params {
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    pub peptide_prior: f64,
    pub max_iter: usize,
    pub tolerance: f64,
    pub damping: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            alpha: 0.1,
            beta: 0.01,
            gamma: 0.5,
            peptide_prior: 0.1,
            max_iter: 200,
            tolerance: 1e-10,
            damping: 0.5,
        }
    }
}

impl Params {
    pub fn validate(&self) -> Result<(), &'static str> {
        if !self.alpha.is_finite() || !(0.0..=1.0).contains(&self.alpha) || self.alpha == 0.0 {
            return Err("alpha must be finite and in (0, 1]");
        }
        if !self.beta.is_finite() || !(0.0..1.0).contains(&self.beta) {
            return Err("beta must be finite and in [0, 1)");
        }
        if !self.gamma.is_finite() || !(0.0..1.0).contains(&self.gamma) {
            return Err("gamma must be finite and in (0, 1)");
        }
        if !self.peptide_prior.is_finite() || !(0.0..1.0).contains(&self.peptide_prior) {
            return Err("peptide-prior must be finite and in (0, 1)");
        }
        if self.max_iter == 0 {
            return Err("max-iter must be greater than zero");
        }
        if !self.tolerance.is_finite() || self.tolerance <= 0.0 {
            return Err("tolerance must be finite and positive");
        }
        if !self.damping.is_finite() || !(0.0..=1.0).contains(&self.damping) || self.damping == 0.0
        {
            return Err("damping must be finite and in (0, 1]");
        }
        Ok(())
    }
}

#[derive(Default, Debug)]
pub struct Diagnostics {
    pub components: usize,
    pub tree_components: usize,
    pub loopy_components: usize,
    pub iterations: usize,
    pub converged: bool,
}

pub struct InferenceResult {
    pub groups: Vec<ProtGroup>,
    pub diagnostics: Diagnostics,
}

struct ProteinCluster {
    proteins: Vec<String>,
    peptides: Vec<usize>,
    edges: Vec<usize>,
    prior: Vec<f64>,
}

struct PeptideFactor {
    edges: Vec<usize>,
    likelihood_ratio: f64,
}

struct Edge {
    cluster: usize,
    factor: usize,
    variable_to_factor: Vec<f64>,
    factor_to_variable: Vec<f64>,
}

struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        while self.parent[x] != root {
            let next = self.parent[x];
            self.parent[x] = root;
            x = next;
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let a = self.find(a);
        let b = self.find(b);
        if a != b {
            self.parent[b] = a;
        }
    }
}

fn normalize(values: &mut [f64]) {
    let sum: f64 = values.iter().sum();
    if sum.is_finite() && sum > 0.0 {
        for value in values {
            *value /= sum;
        }
    } else {
        let uniform = 1.0 / values.len() as f64;
        values.fill(uniform);
    }
}

fn softmax(log_values: &[f64]) -> Vec<f64> {
    let max = log_values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut values: Vec<f64> = log_values
        .iter()
        .map(|value| (*value - max).exp())
        .collect();
    normalize(&mut values);
    values
}

fn binomial_prior(size: usize, gamma: f64) -> Vec<f64> {
    let log_odds = gamma.ln() - (1.0 - gamma).ln();
    let mut logs = Vec::with_capacity(size + 1);
    logs.push(size as f64 * (1.0 - gamma).ln());
    for k in 0..size {
        let next = logs[k] + ((size - k) as f64).ln() - ((k + 1) as f64).ln() + log_odds;
        logs.push(next);
    }
    softmax(&logs)
}

fn peptide_likelihood_ratio(pep: f64, prior_true: f64) -> f64 {
    let posterior_true = (1.0 - pep).clamp(1e-12, 1.0 - 1e-12);
    let prior_true = prior_true.clamp(1e-12, 1.0 - 1e-12);
    let posterior_odds = posterior_true / (1.0 - posterior_true);
    let prior_odds = prior_true / (1.0 - prior_true);
    (posterior_odds / prior_odds).clamp(1e-12, 1e12)
}

/// Infer posterior probabilities for indistinguishable protein groups.
///
/// `entries` contains one best identification per peptide as `(score, PEP,
/// protein ids)`. The discriminant score is deliberately ignored: after PEP
/// calibration it would count the same evidence twice.
pub fn infer(entries: &[(f64, f64, String)], params: &Params) -> InferenceResult {
    debug_assert!(params.validate().is_ok());

    let mut protein_index: HashMap<String, usize> = HashMap::new();
    let mut protein_names = Vec::new();
    let mut peptide_proteins: Vec<Vec<usize>> = Vec::new();
    let mut peptide_peps = Vec::new();

    for (_, pep, raw) in entries {
        let mut proteins = Vec::new();
        for name in split_proteins(raw) {
            let index = if let Some(&index) = protein_index.get(name) {
                index
            } else {
                let index = protein_names.len();
                protein_index.insert(name.to_string(), index);
                protein_names.push(name.to_string());
                index
            };
            proteins.push(index);
        }
        proteins.sort_unstable();
        proteins.dedup();
        if !proteins.is_empty() {
            peptide_proteins.push(proteins);
            peptide_peps.push(pep.clamp(1e-12, 1.0 - 1e-12));
        }
    }

    if protein_names.is_empty() {
        return InferenceResult {
            groups: Vec::new(),
            diagnostics: Diagnostics::default(),
        };
    }

    // Fido's protein clustering transformation: proteins attached to exactly
    // the same observed peptides are exchangeable and become one count-valued
    // variable (0..cluster_size present), with a binomial prior.
    let mut protein_peptides = vec![Vec::new(); protein_names.len()];
    for (peptide, proteins) in peptide_proteins.iter().enumerate() {
        for &protein in proteins {
            protein_peptides[protein].push(peptide);
        }
    }
    let mut cluster_of_neighborhood: HashMap<(bool, Vec<usize>), usize> = HashMap::new();
    let mut clusters: Vec<ProteinCluster> = Vec::new();
    let mut protein_cluster = vec![0; protein_names.len()];
    for protein in 0..protein_names.len() {
        let neighborhood = protein_peptides[protein].clone();
        let key = (
            is_decoy_protein(&protein_names[protein]),
            neighborhood.clone(),
        );
        let cluster = if let Some(&cluster) = cluster_of_neighborhood.get(&key) {
            cluster
        } else {
            let cluster = clusters.len();
            cluster_of_neighborhood.insert(key, cluster);
            clusters.push(ProteinCluster {
                proteins: Vec::new(),
                peptides: neighborhood,
                edges: Vec::new(),
                prior: Vec::new(),
            });
            cluster
        };
        protein_cluster[protein] = cluster;
        clusters[cluster]
            .proteins
            .push(protein_names[protein].clone());
    }
    for cluster in &mut clusters {
        cluster.proteins.sort();
        cluster.prior = binomial_prior(cluster.proteins.len(), params.gamma);
    }

    let mut factors: Vec<PeptideFactor> = peptide_proteins
        .iter()
        .enumerate()
        .map(|(peptide, proteins)| {
            let mut parents: Vec<usize> = proteins.iter().map(|&p| protein_cluster[p]).collect();
            parents.sort_unstable();
            parents.dedup();
            PeptideFactor {
                edges: Vec::with_capacity(parents.len()),
                likelihood_ratio: peptide_likelihood_ratio(
                    peptide_peps[peptide],
                    params.peptide_prior,
                ),
            }
        })
        .collect();

    let mut edges = Vec::new();
    for factor in 0..factors.len() {
        let mut parents: Vec<usize> = peptide_proteins[factor]
            .iter()
            .map(|&protein| protein_cluster[protein])
            .collect();
        parents.sort_unstable();
        parents.dedup();
        for cluster in parents {
            let edge = edges.len();
            let prior = clusters[cluster].prior.clone();
            let uniform = vec![1.0 / prior.len() as f64; prior.len()];
            edges.push(Edge {
                cluster,
                factor,
                variable_to_factor: prior,
                factor_to_variable: uniform,
            });
            clusters[cluster].edges.push(edge);
            factors[factor].edges.push(edge);
        }
    }

    let mut diagnostics = component_diagnostics(&clusters, &factors, &edges);
    let no_emit = 1.0 - params.alpha;
    let mut converged = false;
    for iteration in 1..=params.max_iter {
        let mut max_delta: f64 = 0.0;

        // Peptide noisy-OR factor -> protein-cluster messages. For a count k
        // of present proteins, P(peptide=true)=1-(1-beta)(1-alpha)^k.
        for factor in &factors {
            let moments: Vec<f64> = factor
                .edges
                .iter()
                .map(|&edge| {
                    edges[edge]
                        .variable_to_factor
                        .iter()
                        .enumerate()
                        .map(|(k, &prob)| prob * no_emit.powi(k as i32))
                        .sum::<f64>()
                })
                .collect();
            let all_moments: f64 = moments.iter().product();
            for (position, &edge_index) in factor.edges.iter().enumerate() {
                let other_moments = if moments[position] > 0.0 {
                    all_moments / moments[position]
                } else {
                    moments
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i != position)
                        .map(|(_, value)| *value)
                        .product()
                };
                let states = edges[edge_index].factor_to_variable.len();
                let mut message = Vec::with_capacity(states);
                for k in 0..states {
                    let no_emission = (1.0 - params.beta) * no_emit.powi(k as i32) * other_moments;
                    let emission = (1.0 - no_emission).clamp(0.0, 1.0);
                    let likelihood = 1.0 + (factor.likelihood_ratio - 1.0) * emission;
                    message.push(likelihood.max(1e-300));
                }
                normalize(&mut message);
                for (old, calculated) in
                    edges[edge_index].factor_to_variable.iter_mut().zip(message)
                {
                    let updated = *old + params.damping * (calculated - *old);
                    max_delta = max_delta.max((updated - *old).abs());
                    *old = updated;
                }
                normalize(&mut edges[edge_index].factor_to_variable);
            }
        }

        // Protein-cluster variable -> peptide factor messages. Work in log
        // space because a well-supported protein may have many peptide factors.
        for cluster in &clusters {
            let total_logs: Vec<f64> = cluster
                .prior
                .iter()
                .enumerate()
                .map(|(k, prior)| {
                    prior.max(1e-300).ln()
                        + cluster
                            .edges
                            .iter()
                            .map(|&edge| edges[edge].factor_to_variable[k].max(1e-300).ln())
                            .sum::<f64>()
                })
                .collect();
            for &edge_index in &cluster.edges {
                let message_logs: Vec<f64> = total_logs
                    .iter()
                    .enumerate()
                    .map(|(k, total)| {
                        *total - edges[edge_index].factor_to_variable[k].max(1e-300).ln()
                    })
                    .collect();
                let message = softmax(&message_logs);
                for (old, calculated) in
                    edges[edge_index].variable_to_factor.iter_mut().zip(message)
                {
                    let updated = *old + params.damping * (calculated - *old);
                    max_delta = max_delta.max((updated - *old).abs());
                    *old = updated;
                }
                normalize(&mut edges[edge_index].variable_to_factor);
            }
        }

        diagnostics.iterations = iteration;
        if max_delta < params.tolerance {
            converged = true;
            break;
        }
    }
    diagnostics.converged = converged;

    let mut groups: Vec<ProtGroup> = clusters
        .into_iter()
        .map(|cluster| {
            let logs: Vec<f64> = cluster
                .prior
                .iter()
                .enumerate()
                .map(|(k, prior)| {
                    prior.max(1e-300).ln()
                        + cluster
                            .edges
                            .iter()
                            .map(|&edge| edges[edge].factor_to_variable[k].max(1e-300).ln())
                            .sum::<f64>()
                })
                .collect();
            let posterior = softmax(&logs);
            let pep = posterior[0].clamp(0.0, 1.0);
            let is_decoy = cluster
                .proteins
                .iter()
                .all(|protein| is_decoy_protein(protein));
            ProtGroup {
                proteins: cluster.proteins,
                score: 1.0 - pep,
                qval: 1.0,
                // The Bayesian path does estimate a protein-level posterior, so
                // unlike picked-protein FDR it can fill this in.
                pep: Some(pep),
                n_peptides: cluster.peptides.len(),
                is_decoy,
                picked: true,
            }
        })
        .collect();
    assign_bayesian_qvalues(&mut groups);
    groups.sort_by(|a, b| {
        a.pep
            .partial_cmp(&b.pep)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.proteins.cmp(&b.proteins))
    });

    InferenceResult {
        groups,
        diagnostics,
    }
}

fn assign_bayesian_qvalues(groups: &mut [ProtGroup]) {
    for want_decoy in [false, true] {
        let mut order: Vec<usize> = groups
            .iter()
            .enumerate()
            .filter(|(_, group)| group.is_decoy == want_decoy)
            .map(|(index, _)| index)
            .collect();
        order.sort_by(|&a, &b| {
            groups[a]
                .pep
                .partial_cmp(&groups[b].pep)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut expected_errors = 0.0;
        let mut start = 0;
        while start < order.len() {
            let pep = groups[order[start]].pep.unwrap_or(1.0);
            let mut end = start + 1;
            while end < order.len() && groups[order[end]].pep.unwrap_or(1.0) == pep {
                end += 1;
            }
            expected_errors += pep * (end - start) as f64;
            let qvalue = expected_errors / end as f64;
            for &index in &order[start..end] {
                groups[index].qval = qvalue;
            }
            start = end;
        }
    }
}

fn component_diagnostics(
    clusters: &[ProteinCluster],
    factors: &[PeptideFactor],
    edges: &[Edge],
) -> Diagnostics {
    let nodes = clusters.len() + factors.len();
    let mut union_find = UnionFind::new(nodes);
    for edge in edges {
        union_find.union(edge.cluster, clusters.len() + edge.factor);
    }
    let mut counts: HashMap<usize, (usize, usize)> = HashMap::new();
    for node in 0..nodes {
        let root = union_find.find(node);
        counts.entry(root).or_default().0 += 1;
    }
    for edge in edges {
        let root = union_find.find(edge.cluster);
        counts.entry(root).or_default().1 += 1;
    }
    let tree_components = counts
        .values()
        .filter(|(nodes, edges)| *edges + 1 == *nodes)
        .count();
    Diagnostics {
        components: counts.len(),
        tree_components,
        loopy_components: counts.len() - tree_components,
        ..Diagnostics::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indistinguishable_proteins_form_one_group() {
        let entries = vec![(4.0, 0.001, "A B".to_string())];
        let params = Params::default();
        let result = infer(&entries, &params);
        assert_eq!(result.groups.len(), 1);
        assert_eq!(result.groups[0].proteins, vec!["A", "B"]);
        assert_eq!(result.groups[0].n_peptides, 1);
        assert!(result.groups[0].pep.unwrap() < 0.05);
        assert_eq!(result.diagnostics.tree_components, 1);

        let mut weights = [0.0; 3];
        for (k, weight) in weights.iter_mut().enumerate() {
            let binomial = [1.0, 2.0, 1.0][k]
                * params.gamma.powi(k as i32)
                * (1.0 - params.gamma).powi((2 - k) as i32);
            let emitted = 1.0 - (1.0 - params.beta) * (1.0 - params.alpha).powi(k as i32);
            let likelihood = 0.999 / params.peptide_prior * emitted
                + 0.001 / (1.0 - params.peptide_prior) * (1.0 - emitted);
            *weight = binomial * likelihood;
        }
        let exact_presence = 1.0 - weights[0] / weights.iter().sum::<f64>();
        assert!((result.groups[0].score - exact_presence).abs() < 1e-8);
    }

    #[test]
    fn unique_evidence_resolves_shared_peptide_ambiguity() {
        let entries = vec![
            (5.0, 0.01, "A B".to_string()),
            (8.0, 0.0001, "A".to_string()),
        ];
        let result = infer(&entries, &Params::default());
        assert_eq!(result.groups.len(), 2);
        let a = result
            .groups
            .iter()
            .find(|group| group.proteins == ["A"])
            .unwrap();
        let b = result
            .groups
            .iter()
            .find(|group| group.proteins == ["B"])
            .unwrap();
        assert!(
            a.pep < b.pep,
            "unique evidence should favor A: {:?} vs {:?}",
            a.pep,
            b.pep
        );
        assert_eq!(result.diagnostics.tree_components, 1);
    }

    #[test]
    fn weak_peptide_can_reduce_presence_below_prior() {
        let entries = vec![(0.0, 0.99, "A".to_string())];
        let result = infer(&entries, &Params::default());
        assert!(result.groups[0].score < 0.5);
    }

    #[test]
    fn tree_belief_propagation_matches_brute_force_marginal() {
        let params = Params::default();
        let entries = vec![
            (8.0, 0.001, "A".to_string()),
            (5.0, 0.02, "A B".to_string()),
        ];
        let result = infer(&entries, &params);
        let inferred_a = result
            .groups
            .iter()
            .find(|group| group.proteins == ["A"])
            .map(|group| group.score)
            .unwrap();

        let mut total = 0.0;
        let mut a_present = 0.0;
        for mask in 0..4usize {
            let a = mask & 1 != 0;
            let b = mask & 2 != 0;
            let active = a as usize + b as usize;
            let prior =
                params.gamma.powi(active as i32) * (1.0 - params.gamma).powi((2 - active) as i32);
            let factors = [(0.001, a as usize), (0.02, active)];
            let likelihood: f64 = factors
                .iter()
                .map(|&(pep, count)| {
                    let emitted =
                        1.0 - (1.0 - params.beta) * (1.0 - params.alpha).powi(count as i32);
                    (1.0 - pep) / params.peptide_prior * emitted
                        + pep / (1.0 - params.peptide_prior) * (1.0 - emitted)
                })
                .product();
            let weight = prior * likelihood;
            total += weight;
            if a {
                a_present += weight;
            }
        }
        let brute_a = a_present / total;
        assert!(
            (inferred_a - brute_a).abs() < 1e-8,
            "BP={inferred_a}, brute={brute_a}"
        );
        assert_eq!(result.diagnostics.tree_components, 1);
        assert!(result.diagnostics.converged);
    }

    #[test]
    fn qvalues_are_cumulative_expected_error() {
        let entries = vec![
            (9.0, 0.0001, "A".to_string()),
            (8.0, 0.01, "B".to_string()),
            (1.0, 0.8, "C".to_string()),
        ];
        let result = infer(&entries, &Params::default());
        let targets: Vec<&ProtGroup> = result
            .groups
            .iter()
            .filter(|group| !group.is_decoy)
            .collect();
        assert!(targets.windows(2).all(|pair| pair[0].qval <= pair[1].qval));
        assert!(targets.iter().all(|group| group.qval.is_finite()));
    }

    #[test]
    fn rejects_invalid_parameters() {
        let params = Params {
            gamma: 1.0,
            ..Params::default()
        };
        assert!(params.validate().is_err());
        let params = Params {
            alpha: f64::NAN,
            ..Params::default()
        };
        assert!(params.validate().is_err());
    }
}
