//! Small deterministic CPU MLP used as an optional PSM rescoring model.
//!
//! The network has one tanh hidden layer and a trainable linear skip connection.
//! Initializing the skip connection from Percolator's best-feature direction and
//! the hidden output weights to zero makes its initial score exactly match the
//! linear initialization used by the SVM.

pub struct Network {
    dim: usize,
    hidden: usize,
    skip: Vec<f64>,
    w1: Vec<f64>,
    w2: Vec<f64>,
    bias2: f64,
    moment1: Vec<f64>,
    moment2: Vec<f64>,
    step: u64,
    rng: Rng,
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 * (1.0 / ((1u64 << 53) as f64))
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

impl Network {
    pub fn new(dim: usize, hidden: usize, initial_linear: &[f64], seed: u64) -> Self {
        assert_eq!(dim, initial_linear.len());
        assert!(hidden > 0);
        let parameter_count = dim + hidden * dim + hidden + 1;
        let mut rng = Rng(seed.max(1));
        let scale = (6.0 / (dim + hidden) as f64).sqrt();
        let mut w1 = vec![0.0; hidden * dim];
        for weight in &mut w1 {
            *weight = (2.0 * rng.uniform() - 1.0) * scale;
        }
        Self {
            dim,
            hidden,
            skip: initial_linear.to_vec(),
            w1,
            // A zero output layer makes the initial network exactly reproduce
            // `initial_linear`; gradients reach w1 after the first Adam step.
            w2: vec![0.0; hidden],
            bias2: 0.0,
            moment1: vec![0.0; parameter_count],
            moment2: vec![0.0; parameter_count],
            step: 0,
            rng,
        }
    }

    #[inline]
    pub fn score(&self, row: &[f64]) -> f64 {
        debug_assert_eq!(row.len(), self.dim);
        let mut score = crate::simd::dot(&self.skip, row) + self.bias2;
        for h in 0..self.hidden {
            let start = h * self.dim;
            let activation = crate::simd::dot(&self.w1[start..start + self.dim], row).tanh();
            score += self.w2[h] * activation;
        }
        score
    }

    /// Optimize weighted binary cross-entropy with deterministic mini-batch Adam.
    #[allow(clippy::too_many_arguments)]
    pub fn train(
        &mut self,
        x: &[f64],
        rows: &[usize],
        labels: &[f64],
        weights: &[f64],
        epochs: usize,
        learning_rate: f64,
        l2: f64,
    ) {
        debug_assert_eq!(rows.len(), labels.len());
        debug_assert_eq!(rows.len(), weights.len());
        let mut order: Vec<usize> = (0..rows.len()).collect();
        let mut grad = vec![0.0; self.parameter_count()];
        let mut hidden_value = vec![0.0; self.hidden];

        for _ in 0..epochs {
            for i in (1..order.len()).rev() {
                let j = self.rng.below(i + 1);
                order.swap(i, j);
            }
            for batch in order.chunks(256) {
                grad.fill(0.0);
                let mut total_weight = 0.0;
                for &k in batch {
                    let row = &x[rows[k] * self.dim..(rows[k] + 1) * self.dim];
                    let mut logit = crate::simd::dot(&self.skip, row) + self.bias2;
                    for (h, hidden) in hidden_value.iter_mut().enumerate() {
                        let start = h * self.dim;
                        let value = crate::simd::dot(&self.w1[start..start + self.dim], row).tanh();
                        *hidden = value;
                        logit += self.w2[h] * value;
                    }
                    let target = if labels[k] > 0.0 { 1.0 } else { 0.0 };
                    let probability = sigmoid(logit);
                    let delta = weights[k] * (probability - target);
                    total_weight += weights[k];

                    for j in 0..self.dim {
                        grad[j] += delta * row[j];
                    }
                    let w1_offset = self.dim;
                    let w2_offset = w1_offset + self.hidden * self.dim;
                    for h in 0..self.hidden {
                        let hidden_delta =
                            delta * self.w2[h] * (1.0 - hidden_value[h] * hidden_value[h]);
                        let start = w1_offset + h * self.dim;
                        for j in 0..self.dim {
                            grad[start + j] += hidden_delta * row[j];
                        }
                        grad[w2_offset + h] += delta * hidden_value[h];
                    }
                    grad[w2_offset + self.hidden] += delta;
                }

                let normalizer = total_weight.max(1.0);
                for value in &mut grad {
                    *value /= normalizer;
                }
                self.add_l2_gradient(&mut grad, l2);
                clip_norm(&mut grad, 5.0);
                self.adam_step(&grad, learning_rate);
            }
        }
    }

    fn parameter_count(&self) -> usize {
        self.dim + self.hidden * self.dim + self.hidden + 1
    }

    fn add_l2_gradient(&self, grad: &mut [f64], l2: f64) {
        if l2 == 0.0 {
            return;
        }
        let mut offset = 0;
        for &value in &self.skip {
            grad[offset] += l2 * value;
            offset += 1;
        }
        for &value in &self.w1 {
            grad[offset] += l2 * value;
            offset += 1;
        }
        for &value in &self.w2 {
            grad[offset] += l2 * value;
            offset += 1;
        }
        // Do not regularize the output bias.
    }

    fn adam_step(&mut self, grad: &[f64], learning_rate: f64) {
        const BETA1: f64 = 0.9;
        const BETA2: f64 = 0.999;
        const EPSILON: f64 = 1e-8;
        self.step += 1;
        let correction1 = 1.0 - BETA1.powi(self.step.min(i32::MAX as u64) as i32);
        let correction2 = 1.0 - BETA2.powi(self.step.min(i32::MAX as u64) as i32);
        let mut offset = 0;

        for value in &mut self.skip {
            update_parameter(
                value,
                grad[offset],
                &mut self.moment1[offset],
                &mut self.moment2[offset],
                learning_rate,
                correction1,
                correction2,
                BETA1,
                BETA2,
                EPSILON,
            );
            offset += 1;
        }
        for value in &mut self.w1 {
            update_parameter(
                value,
                grad[offset],
                &mut self.moment1[offset],
                &mut self.moment2[offset],
                learning_rate,
                correction1,
                correction2,
                BETA1,
                BETA2,
                EPSILON,
            );
            offset += 1;
        }
        for value in &mut self.w2 {
            update_parameter(
                value,
                grad[offset],
                &mut self.moment1[offset],
                &mut self.moment2[offset],
                learning_rate,
                correction1,
                correction2,
                BETA1,
                BETA2,
                EPSILON,
            );
            offset += 1;
        }
        update_parameter(
            &mut self.bias2,
            grad[offset],
            &mut self.moment1[offset],
            &mut self.moment2[offset],
            learning_rate,
            correction1,
            correction2,
            BETA1,
            BETA2,
            EPSILON,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn update_parameter(
    value: &mut f64,
    gradient: f64,
    moment1: &mut f64,
    moment2: &mut f64,
    learning_rate: f64,
    correction1: f64,
    correction2: f64,
    beta1: f64,
    beta2: f64,
    epsilon: f64,
) {
    *moment1 = beta1 * *moment1 + (1.0 - beta1) * gradient;
    *moment2 = beta2 * *moment2 + (1.0 - beta2) * gradient * gradient;
    let adjusted1 = *moment1 / correction1;
    let adjusted2 = *moment2 / correction2;
    *value -= learning_rate * adjusted1 / (adjusted2.sqrt() + epsilon);
}

fn sigmoid(value: f64) -> f64 {
    if value >= 0.0 {
        1.0 / (1.0 + (-value).exp())
    } else {
        let exp = value.exp();
        exp / (1.0 + exp)
    }
}

fn clip_norm(gradient: &mut [f64], max_norm: f64) {
    let norm = gradient
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    if norm > max_norm {
        let scale = max_norm / norm;
        for value in gradient {
            *value *= scale;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn initial_score_matches_linear_direction() {
        let linear = vec![0.5, -2.0, 0.25];
        let network = Network::new(3, 4, &linear, 17);
        let row = vec![2.0, 0.5, 1.0];
        assert_eq!(network.score(&row), crate::simd::dot(&linear, &row));
    }

    #[test]
    fn learns_xor() {
        let mut x = Vec::new();
        let mut rows = Vec::new();
        let mut labels = Vec::new();
        for _ in 0..32 {
            for (a, b, label) in [
                (-1.0, -1.0, -1.0),
                (-1.0, 1.0, 1.0),
                (1.0, -1.0, 1.0),
                (1.0, 1.0, -1.0),
            ] {
                rows.push(rows.len());
                x.extend_from_slice(&[a, b, 1.0]);
                labels.push(label);
            }
        }
        let weights = vec![1.0; rows.len()];
        let mut network = Network::new(3, 8, &[0.0; 3], 3);
        network.train(&x, &rows, &labels, &weights, 200, 0.02, 1e-4);
        for (row, &label) in rows.iter().zip(&labels).take(4) {
            assert!(network.score(&x[row * 3..row * 3 + 3]) * label > 1.0);
        }
    }
}
