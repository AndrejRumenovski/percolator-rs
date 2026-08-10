//! L2-regularized L2-loss (squared-hinge) linear SVM, per-sample weighted,
//! minimized in the primal with a truncated-Newton (Newton-CG / TRON-style) solver.
//! This is the same objective family as the reference L2-SVM-MFN routine.
//!
//! Objective:  f(w) = 1/2 ||w||^2 + sum_i C_i * max(0, 1 - y_i (w·x_i))^2
//!
//! Samples are passed as index lists into a shared feature matrix; a constant
//! bias feature is appended by the caller (so w has n_feat+1 entries).

pub struct Problem<'a> {
    pub x: &'a [f64],   // row-major, rows * dim (dim includes the bias column)
    pub dim: usize,
    pub rows: &'a [usize], // indices of samples to train on
    pub y: &'a [f64],   // +1 / -1, aligned with rows
    pub c: &'a [f64],   // per-sample penalty, aligned with rows
}

impl<'a> Problem<'a> {
    #[inline]
    fn xi(&self, k: usize) -> &[f64] {
        let r = self.rows[k];
        &self.x[r * self.dim..(r + 1) * self.dim]
    }

    fn wx(&self, w: &[f64], k: usize) -> f64 {
        crate::simd::dot(w, self.xi(k))
    }

    // value, and cache of (1 - y*wx) for active samples
    fn f_and_active(&self, w: &[f64], z: &mut [f64], active: &mut Vec<usize>) -> f64 {
        let mut f = 0.0;
        for j in 0..self.dim {
            f += 0.5 * w[j] * w[j];
        }
        active.clear();
        for k in 0..self.rows.len() {
            let d = 1.0 - self.y[k] * self.wx(w, k);
            z[k] = d;
            if d > 0.0 {
                f += self.c[k] * d * d;
                active.push(k);
            }
        }
        f
    }

    // gradient g = w - 2 sum_{active} C_k z_k y_k x_k
    fn grad(&self, w: &[f64], z: &[f64], active: &[usize], g: &mut [f64]) {
        g[..self.dim].copy_from_slice(&w[..self.dim]);
        for &k in active {
            let coef = -2.0 * self.c[k] * z[k] * self.y[k];
            crate::simd::axpy(&mut g[..self.dim], coef, self.xi(k));
        }
    }

    // Explicitly form the Hessian H = I + 2 sum_{active} C_k x_k x_k^T (dim x dim, row-major).
    // dim is small (~22), so a single pass + direct solve beats matrix-free CG.
    fn hessian(&self, active: &[usize], h: &mut [f64]) {
        let dim = self.dim;
        for v in h.iter_mut() {
            *v = 0.0;
        }
        for &k in active {
            let xi = self.xi(k);
            let w = 2.0 * self.c[k];
            for a in 0..dim {
                let xa = w * xi[a];
                // H[a, a..dim] += xa * xi[a..dim]  (upper triangle)
                crate::simd::axpy(&mut h[a * dim + a..a * dim + dim], xa, &xi[a..dim]);
            }
        }
        // add I and mirror upper -> lower
        for a in 0..dim {
            h[a * dim + a] += 1.0;
            for b in (a + 1)..dim {
                h[b * dim + a] = h[a * dim + b];
            }
        }
    }
}

/// Solve SPD system H d = rhs by Cholesky (H, rhs consumed; result in `d`). Returns false if not PD.
fn cholesky_solve(h: &mut [f64], rhs: &[f64], d: &mut [f64], dim: usize) -> bool {
    // in-place Cholesky: H = L L^T (lower)
    for j in 0..dim {
        let mut sum = h[j * dim + j];
        for k in 0..j {
            sum -= h[j * dim + k] * h[j * dim + k];
        }
        if sum <= 1e-12 {
            return false;
        }
        let ljj = sum.sqrt();
        h[j * dim + j] = ljj;
        for i in (j + 1)..dim {
            let mut s = h[i * dim + j];
            for k in 0..j {
                s -= h[i * dim + k] * h[j * dim + k];
            }
            h[i * dim + j] = s / ljj;
        }
    }
    // forward solve L y = rhs
    for i in 0..dim {
        let mut s = rhs[i];
        for k in 0..i {
            s -= h[i * dim + k] * d[k];
        }
        d[i] = s / h[i * dim + i];
    }
    // back solve L^T d = y
    for i in (0..dim).rev() {
        let mut s = d[i];
        for k in (i + 1)..dim {
            s -= h[k * dim + i] * d[k];
        }
        d[i] = s / h[i * dim + i];
    }
    true
}

/// Train, warm-started from `w`. `max_newton` outer Newton steps.
pub fn train(p: &Problem, w: &mut [f64], max_newton: usize) {
    let dim = p.dim;
    let n = p.rows.len();
    let mut z = vec![0.0f64; n];
    let mut active: Vec<usize> = Vec::with_capacity(n);
    let mut g = vec![0.0f64; dim];
    let mut d = vec![0.0f64; dim]; // newton direction
    let mut neg_g = vec![0.0f64; dim];
    let mut h = vec![0.0f64; dim * dim];
    let mut w_new = vec![0.0f64; dim];

    let mut f = p.f_and_active(w, &mut z, &mut active);
    for _ in 0..max_newton {
        p.grad(w, &z, &active, &mut g);
        let gnorm2: f64 = g.iter().map(|v| v * v).sum();
        if gnorm2.sqrt() < 1e-5 {
            break;
        }
        // Newton step: form the (small) Hessian explicitly and Cholesky-solve H d = -g.
        p.hessian(&active, &mut h);
        for j in 0..dim {
            neg_g[j] = -g[j];
            d[j] = 0.0;
        }
        if !cholesky_solve(&mut h, &neg_g, &mut d, dim) {
            // fall back to gradient descent direction if not PD (shouldn't happen: H >= I)
            d.copy_from_slice(&neg_g);
        }
        // Backtracking line search on f along d
        let gd: f64 = g.iter().zip(d.iter()).map(|(a, b)| a * b).sum();
        let mut step = 1.0;
        let mut ok = false;
        for _ls in 0..20 {
            for j in 0..dim {
                w_new[j] = w[j] + step * d[j];
            }
            let f_new = p.f_and_active(&w_new, &mut z, &mut active);
            if f_new <= f + 1e-4 * step * gd {
                w.copy_from_slice(&w_new);
                f = f_new;
                ok = true;
                break;
            }
            step *= 0.5;
        }
        if !ok {
            // recompute active set at current w (side effect on z/active) and stop
            let _ = p.f_and_active(w, &mut z, &mut active);
            break;
        }
    }
}
