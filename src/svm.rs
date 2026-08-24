//! L2-regularized L2-loss (squared-hinge) linear SVM, per-sample weighted,
//! minimized in the primal with a truncated-Newton (Newton-CG / TRON-style) solver.
//! This is the same objective family as the reference L2-SVM-MFN routine.
//!
//! Objective:  f(w) = 1/2 ||w||^2 + sum_i C_i * max(0, 1 - y_i (w·x_i))^2
//!
//! Samples are passed as index lists into a shared feature matrix; a constant
//! bias feature is appended by the caller (so w has n_feat+1 entries).

pub struct Problem<'a> {
    pub x: &'a [f64], // row-major, rows * dim (dim includes the bias column)
    pub dim: usize,
    pub rows: &'a [usize], // indices of samples to train on
    pub y: &'a [f64],      // +1 / -1, aligned with rows
    pub c: &'a [f64],      // per-sample penalty, aligned with rows
    /// Optional active-feature mask. Excluded weights remain zero.
    pub feature_mask: Option<&'a [bool]>,
}

impl<'a> Problem<'a> {
    #[inline]
    fn xi(&self, k: usize) -> &[f64] {
        let r = self.rows[k];
        &self.x[r * self.dim..(r + 1) * self.dim]
    }

    fn wx(&self, w: &[f64], k: usize) -> f64 {
        match self.feature_mask {
            None => crate::simd::dot(w, self.xi(k)),
            Some(mask) => w
                .iter()
                .zip(self.xi(k))
                .zip(mask)
                .filter(|(_, active)| **active)
                .map(|((weight, value), _)| weight * value)
                .sum(),
        }
    }

    // value, and cache of (1 - y*wx) for active samples
    #[allow(clippy::needless_range_loop)]
    fn f_and_active(&self, w: &[f64], z: &mut [f64], active: &mut Vec<usize>) -> f64 {
        #[cfg(feature = "profiling")]
        let active_start = std::time::Instant::now();
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
        #[cfg(feature = "profiling")]
        crate::profile::record(
            "svm",
            "active_set_and_margin_scoring",
            active_start.elapsed(),
            Some(self.rows.len() as u64),
            Some(active.len() as u64),
        );
        f
    }

    // gradient g = w - 2 sum_{active} C_k z_k y_k x_k
    fn grad(&self, w: &[f64], z: &[f64], active: &[usize], g: &mut [f64]) {
        #[cfg(feature = "profiling")]
        let gradient_start = std::time::Instant::now();
        g[..self.dim].copy_from_slice(&w[..self.dim]);
        for &k in active {
            let coef = -2.0 * self.c[k] * z[k] * self.y[k];
            match self.feature_mask {
                None => crate::simd::axpy(&mut g[..self.dim], coef, self.xi(k)),
                Some(mask) => {
                    for j in 0..self.dim {
                        if mask[j] {
                            g[j] += coef * self.xi(k)[j];
                        }
                    }
                }
            }
        }
        #[cfg(feature = "profiling")]
        crate::profile::record(
            "svm",
            "gradient_computation",
            gradient_start.elapsed(),
            Some(active.len() as u64),
            None,
        );
    }

    // Explicitly form the Hessian H = I + 2 sum_{active} C_k x_k x_k^T (dim x dim, row-major).
    // dim is small (~22), so a single pass + direct solve beats matrix-free CG.
    fn hessian(&self, active: &[usize], h: &mut [f64]) {
        #[cfg(feature = "profiling")]
        let hessian_start = std::time::Instant::now();
        let dim = self.dim;
        for v in h.iter_mut() {
            *v = 0.0;
        }
        for &k in active {
            let xi = self.xi(k);
            let w = 2.0 * self.c[k];
            for a in 0..dim {
                if self.feature_mask.is_some_and(|mask| !mask[a]) {
                    continue;
                }
                let xa = w * xi[a];
                match self.feature_mask {
                    None => {
                        // H[a, a..dim] += xa * xi[a..dim]  (upper triangle)
                        crate::simd::axpy(&mut h[a * dim + a..a * dim + dim], xa, &xi[a..dim]);
                    }
                    Some(mask) => {
                        for b in a..dim {
                            if mask[b] {
                                h[a * dim + b] += xa * xi[b];
                            }
                        }
                    }
                }
            }
        }
        // add I and mirror upper -> lower
        for a in 0..dim {
            h[a * dim + a] += 1.0;
            for b in (a + 1)..dim {
                h[b * dim + a] = h[a * dim + b];
            }
        }
        #[cfg(feature = "profiling")]
        crate::profile::record(
            "svm",
            "hessian_construction",
            hessian_start.elapsed(),
            Some(active.len() as u64),
            Some((dim * dim * std::mem::size_of::<f64>()) as u64),
        );
    }
}

/// Solve SPD system H d = rhs by Cholesky (H, rhs consumed; result in `d`). Returns false if not PD.
fn cholesky_solve(h: &mut [f64], rhs: &[f64], d: &mut [f64], dim: usize) -> bool {
    #[cfg(feature = "profiling")]
    let factor_start = std::time::Instant::now();
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
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "svm",
        "cholesky_factorization",
        factor_start.elapsed(),
        Some(dim as u64),
        None,
    );
    #[cfg(feature = "profiling")]
    let solve_start = std::time::Instant::now();
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
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "svm",
        "linear_solve",
        solve_start.elapsed(),
        Some(dim as u64),
        None,
    );
    true
}

/// Train, warm-started from `w`. `max_newton` outer Newton steps.
pub fn train(p: &Problem, w: &mut [f64], max_newton: usize, tolerance: f64) {
    #[cfg(feature = "profiling")]
    let _svm_training =
        crate::profile::Scope::with_elements("svm", "svm_training_total", p.rows.len());
    let dim = p.dim;
    let n = p.rows.len();
    #[cfg(feature = "profiling")]
    let allocation_start = std::time::Instant::now();
    let mut z = vec![0.0f64; n];
    let mut active: Vec<usize> = Vec::with_capacity(n);
    let mut g = vec![0.0f64; dim];
    let mut d = vec![0.0f64; dim]; // newton direction
    let mut neg_g = vec![0.0f64; dim];
    let mut h = vec![0.0f64; dim * dim];
    let mut w_new = vec![0.0f64; dim];
    #[cfg(feature = "profiling")]
    {
        crate::profile::record(
            "svm",
            "allocation_and_buffer_initialization",
            allocation_start.elapsed(),
            Some(7),
            Some(
                (n * (std::mem::size_of::<f64>() + std::mem::size_of::<usize>())
                    + (4 * dim + dim * dim) * std::mem::size_of::<f64>()) as u64,
            ),
        );
        crate::profile::allocation_site(
            "svm::train work buffers",
            7,
            (n * (std::mem::size_of::<f64>() + std::mem::size_of::<usize>())
                + (4 * dim + dim * dim) * std::mem::size_of::<f64>()) as u64,
        );
    }

    let mut f = p.f_and_active(w, &mut z, &mut active);
    for newton_iteration in 0..max_newton {
        #[cfg(not(feature = "profiling"))]
        let _ = newton_iteration;
        #[cfg(feature = "profiling")]
        let _newton_context = crate::profile::context(None, None, None, Some(newton_iteration));
        #[cfg(feature = "profiling")]
        let _newton = crate::profile::Scope::new("svm", "newton_iteration_total");
        p.grad(w, &z, &active, &mut g);
        #[cfg(feature = "profiling")]
        let convergence_start = std::time::Instant::now();
        let gnorm2: f64 = g.iter().map(|v| v * v).sum();
        if gnorm2.sqrt() < tolerance {
            #[cfg(feature = "profiling")]
            crate::profile::record(
                "svm",
                "convergence_logic",
                convergence_start.elapsed(),
                Some(dim as u64),
                None,
            );
            break;
        }
        #[cfg(feature = "profiling")]
        crate::profile::record(
            "svm",
            "convergence_logic",
            convergence_start.elapsed(),
            Some(dim as u64),
            None,
        );
        // Newton step: form the (small) Hessian explicitly and Cholesky-solve H d = -g.
        p.hessian(&active, &mut h);
        #[cfg(feature = "profiling")]
        let update_start = std::time::Instant::now();
        for j in 0..dim {
            neg_g[j] = -g[j];
            d[j] = 0.0;
        }
        #[cfg(feature = "profiling")]
        crate::profile::record(
            "svm",
            "solver_buffer_update",
            update_start.elapsed(),
            Some(dim as u64),
            None,
        );
        if !cholesky_solve(&mut h, &neg_g, &mut d, dim) {
            // fall back to gradient descent direction if not PD (shouldn't happen: H >= I)
            d.copy_from_slice(&neg_g);
        }
        // Backtracking line search on f along d
        let gd: f64 = g.iter().zip(d.iter()).map(|(a, b)| a * b).sum();
        let mut step = 1.0;
        let mut ok = false;
        #[cfg(feature = "profiling")]
        let line_search_start = std::time::Instant::now();
        for _ls in 0..20 {
            #[cfg(feature = "profiling")]
            let weight_update_start = std::time::Instant::now();
            for j in 0..dim {
                w_new[j] = w[j] + step * d[j];
            }
            #[cfg(feature = "profiling")]
            crate::profile::record(
                "svm",
                "line_search_weight_update",
                weight_update_start.elapsed(),
                Some(dim as u64),
                None,
            );
            let f_new = p.f_and_active(&w_new, &mut z, &mut active);
            if f_new <= f + 1e-4 * step * gd {
                w.copy_from_slice(&w_new);
                f = f_new;
                ok = true;
                break;
            }
            step *= 0.5;
        }
        #[cfg(feature = "profiling")]
        crate::profile::record(
            "svm",
            "line_search_total",
            line_search_start.elapsed(),
            None,
            None,
        );
        if !ok {
            // recompute active set at current w (side effect on z/active) and stop
            let _ = p.f_and_active(w, &mut z, &mut active);
            break;
        }
    }
}
