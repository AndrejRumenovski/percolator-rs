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
    /// Rows have been copied into `x` in training order.
    pub packed_rows: bool,
    /// Optional active-feature mask. Excluded weights remain zero.
    pub feature_mask: Option<&'a [bool]>,
}

#[derive(Default)]
pub struct Workspace {
    // Raw scores for the current evaluation; gradients recompute their margins.
    z: Vec<f64>,
    row_scores_valid: bool,
    active: Vec<usize>,
    g: Vec<f64>,
    d: Vec<f64>,
    neg_g: Vec<f64>,
    h: Vec<f64>,
    w_new: Vec<f64>,
}

impl Workspace {
    /// Exact row scores at the weights returned by the most recent train call.
    /// Supplied initial scores and masked evaluations are never exposed here.
    pub(crate) fn final_row_scores(&self) -> Option<&[f64]> {
        self.row_scores_valid.then_some(self.z.as_slice())
    }
}

impl<'a> Problem<'a> {
    #[inline]
    fn xi(&self, k: usize) -> &[f64] {
        let r = if self.packed_rows { k } else { self.rows[k] };
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

    // Objective value, raw score cache for every sample, and active membership.
    #[allow(clippy::needless_range_loop)]
    fn f_and_active(&self, w: &[f64], z: &mut [f64], active: &mut Vec<usize>) -> f64 {
        #[cfg(feature = "profiling")]
        let active_start = std::time::Instant::now();
        let mut f = 0.0;
        for j in 0..self.dim {
            f += 0.5 * w[j] * w[j];
        }
        active.clear();
        if self.feature_mask.is_none() && self.dim == 22 {
            let mut first_scalar = 0;
            if self.packed_rows {
                while self.rows.len() - first_scalar >= 4 {
                    let Some(scores) = crate::simd::dot_22x4(
                        w,
                        &self.x[first_scalar * 22..(first_scalar + 4) * 22],
                    ) else {
                        // No batch state has changed. Score this batch and all
                        // remaining rows through the original scalar loop.
                        break;
                    };
                    // Keep margin comparisons and objective additions in the
                    // original row order after computing independent dots.
                    for (lane, score) in scores.into_iter().enumerate() {
                        let k = first_scalar + lane;
                        let d = 1.0 - self.y[k] * score;
                        z[k] = score;
                        if d > 0.0 {
                            f += self.c[k] * d * d;
                            active.push(k);
                        }
                    }
                    first_scalar += 4;
                }
            }
            for k in first_scalar..self.rows.len() {
                let score = crate::simd::dot_22(w, self.xi(k));
                let d = 1.0 - self.y[k] * score;
                z[k] = score;
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
            return f;
        }
        for k in 0..self.rows.len() {
            let score = self.wx(w, k);
            let d = 1.0 - self.y[k] * score;
            z[k] = score;
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

    // gradient g = w - 2 sum_{active} C_k (1 - y_k wx_k) y_k x_k
    fn grad(&self, w: &[f64], z: &[f64], active: &[usize], g: &mut [f64]) {
        #[cfg(feature = "profiling")]
        let gradient_start = std::time::Instant::now();
        g[..self.dim].copy_from_slice(&w[..self.dim]);
        for &k in active {
            let margin = 1.0 - self.y[k] * z[k];
            let coef = -2.0 * self.c[k] * margin * self.y[k];
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
        if dim == 22 && self.feature_mask.is_none() {
            for &k in active {
                hessian_add_22(h, 2.0 * self.c[k], self.xi(k));
            }
        } else {
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

// Fixed upper-triangle ranges let LLVM remove the short AXPY loops and their
// bounds checks for the canonical 21-feature matrix plus bias. Each cell uses
// the same multiply and add as the generic path, in the same active-row order.
#[inline]
fn hessian_add_22(h: &mut [f64], w: f64, xi: &[f64]) {
    let h = &mut h[..22 * 22];
    let xi = &xi[..22];
    macro_rules! accumulate_rows {
        ($($a:literal),* $(,)?) => {
            $(
                crate::simd::axpy(
                    &mut h[$a * 22 + $a..($a + 1) * 22],
                    w * xi[$a],
                    &xi[$a..22],
                );
            )*
        };
    }
    accumulate_rows!(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21);
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
pub fn train(
    p: &Problem,
    w: &mut [f64],
    initial_scores: &[f64],
    max_newton: usize,
    tolerance: f64,
    workspace: &mut Workspace,
) {
    #[cfg(feature = "profiling")]
    let _svm_training =
        crate::profile::Scope::with_elements("svm", "svm_training_total", p.rows.len());
    let dim = p.dim;
    let n = p.rows.len();
    workspace.row_scores_valid = false;
    #[cfg(feature = "profiling")]
    let allocation_start = std::time::Instant::now();
    #[cfg(feature = "profiling")]
    let old_capacities = [
        workspace.z.capacity(),
        workspace.active.capacity(),
        workspace.g.capacity(),
        workspace.d.capacity(),
        workspace.neg_g.capacity(),
        workspace.h.capacity(),
        workspace.w_new.capacity(),
    ];
    workspace.z.resize(n, 0.0);
    workspace.active.clear();
    if workspace.active.capacity() < n {
        workspace.active.reserve(n);
    }
    workspace.g.resize(dim, 0.0);
    workspace.d.resize(dim, 0.0);
    workspace.neg_g.resize(dim, 0.0);
    workspace.h.resize(dim * dim, 0.0);
    workspace.w_new.resize(dim, 0.0);
    let Workspace {
        z,
        row_scores_valid,
        active,
        g,
        d,
        neg_g,
        h,
        w_new,
    } = workspace;
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
            [
                z.capacity(),
                active.capacity(),
                g.capacity(),
                d.capacity(),
                neg_g.capacity(),
                h.capacity(),
                w_new.capacity(),
            ]
            .iter()
            .zip(old_capacities)
            .filter(|(new, old)| **new > *old)
            .count() as u64,
            ((z.capacity() - old_capacities[0]) * std::mem::size_of::<f64>()
                + (active.capacity() - old_capacities[1]) * std::mem::size_of::<usize>()
                + (g.capacity() - old_capacities[2]) * std::mem::size_of::<f64>()
                + (d.capacity() - old_capacities[3]) * std::mem::size_of::<f64>()
                + (neg_g.capacity() - old_capacities[4]) * std::mem::size_of::<f64>()
                + (h.capacity() - old_capacities[5]) * std::mem::size_of::<f64>()
                + (w_new.capacity() - old_capacities[6]) * std::mem::size_of::<f64>())
                as u64,
        );
    }

    debug_assert_eq!(initial_scores.len(), n);
    #[cfg(feature = "profiling")]
    let initial_objective_start = std::time::Instant::now();
    let mut f = if p.feature_mask.is_none() {
        let mut objective = 0.0;
        for &weight in w.iter().take(dim) {
            objective += 0.5 * weight * weight;
        }
        active.clear();
        for k in 0..n {
            let margin = 1.0 - p.y[k] * initial_scores[k];
            z[k] = initial_scores[k];
            if margin > 0.0 {
                objective += p.c[k] * margin * margin;
                active.push(k);
            }
        }
        objective
    } else {
        p.f_and_active(w, z, active)
    };
    #[cfg(feature = "profiling")]
    crate::profile::record(
        "svm",
        "initial_objective_and_active_set",
        initial_objective_start.elapsed(),
        Some(n as u64),
        Some(active.len() as u64),
    );
    for newton_iteration in 0..max_newton {
        #[cfg(not(feature = "profiling"))]
        let _ = newton_iteration;
        #[cfg(feature = "profiling")]
        let _newton_context = crate::profile::context(None, None, None, Some(newton_iteration));
        #[cfg(feature = "profiling")]
        let _newton = crate::profile::Scope::new("svm", "newton_iteration_total");
        p.grad(w, z, active, g);
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
        p.hessian(active, h);
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
        if !cholesky_solve(h, neg_g, d, dim) {
            // fall back to gradient descent direction if not PD (shouldn't happen: H >= I)
            d.copy_from_slice(neg_g);
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
            #[cfg(feature = "profiling")]
            let objective_start = std::time::Instant::now();
            let f_new = p.f_and_active(w_new, z, active);
            #[cfg(feature = "profiling")]
            crate::profile::record(
                "svm",
                "line_search_objective_evaluation",
                objective_start.elapsed(),
                Some(n as u64),
                Some(active.len() as u64),
            );
            if f_new <= f + 1e-4 * step * gd {
                w.copy_from_slice(w_new);
                *row_scores_valid = p.feature_mask.is_none();
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
            // Restore exact scores and active membership at current weights.
            let _ = p.f_and_active(w, z, active);
            *row_scores_valid = p.feature_mask.is_none();
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Problem;

    // Frozen pre-batching objective, including the canonical dot and masked
    // fallback. Equality covers the complete solver-visible return/state.
    #[allow(clippy::needless_range_loop)]
    fn frozen_objective(
        problem: &Problem<'_>,
        w: &[f64],
        z: &mut [f64],
        active: &mut Vec<usize>,
    ) -> f64 {
        let mut f = 0.0;
        for j in 0..problem.dim {
            f += 0.5 * w[j] * w[j];
        }
        active.clear();
        if problem.feature_mask.is_none() && problem.dim == 22 {
            for k in 0..problem.rows.len() {
                let d = 1.0 - problem.y[k] * crate::simd::dot_22(w, problem.xi(k));
                z[k] = d;
                if d > 0.0 {
                    f += problem.c[k] * d * d;
                    active.push(k);
                }
            }
            return f;
        }
        for k in 0..problem.rows.len() {
            let d = 1.0 - problem.y[k] * problem.wx(w, k);
            z[k] = d;
            if d > 0.0 {
                f += problem.c[k] * d * d;
                active.push(k);
            }
        }
        f
    }

    fn check_objective(problem: &Problem<'_>, w: &[f64]) {
        let mut expected_z = vec![f64::from_bits(0x7ff8_0000_0000_0042); problem.rows.len() + 3];
        let mut actual_z = expected_z.clone();
        let mut expected_active = vec![usize::MAX];
        let mut actual_active = expected_active.clone();
        let expected = frozen_objective(problem, w, &mut expected_z, &mut expected_active);
        let actual = problem.f_and_active(w, &mut actual_z, &mut actual_active);
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "objective bits changed"
        );
        assert_eq!(
            actual_active, expected_active,
            "active membership/order changed"
        );
        for row in 0..problem.rows.len() {
            let score = if problem.feature_mask.is_none() && problem.dim == 22 {
                crate::simd::dot_22(w, problem.xi(row))
            } else {
                problem.wx(w, row)
            };
            assert_eq!(
                actual_z[row].to_bits(),
                score.to_bits(),
                "raw score changed for row {row}"
            );
            let margin = 1.0 - problem.y[row] * actual_z[row];
            assert_eq!(
                margin.to_bits(),
                expected_z[row].to_bits(),
                "margin changed for row {row}"
            );
        }
        for (actual, expected) in actual_z[problem.rows.len()..]
            .iter()
            .zip(&expected_z[problem.rows.len()..])
        {
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "score cache overwrote its tail"
            );
        }
    }

    #[test]
    fn batched_objective_matches_frozen_for_remainders_layouts_and_masks() {
        for n in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 17, 64] {
            let physical_rows = n + 7;
            let rows: Vec<_> = (0..n).map(|row| row * 7 % physical_rows).collect();
            let y: Vec<_> = (0..n)
                .map(|row| if row % 3 == 0 { -1.0 } else { 1.0 })
                .collect();
            let c: Vec<_> = (0..n).map(|row| (row % 5) as f64 / 7.0).collect();
            for dim in [0, 1, 4, 21, 22, 23, 32] {
                let w: Vec<_> = (0..dim).map(|j| (j as f64 - 9.0) / 7.0).collect();
                let x: Vec<_> = (0..physical_rows * dim)
                    .map(|j| match j % 29 {
                        0 => 0.0,
                        1 => -0.0,
                        2 => f64::from_bits(1),
                        3 => -f64::from_bits(1),
                        _ => ((j * 17 % 97) as f64 - 48.0) / 23.0,
                    })
                    .collect();
                let all = vec![true; dim];
                let none = vec![false; dim];
                let alternating: Vec<_> = (0..dim).map(|j| j % 2 == 0).collect();
                for packed_rows in [false, true] {
                    for feature_mask in [
                        None,
                        Some(all.as_slice()),
                        Some(none.as_slice()),
                        Some(alternating.as_slice()),
                    ] {
                        let problem = Problem {
                            x: &x,
                            dim,
                            rows: &rows,
                            y: &y,
                            c: &c,
                            packed_rows,
                            feature_mask,
                        };
                        check_objective(&problem, &w);
                    }
                }
            }
        }
    }

    #[test]
    fn batched_objective_preserves_exact_hinge_boundaries() {
        let scores = [
            f64::from_bits(1.0_f64.to_bits() - 1),
            1.0,
            f64::from_bits(1.0_f64.to_bits() + 1),
            -1.0,
            -f64::from_bits(1.0_f64.to_bits() - 1),
            -f64::from_bits(1.0_f64.to_bits() + 1),
            0.0,
            -0.0,
            f64::from_bits(1),
        ];
        let mut x = vec![0.0; scores.len() * 22];
        for (row, &score) in scores.iter().enumerate() {
            x[row * 22] = score;
        }
        let mut w = [0.0; 22];
        w[0] = 1.0;
        let rows: Vec<_> = (0..scores.len()).collect();
        let y = [1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0];
        let c = [1.0; 9];
        let problem = Problem {
            x: &x,
            dim: 22,
            rows: &rows,
            y: &y,
            c: &c,
            packed_rows: true,
            feature_mask: None,
        };
        check_objective(&problem, &w);
        let mut z = vec![0.0; rows.len()];
        let mut active = Vec::new();
        problem.f_and_active(&w, &mut z, &mut active);
        assert_eq!(active, [0, 4, 6, 7, 8]);
    }

    #[test]
    fn batched_objective_preserves_finite_extrema_through_scalar_fallback() {
        const N: usize = 17;
        let rows: Vec<_> = (0..N).collect();
        let reversed: Vec<_> = (0..N).rev().collect();
        let y: Vec<_> = (0..N)
            .map(|row| if row % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let c = [0.5; N];
        let mut w = [0.0; 22];
        w[0] = 2.0;
        w[1] = 2.0;
        // All inputs remain finite. Some products overflow to either infinity,
        // and opposite overflowing products produce the original NaN margin.
        // Place them at each batch boundary and in the final scalar remainder.
        for extreme_row in [0, 3, 4, 7, 8, 12, 15, 16] {
            for pair in [
                [f64::MAX, 0.0],
                [-f64::MAX, 0.0],
                [f64::MAX, -f64::MAX],
                [f64::MAX, f64::MAX],
                [f64::MAX / 4.0, f64::MAX / 4.0],
                [f64::MIN_POSITIVE, -f64::MIN_POSITIVE],
                [f64::from_bits(1), -f64::from_bits(1)],
            ] {
                let mut x = vec![0.0; N * 22];
                for row in 0..N {
                    x[row * 22] = ((row % 5) as f64 - 2.0) / 8.0;
                    x[row * 22 + 1] = ((row % 7) as f64 - 3.0) / 16.0;
                }
                x[extreme_row * 22..extreme_row * 22 + 2].copy_from_slice(&pair);
                for (packed_rows, row_order) in [(true, &rows), (false, &reversed)] {
                    let problem = Problem {
                        x: &x,
                        dim: 22,
                        rows: row_order,
                        y: &y,
                        c: &c,
                        packed_rows,
                        feature_mask: None,
                    };
                    // Bit comparisons remain strict even for internally
                    // overflow-created NaNs from these finite input values.
                    check_objective(&problem, &w);
                }
            }
        }
    }

    #[test]
    fn final_row_scores_cover_inactive_rows_and_reset_for_each_train() {
        const N: usize = 9;
        let rows: Vec<_> = (0..N).collect();
        let y: Vec<_> = (0..N)
            .map(|row| if row % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let c = [1.0; N];
        let mut x = vec![0.0; N * 22];
        for row in 0..N {
            x[row * 22] = y[row] * if row == 0 { 10.0 } else { 1.0 };
        }
        let problem = Problem {
            x: &x,
            dim: 22,
            rows: &rows,
            y: &y,
            c: &c,
            packed_rows: true,
            feature_mask: None,
        };
        let mut workspace = super::Workspace::default();
        let mut w = [0.0; 22];
        super::train(&problem, &mut w, &[0.0; N], 1, 0.0, &mut workspace);
        let scores = workspace
            .final_row_scores()
            .expect("accepted evaluation must expose exact scores");
        for (row, score) in scores.iter().enumerate() {
            assert_eq!(
                score.to_bits(),
                crate::simd::dot_22(&w, problem.xi(row)).to_bits()
            );
        }
        assert!(
            !workspace.active.contains(&0),
            "fixture must include an inactive final row"
        );

        // Even plausible supplied scores are not certified actual dots at entry.
        super::train(&problem, &mut w, &[123.0; N], 0, 0.0, &mut workspace);
        assert!(workspace.final_row_scores().is_none());
        assert!(workspace
            .z
            .iter()
            .all(|score| score.to_bits() == 123.0_f64.to_bits()));
        let initial: Vec<_> = (0..N)
            .map(|row| crate::simd::dot_22(&w, problem.xi(row)))
            .collect();
        super::train(&problem, &mut w, &initial, 1, f64::INFINITY, &mut workspace);
        assert!(
            workspace.final_row_scores().is_none(),
            "immediate convergence must not certify supplied scores"
        );

        w.fill(0.0);
        super::train(&problem, &mut w, &[0.0; N], 1, 0.0, &mut workspace);
        assert!(workspace.final_row_scores().is_some());
        let mask = [true; 22];
        let masked = Problem {
            feature_mask: Some(&mask),
            ..problem
        };
        w.fill(0.0);
        super::train(&masked, &mut w, &[0.0; N], 1, 0.0, &mut workspace);
        assert!(
            workspace.final_row_scores().is_none(),
            "masked evaluations must remain private"
        );
    }

    #[test]
    fn failed_twenty_trial_search_restores_certified_raw_scores() {
        let x = [0.0; 44];
        let rows = [0, 1];
        let y = [1.0, -1.0];
        let c = [1.0, 1.0];
        let problem = Problem {
            x: &x,
            dim: 22,
            rows: &rows,
            y: &y,
            c: &c,
            packed_rows: true,
            feature_mask: None,
        };
        let mut w = [0.0; 22];
        w[0] = 1.0;
        let before = w;
        let mut workspace = super::Workspace::default();
        // These finite supplied scores imply zero hinge loss, although actual
        // features are zero. Every candidate has true loss >= 2 while the
        // initial Armijo threshold is below 0.5, forcing all twenty failures.
        super::train(&problem, &mut w, &y, 1, 0.0, &mut workspace);
        assert_eq!(w.map(f64::to_bits), before.map(f64::to_bits));
        assert_eq!(
            workspace.w_new[0].to_bits(),
            (1.0_f64 - 2.0_f64.powi(-19)).to_bits()
        );
        let scores = workspace
            .final_row_scores()
            .expect("restoration computes actual final scores");
        assert!(scores
            .iter()
            .all(|score| score.to_bits() == 0.0_f64.to_bits()));
        assert_eq!(workspace.active, [0, 1]);
    }

    // Frozen pre-specialization implementation, including SIMD/scalar AXPY
    // tails, feature masks, identity addition, and mirroring.
    fn frozen_hessian(problem: &Problem<'_>, active: &[usize], h: &mut [f64]) {
        let dim = problem.dim;
        for v in h.iter_mut() {
            *v = 0.0;
        }
        for &k in active {
            let xi = problem.xi(k);
            let w = 2.0 * problem.c[k];
            for a in 0..dim {
                if problem.feature_mask.is_some_and(|mask| !mask[a]) {
                    continue;
                }
                let xa = w * xi[a];
                match problem.feature_mask {
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
        for a in 0..dim {
            h[a * dim + a] += 1.0;
            for b in (a + 1)..dim {
                h[b * dim + a] = h[a * dim + b];
            }
        }
    }

    fn check_hessian(problem: &Problem<'_>, active: &[usize]) {
        // The entire output slice is cleared, even beyond the square matrix.
        let mut expected = vec![f64::NAN; problem.dim * problem.dim + 3];
        let mut actual = expected.clone();
        frozen_hessian(problem, active, &mut expected);
        problem.hessian(active, &mut actual);
        for (index, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "Hessian cell {index} differs: dim={}, packed={}, masked={}, active={active:?}",
                problem.dim,
                problem.packed_rows,
                problem.feature_mask.is_some(),
            );
        }
    }

    #[test]
    fn hessian_preserves_bits_for_layouts_masks_and_active_order() {
        const N: usize = 19;
        let rows: Vec<_> = (0..N).map(|row| row * 7 % N).collect();
        let y = vec![1.0; N];
        let c: Vec<_> = (0..N).map(|row| (row % 5) as f64 / 7.0).collect();
        let ascending: Vec<_> = (0..N).collect();
        let descending: Vec<_> = (0..N).rev().collect();
        let permuted: Vec<_> = (0..N).map(|row| row * 11 % N).collect();
        let repeated = [18, 0, 7, 1, 18, 3, 7];
        for dim in [0, 1, 2, 3, 4, 5, 21, 22, 23, 32] {
            let all = vec![true; dim];
            let none = vec![false; dim];
            let alternating: Vec<_> = (0..dim).map(|column| column % 2 == 0).collect();
            let bias_only: Vec<_> = (0..dim).map(|column| column + 1 == dim).collect();
            let mut seed = 0x517cc1b727220a95_u64;
            let x: Vec<_> = (0..N * dim)
                .map(|index| {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    match index % 17 {
                        0 => 0.0,
                        1 => -0.0,
                        _ => ((seed >> 32) as i32 as f64) / 1_000_000_007.0,
                    }
                })
                .collect();
            for packed_rows in [false, true] {
                for feature_mask in [
                    None,
                    Some(all.as_slice()),
                    Some(none.as_slice()),
                    Some(alternating.as_slice()),
                    Some(bias_only.as_slice()),
                ] {
                    let problem = Problem {
                        x: &x,
                        dim,
                        rows: &rows,
                        y: &y,
                        c: &c,
                        packed_rows,
                        feature_mask,
                    };
                    for active in [
                        &[][..],
                        &ascending[..1],
                        ascending.as_slice(),
                        descending.as_slice(),
                        permuted.as_slice(),
                        repeated.as_slice(),
                    ] {
                        check_hessian(&problem, active);
                    }
                }
            }
        }
    }

    #[test]
    fn fixed_hessian_preserves_bits_at_float_boundaries() {
        let values = [
            0.0,
            -0.0,
            f64::from_bits(1),
            -f64::from_bits(1),
            f64::from_bits(f64::MIN_POSITIVE.to_bits() - 1),
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
            f64::from_bits(1.0_f64.to_bits() - 1),
            1.0,
            f64::from_bits(1.0_f64.to_bits() + 1),
            -1.0,
            1e-150,
            -1e-150,
            1e150,
            -1e150,
            f64::MAX,
            -f64::MAX,
        ];
        let rows: Vec<_> = (0..values.len()).rev().collect();
        let y = vec![1.0; rows.len()];
        let ascending: Vec<_> = (0..rows.len()).collect();
        let x: Vec<_> = (0..rows.len() * 22)
            .map(|index| values[(index / 22 + index % 22 * 7) % values.len()])
            .collect();
        for penalty in [0.0, f64::MIN_POSITIVE, 0.1, 0.5, 1.0, 1e150] {
            let c = vec![penalty; rows.len()];
            for packed_rows in [false, true] {
                let problem = Problem {
                    x: &x,
                    dim: 22,
                    rows: &rows,
                    y: &y,
                    c: &c,
                    packed_rows,
                    feature_mask: None,
                };
                // Check single-row products before combinations can overflow,
                // then the original accumulation order in both directions.
                for row in 0..rows.len() {
                    check_hessian(&problem, &[row]);
                }
                check_hessian(&problem, &ascending);
                check_hessian(&problem, &rows);
            }
        }
    }
}
