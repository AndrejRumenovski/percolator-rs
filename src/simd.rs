//! SIMD helpers for the hot inner loops (portable via `wide`, compiles to AVX2/AVX-512
//! with `target-cpu=native`). `axpy` is exact (elementwise); `dot` reassociates the
//! summation across 4 lanes, which can differ from a sequential sum by ~1 ULP.

use wide::f64x4;

/// Dot product a·b. Kept as an exact sequential sum: on the short (~22-element)
/// feature vectors here, 4-lane reassociation gave no measurable speedup while
/// perturbing borderline q<0.01 yields, so we preserve exact summation order.
/// (The vectorized win lives in `axpy`, which is exact/elementwise — used to
/// accumulate the Hessian outer products, the actual SVM hot loop.)
#[inline]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut s = 0.0;
    for i in 0..a.len() {
        s += a[i] * b[i];
    }
    s
}

/// y += alpha * x  (elementwise, exact).
#[inline]
pub fn axpy(y: &mut [f64], alpha: f64, x: &[f64]) {
    debug_assert_eq!(y.len(), x.len());
    let n = y.len();
    let va = f64x4::splat(alpha);
    let mut i = 0;
    while i + 4 <= n {
        let vx = f64x4::from([x[i], x[i + 1], x[i + 2], x[i + 3]]);
        let vy = f64x4::from([y[i], y[i + 1], y[i + 2], y[i + 3]]);
        let r = (va * vx + vy).to_array();
        y[i] = r[0];
        y[i + 1] = r[1];
        y[i + 2] = r[2];
        y[i + 3] = r[3];
        i += 4;
    }
    while i < n {
        y[i] += alpha * x[i];
        i += 1;
    }
}
