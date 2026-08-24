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

/// The canonical PIN matrix has 21 features plus bias. Keeping this fixed-size
/// dot product separate lets LLVM remove the loop without reassociating any
/// floating-point additions.
#[inline]
pub fn dot_22(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), 22);
    debug_assert_eq!(b.len(), 22);
    let mut s = 0.0;
    s += a[0] * b[0];
    s += a[1] * b[1];
    s += a[2] * b[2];
    s += a[3] * b[3];
    s += a[4] * b[4];
    s += a[5] * b[5];
    s += a[6] * b[6];
    s += a[7] * b[7];
    s += a[8] * b[8];
    s += a[9] * b[9];
    s += a[10] * b[10];
    s += a[11] * b[11];
    s += a[12] * b[12];
    s += a[13] * b[13];
    s += a[14] * b[14];
    s += a[15] * b[15];
    s += a[16] * b[16];
    s += a[17] * b[17];
    s += a[18] * b[18];
    s += a[19] * b[19];
    s += a[20] * b[20];
    s += a[21] * b[21];
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_dot_preserves_sequential_result() {
        let mut a = [0.0; 22];
        let mut b = [0.0; 22];
        for i in 0..22 {
            a[i] = (i as f64 - 9.0) / 7.0;
            b[i] = ((i * 17 % 23) as f64 - 11.0) * 1e-3;
        }
        assert_eq!(dot_22(&a, &b).to_bits(), dot(&a, &b).to_bits());
    }
}
