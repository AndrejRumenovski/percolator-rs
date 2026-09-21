//! Exact dot products and elementwise SIMD helpers, portable through `wide`.
//! Dot products retain the scalar left fold; SIMD lanes represent independent
//! rows or independent element updates, never partial sums of one dot product.

use wide::f64x4;

/// Dot product a·b, accumulated left to right from positive zero.
/// Keep each multiply and add separate and preserve their order.
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

/// Score four contiguous 22-column rows, preserving each scalar dot's left fold.
/// Each SIMD lane holds one complete row; objective reductions stay in the caller.
/// Non-finite results defer the batch to the caller's existing scalar loop.
#[inline(always)]
pub(crate) fn dot_22x4(w: &[f64], rows: &[f64]) -> Option<[f64; 4]> {
    debug_assert_eq!(w.len(), 22);
    debug_assert_eq!(rows.len(), 4 * 22);
    let w = &w[..22];
    let rows = &rows[..4 * 22];
    let mut sums = f64x4::splat(0.0);
    macro_rules! columns {
        ($($j:literal),* $(,)?) => {
            $(
                let values = f64x4::from([
                    rows[$j], rows[22 + $j], rows[44 + $j], rows[66 + $j],
                ]);
                sums += f64x4::splat(w[$j]) * values;
            )*
        };
    }
    columns!(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21);
    if sums.is_finite().all() {
        Some(sums.to_array())
    } else {
        None
    }
}

/// y += alpha * x  (elementwise, exact).
// Constant-length Hessian calls must inline so their loop bounds can fold.
#[inline(always)]
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

    fn check_four_rows(w: &[f64; 22], rows: &[f64; 88]) {
        let expected: [f64; 4] =
            std::array::from_fn(|lane| dot_22(w, &rows[lane * 22..(lane + 1) * 22]));
        let actual = dot_22x4(w, rows);
        if expected.iter().all(|score| score.is_finite()) {
            let actual = actual.expect("finite dots must use the SIMD batch");
            for (lane, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "row {lane} changed its scalar dot bits"
                );
            }
        } else {
            assert!(
                actual.is_none(),
                "a non-finite dot must defer the whole batch"
            );
        }
    }

    #[test]
    fn four_row_dot_matches_scalar_for_arbitrary_finite_bits() {
        let mut seed = 0x517cc1b727220a95_u64;
        let mut next_finite = |bounded: bool| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            let mut bits = seed;
            if bounded {
                bits = (bits & 0x800f_ffff_ffff_ffff) | ((992 + ((bits >> 52) & 63)) << 52);
            }
            if bits & 0x7ff0_0000_0000_0000 == 0x7ff0_0000_0000_0000 {
                bits ^= 0x0010_0000_0000_0000;
            }
            f64::from_bits(bits)
        };
        for case in 0..512 {
            let w = std::array::from_fn(|_| next_finite(case % 2 == 0));
            let rows = std::array::from_fn(|_| next_finite(case % 2 == 0));
            check_four_rows(&w, &rows);
        }
    }

    #[test]
    fn four_row_dot_preserves_cancellation_zero_subnormals_and_special_values() {
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
            1e150,
            -1e150,
            f64::MAX,
            -f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8_0000_0000_0042),
            f64::from_bits(0xfff8_0000_0000_1234),
        ];
        for offset in 0..values.len() {
            let w = std::array::from_fn(|j| values[(j + offset) % values.len()]);
            let rows = std::array::from_fn(|j| values[(j * 7 + offset) % values.len()]);
            check_four_rows(&w, &rows);
        }
        let w = [1.0; 22];
        let mut rows = [0.0; 88];
        rows[..6].copy_from_slice(&[1e16, 1.0, -1e16, 1.0, -0.0, 0.0]);
        rows[22..28].copy_from_slice(&[1e16, -1e16, 1.0, -1.0, 0.0, -0.0]);
        rows[44..50].copy_from_slice(&[
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
            f64::from_bits(1),
            f64::from_bits(1),
            -f64::from_bits(1),
            -0.0,
        ]);
        rows[66..].fill(-0.0);
        check_four_rows(&w, &rows);
        assert_eq!(dot_22x4(&w, &rows).unwrap()[0].to_bits(), 1.0_f64.to_bits());
        assert_eq!(dot_22x4(&w, &rows).unwrap()[3].to_bits(), 0.0_f64.to_bits());
    }

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
