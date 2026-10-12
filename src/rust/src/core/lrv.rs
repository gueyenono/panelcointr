//! Long-run variance estimation.

/// Bartlett-kernel long-run variance with bandwidth `m`:
///
/// `(1/n) [ sum v_t^2 + 2 sum_{j=1..m} (1 - j/(m+1)) sum_t v_t v_{t-j} ]`,
///
/// summing over the available pairs, with no demeaning (Stata's `lrvar,
/// weighted nodemean`, used by `xtwest`). Returns NaN for an empty series.
pub fn bartlett_lrv(v: &[f64], m: usize) -> f64 {
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    let mut s: f64 = v.iter().map(|a| a * a).sum();
    let w_den = m as f64 + 1.0;
    for j in 1..=m.min(n - 1) {
        let c: f64 = v[j..].iter().zip(v).map(|(a, b)| a * b).sum();
        s += 2.0 * (1.0 - j as f64 / w_den) * c;
    }
    s / n as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_direct_formula() {
        let v = [0.3, -1.2, 0.8, 2.0, -0.5, 0.1, 0.7];
        let n = v.len() as f64;
        for m in 0..10usize {
            let mut s = 0.0;
            for j in -(m as isize)..=(m as isize) {
                let w = 1.0 - (j.unsigned_abs() as f64) / (m as f64 + 1.0);
                for t in 0..v.len() as isize {
                    let u = t - j;
                    if u >= 0 && (u as usize) < v.len() {
                        s += w * v[t as usize] * v[u as usize];
                    }
                }
            }
            assert!((bartlett_lrv(&v, m) - s / n).abs() < 1e-14);
        }
    }

    #[test]
    fn zero_bandwidth_is_the_uncentred_variance() {
        let v = [1.0, 2.0, 3.0];
        assert!((bartlett_lrv(&v, 0) - 14.0 / 3.0).abs() < 1e-15);
        assert!(bartlett_lrv(&[], 2).is_nan());
    }
}
