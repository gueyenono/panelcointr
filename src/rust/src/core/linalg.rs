//! Small dense linear algebra: Householder QR least squares and Cholesky.
//!
//! Matrices are column-major. The designs are small (n up to a few hundred,
//! k up to a few dozen), where plain loops beat a general library.

/// Rank tolerance of the QR. Columns are scaled to unit norm first, so
/// `|R_jj|` is the norm of column j after removing its projection on the
/// preceding columns, relative to its original norm. Below this value the
/// column counts as collinear: the criterion and default tolerance of
/// `lm.fit`.
pub const RANK_TOL: f64 = 1e-7;

/// Pivot tolerance of the Cholesky factorisation of a Gram matrix of
/// unit-norm columns: a pivot is the squared norm of the column after
/// removing the preceding ones, so this rejects angles below about 1e-5,
/// where squaring the condition number starts to cost accuracy. Callers
/// fall back to the QR.
pub const CHOL_TOL: f64 = 1e-10;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RankDeficient {
    /// Index of the first column found to be collinear with the preceding ones.
    pub column: usize,
}

/// Householder QR of an `n x k` design together with `Q'y`, reusable
/// across fits without reallocating.
#[derive(Clone, Debug, Default)]
pub struct Qr {
    n: usize,
    k: usize,
    // R in the upper triangle, Householder vectors (without their unit
    // leading element) below the diagonal; columns scaled to unit norm.
    a: Vec<f64>,
    scale: Vec<f64>,
    qty: Vec<f64>,
}

impl Qr {
    /// Factorise `z` (`n x k`, column-major) and compute `Q'y`.
    pub fn factor(&mut self, z: &[f64], n: usize, k: usize, y: &[f64]) -> Result<(), RankDeficient> {
        self.n = n;
        self.k = k;
        self.a.clear();
        self.a.extend_from_slice(&z[..n * k]);
        self.scale.clear();
        self.scale.resize(k, 0.0);
        self.qty.clear();
        self.qty.extend_from_slice(&y[..n]);
        if k > n {
            return Err(RankDeficient { column: n });
        }
        for j in 0..k {
            let col = &mut self.a[j * n..(j + 1) * n];
            let s = norm2(col);
            if !(s > 0.0 && s.is_finite()) {
                return Err(RankDeficient { column: j });
            }
            self.scale[j] = s;
            let inv = 1.0 / s;
            col.iter_mut().for_each(|v| *v *= inv);
        }
        for j in 0..k {
            let (head, tail) = self.a.split_at_mut((j + 1) * n);
            let colj = &mut head[j * n..];
            let alpha = colj[j];
            let xnorm = norm2(&colj[j + 1..]);
            let (beta, tau) = if xnorm == 0.0 {
                (alpha, 0.0)
            } else {
                let beta = -alpha.signum() * alpha.hypot(xnorm);
                let inv = 1.0 / (alpha - beta);
                colj[j + 1..].iter_mut().for_each(|v| *v *= inv);
                colj[j] = beta;
                (beta, (beta - alpha) / beta)
            };
            if !(beta.abs() >= RANK_TOL) {
                return Err(RankDeficient { column: j });
            }
            if tau != 0.0 {
                let v = &colj[j + 1..];
                for col in tail.chunks_exact_mut(n) {
                    reflect(col, j, v, tau);
                }
                reflect(&mut self.qty, j, v, tau);
            }
        }
        Ok(())
    }

    /// Least-squares coefficients, on the original scale of the columns.
    pub fn coef(&self, b: &mut Vec<f64>) {
        let (n, k) = (self.n, self.k);
        b.clear();
        b.extend_from_slice(&self.qty[..k]);
        for j in (0..k).rev() {
            let mut s = b[j];
            for (l, bl) in b.iter().enumerate().skip(j + 1) {
                s -= self.a[l * n + j] * bl;
            }
            b[j] = s / self.a[j * n + j];
        }
        for (bj, sj) in b.iter_mut().zip(&self.scale) {
            *bj /= sj;
        }
    }

    /// Residual sum of squares of the regression on the first `m` columns.
    pub fn rss_first(&self, m: usize) -> f64 {
        self.qty[m..self.n].iter().map(|v| v * v).sum()
    }

    /// Residual sum of squares of the full regression.
    pub fn rss(&self) -> f64 {
        self.rss_first(self.k)
    }

    /// `R_jj` on the original scale of the columns.
    pub fn r_diag(&self, j: usize) -> f64 {
        self.a[j * self.n + j] * self.scale[j]
    }

    /// Element `j` of `Q'y`.
    pub fn qty(&self, j: usize) -> f64 {
        self.qty[j]
    }
}

/// Apply `I - tau w w'`, `w = (1, v)`, to rows `j..` of `col`.
#[inline]
fn reflect(col: &mut [f64], j: usize, v: &[f64], tau: f64) {
    let rest = &mut col[j..];
    let (first, tail) = rest.split_first_mut().unwrap();
    let mut s = *first;
    for (vi, ci) in v.iter().zip(tail.iter()) {
        s += vi * ci;
    }
    s *= tau;
    *first -= s;
    for (vi, ci) in v.iter().zip(tail.iter_mut()) {
        *ci -= s * vi;
    }
}

pub fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// In-place lower Cholesky factor of the `m x m` column-major matrix `a`.
/// `L[i, j]` is stored at `a[j * m + i]`; the strict upper triangle is not
/// referenced. Returns `Err(j)` if pivot `j` falls below `CHOL_TOL` times
/// the corresponding diagonal element.
pub fn cholesky(a: &mut [f64], m: usize) -> Result<(), usize> {
    for j in 0..m {
        let diag = a[j * m + j];
        let mut d = diag;
        for l in 0..j {
            let v = a[l * m + j];
            d -= v * v;
        }
        if !(d > CHOL_TOL * diag) {
            return Err(j);
        }
        let ljj = d.sqrt();
        a[j * m + j] = ljj;
        for i in j + 1..m {
            let mut s = a[j * m + i];
            for l in 0..j {
                s -= a[l * m + i] * a[l * m + j];
            }
            a[j * m + i] = s / ljj;
        }
    }
    Ok(())
}

/// Solve `L v = b` in place, with `L` as returned by `cholesky`.
pub fn forward_solve(l: &[f64], m: usize, b: &mut [f64]) {
    for j in 0..m {
        let mut s = b[j];
        for i in 0..j {
            s -= l[i * m + j] * b[i];
        }
        b[j] = s / l[j * m + j];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Small deterministic pseudo-random numbers for tests.
    pub(crate) fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }

    fn naive_ls(z: &[f64], n: usize, k: usize, y: &[f64]) -> (Vec<f64>, f64) {
        // Normal equations solved by Gaussian elimination (well-conditioned designs only).
        let mut a = vec![0.0; k * k];
        let mut b = vec![0.0; k];
        for i in 0..k {
            for j in 0..k {
                a[i * k + j] = (0..n).map(|r| z[i * n + r] * z[j * n + r]).sum();
            }
            b[i] = (0..n).map(|r| z[i * n + r] * y[r]).sum();
        }
        for c in 0..k {
            for r in c + 1..k {
                let f = a[r * k + c] / a[c * k + c];
                for j in c..k {
                    a[r * k + j] -= f * a[c * k + j];
                }
                b[r] -= f * b[c];
            }
        }
        let mut x = vec![0.0; k];
        for c in (0..k).rev() {
            let s: f64 = (c + 1..k).map(|j| a[c * k + j] * x[j]).sum();
            x[c] = (b[c] - s) / a[c * k + c];
        }
        let rss = (0..n)
            .map(|r| {
                let f: f64 = (0..k).map(|j| z[j * n + r] * x[j]).sum();
                (y[r] - f).powi(2)
            })
            .sum();
        (x, rss)
    }

    #[test]
    fn qr_matches_normal_equations() {
        let mut seed = 42;
        let (n, k) = (40, 5);
        let mut z: Vec<f64> = (0..n * k).map(|_| lcg(&mut seed)).collect();
        for r in 0..n {
            z[r] = 1.0; // constant column
            z[n + r] = 1000.0 * (r as f64); // badly scaled trend
        }
        let y: Vec<f64> = (0..n).map(|_| lcg(&mut seed)).collect();
        let mut qr = Qr::default();
        qr.factor(&z, n, k, &y).unwrap();
        let mut b = Vec::new();
        qr.coef(&mut b);
        let (b0, rss0) = naive_ls(&z, n, k, &y);
        for (u, v) in b.iter().zip(&b0) {
            assert!((u - v).abs() <= 1e-9 * v.abs().max(1e-6), "{u} vs {v}");
        }
        assert!((qr.rss() - rss0).abs() <= 1e-10 * rss0);
        // Frisch-Waugh-Lovell through the last column: alpha = c / R_kk.
        assert!((qr.qty(k - 1) / qr.r_diag(k - 1) - b[k - 1]).abs() <= 1e-12 * b[k - 1].abs());
    }

    #[test]
    fn qr_detects_collinearity() {
        let mut seed = 7;
        let n = 20;
        let mut z: Vec<f64> = (0..n * 3).map(|_| lcg(&mut seed)).collect();
        for r in 0..n {
            z[2 * n + r] = 2.0 * z[r] - 3.0 * z[n + r];
        }
        let y = vec![1.0; n];
        let mut qr = Qr::default();
        assert_eq!(qr.factor(&z, n, 3, &y), Err(RankDeficient { column: 2 }));
    }

    #[test]
    fn partial_rss_of_nested_models() {
        let mut seed = 3;
        let (n, k) = (30, 4);
        let z: Vec<f64> = (0..n * k).map(|_| lcg(&mut seed)).collect();
        let y: Vec<f64> = (0..n).map(|_| lcg(&mut seed)).collect();
        let mut qr = Qr::default();
        qr.factor(&z, n, k, &y).unwrap();
        for m in 1..=k {
            let (_, rss) = naive_ls(&z[..n * m], n, m, &y);
            assert!((qr.rss_first(m) - rss).abs() <= 1e-11 * rss);
        }
    }

    #[test]
    fn cholesky_solves_gram_system() {
        let mut seed = 11;
        let (n, k) = (25, 4);
        let z: Vec<f64> = (0..n * k).map(|_| lcg(&mut seed)).collect();
        let y: Vec<f64> = (0..n).map(|_| lcg(&mut seed)).collect();
        let mut g = vec![0.0; k * k];
        let mut v = vec![0.0; k];
        for i in 0..k {
            for j in 0..k {
                g[j * k + i] = (0..n).map(|r| z[i * n + r] * z[j * n + r]).sum();
            }
            v[i] = (0..n).map(|r| z[i * n + r] * y[r]).sum();
        }
        cholesky(&mut g, k).unwrap();
        forward_solve(&g, k, &mut v);
        let yy: f64 = y.iter().map(|a| a * a).sum();
        let (_, rss) = naive_ls(&z, n, k, &y);
        let rss_chol = yy - v.iter().map(|a| a * a).sum::<f64>();
        assert!((rss_chol - rss).abs() <= 1e-10 * rss);
    }
}
