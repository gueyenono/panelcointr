//! Principal components of an `n x N` panel matrix and the factor-number
//! criteria of Bai and Ng (2002).
//!
//! The factors are `sqrt(n)` times the leading eigenvectors of `Y Y'`, so
//! that `f'f / n = I`, the loadings are `f'Y / n` and the residuals
//! `Y - f loadings'` (Bai and Ng 2004). The smaller Gram matrix is
//! factorised.

use super::eigen::sym_eigen;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Criterion {
    Pc1,
    Pc2,
    Pc3,
    Ic1,
    Ic2,
    Ic3,
    Bic3,
}

/// Eigendecomposition of the smaller Gram matrix of `Y` (`n x nn`).
#[derive(Clone, Debug, Default)]
pub struct GramEigen {
    pub n: usize,
    pub nn: usize,
    /// All eigenvalues of `Y'Y` (equivalently `Y Y'`), decreasing, clamped at 0.
    pub values: Vec<f64>,
    vectors: Vec<f64>,
    /// `true` if the `nn x nn` matrix `Y'Y` was factorised.
    by_units: bool,
}

/// Eigenvalues and vectors of the smaller of `Y'Y` and `Y Y'`.
pub fn gram_eigen(y: &[f64], n: usize, nn: usize) -> Result<GramEigen, String> {
    gram_eigen_by(y, n, nn, nn <= n)
}

fn gram_eigen_by(y: &[f64], n: usize, nn: usize, by_units: bool) -> Result<GramEigen, String> {
    let m = if by_units { nn } else { n };
    let mut g = vec![0.0; m * m];
    if by_units {
        for a in 0..nn {
            let ca = &y[a * n..(a + 1) * n];
            for b in 0..=a {
                let cb = &y[b * n..(b + 1) * n];
                let s: f64 = ca.iter().zip(cb).map(|(u, v)| u * v).sum();
                g[b * m + a] = s;
                g[a * m + b] = s;
            }
        }
    } else {
        for a in 0..n {
            for b in 0..=a {
                let s: f64 = (0..nn).map(|i| y[i * n + a] * y[i * n + b]).sum();
                g[b * m + a] = s;
                g[a * m + b] = s;
            }
        }
    }
    let eg = sym_eigen(&g, m)?;
    let total: f64 = y.iter().map(|v| v * v).sum();
    let values: Vec<f64> = eg.values.iter().map(|&v| if v > 0.0 { v } else { 0.0 }).collect();
    let sum: f64 = values.iter().sum();
    if (sum - total).abs() > 1e-8 * total.max(f64::MIN_POSITIVE) {
        return Err("the eigenvalues of the principal components do not add up to the total variation".into());
    }
    Ok(GramEigen { n, nn, values, vectors: eg.vectors, by_units })
}

/// Criterion values for `k = 0..=r_max` factors. `V(k)` is the residual
/// variance `(N n)^{-1} sum_{j > k} mu_j` and `sigma2 = V(r_max)`.
pub fn criteria(values: &[f64], n: usize, nn: usize, r_max: usize, crit: Criterion) -> Vec<f64> {
    let (nf, tf) = (nn as f64, n as f64);
    let nt = nf * tf;
    let total: f64 = values.iter().sum();
    let mut v = Vec::with_capacity(r_max + 1);
    let mut acc = total;
    for k in 0..=r_max {
        if k > 0 {
            acc -= values[k - 1];
        }
        v.push(acc.max(0.0) / nt);
    }
    let sigma2 = v[r_max];
    let c2 = nf.min(tf);
    let ratio = (nf + tf) / nt;
    (0..=r_max)
        .map(|k| {
            let kf = k as f64;
            match crit {
                Criterion::Pc1 => v[k] + kf * sigma2 * ratio * (nt / (nf + tf)).ln(),
                Criterion::Pc2 => v[k] + kf * sigma2 * ratio * c2.ln(),
                Criterion::Pc3 => v[k] + kf * sigma2 * c2.ln() / c2,
                Criterion::Ic1 => v[k].ln() + kf * ratio * (nt / (nf + tf)).ln(),
                Criterion::Ic2 => v[k].ln() + kf * ratio * c2.ln(),
                Criterion::Ic3 => v[k].ln() + kf * c2.ln() / c2,
                Criterion::Bic3 => v[k] + kf * sigma2 * (nf + tf - kf) * nt.ln() / nt,
            }
        })
        .collect()
}

/// Index of the smallest value, the first one among ties.
pub fn argmin(v: &[f64]) -> usize {
    let mut best = 0;
    for (k, &c) in v.iter().enumerate() {
        if c < v[best] {
            best = k;
        }
    }
    best
}

/// Factors, loadings and residuals of `Y` with `r` factors.
#[derive(Clone, Debug, Default)]
pub struct Extracted {
    /// `n x r`, column-major, `f'f / n = I`.
    pub f: Vec<f64>,
    /// `nn x r`, column-major.
    pub loadings: Vec<f64>,
    /// `n x nn`, column-major.
    pub resid: Vec<f64>,
}

pub fn extract(y: &[f64], ge: &GramEigen, r: usize) -> Result<Extracted, String> {
    let (n, nn) = (ge.n, ge.nn);
    if r > n.min(nn) {
        return Err(format!("cannot extract {r} factors from a {n} x {nn} matrix"));
    }
    let sn = (n as f64).sqrt();
    let m = if ge.by_units { nn } else { n };
    let mut f = vec![0.0; n * r];
    for j in 0..r {
        let mu = ge.values[j];
        if !(mu > 0.0) {
            return Err(format!("factor {} has zero variance", j + 1));
        }
        let v = &ge.vectors[j * m..(j + 1) * m];
        let col = &mut f[j * n..(j + 1) * n];
        if ge.by_units {
            // f = sqrt(n) Y v / sqrt(mu).
            let s = sn / mu.sqrt();
            for (i, &vi) in v.iter().enumerate() {
                let yi = &y[i * n..(i + 1) * n];
                for (c, &yv) in col.iter_mut().zip(yi) {
                    *c += yv * vi;
                }
            }
            col.iter_mut().for_each(|c| *c *= s);
        } else {
            for (c, &vi) in col.iter_mut().zip(v) {
                *c = sn * vi;
            }
        }
        // Sign: the element of largest magnitude is positive.
        let mut big = 0;
        for (t, c) in col.iter().enumerate() {
            if c.abs() > col[big].abs() {
                big = t;
            }
        }
        if col[big] < 0.0 {
            col.iter_mut().for_each(|c| *c = -*c);
        }
    }
    // Loadings f'Y / n and residuals.
    let mut loadings = vec![0.0; nn * r];
    let mut resid = y[..n * nn].to_vec();
    for i in 0..nn {
        let yi = &y[i * n..(i + 1) * n];
        for j in 0..r {
            let fj = &f[j * n..(j + 1) * n];
            let l: f64 = fj.iter().zip(yi).map(|(a, b)| a * b).sum::<f64>() / n as f64;
            loadings[j * nn + i] = l;
            for (e, &fv) in resid[i * n..(i + 1) * n].iter_mut().zip(fj) {
                *e -= fv * l;
            }
        }
    }
    Ok(Extracted { f, loadings, resid })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }

    fn factor_panel(n: usize, nn: usize, r: usize, seed: &mut u64) -> Vec<f64> {
        let f: Vec<f64> = (0..n * r).map(|_| lcg(seed)).collect();
        let l: Vec<f64> = (0..nn * r).map(|_| 3.0 * lcg(seed)).collect();
        let mut y = vec![0.0; n * nn];
        for i in 0..nn {
            for t in 0..n {
                y[i * n + t] = (0..r).map(|j| f[j * n + t] * l[j * nn + i]).sum::<f64>() + 0.1 * lcg(seed);
            }
        }
        y
    }

    #[test]
    fn both_routes_agree() {
        let mut seed = 3;
        let (n, nn, r) = (25, 12, 2);
        let y = factor_panel(n, nn, r, &mut seed);
        // Transposed problem: the same matrix seen through the other Gram.
        let ge = gram_eigen(&y, n, nn).unwrap();
        assert!(ge.by_units);
        let a = extract(&y, &ge, r).unwrap();
        let ge2 = gram_eigen_by(&y, n, nn, false).unwrap();
        for (u, v) in ge.values.iter().zip(&ge2.values) {
            assert!((u - v).abs() < 1e-10 * ge.values[0]);
        }
        let b = extract(&y, &ge2, r).unwrap();
        for (u, v) in a.f.iter().zip(&b.f) {
            assert!((u - v).abs() < 1e-10, "{u} vs {v}");
        }
        for (u, v) in a.resid.iter().zip(&b.resid) {
            assert!((u - v).abs() < 1e-10);
        }
        // f'f / n = I.
        for j in 0..r {
            for k in 0..r {
                let s: f64 = (0..n).map(|t| a.f[j * n + t] * a.f[k * n + t]).sum::<f64>() / n as f64;
                assert!((s - if j == k { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn bic3_finds_the_factors() {
        let mut seed = 17;
        let (n, nn, r) = (60, 30, 2);
        let y = factor_panel(n, nn, r, &mut seed);
        let ge = gram_eigen(&y, n, nn).unwrap();
        for crit in [Criterion::Bic3, Criterion::Ic1, Criterion::Pc2] {
            let c = criteria(&ge.values, n, nn, 6, crit);
            assert_eq!(argmin(&c), r);
        }
    }
}
