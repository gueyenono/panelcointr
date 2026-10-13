//! Tests on the estimated common factors (Banerjee and Carrion-i-Silvestre
//! 2015, Section 3; Bai and Ng 2004, Section 2.2): detrending of the
//! cumulated factors, the ADF test when there is one factor and the MQ_c
//! and MQ_f tests of the number of common stochastic trends otherwise.

use super::adf::{self, AdfFit, AdfSpec, AdfWorkspace};
use super::eigen::{gen_sym_eigen, sym_eigen};
use super::linalg::{cholesky, Qr};

/// Orthonormal basis (`n x p`, column-major) of the columns of `x`, by
/// modified Gram-Schmidt applied twice.
pub fn orthonormal_basis(x: &[f64], n: usize, p: usize) -> Result<Vec<f64>, String> {
    let mut q = x[..n * p].to_vec();
    for j in 0..p {
        let norm0 = q[j * n..(j + 1) * n].iter().map(|v| v * v).sum::<f64>().sqrt();
        for _ in 0..2 {
            for l in 0..j {
                let (done, rest) = q.split_at_mut(j * n);
                let ql = &done[l * n..(l + 1) * n];
                let qj = &mut rest[..n];
                let d: f64 = ql.iter().zip(qj.iter()).map(|(a, b)| a * b).sum();
                qj.iter_mut().zip(ql).for_each(|(b, a)| *b -= d * a);
            }
        }
        let qj = &mut q[j * n..(j + 1) * n];
        let norm = qj.iter().map(|v| v * v).sum::<f64>().sqrt();
        if !(norm > 1e-9 * norm0) {
            return Err("the deterministic terms of the factor regression are collinear".into());
        }
        qj.iter_mut().for_each(|v| *v /= norm);
    }
    Ok(q)
}

/// Replace `v` by its residual from the projection on the basis `q`.
pub fn project_out(q: &[f64], n: usize, v: &mut [f64]) {
    for ql in q.chunks_exact(n) {
        let d: f64 = ql.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
        v.iter_mut().zip(ql).for_each(|(b, a)| *b -= d * a);
    }
}

/// Deterministic terms of the factor regression on `n` rows: constant
/// (`trend = false`), or constant and trend plus a trend shift `DT_j` for
/// each break, where `breaks` are row indices of the last period of the old
/// regime.
pub fn deterministics(n: usize, trend: bool, breaks: &[usize]) -> Vec<f64> {
    let mut x = vec![1.0; n];
    if trend {
        x.extend((0..n).map(|s| (s + 1) as f64));
        for &b in breaks {
            x.extend((0..n).map(|s| if s > b { (s - b) as f64 } else { 0.0 }));
        }
    }
    x
}

#[derive(Clone, Debug, Default)]
pub struct MqRow {
    pub q: usize,
    pub mq_c: f64,
    pub mq_f: f64,
    /// VAR order of MQ_f.
    pub p: usize,
}

#[derive(Clone, Debug, Default)]
pub struct FactorTests {
    /// Detrended cumulated factors, `n x r`.
    pub fd: Vec<f64>,
    pub adf: Option<AdfFit>,
    pub mq: Vec<MqRow>,
}

/// Cumulate the factors `f` (`n x r`), detrend them on `det` (`n x p`) and
/// test them: ADF if `r = 1`, MQ for `q = r..1` otherwise.
pub fn factor_tests(
    f: &[f64],
    n: usize,
    r: usize,
    det: &[f64],
    adf_spec: &AdfSpec,
    bandwidth: usize,
    p_max: usize,
) -> Result<FactorTests, String> {
    let p = det.len() / n;
    let q = orthonormal_basis(det, n, p)?;
    let mut fd = vec![0.0; n * r];
    for j in 0..r {
        let col = &mut fd[j * n..(j + 1) * n];
        let mut acc = 0.0;
        for (c, &v) in col.iter_mut().zip(&f[j * n..(j + 1) * n]) {
            acc += v;
            *c = acc;
        }
        project_out(&q, n, col);
    }
    let mut out = FactorTests { fd, adf: None, mq: Vec::new() };
    if r == 1 {
        let mut ws = AdfWorkspace::default();
        out.adf = Some(adf::adf(&out.fd, adf_spec, &mut ws).map_err(|e| format!("factor ADF test: {e}"))?);
    } else if r > 1 {
        for qq in (1..=r).rev() {
            out.mq.push(mq(&out.fd, n, r, qq, bandwidth, p_max)?);
        }
    }
    Ok(out)
}

/// `sum_t a_t b_{t+lag}'` for rows of `n x q` matrices, `t = lo..hi`.
fn cross(a: &[f64], b: &[f64], n: usize, q: usize, rows: std::ops::Range<usize>, lag: usize) -> Vec<f64> {
    let mut s = vec![0.0; q * q];
    for t in rows {
        for i in 0..q {
            let ai = a[i * n + t];
            for j in 0..q {
                s[j * q + i] += ai * b[j * n + t + lag];
            }
        }
    }
    s
}

/// Smallest eigenvalue of `0.5 (num + num') den^{-1}`.
fn smallest(num: &[f64], den: &[f64], q: usize) -> Result<f64, String> {
    let mut sym = vec![0.0; q * q];
    for i in 0..q {
        for j in 0..q {
            sym[j * q + i] = 0.5 * (num[j * q + i] + num[i * q + j]);
        }
    }
    let eg = gen_sym_eigen(&sym, den, q).map_err(|_| "MQ test: the denominator matrix is not positive definite".to_string())?;
    Ok(*eg.values.last().unwrap())
}

/// Least squares of the `q` columns of `y` (`n x q`) on `x` (`n x m`):
/// coefficients (`m x q`) and residuals (`n x q`).
fn mols(x: &[f64], y: &[f64], n: usize, m: usize, q: usize) -> Result<(Vec<f64>, Vec<f64>), String> {
    if m == 0 {
        return Ok((Vec::new(), y.to_vec()));
    }
    let mut qr = Qr::default();
    let mut coef = Vec::with_capacity(m * q);
    let mut resid = Vec::with_capacity(n * q);
    let mut b = Vec::new();
    for j in 0..q {
        let yj = &y[j * n..(j + 1) * n];
        qr.factor(x, n, m, yj).map_err(|_| "MQ test: the VAR regression is rank deficient".to_string())?;
        qr.coef(&mut b);
        coef.extend_from_slice(&b);
        resid.extend((0..n).map(|t| yj[t] - (0..m).map(|c| x[c * n + t] * b[c]).sum::<f64>()));
    }
    Ok((coef, resid))
}

/// `ln det(E'E / n)` of an `n x q` residual matrix.
fn log_det_cov(e: &[f64], n: usize, q: usize) -> Result<f64, String> {
    let mut s = cross(e, e, n, q, 0..n, 0);
    s.iter_mut().for_each(|v| *v /= n as f64);
    cholesky(&mut s, q).map_err(|_| "MQ test: singular VAR residual covariance".to_string())?;
    Ok((0..q).map(|j| 2.0 * s[j * q + j].ln()).sum())
}

/// MQ_c and MQ_f statistics for `H0: r_1 = q` from the detrended factors
/// `fd` (`n x r`).
pub fn mq(fd: &[f64], n: usize, r: usize, q: usize, bandwidth: usize, p_max: usize) -> Result<MqRow, String> {
    let nf = n as f64;
    // beta_perp: eigenvectors of the q largest eigenvalues of n^{-2} sum F F'.
    let mut s = cross(fd, fd, n, r, 0..n, 0);
    s.iter_mut().for_each(|v| *v /= nf * nf);
    let eg = sym_eigen(&s, r)?;
    let mut y = vec![0.0; n * q];
    for j in 0..q {
        let b = eg.vector(j);
        for t in 0..n {
            y[j * n + t] = (0..r).map(|l| fd[l * n + t] * b[l]).sum();
        }
    }

    // MQ_c: VAR(1) without intercept and Bartlett correction.
    let sxx = cross(&y, &y, n, q, 0..n - 1, 0); // sum Y_{t-1} Y_{t-1}'
    let sxy = cross(&y, &y, n, q, 0..n - 1, 1); // sum Y_{t-1} Y_t'
    let ylag: Vec<f64> = (0..q).flat_map(|j| y[j * n..j * n + n - 1].to_vec()).collect();
    let ycur: Vec<f64> = (0..q).flat_map(|j| y[j * n + 1..(j + 1) * n].to_vec()).collect();
    let (_, xi) = mols(&ylag, &ycur, n - 1, q, q)?;
    let ne = n - 1;
    let mut sig1 = vec![0.0; q * q];
    for jl in 1..=bandwidth.min(ne.saturating_sub(1)) {
        let w = 1.0 - jl as f64 / (bandwidth as f64 + 1.0);
        let c = cross(&xi, &xi, ne, q, 0..ne - jl, jl); // sum xi_{t-j} xi_t'
        for (a, b) in sig1.iter_mut().zip(&c) {
            *a += w * b / nf;
        }
    }
    let mut num_c = vec![0.0; q * q];
    for i in 0..q {
        for j in 0..q {
            num_c[j * q + i] = sxy[j * q + i] - nf * sig1[j * q + i];
        }
    }
    let mq_c = nf * (smallest(&num_c, &sxx, q)? - 1.0);

    // MQ_f: VAR(p) in DY by BIC on the common sample, then filter Y.
    let dy: Vec<f64> = (0..q).flat_map(|j| (1..n).map(move |t| (j, t))).map(|(j, t)| y[j * n + t] - y[j * n + t - 1]).collect();
    let nd = n - 1; // rows of dy: index d <-> t = d + 1
    let lagged = |p: usize, d0: usize| -> (Vec<f64>, Vec<f64>, usize) {
        // Regress dy_d on dy_{d-1..d-p} for d = d0..nd.
        let rows = nd - d0;
        let mut x = Vec::with_capacity(rows * p * q);
        for l in 1..=p {
            for j in 0..q {
                x.extend((d0..nd).map(|d| dy[j * nd + d - l]));
            }
        }
        let yv: Vec<f64> = (0..q).flat_map(|j| dy[j * nd + d0..(j + 1) * nd].to_vec()).collect();
        (x, yv, rows)
    };
    let p_cap = p_max.min(nd.saturating_sub(q + 2) / (q + 1));
    let mut best = (0usize, f64::INFINITY);
    for p in 0..=p_cap {
        let (x, yv, rows) = lagged(p, p_cap);
        let (_, e) = mols(&x, &yv, rows, p * q, q)?;
        let bic = log_det_cov(&e, rows, q)? + (p * q * q) as f64 * (rows as f64).ln() / rows as f64;
        if bic < best.1 {
            best = (p, bic);
        }
    }
    let p = best.0;
    let (x, yv, rows) = lagged(p, p);
    let (coef, _) = mols(&x, &yv, rows, p * q, q)?;
    // Filter: ytil_t = Y_t - sum_l Pi_l Y_{t-l}, t = p..n-1.
    let nt = n - p;
    let mut yt = vec![0.0; nt * q];
    for j in 0..q {
        for t in p..n {
            let mut v = y[j * n + t];
            for l in 1..=p {
                for c in 0..q {
                    v -= coef[j * (p * q) + (l - 1) * q + c] * y[c * n + t - l];
                }
            }
            yt[j * nt + t - p] = v;
        }
    }
    let den_f = cross(&yt, &yt, nt, q, 0..nt - 1, 0);
    let num_f = cross(&yt, &yt, nt, q, 0..nt - 1, 1);
    let mq_f = nf * (smallest(&num_f, &den_f, q)? - 1.0);
    Ok(MqRow { q, mq_c, mq_f, p })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basis_projects_out_trend() {
        let n = 30;
        let det = deterministics(n, true, &[14]);
        let q = orthonormal_basis(&det, n, 3).unwrap();
        let mut v: Vec<f64> = (0..n).map(|s| 2.0 + 0.5 * s as f64 + if s > 14 { 3.0 * (s - 14) as f64 } else { 0.0 }).collect();
        project_out(&q, n, &mut v);
        assert!(v.iter().all(|x| x.abs() < 1e-10));
    }
}
