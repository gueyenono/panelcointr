//! Augmented Dickey-Fuller regressions without deterministic terms,
//! `De_s = a0 e_{s-1} + sum_{j=1..k} a_j De_{s-j} + eps_s`, with the lag
//! order fixed or selected by the t-sig rule of Ng and Perron (1995) or by
//! an information criterion.

use super::linalg::Qr;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LagRule {
    /// Drop the last lag while `|t| < cut`, from `k_max` down to `k_min`.
    Tsig(f64),
    Aic,
    Bic,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LagSample {
    /// Compare candidates on the sample of the largest model.
    Common,
    /// Each candidate on its own sample.
    Varying,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdfSpec {
    pub k_min: usize,
    pub k_max: usize,
    pub rule: LagRule,
    pub sample: LagSample,
}

#[derive(Clone, Debug, Default)]
pub struct AdfFit {
    pub k: usize,
    /// Number of observations in the regression.
    pub n: usize,
    /// t-ratio of `a0`.
    pub t: f64,
    pub a0: f64,
    /// Residuals for `s = start..len`.
    pub resid: Vec<f64>,
    pub start: usize,
}

#[derive(Clone, Debug, Default)]
pub struct AdfWorkspace {
    z: Vec<f64>,
    dy: Vec<f64>,
    b: Vec<f64>,
    qr: Qr,
}

/// Fill the ADF design on rows `s = s0..len`: columns `e_{s-1}` (first if
/// `level_first`, otherwise last) and `De_{s-1..s-k}`.
fn fill(e: &[f64], k: usize, s0: usize, level_first: bool, ws: &mut AdfWorkspace) -> usize {
    let n = e.len() - s0;
    let m = k + 1;
    ws.z.clear();
    ws.z.resize(n * m, 0.0);
    ws.dy.clear();
    ws.dy.extend((s0..e.len()).map(|s| e[s] - e[s - 1]));
    let lev = if level_first { 0 } else { k };
    for (r, s) in (s0..e.len()).enumerate() {
        ws.z[lev * n + r] = e[s - 1];
        for j in 1..=k {
            let c = if level_first { j } else { j - 1 };
            ws.z[c * n + r] = e[s - j] - e[s - j - 1];
        }
    }
    n
}

fn too_short(len: usize, k: usize) -> String {
    format!("too few observations ({len}) for an ADF regression with {k} lags")
}

/// t-ratio of the last of the first `m` columns of a factorised design with
/// `n` rows (`Q'y` and `R` from the QR).
fn last_t(qr: &Qr, n: usize, m: usize) -> f64 {
    let df = n as f64 - m as f64;
    let s = (qr.rss_first(m) / df).sqrt();
    qr.qty(m - 1) * qr.r_diag(m - 1).signum() / s
}

fn ic(rss: f64, n: usize, m: usize, rule: LagRule) -> f64 {
    let nf = n as f64;
    let pen = match rule {
        LagRule::Aic => 2.0,
        _ => nf.ln(),
    };
    (rss / nf).ln() + m as f64 * pen / nf
}

/// Select the lag order.
pub fn select(e: &[f64], spec: &AdfSpec, ws: &mut AdfWorkspace) -> Result<usize, String> {
    let (k_min, k_max) = (spec.k_min, spec.k_max);
    if k_min == k_max {
        return Ok(k_max);
    }
    let len = e.len();
    if len < k_max + 3 {
        return Err(too_short(len, k_max));
    }
    match spec.sample {
        LagSample::Common => {
            let s0 = k_max + 1;
            let n = fill(e, k_max, s0, true, ws);
            let m = k_max + 1;
            if n <= m {
                return Err(too_short(len, k_max));
            }
            let AdfWorkspace { z, dy, qr, .. } = ws;
            qr.factor(z, n, m, dy).map_err(|_| "the ADF regression is rank deficient".to_string())?;
            match spec.rule {
                LagRule::Tsig(cut) => {
                    let mut k = k_max;
                    while k > k_min && last_t(qr, n, k + 1).abs() < cut {
                        k -= 1;
                    }
                    Ok(k)
                }
                rule => {
                    let mut best = k_min;
                    let mut best_ic = f64::INFINITY;
                    for k in k_min..=k_max {
                        let v = ic(qr.rss_first(k + 1), n, k + 1, rule);
                        if v < best_ic {
                            best_ic = v;
                            best = k;
                        }
                    }
                    Ok(best)
                }
            }
        }
        LagSample::Varying => {
            let fit_k = |k: usize, ws: &mut AdfWorkspace| -> Result<(usize, f64, f64), String> {
                let n = fill(e, k, k + 1, true, ws);
                if n <= k + 1 {
                    return Err(too_short(len, k));
                }
                let AdfWorkspace { z, dy, qr, .. } = ws;
                qr.factor(z, n, k + 1, dy).map_err(|_| "the ADF regression is rank deficient".to_string())?;
                Ok((n, qr.rss(), if k > 0 { last_t(qr, n, k + 1) } else { 0.0 }))
            };
            match spec.rule {
                LagRule::Tsig(cut) => {
                    let mut k = k_max;
                    while k > k_min {
                        let (_, _, t) = fit_k(k, ws)?;
                        if t.abs() >= cut {
                            break;
                        }
                        k -= 1;
                    }
                    Ok(k)
                }
                rule => {
                    let mut best = k_min;
                    let mut best_ic = f64::INFINITY;
                    for k in k_min..=k_max {
                        let (n, rss, _) = fit_k(k, ws)?;
                        let v = ic(rss, n, k + 1, rule);
                        if v < best_ic {
                            best_ic = v;
                            best = k;
                        }
                    }
                    Ok(best)
                }
            }
        }
    }
}

/// ADF regression with `k` lags on its own sample `s = k+1..len`.
pub fn fit(e: &[f64], k: usize, ws: &mut AdfWorkspace) -> Result<AdfFit, String> {
    let len = e.len();
    if len < k + 3 {
        return Err(too_short(len, k));
    }
    let s0 = k + 1;
    let n = fill(e, k, s0, false, ws);
    let m = k + 1;
    if n <= m {
        return Err(too_short(len, k));
    }
    let AdfWorkspace { z, dy, b, qr } = ws;
    qr.factor(z, n, m, dy).map_err(|_| "the ADF regression is rank deficient".to_string())?;
    qr.coef(b);
    let rss = qr.rss();
    let s = (rss / (n - m) as f64).sqrt();
    let r = qr.r_diag(m - 1);
    let a0 = qr.qty(m - 1) / r;
    let t = a0 * r.abs() / s;
    let resid = (0..n)
        .map(|row| dy[row] - (0..m).map(|c| z[c * n + row] * b[c]).sum::<f64>())
        .collect();
    Ok(AdfFit { k, n, t, a0, resid, start: s0 })
}

/// Select the lag order and fit.
pub fn adf(e: &[f64], spec: &AdfSpec, ws: &mut AdfWorkspace) -> Result<AdfFit, String> {
    let k = select(e, spec, ws)?;
    fit(e, k, ws)
}

/// Dickey-Fuller t-ratio without lags or deterministic terms, in O(n).
pub fn df_t(e: &[f64]) -> f64 {
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for s in 1..e.len() {
        let x = e[s - 1];
        let y = e[s] - x;
        sxx += x * x;
        sxy += x * y;
        syy += y * y;
    }
    let b = sxy / sxx;
    let n = (e.len() - 1) as f64;
    let s2 = (syy - b * sxy) / (n - 1.0);
    b / (s2 / sxx).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn walk(n: usize, seed: &mut u64) -> Vec<f64> {
        let mut e = vec![0.0; n];
        let mut acc = 0.0;
        for v in e.iter_mut() {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            acc += ((*seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5;
            *v = acc;
        }
        e
    }

    #[test]
    fn one_qr_tsig_matches_separate_fits() {
        let mut seed = 21;
        for _ in 0..20 {
            let e = walk(80, &mut seed);
            let k_max = 5;
            let s0 = k_max + 1;
            let mut ws = AdfWorkspace::default();
            let n = fill(&e, k_max, s0, true, &mut ws);
            let AdfWorkspace { z, dy, qr, .. } = &mut ws;
            qr.factor(z, n, k_max + 1, dy).unwrap();
            for k in 1..=k_max {
                let t_one = last_t(qr, n, k + 1);
                // Separate fit of the nested model on the same sample.
                let mut q2 = Qr::default();
                q2.factor(&z[..n * (k + 1)], n, k + 1, dy).unwrap();
                let mut b = Vec::new();
                q2.coef(&mut b);
                let s = (q2.rss() / (n - k - 1) as f64).sqrt();
                let t_sep = q2.qty(k) * q2.r_diag(k).signum() / s;
                assert!((t_one - t_sep).abs() < 1e-11 * t_sep.abs().max(1.0));
                let se = s / q2.r_diag(k).abs();
                assert!((b[k] / se - t_sep).abs() < 1e-9 * t_sep.abs().max(1.0));
            }
        }
    }

    #[test]
    fn df_t_matches_fit() {
        let mut seed = 4;
        let e = walk(60, &mut seed);
        let mut ws = AdfWorkspace::default();
        let f = fit(&e, 0, &mut ws).unwrap();
        assert!((f.t - df_t(&e)).abs() < 1e-11);
    }
}
