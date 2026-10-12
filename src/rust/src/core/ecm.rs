//! Per-unit error-correction regressions: fit, composite error, the
//! semiparametric `α_i(1)` and lag/lead selection.

use super::design::{self, Col};
use super::linalg::{self, Qr, RankDeficient};
use super::lrv::bartlett_lrv;
use super::panel::Unit;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Criterion {
    /// `-2 log L + 2k` of the Gaussian regression (Stata's `estat ic` AIC).
    Aic,
    /// `-2 log L + k log n`.
    Bic,
    /// `log(RSS/n) + 2(p + q + dim d + 1) / (T - p_max - q_max)`, xtwest's
    /// `westerlund` option.
    Westerlund,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IcSample {
    /// All candidates on `S(p_max, q_max)` (Ng and Perron 2005).
    Common,
    /// Each candidate on its own sample `S(p, q)` (xtwest).
    Varying,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IcTies {
    /// Visit `p` then `q` ascending; strict improvement: ties go to the smaller model.
    Smaller,
    /// Visit `p` then `q` descending; strict improvement: ties go to the larger model (xtwest).
    Larger,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LrvDet {
    /// Residuals of `Δy` on `d_t` (paper §3.1).
    Regression,
    /// Demean in Case 3 only (xtwest).
    Implied,
    /// No adjustment (xtwest's `westerlund`).
    None,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum YSample {
    /// All available `Δy_t`, `t = 2..T`.
    All,
    /// The estimation sample `S(p, q)` (xtwest's `westerlund`).
    Estimation,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum USample {
    /// The estimation sample `S(p, q)`.
    Estimation,
    /// `S(p, q)` plus the last `q` periods (xtwest).
    Extended,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GridMethod {
    /// Nested Cholesky factors of the Gram matrix (falls back to QR).
    Gram,
    /// One QR per candidate (reference for tests).
    Qr,
}

/// Options of the per-unit computations.
#[derive(Clone, Debug)]
pub struct EcmOptions {
    /// Number of deterministic terms: 0 (none), 1 (constant), 2 (constant and trend).
    pub dd: usize,
    pub criterion: Criterion,
    pub ic_sample: IcSample,
    pub ic_ties: IcTies,
    pub lrv_det: LrvDet,
    pub y_sample: YSample,
    pub u_sample: USample,
    /// Demean `û` before its long-run variance when `d_t` has a constant.
    pub demean_u: bool,
}

/// Lag and lead orders to consider for one unit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LagSpec {
    pub p_min: usize,
    pub p_max: usize,
    pub q_min: usize,
    pub q_max: usize,
}

impl LagSpec {
    pub fn fixed(p: usize, q: usize) -> Self {
        LagSpec { p_min: p, p_max: p, q_min: q, q_max: q }
    }

    pub fn is_fixed(&self) -> bool {
        self.p_min == self.p_max && self.q_min == self.q_max
    }
}

/// Result of the ECM regression of one unit at given `(p, q)`.
#[derive(Clone, Copy, Debug, Default)]
pub struct UnitFit {
    pub p: usize,
    pub q: usize,
    /// Estimation sample size `T - p - q - 1`.
    pub n: usize,
    /// Number of regressors.
    pub k: usize,
    pub alpha: f64,
    /// `R_kk` of the QR with `y_{t-1}` last: `Σ ỹ²_{t-1} = r_kk²`.
    pub r_kk: f64,
    /// `(Q'Δy)_k`: `Σ ỹ_{t-1} Δỹ_t = r_kk c_k`.
    pub c_k: f64,
    pub rss: f64,
    pub omega_u2: f64,
    pub omega_y2: f64,
    pub alpha1: f64,
}

impl UnitFit {
    /// Regression standard error `sqrt(RSS / (n - k))`.
    pub fn s(&self) -> f64 {
        (self.rss / (self.n - self.k) as f64).sqrt()
    }

    /// Standard error of `α̂`, from `RSS / (n - k)` or, if `df` is false, `RSS / n`.
    pub fn se_alpha(&self, df: bool) -> f64 {
        let d = if df { self.n - self.k } else { self.n };
        (self.rss / d as f64).sqrt() / self.r_kk.abs()
    }
}

/// Reusable buffers for the per-unit computations.
#[derive(Clone, Debug, Default)]
pub struct Workspace {
    pub cols: Vec<Col>,
    pub z: Vec<f64>,
    pub dy: Vec<f64>,
    pub qr: Qr,
    /// Coefficients of the last fit, in the order of `cols`.
    pub b: Vec<f64>,
    v: Vec<f64>,
    order: Vec<Col>,
    idx: Vec<usize>,
    scale: Vec<f64>,
    gram: Vec<f64>,
    g: Vec<f64>,
    sub: Vec<f64>,
    subg: Vec<f64>,
    /// IC grid of the last selection, `ic[(p - p_min) * n_q + (q - q_min)]`.
    pub ic: Vec<f64>,
}

fn model(p: usize, q: usize) -> String {
    format!("{p} lag{} and {q} lead{}", if p == 1 { "" } else { "s" }, if q == 1 { "" } else { "s" })
}

fn rank_message(what: &str, cols: &[Col], e: RankDeficient, p: usize, q: usize, n: usize) -> String {
    if e.column >= cols.len() || cols.len() > n {
        format!(
            "{what} with {} has {} regressors but only {n} observations",
            model(p, q),
            cols.len()
        )
    } else {
        format!(
            "{what} with {} is rank deficient: {} is collinear with the preceding regressors",
            model(p, q),
            cols[e.column].label()
        )
    }
}

fn demean(v: &mut [f64]) {
    if v.is_empty() {
        return;
    }
    let m = v.iter().sum::<f64>() / v.len() as f64;
    v.iter_mut().for_each(|a| *a -= m);
}

/// Residuals of `v` on `(1, τ)` with `τ = t0 + 1, t0 + 2, …` (any affine
/// coding of time gives the same residuals).
fn detrend(v: &mut [f64], t0: usize) {
    let n = v.len();
    if n < 2 {
        demean(v);
        return;
    }
    let tbar = t0 as f64 + (n as f64 + 1.0) / 2.0;
    let vbar = v.iter().sum::<f64>() / n as f64;
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (r, a) in v.iter().enumerate() {
        let d = (t0 + r + 1) as f64 - tbar;
        sxy += d * (a - vbar);
        sxx += d * d;
    }
    let slope = sxy / sxx;
    for (r, a) in v.iter_mut().enumerate() {
        *a -= vbar + slope * ((t0 + r + 1) as f64 - tbar);
    }
}

/// Fit the ECM of one unit at `(p, q)` on `S(p, q)`, with bandwidth `m`
/// for the long-run variances. The coefficients are left in `ws.b`, in the
/// order of `ws.cols`.
pub fn fit(u: &Unit, p: usize, q: usize, m: usize, o: &EcmOptions, ws: &mut Workspace) -> Result<UnitFit, String> {
    let t_len = u.t_len;
    design::fit_columns(o.dd, u.k, p, q, &mut ws.cols);
    let k = ws.cols.len();
    if t_len < p + q + 2 || t_len - p - q - 1 <= k {
        return Err(format!(
            "the ECM regression with {} has {k} regressors but only {} observations",
            model(p, q),
            t_len.saturating_sub(p + q + 1)
        ));
    }
    let n = t_len - p - q - 1;
    let t0 = p + 1;
    design::fill(u, &ws.cols, t0, n, &mut ws.z);
    design::fill_dy(u, t0, n, &mut ws.dy);
    ws.qr
        .factor(&ws.z, n, k, &ws.dy)
        .map_err(|e| rank_message("the ECM regression", &ws.cols, e, p, q, n))?;
    ws.qr.coef(&mut ws.b);
    let rss = ws.qr.rss();
    let r_kk = ws.qr.r_diag(k - 1);
    let c_k = ws.qr.qty(k - 1);
    if !(rss > 0.0) {
        return Err(format!("the ECM regression with {} fits perfectly (zero residuals)", model(p, q)));
    }

    // Composite error û_t = Δy_t - δ'd_t - α y_{t-1} - λ'x_{t-1} - Σ α_j Δy_{t-j}.
    let t_end = match o.u_sample {
        USample::Estimation => t_len - 1 - q,
        USample::Extended => t_len - 1,
    };
    ws.v.clear();
    for t in t0..=t_end {
        let mut f = 0.0;
        for (col, b) in ws.cols.iter().zip(&ws.b) {
            if !matches!(col, Col::Dx(..)) {
                f += b * col.value(u, t);
            }
        }
        ws.v.push(u.dy(t) - f);
    }
    if o.demean_u && o.dd >= 1 {
        demean(&mut ws.v);
    }
    let omega_u2 = bartlett_lrv(&ws.v, m);

    // Δy with the deterministic terms removed.
    let (ys, ye) = match o.y_sample {
        YSample::All => (1, t_len - 1),
        YSample::Estimation => (t0, t_len - 1 - q),
    };
    ws.v.clear();
    ws.v.extend((ys..=ye).map(|t| u.dy(t)));
    match (o.lrv_det, o.dd) {
        (LrvDet::Regression, 1) | (LrvDet::Implied, 2) => demean(&mut ws.v),
        (LrvDet::Regression, 2) => detrend(&mut ws.v, ys),
        _ => {}
    }
    let omega_y2 = bartlett_lrv(&ws.v, m);

    if !(omega_y2.is_finite() && omega_y2 > 0.0) {
        return Err(format!(
            "the long-run variance of the differenced dependent variable is not positive ({omega_y2}) with {}",
            model(p, q)
        ));
    }
    if !(omega_u2.is_finite() && omega_u2 > 0.0) {
        return Err(format!(
            "the long-run variance of the composite error is not positive ({omega_u2}) with {}",
            model(p, q)
        ));
    }
    Ok(UnitFit {
        p,
        q,
        n,
        k,
        alpha: c_k / r_kk,
        r_kk,
        c_k,
        rss,
        omega_u2,
        omega_y2,
        alpha1: (omega_u2 / omega_y2).sqrt(),
    })
}

#[allow(clippy::too_many_arguments)]
fn ic_value(o: &EcmOptions, rss: f64, n: usize, k: usize, p: usize, q: usize, t_len: usize, spec: &LagSpec) -> f64 {
    let nf = n as f64;
    let m2ll = nf * (1.0 + (2.0 * std::f64::consts::PI).ln() + (rss / nf).ln());
    match o.criterion {
        Criterion::Aic => m2ll + 2.0 * k as f64,
        Criterion::Bic => m2ll + k as f64 * nf.ln(),
        Criterion::Westerlund => {
            (rss / nf).ln() + 2.0 * (p + q + o.dd + 1) as f64 / (t_len - spec.p_max - spec.q_max) as f64
        }
    }
}

/// Fill `ws.ic` with the information criterion of every candidate in `spec`.
pub fn ic_grid(u: &Unit, spec: &LagSpec, o: &EcmOptions, ws: &mut Workspace, method: GridMethod) -> Result<(), String> {
    let (np, nq) = (spec.p_max - spec.p_min + 1, spec.q_max - spec.q_min + 1);
    let t_len = u.t_len;
    let kmax = design::n_cols(o.dd, u.k, spec.p_max, spec.q_max);
    if t_len < spec.p_max + spec.q_max + 2 || t_len - spec.p_max - spec.q_max - 1 <= kmax {
        return Err(format!(
            "the ECM regression with {} has {kmax} regressors but only {} observations",
            model(spec.p_max, spec.q_max),
            t_len.saturating_sub(spec.p_max + spec.q_max + 1)
        ));
    }
    ws.ic.clear();
    ws.ic.resize(np * nq, f64::NAN);

    if o.ic_sample == IcSample::Varying || method == GridMethod::Qr {
        for p in spec.p_min..=spec.p_max {
            for q in spec.q_min..=spec.q_max {
                let (t0, n) = match o.ic_sample {
                    IcSample::Varying => (p + 1, t_len - p - q - 1),
                    IcSample::Common => (spec.p_max + 1, t_len - spec.p_max - spec.q_max - 1),
                };
                design::fit_columns(o.dd, u.k, p, q, &mut ws.cols);
                let k = ws.cols.len();
                design::fill(u, &ws.cols, t0, n, &mut ws.z);
                design::fill_dy(u, t0, n, &mut ws.dy);
                ws.qr
                    .factor(&ws.z, n, k, &ws.dy)
                    .map_err(|e| rank_message("the ECM regression", &ws.cols, e, p, q, n))?;
                let rss = ws.qr.rss();
                ws.ic[(p - spec.p_min) * nq + (q - spec.q_min)] = ic_value(o, rss, n, k, p, q, t_len, spec);
            }
        }
        return Ok(());
    }

    // Common sample: one pass over the data for the Gram matrix of the
    // largest design (unit-norm columns), then nested Cholesky factors.
    let t0 = spec.p_max + 1;
    let n = t_len - spec.p_max - spec.q_max - 1;
    design::fit_columns(o.dd, u.k, spec.p_max, spec.q_max, &mut ws.cols);
    design::fill(u, &ws.cols, t0, n, &mut ws.z);
    design::fill_dy(u, t0, n, &mut ws.dy);
    ws.scale.clear();
    for (j, col) in ws.z.chunks_exact_mut(n).enumerate() {
        let s = linalg::norm2(col);
        if !(s > 0.0 && s.is_finite()) {
            return Err(rank_message(
                "the ECM regression",
                &ws.cols,
                RankDeficient { column: j },
                spec.p_max,
                spec.q_max,
                n,
            ));
        }
        col.iter_mut().for_each(|v| *v /= s);
        ws.scale.push(s);
    }
    ws.gram.clear();
    ws.gram.resize(kmax * kmax, 0.0);
    ws.g.clear();
    for i in 0..kmax {
        let zi = &ws.z[i * n..(i + 1) * n];
        for j in 0..=i {
            let zj = &ws.z[j * n..(j + 1) * n];
            let s: f64 = zi.iter().zip(zj).map(|(a, b)| a * b).sum();
            ws.gram[j * kmax + i] = s;
            ws.gram[i * kmax + j] = s;
        }
        ws.g.push(zi.iter().zip(&ws.dy).map(|(a, b)| a * b).sum());
    }
    let yy: f64 = ws.dy.iter().map(|a| a * a).sum();

    // Nest the longer dimension inside the shorter one.
    let outer_q = nq <= np;
    let n_outer = if outer_q { nq } else { np };
    for oi in 0..n_outer {
        let (p_top, q_top) = if outer_q {
            design::nested_by_lags(o.dd, u.k, spec.p_max, spec.q_min + oi, &mut ws.order);
            (spec.p_max, spec.q_min + oi)
        } else {
            design::nested_by_leads(o.dd, u.k, spec.p_min + oi, spec.q_max, &mut ws.order);
            (spec.p_min + oi, spec.q_max)
        };
        let m = ws.order.len();
        ws.idx.clear();
        ws.idx
            .extend(ws.order.iter().map(|c| design::canonical_index(*c, o.dd, u.k, spec.p_max, spec.q_max)));
        ws.sub.clear();
        ws.sub.resize(m * m, 0.0);
        ws.subg.clear();
        for (b, &jb) in ws.idx.iter().enumerate() {
            for (a, &ja) in ws.idx.iter().enumerate() {
                ws.sub[b * m + a] = ws.gram[jb * kmax + ja];
            }
            ws.subg.push(ws.g[jb]);
        }
        let chol_ok = linalg::cholesky(&mut ws.sub, m).is_ok() && {
            linalg::forward_solve(&ws.sub, m, &mut ws.subg);
            yy - ws.subg.iter().map(|a| a * a).sum::<f64>() > 0.0
        };
        if chol_ok {
            // RSS of the leading block of length l: RSS(m) + Σ_{j >= l} v_j².
            let full = yy - ws.subg.iter().map(|a| a * a).sum::<f64>();
            let inner = if outer_q { spec.p_min..=spec.p_max } else { spec.q_min..=spec.q_max };
            for iv in inner {
                let (p, q) = if outer_q { (iv, q_top) } else { (p_top, iv) };
                let l = design::n_cols(o.dd, u.k, p, q);
                let rss = full + ws.subg[l..].iter().map(|a| a * a).sum::<f64>();
                ws.ic[(p - spec.p_min) * nq + (q - spec.q_min)] = ic_value(o, rss, n, l, p, q, t_len, spec);
            }
        } else {
            // Ill-conditioned block: QR of the nested design.
            let (cols, z, dy, qr, order) = (&mut ws.cols, &mut ws.z, &mut ws.dy, &mut ws.qr, &ws.order);
            cols.clear();
            cols.extend_from_slice(order);
            design::fill(u, cols, t0, n, z);
            design::fill_dy(u, t0, n, dy);
            qr.factor(z, n, m, dy)
                .map_err(|e| rank_message("the ECM regression", cols, e, p_top, q_top, n))?;
            let inner = if outer_q { spec.p_min..=spec.p_max } else { spec.q_min..=spec.q_max };
            for iv in inner {
                let (p, q) = if outer_q { (iv, q_top) } else { (p_top, iv) };
                let l = design::n_cols(o.dd, u.k, p, q);
                let rss = qr.rss_first(l);
                ws.ic[(p - spec.p_min) * nq + (q - spec.q_min)] = ic_value(o, rss, n, l, p, q, t_len, spec);
            }
            // The Gram data are still needed by the next outer value.
            design::fit_columns(o.dd, u.k, spec.p_max, spec.q_max, &mut ws.cols);
        }
    }
    Ok(())
}

/// Select `(p, q)` for one unit: the fixed orders, or the minimiser of the
/// information criterion over the grid in `spec`.
pub fn select(u: &Unit, spec: &LagSpec, o: &EcmOptions, ws: &mut Workspace, method: GridMethod) -> Result<(usize, usize), String> {
    if spec.is_fixed() {
        return Ok((spec.p_max, spec.q_max));
    }
    ic_grid(u, spec, o, ws, method)?;
    let nq = spec.q_max - spec.q_min + 1;
    let mut best = f64::INFINITY;
    let mut choice = None;
    let mut visit = |p: usize, q: usize| {
        let ic = ws.ic[(p - spec.p_min) * nq + (q - spec.q_min)];
        if ic < best {
            best = ic;
            choice = Some((p, q));
        }
    };
    match o.ic_ties {
        IcTies::Smaller => {
            for p in spec.p_min..=spec.p_max {
                for q in spec.q_min..=spec.q_max {
                    visit(p, q);
                }
            }
        }
        IcTies::Larger => {
            for p in (spec.p_min..=spec.p_max).rev() {
                for q in (spec.q_min..=spec.q_max).rev() {
                    visit(p, q);
                }
            }
        }
    }
    choice.ok_or_else(|| "the information criterion is not finite for any lag/lead order".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::panel::Panel;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }

    fn random_unit(t: usize, k: usize, seed: &mut u64) -> Panel {
        let mut x = vec![0.0; t * k];
        for j in 0..k {
            let mut acc = 0.0;
            for r in 0..t {
                acc += lcg(seed);
                x[j * t + r] = acc;
            }
        }
        let mut y = vec![0.0; t];
        let mut acc = 0.0;
        for r in 0..t {
            acc += lcg(seed) + 0.3 * if r > 0 { x[r] - x[r - 1] } else { 0.0 };
            y[r] = acc;
        }
        Panel::new(y, x, vec![0, t], k).unwrap()
    }

    fn options(dd: usize, criterion: Criterion) -> EcmOptions {
        EcmOptions {
            dd,
            criterion,
            ic_sample: IcSample::Common,
            ic_ties: IcTies::Smaller,
            lrv_det: LrvDet::Regression,
            y_sample: YSample::All,
            u_sample: USample::Estimation,
            demean_u: true,
        }
    }

    #[test]
    fn gram_grid_matches_qr_grid() {
        let mut seed = 99;
        for &(dd, k) in &[(0usize, 1usize), (1, 2), (2, 1), (2, 3)] {
            let panel = random_unit(80, k, &mut seed);
            let u = panel.unit(0);
            for &(p_min, p_max, q_min, q_max) in &[(0, 3, 0, 2), (1, 2, 0, 4), (0, 4, 1, 1)] {
                let spec = LagSpec { p_min, p_max, q_min, q_max };
                let o = options(dd, Criterion::Aic);
                let mut ws = Workspace::default();
                ic_grid(&u, &spec, &o, &mut ws, GridMethod::Gram).unwrap();
                let a = ws.ic.clone();
                ic_grid(&u, &spec, &o, &mut ws, GridMethod::Qr).unwrap();
                for (x, y) in a.iter().zip(&ws.ic) {
                    assert!((x - y).abs() <= 1e-9 * y.abs().max(1.0), "{x} vs {y}");
                }
            }
        }
    }

    #[test]
    fn fwl_quantities_are_consistent() {
        let mut seed = 5;
        let panel = random_unit(60, 2, &mut seed);
        let u = panel.unit(0);
        let mut ws = Workspace::default();
        let f = fit(&u, 2, 1, 3, &options(2, Criterion::Aic), &mut ws).unwrap();
        let k = ws.cols.len();
        assert!((f.alpha - ws.b[k - 1]).abs() <= 1e-12 * f.alpha.abs().max(1e-12));
        assert_eq!(f.n, 60 - 2 - 1 - 1);
        assert_eq!(f.k, design::n_cols(2, 2, 2, 1));
    }
}
