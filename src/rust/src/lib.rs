//! R bindings. The wrappers only convert and validate inputs; the
//! computations live in `core`, which has no extendr types.

// `!(x > 0.0)` is used on purpose: unlike `x <= 0.0`, it also rejects NaN.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use extendr_api::prelude::*;
use extendr_api::Result;

pub mod core;

use crate::core::adf::{self, AdfSpec, LagRule, LagSample};
use crate::core::banerjee::{self, Context, Factors, Options};
use crate::core::bootstrap::{self, BootState};
use crate::core::breaks::{self, Spec};
use crate::core::eigen;
use crate::core::factors::{self, Criterion as FactorCriterion};
use crate::core::mq;
use crate::core::nulldist::{self, Det, Kind};
use crate::core::rng::Rng;
use crate::core::design;
use crate::core::ecm::{self, Criterion, EcmOptions, GridMethod, IcSample, IcTies, LagSpec, LrvDet, USample, YSample};
use crate::core::linalg::Qr;
use crate::core::lrv::bartlett_lrv;
use crate::core::panel::Panel;
use crate::core::westerlund::{self, PanelLags, SeAlpha, Settings, SigmaDf, TNorm, TbarMode, Weights};

fn err<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::Other(msg.into()))
}

fn setting(list: &List, key: &str) -> Result<Robj> {
    list.iter()
        .find(|(n, _)| *n == key)
        .map(|(_, v)| v)
        .ok_or_else(|| Error::Other(format!("missing setting `{key}`")))
}

fn setting_str(list: &List, key: &str) -> Result<String> {
    setting(list, key)?
        .as_str()
        .map(str::to_string)
        .ok_or_else(|| Error::Other(format!("setting `{key}` must be a string")))
}

fn choose<T: Copy>(list: &List, key: &str, options: &[(&str, T)]) -> Result<T> {
    let v = setting_str(list, key)?;
    options
        .iter()
        .find(|(name, _)| *name == v)
        .map(|(_, t)| *t)
        .ok_or_else(|| Error::Other(format!("invalid value \"{v}\" for setting `{key}`")))
}

fn parse_settings(list: &List) -> Result<Settings> {
    let dd = choose(list, "deterministic", &[("none", 0usize), ("constant", 1), ("trend", 2)])?;
    let demean_u = setting(list, "demean_u")?
        .as_bool()
        .ok_or_else(|| Error::Other("setting `demean_u` must be TRUE or FALSE".into()))?;
    Ok(Settings {
        ecm: EcmOptions {
            dd,
            criterion: choose(
                list,
                "criterion",
                &[("aic", Criterion::Aic), ("bic", Criterion::Bic), ("westerlund", Criterion::Westerlund)],
            )?,
            ic_sample: choose(list, "ic_sample", &[("common", IcSample::Common), ("varying", IcSample::Varying)])?,
            ic_ties: choose(list, "ic_ties", &[("smaller", IcTies::Smaller), ("larger", IcTies::Larger)])?,
            lrv_det: choose(
                list,
                "lrv_deterministics",
                &[("regression", LrvDet::Regression), ("implied", LrvDet::Implied), ("none", LrvDet::None)],
            )?,
            y_sample: choose(list, "lrv_y_sample", &[("all", YSample::All), ("estimation", YSample::Estimation)])?,
            u_sample: choose(list, "u_sample", &[("estimation", USample::Estimation), ("extended", USample::Extended)])?,
            demean_u,
        },
        t_norm: choose(list, "t_norm", &[("effective", TNorm::Effective), ("df", TNorm::Df), ("T", TNorm::T)])?,
        panel_lags: choose(list, "panel_lags", &[("individual", PanelLags::Individual), ("mean", PanelLags::Mean)])?,
        se_alpha: choose(list, "se_alpha", &[("df", SeAlpha::Df), ("n", SeAlpha::N)])?,
        sigma: choose(list, "sigma_df", &[("own", SigmaDf::Own), ("tbar", SigmaDf::Tbar)])?,
        ga_weights: choose(list, "ga_weights", &[("units", Weights::Units), ("rows", Weights::Rows)])?,
        tbar: choose(list, "tbar", &[("units", TbarMode::Units), ("rows", TbarMode::Rows)])?,
    })
}

fn to_usize(v: &[i32], what: &str) -> Result<Vec<usize>> {
    v.iter()
        .map(|&a| {
            if a < 0 || a == i32::MIN {
                err(format!("`{what}` must be non-negative integers"))
            } else {
                Ok(a as usize)
            }
        })
        .collect()
}

fn build_panel(y: &[f64], x: &[f64], k: i32, offsets: &[i32]) -> Result<Panel> {
    if k < 1 {
        return err("at least one regressor is required");
    }
    let offsets = to_usize(offsets, "offsets")?;
    Panel::new(y.to_vec(), x.to_vec(), offsets, k as usize).map_err(Error::Other)
}

fn lag_specs(nu: usize, p_min: &[i32], p_max: &[i32], q_min: &[i32], q_max: &[i32]) -> Result<Vec<LagSpec>> {
    let (a, b, c, d) = (
        to_usize(p_min, "p_min")?,
        to_usize(p_max, "p_max")?,
        to_usize(q_min, "q_min")?,
        to_usize(q_max, "q_max")?,
    );
    if [a.len(), b.len(), c.len(), d.len()].iter().any(|&l| l != nu) {
        return err("lag and lead bounds must have one value per unit");
    }
    (0..nu)
        .map(|i| {
            if a[i] > b[i] || c[i] > d[i] {
                err("lag and lead ranges must satisfy min <= max")
            } else {
                Ok(LagSpec { p_min: a[i], p_max: b[i], q_min: c[i], q_max: d[i] })
            }
        })
        .collect()
}

/// Westerlund (2007) statistics and per-unit results.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn westerlund_rs(
    y: &[f64],
    x: &[f64],
    k: i32,
    offsets: &[i32],
    labels: Vec<String>,
    p_min: &[i32],
    p_max: &[i32],
    q_min: &[i32],
    q_max: &[i32],
    lrwindow: &[i32],
    settings: List,
) -> Result<List> {
    let panel = build_panel(y, x, k, offsets)?;
    let nu = panel.n_units();
    let specs = lag_specs(nu, p_min, p_max, q_min, q_max)?;
    let bw = to_usize(lrwindow, "lrwindow")?;
    if bw.len() != nu {
        return err("`lrwindow` must have one value per unit");
    }
    let s = parse_settings(&settings)?;
    let p_cap = specs.iter().map(|s| s.p_max).max().unwrap_or(0);
    let q_cap = specs.iter().map(|s| s.q_max).max().unwrap_or(0);
    let mut coef = Vec::new();
    let mut ws = ecm::Workspace::default();
    let out = westerlund::compute(&panel, &specs, &bw, &s, &mut ws, &labels, Some((&mut coef, p_cap, q_cap)))
        .map_err(Error::Other)?;
    let gi = |f: fn(&ecm::UnitFit) -> f64| out.group.iter().map(f).collect::<Vec<f64>>();
    let pi = |f: fn(&ecm::UnitFit) -> f64| out.panel.iter().map(f).collect::<Vec<f64>>();
    let mut cols = Vec::new();
    design::fit_columns(s.ecm.dd, panel.k, p_cap, q_cap, &mut cols);
    let coef_names: Vec<String> = cols.iter().map(|c| c.label()).collect();
    Ok(list!(
        statistics = out.stats.to_vec(),
        units = list!(
            lags = gi(|f| f.p as f64),
            leads = gi(|f| f.q as f64),
            n = gi(|f| f.n as f64),
            k = gi(|f| f.k as f64),
            alpha = gi(|f| f.alpha),
            se_alpha = out.se.clone(),
            alpha1 = gi(|f| f.alpha1),
            sigma = gi(|f| f.s()),
            rss = gi(|f| f.rss),
            omega_u2 = gi(|f| f.omega_u2),
            omega_y2 = gi(|f| f.omega_y2),
            t_tilde = out.t_tilde.clone()
        ),
        panel = list!(
            alpha = out.alpha,
            se = out.se_pooled,
            t_bar = out.t_bar,
            p_bar = out.p_bar as f64,
            q_bar = out.q_bar as f64,
            alpha1 = pi(|f| f.alpha1),
            sigma = out.sigma_panel.clone(),
            sum_y2 = pi(|f| f.r_kk * f.r_kk),
            sum_ydy = pi(|f| f.r_kk * f.c_k)
        ),
        coefficients = coef,
        coef_names = coef_names
    ))
}

/// Bootstrap state for `westerlund_boot_rs()`.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn westerlund_boot_prepare_rs(
    y: &[f64],
    x: &[f64],
    k: i32,
    offsets: &[i32],
    labels: Vec<String>,
    start: &[i32],
    lags: &[i32],
    leads: &[i32],
    p_min: &[i32],
    p_max: &[i32],
    q_min: &[i32],
    q_max: &[i32],
    lrwindow: &[i32],
    settings: List,
    reselect: bool,
    burn: i32,
) -> Result<ExternalPtr<BootState>> {
    let panel = build_panel(y, x, k, offsets)?;
    let nu = panel.n_units();
    let specs = lag_specs(nu, p_min, p_max, q_min, q_max)?;
    let bw = to_usize(lrwindow, "lrwindow")?;
    let start = to_usize(start, "start")?;
    let lags = to_usize(lags, "lags")?;
    let leads = to_usize(leads, "leads")?;
    if [bw.len(), start.len(), lags.len(), leads.len()].iter().any(|&l| l != nu) {
        return err("per-unit bootstrap inputs must have one value per unit");
    }
    if burn < 0 {
        return err("`burn` must be non-negative");
    }
    let sel: Vec<(usize, usize)> = lags.into_iter().zip(leads).collect();
    let s = parse_settings(&settings)?;
    let st = bootstrap::prepare(&panel, &start, &sel, &specs, &bw, &s, reselect, burn as usize, &labels)
        .map_err(Error::Other)?;
    Ok(ExternalPtr::new(st))
}

/// Dimensions of a bootstrap state.
/// @noRd
#[extendr]
fn westerlund_boot_info_rs(state: ExternalPtr<BootState>) -> List {
    let st: &BootState = &state;
    list!(n_periods = st.n_periods as f64, len = st.len as f64, t_span = st.t_span as f64)
}

/// Bootstrap statistics for the replicates whose indices are the columns
/// of `idx` (`len x n_rep`, 1-based); returns them row-major.
/// @noRd
#[extendr]
fn westerlund_boot_rs(state: ExternalPtr<BootState>, idx: &[i32], n_rep: i32, n_threads: i32) -> Result<Vec<f64>> {
    if n_rep < 0 || n_threads < 1 {
        return err("`n_rep` must be non-negative and `n_threads` positive");
    }
    let st: &BootState = &state;
    bootstrap::run(st, idx, n_rep as usize, n_threads as usize).map_err(Error::Other)
}

/// For tests: contents of a bootstrap state.
/// @noRd
#[extendr]
fn westerlund_boot_state_rs(state: ExternalPtr<BootState>) -> List {
    let st: &BootState = &state;
    let units: Vec<Robj> = st
        .units
        .iter()
        .map(|b| {
            list!(start = b.start as i32, t_len = b.t_len as i32, p = b.p as i32, q = b.q as i32, a = b.a.clone(), g = b.g.clone())
                .into()
        })
        .collect();
    list!(w = st.w.clone(), n_periods = st.n_periods as i32, burn = st.burn as i32, units = List::from_values(units))
}

/// For tests: the bootstrap panel generated from one replicate's indices.
/// @noRd
#[extendr]
fn westerlund_boot_panel_rs(state: ExternalPtr<BootState>, idx: &[i32]) -> Result<List> {
    let st: &BootState = &state;
    let mut ws = bootstrap::BootWorkspace::new(st);
    bootstrap::generate(st, idx, &mut ws).map_err(Error::Other)?;
    let p = ws.panel();
    Ok(list!(y = p.y.clone(), x = p.x.clone()))
}

/// Kernel for tests: least squares by Householder QR.
/// @noRd
#[extendr]
fn ols_rs(z: &[f64], ncol: i32, y: &[f64]) -> Result<List> {
    let n = y.len();
    let k = ncol.max(0) as usize;
    if z.len() != n * k {
        return err("dimensions of `z` and `y` do not match");
    }
    let mut qr = Qr::default();
    qr.factor(z, n, k, y)
        .map_err(|e| Error::Other(format!("rank deficient at column {}", e.column + 1)))?;
    let mut b = Vec::new();
    qr.coef(&mut b);
    Ok(list!(coef = b, rss = qr.rss(), r_last = qr.r_diag(k - 1), qty_last = qr.qty(k - 1)))
}

/// Kernel for tests: Bartlett long-run variance.
/// @noRd
#[extendr]
fn lrv_rs(v: &[f64], m: i32) -> f64 {
    bartlett_lrv(v, m.max(0) as usize)
}

/// Kernel for tests: information-criterion grid of one unit, by nested
/// Gram-Cholesky (`method = "gram"`) or one QR per candidate (`"qr"`).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn ecm_ic_rs(y: &[f64], x: &[f64], k: i32, p_min: i32, p_max: i32, q_min: i32, q_max: i32, settings: List, method: &str) -> Result<List> {
    let offsets = [0, y.len() as i32];
    let panel = build_panel(y, x, k, &offsets)?;
    let specs = lag_specs(1, &[p_min], &[p_max], &[q_min], &[q_max])?;
    let s = parse_settings(&settings)?;
    let method = match method {
        "gram" => GridMethod::Gram,
        "qr" => GridMethod::Qr,
        _ => return err("`method` must be \"gram\" or \"qr\""),
    };
    let mut ws = ecm::Workspace::default();
    let u = panel.unit(0);
    let (p, q) = if specs[0].is_fixed() {
        ecm::ic_grid(&u, &specs[0], &s.ecm, &mut ws, method).map_err(Error::Other)?;
        (specs[0].p_max, specs[0].q_max)
    } else {
        ecm::select(&u, &specs[0], &s.ecm, &mut ws, method).map_err(Error::Other)?
    };
    Ok(list!(ic = ws.ic.clone(), lags = p as i32, leads = q as i32))
}

// ---------------------------------------------------------------------------
// Banerjee and Carrion-i-Silvestre (2015)

fn setting_int(list: &List, key: &str) -> Result<usize> {
    let v = setting(list, key)?;
    let x = v
        .as_integer()
        .map(f64::from)
        .or_else(|| v.as_real())
        .ok_or_else(|| Error::Other(format!("setting `{key}` must be a number")))?;
    if !(x >= 0.0) || x.fract() != 0.0 {
        return err(format!("setting `{key}` must be a non-negative whole number"));
    }
    Ok(x as usize)
}

fn setting_real(list: &List, key: &str) -> Result<f64> {
    let v = setting(list, key)?;
    v.as_real()
        .or_else(|| v.as_integer().map(f64::from))
        .ok_or_else(|| Error::Other(format!("setting `{key}` must be a number")))
}

fn setting_bool(list: &List, key: &str) -> Result<bool> {
    setting(list, key)?
        .as_bool()
        .ok_or_else(|| Error::Other(format!("setting `{key}` must be TRUE or FALSE")))
}

const CRITERIA: [(&str, FactorCriterion); 7] = [
    ("bic3", FactorCriterion::Bic3),
    ("ic1", FactorCriterion::Ic1),
    ("ic2", FactorCriterion::Ic2),
    ("ic3", FactorCriterion::Ic3),
    ("pc1", FactorCriterion::Pc1),
    ("pc2", FactorCriterion::Pc2),
    ("pc3", FactorCriterion::Pc3),
];

fn parse_adf(list: &List) -> Result<AdfSpec> {
    let rule = match setting_str(list, "lag_criterion")?.as_str() {
        "tsig" => LagRule::Tsig(setting_real(list, "tsig_cut")?),
        "aic" => LagRule::Aic,
        "bic" => LagRule::Bic,
        v => return err(format!("invalid value \"{v}\" for setting `lag_criterion`")),
    };
    let k_min = setting_int(list, "lag_min")?;
    let k_max = setting_int(list, "lag_max")?;
    if k_min > k_max {
        return err("`lag_min` must not exceed `lag_max`");
    }
    Ok(AdfSpec {
        k_min,
        k_max,
        rule,
        sample: choose(list, "lag_sample", &[("common", LagSample::Common), ("varying", LagSample::Varying)])?,
    })
}

fn parse_banerjee(list: &List, k: usize) -> Result<Options> {
    let model = setting_int(list, "model")?;
    if !(1..=6).contains(&model) {
        return err("`model` must be between 1 and 6");
    }
    let factors = match setting_str(list, "factors")?.as_str() {
        "fixed" => Factors::Fixed(setting_int(list, "r_fixed")?),
        _ => Factors::Select(choose(list, "factors", &CRITERIA)?),
    };
    Ok(Options {
        spec: Spec { model: model as u8, k, dols_p: setting_int(list, "dols_p")?, dols_q: setting_int(list, "dols_q")? },
        factors,
        r_max: setting_int(list, "max_factors")?,
        adf: parse_adf(list)?,
        pca_scale: setting_bool(list, "pca_scale")?,
        mq_bandwidth: setting_int(list, "mq_bandwidth")?,
        mq_max_var: setting_int(list, "mq_max_var")?,
    })
}

/// Break positions from an `N x m_max` column-major matrix, -1 padded.
fn break_positions(pos: &[i32], nu: usize, m_max: usize) -> Result<Vec<Vec<usize>>> {
    if pos.len() != nu * m_max {
        return err("`pos` must have one row per unit");
    }
    (0..nu)
        .map(|i| {
            let mut v = Vec::new();
            for j in 0..m_max {
                let p = pos[j * nu + i];
                if p >= 0 {
                    v.push(p as usize);
                }
            }
            if v.windows(2).any(|w| w[1] <= w[0]) {
                return err("break positions must be increasing");
            }
            Ok(v)
        })
        .collect()
}

/// Banerjee and Carrion-i-Silvestre (2015) procedure for known break
/// positions (0-based, `N x m_max`, -1 padded).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn banerjee_rs(y: &[f64], x: &[f64], k: i32, offsets: &[i32], labels: Vec<String>, times: Vec<String>, pos: &[i32], m_max: i32, settings: List) -> Result<List> {
    let panel = build_panel(y, x, k, offsets)?;
    let nu = panel.n_units();
    let opts = parse_banerjee(&settings, panel.k)?;
    let pos = break_positions(pos, nu, m_max.max(0) as usize)?;
    let ctx = Context { panel: &panel, labels: &labels, times: &times };
    let res = banerjee::run_known(&ctx, &opts, &pos, true).map_err(Error::Other)?;
    let n = res.n;
    let ne = n - 1;
    // ADF residuals on the rows s = 1..n-1 of the cumulated series.
    let mut resid = vec![f64::NAN; ne * nu];
    for (i, f) in res.adf.iter().enumerate() {
        for (j, v) in f.resid.iter().enumerate() {
            resid[i * ne + f.start - 1 + j] = *v;
        }
    }
    let fs = &res.factor;
    let empty = || -> Robj { ().into() };
    let (fd, adf_f, mq_tab): (Robj, Robj, Robj) = match &res.tests {
        Some(t) => (
            t.fd.clone().into(),
            match &t.adf {
                Some(a) => list!(t = a.t, lags = a.k as i32, n = a.n as i32).into(),
                None => empty(),
            },
            if t.mq.is_empty() {
                empty()
            } else {
                list!(
                    q = t.mq.iter().map(|m| m.q as i32).collect::<Vec<i32>>(),
                    mq_c = t.mq.iter().map(|m| m.mq_c).collect::<Vec<f64>>(),
                    mq_f = t.mq.iter().map(|m| m.mq_f).collect::<Vec<f64>>(),
                    var_order = t.mq.iter().map(|m| m.p as i32).collect::<Vec<i32>>()
                )
                .into()
            },
        ),
        None => (empty(), empty(), empty()),
    };
    Ok(list!(
        t = res.adf.iter().map(|f| f.t).collect::<Vec<f64>>(),
        lags = res.adf.iter().map(|f| f.k as i32).collect::<Vec<i32>>(),
        nobs = res.adf.iter().map(|f| f.n as i32).collect::<Vec<i32>>(),
        tbar = res.tbar,
        r0 = res.r0 as i32,
        n = n as i32,
        eigenvalues = fs.values.clone(),
        criteria = fs.criteria.clone(),
        r = fs.r as i32,
        f = fs.ext.f.clone(),
        loadings = fs.ext.loadings.clone(),
        resid = resid,
        fd = fd,
        factor_adf = adf_f,
        mq = mq_tab
    ))
}

/// Iterative estimation of unit-specific break dates.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn banerjee_breaks_rs(y: &[f64], x: &[f64], k: i32, offsets: &[i32], labels: Vec<String>, times: Vec<String>, m: i32, h: i32, settings: List, n_threads: i32) -> Result<List> {
    let panel = build_panel(y, x, k, offsets)?;
    let opts = parse_banerjee(&settings, panel.k)?;
    if !(1..=2).contains(&m) || h < 1 || n_threads < 1 {
        return err("`m` must be 1 or 2, `h` and `n_threads` positive");
    }
    let iter_max = setting_int(&settings, "iter_max")?;
    let iter_tol = setting_real(&settings, "iter_tol")?;
    let ctx = Context { panel: &panel, labels: &labels, times: &times };
    let res = banerjee::estimate_breaks(&ctx, &opts, m as usize, h as usize, iter_max, iter_tol, n_threads as usize)
        .map_err(Error::Other)?;
    let nu = panel.n_units();
    let mut pos = vec![0i32; nu * m as usize];
    for (i, p) in res.pos.iter().enumerate() {
        for (j, v) in p.iter().enumerate() {
            pos[j * nu + i] = *v as i32;
        }
    }
    Ok(list!(pos = pos, ssr = res.ssr, converged = res.converged, r = res.r as i32))
}

/// t-bar and the number of factors for every admissible common break vector.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn banerjee_profile_rs(y: &[f64], x: &[f64], k: i32, offsets: &[i32], labels: Vec<String>, times: Vec<String>, m: i32, h: i32, settings: List, n_threads: i32) -> Result<List> {
    let panel = build_panel(y, x, k, offsets)?;
    let opts = parse_banerjee(&settings, panel.k)?;
    if !(1..=2).contains(&m) || h < 1 || n_threads < 1 {
        return err("`m` must be 1 or 2, `h` and `n_threads` positive");
    }
    let ctx = Context { panel: &panel, labels: &labels, times: &times };
    let rows = banerjee::profile(&ctx, &opts, m as usize, h as usize, n_threads as usize).map_err(Error::Other)?;
    let nc = rows.len();
    let mut pos = vec![0i32; nc * m as usize];
    for (c, row) in rows.iter().enumerate() {
        for (j, v) in row.pos.iter().enumerate() {
            pos[j * nc + c] = *v as i32;
        }
    }
    Ok(list!(
        pos = pos,
        tbar = rows.iter().map(|r| r.tbar).collect::<Vec<f64>>(),
        r = rows.iter().map(|r| r.r as i32).collect::<Vec<i32>>()
    ))
}

/// Simulated null distribution: `kind` is "idio", "factor_df" or "mq",
/// `det` "c", "tau" or "gamma"; returns a `reps x ncol` matrix as a vector.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn banerjee_null_rs(kind: &str, det: &str, lambda: &[f64], steps: i32, reps: i32, q: i32, seed: f64, n_threads: i32) -> Result<Vec<f64>> {
    let kind = match kind {
        "idio" => Kind::Idio,
        "factor_df" => Kind::FactorDf,
        "mq" if q >= 1 => Kind::Mq(q as usize),
        _ => return err("invalid `kind` or `q`"),
    };
    let det = match det {
        "c" => Det::C,
        "tau" => Det::Tau,
        "gamma" => Det::Gamma,
        _ => return err("`det` must be \"c\", \"tau\" or \"gamma\""),
    };
    if steps < 10 || reps < 1 || n_threads < 1 || !(seed >= 0.0) {
        return err("invalid simulation size");
    }
    nulldist::simulate(kind, det, lambda, steps as usize, reps as usize, seed as u64, n_threads as usize).map_err(Error::Other)
}

/// Kernel for tests: symmetric eigendecomposition (decreasing).
/// @noRd
#[extendr]
fn eigen_sym_rs(a: &[f64], n: i32) -> Result<List> {
    let eg = eigen::sym_eigen(a, n.max(0) as usize).map_err(Error::Other)?;
    Ok(list!(values = eg.values, vectors = eg.vectors))
}

/// Kernel for tests: generalised symmetric eigenproblem `A v = mu B v`.
/// @noRd
#[extendr]
fn eigen_gen_rs(a: &[f64], b: &[f64], n: i32) -> Result<List> {
    let eg = eigen::gen_sym_eigen(a, b, n.max(0) as usize).map_err(Error::Other)?;
    Ok(list!(values = eg.values, vectors = eg.vectors))
}

/// Kernel for tests: principal components of an `n x N` matrix.
/// @noRd
#[extendr]
fn pca_rs(y: &[f64], n: i32, criterion: &str, r_max: i32, r: i32) -> Result<List> {
    let n = n.max(1) as usize;
    let nn = y.len() / n;
    let crit = CRITERIA.iter().find(|(c, _)| *c == criterion).map(|(_, c)| *c).ok_or_else(|| Error::Other("unknown criterion".into()))?;
    let ge = factors::gram_eigen(y, n, nn).map_err(Error::Other)?;
    let crits = factors::criteria(&ge.values, n, nn, r_max.max(0) as usize, crit);
    let r = if r < 0 { factors::argmin(&crits) } else { r as usize };
    let ext = factors::extract(y, &ge, r).map_err(Error::Other)?;
    Ok(list!(values = ge.values, criteria = crits, r = r as i32, f = ext.f, loadings = ext.loadings, resid = ext.resid))
}

/// Kernel for tests: ADF regression without deterministic terms.
/// @noRd
#[extendr]
fn adf_rs(e: &[f64], settings: List) -> Result<List> {
    let spec = parse_adf(&settings)?;
    let mut ws = adf::AdfWorkspace::default();
    let f = adf::adf(e, &spec, &mut ws).map_err(Error::Other)?;
    Ok(list!(t = f.t, lags = f.k as i32, n = f.n as i32, resid = f.resid, start = f.start as i32))
}

/// Kernel for tests: break search of one unit on a target (`NULL`: `Dy`).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn break_search_rs(y: &[f64], x: &[f64], k: i32, model: i32, m: i32, h: i32, dols_p: i32, dols_q: i32, target: Nullable<Vec<f64>>) -> Result<List> {
    let panel = build_panel(y, x, k, &[0, y.len() as i32])?;
    let spec = Spec { model: model as u8, k: panel.k, dols_p: dols_p.max(0) as usize, dols_q: dols_q.max(0) as usize };
    let u = panel.unit(0);
    let (r0, n) = spec.rows(u.t_len).ok_or_else(|| Error::Other("too few periods".into()))?;
    let target = match target {
        Nullable::NotNull(t) => t,
        Nullable::Null => {
            let mut dy = Vec::new();
            breaks::fill_dy(&u, r0, n, &mut dy);
            dy
        }
    };
    if target.len() != n {
        return err("`target` has the wrong length");
    }
    let cands = breaks::candidates(u.t_len, m.max(0) as usize, h.max(1) as usize, r0, n);
    let (pos, ssr) = breaks::search(&u, &spec, &cands, r0, n, &target).ok_or_else(|| Error::Other("no admissible break dates".into()))?;
    Ok(list!(pos = pos.iter().map(|&p| p as i32).collect::<Vec<i32>>(), ssr = ssr, n_candidates = cands.len() as i32))
}

/// Kernel for tests: factor detrending, factor ADF and MQ statistics.
/// @noRd
#[extendr]
fn factor_tests_rs(f: &[f64], n: i32, det: &[f64], settings: List) -> Result<List> {
    let n = n.max(1) as usize;
    let r = f.len() / n;
    let spec = parse_adf(&settings)?;
    let t = mq::factor_tests(f, n, r, det, &spec, setting_int(&settings, "mq_bandwidth")?, setting_int(&settings, "mq_max_var")?)
        .map_err(Error::Other)?;
    Ok(list!(
        fd = t.fd,
        adf_t = t.adf.as_ref().map_or(f64::NAN, |a| a.t),
        adf_lags = t.adf.as_ref().map_or(-1, |a| a.k as i32),
        mq_c = t.mq.iter().map(|m| m.mq_c).collect::<Vec<f64>>(),
        mq_f = t.mq.iter().map(|m| m.mq_f).collect::<Vec<f64>>(),
        var_order = t.mq.iter().map(|m| m.p as i32).collect::<Vec<i32>>()
    ))
}

/// Kernel for tests: raw xoshiro256++ outputs (as doubles of the top 53
/// bits) and normal draws of one stream.
/// @noRd
#[extendr]
fn rng_rs(seed: f64, stream: f64, n: i32) -> List {
    let mut g = Rng::stream(seed as u64, stream as u64);
    let u: Vec<f64> = (0..n.max(0)).map(|_| g.uniform()).collect();
    let z: Vec<f64> = (0..n.max(0)).map(|_| g.normal()).collect();
    list!(uniform = u, normal = z)
}

extendr_module! {
    mod panelcointr;
    fn westerlund_rs;
    fn westerlund_boot_prepare_rs;
    fn westerlund_boot_info_rs;
    fn westerlund_boot_rs;
    fn westerlund_boot_state_rs;
    fn westerlund_boot_panel_rs;
    fn ols_rs;
    fn lrv_rs;
    fn ecm_ic_rs;
    fn banerjee_rs;
    fn banerjee_breaks_rs;
    fn banerjee_profile_rs;
    fn banerjee_null_rs;
    fn eigen_sym_rs;
    fn eigen_gen_rs;
    fn pca_rs;
    fn adf_rs;
    fn break_search_rs;
    fn factor_tests_rs;
    fn rng_rs;
}
