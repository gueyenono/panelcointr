//! R bindings. The wrappers only convert and validate inputs; the
//! computations live in `core`, which has no extendr types.

// `!(x > 0.0)` is used on purpose: unlike `x <= 0.0`, it also rejects NaN.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use extendr_api::prelude::*;
use extendr_api::Result;

pub mod core;

use crate::core::bootstrap::{self, BootState};
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
}
