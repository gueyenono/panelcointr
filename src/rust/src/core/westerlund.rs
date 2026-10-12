//! Westerlund (2007) error-correction panel cointegration statistics
//! Gt, Ga, Pt and Pa, including the variants needed to reproduce `xtwest`.

use super::design;
use super::ecm::{self, EcmOptions, GridMethod, LagSpec, UnitFit, Workspace};
use super::panel::Panel;

/// Normalising sample size `T̃_i` in Ga and Pa.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TNorm {
    /// Effective number of observations `n_i` (paper's recommendation).
    Effective,
    /// Residual degrees of freedom `n_i - k_i` (xtwest).
    Df,
    /// Number of time periods `T_i`.
    T,
}

/// Lag and lead orders of the panel statistics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PanelLags {
    /// Each unit's selected `(p_i, q_i)` (paper §3.2).
    Individual,
    /// `floor(mean p_i)`, `floor(mean q_i)` for every unit (xtwest).
    Mean,
}

/// Variance estimator in `SE(α̂_i)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeAlpha {
    /// `RSS / (n - k)`.
    Df,
    /// `RSS / n` (xtwest's `westerlund`).
    N,
}

/// Degrees of freedom of `σ̂_i` in the panel statistics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SigmaDf {
    /// Each unit's `n_i - k_i`.
    Own,
    /// The panel normaliser `T̄` for every unit (xtwest).
    Tbar,
}

/// Averaging of the group-mean statistics' unit terms in Ga.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Weights {
    /// Simple mean over units.
    Units,
    /// Weighted by `T_i`, xtwest's mean over rows.
    Rows,
}

/// How `T̄` in Pa is formed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TbarMode {
    /// Mean over units of `T̃_i` from the panel regressions.
    Units,
    /// xtwest: `T̃` evaluated at the mean lags and leads and at the row
    /// mean of `T_i`, `Σ T_i² / Σ T_i`.
    Rows,
}

#[derive(Clone, Debug)]
pub struct Settings {
    pub ecm: EcmOptions,
    pub t_norm: TNorm,
    pub panel_lags: PanelLags,
    pub se_alpha: SeAlpha,
    pub sigma: SigmaDf,
    pub ga_weights: Weights,
    pub tbar: TbarMode,
}

/// Everything computed for one panel.
#[derive(Clone, Debug)]
pub struct Output {
    /// Gt, Ga, Pt, Pa.
    pub stats: [f64; 4],
    /// Group regressions at the selected orders.
    pub group: Vec<UnitFit>,
    /// Panel regressions (equal to `group` with individual panel lags).
    pub panel: Vec<UnitFit>,
    /// `T̃_i` of the group regressions.
    pub t_tilde: Vec<f64>,
    /// `SE(α̂_i)` of the group regressions.
    pub se: Vec<f64>,
    /// `σ̂_i` of the panel statistics.
    pub sigma_panel: Vec<f64>,
    pub p_bar: usize,
    pub q_bar: usize,
    /// Pooled `α̂` and its standard error.
    pub alpha: f64,
    pub se_pooled: f64,
    pub t_bar: f64,
}

fn t_tilde(mode: TNorm, f: &UnitFit, t_len: usize) -> f64 {
    match mode {
        TNorm::Effective => f.n as f64,
        TNorm::Df => (f.n - f.k) as f64,
        TNorm::T => t_len as f64,
    }
}

/// Compute the four statistics. `specs` and `bw` (bandwidths) are per unit,
/// `labels` names the units in error messages. If `coef` is given, the
/// group regressions' coefficients are written to it, unit by unit, in the
/// layout `design::fit_columns(dd, k, p_cap, q_cap)` (NaN where absent).
pub fn compute(
    panel: &Panel,
    specs: &[LagSpec],
    bw: &[usize],
    s: &Settings,
    ws: &mut Workspace,
    labels: &[String],
    mut coef: Option<(&mut Vec<f64>, usize, usize)>,
) -> Result<Output, String> {
    let nu = panel.n_units();
    let o = &s.ecm;
    let k = panel.k;
    let label = |i: usize| labels.get(i).cloned().unwrap_or_else(|| format!("#{}", i + 1));

    if let Some((c, p_cap, q_cap)) = coef.as_mut() {
        let width = design::n_cols(o.dd, k, *p_cap, *q_cap);
        c.clear();
        c.resize(nu * width, f64::NAN);
    }

    let mut group = Vec::with_capacity(nu);
    for i in 0..nu {
        let u = panel.unit(i);
        let (p, q) = ecm::select(&u, &specs[i], o, ws, GridMethod::Gram).map_err(|e| format!("unit {}: {e}", label(i)))?;
        let f = ecm::fit(&u, p, q, bw[i], o, ws).map_err(|e| format!("unit {}: {e}", label(i)))?;
        if let Some((c, p_cap, q_cap)) = coef.as_mut() {
            let width = design::n_cols(o.dd, k, *p_cap, *q_cap);
            for (col, b) in ws.cols.iter().zip(&ws.b) {
                c[i * width + design::canonical_index(*col, o.dd, k, *p_cap, *q_cap)] = *b;
            }
        }
        group.push(f);
    }

    // Group-mean statistics.
    let se: Vec<f64> = group.iter().map(|f| f.se_alpha(s.se_alpha == SeAlpha::Df)).collect();
    let t_tilde_g: Vec<f64> = (0..nu).map(|i| t_tilde(s.t_norm, &group[i], panel.t_len(i))).collect();
    let gt = group.iter().zip(&se).map(|(f, e)| f.alpha / e).sum::<f64>() / nu as f64;
    let (mut ga_num, mut ga_den) = (0.0, 0.0);
    for i in 0..nu {
        let w = match s.ga_weights {
            Weights::Units => 1.0,
            Weights::Rows => panel.t_len(i) as f64,
        };
        ga_num += w * t_tilde_g[i] * group[i].alpha / group[i].alpha1;
        ga_den += w;
    }
    let ga = ga_num / ga_den;

    // Panel statistics.
    let p_bar = group.iter().map(|f| f.p).sum::<usize>() / nu;
    let q_bar = group.iter().map(|f| f.q).sum::<usize>() / nu;
    let pan: Vec<UnitFit> = match s.panel_lags {
        PanelLags::Individual => group.clone(),
        PanelLags::Mean => {
            let mut v = Vec::with_capacity(nu);
            for (i, &m) in bw.iter().enumerate().take(nu) {
                let u = panel.unit(i);
                v.push(ecm::fit(&u, p_bar, q_bar, m, o, ws).map_err(|e| format!("unit {}: {e}", label(i)))?);
            }
            v
        }
    };
    let t_bar = match s.tbar {
        TbarMode::Units => (0..nu).map(|i| t_tilde(s.t_norm, &pan[i], panel.t_len(i))).sum::<f64>() / nu as f64,
        TbarMode::Rows => {
            let (s1, s2) = (0..nu).fold((0.0, 0.0), |(a, b), i| {
                let t = panel.t_len(i) as f64;
                (a + t, b + t * t)
            });
            let t_rows = s2 / s1;
            match s.t_norm {
                TNorm::T => t_rows,
                TNorm::Effective => t_rows - (p_bar + q_bar + 1) as f64,
                TNorm::Df => t_rows - (p_bar + q_bar + 1) as f64 - design::n_cols(o.dd, k, p_bar, q_bar) as f64,
            }
        }
    };
    let (mut num, mut den) = (0.0, 0.0);
    for f in &pan {
        num += f.r_kk * f.c_k / f.alpha1;
        den += f.r_kk * f.r_kk;
    }
    let alpha = num / den;
    let sigma_panel: Vec<f64> = pan
        .iter()
        .map(|f| match s.sigma {
            SigmaDf::Own => f.s(),
            SigmaDf::Tbar => (f.rss / t_bar).sqrt(),
        })
        .collect();
    let s2 = pan.iter().zip(&sigma_panel).map(|(f, sg)| (sg / f.alpha1).powi(2)).sum::<f64>() / nu as f64;
    let se_pooled = s2.sqrt() / den.sqrt();
    let pt = alpha / se_pooled;
    let pa = t_bar * alpha;

    Ok(Output {
        stats: [gt, ga, pt, pa],
        group,
        panel: pan,
        t_tilde: t_tilde_g,
        se,
        sigma_panel,
        p_bar,
        q_bar,
        alpha,
        se_pooled,
        t_bar,
    })
}
