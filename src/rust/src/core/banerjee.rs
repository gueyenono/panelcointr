//! Panel cointegration tests of Banerjee and Carrion-i-Silvestre (2015):
//! per-unit projection in first differences, principal components, ADF
//! tests of the cumulated idiosyncratic components and tests of the common
//! factors, for known break dates; the iterative estimation of
//! heterogeneous break dates (Section 4.2.1); and the profile over common
//! break dates behind the infimum statistic (Section 4.2.2).

use super::adf::{self, AdfFit, AdfSpec, AdfWorkspace};
use super::breaks::{self, Spec};
use super::factors::{self, Criterion, Extracted};
use super::linalg::Qr;
use super::mq::{self, FactorTests};
use super::panel::Panel;
use super::parallel::par_map;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Factors {
    Fixed(usize),
    Select(Criterion),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Options {
    pub spec: Spec,
    pub factors: Factors,
    pub r_max: usize,
    pub adf: AdfSpec,
    pub pca_scale: bool,
    pub mq_bandwidth: usize,
    pub mq_max_var: usize,
}

/// The panel with the labels used in error messages.
pub struct Context<'a> {
    pub panel: &'a Panel,
    /// Unit identifiers.
    pub labels: &'a [String],
    /// Time values of the common calendar.
    pub times: &'a [String],
}

impl Context<'_> {
    fn unit(&self, i: usize) -> String {
        self.labels.get(i).cloned().unwrap_or_else(|| format!("#{}", i + 1))
    }

    fn dates(&self, pos: &[usize]) -> String {
        if pos.is_empty() {
            return "no breaks".into();
        }
        let d: Vec<String> = pos.iter().map(|&p| self.times.get(p).cloned().unwrap_or_else(|| format!("period {}", p + 1))).collect();
        format!("breaks after {}", d.join(", "))
    }

    /// Common number of periods and the first-difference sample rows.
    pub fn sample(&self, spec: &Spec) -> Result<(usize, usize, usize), String> {
        let p = self.panel;
        let t_len = p.t_len(0);
        if (1..p.n_units()).any(|i| p.t_len(i) != t_len) {
            return Err("the panel must be balanced".into());
        }
        let (r0, n) = spec.rows(t_len).ok_or_else(|| "too few periods for the first-difference regression".to_string())?;
        Ok((t_len, r0, n))
    }
}

/// Residuals of the units' first-difference regressions, `n x N`.
pub struct Projected {
    pub r0: usize,
    pub n: usize,
    pub y: Vec<f64>,
}

pub fn project(ctx: &Context, spec: &Spec, pos: &[Vec<usize>]) -> Result<Projected, String> {
    let (_, r0, n) = ctx.sample(spec)?;
    let nn = ctx.panel.n_units();
    if pos.len() != nn {
        return Err("one vector of break positions per unit is needed".into());
    }
    let mut y = Vec::with_capacity(n * nn);
    let (mut qr, mut z, mut dy) = (Qr::default(), Vec::new(), Vec::new());
    for (i, pi) in pos.iter().enumerate().take(nn) {
        let u = ctx.panel.unit(i);
        breaks::fill_dy(&u, r0, n, &mut dy);
        let e = breaks::residuals(&u, spec, pi, r0, n, &dy, &mut qr, &mut z).map_err(|_| {
            format!(
                "unit {}: the first-difference regression is rank deficient ({}); \
                 a regime may be too short or a regressor constant within a regime",
                ctx.unit(i),
                ctx.dates(pi)
            )
        })?;
        y.extend(e);
    }
    Ok(Projected { r0, n, y })
}

#[derive(Clone, Debug, Default)]
pub struct FactorStep {
    pub values: Vec<f64>,
    pub criteria: Vec<f64>,
    pub r: usize,
    pub ext: Extracted,
}

/// Principal components of `y` (`n x nn`) with the number of factors fixed
/// or selected; with `pca_scale`, on the standardised columns, with the
/// residuals and loadings scaled back.
pub fn factor_step(y: &[f64], n: usize, nn: usize, opts: &Options) -> Result<FactorStep, String> {
    let sd: Option<Vec<f64>> = if opts.pca_scale {
        let sd: Vec<f64> = (0..nn)
            .map(|i| {
                let c = &y[i * n..(i + 1) * n];
                let m = c.iter().sum::<f64>() / n as f64;
                (c.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (n as f64 - 1.0)).sqrt()
            })
            .collect();
        if sd.iter().any(|s| !(*s > 0.0)) {
            return Err("a unit has no variation left after the first-difference regression".into());
        }
        Some(sd)
    } else {
        None
    };
    let scaled: Vec<f64>;
    let ys = match &sd {
        Some(sd) => {
            scaled = (0..nn).flat_map(|i| y[i * n..(i + 1) * n].iter().map(move |v| v / sd[i])).collect();
            &scaled[..]
        }
        None => y,
    };
    let ge = factors::gram_eigen(ys, n, nn)?;
    let crit = match opts.factors {
        Factors::Select(c) => c,
        Factors::Fixed(_) => Criterion::Bic3,
    };
    let criteria = factors::criteria(&ge.values, n, nn, opts.r_max, crit);
    let r = match opts.factors {
        Factors::Select(_) => factors::argmin(&criteria),
        Factors::Fixed(r) => r,
    };
    let mut ext = factors::extract(ys, &ge, r)?;
    if let Some(sd) = &sd {
        for (i, s) in sd.iter().enumerate() {
            ext.resid[i * n..(i + 1) * n].iter_mut().for_each(|v| *v *= s);
            for j in 0..r {
                ext.loadings[j * nn + i] *= s;
            }
        }
    }
    Ok(FactorStep { values: ge.values, criteria, r, ext })
}

/// ADF tests of the cumulated idiosyncratic components (`z`: `n x nn`).
pub fn idiosyncratic(ctx: &Context, z: &[f64], n: usize, adf_spec: &AdfSpec) -> Result<Vec<AdfFit>, String> {
    let nn = z.len() / n;
    let mut ws = AdfWorkspace::default();
    let mut e = vec![0.0; n];
    (0..nn)
        .map(|i| {
            let mut acc = 0.0;
            for (ev, zv) in e.iter_mut().zip(&z[i * n..(i + 1) * n]) {
                acc += zv;
                *ev = acc;
            }
            adf::adf(&e, adf_spec, &mut ws).map_err(|m| format!("unit {}: {m}", ctx.unit(i)))
        })
        .collect()
}

#[derive(Clone, Debug, Default)]
pub struct KnownResult {
    pub r0: usize,
    pub n: usize,
    pub adf: Vec<AdfFit>,
    pub tbar: f64,
    pub factor: FactorStep,
    pub tests: Option<FactorTests>,
}

/// The whole procedure for given break positions (one vector per unit).
/// Models 3 and 6 detrend the factors on the breaks of the first unit,
/// which are common.
pub fn run_known(ctx: &Context, opts: &Options, pos: &[Vec<usize>], factor_tests: bool) -> Result<KnownResult, String> {
    let nn = ctx.panel.n_units();
    let pr = project(ctx, &opts.spec, pos)?;
    let fs = factor_step(&pr.y, pr.n, nn, opts)?;
    let adf = idiosyncratic(ctx, &fs.ext.resid, pr.n, &opts.adf)?;
    let tbar = adf.iter().map(|f| f.t).sum::<f64>() / nn as f64;
    let tests = if factor_tests && fs.r > 0 {
        let spec = &opts.spec;
        let det = match spec.model {
            1 | 4 => mq::deterministics(pr.n, false, &[]),
            2 | 5 => mq::deterministics(pr.n, true, &[]),
            _ => {
                let b: Vec<usize> = pos[0].iter().map(|&p| p - pr.r0).collect();
                mq::deterministics(pr.n, true, &b)
            }
        };
        Some(mq::factor_tests(&fs.ext.f, pr.n, fs.r, &det, &opts.adf, opts.mq_bandwidth, opts.mq_max_var)?)
    } else {
        None
    };
    Ok(KnownResult { r0: pr.r0, n: pr.n, adf, tbar, factor: fs, tests })
}

/// Break positions of a unit and its residuals.
type UnitBreaks = (Vec<usize>, Vec<f64>);

#[derive(Clone, Debug, Default)]
pub struct HeteroResult {
    pub pos: Vec<Vec<usize>>,
    /// Sum of squared residuals after the principal components, per iteration.
    pub ssr: Vec<f64>,
    pub converged: bool,
    pub r: usize,
}

/// Iterative estimation of unit-specific break dates (Section 4.2.1):
/// break dates and slopes unit by unit, then principal components of the
/// residuals with `r` fixed after the first round, then break dates again
/// on the data net of the common component, until the relative decrease of
/// the sum of squared residuals falls below `iter_tol`.
pub fn estimate_breaks(ctx: &Context, opts: &Options, m: usize, h: usize, iter_max: usize, iter_tol: f64, n_threads: usize) -> Result<HeteroResult, String> {
    let spec = &opts.spec;
    let (t_len, r0, n) = ctx.sample(spec)?;
    let nn = ctx.panel.n_units();
    let cands = breaks::candidates(t_len, m, h, r0, n);
    if cands.is_empty() {
        return Err(format!("no admissible dates for {m} break(s) with a minimum regime length of {h} periods"));
    }
    // Search every unit given the common component `c` (n x nn).
    let search_all = |c: Option<&[f64]>| -> Result<(Vec<Vec<usize>>, Vec<f64>), String> {
        let res: Vec<Result<UnitBreaks, String>> = par_map(nn, n_threads, |i| {
            let u = ctx.panel.unit(i);
            let mut dy = Vec::new();
            breaks::fill_dy(&u, r0, n, &mut dy);
            let ci = c.map(|c| &c[i * n..(i + 1) * n]);
            let target: Vec<f64> = match ci {
                Some(ci) => dy.iter().zip(ci).map(|(a, b)| a - b).collect(),
                None => dy,
            };
            let (pos, _) = breaks::search(&u, spec, &cands, r0, n, &target)
                .ok_or_else(|| format!("unit {}: no admissible break dates give a full-rank regression", ctx.unit(i)))?;
            let (mut qr, mut z) = (Qr::default(), Vec::new());
            let mut e = breaks::residuals(&u, spec, &pos, r0, n, &target, &mut qr, &mut z)
                .map_err(|_| format!("unit {}: the first-difference regression is rank deficient ({})", ctx.unit(i), ctx.dates(&pos)))?;
            if let Some(ci) = ci {
                e.iter_mut().zip(ci).for_each(|(a, b)| *a += b);
            }
            Ok((pos, e))
        });
        let mut pos = Vec::with_capacity(nn);
        let mut u = Vec::with_capacity(n * nn);
        for r in res {
            let (p, e) = r?;
            pos.push(p);
            u.extend(e);
        }
        Ok((pos, u))
    };

    let (mut pos, mut u) = search_all(None)?;
    let ge = factors::gram_eigen(&u, n, nn)?;
    let r = match opts.factors {
        Factors::Fixed(r) => r,
        Factors::Select(c) => factors::argmin(&factors::criteria(&ge.values, n, nn, opts.r_max, c)),
    };
    let mut ssr = Vec::new();
    let mut converged = false;
    for it in 0..iter_max.max(1) {
        let ge = factors::gram_eigen(&u, n, nn)?;
        let s: f64 = ge.values[r..].iter().sum();
        ssr.push(s);
        if it > 0 {
            let prev = ssr[it - 1];
            if !(prev > 0.0) || (prev - s) / prev < iter_tol {
                converged = true;
                break;
            }
        }
        if r == 0 {
            converged = true;
            break;
        }
        let ext = factors::extract(&u, &ge, r)?;
        let c: Vec<f64> = u.iter().zip(&ext.resid).map(|(a, b)| a - b).collect();
        let next = search_all(Some(&c))?;
        pos = next.0;
        u = next.1;
    }
    Ok(HeteroResult { pos, ssr, converged, r })
}

#[derive(Clone, Debug, Default)]
pub struct ProfileRow {
    pub pos: Vec<usize>,
    /// NaN if the procedure fails at these dates.
    pub tbar: f64,
    pub r: usize,
}

/// `t-bar` and the number of factors for every admissible common break
/// vector (Section 4.2.2).
pub fn profile(ctx: &Context, opts: &Options, m: usize, h: usize, n_threads: usize) -> Result<Vec<ProfileRow>, String> {
    let (t_len, r0, n) = ctx.sample(&opts.spec)?;
    let nn = ctx.panel.n_units();
    let cands = breaks::candidates(t_len, m, h, r0, n);
    if cands.is_empty() {
        return Err(format!("no admissible dates for {m} break(s) with a minimum regime length of {h} periods"));
    }
    Ok(par_map(cands.len(), n_threads, |c| {
        let pos = vec![cands[c].clone(); nn];
        match run_known(ctx, opts, &pos, false) {
            Ok(k) => ProfileRow { pos: cands[c].clone(), tbar: k.tbar, r: k.factor.r },
            Err(_) => ProfileRow { pos: cands[c].clone(), tbar: f64::NAN, r: 0 },
        }
    }))
}
