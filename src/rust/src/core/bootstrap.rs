//! Bootstrap of Westerlund (2007, §5), extended to leads: null-imposed
//! first stage, resampling of common time indices (drawn in R), and the
//! replicate loop, parallel over replicates.

use super::design::{self, Col};
use super::ecm::{LagSpec, Workspace};
use super::linalg::Qr;
use super::panel::Panel;
use super::westerlund::{self, Settings};

/// First-stage coefficients of one unit and its place on the calendar.
#[derive(Clone, Debug)]
pub struct BootUnit {
    /// Calendar position of the unit's first period (0 = first period of the panel).
    pub start: usize,
    pub t_len: usize,
    pub p: usize,
    pub q: usize,
    /// `α̂_j`, `j = 1..p`.
    pub a: Vec<f64>,
    /// `γ̂_{l,j}` at `(l + q) * k + j` for `l = -q..p`.
    pub g: Vec<f64>,
}

/// Everything a replicate needs, built once and kept in R as an external pointer.
#[derive(Clone, Debug)]
pub struct BootState {
    pub units: Vec<BootUnit>,
    pub k: usize,
    /// Centred `w_t = (ê_1t..ê_Nt, Δx_1t'..Δx_Nt')`, `n_periods x (N + N k)`, column-major.
    pub w: Vec<f64>,
    pub n_periods: usize,
    pub burn: usize,
    /// Length of the calendar range of the panel.
    pub t_span: usize,
    /// Indices drawn per replicate: `burn + t_span + max_i q_i`.
    pub len: usize,
    pub offsets: Vec<usize>,
    pub specs: Vec<LagSpec>,
    pub bw: Vec<usize>,
    pub settings: Settings,
    pub labels: Vec<String>,
}

/// First-stage coefficients `α̂_j` and `γ̂_{l,j}` and residuals of one unit.
type FirstStage = (Vec<f64>, Vec<f64>, Vec<f64>);

/// Null-imposed regression of `Δy_t` on `Δy_{t-1..t-p}` and `Δx_{t+q..t-p}`
/// (plus a constant in Case 3) on `S(p, q)`. Returns the coefficients and
/// the residuals.
fn first_stage(panel: &Panel, i: usize, p: usize, q: usize, constant: bool) -> Result<FirstStage, String> {
    let u = panel.unit(i);
    let k = panel.k;
    let mut cols = Vec::new();
    if constant {
        cols.push(Col::Const);
    }
    cols.extend((1..=p).map(Col::DyLag));
    for l in -(q as isize)..=(p as isize) {
        cols.extend((0..k).map(|j| Col::Dx(j, l)));
    }
    let m = cols.len();
    let t0 = p + 1;
    if u.t_len < p + q + 2 || u.t_len - p - q - 1 <= m {
        return Err(format!("too few observations for the first-stage regression with {p} lags and {q} leads"));
    }
    let n = u.t_len - p - q - 1;
    let (mut z, mut dy, mut b) = (Vec::new(), Vec::new(), Vec::new());
    design::fill(&u, &cols, t0, n, &mut z);
    design::fill_dy(&u, t0, n, &mut dy);
    let mut qr = Qr::default();
    qr.factor(&z, n, m, &dy).map_err(|e| {
        format!(
            "the first-stage bootstrap regression with {p} lags and {q} leads is rank deficient ({} is collinear with the preceding regressors)",
            cols[e.column.min(m - 1)].label()
        )
    })?;
    qr.coef(&mut b);
    let resid: Vec<f64> = (0..n)
        .map(|r| dy[r] - (0..m).map(|c| z[c * n + r] * b[c]).sum::<f64>())
        .collect();
    let off = usize::from(constant);
    let a = b[off..off + p].to_vec();
    let g = b[off + p..].to_vec();
    Ok((a, g, resid))
}

/// Build the bootstrap state. `start` holds 0-based calendar positions,
/// `sel` the `(p_i, q_i)` selected on the data. With `reselect`, the
/// replicates select their own orders over `specs`; otherwise they use `sel`.
#[allow(clippy::too_many_arguments)]
pub fn prepare(
    panel: &Panel,
    start: &[usize],
    sel: &[(usize, usize)],
    specs: &[LagSpec],
    bw: &[usize],
    settings: &Settings,
    reselect: bool,
    burn: usize,
    labels: &[String],
) -> Result<BootState, String> {
    let nu = panel.n_units();
    let k = panel.k;
    let constant = settings.ecm.dd == 2;
    let mut units = Vec::with_capacity(nu);
    let mut resid = Vec::with_capacity(nu);
    for i in 0..nu {
        let (p, q) = sel[i];
        let (a, g, e) = first_stage(panel, i, p, q, constant).map_err(|e| {
            format!("unit {}: {e}", labels.get(i).cloned().unwrap_or_else(|| format!("#{}", i + 1)))
        })?;
        units.push(BootUnit { start: start[i], t_len: panel.t_len(i), p, q, a, g });
        resid.push(e);
    }
    let t_span = units.iter().map(|b| b.start + b.t_len).max().unwrap_or(0);
    // Calendar periods where every unit has a first-stage residual.
    let lo = units.iter().map(|b| b.start + b.p + 1).max().unwrap_or(0);
    let hi = units.iter().map(|b| b.start + b.t_len - 1 - b.q).min().unwrap_or(0);
    let n_periods = if hi >= lo { hi - lo + 1 } else { 0 };
    let width = nu + nu * k;
    let mut w = vec![0.0; n_periods * width];
    for r in 0..n_periods {
        let cal = lo + r;
        for (i, b) in units.iter().enumerate() {
            let t = cal - b.start;
            w[i * n_periods + r] = resid[i][t - b.p - 1];
            let u = panel.unit(i);
            for j in 0..k {
                w[(nu + i * k + j) * n_periods + r] = u.dx(j, t);
            }
        }
    }
    if n_periods > 0 {
        for col in w.chunks_exact_mut(n_periods) {
            let m = col.iter().sum::<f64>() / n_periods as f64;
            col.iter_mut().for_each(|v| *v -= m);
        }
    }
    let q_need = units.iter().map(|b| b.q).max().unwrap_or(0);
    let specs = if reselect {
        specs.to_vec()
    } else {
        sel.iter().map(|&(p, q)| LagSpec::fixed(p, q)).collect()
    };
    Ok(BootState {
        units,
        k,
        w,
        n_periods,
        burn,
        t_span,
        len: burn + t_span + q_need,
        offsets: panel.offsets.clone(),
        specs,
        bw: bw.to_vec(),
        settings: settings.clone(),
        labels: labels.to_vec(),
    })
}

/// Per-thread buffers.
pub struct BootWorkspace {
    panel: Panel,
    s: Vec<usize>,
    dys: Vec<f64>,
    dxs: Vec<f64>,
    ecm: Workspace,
}

impl BootWorkspace {
    pub fn new(st: &BootState) -> Self {
        let rows = *st.offsets.last().unwrap_or(&0);
        BootWorkspace {
            panel: Panel { y: vec![0.0; rows], x: vec![0.0; rows * st.k], offsets: st.offsets.clone(), k: st.k },
            s: Vec::with_capacity(st.len),
            dys: vec![0.0; st.len],
            dxs: vec![0.0; st.len * st.k],
            ecm: Workspace::default(),
        }
    }

    /// The panel of the last generated replicate.
    pub fn panel(&self) -> &Panel {
        &self.panel
    }
}

/// Generate the bootstrap panel for one replicate from its `len` indices
/// (1-based, as drawn in R) into `ws.panel`.
pub fn generate(st: &BootState, idx: &[i32], ws: &mut BootWorkspace) -> Result<(), String> {
    let (l, k, np, nu) = (st.len, st.k, st.n_periods, st.units.len());
    if idx.len() != l {
        return Err(format!("expected {l} indices per replicate, got {}", idx.len()));
    }
    ws.s.clear();
    for &v in idx {
        if v < 1 || v as usize > np {
            return Err(format!("bootstrap index {v} outside 1..{np}"));
        }
        ws.s.push(v as usize - 1);
    }
    let rows = ws.panel.n_rows();
    for (i, b) in st.units.iter().enumerate() {
        // Periods needed: Δy* up to the unit's last calendar period, Δx* also
        // over the q leads after it.
        let need = st.burn + b.start + b.t_len;
        let need_x = (need + b.q).min(l);
        for j in 0..k {
            let src = &st.w[(nu + i * k + j) * np..(nu + i * k + j + 1) * np];
            let dst = &mut ws.dxs[j * l..j * l + need_x];
            for (d, &s) in dst.iter_mut().zip(&ws.s) {
                *d = src[s];
            }
        }
        let e = &st.w[i * np..(i + 1) * np];
        let (p, q) = (b.p as isize, b.q as isize);
        for tau in 0..need {
            let mut u = e[ws.s[tau]];
            for lag in -q..=p {
                let src = tau as isize - lag;
                if src < 0 {
                    continue;
                }
                let gl = &b.g[(lag + q) as usize * k..(lag + q + 1) as usize * k];
                for (j, g) in gl.iter().enumerate() {
                    u += g * ws.dxs[j * l + src as usize];
                }
            }
            let mut dy = u;
            for (j, a) in b.a.iter().enumerate() {
                if tau > j {
                    dy += a * ws.dys[tau - j - 1];
                }
            }
            ws.dys[tau] = dy;
        }
        // Levels: cumulative sums over the unit's own window, started at zero.
        let off = st.offsets[i];
        let base = st.burn + b.start;
        let mut acc = 0.0;
        for t in 0..b.t_len {
            acc += ws.dys[base + t];
            ws.panel.y[off + t] = acc;
        }
        for j in 0..k {
            let mut acc = 0.0;
            for t in 0..b.t_len {
                acc += ws.dxs[j * l + base + t];
                ws.panel.x[j * rows + off + t] = acc;
            }
        }
    }
    Ok(())
}

/// Statistics of one replicate.
pub fn replicate(st: &BootState, idx: &[i32], ws: &mut BootWorkspace) -> Result<[f64; 4], String> {
    generate(st, idx, ws)?;
    let out = westerlund::compute(&ws.panel, &st.specs, &st.bw, &st.settings, &mut ws.ecm, &st.labels, None)?;
    Ok(out.stats)
}

/// Run `n_rep` replicates whose indices are the consecutive blocks of `len`
/// values in `idx`, on `n_threads` threads. Returns the statistics row by
/// row (`n_rep x 4`, row-major). Each replicate depends only on its own
/// indices, so the result does not depend on the number of threads.
pub fn run(st: &BootState, idx: &[i32], n_rep: usize, n_threads: usize) -> Result<Vec<f64>, String> {
    if idx.len() != n_rep * st.len {
        return Err(format!("expected {} bootstrap indices, got {}", n_rep * st.len, idx.len()));
    }
    let mut out = vec![f64::NAN; n_rep * 4];
    if n_rep == 0 {
        return Ok(out);
    }
    let work = |first: usize, rows: &mut [f64]| -> Result<(), String> {
        let mut ws = BootWorkspace::new(st);
        for (r, row) in rows.chunks_exact_mut(4).enumerate() {
            let b = first + r;
            let stats = replicate(st, &idx[b * st.len..(b + 1) * st.len], &mut ws)
                .map_err(|e| format!("bootstrap replicate {}: {e}", b + 1))?;
            row.copy_from_slice(&stats);
        }
        Ok(())
    };
    let n_threads = n_threads.clamp(1, n_rep);
    if n_threads == 1 {
        work(0, &mut out)?;
        return Ok(out);
    }
    let per = (n_rep + n_threads - 1) / n_threads;
    let results: Vec<Result<(), String>> = std::thread::scope(|scope| {
        let handles: Vec<_> = out
            .chunks_mut(per * 4)
            .enumerate()
            .map(|(t, rows)| {
                let work = &work;
                scope.spawn(move || work(t * per, rows))
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|p| std::panic::resume_unwind(p)))
            .collect()
    });
    for r in results {
        r?;
    }
    Ok(out)
}
