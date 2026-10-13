//! Design matrices of the conditional error-correction model
//!
//! `Δy_t = δ'd_t + α y_{t-1} + λ'x_{t-1} + Σ_{j=1..p} α_j Δy_{t-j}
//!         + Σ_{j=-q..p} γ_j'Δx_{t-j} + e_t`
//!
//! for one unit, with time `t` 0-based and the estimation sample
//! `p + 1 <= t <= T - 1 - q` (the paper's `p + 2 <= t <= T - q`).

use super::panel::Unit;

/// One regressor of the ECM.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Col {
    Const,
    /// Linear trend `t = 1..T` within the unit.
    Trend,
    /// `x_{j,t-1}`.
    XLag(usize),
    /// `y_{t-1}`.
    YLag,
    /// `Δy_{t-j}`, `j >= 1`.
    DyLag(usize),
    /// `Δx_{j,t-l}`: a lag for `l > 0`, a lead for `l < 0`.
    Dx(usize, isize),
}

impl Col {
    #[inline]
    pub fn value(self, u: &Unit, t: usize) -> f64 {
        match self {
            Col::Const => 1.0,
            Col::Trend => (t + 1) as f64,
            Col::XLag(j) => u.x(j, t - 1),
            Col::YLag => u.y(t - 1),
            Col::DyLag(j) => u.dy(t - j),
            Col::Dx(j, l) => u.dx(j, (t as isize - l) as usize),
        }
    }

    /// Label used in error messages; regressors are numbered from 1.
    pub fn label(self) -> String {
        let shift = |l: isize| match l {
            0 => "t".to_string(),
            l if l > 0 => format!("t-{l}"),
            l => format!("t+{}", -l),
        };
        match self {
            Col::Const => "constant".into(),
            Col::Trend => "trend".into(),
            Col::XLag(j) => format!("x{}(t-1)", j + 1),
            Col::YLag => "y(t-1)".into(),
            Col::DyLag(j) => format!("dy(t-{j})"),
            Col::Dx(j, l) => format!("dx{}({})", j + 1, shift(l)),
        }
    }
}

/// Number of ECM regressors with `dd` deterministic terms and `k` regressors.
pub fn n_cols(dd: usize, k: usize, p: usize, q: usize) -> usize {
    dd + k + p + k * (p + q + 1) + 1
}

fn push_deterministics(dd: usize, out: &mut Vec<Col>) {
    if dd >= 1 {
        out.push(Col::Const);
    }
    if dd >= 2 {
        out.push(Col::Trend);
    }
}

/// Regressors of the ECM for `(p, q)` in the fitting order
/// `[d_t, x_{t-1}, Δy_{t-1..t-p}, Δx_{t+q..t-p}, y_{t-1}]`, with `y_{t-1}` last
/// so that the last diagonal element of the QR gives the projection on it.
pub fn fit_columns(dd: usize, k: usize, p: usize, q: usize, out: &mut Vec<Col>) {
    out.clear();
    push_deterministics(dd, out);
    out.extend((0..k).map(Col::XLag));
    out.extend((1..=p).map(Col::DyLag));
    for l in -(q as isize)..=(p as isize) {
        out.extend((0..k).map(|j| Col::Dx(j, l)));
    }
    out.push(Col::YLag);
}

/// Position of `col` in `fit_columns(dd, k, p_max, q_max)`.
pub fn canonical_index(col: Col, dd: usize, k: usize, p_max: usize, q_max: usize) -> usize {
    match col {
        Col::Const => 0,
        Col::Trend => 1,
        Col::XLag(j) => dd + j,
        Col::DyLag(j) => dd + k + j - 1,
        Col::Dx(j, l) => dd + k + p_max + (l + q_max as isize) as usize * k + j,
        Col::YLag => dd + k + p_max + k * (p_max + q_max + 1),
    }
}

/// Column order in which every model `(p, q_fixed)`, `p <= p_max`, is a
/// leading block: `[d, x_{t-1}, y_{t-1}, Δx_t, Δx_{t+1..t+q}, (Δy_{t-1}, Δx_{t-1}), …]`.
/// Model `(p, q)` uses the first `nested_len(dd, k, p, q)` columns.
pub fn nested_by_lags(dd: usize, k: usize, p_max: usize, q: usize, out: &mut Vec<Col>) {
    out.clear();
    push_deterministics(dd, out);
    out.extend((0..k).map(Col::XLag));
    out.push(Col::YLag);
    for l in 0..=q as isize {
        out.extend((0..k).map(|j| Col::Dx(j, -l)));
    }
    for lag in 1..=p_max {
        out.push(Col::DyLag(lag));
        out.extend((0..k).map(|j| Col::Dx(j, lag as isize)));
    }
}

/// Column order in which every model `(p_fixed, q)`, `q <= q_max`, is a
/// leading block: `[d, x_{t-1}, y_{t-1}, Δx_t, (Δy_{t-1}, Δx_{t-1}), …, Δx_{t+1}, …]`.
pub fn nested_by_leads(dd: usize, k: usize, p: usize, q_max: usize, out: &mut Vec<Col>) {
    out.clear();
    push_deterministics(dd, out);
    out.extend((0..k).map(Col::XLag));
    out.push(Col::YLag);
    out.extend((0..k).map(|j| Col::Dx(j, 0)));
    for lag in 1..=p {
        out.push(Col::DyLag(lag));
        out.extend((0..k).map(|j| Col::Dx(j, lag as isize)));
    }
    for l in 1..=q_max as isize {
        out.extend((0..k).map(|j| Col::Dx(j, -l)));
    }
}

/// Fill the column-major `n x cols.len()` design for times `t0..t0 + n`.
pub fn fill(u: &Unit, cols: &[Col], t0: usize, n: usize, z: &mut Vec<f64>) {
    z.clear();
    z.resize(n * cols.len(), 0.0);
    for (col, dst) in cols.iter().zip(z.chunks_exact_mut(n.max(1))) {
        match *col {
            Col::Const => dst.fill(1.0),
            Col::Trend => dst.iter_mut().enumerate().for_each(|(r, v)| *v = (t0 + r + 1) as f64),
            Col::XLag(j) => dst.iter_mut().enumerate().for_each(|(r, v)| *v = u.x(j, t0 + r - 1)),
            Col::YLag => dst.iter_mut().enumerate().for_each(|(r, v)| *v = u.y(t0 + r - 1)),
            Col::DyLag(j) => dst.iter_mut().enumerate().for_each(|(r, v)| *v = u.dy(t0 + r - j)),
            Col::Dx(j, l) => dst
                .iter_mut()
                .enumerate()
                .for_each(|(r, v)| *v = u.dx(j, ((t0 + r) as isize - l) as usize)),
        }
    }
}

/// Fill `Δy_t` for times `t0..t0 + n`.
pub fn fill_dy(u: &Unit, t0: usize, n: usize, dy: &mut Vec<f64>) {
    dy.clear();
    dy.extend((t0..t0 + n).map(|t| u.dy(t)));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layouts_agree_on_columns() {
        let (dd, k, p_max, q_max) = (2, 2, 3, 2);
        let mut full = Vec::new();
        fit_columns(dd, k, p_max, q_max, &mut full);
        assert_eq!(full.len(), n_cols(dd, k, p_max, q_max));
        for (i, c) in full.iter().enumerate() {
            assert_eq!(canonical_index(*c, dd, k, p_max, q_max), i);
        }
        let mut nested = Vec::new();
        for q in 0..=q_max {
            nested_by_lags(dd, k, p_max, q, &mut nested);
            for p in 0..=p_max {
                let len = n_cols(dd, k, p, q);
                let mut a: Vec<usize> = nested[..len]
                    .iter()
                    .map(|c| canonical_index(*c, dd, k, p_max, q_max))
                    .collect();
                let mut fit = Vec::new();
                fit_columns(dd, k, p, q, &mut fit);
                let mut b: Vec<usize> = fit.iter().map(|c| canonical_index(*c, dd, k, p_max, q_max)).collect();
                a.sort_unstable();
                b.sort_unstable();
                assert_eq!(a, b, "lags nested, p = {p}, q = {q}");
            }
        }
        for p in 0..=p_max {
            nested_by_leads(dd, k, p, q_max, &mut nested);
            for q in 0..=q_max {
                let len = n_cols(dd, k, p, q);
                let mut a: Vec<usize> = nested[..len]
                    .iter()
                    .map(|c| canonical_index(*c, dd, k, p_max, q_max))
                    .collect();
                let mut fit = Vec::new();
                fit_columns(dd, k, p, q, &mut fit);
                let mut b: Vec<usize> = fit.iter().map(|c| canonical_index(*c, dd, k, p_max, q_max)).collect();
                a.sort_unstable();
                b.sort_unstable();
                assert_eq!(a, b, "leads nested, p = {p}, q = {q}");
            }
        }
    }
}
