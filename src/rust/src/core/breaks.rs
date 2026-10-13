//! First-difference designs of the six models of Banerjee and
//! Carrion-i-Silvestre (2015, Section 2) and the search for break dates
//! that minimise a unit's sum of squared residuals.
//!
//! A break position `p` is the 0-based calendar position of the last period
//! of the old regime: `DU_t = 1(t > p)` and the impulse is at `t = p + 1`.

use super::linalg::{cholesky, forward_solve, Qr};
use super::panel::Unit;

/// Model (1 to 6) and the leads and lags of the DOLS terms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Spec {
    pub model: u8,
    pub k: usize,
    /// Lags `p` and leads `q` of `D^2 x`; no DOLS terms when both are 0.
    pub dols_p: usize,
    pub dols_q: usize,
}

impl Spec {
    pub fn constant(&self) -> bool {
        matches!(self.model, 2 | 3 | 5 | 6)
    }

    /// Regime constants `DU_j` in first differences (trend shifts in levels).
    pub fn du(&self) -> bool {
        matches!(self.model, 3 | 6)
    }

    /// Shifts in the cointegrating vector.
    pub fn shift(&self) -> bool {
        self.model >= 4
    }

    pub fn dols(&self) -> bool {
        self.dols_p + self.dols_q > 0
    }

    fn n_dols(&self) -> usize {
        if self.dols() {
            self.dols_p + self.dols_q + 1
        } else {
            0
        }
    }

    /// First calendar row and number of rows of the first-difference sample
    /// of a unit with `t_len` periods (`None` if empty).
    pub fn rows(&self, t_len: usize) -> Option<(usize, usize)> {
        let (r0, end) = if self.dols() { (2 + self.dols_p, t_len.checked_sub(self.dols_q)?) } else { (1, t_len) };
        if end > r0 {
            Some((r0, end - r0))
        } else {
            None
        }
    }

    /// Columns common to all regimes: constant, `Dx`, DOLS terms.
    pub fn n_base(&self) -> usize {
        usize::from(self.constant()) + self.k * (1 + self.n_dols())
    }

    /// Number of columns with `m` breaks.
    pub fn ncol(&self, m: usize) -> usize {
        let per = 1 + usize::from(self.du()) + if self.shift() { self.k } else { 0 };
        self.n_base() + m * per
    }
}

/// Fill the first-difference design (`n x ncol`, column-major) on calendar
/// rows `r0..r0 + n`. Columns: constant, impulses, `DU_j`, `Dx`, `Dx DU_j`,
/// `D^2 x_{t-j}` for `j = -q..p`.
pub fn fill_design(u: &Unit, spec: &Spec, pos: &[usize], r0: usize, n: usize, z: &mut Vec<f64>) -> usize {
    let m = spec.ncol(pos.len());
    z.clear();
    z.resize(n * m, 0.0);
    let mut c = 0;
    if spec.constant() {
        z[..n].iter_mut().for_each(|v| *v = 1.0);
        c += 1;
    }
    for &p in pos {
        if p + 1 >= r0 && p + 1 < r0 + n {
            z[c * n + p + 1 - r0] = 1.0;
        }
        c += 1;
    }
    if spec.du() {
        for &p in pos {
            for r in 0..n {
                z[c * n + r] = if r0 + r > p { 1.0 } else { 0.0 };
            }
            c += 1;
        }
    }
    for j in 0..spec.k {
        for r in 0..n {
            z[c * n + r] = u.dx(j, r0 + r);
        }
        c += 1;
    }
    if spec.shift() {
        for &p in pos {
            for j in 0..spec.k {
                for r in 0..n {
                    let t = r0 + r;
                    z[c * n + r] = if t > p { u.dx(j, t) } else { 0.0 };
                }
                c += 1;
            }
        }
    }
    if spec.dols() {
        let (p, q) = (spec.dols_p as isize, spec.dols_q as isize);
        for l in -q..=p {
            for j in 0..spec.k {
                for r in 0..n {
                    let t = (r0 + r) as isize - l;
                    let t = t as usize;
                    z[c * n + r] = u.dx(j, t) - u.dx(j, t - 1);
                }
                c += 1;
            }
        }
    }
    debug_assert_eq!(c, m);
    m
}

/// `Dy_t` on rows `r0..r0 + n`.
pub fn fill_dy(u: &Unit, r0: usize, n: usize, dy: &mut Vec<f64>) {
    dy.clear();
    dy.extend((r0..r0 + n).map(|t| u.dy(t)));
}

/// Residuals of `target` on the design with breaks at `pos`.
#[allow(clippy::too_many_arguments)]
pub fn residuals(u: &Unit, spec: &Spec, pos: &[usize], r0: usize, n: usize, target: &[f64], qr: &mut Qr, z: &mut Vec<f64>) -> Result<Vec<f64>, usize> {
    let m = fill_design(u, spec, pos, r0, n, z);
    qr.factor(z, n, m, target).map_err(|e| e.column)?;
    let mut b = Vec::new();
    qr.coef(&mut b);
    Ok((0..n)
        .map(|r| target[r] - (0..m).map(|c| z[c * n + r] * b[c]).sum::<f64>())
        .collect())
}

/// Admissible break positions with `m` breaks: every regime spans at least
/// `h` periods and every impulse falls in the sample rows `r0..r0 + n`.
/// Lexicographic order.
pub fn candidates(t_len: usize, m: usize, h: usize, r0: usize, n: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    if h == 0 || t_len < (m + 1) * h {
        return out;
    }
    let ok = |p: usize| p + 1 >= r0 && p + 1 < r0 + n;
    match m {
        1 => {
            for p in h - 1..=t_len - 1 - h {
                if ok(p) {
                    out.push(vec![p]);
                }
            }
        }
        2 => {
            for p1 in h - 1..=t_len - 1 - 2 * h {
                if !ok(p1) {
                    continue;
                }
                for p2 in p1 + h..=t_len - 1 - h {
                    if ok(p2) {
                        out.push(vec![p1, p2]);
                    }
                }
            }
        }
        _ => {}
    }
    out
}

/// Prefix sums of `a_t a_t'`, `a_t = (w_t, target_t)`, over the sample rows
/// of one unit, where `w_t` holds the columns common to all regimes.
#[derive(Clone, Debug, Default)]
pub struct Gram {
    n: usize,
    a: usize,
    /// Positions of `Dx` within `w`.
    xi: Vec<usize>,
    rows: Vec<f64>,
    prefix: Vec<f64>,
    r0: usize,
}

impl Gram {
    pub fn new(u: &Unit, spec: &Spec, r0: usize, n: usize, target: &[f64]) -> Self {
        // w = (constant, Dx, DOLS terms): the design without breaks.
        let mut z = Vec::new();
        let k0 = fill_design(u, spec, &[], r0, n, &mut z);
        let a = k0 + 1;
        let mut rows = vec![0.0; n * a];
        for r in 0..n {
            for c in 0..k0 {
                rows[r * a + c] = z[c * n + r];
            }
            rows[r * a + k0] = target[r];
        }
        let mut prefix = vec![0.0; (n + 1) * a * a];
        for r in 0..n {
            let (done, next) = prefix.split_at_mut((r + 1) * a * a);
            let prev = &done[r * a * a..];
            let next = &mut next[..a * a];
            let row = &rows[r * a..(r + 1) * a];
            for i in 0..a {
                for j in 0..a {
                    next[i * a + j] = prev[i * a + j] + row[i] * row[j];
                }
            }
        }
        let first_x = usize::from(spec.constant());
        Gram { n, a, xi: (first_x..first_x + spec.k).collect(), rows, prefix, r0 }
    }

    /// Sum of `a_i a_j` over rows `lo..n`, excluding the deleted rows.
    #[inline]
    fn sum(&self, lo: usize, i: usize, j: usize, del: &[usize]) -> f64 {
        let a = self.a;
        let mut s = self.prefix[self.n * a * a + i * a + j] - self.prefix[lo * a * a + i * a + j];
        for &d in del {
            if d >= lo {
                s -= self.rows[d * a + i] * self.rows[d * a + j];
            }
        }
        s
    }

    /// Sum of squared residuals of the target on the design with breaks at
    /// `pos`, from the normal equations; `None` if a Cholesky pivot is too
    /// small (the caller falls back to the QR).
    pub fn ssr(&self, spec: &Spec, pos: &[usize], ws: &mut GramWork) -> Option<f64> {
        let k0 = self.a - 1;
        let GramWork { del, cols, g, v, scale } = ws;
        del.clear();
        del.extend(pos.iter().map(|&p| p + 1 - self.r0));
        // Column c of the reduced design: (column of a, first row).
        cols.clear();
        cols.extend((0..k0).map(|c| (c, 0)));
        if spec.shift() {
            for &d in del.iter() {
                cols.extend(self.xi.iter().map(|&x| (x, d)));
            }
        }
        let m = cols.len();
        g.clear();
        g.resize(m * m, 0.0);
        v.clear();
        v.resize(m, 0.0);
        for (ci, &(ai, li)) in cols.iter().enumerate() {
            for (cj, &(aj, lj)) in cols.iter().enumerate().take(ci + 1) {
                let s = self.sum(li.max(lj), ai, aj, del);
                g[cj * m + ci] = s;
                g[ci * m + cj] = s;
            }
            v[ci] = self.sum(li, ai, k0, del);
        }
        let yy = self.sum(0, k0, k0, del);
        scale.clear();
        scale.extend((0..m).map(|c| g[c * m + c].sqrt()));
        if scale.iter().any(|s| !(*s > 0.0)) {
            return None;
        }
        for i in 0..m {
            for j in 0..m {
                g[j * m + i] /= scale[i] * scale[j];
            }
            v[i] /= scale[i];
        }
        cholesky(g, m).ok()?;
        forward_solve(g, m, v);
        let fit: f64 = v.iter().map(|x| x * x).sum();
        Some((yy - fit).max(0.0))
    }
}

/// Buffers of `Gram::ssr`.
#[derive(Clone, Debug, Default)]
pub struct GramWork {
    del: Vec<usize>,
    cols: Vec<(usize, usize)>,
    g: Vec<f64>,
    v: Vec<f64>,
    scale: Vec<f64>,
}

/// Break positions minimising the sum of squared residuals of `target` over
/// the candidates (the earliest among ties), and that minimum.
pub fn search(u: &Unit, spec: &Spec, cands: &[Vec<usize>], r0: usize, n: usize, target: &[f64]) -> Option<(Vec<usize>, f64)> {
    let gram = Gram::new(u, spec, r0, n, target);
    let mut ws = GramWork::default();
    let mut z = Vec::new();
    let mut qr = Qr::default();
    let mut best: Option<(usize, f64)> = None;
    for (ci, pos) in cands.iter().enumerate() {
        let ssr = match gram.ssr(spec, pos, &mut ws) {
            Some(s) => Some(s),
            None => {
                let m = fill_design(u, spec, pos, r0, n, &mut z);
                qr.factor(&z, n, m, target).ok().map(|_| qr.rss())
            }
        };
        if let Some(s) = ssr {
            if best.map_or(true, |(_, b)| s < b) {
                best = Some((ci, s));
            }
        }
    }
    best.map(|(ci, s)| (cands[ci].clone(), s))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::panel::Panel;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }

    fn unit_panel(t: usize, k: usize, seed: &mut u64) -> Panel {
        let mut x = vec![0.0; t * k];
        for j in 0..k {
            let mut acc = 0.0;
            for r in 0..t {
                acc += lcg(seed);
                x[j * t + r] = acc;
            }
        }
        let mut acc = 0.0;
        let y: Vec<f64> = (0..t)
            .map(|r| {
                acc += lcg(seed);
                acc + 0.7 * x[r] + if r > t / 2 { 0.8 } else { 0.0 }
            })
            .collect();
        Panel::new(y, x, vec![0, t], k).unwrap()
    }

    #[test]
    fn gram_search_matches_brute_force_qr() {
        let mut seed = 12;
        for model in [1u8, 2, 4, 5] {
            for (dp, dq) in [(0, 0), (1, 1)] {
                for m in [1usize, 2] {
                    let panel = unit_panel(40, 2, &mut seed);
                    let u = panel.unit(0);
                    let spec = Spec { model, k: 2, dols_p: dp, dols_q: dq };
                    let (r0, n) = spec.rows(40).unwrap();
                    let mut dy = Vec::new();
                    fill_dy(&u, r0, n, &mut dy);
                    let cands = candidates(40, m, 6, r0, n);
                    assert!(!cands.is_empty());
                    let gram = Gram::new(&u, &spec, r0, n, &dy);
                    let mut ws = GramWork::default();
                    let mut z = Vec::new();
                    let mut qr = Qr::default();
                    for pos in &cands {
                        let a = gram.ssr(&spec, pos, &mut ws).unwrap();
                        let mm = fill_design(&u, &spec, pos, r0, n, &mut z);
                        qr.factor(&z, n, mm, &dy).unwrap();
                        let b = qr.rss();
                        assert!((a - b).abs() < 1e-10 * b, "model {model} m {m}: {a} vs {b}");
                    }
                    let (pos, _) = search(&u, &spec, &cands, r0, n, &dy).unwrap();
                    assert!(cands.contains(&pos));
                }
            }
        }
    }

    #[test]
    fn candidate_sets() {
        // T = 20, h = 3: one break at p = 2..=16.
        let c = candidates(20, 1, 3, 1, 19);
        assert_eq!(c.first().unwrap(), &vec![2]);
        assert_eq!(c.last().unwrap(), &vec![16]);
        let c2 = candidates(20, 2, 3, 1, 19);
        assert!(c2.iter().all(|p| p[1] - p[0] >= 3 && p[0] >= 2 && p[1] <= 16));
        assert_eq!(c2.len(), (2..=13).map(|p1| 16 - (p1 + 3) + 1).sum::<usize>());
    }
}
