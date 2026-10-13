//! Null distributions of the statistics of Banerjee and Carrion-i-Silvestre
//! (2015, Theorem 1), simulated as the paper does: the procedure applied to
//! Gaussian random walks of `steps` periods (the paper uses 1000).
//!
//! * `Idio`: the ADF t-ratio (no lags, no deterministic terms) of the
//!   cumulated first differences after the model's first-difference
//!   projection: nothing (c), demeaning (tau), or demeaning within regimes
//!   with the impulse periods removed (gamma).
//! * `FactorDf`: the same series detrended on the constant (c), constant and
//!   trend (tau) or constant, trend and trend shifts (gamma), then the DF
//!   t-ratio.
//! * `Mq`: `T (v - 1)` for `q = 1..q_max` independent walks detrended as in
//!   `FactorDf`, where `v` is the smallest eigenvalue of
//!   `0.5 sum (Y_t Y_{t-1}' + Y_{t-1} Y_t') (sum Y_{t-1} Y_{t-1}')^{-1}`.
//!
//! Replicates come in blocks of `BLOCK`, block `b` drawing from stream `b`
//! of the seed, so the results do not depend on the number of threads.

use super::adf::df_t;
use super::mq::{deterministics, orthonormal_basis, project_out};
use super::parallel::par_map;
use super::rng::Rng;

pub const BLOCK: usize = 1000;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Idio,
    FactorDf,
    Mq(usize),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Det {
    C,
    Tau,
    Gamma,
}

/// Fixed design of a simulation.
struct Design {
    n: usize,
    det: Det,
    /// Break positions on the rows of the differences (row `s` is period
    /// `s + 1`): the impulse row and the regime of each row.
    impulse: Vec<bool>,
    regime: Vec<usize>,
    n_regimes: usize,
    /// Orthonormal basis of the factor detrending regressors.
    basis: Vec<f64>,
}

impl Design {
    fn new(steps: usize, det: Det, lambda: &[f64]) -> Result<Self, String> {
        if steps < 10 {
            return Err("at least 10 steps are needed".into());
        }
        let n = steps - 1;
        // Break: the last period of the old regime is round(lambda T) (1-based),
        // i.e. 0-based period p = round(lambda T) - 1 and difference row p.
        let mut pos: Vec<usize> = Vec::new();
        if det == Det::Gamma {
            for &l in lambda {
                let tb = (l * steps as f64).round() as usize;
                if tb < 2 || tb + 2 > steps {
                    return Err(format!("break fraction {l} is too close to the ends of the sample"));
                }
                pos.push(tb - 1);
            }
            if pos.windows(2).any(|w| w[1] <= w[0] + 1) {
                return Err("break fractions must be increasing and separated".into());
            }
        }
        // Row s of the differences is period t = s + 1.
        let mut impulse = vec![false; n];
        let mut regime = vec![0; n];
        for s in 0..n {
            let t = s + 1;
            regime[s] = pos.iter().filter(|&&p| t > p).count();
            impulse[s] = pos.iter().any(|&p| t == p + 1);
        }
        let n_regimes = pos.len() + 1;
        // Factor detrending on the n cumulated values (row s), with DT at the
        // row of the last old period: row p - 1 for period p.
        let rows_b: Vec<usize> = pos.iter().map(|&p| p - 1).collect();
        let det_x = match det {
            Det::C => deterministics(n, false, &[]),
            Det::Tau => deterministics(n, true, &[]),
            Det::Gamma => deterministics(n, true, &rows_b),
        };
        let basis = orthonormal_basis(&det_x, n, det_x.len() / n)?;
        Ok(Design { n, det, impulse, regime, n_regimes, basis })
    }

    /// First-difference projection and cumulation, in place.
    fn cumulate(&self, v: &mut [f64], sums: &mut Vec<(f64, usize)>) {
        match self.det {
            Det::C => {}
            Det::Tau => {
                let m = v.iter().sum::<f64>() / v.len() as f64;
                v.iter_mut().for_each(|x| *x -= m);
            }
            Det::Gamma => {
                sums.clear();
                sums.resize(self.n_regimes, (0.0, 0));
                for s in 0..self.n {
                    if !self.impulse[s] {
                        let e = &mut sums[self.regime[s]];
                        e.0 += v[s];
                        e.1 += 1;
                    }
                }
                for s in 0..self.n {
                    v[s] = if self.impulse[s] {
                        0.0
                    } else {
                        let (a, c) = sums[self.regime[s]];
                        v[s] - a / c as f64
                    };
                }
            }
        }
        let mut acc = 0.0;
        for x in v.iter_mut() {
            acc += *x;
            *x = acc;
        }
    }
}

fn mq_raw(y: &[f64], n: usize, q: usize) -> Result<f64, String> {
    let mut num = vec![0.0; q * q];
    let mut den = vec![0.0; q * q];
    for t in 1..n {
        for i in 0..q {
            let yi0 = y[i * n + t - 1];
            let yi1 = y[i * n + t];
            for j in 0..q {
                let yj0 = y[j * n + t - 1];
                let yj1 = y[j * n + t];
                num[j * q + i] += 0.5 * (yi1 * yj0 + yi0 * yj1);
                den[j * q + i] += yi0 * yj0;
            }
        }
    }
    let eg = super::eigen::gen_sym_eigen(&num, &den, q)?;
    Ok(n as f64 * (eg.values.last().unwrap() - 1.0))
}

/// Simulate `reps` draws. Returns a column-major `reps x ncol` matrix, with
/// `ncol = q_max` for `Mq` and 1 otherwise.
pub fn simulate(kind: Kind, det: Det, lambda: &[f64], steps: usize, reps: usize, seed: u64, n_threads: usize) -> Result<Vec<f64>, String> {
    let design = Design::new(steps, det, lambda)?;
    let ncol = match kind {
        Kind::Mq(q) => q,
        _ => 1,
    };
    if ncol == 0 {
        return Err("q_max must be positive".into());
    }
    let n_blocks = (reps + BLOCK - 1) / BLOCK;
    let blocks: Vec<Result<Vec<f64>, String>> = par_map(n_blocks, n_threads, |b| {
        let mut rng = Rng::stream(seed, b as u64);
        let count = BLOCK.min(reps - b * BLOCK);
        let n = design.n;
        let mut out = vec![0.0; count * ncol];
        let mut v = vec![0.0; n * ncol];
        let mut sums = Vec::new();
        for i in 0..count {
            for c in 0..ncol {
                let col = &mut v[c * n..(c + 1) * n];
                col.iter_mut().for_each(|x| *x = rng.normal());
                design.cumulate(col, &mut sums);
                if kind != Kind::Idio {
                    project_out(&design.basis, n, col);
                }
            }
            match kind {
                Kind::Idio | Kind::FactorDf => out[i] = df_t(&v),
                Kind::Mq(q) => {
                    for qq in 1..=q {
                        out[(qq - 1) * count + i] = mq_raw(&v[..n * qq], n, qq)?;
                    }
                }
            }
        }
        Ok(out)
    });
    let mut res = vec![0.0; reps * ncol];
    for (b, blk) in blocks.into_iter().enumerate() {
        let blk = blk?;
        let count = BLOCK.min(reps - b * BLOCK);
        for c in 0..ncol {
            res[c * reps + b * BLOCK..c * reps + b * BLOCK + count].copy_from_slice(&blk[c * count..(c + 1) * count]);
        }
    }
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mean_var(x: &[f64]) -> (f64, f64) {
        let n = x.len() as f64;
        let m = x.iter().sum::<f64>() / n;
        (m, x.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (n - 1.0))
    }

    #[test]
    fn idiosyncratic_moments_close_to_remark_1() {
        // 4000 draws of 400 steps: standard errors about 0.016 (c) and 0.009 (tau).
        let c = simulate(Kind::Idio, Det::C, &[], 400, 4000, 1, 2).unwrap();
        let (m, v) = mean_var(&c);
        assert!((m + 0.424).abs() < 0.06 && (v - 0.964).abs() < 0.1, "{m} {v}");
        let t = simulate(Kind::Idio, Det::Tau, &[], 400, 4000, 1, 2).unwrap();
        let (m, v) = mean_var(&t);
        assert!((m + 1.535).abs() < 0.04 && (v - 0.341).abs() < 0.04, "{m} {v}");
        let g = simulate(Kind::Idio, Det::Gamma, &[0.5], 400, 4000, 1, 2).unwrap();
        let (m, _) = mean_var(&g);
        assert!((m + 1.972).abs() < 0.04, "{m}");
    }

    #[test]
    fn independent_of_threads() {
        let a = simulate(Kind::Mq(2), Det::Tau, &[], 50, 2500, 3, 1).unwrap();
        let b = simulate(Kind::Mq(2), Det::Tau, &[], 50, 2500, 3, 3).unwrap();
        assert_eq!(a, b);
    }
}
