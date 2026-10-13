//! Symmetric eigendecomposition: Householder reduction to tridiagonal form
//! followed by the implicit QL algorithm, as EISPACK's `tred2` and `tql2`
//! (in the formulation of the public-domain JAMA library). Also the
//! generalised problem `A v = mu B v` with `B` symmetric positive definite,
//! through the Cholesky factor of `B`.

use super::linalg::cholesky;

/// Eigenvalues in decreasing order and the eigenvectors as the columns of an
/// `n x n` column-major matrix.
#[derive(Clone, Debug, Default)]
pub struct SymEigen {
    pub n: usize,
    pub values: Vec<f64>,
    pub vectors: Vec<f64>,
}

impl SymEigen {
    pub fn vector(&self, j: usize) -> &[f64] {
        &self.vectors[j * self.n..(j + 1) * self.n]
    }
}

/// Eigendecomposition of the symmetric `n x n` matrix `a` (column-major;
/// only symmetry is assumed, both triangles are read).
pub fn sym_eigen(a: &[f64], n: usize) -> Result<SymEigen, String> {
    if a.len() != n * n {
        return Err("the matrix is not square".into());
    }
    if a.iter().any(|v| !v.is_finite()) {
        return Err("the matrix has non-finite elements".into());
    }
    if n == 0 {
        return Ok(SymEigen::default());
    }
    // v[i][j] of JAMA is v[i * n + j]; symmetrise to be safe.
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            v[i * n + j] = 0.5 * (a[j * n + i] + a[i * n + j]);
        }
    }
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    tred2(n, &mut v, &mut d, &mut e);
    tql2(n, &mut v, &mut d, &mut e)?;
    // tql2 sorts increasingly; return decreasing, vectors as columns.
    let mut values = Vec::with_capacity(n);
    let mut vectors = Vec::with_capacity(n * n);
    for j in (0..n).rev() {
        values.push(d[j]);
        vectors.extend((0..n).map(|i| v[i * n + j]));
    }
    Ok(SymEigen { n, values, vectors })
}

fn tred2(n: usize, v: &mut [f64], d: &mut [f64], e: &mut [f64]) {
    for j in 0..n {
        d[j] = v[(n - 1) * n + j];
    }
    for i in (1..n).rev() {
        let mut scale = 0.0;
        let mut h = 0.0;
        for dk in d.iter().take(i) {
            scale += dk.abs();
        }
        if scale == 0.0 {
            e[i] = d[i - 1];
            for j in 0..i {
                d[j] = v[(i - 1) * n + j];
                v[i * n + j] = 0.0;
                v[j * n + i] = 0.0;
            }
        } else {
            for dk in d.iter_mut().take(i) {
                *dk /= scale;
                h += *dk * *dk;
            }
            let mut f = d[i - 1];
            let mut g = h.sqrt();
            if f > 0.0 {
                g = -g;
            }
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for ej in e.iter_mut().take(i) {
                *ej = 0.0;
            }
            for j in 0..i {
                f = d[j];
                v[j * n + i] = f;
                g = e[j] + v[j * n + j] * f;
                for k in j + 1..i {
                    g += v[k * n + j] * d[k];
                    e[k] += v[k * n + j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for j in 0..i {
                e[j] /= h;
                f += e[j] * d[j];
            }
            let hh = f / (h + h);
            for j in 0..i {
                e[j] -= hh * d[j];
            }
            for j in 0..i {
                f = d[j];
                g = e[j];
                for k in j..i {
                    v[k * n + j] -= f * e[k] + g * d[k];
                }
                d[j] = v[(i - 1) * n + j];
                v[i * n + j] = 0.0;
            }
        }
        d[i] = h;
    }
    for i in 0..n - 1 {
        v[(n - 1) * n + i] = v[i * n + i];
        v[i * n + i] = 1.0;
        let h = d[i + 1];
        if h != 0.0 {
            for k in 0..=i {
                d[k] = v[k * n + i + 1] / h;
            }
            for j in 0..=i {
                let mut g = 0.0;
                for k in 0..=i {
                    g += v[k * n + i + 1] * v[k * n + j];
                }
                for k in 0..=i {
                    v[k * n + j] -= g * d[k];
                }
            }
        }
        for k in 0..=i {
            v[k * n + i + 1] = 0.0;
        }
    }
    for j in 0..n {
        d[j] = v[(n - 1) * n + j];
        v[(n - 1) * n + j] = 0.0;
    }
    v[(n - 1) * n + n - 1] = 1.0;
    e[0] = 0.0;
}

fn tql2(n: usize, v: &mut [f64], d: &mut [f64], e: &mut [f64]) -> Result<(), String> {
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    let mut f = 0.0;
    let mut tst1: f64 = 0.0;
    let eps = f64::EPSILON;
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n {
            if e[m].abs() <= eps * tst1 {
                break;
            }
            m += 1;
        }
        if m > l {
            let mut iter = 0;
            loop {
                iter += 1;
                if iter > 60 {
                    return Err("the eigenvalue iteration did not converge".into());
                }
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for di in d.iter_mut().skip(l + 2) {
                    *di -= h;
                }
                f += h;
                p = d[m];
                let mut c = 1.0;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0;
                let mut s2 = 0.0;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for k in 0..n {
                        h = v[k * n + i + 1];
                        v[k * n + i + 1] = s * v[k * n + i] + c * h;
                        v[k * n + i] = c * v[k * n + i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if !(e[l].abs() > eps * tst1) {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }
    // Selection sort into increasing order.
    for i in 0..n.saturating_sub(1) {
        let mut k = i;
        let mut p = d[i];
        for (j, &dj) in d.iter().enumerate().skip(i + 1) {
            if dj < p {
                k = j;
                p = dj;
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for j in 0..n {
                v.swap(j * n + i, j * n + k);
            }
        }
    }
    Ok(())
}

/// Eigenvalues (decreasing) and `B`-orthonormal eigenvectors of the
/// generalised problem `A v = mu B v`, `A` symmetric and `B` symmetric
/// positive definite: with `B = L L'`, the eigenvalues of
/// `L^{-1} A L^{-T}`. They are also the eigenvalues of `A B^{-1}`.
pub fn gen_sym_eigen(a: &[f64], b: &[f64], n: usize) -> Result<SymEigen, String> {
    if a.len() != n * n || b.len() != n * n {
        return Err("the matrices are not square".into());
    }
    let mut l = b.to_vec();
    cholesky(&mut l, n).map_err(|_| "the matrix B is not positive definite".to_string())?;
    // C = L^{-1} A L^{-T}: solve column by column.
    let mut w = a.to_vec(); // W = L^{-1} A
    for j in 0..n {
        lower_solve(&l, n, &mut w[j * n..(j + 1) * n]);
    }
    // C = L^{-1} W' (W' = A L^{-T} transposed back since A is symmetric).
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            c[j * n + i] = w[i * n + j];
        }
    }
    for j in 0..n {
        lower_solve(&l, n, &mut c[j * n..(j + 1) * n]);
    }
    let mut eg = sym_eigen(&c, n)?;
    // v = L^{-T} w.
    for j in 0..n {
        upper_solve_lt(&l, n, &mut eg.vectors[j * n..(j + 1) * n]);
    }
    Ok(eg)
}

/// Solve `L x = b` in place (`L` lower triangular, column-major in `l`).
fn lower_solve(l: &[f64], n: usize, b: &mut [f64]) {
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[k * n + i] * b[k];
        }
        b[i] = s / l[i * n + i];
    }
}

/// Solve `L' x = b` in place.
fn upper_solve_lt(l: &[f64], n: usize, b: &mut [f64]) {
    for i in (0..n).rev() {
        let mut s = b[i];
        for k in i + 1..n {
            s -= l[i * n + k] * b[k];
        }
        b[i] = s / l[i * n + i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }

    fn random_sym(n: usize, seed: &mut u64) -> Vec<f64> {
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = lcg(seed);
                a[j * n + i] = v;
                a[i * n + j] = v;
            }
        }
        a
    }

    fn check_decomposition(a: &[f64], n: usize, eg: &SymEigen, tol: f64) {
        let scale = a.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
        for j in 0..n {
            let v = eg.vector(j);
            for i in 0..n {
                let av: f64 = (0..n).map(|k| a[k * n + i] * v[k]).sum();
                assert!((av - eg.values[j] * v[i]).abs() <= tol * scale, "A v != mu v");
            }
            for k in 0..n {
                let dot: f64 = (0..n).map(|i| v[i] * eg.vector(k)[i]).sum();
                let target = if j == k { 1.0 } else { 0.0 };
                assert!((dot - target).abs() <= tol, "not orthonormal");
            }
        }
        assert!(eg.values.windows(2).all(|w| w[0] >= w[1]));
    }

    #[test]
    fn random_matrices() {
        let mut seed = 5;
        for &n in &[1usize, 2, 3, 7, 30, 120] {
            let a = random_sym(n, &mut seed);
            let eg = sym_eigen(&a, n).unwrap();
            check_decomposition(&a, n, &eg, 1e-11);
        }
    }

    #[test]
    fn repeated_and_zero_eigenvalues() {
        // Rank-2 Gram matrix with a repeated non-zero eigenvalue.
        let n = 6;
        let mut a = vec![0.0; n * n];
        let u = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let w = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        for i in 0..n {
            for j in 0..n {
                a[j * n + i] = 3.0 * (u[i] * u[j] + w[i] * w[j]);
            }
        }
        let eg = sym_eigen(&a, n).unwrap();
        assert!((eg.values[0] - 3.0).abs() < 1e-14 && (eg.values[1] - 3.0).abs() < 1e-14);
        assert!(eg.values[2..].iter().all(|v| v.abs() < 1e-14));
        check_decomposition(&a, n, &eg, 1e-12);
        let id: Vec<f64> = (0..16).map(|k| if k % 5 == 0 { 2.0 } else { 0.0 }).collect();
        let eg = sym_eigen(&id, 4).unwrap();
        assert!(eg.values.iter().all(|v| (v - 2.0).abs() < 1e-15));
    }

    #[test]
    fn generalised_problem() {
        let mut seed = 9;
        let n = 5;
        let a = random_sym(n, &mut seed);
        // B = M M' + I, positive definite.
        let m: Vec<f64> = (0..n * n).map(|_| lcg(&mut seed)).collect();
        let mut b = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                b[j * n + i] = (0..n).map(|k| m[k * n + i] * m[k * n + j]).sum::<f64>() + if i == j { 1.0 } else { 0.0 };
            }
        }
        let eg = gen_sym_eigen(&a, &b, n).unwrap();
        for j in 0..n {
            let v = eg.vector(j);
            for i in 0..n {
                let av: f64 = (0..n).map(|k| a[k * n + i] * v[k]).sum();
                let bv: f64 = (0..n).map(|k| b[k * n + i] * v[k]).sum();
                assert!((av - eg.values[j] * bv).abs() < 1e-12);
            }
        }
    }
}
