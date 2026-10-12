//! Stacked panel layout shared by the tests.

/// A panel stored unit by unit. `y` is stacked over units, `x` is
/// column-major with `n_rows()` rows and `k` columns, and the rows of unit
/// `i` are `offsets[i]..offsets[i + 1]`.
#[derive(Clone, Debug, Default)]
pub struct Panel {
    pub y: Vec<f64>,
    pub x: Vec<f64>,
    pub offsets: Vec<usize>,
    pub k: usize,
}

impl Panel {
    pub fn new(y: Vec<f64>, x: Vec<f64>, offsets: Vec<usize>, k: usize) -> Result<Self, String> {
        if offsets.len() < 2 {
            return Err("the panel has no units".into());
        }
        if offsets[0] != 0 || offsets.windows(2).any(|w| w[1] < w[0]) {
            return Err("unit offsets must start at 0 and be non-decreasing".into());
        }
        let n = *offsets.last().unwrap();
        if y.len() != n || x.len() != n * k {
            return Err(format!(
                "panel dimensions do not match: {} rows in the offsets, length(y) = {}, length(X) = {}",
                n,
                y.len(),
                x.len()
            ));
        }
        if y.iter().chain(x.iter()).any(|v| !v.is_finite()) {
            return Err("the panel contains non-finite values".into());
        }
        Ok(Panel { y, x, offsets, k })
    }

    pub fn n_units(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn n_rows(&self) -> usize {
        self.y.len()
    }

    pub fn t_len(&self, i: usize) -> usize {
        self.offsets[i + 1] - self.offsets[i]
    }

    pub fn unit(&self, i: usize) -> Unit<'_> {
        let (a, b) = (self.offsets[i], self.offsets[i + 1]);
        Unit {
            y: &self.y[a..b],
            x: &self.x,
            stride: self.n_rows(),
            offset: a,
            t_len: b - a,
            k: self.k,
        }
    }
}

/// View of one unit; time `t` runs over `0..t_len`.
#[derive(Clone, Copy, Debug)]
pub struct Unit<'a> {
    y: &'a [f64],
    x: &'a [f64],
    stride: usize,
    offset: usize,
    pub t_len: usize,
    pub k: usize,
}

impl<'a> Unit<'a> {
    #[inline]
    pub fn y(&self, t: usize) -> f64 {
        self.y[t]
    }

    #[inline]
    pub fn x(&self, j: usize, t: usize) -> f64 {
        self.x[j * self.stride + self.offset + t]
    }

    /// `y_t - y_{t-1}`, for `t >= 1`.
    #[inline]
    pub fn dy(&self, t: usize) -> f64 {
        self.y[t] - self.y[t - 1]
    }

    /// `x_{j,t} - x_{j,t-1}`, for `t >= 1`.
    #[inline]
    pub fn dx(&self, j: usize, t: usize) -> f64 {
        self.x(j, t) - self.x(j, t - 1)
    }
}
