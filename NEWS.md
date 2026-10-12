# panelcointr 0.0.0.9000

* New `pcoint_westerlund()`: the error-correction panel cointegration tests
  Gt, Ga, Pt and Pa of Westerlund (2007), with asymptotic p-values from the
  paper's Table 1 and bootstrap p-values robust to cross-sectional
  dependence. Lags and leads can be selected per unit by AIC or BIC. The
  computations, including the bootstrap, are written in Rust; the bootstrap
  runs in parallel and is reproducible with `set.seed()` for any number of
  threads.
* `compat = "xtwest"` and `compat = "xtwest_westerlund"` reproduce the
  Stata command `xtwest` and its `westerlund` option;
  `westerlund_control()` switches each difference from the paper separately.
* New `hce` data: health care expenditure and GDP in 20 OECD countries,
  1970-2001 (Westerlund 2007; Persyn and Westerlund 2008).
* Methods: `print()`, `summary()` (with the mean-group error-correction
  model), `plot()` (bootstrap distributions), `as.data.frame()`, and
  `tidy()` and `glance()` for the generics package.
