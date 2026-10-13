# panelcointr 0.0.0.9000

* New `pcoint_banerjee()`: the panel cointegration tests of Banerjee and
  Carrion-i-Silvestre (2015) with common factors and structural breaks. Six
  models (breaks in the level, the trend and the cointegrating vector);
  breaks absent, known (common or unit-specific) or estimated (one or two
  per unit, heterogeneous by the iterative least-squares procedure or
  common by the infimum statistic); factors by principal components with
  the Bai and Ng (2002) criteria; ADF lags by the t-sig rule; DOLS terms for
  endogenous regressors. The common factors are tested by an ADF test or
  the MQ tests of Bai and Ng (2004), and the CD statistic checks the
  independence of the idiosyncratic errors. Moments and critical values
  come from the papers' tables, or are simulated for break fractions they do
  not cover. The computations, including the simulations, are in Rust and
  run in parallel with results independent of the number of threads.
* New `banerjee_control()` for the fine-grained settings.

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
