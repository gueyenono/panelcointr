# Validation results for `pcoint_westerlund()`

Run on 2026-10-11 with panelcointr 0.0.0.9000 (release build), R 4.6.0 and
rustc 1.95.0, on 8 cores. The scripts in this directory and in `inst/bench/`
reproduce these results; seeds are set inside them. Rejection rates are in
percent at the 5% level.

## T7(c): Table 1 moments by simulation

`Rscript inst/validation/westerlund-table1.R 10000 1000` (the paper's 10,000
draws of random walks of length 1,000).

All 72 means (3 cases x 6 values of K x 4 statistics) agree with the literal
in `R/westerlund-moments.R` within Monte Carlo error: the mean |z| of the
differences is 0.94 and the largest is 2.7. The variances agree to a few
percent (median absolute relative difference 2%). The variances of Pt and
Pa come from the delta method and are the noisiest; the largest difference
is 20% (Case 1, K = 4).

## T7(a): size of the asymptotic tests (paper's Table 2)

`Rscript inst/validation/westerlund-size.R 2000` (Case 1, K = 1; "a": lags and
leads by AIC up to 2 (T/100)^(2/9); "b": fixed at that value). The Monte
Carlo standard error is 0.5 points. N = 20, T = 100:

| Design | Source | Gt a | Gt b | Ga a | Ga b | Pt a | Pt b | Pa a | Pa b |
|---|---|---|---|---|---|---|---|---|---|
| phi = 0 | panelcointr | 6.0 | 4.2 | 0.8 | 1.4 | 7.6 | 6.4 | 4.0 | 4.7 |
| | paper | 10.8 | 6.8 | 9.4 | 7.2 | 7.3 | 9.3 | 11.3 | 11.7 |
| phi = -0.4 | panelcointr | 32.3 | 12.1 | 31.9 | 18.6 | 21.7 | 11.8 | 25.2 | 16.7 |
| | paper | 20.7 | 21.2 | 36.4 | 42.4 | 17.8 | 19.8 | 31.6 | 34.8 |
| phi = 0.4 | panelcointr | 7.4 | 3.3 | 1.0 | 0.6 | 9.5 | 6.2 | 4.9 | 3.7 |
| | paper | 12.2 | 4.2 | 5.7 | 2.2 | 6.9 | 6.0 | 7.9 | 7.5 |
| gamma = 0.8 | panelcointr | 6.1 | 4.7 | 1.9 | 2.4 | 8.3 | 7.4 | 5.2 | 5.6 |
| | paper | 9.7 | 5.5 | 10.2 | 7.0 | 7.0 | 10.2 | 10.2 | 12.9 |

The distortions follow the paper's pattern (large with negative moving-average
errors, small with leads of the regressor), but Ga, and to a lesser extent
Pa, reject less often than in the paper. This comes from the initial level
of the data. In Case 1 (no deterministic terms) the statistics are not
invariant to it, and after the 50-period burn-in of the design the data
start far from zero (standard deviation sqrt(50)), whereas Table 1 assumes a
zero initial value. With fixed orders, N = 20, T = 100, phi = gamma = 0 and
1,000 replications:

| Burn-in | Deterministic | Gt | Ga | Pt | Pa |
|---|---|---|---|---|---|
| 50 | none | 4.5 | 1.4 | 7.0 | 5.1 |
| 0 | none | 4.0 | 6.4 | 6.7 | 12.0 |
| 50 | constant | 3.7 | 9.6 | 5.5 | 15.9 |
| 0 | constant | 3.9 | 10.2 | 5.3 | 14.4 |

With zero initial values the Case 1 sizes match the paper's (Ga 7.2, Pa 11.7
in Table 2), which suggests that its simulations started the observed sample
at zero; with a constant the initial level is irrelevant, as it should be.
`Rscript inst/validation/westerlund-size.R 2000 0` runs Table 2 without burn-in.

## T7(b): bootstrap under cross-sectional dependence (paper's Table 5)

`Rscript inst/validation/westerlund-bootstrap-size.R 300 4 [none|constant]`
(alpha = 0, factor loadings drawn from U(0, 2), which the paper does not
specify; lags and leads fixed at 2; 200 bootstrap replications; 300 Monte
Carlo replications, standard error 1.3 points). A star marks the bootstrap test.

Case 2 (constant), where the initial level does not matter:

| N | T | b | Gt | Gt* | Ga | Ga* | Pt | Pt* | Pa | Pa* |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 100 | 0 | 4.0 | 4.3 | 8.7 | 3.0 | 7.3 | 4.3 | 16.7 | 5.0 |
| 20 | 100 | 0 | 3.0 | 4.3 | 12.0 | 5.0 | 5.0 | 4.0 | 15.7 | 4.3 |
| 10 | 200 | 0 | 5.0 | 5.3 | 7.7 | 5.3 | 4.7 | 4.0 | 10.7 | 3.7 |
| 20 | 200 | 0 | 4.7 | 5.0 | 7.0 | 3.3 | 5.3 | 2.3 | 10.7 | 5.0 |
| 10 | 100 | 1 | 16.3 | 8.0 | 52.0 | 9.7 | 21.0 | 5.3 | 53.7 | 6.0 |
| 20 | 100 | 1 | 26.7 | 10.7 | 76.7 | 14.0 | 27.3 | 3.7 | 72.7 | 5.3 |
| 10 | 200 | 1 | 21.3 | 12.0 | 54.0 | 15.7 | 23.0 | 6.0 | 54.0 | 9.3 |
| 20 | 200 | 1 | 31.7 | 8.3 | 73.7 | 16.0 | 28.3 | 3.7 | 67.3 | 8.0 |

With cross-sectional dependence (b = 1) the asymptotic tests reject far too
often and the bootstrap brings the size close to the nominal level, with
some remaining excess for Ga*, as in the paper's Table 5 (Ga* 8.3 to 22.7).
In the paper's Case 1 design the bootstrap tests are conservative (Pa* 0.3
to 4.3; asymptotic tests 12 to 36 with b = 1): the bootstrap starts its
samples at zero, as the paper prescribes, while the simulated data do not
(see T7(a)).

## Section 9: performance

`Rscript inst/bench/westerlund.R 1 2 4` (medians):

| Workload | Target | Result |
|---|---|---|
| SJ example, AIC lags 1-3, leads 0-3, no bootstrap | < 5 ms | 3.7 ms (pure-R oracle: 78 ms) |
| SJ example, bootstrap 800, fixed lags | < 0.5 s | 0.067 s |
| SJ example, bootstrap 800, AIC re-selection (3 x 4 grid) | < 2 s | 0.165 s |
| N = 100, T = 200, K = 2, 5 x 5 grid, B = 999, re-selection | < 60 s | 10.3 s (1 thread), 5.3 s (2), 2.9 s (4) |

## For information: CRAN package Westerlund

`Westerlund::westerlund_test()` (version 0.1.3, also a port of `xtwest`)
gives the same statistics as `compat = "xtwest"` and
`compat = "xtwest_westerlund"` on the three runs of Persyn and Westerlund
(2008), to within 2e-12. It is not a dependency.

# Validation results for `pcoint_banerjee()`

Run on 2026-10-12 with panelcointr 0.0.0.9000 (release build), R 4.6.0 and
rustc 1.95.0, on 8 cores (7 threads). Numbering as in `plan/banerjee.md`,
Section 8. T2 to T6 are in `tests/testthat/test-banerjee*.R`.

## T1: published tables by simulation

`Rscript inst/validation/banerjee-tables.R 100000` (the papers' design:
random walks of 1,000 steps, 100,000 draws; 105 s). The z-scores compare our
draws with the printed values, allowing for the papers' own Monte Carlo error
(assumed 100,000 draws) and their rounding to three decimals.

| Table | Entries | max abs z | share abs z > 2 |
|---|---|---|---|
| BCS Remark 1 (c, tau moments) | 4 | 2.1 | 0.25 |
| BCS Table II, one break | 18 | 4.8 | 0.33 |
| BCS Table II, two breaks | 72 | 3.0 | 0.17 |
| BCS Table I (MQ, one break) | 162 | 5.1 | 0.31 |
| WP Table 9 (finite T) | 110 | 12.2 | 0.53 |
| DF quantiles vs MacKinnon (1996) | 6 | 2.0 | 0.17 |
| Bai and Ng (2004) Table I | 36 | 28.4 | 0.72 |

- Remark 1 and Table II are reproduced: the moments agree to the third
  decimal or within a few units of it, so the reading of the procedure in
  the plan is confirmed. The literals are used as printed.
- BCS Table I agrees within about 0.5 everywhere (the largest z-scores are
  for q = 5 and 6, where the quantiles are flattest).
- WP Table 9 (only used with `moments = "finite_sample"`) agrees within
  0.005 at T >= 250. At T = 50 and 100 our means are up to 0.022 higher (less
  negative) than the working paper's, e.g. -2.000 against -2.022 for gamma at
  lambda = 0.5 and T = 50; the working paper's finite-sample design is not
  described precisely enough to tell why.
- Bai and Ng (2004) Table I: our critical values agree within about 1 for
  q <= 4 but are 1.5 to 3 more negative for q = 5 and 6 in the detrended case
  (e.g. 5%, q = 5: -57.39 against -55.82). Their 10% value for q = 5,
  detrended, printed as -55.286, is a misprint: it is barely above the 5%
  value, and we simulate -52.61. The table is kept as printed (it only
  matters with `mq_level = 0.10` and five or more factors in Models 2 and 5).
- `R/sysdata.rda` holds quantiles of the c and tau individual statistics and
  of DF_mu and DF_tau from 10^6 draws (`data-raw/banerjee-tables.R`, 12 s).
  Their 5% quantiles are -1.943, -2.619, -2.865 and -3.418; MacKinnon's
  asymptotic values for the Dickey-Fuller cases are -1.941, -2.862 and
  -3.410.

## T7: size in the paper's Monte Carlo design (Section 5)

`Rscript inst/validation/banerjee-size.R 300` (N = 40, one factor, rho = 1,
phi = 1, a break at lambda = 0.5 with theta = 3 and, in Models 4 to 6,
delta_2 = 3; 102 s). Rejection rates at 5%; Monte Carlo standard error about
1.3 points.

| Model | Breaks | T | sigma2_F | Z, N(0,1) | Z*, Table III | factor ADF | mean r |
|---|---|---|---|---|---|---|---|
| 1 | none | 100 | 1 | 3.3 | | 4.3 | 1.00 |
| 2 | none | 100 | 1 | 9.0 | | 8.3 | 1.00 |
| 3 | known, 50 | 100 | 1 | 12.7 | | 7.7 | 1.00 |
| 6 | known, 50 | 100 | 1 | 11.0 | | 7.3 | 1.00 |
| 5 | known, 50 | 100 | 1 | 9.7 | | 7.7 | 1.00 |
| 2 | 1, heterogeneous | 100 | 1 | 8.7 | | 7.7 | 1.00 |
| 2 | 1, homogeneous | 50 | 0.5 | 33.0 | 22.0 | 10.7 | 1.00 |
| 2 | 1, homogeneous | 100 | 0.5 | 20.3 | 8.3 | 6.0 | 1.00 |
| 5 | 1, homogeneous | 50 | 0.5 | 94.0 | 80.0 | 20.0 | 1.99 |
| 5 | 1, homogeneous | 100 | 0.5 | 90.7 | 74.3 | | 2.00 |
| 2 | 1, homogeneous | 100 | 10 | 21.3 | 8.0 | 9.1 | 1.01 |
| 5 | 1, homogeneous | 100 | 10 | 87.0 | 69.7 | | 2.00 |

Findings (400 further replications isolated the causes):

- **Known or no breaks, and heterogeneous breaks.** The statistic is
  calibrated: with the ADF lag fixed at 0 (the true order), the size of Z in
  Model 2 is 5.5% (mean z 0.00, sd 1.01), and heterogeneous breaks in Model 5
  give 5.0%. The excess at the defaults (9% in Model 2) comes from the t-sig
  rule with k_max = 4 at T = 100, which shifts the mean of Z by -0.16; BIC
  lags give 6.8%. In Models 3 and 6 the asymptotic moments add about 0.13 to
  the shift at T = 100 (WP Table 9's T = 100 mean is 0.026 lower), so
  `moments = "finite_sample"` helps there.
- **Homogeneous breaks, Model 2.** Against the Table III critical values the
  size is 8% at T = 100 (the table is for N = 100; we use N = 40), close to
  the paper's Table IV. The N(0, 1) p-value of the infimum over-rejects
  badly, as the print method warns.
- **Homogeneous breaks, Model 5: does not reproduce the paper.** Table IV
  reports 6.5% to 9% for Z*_tau; we find 70% to 94%. At the true date Z is
  centred (mean 0.11), but at misplaced dates it is shifted down by about 1.3:
  a misplaced break in the cointegrating vector leaves `delta_2 Delta x` in
  part of the residuals, a variance shift in the cumulated errors that moves
  the ADF distribution to the left (a level shift alone moves it by only
  0.1). The infimum over the candidate dates then picks a wrong date (the
  true date in about 10% of the samples), and BIC3 finds a second factor in
  the common break spike. Fixing r = 1 and the lags does not remove it (67%).
  So the literal eq. 15 is not robust to large shifts of the cointegrating
  vector; known breaks or heterogeneous estimation (which minimises the SSR
  and finds the true date) have correct size. Documented in
  `?pcoint_banerjee`.

## Section 9: performance

`Rscript inst/bench/banerjee.R` (medians):

| Workload | Target | 1 thread | 2 threads | 4 threads |
|---|---|---|---|---|
| Known breaks, Model 5, m = 2, N = 63, T = 123, K = 2 | < 10 ms | 7.3 ms | | |
| Heterogeneous, same panel, m = 2 | < 1 s | 0.86 s | 0.48 s | 0.27 s |
| Homogeneous, Model 2, m = 1, N = 40, T = 250 | < 0.5 s | 0.18 s | | |
| Homogeneous, Model 5, m = 2, N = 40, T = 100 | < 3 s | 0.98 s | 0.51 s | 0.26 s |
| Gamma moments, 100,000 x 1000 | < 2 s | 1.75 s | 0.91 s | 0.46 s |

## T8: CRAN `xtbreakcoint` 1.0.4 (information only)

`Rscript inst/validation/banerjee-xtbreakcoint.R` (outputs only; the package
is GPL-3 and none of its code is used).

- **One unit, no factors, fixed lags (k = 0, 1, 2; 60 t-ratios).** Its
  "constant" and "trend" models both project the first differences on a
  constant and `Delta x`, the projection of our Model 2, and its ADF
  regression (own sample per lag, `RSS / (n - k - 1)`, cumulation without a
  zero start) is ours: the t-ratios agree to 9e-15.
- **N = 20, T = 60, heterogeneous slopes, no factors, k = 1.** The unit
  t-ratios hardly correlate (-0.25 to 0.44): `xtbreakcoint` pools the
  cointegrating vector over units (Appendix A, X2), so each unit keeps
  `(b_i - b) Delta x` and its own drift in the residuals. Its t-bar is the
  same for "constant" and "trend" (the same projection), around -0.5, close
  to the mean of the c distribution; in the trend model it is standardised
  with tau moments, which gives Z = +7.3 to +10.0 on data generated under
  H0 with a trend, where our Model 2 gives -1.0 to 0.9.
- **One common factor, both at their defaults.** Both find r = 1 (IC_p2 and
  BIC3). `xtbreakcoint` selects fewer lags (BIC: mean 0.05 to 0.15; t-sig:
  0.25 to 1.05). Z in the trend model is again +9.6 to +11.8 for
  `xtbreakcoint` and -0.5 to 0.8 for ours.

## Not done

- T9 (golden test on the PWT 6.3 panel of BCS 2017, Table 16): the data have
  not been obtained.
