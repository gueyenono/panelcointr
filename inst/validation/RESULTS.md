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
