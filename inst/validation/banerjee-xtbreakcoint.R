# T8 of plan/banerjee.md (information only): pcoint_banerjee() against CRAN
# xtbreakcoint 1.0.4, which departs from Banerjee and Carrion-i-Silvestre
# (2015) in most steps (plan, Appendix A). xtbreakcoint is GPL-3: only its
# outputs are compared, no code is used.
#
# 1. One unit, no factors, fixed lags: xtbreakcoint's "constant" and "trend"
#    models both project the first differences on a constant and Delta x,
#    which is the projection of our Model 2, and its ADF regression is ours.
#    The t-ratios should agree to rounding error.
# 2. A panel with heterogeneous slopes, no factors, fixed lags: xtbreakcoint
#    pools the cointegrating vector over units (Appendix A, X2) and has a
#    constant in Model 1, so t-bar differs.
# 3. Both packages at their defaults: different factor criteria (IC_p2 vs
#    BIC3), lag rules (BIC vs t-sig) and moments (T = 100 vs asymptotic).
#
# Usage: Rscript inst/validation/banerjee-xtbreakcoint.R (needs xtbreakcoint;
# a release build of panelcointr is faster but not required).

library(panelcointr)
library(xtbreakcoint)

sim <- function(N, T, seed, het = TRUE, r = 1) {
  set.seed(seed)
  F <- matrix(apply(matrix(rnorm(T * r), T), 2, cumsum), T)
  do.call(rbind, lapply(seq_len(N), function(i) {
    x <- cumsum(rnorm(T))
    b <- if (het) runif(1, 0.5, 1.5) else 1
    y <- 1 + 0.2 * seq_len(T) + b * x + drop(F %*% rnorm(r, 1, 1)) + cumsum(rnorm(T))
    data.frame(id = i, time = seq_len(T), y = y, x = x)
  }))
}
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

cat("1. One unit, no factors, fixed lags: individual t-ratios\n")
one <- do.call(rbind, lapply(1:20, function(s) {
  d <- sim(1, 60, s, r = 0)
  do.call(rbind, lapply(0:2, function(k) {
    ours <- pcoint_banerjee(y ~ x, data = d, index = c("id", "time"), model = 2,
                            factors = 0, lags = k)
    xc <- quiet(xtbreakcoint(y ~ x, data = d, id = "id", time = "time", model = "constant",
                             max_factors = 0, max_lag = k, lag_method = "fixed"))
    xt <- quiet(xtbreakcoint(y ~ x, data = d, id = "id", time = "time", model = "trend",
                             max_factors = 0, max_lag = k, lag_method = "fixed"))
    data.frame(seed = s, k = k, ours = ours$units$t_stat, constant = xc$adf_stats,
               trend = xt$adf_stats)
  }))
}))
cat(sprintf("   %d t-ratios; largest |difference|: %.2e (\"constant\"), %.2e (\"trend\")\n",
            nrow(one), max(abs(one$ours - one$constant)), max(abs(one$ours - one$trend))))

cat("\n2. N = 20, T = 60, heterogeneous slopes, no factors, fixed lags: t-bar\n")
two <- do.call(rbind, lapply(1:3, function(s) {
  d <- sim(20, 60, 100 + s, r = 0)
  do.call(rbind, lapply(1:2, function(m) {
    k <- 1
    ours <- pcoint_banerjee(y ~ x, data = d, index = c("id", "time"), model = m,
                            factors = 0, lags = k)
    xb <- quiet(xtbreakcoint(y ~ x, data = d, id = "id", time = "time", model = m,
                             max_factors = 0, max_lag = k, lag_method = "fixed"))
    data.frame(seed = s, model = m, tbar_ours = ours$statistics$value, tbar_xtb = xb$tbar,
               cor_t = cor(ours$units$t_stat, xb$adf_stats),
               z_ours = ours$statistics$z, z_xtb = xb$Z_t)
  }))
}))
print(two, digits = 3, row.names = FALSE)

cat("\n3. N = 20, T = 60, one common factor, both packages at their defaults\n")
three <- do.call(rbind, lapply(1:3, function(s) {
  d <- sim(20, 60, 200 + s)
  do.call(rbind, lapply(1:2, function(m) {
    ours <- pcoint_banerjee(y ~ x, data = d, index = c("id", "time"), model = m)
    xb <- quiet(xtbreakcoint(y ~ x, data = d, id = "id", time = "time", model = m))
    data.frame(seed = s, model = m, r_ours = ours$factors$r, r_xtb = xb$n_factors,
               lags_ours = mean(ours$units$lags), lags_xtb = mean(xb$lag_orders),
               tbar_ours = ours$statistics$value, tbar_xtb = xb$tbar,
               z_ours = ours$statistics$z, z_xtb = xb$Z_t)
  }))
}))
print(three, digits = 3, row.names = FALSE)
