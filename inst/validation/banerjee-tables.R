# T1 of plan/banerjee.md: the published moments and critical values of
# Banerjee and Carrion-i-Silvestre (2015) and Bai and Ng (2004), re-simulated
# with the package's engine in the papers' design (random walks of 1000
# steps, 100,000 replications), and compared with the literals in
# R/banerjee-tables.R through z-scores of the differences.
#
# Needs a release build (R CMD INSTALL .). Usage:
#   Rscript inst/validation/banerjee-tables.R [reps] [threads]
# Prints the largest |z| per table; |z| above 4 deserves a look.

library(panelcointr)
args <- commandArgs(TRUE)
reps <- if (length(args)) as.integer(args[1]) else 100000L
threads <- if (length(args) > 1) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)
ns <- asNamespace("panelcointr")
null <- function(kind, det, lambda = numeric(), q = 1L, steps = 1000L, seed = 1L) {
  matrix(ns$banerjee_null_rs(kind, det, lambda, steps, reps, q, seed, threads), nrow = reps)
}

# z-score of a simulated mean or variance against a printed value that came
# from an independent simulation of `n_paper` draws, rounded to 3 decimals.
z_mean <- function(d, paper, n_paper = 1e5) {
  se <- sqrt(var(d) / length(d) + var(d) / n_paper + (0.0005^2) / 3)
  (mean(d) - paper) / se
}
z_var <- function(d, paper, n_paper = 1e5) {
  v <- var(d)
  m4 <- mean((d - mean(d))^4)
  se <- sqrt((m4 - v^2) / length(d) + (m4 - v^2) / n_paper + (0.0005^2) / 3)
  (v - paper) / se
}
# z-score of a simulated quantile, its standard error from the binomial
# variance of the empirical distribution function and a difference quotient
# of the quantile function.
z_quant <- function(d, p, paper, n_paper = 1e5) {
  q <- quantile(d, c(p - 0.005, p, p + 0.005), names = FALSE, type = 8)
  slope <- (q[3] - q[1]) / 0.01
  se <- slope * sqrt(p * (1 - p)) * sqrt(1 / length(d) + 1 / n_paper)
  (q[2] - paper) / sqrt(se^2 + (0.0005^2) / 3)
}

report <- function(name, z) {
  z <- unlist(z)
  cat(sprintf("%-36s %4d entries   max |z| = %.2f   share |z| > 2: %.2f\n",
              name, length(z), max(abs(z)), mean(abs(z) > 2)))
  invisible(z)
}

t0 <- proc.time()[["elapsed"]]
grid <- ns$bcs_lambda_grid

# (a) Remark 1.
dc <- null("idio", "c")
dt <- null("idio", "tau")
report("Remark 1 (c, tau)", list(
  z_mean(dc, ns$bcs_remark1$c[["mean"]]), z_var(dc, ns$bcs_remark1$c[["var"]]),
  z_mean(dt, ns$bcs_remark1$tau[["mean"]]), z_var(dt, ns$bcs_remark1$tau[["var"]])))

# (b) Table II, one and two breaks.
z1 <- lapply(seq_along(grid), function(i) {
  d <- null("idio", "gamma", grid[i])
  c(z_mean(d, ns$bcs_table2_one["mean", i]), z_var(d, ns$bcs_table2_one["var", i]))
})
report("Table II, one break", z1)
tab2 <- ns$bcs_table2_two
z2 <- list()
for (i in 1:8) for (j in (i + 1):9) {
  d <- null("idio", "gamma", grid[c(i, j)])
  z2[[length(z2) + 1]] <- c(z_mean(d, tab2[i, j]), z_var(d, tab2[j, i]))
}
report("Table II, two breaks", z2)

# (c) Table I: MQ with one trend break, q = 1..6, 1/5/10%.
lev <- c(0.01, 0.05, 0.10)
zt1 <- lapply(seq_along(grid), function(l) {
  d <- null("mq", "gamma", grid[l], q = 6L)
  outer(1:6, 1:3, Vectorize(function(q, k) z_quant(d[, q], lev[k], ns$bcs_table1[q, k, l])))
})
report("Table I (MQ, one break)", zt1)

# (d) WP Table 9: finite-sample moments.
wp <- ns$bcs_wp_table9
zw <- list()
for (k in seq_along(wp$T)) {
  TT <- wp$T[k]
  d <- null("idio", "c", steps = TT)
  zw[[length(zw) + 1]] <- c(z_mean(d, wp$c["mean", k]), z_var(d, wp$c["var", k]))
  d <- null("idio", "tau", steps = TT)
  zw[[length(zw) + 1]] <- c(z_mean(d, wp$tau["mean", k]), z_var(d, wp$tau["var", k]))
  for (i in seq_along(grid)) {
    d <- null("idio", "gamma", grid[i], steps = TT)
    zw[[length(zw) + 1]] <- c(z_mean(d, wp$gamma_mean[k, i]), z_var(d, wp$gamma_var[k, i]))
  }
}
zw <- report("WP Table 9 (finite T)", zw)

# (e) Asymptotic Dickey-Fuller quantiles of the factor ADF test against
# MacKinnon (1996, Table 1, T = infinity): 1%, 5%, 10%.
mk <- list(c = c(-3.43035, -2.86154, -2.56677), tau = c(-3.95877, -3.41049, -3.12705))
ze <- lapply(c("c", "tau"), function(j) {
  d <- null("factor_df", j)
  mapply(function(p, v) z_quant(d, p, v, n_paper = Inf), lev, mk[[j]])
})
report("DF quantiles vs MacKinnon (1996)", ze)
# Their finite-T (1000 steps) bias: these are not exactly the asymptotic values.

# (f) Bai and Ng (2004) Table I, demeaned and detrended.
sims <- list()
zf <- lapply(c("c", "tau"), function(j) {
  d <- null("mq", j, q = 6L)
  sims[[j]] <<- t(apply(d, 2, quantile, lev, type = 8))
  z <- outer(1:6, 1:3, Vectorize(function(q, k) z_quant(d[, q], lev[k], ns$bai_ng_table1[[j]][q, k])))
  dimnames(z) <- list(q = 1:6, level = lev)
  z
})
names(zf) <- c("c", "tau")
report("Bai and Ng (2004) Table I", zf)
for (j in c("c", "tau")) {
  cat("\nMQ critical values,", j, "- simulated (left) and Bai and Ng (2004) Table I (right)\n")
  print(round(cbind(sims[[j]], ns$bai_ng_table1[[j]]), 2))
}

cat(sprintf("\n%d replications, %d threads, %.0f s\n", reps, threads,
            proc.time()[["elapsed"]] - t0))
