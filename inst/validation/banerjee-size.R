# T7 of plan/banerjee.md: size of the Banerjee and Carrion-i-Silvestre (2015)
# tests in the paper's Monte Carlo design (Section 5): N = 40 units,
#   y_it = mu + beta t + theta DU_t + delta_1 x_it + delta_2 DU_t x_it + u_it,
#   u_it = pi_i' F_t + e_it,  F_jt = phi F_j,t-1 + sigma_F w_jt,
#   e_it = rho e_i,t-1 + eps_it,  Delta x_it = v_it,
# with mu = 1, beta = 0.3, theta = 3, delta_1 = 1, delta_2 = 3 (Models 4-6,
# else 0), pi_ij ~ N(1, 1) and a break at lambda = 0.5. Under H0, rho = 1.
# Rejection rates at 5%: Z against the N(0, 1) quantile, Z* also against the
# paper's Table III, and the factor ADF test (t_Fd) when r = 1.
#
# Needs a release build (R CMD INSTALL .). Usage:
#   Rscript inst/validation/banerjee-size.R [reps] [threads]

library(panelcointr)
args <- commandArgs(TRUE)
reps <- if (length(args)) as.integer(args[1]) else 300L
threads <- if (length(args) > 1) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)

dgp <- function(N, T, model, r, sigma2_F, phi, rho, seed) {
  set.seed(seed)
  tb <- T / 2
  DU <- as.numeric(seq_len(T) > tb)
  F <- apply(matrix(rnorm(T * r, sd = sqrt(sigma2_F)), T), 2, function(w) stats::filter(w, phi, "recursive"))
  F <- matrix(F, T)
  do.call(rbind, lapply(seq_len(N), function(i) {
    x <- cumsum(rnorm(T))
    e <- as.numeric(stats::filter(rnorm(T), rho, "recursive"))
    trend <- if (model %in% c(2, 3, 5, 6)) 0.3 * seq_len(T) else 0
    d2 <- if (model >= 4) 3 else 0
    y <- 1 + trend + 3 * DU + x + d2 * DU * x + drop(F %*% rnorm(r, 1, 1)) + e
    data.frame(id = i, time = seq_len(T), y = y, x = x)
  }))
}

run <- function(model, breaks, break_type, T, r, sigma2_F, phi = 1, rho = 1, N = 40) {
  res <- t(vapply(seq_len(reps), function(b) {
    d <- dgp(N, T, model, r, sigma2_F, phi, rho, seed = 1000 * T + b)
    a <- list(y ~ x, data = d, model = model, breaks = breaks, n_threads = threads)
    if (!is.null(break_type)) a$break_type <- break_type
    fit <- suppressWarnings(do.call(pcoint_banerjee, a))
    t3 <- fit$breaks$table3
    c(z = fit$statistics$z < qnorm(0.05),
      z_t3 = if (is.null(t3)) NA else fit$statistics$z < t3[["5%"]],
      tf = if (fit$factors$r == 1) fit$factors$adf$p_value < 0.05 else NA,
      r = fit$factors$r)
  }, numeric(4)))
  data.frame(model = model, breaks = if (is.null(breaks)) "none" else
               paste(breaks, if (is.null(break_type)) "known" else break_type),
             T = T, r = r, sigma2_F = sigma2_F,
             size_Z = mean(res[, "z"]), size_Z_table3 = mean(res[, "z_t3"]),
             size_tF = mean(res[, "tf"], na.rm = TRUE), mean_r = mean(res[, "r"]))
}

t0 <- proc.time()[["elapsed"]]
out <- rbind(
  run(1, NULL, NULL, 100, 1, 1),
  run(2, NULL, NULL, 100, 1, 1),
  run(3, 50, NULL, 100, 1, 1),
  run(6, 50, NULL, 100, 1, 1),
  run(5, 50, NULL, 100, 1, 1),
  run(2, 1, "heterogeneous", 100, 1, 1),
  run(2, 1, "homogeneous", 50, 1, 0.5),
  run(2, 1, "homogeneous", 100, 1, 0.5),
  run(5, 1, "homogeneous", 50, 1, 0.5),
  run(5, 1, "homogeneous", 100, 1, 0.5),
  run(2, 1, "homogeneous", 100, 1, 10),
  run(5, 1, "homogeneous", 100, 1, 10)
)
print(out, digits = 3, row.names = FALSE)
cat(sprintf("\n%d replications, %d threads, %.0f s\n", reps, threads, proc.time()[["elapsed"]] - t0))
