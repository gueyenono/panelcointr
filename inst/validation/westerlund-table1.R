# T7(c): regenerate the asymptotic moments of Westerlund (2007, Table 1) by
# simulation, as in the paper's Remark 2: draws of K + 1 independent random
# walks of length T stand in for the Brownian motions V and W, and
#   C  = (int U^2, int U dV)',  C~ = (int U dV / int U^2, int U dV / sqrt(int U^2))',
# with U the residual of V on (d, W). Gt and Ga have mean and variance of the
# two elements of C~; Pt and Pa follow by the delta method from the mean and
# covariance of C.
#
# Not run by R CMD check. Usage (from the package root):
#   Rscript inst/validation/westerlund-table1.R [draws] [T]
# The paper uses 10000 draws and T = 1000. Prints the simulated moments, the
# Table 1 values and the difference in Monte Carlo standard errors.

args <- commandArgs(TRUE)
draws <- if (length(args) >= 1) as.integer(args[1]) else 10000L
T <- if (length(args) >= 2) as.integer(args[2]) else 1000L
set.seed(20071)

table1 <- get("westerlund_table1", envir = asNamespace("panelcointr"))

simulate_case <- function(case, K) {
  C <- matrix(NA_real_, draws, 2)
  Ct <- matrix(NA_real_, draws, 2)
  r <- (0:(T - 1)) / T # time at the start of each step
  D <- switch(case, NULL, rep(1, T), cbind(1, r))
  for (b in seq_len(draws)) {
    dV <- stats::rnorm(T) / sqrt(T)
    dW <- matrix(stats::rnorm(T * K), T) / sqrt(T)
    # Paths at the start of each step (V_0 = W_0 = 0): the integrands of the
    # stochastic integrals are evaluated before the increment.
    V <- c(0, cumsum(dV)[-T])
    W <- rbind(0, apply(dW, 2, cumsum)[-T, , drop = FALSE])
    U <- stats::lm.fit(cbind(D, W), V)$residuals
    iU2 <- mean(U^2)
    iUdV <- sum(U * dV)
    C[b, ] <- c(iU2, iUdV)
    Ct[b, ] <- c(iUdV / iU2, iUdV / sqrt(iU2))
  }
  Th <- colMeans(C)
  S <- stats::cov(C)
  phi_t <- c(-Th[2] / (2 * Th[1]^1.5), 1 / sqrt(Th[1]))
  phi_a <- c(-Th[2] / Th[1]^2, 1 / Th[1])
  # Monte Carlo standard errors of the means (delta method for Pt and Pa).
  c(mean_Gt = mean(Ct[, 2]), mean_Ga = mean(Ct[, 1]),
    mean_Pt = Th[2] / sqrt(Th[1]), mean_Pa = Th[2] / Th[1],
    var_Gt = stats::var(Ct[, 2]), var_Ga = stats::var(Ct[, 1]),
    var_Pt = drop(t(phi_t) %*% S %*% phi_t), var_Pa = drop(t(phi_a) %*% S %*% phi_a),
    se_Gt = stats::sd(Ct[, 2]) / sqrt(draws), se_Ga = stats::sd(Ct[, 1]) / sqrt(draws),
    se_Pt = sqrt(drop(t(phi_t) %*% S %*% phi_t) / draws),
    se_Pa = sqrt(drop(t(phi_a) %*% S %*% phi_a) / draws))
}

rows <- list()
for (case in 1:3) for (K in 1:6) {
  sim <- simulate_case(case, K)
  for (s in c("Gt", "Ga", "Pt", "Pa")) {
    rows[[length(rows) + 1]] <- data.frame(
      case = case, K = K, statistic = s,
      sim_mean = sim[[paste0("mean_", s)]], table_mean = table1$mean[[s]][case, K],
      z_mean = (sim[[paste0("mean_", s)]] - table1$mean[[s]][case, K]) / sim[[paste0("se_", s)]],
      sim_var = sim[[paste0("var_", s)]], table_var = table1$var[[s]][case, K]
    )
  }
  cat("case", case, "K", K, "done\n")
}
res <- do.call(rbind, rows)
res$rel_var <- res$sim_var / res$table_var - 1
options(width = 120)
print(res, digits = 4, row.names = FALSE)
cat("\nMean |z| of the means:", round(mean(abs(res$z_mean)), 2),
    "; max |relative difference| of the variances:", round(max(abs(res$rel_var)), 3), "\n")
