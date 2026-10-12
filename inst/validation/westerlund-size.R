# T7(a): size of the asymptotic tests in the designs of Westerlund (2007,
# Table 2). Case 1, one regressor, DGP (11)-(12) under the null (alpha = 0):
#   Delta y_it = gamma (Delta x_i,t+1 + Delta x_it + Delta x_i,t-1) + e_it,
#   e_it = u_it + phi u_i,t-1,   Delta x_it = v_it,   u, v ~ N(0, 1),
# generated with T + 50 periods from zero, discarding the first 50. The tests
# use the bandwidth 4 (T/100)^(2/9) and either lags and leads chosen by AIC
# up to 2 (T/100)^(2/9) ("a") or fixed at that value ("b"); 2 (T/100)^(2/9)
# is rounded down. Rejection at the 5% level with the asymptotic p-values.
#
# In Case 1 the statistics depend on the initial level of the data, which
# after the burn-in is far from zero (standard deviation sqrt(50)), whereas
# the moments of Table 1 assume a zero initial value. The paper's Table 2 is
# matched closely with no burn-in (second argument 0), which suggests that
# its simulations started the observed sample at zero.
#
# Not run by R CMD check. Usage (from the package root):
#   Rscript inst/validation/westerlund-size.R [replications] [burn]
# The paper uses 2000 replications.

library(panelcointr)
args <- commandArgs(TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 2000L
burn <- if (length(args) >= 2) as.integer(args[2]) else 50L
set.seed(20072)

simulate <- function(N, T, phi, gamma) {
  TT <- T + burn
  do.call(rbind, lapply(seq_len(N), function(i) {
    dx <- stats::rnorm(TT + 2) # Delta x at times 0..TT+1
    u <- stats::rnorm(TT + 1)  # u at times 0..TT
    e <- u[-1] + phi * u[-(TT + 1)]
    dy <- gamma * (dx[3:(TT + 2)] + dx[2:(TT + 1)] + dx[1:TT]) + e
    keep <- (burn + 1):TT
    data.frame(id = i, t = seq_len(T), y = cumsum(dy)[keep], x = cumsum(dx[2:(TT + 1)])[keep])
  }))
}

reject <- function(d, rule, T) {
  m <- floor(2 * (T / 100)^(2 / 9))
  orders <- if (rule == "a") c(0, m) else m
  fit <- pcoint_westerlund(y ~ x, data = d, index = c("id", "t"), deterministic = "none",
                           lags = orders, leads = orders)
  fit$statistics$p_value < 0.05
}

designs <- rbind(
  expand.grid(N = c(10, 20), T = c(100, 200), phi = c(0, -0.4, 0.4), gamma = 0),
  expand.grid(N = c(10, 20), T = c(100, 200), phi = 0, gamma = c(0.4, 0.8))
)
# Westerlund (2007, Table 2), in percent: Gt a, Gt b, Ga a, Ga b, Pt a, Pt b, Pa a, Pa b.
paper <- rbind(
  c(9.8, 6.6, 8.6, 7.1, 7.9, 9.8, 13.1, 14.4), c(10.8, 6.8, 9.4, 7.2, 7.3, 9.3, 11.3, 11.7),
  c(7.7, 6.9, 8.7, 6.8, 7.2, 8.7, 13.8, 13.9), c(7.5, 5.0, 7.0, 5.3, 7.0, 6.8, 9.5, 11.2),
  c(16.3, 14.6, 28.9, 28.1, 14.3, 14.4, 26.5, 29.9), c(20.7, 21.2, 36.4, 42.4, 17.8, 19.8, 31.6, 34.8),
  c(12.5, 14.2, 24.9, 27.8, 11.6, 13.4, 25.2, 29.3), c(14.1, 21.1, 34.4, 41.8, 15.6, 17.3, 31.4, 33.9),
  c(8.3, 3.9, 5.4, 3.1, 6.3, 6.5, 9.3, 9.8), c(12.2, 4.2, 5.7, 2.2, 6.9, 6.0, 7.9, 7.5),
  c(6.7, 3.4, 5.3, 3.6, 6.3, 5.7, 10.2, 7.2), c(7.3, 3.2, 4.9, 1.5, 5.5, 4.8, 7.0, 6.2),
  c(10.6, 5.6, 10.0, 8.0, 10.4, 8.3, 15.6, 13.3), c(11.7, 6.9, 12.0, 7.7, 7.3, 9.7, 11.0, 12.5),
  c(8.4, 6.0, 9.1, 7.5, 8.8, 7.4, 14.3, 12.8), c(7.7, 5.7, 8.4, 6.3, 7.4, 7.9, 11.8, 12.0),
  c(9.2, 6.7, 10.9, 7.6, 8.5, 9.8, 13.4, 16.2), c(9.7, 5.5, 10.2, 7.0, 7.0, 10.2, 10.2, 12.9),
  c(8.1, 7.7, 8.1, 11.2, 8.0, 8.3, 13.5, 15.8), c(8.8, 5.0, 9.5, 5.9, 9.1, 7.4, 12.6, 11.1)
)
colnames(paper) <- c("Gt_a", "Gt_b", "Ga_a", "Ga_b", "Pt_a", "Pt_b", "Pa_a", "Pa_b")
# Rows above: five scenarios (phi = 0, -0.4, 0.4; gamma = 0.4, 0.8), each with
# (N, T) = (10, 100), (20, 100), (10, 200), (20, 200).
blocks <- split(seq_len(nrow(paper)), rep(1:5, each = 4))

out <- list()
for (s in seq_len(nrow(designs))) {
  g <- designs[s, ]
  ra <- rb <- matrix(FALSE, reps, 4)
  for (r in seq_len(reps)) {
    d <- simulate(g$N, g$T, g$phi, g$gamma)
    ra[r, ] <- reject(d, "a", g$T)
    rb[r, ] <- reject(d, "b", g$T)
  }
  ours <- 100 * c(colMeans(ra), colMeans(rb))[c(1, 5, 2, 6, 3, 7, 4, 8)]
  scen <- if (g$gamma == 0) match(g$phi, c(0, -0.4, 0.4)) else 3 + match(g$gamma, c(0.4, 0.8))
  pos <- (g$N == 20) + 2 * (g$T == 200) + 1
  out[[s]] <- data.frame(g, source = c("panelcointr", "paper"),
                         rbind(ours, paper[blocks[[scen]][pos], ]), row.names = NULL)
  cat(sprintf("N = %d, T = %d, phi = %.1f, gamma = %.1f done\n", g$N, g$T, g$phi, g$gamma))
}
res <- do.call(rbind, out)
cat("Burn-in:", burn, "periods\n")
options(width = 140)
print(res, digits = 3, row.names = FALSE)
cat("\nMonte Carlo standard error of a 5% rejection rate with", reps, "replications:",
    round(100 * sqrt(0.05 * 0.95 / reps), 2), "percentage points\n")
