# T7(b): size of the asymptotic and bootstrap tests under cross-sectional
# dependence, in the design of Westerlund (2007, Table 5). Case 1, one
# regressor, alpha = gamma = 0:
#   Delta y_it = e_it = b lambda_i Delta f_t + u_it,  Delta x_it = v_it,
#   f_t, u_it, v_it ~ N(0, 1),
# T + 50 periods from zero, the first 50 discarded. Lags and leads fixed at
# 2 (T/100)^(2/9) (rounded down), bandwidth 4 (T/100)^(2/9), 200 bootstrap
# replications, 5% level. The paper does not give the distribution of the
# loadings lambda_i; this script draws them from U(0, 2) (an assumption).
#
# In Case 1 the statistics depend on the initial level of the data: after the
# burn-in, y and x start far from zero (standard deviation sqrt(50)), while
# the bootstrap starts its samples at zero, as the paper prescribes. The
# optional third argument runs the same design with a constant (Case 2),
# where the initial level is irrelevant.
#
# Not run by R CMD check. Usage (from the package root):
#   Rscript inst/validation/westerlund-bootstrap-size.R [replications] [threads] [deterministic]
# The paper's Table 5 is based on 300 replications and Case 1 ("none").

library(panelcointr)
args <- commandArgs(TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 300L
threads <- if (length(args) >= 2) as.integer(args[2]) else 1L
deterministic <- if (length(args) >= 3) args[3] else "none"
set.seed(20075)

simulate <- function(N, T, b, burn = 50L) {
  TT <- T + burn
  lambda <- stats::runif(N, 0, 2)
  df <- diff(c(0, stats::rnorm(TT))) # Delta f_t, f_0 = 0
  do.call(rbind, lapply(seq_len(N), function(i) {
    dy <- b * lambda[i] * df + stats::rnorm(TT)
    dx <- stats::rnorm(TT)
    keep <- (burn + 1):TT
    data.frame(id = i, t = seq_len(T), y = cumsum(dy)[keep], x = cumsum(dx)[keep])
  }))
}

designs <- expand.grid(N = c(10, 20), T = c(100, 200), b = c(0, 1))
# Westerlund (2007, Table 5), alpha = 0, in percent: Gt, Gt*, Ga, Ga*, Pt, Pt*, Pa, Pa*.
paper <- rbind(
  c(6.0, 4.0, 6.7, 3.7, 10.0, 6.3, 16.7, 7.7), c(5.3, 3.0, 5.7, 3.3, 6.3, 4.7, 12.0, 6.0),
  c(7.0, 4.7, 6.7, 5.0, 11.3, 7.0, 14.3, 7.7), c(7.0, 5.7, 8.0, 4.3, 8.7, 7.3, 12.7, 9.0),
  c(26.3, 8.0, 48.3, 12.3, 27.7, 7.0, 48.3, 10.0), c(32.0, 12.7, 56.3, 20.3, 30.0, 8.7, 55.0, 13.0),
  c(40.7, 6.0, 65.0, 8.3, 43.3, 5.0, 60.7, 8.7), c(44.3, 9.3, 72.7, 22.7, 40.3, 4.7, 59.3, 7.7)
)
colnames(paper) <- c("Gt", "Gt_boot", "Ga", "Ga_boot", "Pt", "Pt_boot", "Pa", "Pa_boot")
# Paper rows: within b, (N, T) = (10, 100), (10, 200), (20, 100), (20, 200).
paper_row <- function(N, T, b) 4 * b + 2 * (N == 20) + (T == 200) + 1

out <- list()
for (s in seq_len(nrow(designs))) {
  g <- designs[s, ]
  m <- floor(2 * (g$T / 100)^(2 / 9))
  rej <- matrix(NA, reps, 8)
  for (r in seq_len(reps)) {
    d <- simulate(g$N, g$T, g$b)
    fit <- pcoint_westerlund(y ~ x, data = d, index = c("id", "t"), deterministic = deterministic,
                             lags = m, leads = m, bootstrap = 200, n_threads = threads,
                             control = westerlund_control(keep_boot = FALSE))
    rej[r, ] <- c(rbind(fit$statistics$p_value, fit$statistics$p_boot)) < 0.05
  }
  out[[s]] <- data.frame(g, source = c("panelcointr", "paper"),
                         rbind(100 * colMeans(rej), paper[paper_row(g$N, g$T, g$b), ]),
                         row.names = NULL)
  cat(sprintf("N = %d, T = %d, b = %d done\n", g$N, g$T, g$b))
}
res <- do.call(rbind, out)
names(res)[-(1:4)] <- colnames(paper)
cat("Deterministic terms:", deterministic, "(the paper's values are for Case 1)\n")
options(width = 140)
print(res, digits = 3, row.names = FALSE)
cat("\nMonte Carlo standard error of a 5% rejection rate with", reps, "replications:",
    round(100 * sqrt(0.05 * 0.95 / reps), 2), "percentage points\n")
