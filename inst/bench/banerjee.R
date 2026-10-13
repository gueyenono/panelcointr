# Performance targets of plan/banerjee.md, Section 9. Not run by R CMD check.
# Needs a release build (R CMD INSTALL .) and the bench package. Usage (from
# the package root):
#   Rscript inst/bench/banerjee.R [threads for the parallel runs, default 1 2 4]

library(panelcointr)
args <- as.integer(commandArgs(TRUE))
threads <- if (length(args)) args else c(1L, 2L, 4L)

# Panel of the paper's application size: N = 63, T = 123, K = 2, one factor.
sim <- function(N, T, K, seed) {
  set.seed(seed)
  f <- cumsum(rnorm(T))
  do.call(rbind, lapply(seq_len(N), function(i) {
    X <- apply(matrix(rnorm(T * K), T), 2, cumsum)
    colnames(X) <- paste0("X", seq_len(K))
    y <- drop(X %*% rep(1, K)) + rnorm(1) * f + cumsum(rnorm(T)) + 2 * (seq_len(T) > T / 2)
    data.frame(id = i, time = 1900 + seq_len(T), y = y, X)
  }))
}
fml <- function(K) reformulate(paste0("X", seq_len(K)), "y")
show <- function(b) print(b[, c("expression", "min", "median")])

big <- sim(63, 123, 2, 1)
cat("Known breaks, Model 5, m = 2, N = 63, T = 123, K = 2 (target < 10 ms)\n")
show(bench::mark(pcoint_banerjee(fml(2), big, model = 5, breaks = c(1950, 1990)),
                 iterations = 20, check = FALSE))

cat("\nHeterogeneous estimated breaks, same panel, m = 2 (target < 1 s)\n")
for (th in threads) {
  b <- bench::mark(pcoint_banerjee(fml(2), big, model = 5, breaks = 2, n_threads = th),
                   iterations = 3, check = FALSE)
  cat("  threads =", th, ": median", format(b$median), "\n")
}

d40 <- sim(40, 250, 1, 2)
cat("\nHomogeneous, Model 2, m = 1, N = 40, T = 250 (target < 0.5 s)\n")
show(bench::mark(pcoint_banerjee(fml(1), d40, model = 2, breaks = 1, break_type = "homogeneous"),
                 iterations = 3, check = FALSE))

d100 <- sim(40, 100, 1, 3)
cat("\nHomogeneous, Model 5, m = 2, N = 40, T = 100 (target < 3 s; threads)\n")
for (th in threads) {
  b <- bench::mark(pcoint_banerjee(fml(1), d100, model = 5, breaks = 2, break_type = "homogeneous",
                                   n_threads = th), iterations = 2, check = FALSE)
  cat("  threads =", th, ": median", format(b$median), "\n")
}

cat("\nGamma moments by simulation, 100,000 x 1000 (target < 2 s)\n")
null <- getFromNamespace("banerjee_null_rs", "panelcointr")
for (th in threads) {
  b <- bench::mark(null("idio", "gamma", 0.37, 1000L, 100000L, 1L, 1, th), iterations = 2,
                   check = FALSE)
  cat("  threads =", th, ": median", format(b$median), "\n")
}
