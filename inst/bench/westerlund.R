# Performance targets of plan/westerlund.md, Section 9. Not run by R CMD check.
# Needs a release build (R CMD INSTALL .), the bench package, and the package
# source tree for the pure-R oracle. Usage (from the package root):
#   Rscript inst/bench/westerlund.R [threads for the large panel, default 1 2 4]

library(panelcointr)
args <- as.integer(commandArgs(TRUE))
threads <- if (length(args)) args else c(1L, 2L, 4L)
source("tests/testthat/helper-reference.R")

sj <- function(...) {
  pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
                    deterministic = "trend", ...)
}

cat("SJ example, AIC lags 1-3, leads 0-3, no bootstrap (target < 5 ms)\n")
print(bench::mark(rust = sj(lags = c(1, 3), leads = c(0, 3)), iterations = 50,
                  check = FALSE)[, c("expression", "min", "median", "mem_alloc")])

cat("\nPure-R oracle on the same problem\n")
units <- lapply(split(hce, hce$ctr), function(d) list(y = d$loghex, X = as.matrix(d["loggdp"])))
s <- sj(lags = c(1, 3), leads = c(0, 3))$settings
print(bench::mark(oracle = ref_westerlund(units, s, rep(3L, 20)), iterations = 3,
                  check = FALSE)[, c("expression", "min", "median")])

cat("\nSJ example, bootstrap 800, fixed lags (target < 0.5 s)\n")
print(bench::mark(
  sj(lags = 1, leads = 1, bootstrap = 800, control = westerlund_control(keep_boot = FALSE)),
  iterations = 5, check = FALSE)[, c("min", "median")])

cat("\nSJ example, bootstrap 800, AIC re-selection over a 3 x 4 grid (target < 2 s)\n")
print(bench::mark(
  sj(lags = c(1, 3), leads = c(0, 3), bootstrap = 800,
     control = westerlund_control(keep_boot = FALSE)),
  iterations = 3, check = FALSE)[, c("min", "median")])

cat("\nN = 100, T = 200, K = 2, grid 5 x 5, B = 999, re-selection (target < 60 s)\n")
set.seed(1)
big <- do.call(rbind, lapply(1:100, function(i) {
  x <- apply(matrix(rnorm(400), 200), 2, cumsum)
  data.frame(id = i, t = 1:200, y = cumsum(rnorm(200)) + 0.5 * x[, 1], x1 = x[, 1], x2 = x[, 2])
}))
times <- vapply(threads, function(n) {
  set.seed(2)
  system.time(pcoint_westerlund(y ~ x1 + x2, big, index = c("id", "t"), lags = c(0, 4),
                                leads = c(0, 4), bootstrap = 999, n_threads = n,
                                control = westerlund_control(keep_boot = FALSE)))[["elapsed"]]
}, numeric(1))
print(data.frame(threads = threads, seconds = times, speedup = times[1] / times))
