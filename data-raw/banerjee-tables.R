# Quantile tables of the null distributions used by pcoint_banerjee() for
# the unit critical values and the factor ADF p-values (plan/banerjee.md,
# B19), simulated with the package's own engine as the paper does
# (Gaussian random walks of 1000 steps), with 10^6 replications each:
#   idio_c      ADF t-ratio without deterministic terms (Models 1 and 4)
#   idio_tau    the same after demeaning the differences (Models 2 and 5)
#   factor_c    DF t-ratio of the demeaned walk (DF_mu; factor test, 1 and 4)
#   factor_tau  DF t-ratio of the detrended walk (DF_tau; factor test, 2 and 5)
# Needs a release build: R CMD INSTALL . first, then run from the package
# root with Rscript data-raw/banerjee-tables.R. Writes R/sysdata.rda.

library(panelcointr)
null <- getFromNamespace("banerjee_null_rs", "panelcointr")
seed <- getFromNamespace("bcs_sim_seed", "panelcointr")
steps <- 1000L
reps <- 1000000L
threads <- max(1L, parallel::detectCores() - 1L)
probs <- c(0.0005, 0.001, 0.0025, 0.005, 0.0075, seq(0.01, 0.99, by = 0.005), 0.9925,
           0.995, 0.9975, 0.999, 0.9995)
draw <- function(kind, det) {
  t0 <- proc.time()[["elapsed"]]
  v <- null(kind, det, numeric(), steps, reps, 1L, seed, threads)
  message(kind, " ", det, ": ", round(proc.time()[["elapsed"]] - t0, 1), " s; mean ",
          round(mean(v), 4), ", variance ", round(var(v), 4))
  unname(quantile(v, probs, type = 8))
}
bcs_quantiles <- list(
  probs = probs,
  idio_c = draw("idio", "c"),
  idio_tau = draw("idio", "tau"),
  factor_c = draw("factor_df", "c"),
  factor_tau = draw("factor_df", "tau"),
  steps = steps,
  reps = reps
)
save(bcs_quantiles, file = "R/sysdata.rda", compress = "xz", version = 2)
