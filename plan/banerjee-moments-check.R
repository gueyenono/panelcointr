# Simulation check for plan/banerjee.md (section 0): the procedure as read in
# section 2 reproduces the null moments of Banerjee and Carrion-i-Silvestre
# (2015), Theorem 2 Remark 1 and Table II.
# Usage: Rscript plan/banerjee-moments-check.R [reps] [T]   (defaults 20000, 1000)
#
args <- commandArgs(TRUE)
reps <- if (length(args)) as.integer(args[1]) else 20000L
TT <- if (length(args) > 1) as.integer(args[2]) else 1000L
set.seed(20151)

tstat <- function(z) {               # z: projected differences, t = 2..T
  e <- cumsum(z)                      # e_2..e_T
  y <- diff(e); x <- e[-length(e)]
  b <- sum(x * y) / sum(x^2)
  r <- y - b * x
  s2 <- sum(r^2) / (length(y) - 1)
  b / sqrt(s2 / sum(x^2))
}

run <- function(cols) {               # cols: n x k matrix of deterministic columns (or NULL)
  n <- TT - 1
  P <- if (is.null(cols)) NULL else qr(cols)
  out <- numeric(reps)
  for (b in seq_len(reps)) {
    de <- rnorm(n)                    # Delta e under the null (random walk)
    z <- if (is.null(P)) de else qr.resid(P, de)
    out[b] <- tstat(z)
  }
  c(mean = mean(out), var = var(out))
}

n <- TT - 1
tt <- 2:TT                            # calendar index of the differenced rows
imp <- function(Tb) as.numeric(tt == Tb + 1)
du <- function(Tb) as.numeric(tt > Tb)
res <- list(
  c_nobreak         = run(NULL),
  c_level_0.5       = run(cbind(imp(TT / 2))),
  tau_nobreak       = run(cbind(rep(1, n))),
  tau_level_0.5     = run(cbind(1, imp(TT / 2))),
  gamma_0.5         = run(cbind(1, imp(TT / 2), du(TT / 2))),
  gamma_0.2         = run(cbind(1, imp(TT / 5), du(TT / 5))),
  gamma_0.3_0.7     = run(cbind(1, imp(0.3 * TT), imp(0.7 * TT), du(0.3 * TT), du(0.7 * TT)))
)
print(round(do.call(rbind, res), 3))
cat("Paper: c -0.424/0.964; tau -1.535/0.341; gamma(0.5) -1.972/0.339;",
    "gamma(0.2) -1.809/0.392; gamma(0.3,0.7) -2.315/0.337\n")
