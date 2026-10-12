# Pure-R prototype of xtwest (Persyn & Westerlund 2008, xtwest.ado v1.5).
#
# This is research code, not package code. It encodes the reading of
# xtwest.ado described in plan/westerlund.md and reproduces the published
# Stata Journal output (st0146) for the default (non-"westerlund") mode:
#
#   lags(1)   leads(1)   lrwindow(3) constant trend -> -2.681 -10.927 -12.035 -10.524  (exact)
#   lags(1 2) leads(0 1) lrwindow(2) constant trend -> -2.736 -11.254 -12.859 -11.773  (exact, AIC)
#
#   lags(1 3) leads(0 3) lrwindow(3) constant trend westerlund
#                                                   -> -4.082 -27.702 -12.969 -22.470  (exact, incl. z)
#
# It is the starting point for tests/testthat/helper-reference.R (the pure-R
# oracle the Rust engine is tested against).
#
# Usage:
#   d <- haven::read_dta("xtwestdata.dta")   # http://fmwww.bc.edu/repec/bocode/x/xtwestdata.dta
#   xtwest_proto(d, "loghex", "loggdp", id = "ctr", time = "year",
#                lags = 1, leads = 1, lrwindow = 3, constant = TRUE, trend = TRUE)

lag_vec <- function(v, k) { # k > 0: lag, k < 0: lead; NA padded
  n <- length(v)
  if (k == 0) return(v)
  if (k > 0) c(rep(NA, k), v[seq_len(n - k)]) else c(v[(1 - k):n], rep(NA, -k))
}

# xtwest's `lrvar ..., weighted nodemean`: Bartlett kernel, sums over available
# pairs, divided by the number of non-missing observations.
lrvar_bartlett <- function(v, M) {
  n <- sum(!is.na(v))
  v0 <- ifelse(is.na(v), 0, v)
  J <- sum(v0 * v0)
  for (l in seq_len(M)) J <- J + 2 * (1 - l / (M + 1)) * sum(v0 * lag_vec(v0, l), na.rm = TRUE)
  J / n
}

ecm_design <- function(y, X, tt, p, q, constant, trend) {
  dy <- c(NA, diff(y))
  dX <- apply(X, 2, function(x) c(NA, diff(x)))
  Z <- list()
  if (constant) Z$cons <- rep(1, length(y))
  if (trend) Z$trend <- tt
  Z$ly <- lag_vec(y, 1)
  for (k in seq_len(ncol(X))) Z[[paste0("lx", k)]] <- lag_vec(X[, k], 1)
  for (j in seq_len(p)) Z[[paste0("ldy", j)]] <- lag_vec(dy, j)
  for (k in seq_len(ncol(X))) for (j in -q:p) Z[[paste0("dx", k, "_", j)]] <- lag_vec(dX[, k], j)
  list(dy = dy, Z = do.call(cbind, Z))
}

fit_ecm <- function(y, X, tt, p, q, constant, trend) {
  des <- ecm_design(y, X, tt, p, q, constant, trend)
  ok <- stats::complete.cases(des$dy, des$Z)
  fit <- stats::lm.fit(des$Z[ok, , drop = FALSE], des$dy[ok])
  list(des = des, ok = ok, fit = fit, n = sum(ok), k = ncol(des$Z), rss = sum(fit$residuals^2))
}

xtwest_proto <- function(data, yvar, xvars, id, time, lags, leads = 0, lrwindow = 2,
                         constant = FALSE, trend = FALSE, westerlund = FALSE) {
  data <- data[order(data[[id]], data[[time]]), ]
  K <- length(xvars)
  pmin <- min(lags); pmax <- max(lags); qmin <- min(leads); qmax <- max(leads)
  auto <- pmin != pmax || qmin != qmax
  units <- lapply(split(data, data[[id]]), function(di) {
    list(y = di[[yvar]], X = as.matrix(di[, xvars]), tt = di[[time]], Ti = nrow(di))
  })

  # Lag/lead selection: xtwest loops max -> min and keeps the first strict minimum.
  # Each candidate is fitted on its own sample (sample size varies with p, q).
  sel <- t(vapply(units, function(u) {
    if (!auto) return(c(pmax, qmax))
    best <- Inf; bp <- NA; bq <- NA
    for (p in pmax:pmin) for (q in qmax:qmin) {
      f <- fit_ecm(u$y, u$X, u$tt, p, q, constant, trend)
      ic <- if (westerlund) {
        log(f$rss / (u$Ti - p - q - 1)) + 2 * (p + q + constant + trend + 1) / (u$Ti - pmax - qmax)
      } else { # Stata `estat ic` AIC
        f$n * (1 + log(2 * pi) + log(f$rss / f$n)) + 2 * f$k
      }
      if (ic < best) { best <- ic; bp <- p; bq <- q }
    }
    c(bp, bq)
  }, numeric(2)))

  mask_to_sample <- function(v, dy, p, q) ifelse(!is.na(lag_vec(dy, p)) & !is.na(lag_vec(dy, -q)), v, NA)

  unit_stats <- function(u, p, q) {
    f <- fit_ecm(u$y, u$X, u$tt, p, q, constant, trend)
    nm <- colnames(f$des$Z); b <- f$fit$coefficients; j <- which(nm == "ly")
    piv <- order(f$fit$qr$pivot)
    XtXi <- chol2inv(qr.R(f$fit$qr))[piv, piv]
    s2 <- f$rss / (f$n - f$k)
    dy <- f$des$dy
    dyt <- if (westerlund) mask_to_sample(dy, dy, p, q)
           else if (constant && trend) dy - mean(dy, na.rm = TRUE)
           else dy
    keep <- grepl("^(cons|trend|ly|lx|ldy)", nm)
    uu <- as.vector(dy - f$des$Z[, keep, drop = FALSE] %*% b[keep]) # = gamma(L)'dx + e
    if (westerlund) uu <- mask_to_sample(uu, dy, p, q)
    a1 <- sqrt(lrvar_bartlett(uu, lrwindow) / lrvar_bartlett(dyt, lrwindow))
    Zp <- f$des$Z[f$ok, nm != "ly", drop = FALSE]
    list(ai = unname(b[j]), se = sqrt(s2 * XtXi[j, j]), a1 = a1, n = f$n, df = f$n - f$k, rss = f$rss,
         dyres = stats::lm.fit(Zp, dy[f$ok])$residuals,
         yres = stats::lm.fit(Zp, f$des$Z[f$ok, "ly"])$residuals)
  }

  g <- function(L, f) vapply(L, `[[`, numeric(1), f)
  G <- Map(unit_stats, units, sel[, 1], sel[, 2])
  se <- g(G, "se"); if (westerlund) se <- se * sqrt(g(G, "df") / g(G, "n"))
  tn <- if (westerlund) g(G, "n") else g(G, "df")
  Gt <- mean(g(G, "ai") / se)
  Ga <- mean(tn * g(G, "ai") / g(G, "a1")) # balanced only; xtwest weights by rows otherwise

  # Panel statistics: xtwest uses floor(mean selected lag/lead) for every unit.
  mp <- floor(mean(sel[, 1])); mq <- floor(mean(sel[, 2]))
  P <- lapply(units, unit_stats, p = mp, q = mq)
  bot <- sum(vapply(P, function(r) sum(r$yres^2), numeric(1)))
  alpha <- sum(vapply(P, function(r) sum(r$yres * r$dyres) / r$a1, numeric(1))) / bot
  meanT <- mean(vapply(units, `[[`, numeric(1), "Ti"))
  tnp <- if (westerlund) meanT - mp - mq - 1
         else meanT - mp - mq - 1 - (constant + trend) - 1 - K - mp - K * (mp + mq + 1)
  S2 <- mean((sqrt(g(P, "rss") / tnp) / g(P, "a1"))^2)
  stats <- c(Gt = Gt, Ga = Ga, Pt = alpha / (sqrt(S2) / sqrt(bot)), Pa = tnp * alpha)

  case <- constant + trend + 1
  m <- if (westerlund) westerlund_moments_legacy(trend) else westerlund_moments_table1(case, K)
  N <- length(units)
  z <- c(Gt = sqrt(N) * (stats[["Gt"]] - m$mean[1]) / sqrt(m$var[1]),
         Ga = sqrt(N) * (stats[["Ga"]] - m$mean[2]) / sqrt(m$var[2]),
         Pt = (stats[["Pt"]] - sqrt(N) * m$mean[3]) / sqrt(m$var[3]),
         Pa = sqrt(N) * (stats[["Pa"]] - m$mean[4]) / sqrt(m$var[4]))
  list(stats = stats, z = z, p = stats::pnorm(z), lags = sel[, 1], leads = sel[, 2])
}

# Westerlund (2007) Table 1 (identical to the matrices hard-coded in xtwest.ado).
# Rows: case 1 (none), 2 (constant), 3 (constant + trend); columns: K = 1..6.
westerlund_moments_table1 <- function(case, K) {
  mean_gt <- rbind(c(-0.9763, -1.3816, -1.7093, -1.9789, -2.1985, -2.4262),
                   c(-1.7776, -2.0349, -2.2332, -2.4453, -2.6462, -2.8358),
                   c(-2.3664, -2.5284, -2.7040, -2.8639, -3.0146, -3.1710))
  mean_ga <- rbind(c(-3.8022, -5.8239, -7.8108, -9.8791, -11.7239, -13.8581),
                   c(-7.1423, -9.1249, -10.9667, -12.9561, -14.9752, -17.0673),
                   c(-12.0116, -13.6324, -15.5262, -17.3648, -19.2533, -21.2479))
  mean_pt <- rbind(c(-0.5105, -0.9370, -1.3169, -1.6167, -1.8815, -2.1256),
                   c(-1.4476, -1.7131, -1.9206, -2.1484, -2.3730, -2.5765),
                   c(-2.1124, -2.2876, -2.4633, -2.6275, -2.7858, -2.9537))
  mean_pa <- rbind(c(-1.0263, -2.4988, -4.2699, -6.1141, -8.0317, -10.0074),
                   c(-4.2303, -5.8650, -7.4599, -9.3057, -11.3152, -13.3180),
                   c(-8.9326, -10.4874, -12.1672, -13.8889, -15.6815, -17.6515))
  var_gt <- rbind(c(1.0823, 1.0981, 1.0489, 1.0576, 1.0351, 1.0409),
                  c(0.8071, 0.8481, 0.8886, 0.9119, 0.9083, 0.9236),
                  c(0.6603, 0.7070, 0.7586, 0.8228, 0.8477, 0.8599))
  var_ga <- rbind(c(20.6868, 29.9016, 39.0109, 50.5741, 58.9595, 69.5967),
                  c(29.6336, 39.3428, 49.4880, 58.7035, 67.9499, 79.1093),
                  c(46.2420, 53.7428, 64.5591, 74.7403, 84.7990, 94.0024))
  var_pt <- rbind(c(1.3624, 1.7657, 1.7177, 1.6051, 1.4935, 1.4244),
                  c(0.9885, 1.0663, 1.1168, 1.1735, 1.1684, 1.1589),
                  c(0.7649, 0.8137, 0.8857, 0.9985, 0.9918, 0.9898))
  var_pa <- rbind(c(8.3827, 24.0223, 39.8827, 53.4518, 63.2406, 76.6757),
                  c(19.7090, 31.2637, 42.9975, 57.4844, 69.4374, 81.0384),
                  c(37.5948, 45.6890, 57.9985, 74.1258, 81.3934, 91.2392))
  list(mean = c(mean_gt[case, K], mean_ga[case, K], mean_pt[case, K], mean_pa[case, K]),
       var = c(var_gt[case, K], var_ga[case, K], var_pt[case, K], var_pa[case, K]))
}

# Moments xtwest uses under its `westerlund` option (K = 1 only). They differ
# from Table 1 and presumably come from Westerlund's original GAUSS code.
westerlund_moments_legacy <- function(trend) {
  if (trend) list(mean = c(-2.356, -11.8978, -2.1128, -8.9536), var = c(0.6450, 44.2471, 0.7371, 35.6802))
  else list(mean = c(-1.793, -7.2014, -1.4746, -4.3559), var = c(0.7904, 29.3677, 1.0262, 21.0535))
}
