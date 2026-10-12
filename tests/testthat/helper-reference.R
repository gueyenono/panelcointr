# Pure-R reference implementation (oracle) of pcoint_westerlund(), ported
# from plan/westerlund-prototype.R and extended to every setting of the
# package. It is deliberately naive and independent of the Rust engine:
# lm.fit() on NA-padded lag matrices, one regression per candidate model,
# explicit Frisch-Waugh-Lovell regressions for the panel statistics, and the
# long-run variance from its double-sum definition.

ref_lag <- function(v, k) { # k > 0: lag, k < 0: lead; NA padded
  n <- length(v)
  if (k == 0) return(v)
  if (abs(k) >= n) return(rep(NA_real_, n))
  if (k > 0) c(rep(NA_real_, k), v[seq_len(n - k)]) else c(v[(1 - k):n], rep(NA_real_, -k))
}

ref_lrv <- function(v, M) {
  n <- length(v)
  s <- 0
  for (j in -M:M) {
    w <- 1 - abs(j) / (M + 1)
    for (t in seq_len(n)) {
      u <- t - j
      if (u >= 1 && u <= n) s <- s + w * v[t] * v[u]
    }
  }
  s / n
}

ref_design <- function(y, X, p, q, dd) {
  T <- length(y)
  dy <- c(NA, diff(y))
  dX <- rbind(NA, diff(X))
  Z <- list()
  if (dd >= 1) Z$const <- rep(1, T)
  if (dd >= 2) Z$trend <- seq_len(T)
  for (k in seq_len(ncol(X))) Z[[paste0("lx", k)]] <- ref_lag(X[, k], 1)
  for (j in seq_len(p)) Z[[paste0("ldy", j)]] <- ref_lag(dy, j)
  for (j in -q:p) for (k in seq_len(ncol(X))) Z[[paste0("dx", k, "_", j)]] <- ref_lag(dX[, k], j)
  Z$ly <- ref_lag(y, 1)
  list(dy = dy, Z = do.call(cbind, Z))
}

ref_fit <- function(y, X, p, q, dd, rows) {
  d <- ref_design(y, X, p, q, dd)
  Z <- d$Z[rows, , drop = FALSE]
  dy <- d$dy[rows]
  stopifnot(!anyNA(Z), !anyNA(dy))
  f <- stats::lm.fit(Z, dy)
  stopifnot(f$rank == ncol(Z))
  list(d = d, Z = Z, dy = dy, coef = f$coefficients, rss = sum(f$residuals^2),
       n = length(rows), k = ncol(Z))
}

ref_ic <- function(f, p, q, T, lags, leads, dd, criterion) {
  m2ll <- f$n * (1 + log(2 * pi) + log(f$rss / f$n))
  switch(criterion,
    aic = m2ll + 2 * f$k,
    bic = m2ll + f$k * log(f$n),
    westerlund = log(f$rss / f$n) + 2 * (p + q + dd + 1) / (T - lags[2] - leads[2])
  )
}

ref_select <- function(y, X, lags, leads, dd, s) {
  if (lags[1] == lags[2] && leads[1] == leads[2]) return(list(p = lags[2], q = leads[2], grid = NULL))
  T <- length(y)
  grid <- expand.grid(q = leads[1]:leads[2], p = lags[1]:lags[2])[, c("p", "q")]
  grid$ic <- NA_real_
  for (g in seq_len(nrow(grid))) {
    p <- grid$p[g]
    q <- grid$q[g]
    rows <- if (s$ic_sample == "common") (lags[2] + 2):(T - leads[2]) else (p + 2):(T - q)
    grid$ic[g] <- ref_ic(ref_fit(y, X, p, q, dd, rows), p, q, T, lags, leads, dd, s$criterion)
  }
  visit <- if (s$ic_ties == "smaller") order(grid$p, grid$q) else order(-grid$p, -grid$q)
  best <- Inf
  choice <- NA
  for (g in visit) {
    if (grid$ic[g] < best) {
      best <- grid$ic[g]
      choice <- g
    }
  }
  list(p = grid$p[choice], q = grid$q[choice], grid = grid)
}

ref_unit <- function(y, X, p, q, dd, M, s) {
  T <- length(y)
  rows <- (p + 2):(T - q)
  f <- ref_fit(y, X, p, q, dd, rows)
  b <- f$coef
  df_den <- if (s$se_alpha == "df") f$n - f$k else f$n
  se <- sqrt(f$rss / df_den * solve(crossprod(f$Z))["ly", "ly"])

  # Composite error: Delta y minus the deterministic, level and lagged Delta y terms.
  keep <- !grepl("^dx", names(b))
  u_rows <- if (s$u_sample == "estimation") rows else (p + 2):T
  u <- f$d$dy[u_rows] - drop(f$d$Z[u_rows, keep, drop = FALSE] %*% b[keep])
  if (s$demean_u && dd >= 1) u <- u - mean(u)

  y_rows <- if (s$lrv_y_sample == "all") 2:T else rows
  v <- f$d$dy[y_rows]
  if (s$lrv_deterministics == "regression") {
    if (dd == 1) v <- v - mean(v)
    if (dd == 2) v <- stats::lm.fit(cbind(1, y_rows), v)$residuals
  } else if (s$lrv_deterministics == "implied" && dd == 2) {
    v <- v - mean(v)
  }
  a1 <- sqrt(ref_lrv(u, M) / ref_lrv(v, M))

  others <- colnames(f$Z) != "ly"
  yres <- stats::lm.fit(f$Z[, others, drop = FALSE], f$Z[, "ly"])$residuals
  dyres <- stats::lm.fit(f$Z[, others, drop = FALSE], f$dy)$residuals
  list(p = p, q = q, n = f$n, k = f$k, alpha = unname(b["ly"]), se = se, a1 = a1,
       rss = f$rss, s = sqrt(f$rss / (f$n - f$k)), sumy2 = sum(yres^2),
       sumydy = sum(yres * dyres))
}

# units: list of list(y, X); s: resolved settings (fit$settings);
# M: bandwidth per unit.
ref_westerlund <- function(units, s, M) {
  dd <- match(s$deterministic, c("none", "constant", "trend")) - 1
  N <- length(units)
  K <- ncol(units[[1]]$X)
  T_i <- vapply(units, function(u) length(u$y), numeric(1))
  sel <- lapply(units, function(u) ref_select(u$y, u$X, s$lags, s$leads, dd, s))
  G <- Map(function(u, sl, m) ref_unit(u$y, u$X, sl$p, sl$q, dd, m, s), units, sel, M)
  g <- function(L, f) unname(vapply(L, function(e) e[[f]], numeric(1)))
  T_i <- unname(T_i)
  ttil <- function(L) switch(s$t_norm, effective = g(L, "n"), df = g(L, "n") - g(L, "k"), T = T_i)

  Gt <- mean(g(G, "alpha") / g(G, "se"))
  w <- if (s$ga_weights == "rows") T_i else rep(1, N)
  Ga <- sum(w * ttil(G) * g(G, "alpha") / g(G, "a1")) / sum(w)

  pbar <- floor(mean(g(G, "p")))
  qbar <- floor(mean(g(G, "q")))
  P <- if (s$panel_lags == "individual") {
    G
  } else {
    Map(function(u, m) ref_unit(u$y, u$X, pbar, qbar, dd, m, s), units, M)
  }
  Tbar <- if (s$tbar == "units") {
    mean(ttil(P))
  } else {
    Tr <- sum(T_i^2) / sum(T_i)
    switch(s$t_norm,
      T = Tr,
      effective = Tr - pbar - qbar - 1,
      df = Tr - pbar - qbar - 1 - (dd + K + pbar + K * (pbar + qbar + 1) + 1)
    )
  }
  alpha <- sum(g(P, "sumydy") / g(P, "a1")) / sum(g(P, "sumy2"))
  sig <- if (s$sigma_df == "own") g(P, "s") else sqrt(g(P, "rss") / Tbar)
  se <- sqrt(mean((sig / g(P, "a1"))^2)) / sqrt(sum(g(P, "sumy2")))
  list(
    stats = c(Gt = Gt, Ga = Ga, Pt = alpha / se, Pa = Tbar * alpha),
    lags = g(G, "p"), leads = g(G, "q"), alpha = g(G, "alpha"), se = g(G, "se"),
    alpha1 = g(G, "a1"), sel = sel
  )
}

# Split a long data frame (columns id, time, y, x1..xK, sorted or not) into
# the units of the oracle, dropping nothing.
ref_units <- function(df, y = "y", x = grep("^x", names(df), value = TRUE)) {
  df <- df[order(df$id, df$time), ]
  lapply(split(df, df$id), function(d) list(y = d[[y]], X = as.matrix(d[x])))
}

# TRUE if some unit's two best information criteria are within `tol`
# (relative): the Gram-Cholesky selection could then differ from the QR one.
ref_near_tie <- function(ref, tol = 1e-8) {
  any(vapply(ref$sel, function(sl) {
    if (is.null(sl$grid)) return(FALSE)
    ic <- sort(sl$grid$ic)
    abs(ic[2] - ic[1]) <= tol * max(1, abs(ic[1]))
  }, logical(1)))
}

# Bootstrap reference (paper Section 5, extended to leads).

# Null-imposed first stage of one unit on S(p, q).
ref_boot_first_stage <- function(y, X, p, q, constant) {
  T <- length(y)
  dy <- c(NA, diff(y))
  dX <- rbind(NA, diff(X))
  Z <- list()
  if (constant) Z$const <- rep(1, T)
  for (j in seq_len(p)) Z[[paste0("ldy", j)]] <- ref_lag(dy, j)
  for (j in -q:p) for (k in seq_len(ncol(X))) Z[[paste0("dx", k, "_", j)]] <- ref_lag(dX[, k], j)
  Z <- do.call(cbind, Z)
  rows <- (p + 2):(T - q)
  f <- stats::lm.fit(Z[rows, , drop = FALSE], dy[rows])
  b <- f$coefficients
  list(a = unname(b[grepl("^ldy", names(b))]), g = unname(b[grepl("^dx", names(b))]),
       e = unname(f$residuals), dX = dX)
}

# Centred joint residual matrix over the common calendar periods. `start`:
# 0-based calendar position of each unit's first period.
ref_boot_w <- function(units, start, p, q, constant) {
  fs <- Map(function(u, pi, qi) ref_boot_first_stage(u$y, u$X, pi, qi, constant), units, p, q)
  T_i <- vapply(units, function(u) length(u$y), numeric(1))
  lo <- max(start + p + 1)
  hi <- min(start + T_i - 1 - q)
  cal <- lo:hi
  e <- vapply(seq_along(units), function(i) fs[[i]]$e[cal - start[i] - p[i]], numeric(length(cal)))
  dx <- do.call(cbind, lapply(seq_along(units), function(i) fs[[i]]$dX[cal - start[i] + 1, , drop = FALSE]))
  w <- cbind(e, dx)
  list(w = sweep(w, 2, colMeans(w)), first_stage = fs)
}

# Bootstrap panel for one replicate, from the state returned by
# westerlund_boot_state_rs() and 1-based indices `idx`.
ref_boot_panel <- function(st, idx, K) {
  N <- length(st$units)
  w <- matrix(st$w, st$n_periods)
  L <- length(idx)
  ys <- list()
  xs <- list()
  for (i in seq_len(N)) {
    b <- st$units[[i]]
    dx <- w[idx, N + (i - 1) * K + seq_len(K), drop = FALSE]
    e <- w[idx, i]
    need <- st$burn + b$start + b$t_len
    dy <- numeric(need)
    for (tau in seq_len(need)) {
      u <- e[tau]
      for (l in -b$q:b$p) {
        src <- tau - l
        if (src >= 1) u <- u + sum(b$g[(l + b$q) * K + seq_len(K)] * dx[src, ])
      }
      for (j in seq_len(b$p)) if (tau - j >= 1) u <- u + b$a[j] * dy[tau - j]
      dy[tau] <- u
    }
    win <- st$burn + b$start + seq_len(b$t_len)
    ys[[i]] <- cumsum(dy[win])
    xs[[i]] <- apply(dx[win, , drop = FALSE], 2, cumsum)
  }
  list(y = unlist(ys), x = as.vector(do.call(rbind, xs)))
}
