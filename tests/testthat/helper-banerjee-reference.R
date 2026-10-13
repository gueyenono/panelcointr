# Pure-R oracle of the Banerjee and Carrion-i-Silvestre (2015) engine, written
# from the paper rather than from the Rust code: lm.fit projections, eigen()
# of the n x n Gram matrix, one lm.fit per ADF lag order, brute-force break
# searches and the MQ formulas of Bai and Ng (2004) with solve() and eigen().
# Break positions are 1-based indices of the last period of the old regime.

ref_bcs_units <- function(df, K) {
  lapply(split(df, df$id), function(d) {
    d <- d[order(d$time), ]
    list(y = d$y, X = as.matrix(d[paste0("x", seq_len(K))]))
  })
}

# First-difference sample rows (1-based periods) for T periods.
ref_rows <- function(T, dols) if (sum(dols) > 0) (3 + dols[1]):(T - dols[2]) else 2:T

ref_fd_design <- function(y, X, model, tb, dols) {
  T <- length(y)
  rows <- ref_rows(T, dols)
  dX <- rbind(NA, diff(X))
  cols <- list()
  if (model %in% c(2, 3, 5, 6)) cols$const <- rep(1, length(rows))
  for (b in seq_along(tb)) cols[[paste0("imp", b)]] <- as.numeric(rows == tb[b] + 1)
  if (model %in% c(3, 6)) for (b in seq_along(tb)) cols[[paste0("du", b)]] <- as.numeric(rows > tb[b])
  for (k in seq_len(ncol(X))) cols[[paste0("dx", k)]] <- dX[rows, k]
  if (model >= 4) {
    for (b in seq_along(tb)) for (k in seq_len(ncol(X))) {
      cols[[paste0("dxdu", b, "_", k)]] <- dX[rows, k] * (rows > tb[b])
    }
  }
  if (sum(dols) > 0) {
    d2X <- rbind(NA, diff(dX))
    for (l in -dols[2]:dols[1]) for (k in seq_len(ncol(X))) {
      cols[[paste0("d2x", l, "_", k)]] <- d2X[rows - l, k]
    }
  }
  list(Z = do.call(cbind, cols), dy = diff(y)[rows - 1], rows = rows)
}

ref_resid <- function(Z, v) {
  if (is.null(Z) || !ncol(Z)) return(v)
  fit <- lm.fit(Z, v)
  stopifnot(fit$rank == ncol(Z))
  fit$residuals
}

ref_project <- function(units, model, tb_list, dols) {
  sapply(seq_along(units), function(i) {
    d <- ref_fd_design(units[[i]]$y, units[[i]]$X, model, tb_list[[i]], dols)
    ref_resid(d$Z, d$dy)
  })
}

ref_criteria <- function(values, n, N, rmax, crit) {
  NT <- N * n
  V <- vapply(0:rmax, function(k) (sum(values) - sum(values[seq_len(k)])) / NT, numeric(1))
  s2 <- V[rmax + 1]
  C2 <- min(N, n)
  k <- 0:rmax
  switch(crit,
    pc1 = V + k * s2 * (N + n) / NT * log(NT / (N + n)),
    pc2 = V + k * s2 * (N + n) / NT * log(C2),
    pc3 = V + k * s2 * log(C2) / C2,
    ic1 = log(V) + k * (N + n) / NT * log(NT / (N + n)),
    ic2 = log(V) + k * (N + n) / NT * log(C2),
    ic3 = log(V) + k * log(C2) / C2,
    bic3 = V + k * s2 * (N + n - k) * log(NT) / NT)
}

# Principal components through the n x n Gram matrix (the Rust engine uses
# the smaller one).
ref_pca <- function(Y, r = NULL, crit = "bic3", rmax = 6, scale = FALSE) {
  n <- nrow(Y)
  N <- ncol(Y)
  sdv <- if (scale) apply(Y, 2, sd) else rep(1, N)
  Ys <- sweep(Y, 2, sdv, "/")
  e <- eigen(tcrossprod(Ys), symmetric = TRUE)
  values <- pmax(e$values[seq_len(min(n, N))], 0)
  crits <- ref_criteria(values, n, N, rmax, crit)
  if (is.null(r)) r <- which.min(crits) - 1L
  f <- sqrt(n) * e$vectors[, seq_len(r), drop = FALSE]
  for (j in seq_len(r)) if (f[which.max(abs(f[, j])), j] < 0) f[, j] <- -f[, j]
  L <- crossprod(Ys, f) / n
  resid <- (Ys - f %*% t(L))
  list(values = values, criteria = crits, r = r, f = f,
       loadings = L * sdv, resid = sweep(resid, 2, sdv, "*"))
}

# ADF without deterministic terms: t-ratio of e_{s-1} with k lags of De on
# rows s = s0..n (1-based).
ref_adf_fit <- function(e, k, s0 = k + 2) {
  n <- length(e)
  s <- s0:n
  de <- c(NA, diff(e))
  Z <- do.call(cbind, c(list(e[s - 1]), lapply(seq_len(k), function(j) de[s - j])))
  fit <- lm.fit(Z, de[s])
  rss <- sum(fit$residuals^2)
  df <- length(s) - ncol(Z)
  se <- sqrt(rss / df * diag(chol2inv(qr.R(fit$qr)))[order(fit$qr$pivot)])
  list(coef = fit$coefficients, se = se, rss = rss, n = length(s), resid = fit$residuals,
       t0 = fit$coefficients[1] / se[1], tlast = fit$coefficients[k + 1] / se[k + 1])
}

ref_adf <- function(e, kmin, kmax, rule = "tsig", cut = qnorm(0.95), sample = "common") {
  s0 <- if (sample == "common") kmax + 2 else NULL
  fitk <- function(k) ref_adf_fit(e, k, if (is.null(s0)) k + 2 else s0)
  k <- kmax
  if (kmin < kmax) {
    if (rule == "tsig") {
      while (k > kmin && abs(fitk(k)$tlast) < cut) k <- k - 1
    } else {
      ic <- vapply(kmin:kmax, function(kk) {
        f <- fitk(kk)
        pen <- if (rule == "aic") 2 else log(f$n)
        log(f$rss / f$n) + (kk + 1) * pen / f$n
      }, numeric(1))
      k <- (kmin:kmax)[which.min(ic)]
    }
  }
  f <- ref_adf_fit(e, k)
  list(t = unname(f$t0), k = k, n = f$n, resid = f$resid)
}

ref_detrend_basis <- function(n, j, tb_rows = integer()) {
  s <- seq_len(n)
  switch(j,
    c = cbind(rep(1, n)),
    tau = cbind(1, s),
    gamma = cbind(1, s, sapply(tb_rows, function(b) pmax(s - b, 0))))
}

ref_mq <- function(Fd, q, J, pmax) {
  n <- nrow(Fd)
  ev <- eigen(crossprod(Fd) / n^2, symmetric = TRUE)
  Y <- Fd %*% ev$vectors[, seq_len(q), drop = FALSE]
  Y0 <- Y[-n, , drop = FALSE]
  Y1 <- Y[-1, , drop = FALSE]
  xi <- Y1 - Y0 %*% solve(crossprod(Y0), crossprod(Y0, Y1))
  ne <- n - 1
  S1 <- matrix(0, q, q)
  for (j in seq_len(min(J, ne - 1))) {
    S1 <- S1 + (1 - j / (J + 1)) * crossprod(xi[1:(ne - j), , drop = FALSE],
                                             xi[(1 + j):ne, , drop = FALSE]) / n
  }
  num <- 0.5 * (crossprod(Y1, Y0) + crossprod(Y0, Y1) - n * (S1 + t(S1)))
  mq_c <- n * (min(Re(eigen(num %*% solve(crossprod(Y0)))$values)) - 1)
  dY <- diff(Y)
  nd <- n - 1
  pcap <- min(pmax, max(nd - q - 2, 0) %/% (q + 1))
  design <- function(p, d) do.call(cbind, lapply(seq_len(p), function(l) dY[d - l, , drop = FALSE]))
  bic <- vapply(0:pcap, function(p) {
    d <- (pcap + 1):nd
    E <- if (p == 0) dY[d, , drop = FALSE] else qr.resid(qr(design(p, d)), dY[d, , drop = FALSE])
    log(det(crossprod(E) / length(d))) + p * q^2 * log(length(d)) / length(d)
  }, numeric(1))
  p <- (0:pcap)[which.min(bic)]
  yt <- Y[(p + 1):n, , drop = FALSE]
  if (p > 0) {
    d <- (p + 1):nd
    B <- qr.coef(qr(design(p, d)), dY[d, , drop = FALSE])
    for (l in seq_len(p)) {
      yt <- yt - Y[(p + 1 - l):(n - l), , drop = FALSE] %*% B[((l - 1) * q + 1):(l * q), , drop = FALSE]
    }
  }
  m <- nrow(yt)
  a <- yt[-m, , drop = FALSE]
  b <- yt[-1, , drop = FALSE]
  mq_f <- n * (min(Re(eigen(0.5 * (crossprod(b, a) + crossprod(a, b)) %*% solve(crossprod(a)))$values)) - 1)
  list(mq_c = mq_c, mq_f = mq_f, p = p)
}

ref_factor_tests <- function(f, det, adf_args, J, pmax) {
  Fd <- apply(apply(f, 2, cumsum), 2, function(v) ref_resid(det, v))
  Fd <- matrix(Fd, nrow(f))
  r <- ncol(f)
  out <- list(fd = Fd)
  if (r == 1) out$adf <- do.call(ref_adf, c(list(Fd[, 1]), adf_args))
  if (r > 1) out$mq <- lapply(r:1, function(q) ref_mq(Fd, q, J, pmax))
  out
}

# The whole procedure for known break positions (list of 1-based vectors).
# `s` is the settings list of a fitted pcoint_banerjee object.
ref_bcs <- function(units, s, tb_list) {
  dols <- c(s$dols_p, s$dols_q)
  Y <- ref_project(units, s$model, tb_list, dols)
  pc <- ref_pca(Y, r = if (s$factors == "fixed") s$r_fixed else NULL,
                crit = if (s$factors == "fixed") "bic3" else s$factors,
                rmax = s$max_factors, scale = s$pca_scale)
  adf_args <- list(kmin = s$lag_min, kmax = s$lag_max, rule = s$lag_criterion, cut = s$tsig_cut,
                   sample = s$lag_sample)
  adf <- lapply(seq_len(ncol(Y)), function(i) do.call(ref_adf, c(list(cumsum(pc$resid[, i])), adf_args)))
  t <- vapply(adf, `[[`, numeric(1), "t")
  out <- list(Y = Y, pca = pc, t = t, lags = vapply(adf, `[[`, numeric(1), "k"),
              n = vapply(adf, `[[`, numeric(1), "n"), tbar = mean(t), adf = adf)
  if (pc$r > 0) {
    T <- length(units[[1]]$y)
    r0 <- ref_rows(T, dols)[1]
    j <- c("c", "tau", "gamma", "c", "tau", "gamma")[s$model]
    det <- ref_detrend_basis(nrow(Y), j, tb_list[[1]] - r0 + 1)
    out$factor <- ref_factor_tests(pc$f, det, adf_args, s$mq_bandwidth, s$mq_max_var)
  }
  out
}

# Admissible positions (1-based last period of the old regime) with m breaks.
ref_candidates <- function(T, m, h, dols) {
  rows <- ref_rows(T, dols)
  ok <- function(tb) all((tb + 1) %in% rows)
  if (m == 1) {
    c1 <- h:(T - h)
    return(lapply(c1[vapply(c1, ok, logical(1))], identity))
  }
  out <- list()
  for (a in h:(T - 2 * h)) for (b in (a + h):(T - h)) if (ok(c(a, b))) out[[length(out) + 1]] <- c(a, b)
  out
}

ref_break_search <- function(unit, model, m, h, dols, target = NULL) {
  T <- length(unit$y)
  cands <- ref_candidates(T, m, h, dols)
  ssr <- vapply(cands, function(tb) {
    d <- ref_fd_design(unit$y, unit$X, model, tb, dols)
    v <- if (is.null(target)) d$dy else target
    fit <- lm.fit(d$Z, v)
    if (fit$rank < ncol(d$Z)) Inf else sum(fit$residuals^2)
  }, numeric(1))
  best <- which.min(ssr)
  list(tb = cands[[best]], ssr = ssr[best], all = ssr, cands = cands)
}

ref_hetero <- function(units, s, m, h, r_fixed = NULL, crit = "bic3", rmax = 6) {
  dols <- c(s$dols_p, s$dols_q)
  N <- length(units)
  search <- function(C) {
    res <- lapply(seq_len(N), function(i) {
      d0 <- ref_fd_design(units[[i]]$y, units[[i]]$X, s$model, integer(), dols)
      target <- d0$dy - if (is.null(C)) 0 else C[, i]
      b <- ref_break_search(units[[i]], s$model, m, h, dols, target)
      d <- ref_fd_design(units[[i]]$y, units[[i]]$X, s$model, b$tb, dols)
      list(tb = b$tb, u = ref_resid(d$Z, target) + if (is.null(C)) 0 else C[, i])
    })
    list(tb = lapply(res, `[[`, "tb"), U = sapply(res, `[[`, "u"))
  }
  st <- search(NULL)
  r <- if (is.null(r_fixed)) ref_pca(st$U, crit = crit, rmax = rmax)$r else r_fixed
  ssr <- numeric()
  for (it in seq_len(s$iter_max)) {
    pc <- ref_pca(st$U, r = r, rmax = min(rmax, ncol(st$U) - 1))
    ssr[it] <- sum(pc$resid^2)
    if (it > 1 && (ssr[it - 1] - ssr[it]) / ssr[it - 1] < s$iter_tol) break
    if (r == 0) break
    st <- search(st$U - pc$resid)
  }
  list(tb = st$tb, ssr = ssr, r = r)
}

# Near-ties that make a discrete choice fragile: |t| near the t-sig cut-off.
ref_near_tie_adf <- function(e, kmin, kmax, cut, tol = 1e-6) {
  if (kmin >= kmax) return(FALSE)
  any(vapply(kmin:kmax, function(k) {
    if (k == 0) return(FALSE)
    abs(abs(ref_adf_fit(e, k, kmax + 2)$tlast) - cut) < tol
  }, logical(1)))
}
