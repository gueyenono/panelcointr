# Small simulations keep the tests fast (Models 3 and 6 simulate their
# critical values).
bcs_ctl <- function(...) banerjee_control(sim_reps = 1000L, sim_steps = 100L, ...)

bcs_fit <- function(df, K, ...) {
  pcoint_banerjee(xfml(K), data = df, index = c("id", "time"), ...)
}

# Break positions of a fit as a list of 1-based vectors.
fit_tb <- function(fit) {
  p <- fit$breaks$positions
  lapply(seq_len(fit$N), function(i) if (is.null(p)) integer() else as.integer(stats::na.omit(p[i, ])))
}

# Discrete choices close to a tie make the comparison fragile.
bcs_near_tie <- function(ref, s) {
  crit <- ref$pca$criteria
  if (s$factors != "fixed" && sum(abs(crit - min(crit)) < 1e-8 * max(1, abs(min(crit)))) > 1) {
    return(TRUE)
  }
  if (s$lag_criterion == "tsig") {
    for (i in seq_len(ncol(ref$Y))) {
      if (ref_near_tie_adf(cumsum(ref$pca$resid[, i]), s$lag_min, s$lag_max, s$tsig_cut)) return(TRUE)
    }
  }
  FALSE
}

# T2: the Rust engine against the pure-R oracle.
bcs_compare <- function(df, K, args, tol = 1e-10) {
  fit <- do.call(bcs_fit, c(list(df, K, control = bcs_ctl()), args))
  ref <- ref_bcs(ref_bcs_units(df, K), fit$settings, fit_tb(fit))
  if (bcs_near_tie(ref, fit$settings)) return(invisible(FALSE))
  expect_equal(fit$factors$r, ref$pca$r)
  expect_equal(fit$units$lags, ref$lags)
  expect_equal(fit$units$n, ref$n)
  expect_equal(fit$units$t_stat, ref$t, tolerance = tol)
  expect_equal(fit$statistics$value, ref$tbar, tolerance = tol)
  expect_equal(unname(fit$factors$eigenvalues[seq_along(ref$pca$values)]), ref$pca$values,
               tolerance = tol)
  expect_equal(unname(fit$factors$criteria), ref$pca$criteria, tolerance = tol)
  if (fit$factors$r > 0) {
    expect_equal(unname(fit$factors$f), ref$pca$f, tolerance = 1e-8)
    expect_equal(unname(fit$factors$F_detrended), ref$factor$fd, tolerance = 1e-8)
    if (fit$factors$r == 1) {
      expect_equal(fit$factors$adf$statistic, ref$factor$adf$t, tolerance = 1e-8)
    } else {
      expect_equal(fit$factors$mq$MQ_c, vapply(ref$factor$mq, `[[`, numeric(1), "mq_c"), tolerance = 1e-7)
      expect_equal(fit$factors$mq$MQ_f, vapply(ref$factor$mq, `[[`, numeric(1), "mq_f"), tolerance = 1e-7)
    }
  }
  invisible(TRUE)
}

test_that("known and absent breaks match the oracle in every model (T2)", {
  configs <- list(
    list(model = 1),
    list(model = 2, factors = 2, lags = 1),
    list(model = 1, dols = c(1, 1), factors = 1),
    list(model = 2, factors = "ic2", lags = c(0, 3)),
    list(model = 1, breaks = 15),
    list(model = 2, breaks = c(12, 26), factors = 0),
    list(model = 3, breaks = 20, factors = 1),
    list(model = 3, breaks = c(12, 28), lags = 0),
    list(model = 4, breaks = 20),
    list(model = 5, breaks = c(14, 27), factors = 2, dols = 1),
    list(model = 6, breaks = 20, factors = 2),
    list(model = 6, breaks = c(13, 27), dols = c(0, 1)),
    list(model = 4, breaks = list(`1` = 15, `2` = c(10, 25), `3` = 20, `4` = integer(),
                                  `5` = 30, `6` = 18)),
    list(model = 5, breaks = data.frame(id = 1:6, b = c(12, 14, 16, 18, 20, 22)), factors = 1),
    list(model = 2, breaks = list(`1` = 15, `2` = 16, `3` = 20, `4` = 21, `5` = 30, `6` = 18),
         lags = 2)
  )
  for (K in 1:3) {
    for (cfg in configs) {
      df <- sim_panel(N = 6, T = 40, K = K, seed = K * 100 + cfg$model)
      bcs_compare(df, K, cfg)
    }
  }
})

test_that("the lag rules, lag samples and pca_scale match the oracle (T2)", {
  df <- sim_panel(N = 7, T = 45, K = 2, seed = 9)
  for (crit in c("tsig", "aic", "bic")) for (smp in c("common", "varying")) {
    fit <- bcs_fit(df, 2, model = 2, control = bcs_ctl(lag_criterion = crit, lag_sample = smp))
    ref <- ref_bcs(ref_bcs_units(df, 2), fit$settings, fit_tb(fit))
    if (bcs_near_tie(ref, fit$settings)) next
    expect_equal(fit$units$lags, ref$lags)
    expect_equal(fit$units$t_stat, ref$t, tolerance = 1e-10)
  }
  fit <- bcs_fit(df, 2, model = 1, factors = 2, control = bcs_ctl(pca_scale = TRUE))
  ref <- ref_bcs(ref_bcs_units(df, 2), fit$settings, fit_tb(fit))
  expect_equal(fit$units$t_stat, ref$t, tolerance = 1e-10)
  expect_equal(unname(fit$factors$f), ref$pca$f, tolerance = 1e-8)
  expect_equal(unname(fit$factors$loadings), unname(ref$pca$loadings), tolerance = 1e-8)
})

test_that("estimated heterogeneous breaks match the naive iteration (T2)", {
  for (cfg in list(list(model = 1, m = 1, K = 1), list(model = 2, m = 2, K = 1),
                   list(model = 4, m = 1, K = 2), list(model = 5, m = 1, K = 1, factors = 1))) {
    df <- sim_panel(N = 5, T = 36, K = cfg$K, seed = cfg$model + 10 * cfg$m)
    # With N = 5, BIC3 picks r = 4 and the iteration converges slowly; a
    # few iterations suffice to compare the paths.
    fit <- suppressWarnings(bcs_fit(df, cfg$K, model = cfg$model, breaks = cfg$m,
                                    factors = if (is.null(cfg$factors)) "bic3" else cfg$factors,
                                    control = bcs_ctl(iter_tol = 1e-8, iter_max = 8)))
    ref <- ref_hetero(ref_bcs_units(df, cfg$K), fit$settings, cfg$m, fit$settings$h,
                      r_fixed = cfg$factors, rmax = fit$settings$max_factors)
    expect_equal(fit_tb(fit), unname(ref$tb))
    expect_equal(fit$breaks$ssr, ref$ssr, tolerance = 1e-8)
    expect_equal(fit$breaks$r_iteration, ref$r)
    # The statistic is the known-break procedure at the estimated dates.
    bcs_compare(df, cfg$K, list(model = cfg$model, breaks = data.frame(
      id = rep(fit$breaks$dates$id, cfg$m), d = unlist(fit$breaks$dates[-1L]))))
  }
})

test_that("the homogeneous profile matches the oracle on every candidate (T2)", {
  df <- sim_panel(N = 5, T = 34, K = 1, seed = 77)
  for (model in c(1, 5)) {
    fit <- bcs_fit(df, 1, model = model, breaks = 1, break_type = "homogeneous",
                   factors = 1, lags = 1, control = bcs_ctl())
    units <- ref_bcs_units(df, 1)
    prof <- fit$breaks$profile
    tb <- match(prof$break_date, sort(unique(df$time)))
    tbar <- vapply(tb, function(b) ref_bcs(units, fit$settings, rep(list(b), 5))$tbar, numeric(1))
    expect_equal(prof$t_bar, tbar, tolerance = 1e-10)
    expect_equal(fit$statistics$value, min(tbar), tolerance = 1e-10)
    expect_equal(fit$statistics$statistic, paste0("Zinf_", c("c", "tau")[1 + (model == 5)]))
    expect_equal(fit$breaks$lambda[1, 1], tb[which.min(tbar)] / 34)
    expect_equal(min(prof$lambda) * 34, max(floor(0.15 * 34), if (model == 5) 3 else 2))
  }
})

test_that("Z standardises t-bar with the moments (Theorem 2)", {
  df <- sim_panel(N = 8, T = 40, K = 1, seed = 2)
  fit <- bcs_fit(df, 1, model = 2)
  z <- sqrt(8) * (fit$statistics$value - -1.535) / sqrt(0.341)
  expect_equal(fit$statistics$z, z, tolerance = 1e-14)
  expect_equal(fit$statistics$p_value, pnorm(z))
  fit3 <- bcs_fit(df, 1, model = 3, breaks = c(12, 28), control = bcs_ctl())
  expect_equal(fit3$moments$mean, -2.315)
  expect_equal(fit3$moments$var, 0.337)
  expect_equal(fit3$moments$source, "Table II")
  fit4 <- bcs_fit(df, 1, model = 6, breaks = 20, control = bcs_ctl(moments = "finite_sample"))
  expect_equal(fit4$moments$source, "WP Table 9")
  # T = 40 is below the table: the T = 50 row at lambda = 0.5.
  expect_equal(fit4$moments$mean, -2.022)
  expect_equal(fit4$moments$var, 0.383)
  fit6 <- bcs_fit(sim_panel(N = 8, T = 75, K = 1, seed = 2), 1, model = 2,
                  control = bcs_ctl(moments = "finite_sample"))
  # Linear in 1/T between T = 50 and T = 100.
  expect_equal(fit6$moments$mean, -1.549 + (-1.541 + 1.549) * (1 / 75 - 1 / 50) / (1 / 100 - 1 / 50),
               tolerance = 1e-12)
  fit5 <- bcs_fit(df, 1, model = 3, breaks = 17, control = bcs_ctl())
  expect_equal(fit5$moments$source, "simulated")
})

test_that("the statistics are invariant to scale, units, constants, trends and breaks (T4)", {
  df <- sim_panel(N = 6, T = 40, K = 2, seed = 21)
  same <- function(a, b, tol = 1e-10) {
    expect_equal(a$statistics$value, b$statistics$value, tolerance = tol)
    expect_equal(a$factors$r, b$factors$r)
    expect_equal(a$units$lags, b$units$lags)
  }
  for (args in list(list(model = 2, factors = 1), list(model = 6, breaks = 20, factors = 2),
                    list(model = 4, breaks = c(14, 26)), list(model = 5, breaks = 1, factors = 1),
                    list(model = 1, breaks = 1, break_type = "homogeneous", factors = 1))) {
    base <- do.call(bcs_fit, c(list(df, 2, control = bcs_ctl()), args))
    d2 <- df
    d2$y <- 3 * d2$y
    d2$x1 <- 0.5 * d2$x1
    same(base, do.call(bcs_fit, c(list(d2, 2, control = bcs_ctl()), args)))
    d3 <- df
    d3$y <- d3$y + 5 * d3$id
    same(base, do.call(bcs_fit, c(list(d3, 2, control = bcs_ctl()), args)))
    perm <- df[order(-df$id, df$time), ]
    p <- do.call(bcs_fit, c(list(perm, 2, control = bcs_ctl()), args))
    expect_equal(p$statistics$value, base$statistics$value, tolerance = 1e-10)
    expect_equal(p$units$t_stat[match(base$units$id, p$units$id)], base$units$t_stat, tolerance = 1e-10)
    if (args$model %in% c(2, 3, 5, 6)) {
      d4 <- df
      d4$y <- d4$y + 0.3 * d4$id * d4$time
      same(base, do.call(bcs_fit, c(list(d4, 2, control = bcs_ctl()), args)))
    }
    if (!is.null(args$breaks) && args$breaks[1] > 2) {
      d5 <- df
      d5$y <- d5$y + d5$id * (d5$time > args$breaks[1])
      same(base, do.call(bcs_fit, c(list(d5, 2, control = bcs_ctl()), args)))
    }
    if (args$model %in% c(3, 6)) {
      d6 <- df
      d6$y <- d6$y + 0.2 * d6$id * pmax(d6$time - args$breaks[1], 0)
      same(base, do.call(bcs_fit, c(list(d6, 2, control = bcs_ctl()), args)))
    }
    two <- do.call(bcs_fit, c(list(df, 2, control = bcs_ctl(), n_threads = 2), args))
    expect_identical(two$statistics, base$statistics)
    expect_identical(two$breaks$dates, base$breaks$dates)
  }
})

test_that("results are deterministic and leave the RNG alone (T5)", {
  df <- sim_panel(N = 6, T = 40, K = 1, seed = 5)
  set.seed(99)
  before <- .Random.seed
  a <- bcs_fit(df, 1, model = 6, breaks = 17, control = bcs_ctl())
  expect_identical(.Random.seed, before)
  rm(list = ls(bcs_null_cache), envir = bcs_null_cache)
  b <- bcs_fit(df, 1, model = 6, breaks = 17, control = bcs_ctl(), n_threads = 3)
  expect_identical(a$statistics, b$statistics)
  expect_identical(a$units, b$units)
  expect_identical(a$moments, b$moments)
})

test_that("break dates are time values and the last period of the old regime", {
  df <- sim_panel(N = 4, T = 30, K = 1, seed = 3)
  df$time <- df$time + 1990
  fit <- bcs_fit(df, 1, model = 4, breaks = 2005)
  expect_equal(fit$breaks$positions[1, 1], 15L)
  expect_equal(fit$units$break_date, rep(2005, 4))
  expect_equal(fit$units$lambda, rep(0.5, 4))
  dd <- df
  dd$time <- as.Date("2000-01-01") + df$time - 1991
  fd <- bcs_fit(dd, 1, model = 4, breaks = as.Date("2000-01-15"))
  expect_equal(fd$statistics$value, fit$statistics$value, tolerance = 1e-12)
})

test_that("invalid inputs give informative errors (T6)", {
  df <- sim_panel(N = 4, T = 30, K = 1, seed = 1)
  expect_error(bcs_fit(df, 1), "`model` is required")
  expect_error(bcs_fit(df, 1, model = 7), "between 1 and 6")
  expect_error(bcs_fit(df, 1, model = 4), "Model 4 needs breaks")
  expect_error(bcs_fit(df, 1, model = 3, breaks = 1), "known common breaks only")
  expect_error(bcs_fit(df, 1, model = 6, breaks = list(`1` = 10, `2` = 10, `3` = 10, `4` = 10)),
               "common to all units")
  expect_error(bcs_fit(df, 1, model = 2, breaks = 99), "1 or 2 breaks")
  expect_error(bcs_fit(df, 1, model = 2, breaks = c(5, 99)), "not among the time values")
  # A single 1 or 2 is a number of breaks; 3 is a time value of this panel.
  expect_equal(bcs_fit(df, 1, model = 2, breaks = 3)$breaks$type, "known")
  expect_error(bcs_fit(df, 1, model = 2, breaks = 30), "no period of the new regime")
  expect_error(bcs_fit(df, 1, model = 2, breaks = 15, break_type = "homogeneous"),
               "applies only to estimated breaks")
  expect_error(bcs_fit(df, 1, model = 2, max_factors = 4), "must be below min\\(N, n\\)")
  expect_error(bcs_fit(df, 1, model = 2, factors = "bic9"), "must be a number or one of")
  expect_error(bcs_fit(df, 1, model = 4, breaks = 29), "rank deficient")
  expect_error(bcs_fit(df, 1, model = 4, breaks = list(`1` = 10)), "missing for units 2, 3, 4")
  ub <- sim_panel(N = 4, T = 30, K = 1, balanced = FALSE, seed = 1)
  expect_error(bcs_fit(ub, 1, model = 1), "balanced panel")
  expect_error(bcs_fit(df, 1, model = 1, control = list()), "banerjee_control")
  expect_error(bcs_fit(df, 1, model = 1, lags = 14), "Too few periods for ADF")
})
