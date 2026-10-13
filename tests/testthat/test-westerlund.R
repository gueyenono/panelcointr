# T1: published xtwest output (Persyn and Westerlund 2008, Stata Journal 8(2))
# on xtwestdata (`hce`), loghex ~ loggdp with constant and trend. Absolute
# tolerance 1e-3: the printed precision (Stata also stores temporaries as
# floats).
sj_runs <- list(
  run1 = list(args = list(lags = 1, leads = 1, lrwindow = 3, compat = "xtwest"),
              stats = c(-2.681, -10.927, -12.035, -10.524), lags = 1, leads = 1),
  run2 = list(args = list(lags = c(1, 2), leads = c(0, 1), lrwindow = 2, compat = "xtwest"),
              stats = c(-2.736, -11.254, -12.859, -11.773), lags = 1.15, leads = 0.05),
  run3 = list(args = list(lags = c(1, 3), leads = c(0, 3), lrwindow = 3,
                          compat = "xtwest_westerlund"),
              stats = c(-4.082, -27.702, -12.969, -22.470),
              z = c(-9.613, -10.626, -4.100, -10.119), lags = 2.8, leads = 1.65)
)

test_that("compat modes reproduce the published xtwest output (T1)", {
  for (run in sj_runs) {
    fit <- do.call(hce_fit, run$args)
    expect_lt(max(abs(fit$statistics$value - run$stats)), 1e-3)
    expect_equal(mean(fit$units$lags), run$lags, tolerance = 1e-12)
    expect_equal(mean(fit$units$leads), run$leads, tolerance = 1e-12)
    if (!is.null(run$z)) expect_lt(max(abs(fit$statistics$z - run$z)), 1e-3)
  }
})

test_that("z-values use Table 1 (Theorem 1 standardisation)", {
  fit <- do.call(hce_fit, sj_runs$run1$args)
  m <- westerlund_table1
  N <- 20
  st <- fit$statistics$value
  z <- c(sqrt(N) * (st[1] - m$mean$Gt[3, 1]) / sqrt(m$var$Gt[3, 1]),
         sqrt(N) * (st[2] - m$mean$Ga[3, 1]) / sqrt(m$var$Ga[3, 1]),
         (st[3] - sqrt(N) * m$mean$Pt[3, 1]) / sqrt(m$var$Pt[3, 1]),
         sqrt(N) * (st[4] - m$mean$Pa[3, 1]) / sqrt(m$var$Pa[3, 1]))
  expect_equal(fit$statistics$z, z, tolerance = 1e-14)
  expect_equal(fit$statistics$p_value, pnorm(z), tolerance = 1e-14)
  # Published statistics of run 1 give z-values -1.731, 0.713, -2.959, -1.160.
  expect_lt(max(abs(fit$statistics$z - c(-1.731, 0.713, -2.959, -1.160))), 1e-3)
})

test_that("the oracle reproduces the published xtwest output (T1 on the oracle)", {
  units <- lapply(split(hce, hce$ctr), function(d) list(y = d$loghex, X = as.matrix(d["loggdp"])))
  for (run in sj_runs) {
    s <- do.call(hce_fit, run$args)$settings
    ref <- ref_westerlund(units, s, rep(run$args$lrwindow, 20))
    expect_lt(max(abs(unname(ref$stats) - run$stats)), 1e-3)
    expect_equal(mean(ref$lags), run$lags)
    expect_equal(mean(ref$leads), run$leads)
  }
})

# T2: the Rust engine against the pure-R oracle on random panels.
compare_with_oracle <- function(df, K, args, tol = 1e-10) {
  fit <- do.call(pcoint_westerlund, c(list(xfml(K), data = df, index = c("id", "time")), args))
  ref <- ref_westerlund(ref_units(df), fit$settings, fit$units$lrwindow)
  if (ref_near_tie(ref)) return(FALSE)
  expect_equal(fit$units$lags, ref$lags)
  expect_equal(fit$units$leads, ref$leads)
  expect_equal(fit$statistics$value, unname(ref$stats), tolerance = tol)
  expect_equal(fit$units$alpha, ref$alpha, tolerance = tol)
  expect_equal(fit$units$se_alpha, ref$se, tolerance = tol)
  expect_equal(fit$units$alpha1, ref$alpha1, tolerance = tol)
  TRUE
}

test_that("Rust statistics match the pure-R oracle (T2)", {
  seed <- 100
  for (det in c("none", "constant", "trend")) {
    for (K in 1:3) {
      for (balanced in c(TRUE, FALSE)) {
        for (compat in c("none", "xtwest", "xtwest_westerlund")) {
          if (compat == "xtwest_westerlund" && (det == "none" || K > 1)) next
          for (ic in c(FALSE, TRUE)) {
            args <- list(deterministic = det, compat = compat,
                         lags = if (ic) c(0, 2) else 1, leads = if (ic) c(0, 1) else 1)
            repeat {
              seed <- seed + 1
              df <- sim_panel(N = 5, T = 45, K = K, balanced = balanced, seed = seed,
                              drift = if (seed %% 2) 0.2 else 0)
              if (compare_with_oracle(df, K, args)) break
            }
          }
        }
      }
    }
  }
})

test_that("every control setting matches the oracle (T2)", {
  df <- sim_panel(N = 6, T = 50, K = 2, balanced = FALSE, seed = 7, drift = 0.1)
  base <- list(deterministic = "trend", lags = c(0, 2), leads = c(0, 2))
  variants <- list(
    list(criterion = "bic"),
    list(control = westerlund_control(t_norm = "df")),
    list(control = westerlund_control(t_norm = "T")),
    list(control = westerlund_control(panel_lags = "mean")),
    list(control = westerlund_control(ic_sample = "varying")),
    list(control = westerlund_control(lrv_deterministics = "implied")),
    list(control = westerlund_control(lrv_deterministics = "none", demean_u = FALSE)),
    list(lrwindow = 0),
    list(lrwindow = 6),
    list(compat = "xtwest", control = westerlund_control(ic_sample = "common", t_norm = "T"))
  )
  for (v in variants) {
    expect_true(compare_with_oracle(df, 2, utils::modifyList(base, v)))
  }
})

test_that("Gram-Cholesky selection matches brute-force QR (T3)", {
  for (seed in 1:6) {
    df <- sim_panel(N = 1, T = 70, K = 2, seed = seed)
    for (crit in c("aic", "bic", "westerlund")) {
      s <- list(deterministic = c("none", "constant", "trend")[seed %% 3 + 1], criterion = crit,
                ic_sample = "common", ic_ties = "smaller", lrv_deterministics = "regression",
                lrv_y_sample = "all", u_sample = "estimation", demean_u = TRUE,
                t_norm = "effective", panel_lags = "individual", se_alpha = "df",
                sigma_df = "own", ga_weights = "units", tbar = "units")
      X <- as.vector(as.matrix(df[c("x1", "x2")]))
      a <- ecm_ic_rs(df$y, X, 2L, 0L, 4L, 0L, 3L, s, "gram")
      b <- ecm_ic_rs(df$y, X, 2L, 0L, 4L, 0L, 3L, s, "qr")
      expect_equal(a$ic, b$ic, tolerance = 1e-9)
      expect_identical(c(a$lags, a$leads), c(b$lags, b$leads))
    }
  }
})

# T4: invariances.
invariance_fits <- function(df, K, args, transform) {
  f <- function(d) {
    do.call(pcoint_westerlund, c(list(xfml(K), data = d, index = c("id", "time")), args))
  }
  list(f(df), f(transform(df)))
}

# `rows`: row of the second fit's unit table holding each unit of the first.
expect_same_fit <- function(pair, tol = 1e-10, rows = seq_len(nrow(pair[[1]]$units))) {
  expect_equal(pair[[2]]$statistics$value, pair[[1]]$statistics$value, tolerance = tol)
  expect_identical(pair[[2]]$units$lags[rows], pair[[1]]$units$lags)
  expect_identical(pair[[2]]$units$leads[rows], pair[[1]]$units$leads)
}

test_that("statistics are invariant to rescaling, unit order and deterministics (T4)", {
  df <- sim_panel(N = 6, T = 50, K = 2, balanced = FALSE, seed = 11, drift = 0.15)
  for (lag_args in list(list(lags = 1, leads = 1), list(lags = c(0, 2), leads = c(0, 2)))) {
    for (det in c("constant", "trend")) {
      args <- c(lag_args, deterministic = det)
      # Rescaling y and each x separately.
      expect_same_fit(invariance_fits(df, 2, args, function(d) {
        d$y <- 7.5 * d$y
        d$x1 <- 0.01 * d$x1
        d$x2 <- -3 * d$x2
        d
      }))
      # Permuting the units (relabelling the ids) and shuffling the rows.
      perm <- c(4, 6, 1, 3, 5, 2)
      expect_same_fit(invariance_fits(df, 2, args, function(d) {
        d$id <- perm[d$id]
        d[sample.int(nrow(d)), ]
      }), rows = perm)
      # Adding constants to y and x.
      expect_same_fit(invariance_fits(df, 2, args, function(d) {
        d$y <- d$y + 100
        d$x1 <- d$x1 - 50
        d
      }))
      # Recoding time (t versus calendar years).
      expect_same_fit(invariance_fits(df, 2, args, function(d) {
        d$time <- d$time + 1969L
        d
      }))
    }
    # Linear trends in y and x in Case 3 (needs demeaning of u, D19).
    expect_same_fit(invariance_fits(df, 2, c(lag_args, deterministic = "trend"), function(d) {
      d$y <- d$y + 0.8 * d$time
      d$x1 <- d$x1 - 0.3 * d$time
      d$x2 <- d$x2 + 2 * d$time
      d
    }))
  }
})

test_that("varying-sample selection depends on the units of y (D7 regression test)", {
  # Persyn and Westerlund (2008) run 2 with loghex multiplied by 3: xtwest's
  # mean lags and leads move from 1.15/0.05 to 1.20/0.20.
  hce3 <- transform(hce, loghex = 3 * loghex)
  args <- list(lags = c(1, 2), leads = c(0, 1), lrwindow = 2, compat = "xtwest")
  a <- do.call(hce_fit, args)
  b <- do.call(pcoint_westerlund, c(list(loghex ~ loggdp, data = hce3, index = c("ctr", "year"),
                                         deterministic = "trend"), args))
  expect_equal(c(mean(b$units$lags), mean(b$units$leads)), c(1.20, 0.20))
  expect_lt(abs(b$statistics$value[1] - (-2.843)), 1e-3)
  expect_false(isTRUE(all.equal(a$statistics$value, b$statistics$value)))
  # On the common sample (the default) the selection is scale invariant.
  args$control <- westerlund_control(ic_sample = "common")
  a <- do.call(hce_fit, args)
  b <- do.call(pcoint_westerlund, c(list(loghex ~ loggdp, data = hce3, index = c("ctr", "year"),
                                         deterministic = "trend"), args))
  expect_equal(b$statistics$value, a$statistics$value, tolerance = 1e-10)
})

test_that("the default procedure differs from xtwest as documented", {
  a <- hce_fit(lags = 1, leads = 1, lrwindow = 3)
  b <- hce_fit(lags = 1, leads = 1, lrwindow = 3, compat = "xtwest")
  expect_identical(a$settings$t_norm, "effective")
  expect_identical(b$settings$t_norm, "df")
  # Gt only differs through nothing at fixed lags in a balanced panel.
  expect_equal(a$statistics$value[1], b$statistics$value[1], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(a$statistics$value[2], b$statistics$value[2])))
})

test_that("results are returned in the documented structure", {
  fit <- hce_fit(lags = c(1, 3), leads = c(0, 3))
  expect_s3_class(fit, c("pcoint_westerlund", "pcoint_test"))
  expect_named(fit$statistics, c("statistic", "value", "z", "p_value", "p_boot", "alternative"))
  expect_identical(fit$statistics$statistic, c("Gt", "Ga", "Pt", "Pa"))
  expect_true(all(is.na(fit$statistics$p_boot)))
  expect_named(fit$units, c("id", "T", "n", "k", "lags", "leads", "lrwindow", "alpha",
                            "se_alpha", "t_alpha", "alpha1", "sigma"))
  expect_equal(nrow(fit$units), 20)
  expect_equal(fit$units$n, fit$units$T - fit$units$lags - fit$units$leads - 1L)
  expect_equal(fit$units$lrwindow, rep(3L, 20)) # floor(4 * 0.32^(2/9)) = 3
  expect_equal(fit$units$t_alpha, fit$units$alpha / fit$units$se_alpha)
  expect_equal(dim(fit$coefficients), c(20L, 2 + 1 + 3 + 7 + 1))
  expect_equal(unname(fit$coefficients[, "L.loghex"]), fit$units$alpha, tolerance = 1e-12)
  expect_equal(sum(is.na(fit$coefficients[, "L3D.loghex"])), sum(fit$units$lags < 3))
  expect_equal(fit$N, 20L)
  expect_equal(fit$K, 1L)
  expect_equal(fit$case, 3L)
  expect_identical(fit$settings$criterion, "aic")
  # The panel statistics use the individual orders by default.
  expect_true(is.na(fit$panel$lags))
  fx <- hce_fit(lags = c(1, 3), leads = c(0, 3), compat = "xtwest")
  expect_equal(fx$panel$lags, floor(mean(fx$units$lags)))
})

test_that("coefficients agree with lm() for a unit", {
  fit <- hce_fit(lags = 2, leads = 1)
  d <- hce[hce$ctr == 5, ]
  des <- ref_design(d$loghex, as.matrix(d["loggdp"]), 2, 1, 2)
  rows <- 4:31
  b <- stats::lm.fit(des$Z[rows, ], des$dy[rows])$coefficients
  expect_equal(unname(fit$coefficients["5", ]), unname(b[c(1, 2, 3, 4, 5, 6:9, 10)]),
               tolerance = 1e-10)
})

test_that("more than six regressors give NA asymptotics with a warning", {
  df <- sim_panel(N = 4, T = 60, K = 7, seed = 3)
  expect_warning(fit <- pcoint_westerlund(xfml(7), df, index = c("id", "time"), lags = 0),
                 "No asymptotic moments")
  expect_true(all(is.na(fit$statistics$z)))
  expect_true(all(is.na(fit$statistics$p_value)))
  expect_true(all(is.finite(fit$statistics$value)))
})

test_that("argument errors are informative", {
  expect_error(pcoint_westerlund(loghex ~ loggdp, hce, index = c("ctr", "year")),
               "`lags` is required")
  expect_error(hce_fit(lags = -1), "non-negative whole number")
  expect_error(hce_fit(lags = c(1, 2, 3)), "range")
  expect_error(hce_fit(lags = 1, lrwindow = 1.5), "lrwindow")
  expect_error(hce_fit(lags = 1, bootstrap = -1), "bootstrap")
  expect_error(hce_fit(lags = 1, n_threads = 0), "n_threads")
  expect_error(hce_fit(lags = 1, control = list()), "westerlund_control")
  expect_error(pcoint_westerlund(loghex ~ loggdp, hce, index = c("ctr", "year"), lags = 1,
                                 compat = "xtwest_westerlund", deterministic = "none"),
               "needs a constant")
  expect_warning(hce_fit(lags = 1, compat = "xtwest_westerlund", criterion = "bic"),
                 "criterion")
  df <- sim_panel(N = 3, T = 30, K = 2, seed = 1)
  expect_error(pcoint_westerlund(xfml(2), df, index = c("id", "time"), lags = 1,
                                 compat = "xtwest_westerlund"), "single regressor")
})

test_that("reversed ranges are sorted, as in xtwest", {
  a <- hce_fit(lags = c(3, 1), leads = c(2, 0))
  b <- hce_fit(lags = c(1, 3), leads = c(0, 2))
  expect_identical(a$statistics, b$statistics)
})

test_that("rank deficiency is reported with the unit and the model", {
  df <- sim_panel(N = 4, T = 40, K = 2, seed = 5)
  df$x2[df$id == 3] <- 2 * df$x1[df$id == 3]
  expect_error(pcoint_westerlund(xfml(2), df, index = c("id", "time"), lags = 1),
               "unit 3: .*rank deficient")
})
