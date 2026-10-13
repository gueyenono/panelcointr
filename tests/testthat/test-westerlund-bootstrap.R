# T5: bootstrap mechanics.

boot_state <- function(df, K, fit, reselect = TRUE, burn = 7L) {
  panel <- prepare_panel(xfml(K), df, c("id", "time"))
  N <- length(panel$T)
  s <- fit$settings
  ones <- rep(1L, N)
  state <- westerlund_boot_prepare_rs(
    panel$y, as.vector(panel$X), K, panel$offsets, as.character(panel$ids), panel$start,
    fit$units$lags, fit$units$leads, s$lags[1] * ones, s$lags[2] * ones,
    s$leads[1] * ones, s$leads[2] * ones, fit$units$lrwindow,
    westerlund_engine_settings(s), reselect, burn
  )
  list(state = state, panel = panel)
}

test_that("first stage and centred residuals match R (unbalanced panel)", {
  df <- sim_panel(N = 4, T = 40, K = 2, balanced = FALSE, seed = 21)
  for (det in c("constant", "trend")) {
    fit <- pcoint_westerlund(xfml(2), df, index = c("id", "time"), deterministic = det,
                             lags = c(0, 2), leads = c(0, 2))
    b <- boot_state(df, 2, fit)
    st <- westerlund_boot_state_rs(b$state)
    ref <- ref_boot_w(ref_units(df), b$panel$start, fit$units$lags, fit$units$leads,
                      constant = det == "trend")
    expect_equal(st$n_periods, nrow(ref$w))
    expect_equal(matrix(st$w, st$n_periods), unname(ref$w), tolerance = 1e-10)
    for (i in 1:4) {
      expect_equal(st$units[[i]]$a, ref$first_stage[[i]]$a, tolerance = 1e-10)
      expect_equal(st$units[[i]]$g, ref$first_stage[[i]]$g, tolerance = 1e-10)
    }
  }
})

test_that("bootstrap panels follow the null DGP on each unit's calendar window", {
  df <- sim_panel(N = 4, T = 40, K = 2, balanced = FALSE, seed = 22)
  fit <- pcoint_westerlund(xfml(2), df, index = c("id", "time"), deterministic = "trend",
                           lags = c(1, 2), leads = c(1, 2))
  b <- boot_state(df, 2, fit, burn = 7L)
  st <- westerlund_boot_state_rs(b$state)
  info <- westerlund_boot_info_rs(b$state)
  expect_equal(info$t_span, max(b$panel$start + b$panel$T))
  expect_equal(info$len, 7 + info$t_span + max(fit$units$leads))
  set.seed(5)
  idx <- sample.int(st$n_periods, info$len, replace = TRUE)
  gen <- westerlund_boot_panel_rs(b$state, idx)
  ref <- ref_boot_panel(st, idx, 2)
  expect_equal(gen$y, ref$y, tolerance = 1e-12)
  expect_equal(gen$x, ref$x, tolerance = 1e-12)
  expect_error(westerlund_boot_panel_rs(b$state, idx[-1]), "indices per replicate")
  expect_error(westerlund_boot_panel_rs(b$state, replace(idx, 3, st$n_periods + 1L)),
               "outside")
})

test_that("results are reproducible and do not depend on threads or chunks", {
  run <- function(...) {
    set.seed(42)
    hce_fit(lags = c(1, 2), leads = c(0, 1), bootstrap = 60, ...)
  }
  a <- run()
  b <- run()
  c2 <- run(n_threads = 2)
  d <- run(control = westerlund_control(boot_chunk = 7))
  expect_identical(a$boot$stats, b$boot$stats)
  expect_identical(a$boot$stats, c2$boot$stats)
  expect_identical(a$boot$stats, d$boot$stats)
  expect_identical(a$statistics, d$statistics)
  expect_false(anyNA(a$boot$stats))
  expect_identical(a$boot$rng$seed, {
    set.seed(42)
    .Random.seed
  })
})

test_that("bootstrap p-values and critical values follow their definitions", {
  set.seed(3)
  fit <- hce_fit(lags = 1, leads = 1, bootstrap = 49)
  S <- fit$statistics$value
  below <- colSums(sweep(fit$boot$stats, 2, S, `<=`))
  expect_equal(fit$statistics$p_boot, unname((1 + below) / 50))
  expect_equal(fit$boot$p_value, fit$statistics$p_boot)
  q <- apply(fit$boot$stats, 2, quantile, probs = c(0.01, 0.05, 0.1), type = 1, names = FALSE)
  expect_equal(unname(fit$boot$critical), unname(q))
  expect_identical(dimnames(fit$boot$critical), list(c("1%", "5%", "10%"), c("Gt", "Ga", "Pt", "Pa")))

  set.seed(3)
  fx <- hce_fit(lags = 1, leads = 1, bootstrap = 49,
                control = westerlund_control(boot_pvalue = "xtwest", keep_boot = FALSE))
  expect_equal(fx$statistics$p_boot, unname(below / 49))
  expect_null(fx$boot$stats)
})

test_that("boot_reselect = FALSE keeps the selected orders", {
  set.seed(9)
  a <- hce_fit(lags = c(1, 3), leads = c(0, 2), bootstrap = 30)
  set.seed(9)
  b <- hce_fit(lags = c(1, 3), leads = c(0, 2), bootstrap = 30,
               control = westerlund_control(boot_reselect = FALSE))
  expect_false(isTRUE(all.equal(a$boot$stats, b$boot$stats)))
  # With fixed orders, re-selection is a no-op.
  set.seed(9)
  c1 <- hce_fit(lags = 2, leads = 1, bootstrap = 30)
  set.seed(9)
  c2 <- hce_fit(lags = 2, leads = 1, bootstrap = 30,
                control = westerlund_control(boot_reselect = FALSE))
  expect_identical(c1$boot$stats, c2$boot$stats)
})

test_that("the bootstrap needs enough common periods", {
  mk <- function(start2) {
    d1 <- sim_panel(N = 1, T = 40, seed = 1)
    d2 <- sim_panel(N = 1, T = 40, seed = 2)
    d2$id <- 2
    d2$time <- d2$time + start2 - 1
    rbind(d1, d2)
  }
  expect_error(pcoint_westerlund(y ~ x1, mk(33), index = c("id", "time"), lags = 1,
                                 bootstrap = 9), "calendar periods")
  expect_warning(pcoint_westerlund(y ~ x1, mk(25), index = c("id", "time"), lags = 1,
                                   bootstrap = 9), "common to all units")
  expect_error(pcoint_westerlund(y ~ x1, mk(25), index = c("id", "time"), lags = 1,
                                 bootstrap = 9,
                                 control = westerlund_control(min_boot_periods = 30)),
               "at least 30")
})

test_that("the bootstrap works beyond the tabulated moments", {
  df <- sim_panel(N = 3, T = 60, K = 7, seed = 4)
  set.seed(1)
  expect_warning(fit <- pcoint_westerlund(xfml(7), df, index = c("id", "time"), lags = 0,
                                          bootstrap = 19), "No asymptotic moments")
  expect_true(all(is.na(fit$statistics$p_value)))
  expect_true(all(fit$statistics$p_boot > 0 & fit$statistics$p_boot <= 1))
})
