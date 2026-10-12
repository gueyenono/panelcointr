test_that("print() shows the xtwest-style table", {
  fit <- hce_fit(lags = c(1, 3), leads = c(0, 3))
  out <- capture.output(print(fit))
  expect_true(any(grepl("Westerlund \\(2007\\)", out)))
  expect_true(any(grepl("Lags: 1 to 3 by AIC \\(mean 2.05\\)", out)))
  expect_true(any(grepl("^ +Gt +-2.714", out)))
  expect_false(any(grepl("Robust P", out)))
  expect_identical(print(fit), fit) |> capture.output() |> invisible()
})

test_that("summary() reports the mean-group model", {
  fit <- hce_fit(lags = 1, leads = 1)
  s <- summary(fit)
  expect_s3_class(s, "summary.pcoint_westerlund")
  b <- fit$coefficients
  lr <- -b[, "L.loggdp"] / b[, "L.loghex"]
  expect_equal(s$mean_group$long_run["loggdp", "estimate"], mean(lr))
  expect_equal(s$mean_group$long_run["loggdp", "std_error"], sd(lr) / sqrt(20))
  expect_equal(s$mean_group$short_run["ec", "estimate"], mean(fit$units$alpha))
  expect_true(all(c("D.loggdp", "F1D.loggdp", "L1D.loghex") %in% rownames(s$mean_group$short_run)))
  out <- capture.output(print(s))
  expect_true(any(grepl("Mean-group error-correction model", out)))
  # With selected orders, only terms present for every unit are averaged.
  s2 <- summary(hce_fit(lags = c(1, 3), leads = c(0, 3)))
  expect_false("L3D.loghex" %in% rownames(s2$mean_group$short_run))
  expect_true("L1D.loghex" %in% rownames(s2$mean_group$short_run))
})

test_that("bootstrap output is printed and plotted", {
  set.seed(1)
  fit <- hce_fit(lags = 1, leads = 1, bootstrap = 29)
  out <- capture.output(print(fit))
  expect_true(any(grepl("Robust P", out)))
  expect_true(any(grepl("Bootstrap: 29 replications", out)))
  pdf(NULL)
  on.exit(dev.off())
  expect_identical(plot(fit), fit)
  expect_error(plot(hce_fit(lags = 1)), "No bootstrap distribution")
})

test_that("as.data.frame(), tidy() and glance() return the statistics", {
  fit <- hce_fit(lags = 1, leads = 1)
  expect_identical(as.data.frame(fit), fit$statistics)
  skip_if_not_installed("generics")
  td <- generics::tidy(fit)
  expect_named(td, c("test", "statistic", "z", "p.value", "p.value.boot", "alternative"))
  expect_equal(td$statistic, fit$statistics$value)
  gl <- generics::glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$n_units, 20L)
  expect_equal(gl$n_boot, 0L)
})
