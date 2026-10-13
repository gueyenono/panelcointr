test_that("print, summary, plot and the tidiers run", {
  ctl <- banerjee_control(sim_reps = 1000L, sim_steps = 100L)
  df <- sim_panel(N = 6, T = 40, K = 1, seed = 4)
  fits <- list(
    pcoint_banerjee(y ~ x1, data = df, model = 2),
    pcoint_banerjee(y ~ x1, data = df, model = 6, breaks = c(12, 28), factors = 2, control = ctl),
    pcoint_banerjee(y ~ x1, data = df, model = 1, breaks = 1, break_type = "homogeneous",
                    factors = 1),
    suppressWarnings(pcoint_banerjee(y ~ x1, data = df, model = 4, breaks = 1, factors = 1,
                                     control = banerjee_control(iter_max = 5L)))
  )
  pdf(NULL)
  on.exit(grDevices::dev.off())
  for (fit in fits) {
    expect_output(print(fit), "Banerjee and Carrion-i-Silvestre")
    expect_output(print(summary(fit)), "Units")
    expect_s3_class(as.data.frame(fit), "data.frame")
    expect_identical(plot(fit), fit)
  }
  expect_output(print(fits[[3]]), "Table III")
  expect_output(print(fits[[2]]), "MQ_c")
  if (requireNamespace("generics", quietly = TRUE)) {
    expect_equal(nrow(generics::tidy(fits[[1]])), 1L)
    expect_equal(generics::glance(fits[[2]])$n_factors, 2L)
  }
})
