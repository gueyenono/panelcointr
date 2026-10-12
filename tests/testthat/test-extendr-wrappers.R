# T3: Rust kernels against base R.

test_that("QR least squares matches lm.fit", {
  set.seed(1)
  n <- 50
  Z <- cbind(1, seq_len(n) * 1000, matrix(rnorm(n * 4), n))
  Z[, 6] <- Z[, 5] + 1e-4 * rnorm(n) # near-collinear pair
  y <- rnorm(n)
  r <- ols_rs(as.vector(Z), ncol(Z), y)
  f <- lm.fit(Z, y)
  expect_equal(r$coef, unname(f$coefficients), tolerance = 1e-8)
  expect_equal(r$rss, sum(f$residuals^2), tolerance = 1e-12)
  # Frisch-Waugh-Lovell through the last column.
  yres <- lm.fit(Z[, -6], Z[, 6])$residuals
  dyres <- lm.fit(Z[, -6], y)$residuals
  expect_equal(r$r_last^2, sum(yres^2), tolerance = 1e-10)
  expect_equal(r$r_last * r$qty_last, sum(yres * dyres), tolerance = 1e-10)

  Zw <- cbind(1, matrix(rnorm(n * 3), n))
  expect_equal(ols_rs(as.vector(Zw), 4L, y)$coef, unname(lm.fit(Zw, y)$coefficients),
               tolerance = 1e-12)
  expect_error(ols_rs(as.vector(cbind(Zw, Zw[, 2] - Zw[, 3])), 5L, y), "rank deficient at column 5")
})

test_that("the Bartlett long-run variance matches its definition", {
  set.seed(2)
  v <- rnorm(37)
  for (M in c(0, 1, 3, 10, 40)) {
    expect_equal(lrv_rs(v, M), ref_lrv(v, M), tolerance = 1e-12)
  }
  expect_equal(lrv_rs(v, 0), mean(v^2))
})
