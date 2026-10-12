# T6: input validation and panel preparation.

test_that("panels are sorted and stacked", {
  df <- sim_panel(N = 3, T = 10, K = 2, balanced = FALSE, seed = 1)
  shuffled <- df[sample.int(nrow(df)), c(4, 3, 1, 5, 2)]
  p <- prepare_panel(xfml(2), shuffled, index = c("id", "time"))
  expect_equal(p$ids, 1:3)
  expect_equal(p$y, df$y)
  expect_equal(unname(p$X), unname(as.matrix(df[c("x1", "x2")])))
  expect_equal(p$offsets, c(0L, cumsum(as.integer(table(df$id)))))
  starts <- tapply(df$time, df$id, min)
  expect_equal(p$start, as.integer(starts - min(starts)))
  expect_identical(p$x_names, c("x1", "x2"))
  expect_identical(p$y_name, "y")
})

test_that("the first two columns are the default index", {
  df <- sim_panel(N = 3, T = 30, seed = 2)
  a <- pcoint_westerlund(y ~ x1, df, lags = 1)
  b <- pcoint_westerlund(y ~ x1, df, index = c("id", "time"), lags = 1)
  expect_identical(a$statistics, b$statistics)
})

test_that("tibbles, data.tables and pdata.frames are accepted", {
  df <- sim_panel(N = 3, T = 30, seed = 3)
  ref <- pcoint_westerlund(y ~ x1, df, index = c("id", "time"), lags = 1)$statistics
  skip_if_not_installed("tibble")
  expect_identical(pcoint_westerlund(y ~ x1, tibble::as_tibble(df), index = c("id", "time"),
                                     lags = 1)$statistics, ref)
  skip_if_not_installed("data.table")
  expect_identical(pcoint_westerlund(y ~ x1, data.table::as.data.table(df),
                                     index = c("id", "time"), lags = 1)$statistics, ref)
  skip_if_not_installed("plm")
  pdf <- plm::pdata.frame(df[c("time", "x1", "y", "id")], index = c("id", "time"))
  expect_equal(pcoint_westerlund(y ~ x1, pdf, lags = 1)$statistics, ref)
})

test_that("dates, factors and character time are mapped to ranks", {
  df <- sim_panel(N = 3, T = 30, seed = 4)
  ref <- pcoint_westerlund(y ~ x1, df, index = c("id", "time"), lags = 1)$statistics
  dated <- transform(df, time = as.Date("2000-01-01") + 30 * (time - 1))
  expect_identical(pcoint_westerlund(y ~ x1, dated, index = c("id", "time"), lags = 1)$statistics, ref)
  fct <- transform(df, time = factor(sprintf("t%02d", time)))
  expect_identical(pcoint_westerlund(y ~ x1, fct, index = c("id", "time"), lags = 1)$statistics, ref)
  chr <- transform(df, time = sprintf("t%02d", time), id = paste0("unit", id))
  fit <- pcoint_westerlund(y ~ x1, chr, index = c("id", "time"), lags = 1)
  expect_identical(fit$statistics, ref)
  expect_identical(fit$units$id, paste0("unit", 1:3))
})

test_that("leading and trailing missing values are dropped", {
  df <- sim_panel(N = 3, T = 30, seed = 5)
  df$y[df$id == 2 & df$time %in% 1:2] <- NA
  df$x1[df$id == 3 & df$time == 30] <- NA
  p <- prepare_panel(y ~ x1, df, c("id", "time"))
  expect_equal(p$T, c(30L, 28L, 29L))
  expect_equal(p$start, c(0L, 2L, 0L))
  fit <- pcoint_westerlund(y ~ x1, df, index = c("id", "time"), lags = 1)
  expect_equal(fit$units$T, c(30L, 28L, 29L))
})

test_that("interior missing values and gaps are errors naming the units", {
  df <- sim_panel(N = 4, T = 30, seed = 6)
  bad <- df
  bad$y[bad$id == 2 & bad$time == 10] <- NA
  bad$x1[bad$id == 4 & bad$time == 20] <- NA
  expect_error(prepare_panel(y ~ x1, bad, c("id", "time")),
               "Missing values inside the sample of units 2, 4")
  gap <- df[!(df$id == 3 & df$time == 15), ]
  expect_error(prepare_panel(y ~ x1, gap, c("id", "time")), "Gaps .* within unit 3")
  dup <- rbind(df, df[df$id == 1 & df$time == 5, ])
  expect_error(prepare_panel(y ~ x1, dup, c("id", "time")), "Duplicated time periods within unit: 1")
})

test_that("units too short for the requested orders are listed", {
  df <- sim_panel(N = 4, T = 30, seed = 7)
  df <- df[!(df$id %in% c(2, 4) & df$time > 12), ]
  expect_error(pcoint_westerlund(y ~ x1, df, index = c("id", "time"), lags = c(1, 3), leads = 1,
                                 deterministic = "trend"),
               "needs at least 18 periods, but 2 units have fewer: 2 \\(T = 12\\), 4 \\(T = 12\\)")
})

test_that("invalid inputs are rejected", {
  df <- sim_panel(N = 3, T = 30, seed = 8)
  expect_error(prepare_panel(~ x1, df, c("id", "time")), "two-sided")
  expect_error(prepare_panel(y ~ x1, as.matrix(df), c("id", "time")), "data frame")
  expect_error(prepare_panel(y ~ x1, df, c("id", "period")), "not found in `data`: period")
  expect_error(prepare_panel(y ~ x1, df, "id"), "`index` must name")
  expect_error(prepare_panel(y ~ 1, df, c("id", "time")), "at least one regressor")
  inf <- df
  inf$x1[5] <- Inf
  expect_error(prepare_panel(y ~ x1, inf, c("id", "time")), "infinite")
  chr <- transform(df, x1 = as.character(x1))
  expect_error(prepare_panel(y ~ x1, chr, c("id", "time")), "regressor `x1` must be a numeric")
  frac <- transform(df, time = time + 0.5 * (id == 2))
  expect_error(prepare_panel(y ~ x1, frac, c("id", "time")), "whole numbers")
  na_id <- df
  na_id$id[3] <- NA
  expect_error(prepare_panel(y ~ x1, na_id, c("id", "time")), "id variable `id` has missing")
})

test_that("units without complete observations are dropped with a warning", {
  df <- sim_panel(N = 4, T = 30, seed = 9)
  df$y[df$id == 3] <- NA
  expect_warning(p <- prepare_panel(y ~ x1, df, c("id", "time")), "Dropping unit .*: 3")
  expect_equal(p$ids, c(1, 2, 4))
})

test_that("`.` in the formula excludes the index columns", {
  df <- sim_panel(N = 3, T = 30, K = 2, seed = 10)
  p <- prepare_panel(y ~ ., df, c("id", "time"))
  expect_identical(p$x_names, c("x1", "x2"))
})
