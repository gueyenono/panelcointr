# P0: spot checks of the transcribed literals against the papers.
test_that("the published tables are transcribed as printed", {
  expect_equal(bcs_remark1$c, c(mean = -0.424, var = 0.964))
  expect_equal(bcs_remark1$tau, c(mean = -1.535, var = 0.341))
  expect_equal(bcs_table2_one[, "0.5" == format(bcs_lambda_grid)], c(mean = -1.972, var = 0.339))
  # The paper's example: (0.3, 0.7) has mean -2.315 and variance 0.337.
  expect_equal(bcs_table2_two[3, 7], -2.315)
  expect_equal(bcs_table2_two[7, 3], 0.337)
  expect_true(all(bcs_table2_two[upper.tri(bcs_table2_two)] < -1.8))
  expect_true(all(bcs_table2_two[lower.tri(bcs_table2_two)] > 0.3))
  # Symmetry of the two-break moments under lambda -> 1 - lambda, to the
  # Monte Carlo error of the paper.
  m <- bcs_table2_two
  for (i in 1:8) for (j in (i + 1):9) {
    expect_lt(abs(m[i, j] - m[10 - j, 10 - i]), 0.01)
  }
  expect_equal(bcs_table1["1", , "0.1"], c(`0.01` = -32.163, `0.05` = -23.629, `0.1` = -19.865))
  expect_equal(bcs_table1["6", , "0.9"], c(`0.01` = -82.664, `0.05` = -71.518, `0.1` = -66.449))
  expect_true(all(bcs_table1[, "0.01", ] < bcs_table1[, "0.05", ]))
  expect_true(all(bcs_table1[, "0.05", ] < bcs_table1[, "0.1", ]))
  expect_equal(bai_ng_table1$c["1", ], c(`0.01` = -20.151, `0.05` = -13.730, `0.1` = -11.022))
  expect_equal(bai_ng_table1$tau["6", ], c(`0.01` = -74.251, `0.05` = -64.393, `0.1` = -59.555))
  expect_equal(bcs_table3$c[2, ], c(-2.824, -2.402, -2.113, -1.759))
  expect_equal(bcs_table3$tau[3, ], c(-2.619, -2.269, -1.931, -1.506))
  # WP Table 9: its T = 1000 rows are the JAE's asymptotic moments.
  expect_equal(bcs_wp_table9$c[, 5], c(mean = -0.424, var = 0.964))
  expect_equal(bcs_wp_table9$gamma_mean[5, ], unname(bcs_table2_one["mean", ]))
})

test_that("moments come from the tables on the grid and are simulated off it", {
  ctl <- banerjee_control(sim_reps = 2000L, sim_steps = 200L)
  expect_equal(banerjee_moments("gamma", 0.2, 100, ctl)$mean, -1.809)
  expect_equal(banerjee_moments("gamma", c(0.1, 0.9), 100, ctl)$var, 0.442)
  sim <- banerjee_moments("gamma", 0.25, 100, ctl)
  expect_equal(sim$source, "simulated")
  expect_lt(abs(sim$mean - (-1.809 - 1.902) / 2), 0.08)
  expect_true(on_grid(c(0.3, 0.7)))
  expect_false(on_grid(0.35))
  expect_false(on_grid(1))
})

test_that("the simulated quantile tables reproduce the Dickey-Fuller critical values", {
  q <- bcs_quantiles
  at <- function(v, p) v[which.min(abs(q$probs - p))]
  # MacKinnon (1996): no constant -1.941, constant -2.862, trend -3.410 at 5%.
  expect_lt(abs(at(q$idio_c, 0.05) + 1.941), 0.01)
  expect_lt(abs(at(q$factor_c, 0.05) + 2.862), 0.01)
  expect_lt(abs(at(q$factor_tau, 0.05) + 3.410), 0.015)
  expect_true(all(diff(q$idio_tau) > 0))
})
