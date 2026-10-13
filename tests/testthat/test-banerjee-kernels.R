# T3: the Rust kernels of pcoint_banerjee() against direct R computations.

rand_sym <- function(n, seed) {
  set.seed(seed)
  a <- matrix(rnorm(n * n), n)
  a + t(a)
}

test_that("the symmetric eigensolver matches eigen()", {
  for (n in c(1, 2, 5, 40, 150)) {
    a <- rand_sym(n, n)
    e <- eigen_sym_rs(as.vector(a), n)
    ref <- eigen(a, symmetric = TRUE)
    expect_equal(e$values, ref$values, tolerance = 1e-12)
    V <- matrix(e$vectors, n)
    expect_equal(crossprod(V), diag(n), tolerance = 1e-12)
    expect_equal(a %*% V, V %*% diag(e$values, n), tolerance = 1e-11)
  }
  # Clustered and repeated eigenvalues.
  set.seed(4)
  Q <- qr.Q(qr(matrix(rnorm(64), 8)))
  a <- Q %*% diag(c(5, 5, 5, 1 + 1e-10, 1, 0, 0, -2)) %*% t(Q)
  e <- eigen_sym_rs(as.vector(a), 8)
  expect_equal(e$values, c(5, 5, 5, 1 + 1e-10, 1, 0, 0, -2), tolerance = 1e-12)
  V <- matrix(e$vectors, 8)
  expect_equal(crossprod(V), diag(8), tolerance = 1e-12)
})

test_that("the generalised eigensolver matches eigen(A B^-1)", {
  set.seed(8)
  n <- 4
  A <- rand_sym(n, 3)
  M <- matrix(rnorm(n * n), n)
  B <- crossprod(M) + diag(n)
  e <- eigen_gen_rs(as.vector(A), as.vector(B), n)
  expect_equal(e$values, sort(Re(eigen(A %*% solve(B))$values), decreasing = TRUE), tolerance = 1e-12)
})

test_that("principal components match svd() on rotation-invariant quantities", {
  set.seed(5)
  for (dims in list(c(40, 12), c(15, 30))) {
    n <- dims[1]; N <- dims[2]
    Y <- matrix(rnorm(n * 2), n) %*% matrix(rnorm(2 * N), 2) + matrix(rnorm(n * N, sd = 0.3), n)
    p <- pca_rs(as.vector(Y), n, "bic3", 6L, 2L)
    sv <- svd(Y)
    expect_equal(p$values[1:min(n, N)], sv$d^2, tolerance = 1e-11)
    f <- matrix(p$f, n)
    expect_equal(crossprod(f) / n, diag(2), tolerance = 1e-12)
    common <- sv$u[, 1:2] %*% diag(sv$d[1:2]) %*% t(sv$v[, 1:2])
    expect_equal(matrix(p$resid, n), Y - common, tolerance = 1e-10)
    expect_equal(p$r, 2L)
    ref <- ref_pca(Y, r = 2)
    expect_equal(f, ref$f, tolerance = 1e-10)
    expect_equal(p$criteria, ref_criteria(ref$values, n, N, 6, "bic3"), tolerance = 1e-12)
    for (crit in c("ic1", "ic2", "ic3", "pc1", "pc2", "pc3")) {
      expect_equal(pca_rs(as.vector(Y), n, crit, 6L, -1L)$criteria,
                   ref_criteria(ref$values, n, N, 6, crit), tolerance = 1e-12)
    }
  }
})

adf_settings <- function(kmin, kmax, rule = "tsig", sample = "common") {
  list(lag_min = kmin, lag_max = kmax, lag_criterion = rule, tsig_cut = qnorm(0.95),
       lag_sample = sample)
}

test_that("the ADF regression and its lag selection match lm.fit", {
  for (seed in 1:30) {
    set.seed(seed)
    e <- cumsum(as.numeric(stats::filter(rnorm(80), 0.5, method = "recursive")))
    for (rule in c("tsig", "aic", "bic")) for (sample in c("common", "varying")) {
      if (rule == "tsig" && ref_near_tie_adf(e, 0, 5, qnorm(0.95))) next
      a <- adf_rs(e, adf_settings(0L, 5L, rule, sample))
      r <- ref_adf(e, 0, 5, rule, qnorm(0.95), sample)
      expect_equal(a$lags, r$k)
      expect_equal(a$t, r$t, tolerance = 1e-10)
      expect_equal(a$n, r$n)
      expect_equal(a$resid, unname(r$resid), tolerance = 1e-10)
    }
  }
  e <- cumsum(rnorm(50))
  expect_equal(adf_rs(e, adf_settings(2L, 2L))$lags, 2L)
})

test_that("the Gram break search matches brute force", {
  for (seed in 1:6) {
    df <- sim_panel(N = 1, T = 45, K = 2, seed = seed)
    u <- ref_bcs_units(df, 2)[[1]]
    for (model in c(1, 2, 4, 5)) for (m in 1:2) for (dols in list(c(0, 0), c(1, 2))) {
      b <- break_search_rs(u$y, as.vector(u$X), 2L, model, m, 6L, dols[1], dols[2], NULL)
      r <- ref_break_search(u, model, m, 6, dols)
      expect_equal(b$pos + 1L, r$tb)
      expect_equal(b$ssr, r$ssr, tolerance = 1e-9)
      expect_equal(b$n_candidates, length(r$cands))
    }
  }
})

test_that("factor detrending, the factor ADF and the MQ tests match R", {
  set.seed(11)
  n <- 80
  for (j in c("c", "tau", "gamma")) for (r in 1:3) {
    f <- matrix(rnorm(n * r), n)
    det <- ref_detrend_basis(n, j, if (j == "gamma") c(30, 55) else integer())
    s <- c(adf_settings(0L, 4L), list(mq_bandwidth = 4L, mq_max_var = 4L))
    a <- factor_tests_rs(as.vector(f), n, as.vector(det), s)
    ref <- ref_factor_tests(f, det, list(kmin = 0, kmax = 4, rule = "tsig", cut = qnorm(0.95)), 4, 4)
    expect_equal(matrix(a$fd, n), ref$fd, tolerance = 1e-10)
    if (r == 1) {
      expect_equal(a$adf_t, ref$adf$t, tolerance = 1e-10)
      expect_equal(a$adf_lags, ref$adf$k)
    } else {
      expect_equal(a$mq_c, vapply(ref$mq, `[[`, numeric(1), "mq_c"), tolerance = 1e-9)
      expect_equal(a$mq_f, vapply(ref$mq, `[[`, numeric(1), "mq_f"), tolerance = 1e-9)
      expect_equal(a$var_order, vapply(ref$mq, `[[`, numeric(1), "p"))
    }
  }
})

test_that("the random number streams are uniform and normal and reproducible", {
  a <- rng_rs(1, 0, 20000)
  expect_identical(a, rng_rs(1, 0, 20000))
  expect_false(identical(a$uniform, rng_rs(1, 1, 20000)$uniform))
  expect_true(all(a$uniform >= 0 & a$uniform < 1))
  expect_gt(stats::ks.test(a$uniform, "punif")$p.value, 1e-4)
  expect_gt(stats::ks.test(a$normal, "pnorm")$p.value, 1e-4)
})

test_that("the null simulations do not depend on the number of threads", {
  a <- banerjee_null_rs("idio", "gamma", c(0.3, 0.7), 100L, 2500L, 1L, 1, 1L)
  b <- banerjee_null_rs("idio", "gamma", c(0.3, 0.7), 100L, 2500L, 1L, 1, 3L)
  expect_identical(a, b)
  m <- banerjee_null_rs("mq", "tau", numeric(), 60L, 1200L, 3L, 2, 2L)
  expect_length(m, 3600L)
  expect_identical(m, banerjee_null_rs("mq", "tau", numeric(), 60L, 1200L, 3L, 2, 1L))
})

test_that("the CD statistic matches plm::pcdtest", {
  skip_if_not_installed("plm")
  set.seed(3)
  e <- matrix(rnorm(40 * 6), 40) + rnorm(40)
  ours <- cd_test(e)
  pd <- data.frame(id = rep(1:6, each = 40), time = rep(1:40, 6), e = as.vector(e))
  pd <- plm::pdata.frame(pd, index = c("id", "time"))
  ref <- plm::pcdtest(pd$e, test = "cd")
  expect_equal(ours$cd, unname(ref$statistic), tolerance = 1e-10)
  expect_equal(ours$cd_p_value, ref$p.value, tolerance = 1e-10)
})
