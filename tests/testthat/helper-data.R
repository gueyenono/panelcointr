# Simulated panels for the tests: N units of I(1) series in long format
# (id, time, y, x1..xK), balanced or with unit-specific start and length.
sim_panel <- function(N = 6, T = 40, K = 1, balanced = TRUE, seed = 1,
                      coint = FALSE, drift = 0) {
  withr_seed <- function(code) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
    code
  }
  withr_seed({
    rows <- vector("list", N)
    for (i in seq_len(N)) {
      Ti <- if (balanced) T else sample((T - 10):(T + 4), 1)
      start <- if (balanced) 1L else sample(1:6, 1)
      e <- matrix(stats::rnorm(Ti * K), Ti, K)
      X <- apply(e, 2, cumsum) + outer(seq_len(Ti), rep(drift, K))
      v <- as.numeric(stats::filter(stats::rnorm(Ti), 0.4, method = "recursive"))
      y <- if (coint) {
        drop(X %*% stats::runif(K, 0.5, 1.5)) +
          as.numeric(stats::filter(stats::rnorm(Ti), 0.5, method = "recursive"))
      } else {
        cumsum(v + 0.5 * e[, 1]) + drift * seq_len(Ti)
      }
      rows[[i]] <- data.frame(id = i, time = start:(start + Ti - 1L), y = y, X)
    }
    df <- do.call(rbind, rows)
    names(df) <- c("id", "time", "y", paste0("x", seq_len(K)))
    df
  })
}

xfml <- function(K) stats::reformulate(paste0("x", seq_len(K)), response = "y")

hce_fit <- function(...) {
  pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
                    deterministic = "trend", ...)
}
