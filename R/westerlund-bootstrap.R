# Bootstrap of the Westerlund (2007) statistics; see the Bootstrap section of
# ?pcoint_westerlund. Rust prepares the state and computes the replicates;
# R draws the indices, so that set.seed() governs the results whatever the
# number of threads, and loops over chunks so the user can interrupt.
westerlund_bootstrap <- function(panel, labels, sel_lags, sel_leads, lags, leads, M, engine,
                                 s, B, n_threads, stats) {
  N <- length(panel$T)
  ones <- rep(1L, N)
  state <- westerlund_boot_prepare_rs(
    panel$y, as.vector(panel$X), ncol(panel$X), panel$offsets, labels, panel$start,
    as.integer(sel_lags), as.integer(sel_leads),
    lags[1L] * ones, lags[2L] * ones, leads[1L] * ones, leads[2L] * ones, M, engine,
    s$boot_reselect, s$boot_burn
  )
  info <- westerlund_boot_info_rs(state)
  n_periods <- as.integer(info$n_periods)
  len <- as.integer(info$len)
  if (n_periods < s$min_boot_periods) {
    abort("The bootstrap resamples the calendar periods in which every unit has a ",
          "first-stage residual, and this panel has ", n_periods, " (at least ",
          s$min_boot_periods, " are needed; see `min_boot_periods` in ",
          "`westerlund_control()`). The units' samples overlap too little.")
  }
  med <- stats::median(panel$T)
  if (n_periods < med / 2) {
    warn("Only ", n_periods, " calendar periods are common to all units (median T = ", med,
         "); the bootstrap resamples from these periods only.")
  }

  rng <- list(
    kind = RNGkind(),
    seed = if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv())
    }
  )
  out <- matrix(NA_real_, B, 4L, dimnames = list(NULL, names(stats)))
  show_progress <- interactive() && B > s$boot_chunk && isTRUE(getOption("panelcointr.progress", TRUE))
  pb <- if (show_progress) utils::txtProgressBar(0, B, style = 3)
  done <- 0L
  while (done < B) {
    m <- min(s$boot_chunk, B - done)
    # Replicate b uses the b-th block of `len` draws, whatever the chunk size.
    idx <- sample.int(n_periods, m * len, replace = TRUE)
    out[done + seq_len(m), ] <- matrix(westerlund_boot_rs(state, idx, m, n_threads), m, 4L,
                                       byrow = TRUE)
    done <- done + m
    if (!is.null(pb)) utils::setTxtProgressBar(pb, done)
  }
  if (!is.null(pb)) close(pb)

  below <- colSums(out <= matrix(stats, B, 4L, byrow = TRUE))
  p_value <- if (s$boot_pvalue == "plus_one") (1 + below) / (B + 1) else below / B
  critical <- apply(out, 2L, stats::quantile, probs = c(0.01, 0.05, 0.10), type = 1,
                    names = FALSE)
  critical <- matrix(critical, 3L, 4L, dimnames = list(c("1%", "5%", "10%"), names(stats)))
  list(
    B = B,
    p_value = unname(p_value),
    critical = critical,
    stats = if (isTRUE(s$keep_boot)) out,
    n_periods = n_periods,
    draws_per_replicate = len,
    rng = rng
  )
}
