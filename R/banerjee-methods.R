banerjee_model_label <- function(model) {
  c("constant, level shifts",
    "constant and trend, level shifts",
    "constant and trend, level and trend shifts",
    "constant, level shifts; shifts in the cointegrating vector",
    "constant and trend, level shifts; shifts in the cointegrating vector",
    "constant and trend, level and trend shifts; shifts in the cointegrating vector")[model]
}

banerjee_breaks_label <- function(x) {
  b <- x$breaks
  switch(b$type,
    none = "none",
    known = if (isTRUE(b$common)) {
      paste0("known, common: after ", paste(format(unlist(b$dates[1L, -1L])), collapse = ", "))
    } else {
      paste0("known, unit-specific (up to ", b$m, " per unit)")
    },
    heterogeneous = sprintf("%s estimated per unit (heterogeneous; %s after %d iteration%s)",
                            b$m, if (b$converged) "converged" else "not converged",
                            b$iterations, if (b$iterations != 1L) "s" else ""),
    homogeneous = paste0(b$m, " estimated, common (homogeneous): after ",
                         paste(format(unlist(b$dates[1L, -1L])), collapse = ", "))
  )
}

banerjee_lags_label <- function(x) {
  s <- x$settings
  l <- x$units$lags
  if (s$lags[1L] == s$lags[2L]) return(as.character(s$lags[1L]))
  crit <- c(tsig = "t-sig", aic = "AIC", bic = "BIC")[[s$lag_criterion]]
  sprintf("%d to %d by %s (mean %.2f)", s$lags[1L], s$lags[2L], crit, mean(l))
}

banerjee_factor_line <- function(x, digits) {
  f <- x$factors
  how <- if (f$criterion == "fixed") "fixed" else paste0(toupper(f$criterion), ", at most ",
                                                         x$settings$max_factors)
  cat(sprintf("Common factors: r = %d (%s)\n", f$r, how))
  lev <- paste0(100 * x$settings$mq_level, "%")
  if (f$r == 1L) {
    a <- f$adf
    cat(sprintf("  ADF on the factor: %s (5%% critical value %s, p = %s, %d lag%s): %s\n",
                formatC(a$statistic, format = "f", digits = digits),
                formatC(a$critical[["5%"]], format = "f", digits = digits),
                formatC(a$p_value, format = "f", digits = digits), a$lags,
                if (a$lags != 1L) "s" else "",
                if (f$r1_c == 0L) "stationary" else "I(1)"))
  } else if (f$r > 1L) {
    cat(sprintf("  I(1) common factors at %s: r1 = %d by MQ_c, %d by MQ_f\n", lev, f$r1_c, f$r1_f))
  }
}

#' Print, summarise and plot Banerjee and Carrion-i-Silvestre (2015) tests
#'
#' `print()` shows the panel statistic, the common-factor tests, the CD
#' diagnostic and, for estimated common breaks, the finite-sample critical
#' values of the paper's Table III. `summary()` adds the unit-level results,
#' the factor-number criterion, the MQ sequence and the distribution of the
#' break dates. `plot()` draws the profile of the statistic over the
#' candidate break dates (estimated common breaks) or the estimated common
#' factors with the break dates.
#'
#' @param x,object A `pcoint_banerjee` object.
#' @param digits Number of decimals.
#' @param ... Ignored, or passed to [graphics::matplot()] by `plot()`.
#' @return `print()` and `plot()` return `x` invisibly; `summary()` returns
#'   an object of class `summary.pcoint_banerjee`.
#' @seealso [pcoint_banerjee()]
#' @name banerjee-methods
#' @examples
#' set.seed(2)
#' d <- do.call(rbind, lapply(1:10, function(i) {
#'   x <- cumsum(rnorm(50))
#'   data.frame(id = i, time = 1:50, y = x + cumsum(rnorm(50)), x = x)
#' }))
#' fit <- pcoint_banerjee(y ~ x, data = d, model = 2)
#' fit
#' summary(fit)
NULL

#' @rdname banerjee-methods
#' @export
print.pcoint_banerjee <- function(x, digits = 3L, ...) {
  cat("\nBanerjee and Carrion-i-Silvestre (2015) panel cointegration test\n")
  cat("H0: no cointegration (unit root in every idiosyncratic error)\n\n")
  cat(sprintf("%s, %s (%s); T = %d\n", plural(x$N, "unit"), plural(x$K, "regressor"),
              paste(x$x_names, collapse = ", "), x$T))
  cat(sprintf("Model %d: %s\n", x$model, banerjee_model_label(x$model)))
  cat("Breaks:", banerjee_breaks_label(x), "\n")
  cat("ADF lags:", banerjee_lags_label(x), "\n")
  if (sum(x$settings$dols)) {
    cat(sprintf("DOLS: %d lag(s) and %d lead(s) of D2.x\n", x$settings$dols[1L], x$settings$dols[2L]))
  }
  cat(sprintf("Moments: mean %.3f, variance %.3f (%s)\n\n", x$moments$mean, x$moments$var,
              x$moments$source))
  print(format_statistics(x$statistics, digits), row.names = FALSE, right = TRUE)
  b <- x$breaks
  if (!is.null(b$table3)) {
    cat(sprintf("\nCritical values of Z* (Table III, N = 100, T = %d): %s\n",
                bcs_table3$T[which.min(abs(bcs_table3$T - x$T))],
                paste(names(b$table3), formatC(b$table3, format = "f", digits = 3), collapse = ", ")))
    cat("The N(0, 1) p-value of the infimum is anti-conservative in finite samples.\n")
  }
  cat("\n")
  banerjee_factor_line(x, digits)
  d <- x$diagnostics
  if (!is.na(d$cd)) {
    cat(sprintf("CD test of the ADF residuals: %s (p = %s, %d periods)\n",
                formatC(d$cd, format = "f", digits = digits),
                formatC(d$cd_p_value, format = "f", digits = digits), d$cd_periods))
  }
  cat(sprintf("Individual ADF rejections at 5%%: %.1f%%\n", 100 * d$share_rejected))
  r1 <- max(x$factors$r1_c, x$factors$r1_f)
  if (r1 > 0L && x$statistics$p_value < 0.05) {
    cat("\nThe idiosyncratic errors look stationary but some common factors are I(1): the\n",
        "variables cointegrate only together with the non-stationary factors.\n", sep = "")
  }
  invisible(x)
}

#' @rdname banerjee-methods
#' @export
summary.pcoint_banerjee <- function(object, ...) {
  b <- object$breaks
  dist <- NULL
  if (!is.null(b$dates) && !isTRUE(b$common)) {
    dist <- lapply(b$dates[-1L], function(v) table(v, useNA = "no"))
  }
  structure(list(object = object, lags = table(lags = object$units$lags), break_dates = dist),
            class = "summary.pcoint_banerjee")
}

#' @rdname banerjee-methods
#' @export
print.summary.pcoint_banerjee <- function(x, digits = 3L, ...) {
  obj <- x$object
  print(obj, digits = digits)
  f <- obj$factors
  cat("\nFactor-number criterion", if (f$criterion != "fixed") paste0("(", toupper(f$criterion), ")"),
      "by number of factors\n")
  print(signif(f$criteria, 5))
  if (!is.null(f$mq)) {
    cat("\nMQ tests of H0: r1 = q (critical values at ", 100 * obj$settings$mq_level, "%, ",
        f$mq_source, ")\n", sep = "")
    m <- f$mq
    m[c("MQ_c", "MQ_f", "cv")] <- lapply(m[c("MQ_c", "MQ_f", "cv")], round, digits = digits)
    print(m, row.names = FALSE)
  }
  cat("\nSelected ADF lags (number of units)\n")
  print(x$lags)
  if (!is.null(x$break_dates)) {
    cat("\nEstimated or given break dates (number of units)\n")
    for (n in names(x$break_dates)) {
      cat(n, ":\n", sep = "")
      print(x$break_dates[[n]])
    }
  }
  cat("\nUnits\n")
  u <- obj$units
  num <- vapply(u, is.double, logical(1))
  u[num] <- lapply(u[num], round, digits = digits)
  print(u, row.names = FALSE)
  invisible(x)
}

#' @rdname banerjee-methods
#' @export
plot.pcoint_banerjee <- function(x, ...) {
  b <- x$breaks
  if (!is.null(b$profile) && b$m == 1L) {
    p <- b$profile
    graphics::plot(seq_len(nrow(p)), p$z, type = "l", xaxt = "n",
                   xlab = "Common break date", ylab = "Z(lambda)",
                   main = "Profile of the panel statistic", ...)
    at <- unique(round(seq(1, nrow(p), length.out = min(8L, nrow(p)))))
    graphics::axis(1, at = at, labels = format(p[[1L]][at]))
    best <- which.min(p$z)
    graphics::points(best, p$z[best], pch = 19)
    if (!is.null(b$table3)) graphics::abline(h = b$table3[["5%"]], lty = 2)
    return(invisible(x))
  }
  if (!is.null(b$profile)) {
    p <- b$profile
    graphics::plot(p$lambda_1, p$lambda_2, cex = 0.4 + 2 * (max(p$z, na.rm = TRUE) - p$z) /
                     diff(range(p$z, na.rm = TRUE)),
                   xlab = "lambda_1", ylab = "lambda_2", main = "Z(lambda): larger points are smaller")
    best <- which.min(p$z)
    graphics::points(p$lambda_1[best], p$lambda_2[best], pch = 19, col = "red")
    return(invisible(x))
  }
  f <- x$factors
  if (f$r == 0L) {
    abort("No common factors to plot (r = 0).")
  }
  graphics::matplot(f$F, type = "l", lty = 1, xaxt = "n", xlab = "", ylab = "Cumulated factor",
                    main = "Estimated common factors", ...)
  rn <- rownames(f$F)
  at <- unique(round(seq(1, length(rn), length.out = min(8L, length(rn)))))
  graphics::axis(1, at = at, labels = rn[at])
  if (isTRUE(b$common) && !is.null(b$positions)) {
    graphics::abline(v = b$positions[1L, ] - (x$T - nrow(f$F)) + 0.5, lty = 2)
  }
  invisible(x)
}

#' @rdname westerlund-tidiers
#' @exportS3Method generics::tidy
tidy.pcoint_banerjee <- function(x, ...) {
  s <- x$statistics
  data.frame(test = s$statistic, statistic = s$value, z = s$z, p.value = s$p_value,
             alternative = s$alternative)
}

#' @rdname westerlund-tidiers
#' @exportS3Method generics::glance
glance.pcoint_banerjee <- function(x, ...) {
  data.frame(
    method = "Banerjee and Carrion-i-Silvestre (2015) panel cointegration test",
    n_units = x$N, n_regressors = x$K, n_periods = x$T, model = x$model,
    breaks = x$breaks$type, n_breaks = x$breaks$m, n_factors = x$factors$r,
    r1_c = x$factors$r1_c, r1_f = x$factors$r1_f, mean_lags = mean(x$units$lags),
    cd = x$diagnostics$cd, cd.p.value = x$diagnostics$cd_p_value
  )
}
