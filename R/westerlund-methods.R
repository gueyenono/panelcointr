plural <- function(n, word) paste0(n, " ", word, if (n != 1L) "s")

order_label <- function(range, chosen, criterion) {
  if (range[1L] == range[2L]) return(as.character(range[1L]))
  crit <- c(aic = "AIC", bic = "BIC", westerlund = "xtwest's westerlund criterion")[[criterion]]
  sprintf("%d to %d by %s (mean %.2f)", range[1L], range[2L], crit, mean(chosen))
}

westerlund_header <- function(x) {
  s <- x$settings
  u <- x$units
  det <- c(none = "no deterministic terms", constant = "constant",
           trend = "constant and trend")[[x$deterministic]]
  Trange <- range(u$T)
  Tlab <- if (Trange[1L] == Trange[2L]) Trange[1L] else paste(Trange, collapse = " to ")
  M <- range(u$lrwindow)
  Mlab <- if (M[1L] == M[2L]) M[1L] else paste(M, collapse = " to ")
  cat("\nWesterlund (2007) error-correction panel cointegration tests\n")
  cat("H0: no cointegration\n\n")
  cat(sprintf("%s, %s (%s); T = %s; deterministic: %s\n", plural(x$N, "unit"),
              plural(x$K, "regressor"), paste(x$x_names, collapse = ", "), Tlab, det))
  cat("Lags:", order_label(s$lags, u$lags, s$criterion), "\n")
  cat("Leads:", order_label(s$leads, u$leads, s$criterion), "\n")
  cat("Bartlett bandwidth:", Mlab, if (is.null(s$lrwindow)) "(4 (T/100)^(2/9) rule)", "\n")
  if (s$compat != "none") cat("Compatibility:", s$compat, "\n")
  if (!is.null(x$boot)) cat("Bootstrap:", plural(x$boot$B, "replication"), "\n")
  cat("\n")
}

#' Print, summarise and plot Westerlund (2007) test results
#'
#' `print()` shows the statistics in the layout of Stata's `xtwest`: value,
#' z-value, asymptotic p-value and, with a bootstrap, the robust
#' (bootstrap) p-value. `summary()` adds the bootstrap critical values, the
#' distribution of the selected lags and leads, the unit-level results and
#' the mean-group error-correction model. `plot()` draws the bootstrap
#' distributions of the four statistics with their observed values.
#'
#' The mean-group model averages the units' coefficients (Pesaran and
#' Smith 1995), with standard errors from their dispersion across units.
#' Its long-run part reports, for each regressor and deterministic term,
#' the average of \eqn{-\hat\lambda_i / \hat\alpha_i}{-lambda_i / alpha_i}
#' (\eqn{-\hat\delta_i / \hat\alpha_i}{-delta_i / alpha_i} for the
#' deterministic terms, as `xtwest`'s `mg` option does); its short-run
#' part reports the error-correction coefficient and the differenced terms
#' present in every unit's model.
#'
#' @param x,object A `pcoint_westerlund` object.
#' @param digits Number of decimals.
#' @param which Statistics to plot.
#' @param ... Further arguments, passed to [graphics::hist()] by `plot()`
#'   and ignored otherwise.
#' @return `print()` and `plot()` return `x` invisibly; `summary()` returns
#'   an object of class `summary.pcoint_westerlund`.
#' @seealso [pcoint_westerlund()]
#' @name westerlund-methods
#' @examples
#' fit <- pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
#'                          deterministic = "trend", lags = c(1, 3), leads = c(0, 3))
#' fit
#' summary(fit)
NULL

#' @rdname westerlund-methods
#' @export
print.pcoint_westerlund <- function(x, digits = 3L, ...) {
  westerlund_header(x)
  print(format_statistics(x$statistics, digits), row.names = FALSE, right = TRUE)
  cat("\nGt, Ga: H1 alpha_i < 0 for some i.  Pt, Pa: H1 alpha_i = alpha < 0 for all i.\n")
  invisible(x)
}

# Mean-group estimates: average of the units' coefficients, standard
# errors from their dispersion (xtwest's `mg` option).
mean_group <- function(b) {
  b <- b[, colSums(is.na(b)) == 0L, drop = FALSE]
  N <- nrow(b)
  est <- colMeans(b)
  se <- if (N > 1L) apply(b, 2L, stats::sd) / sqrt(N) else rep(NA_real_, ncol(b))
  z <- est / se
  data.frame(estimate = est, std_error = se, z = z, p_value = 2 * stats::pnorm(-abs(z)),
             row.names = colnames(b))
}

#' @rdname westerlund-methods
#' @export
summary.pcoint_westerlund <- function(object, ...) {
  b <- object$coefficients
  ec <- paste0("L.", object$y_name)
  alpha <- b[, ec]
  lr_cols <- c(intersect(c("constant", "trend"), colnames(b)), paste0("L.", object$x_names))
  long_run <- -b[, lr_cols, drop = FALSE] / alpha
  colnames(long_run) <- c(intersect(c("constant", "trend"), colnames(b)), object$x_names)
  sr <- b[, setdiff(colnames(b), lr_cols), drop = FALSE]
  sr <- sr[, c(ec, setdiff(colnames(sr), ec)), drop = FALSE]
  colnames(sr)[1L] <- "ec"
  structure(list(
    object = object,
    lags = table(lags = object$units$lags),
    leads = table(leads = object$units$leads),
    mean_group = list(long_run = mean_group(long_run), short_run = mean_group(sr))
  ), class = "summary.pcoint_westerlund")
}

#' @rdname westerlund-methods
#' @export
print.summary.pcoint_westerlund <- function(x, digits = 3L, ...) {
  obj <- x$object
  print(obj, digits = digits)
  if (!is.null(obj$boot)) {
    cat("\nBootstrap critical values\n")
    print(round(obj$boot$critical, digits))
  }
  cat("\nSelected lags and leads (number of units)\n")
  print(x$lags)
  print(x$leads)
  cat("\nMean-group error-correction model\n")
  cat("Long run (-coefficient / alpha_i):\n")
  print(round(x$mean_group$long_run, digits))
  cat("Short run:\n")
  print(round(x$mean_group$short_run, digits))
  cat("\nUnits\n")
  u <- obj$units
  num <- vapply(u, is.double, logical(1))
  u[num] <- lapply(u[num], round, digits = digits)
  print(u, row.names = FALSE)
  invisible(x)
}

#' @rdname westerlund-methods
#' @export
plot.pcoint_westerlund <- function(x, which = c("Gt", "Ga", "Pt", "Pa"), ...) {
  if (is.null(x$boot) || is.null(x$boot$stats)) {
    abort("No bootstrap distribution to plot: rerun with `bootstrap > 0` ",
          "(and `keep_boot = TRUE`).")
  }
  which <- match.arg(which, several.ok = TRUE)
  old <- graphics::par(mfrow = grDevices::n2mfrow(length(which)))
  on.exit(graphics::par(old))
  for (s in which) {
    obs <- x$statistics$value[x$statistics$statistic == s]
    draws <- x$boot$stats[, s]
    graphics::hist(draws, breaks = 30, main = s, xlab = "Bootstrap distribution",
                   xlim = range(draws, obs), border = "grey40", col = "grey85", ...)
    graphics::abline(v = obs, lwd = 2)
    graphics::abline(v = x$boot$critical["5%", s], lty = 2)
  }
  graphics::legend("topright", c("observed", "5% critical value"), lty = c(1, 2),
                   lwd = c(2, 1), bty = "n", cex = 0.8)
  invisible(x)
}

#' Tidy and glance methods for panel cointegration tests
#'
#' Methods for the [generics::tidy()] and [generics::glance()] generics
#' (from the broom family), available when the generics package is
#' installed.
#'
#' @param x A `pcoint_westerlund` object.
#' @param ... Ignored.
#' @return `tidy()`: a data frame with one row per statistic and columns
#'   `test`, `statistic`, `z`, `p.value`, `p.value.boot` and `alternative`.
#'   `glance()`: a one-row data frame describing the test.
#' @name westerlund-tidiers
#' @examples
#' fit <- pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
#'                          deterministic = "trend", lags = 1, leads = 1)
#' if (requireNamespace("generics", quietly = TRUE)) {
#'   generics::tidy(fit)
#'   generics::glance(fit)
#' }
NULL

#' @rdname westerlund-tidiers
#' @exportS3Method generics::tidy
tidy.pcoint_westerlund <- function(x, ...) {
  s <- x$statistics
  data.frame(test = s$statistic, statistic = s$value, z = s$z, p.value = s$p_value,
             p.value.boot = s$p_boot, alternative = s$alternative)
}

#' @rdname westerlund-tidiers
#' @exportS3Method generics::glance
glance.pcoint_westerlund <- function(x, ...) {
  data.frame(
    method = "Westerlund (2007) ECM panel cointegration tests",
    n_units = x$N, n_regressors = x$K, deterministic = x$deterministic,
    mean_lags = mean(x$units$lags), mean_leads = mean(x$units$leads),
    n_boot = if (is.null(x$boot)) 0L else x$boot$B,
    compat = x$settings$compat
  )
}
