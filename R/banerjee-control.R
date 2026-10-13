#' Control settings for the Banerjee and Carrion-i-Silvestre (2015) tests
#'
#' Fine-grained choices for [pcoint_banerjee()]. The defaults follow the
#' paper; the settings actually used are stored in the `settings` element of
#' the result.
#'
#' @param lag_criterion Rule for the ADF lag orders: `"tsig"` (default), the
#'   general-to-specific rule of Ng and Perron (1995) used by the paper,
#'   which starts at the largest order and drops the last lag while its
#'   t-ratio is not significant at level `tsig_level` (two-sided); or
#'   `"aic"` / `"bic"`.
#' @param tsig_level Significance level of the t-sig rule (default 0.10, a
#'   cut-off of 1.645).
#' @param lag_sample Sample on which lag orders are compared: `"common"`,
#'   the sample of the largest order for every candidate (default); or
#'   `"varying"`, each candidate's own sample. The chosen order is always
#'   re-estimated on its own sample.
#' @param max_lags Largest ADF lag order when `lags = NULL` in
#'   [pcoint_banerjee()]. `NULL` uses
#'   \eqn{4 \lceil (\min(N, T)/100)^{1/4} \rceil}{4 ceiling((min(N, T)/100)^(1/4))},
#'   the rule of the paper's simulations (4 for any panel with
#'   \eqn{\min(N, T) \le 100}{min(N, T) <= 100}).
#' @param pca_scale Divide each unit's first-difference residuals by their
#'   standard deviation before the principal components (the "transformed"
#'   variant of Banerjee and Carrion-i-Silvestre 2017). Only the factor
#'   estimates and the number of factors change; the idiosyncratic ADF
#'   statistics are scale invariant. The principal components inside the
#'   estimation of heterogeneous break dates are always unscaled, so that
#'   every iteration lowers the sum of squared residuals.
#' @param moments Moments that standardise the panel statistic:
#'   `"asymptotic"` (default; the paper's Remark 1 and Table II), or
#'   `"finite_sample"`, the T-specific moments of the working-paper version
#'   (BCS 2006, Table 9), interpolated linearly in \eqn{1/T}{1/T} and
#'   available for Models 1, 2, 4 and 5 and one break on the 0.1 grid in
#'   Models 3 and 6.
#' @param sim_reps,sim_steps Replications and length of the random walks
#'   used to simulate null distributions that are not tabulated (Models 3
#'   and 6). Defaults: the paper's 100,000 and 1000.
#' @param iter_max,iter_tol Maximum number of iterations and relative
#'   tolerance on the sum of squared residuals in the estimation of
#'   heterogeneous break dates.
#' @param mq_level Level of the sequential MQ tests: 0.01, 0.05 (default)
#'   or 0.10.
#' @param mq_bandwidth Bartlett bandwidth \eqn{J}{J} of the MQ_c test. `NULL`
#'   uses \eqn{4 \lceil (\min(N, T)/100)^{1/4} \rceil}{4 ceiling((min(N, T)/100)^(1/4))}.
#' @param mq_max_var Largest VAR order of the MQ_f test, selected by BIC.
#'   `NULL` uses \eqn{4 \lceil (T/100)^{1/4} \rceil}{4 ceiling((T/100)^(1/4))}.
#' @param min_df Minimum residual degrees of freedom of every unit's
#'   first-difference and ADF regressions.
#'
#' @return An object of class `banerjee_control`.
#' @seealso [pcoint_banerjee()]
#' @export
#' @examples
#' banerjee_control(lag_criterion = "bic", pca_scale = TRUE)
banerjee_control <- function(lag_criterion = c("tsig", "bic", "aic"), tsig_level = 0.10,
                             lag_sample = c("common", "varying"), max_lags = NULL,
                             pca_scale = FALSE, moments = c("asymptotic", "finite_sample"),
                             sim_reps = 100000L, sim_steps = 1000L, iter_max = 100L,
                             iter_tol = 1e-6, mq_level = 0.05, mq_bandwidth = NULL,
                             mq_max_var = NULL, min_df = 1L) {
  lag_criterion <- match.arg(lag_criterion)
  lag_sample <- match.arg(lag_sample)
  moments <- match.arg(moments)
  count <- function(x, name, min, null_ok = FALSE) {
    if (null_ok && is.null(x)) return(NULL)
    if (length(x) != 1L || !is_whole(x) || x < min) {
      abort("`", name, "` must be a whole number of at least ", min,
            if (null_ok) " (or NULL)", ".")
    }
    as.integer(round(x))
  }
  prob <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
      abort("`", name, "` must be a number between 0 and 1.")
    }
    x
  }
  if (!is.logical(pca_scale) || length(pca_scale) != 1L || is.na(pca_scale)) {
    abort("`pca_scale` must be TRUE or FALSE.")
  }
  if (!is.numeric(mq_level) || length(mq_level) != 1L ||
      !any(abs(mq_level - c(0.01, 0.05, 0.10)) < 1e-12)) {
    abort("`mq_level` must be 0.01, 0.05 or 0.10, the levels of the published tables.")
  }
  structure(list(
    lag_criterion = lag_criterion,
    tsig_level = prob(tsig_level, "tsig_level"),
    lag_sample = lag_sample,
    max_lags = count(max_lags, "max_lags", 0L, null_ok = TRUE),
    pca_scale = pca_scale,
    moments = moments,
    sim_reps = count(sim_reps, "sim_reps", 100L),
    sim_steps = count(sim_steps, "sim_steps", 50L),
    iter_max = count(iter_max, "iter_max", 1L),
    iter_tol = prob(iter_tol, "iter_tol"),
    mq_level = mq_level,
    mq_bandwidth = count(mq_bandwidth, "mq_bandwidth", 0L, null_ok = TRUE),
    mq_max_var = count(mq_max_var, "mq_max_var", 0L, null_ok = TRUE),
    min_df = count(min_df, "min_df", 1L)
  ), class = "banerjee_control")
}

#' @export
print.banerjee_control <- function(x, ...) {
  cat("<banerjee_control>\n")
  for (n in names(x)) {
    v <- x[[n]]
    cat("  ", n, ": ", if (is.null(v)) "NULL (default rule)" else format(v), "\n", sep = "")
  }
  invisible(x)
}
