#' Westerlund (2007) error-correction panel cointegration tests
#'
#' Computes the four error-correction-based panel cointegration tests of
#' Westerlund (2007): the group-mean tests Gt and Ga and the panel tests Pt
#' and Pa, with asymptotic p-values from the paper's Table 1 and,
#' optionally, bootstrap p-values that are robust to cross-sectional
#' dependence.
#'
#' For each unit \eqn{i} the conditional error-correction model
#' \deqn{\Delta y_{it} = \delta_i' d_t + \alpha_i y_{i,t-1} + \lambda_i' x_{i,t-1}
#'   + \sum_{j=1}^{p_i} \alpha_{ij} \Delta y_{i,t-j}
#'   + \sum_{j=-q_i}^{p_i} \gamma_{ij}' \Delta x_{i,t-j} + e_{it}}{Delta y_it =
#'   delta_i' d_t + alpha_i y_i,t-1 + lambda_i' x_i,t-1 + sum alpha_ij Delta y_i,t-j
#'   + sum gamma_ij' Delta x_i,t-j + e_it}
#' is estimated by least squares. The null hypothesis is no error
#' correction, \eqn{\alpha_i = 0} for all \eqn{i}, that is, no
#' cointegration. The group-mean statistics test it against
#' \eqn{\alpha_i < 0} for at least one unit, the panel statistics against
#' \eqn{\alpha_i = \alpha < 0} for all units. All four statistics diverge
#' to minus infinity under the alternative, so the tests are left-tailed.
#'
#' Gt is the average of the units' t-ratios of \eqn{\hat\alpha_i} and Ga
#' the average of \eqn{\tilde T_i \hat\alpha_i / \hat\alpha_i(1)}, where
#' \eqn{\hat\alpha_i(1) = \hat\omega_{ui} / \hat\omega_{yi}} is the
#' semiparametric estimate built from Bartlett long-run variances of the
#' composite error and of \eqn{\Delta y}. Pt and Pa are the t-ratio and the
#' normalised pooled estimate of a common \eqn{\alpha}.
#'
#' @section Lags and leads:
#' `lags` and `leads` are each a single order, used for every unit, or a
#' range `c(min, max)`. With a range, each unit's orders minimise the
#' information criterion `criterion` over the grid. By default every
#' candidate is compared on the sample of the largest model (Ng and Perron
#' 2005), ties going to the smaller model; the chosen model is then
#' re-estimated on its own sample. Leads correct for regressors that are
#' weakly but not strictly exogenous (paper, Remark 5).
#'
#' @section Asymptotic p-values:
#' The statistics are standardised with the moments of Westerlund (2007,
#' Table 1) for the deterministic case and the number of regressors `K`.
#' The table covers `K` up to 6; with more regressors the z-values and
#' asymptotic p-values are `NA` (with a warning) and only the bootstrap is
#' available.
#'
#' @section Bootstrap:
#' With `bootstrap = B > 0`, p-values robust to cross-sectional dependence
#' are computed as in Westerlund (2007, Section 5), extended to leads:
#' 1. For each unit, `Delta y` is regressed on its own lags and on the
#'    leads and lags of `Delta x` at the unit's selected orders, imposing
#'    the null (a constant is added in the trend case).
#' 2. The residuals and `Delta x` of all units are stacked, period by
#'    period, over the calendar periods in which every unit has a residual,
#'    and centred.
#' 3. Each replication resamples whole periods with replacement (the same
#'    period for all units), which preserves the cross-sectional
#'    dependence; the indices are drawn with R's random number generator,
#'    so [set.seed()] makes the results reproducible for any `n_threads`.
#' 4. Bootstrap `Delta y` and `Delta x` are generated from the first-stage
#'    coefficients, starting at zero and discarding `boot_burn` periods
#'    (see [westerlund_control()]); the levels are their cumulative sums
#'    over each unit's own window, so unbalanced panels keep their pattern.
#' 5. The four statistics are recomputed with the same settings,
#'    re-selecting the lags and leads unless `boot_reselect = FALSE`.
#'
#' The bootstrap p-value of a statistic \eqn{S} is
#' \eqn{(1 + \#\{S^*_b \le S\}) / (B + 1)}. The 1%, 5% and 10% bootstrap
#' critical values are also returned. Only the bootstrap loop runs in
#' parallel, over `n_threads` threads.
#'
#' @section Compatibility with Stata's xtwest:
#' The defaults follow the paper. `compat = "xtwest"` reproduces the
#' statistics of the Stata command `xtwest` (Persyn and Westerlund 2008),
#' and `compat = "xtwest_westerlund"` those of its `westerlund` option
#' (which requires a constant and a single regressor, as in Stata). The
#' differences are listed in [westerlund_control()]; with `compat` set,
#' `lrwindow` defaults to 2 as in `xtwest`. The `deterministic` argument is
#' not changed by `compat` (`xtwest`'s default is `"none"`). The bootstrap
#' is always the algorithm above: it does not reproduce `xtwest`'s
#' bootstrap draws.
#'
#' @param formula A formula `y ~ x1 + x2` with the levels of the I(1)
#'   variables. The intercept is ignored: use `deterministic`.
#' @param data A data frame (including tibbles, data.tables and
#'   `plm::pdata.frame`s) in long format.
#' @param index The names of the unit and time variables. `NULL` uses the
#'   index of a `pdata.frame`, and otherwise the first two columns of
#'   `data`. Numeric time must step by one within each unit; dates,
#'   factors and character time are mapped to their rank among the time
#'   values present. Leading and trailing missing values of a unit are
#'   dropped; missing values or gaps inside a unit's sample are an error.
#' @param deterministic Deterministic terms of the model: `"constant"`
#'   (Case 2, the default), `"trend"` (constant and linear trend, Case 3)
#'   or `"none"` (Case 1). Case 1 is inconsistent against cointegration
#'   with a constant (paper, Remark 4).
#' @param lags Number of lags of the differences: a single order or a range
#'   `c(min, max)` searched by `criterion`. Required, because the results
#'   can be sensitive to it.
#' @param leads Number of leads of `Delta x`: a single order or a range
#'   `c(min, max)`. Default 0.
#' @param lrwindow Bandwidth of the Bartlett kernel in the long-run
#'   variances. `NULL` (default) uses \eqn{\lfloor 4 (T_i/100)^{2/9} \rfloor}{floor(4 (T_i/100)^(2/9))}
#'   for each unit, the rule of the paper's simulations (2 with `compat`).
#' @param criterion Information criterion for ranges of `lags` or `leads`:
#'   `"aic"` or `"bic"`.
#' @param bootstrap Number of bootstrap replications; 0 (default) for
#'   asymptotic p-values only.
#' @param compat `"none"` (default) for the paper's procedure, `"xtwest"`
#'   or `"xtwest_westerlund"` to reproduce Stata's `xtwest`.
#' @param control Further settings from [westerlund_control()].
#' @param n_threads Number of threads for the bootstrap; the default is
#'   `getOption("panelcointr.threads", 1L)`.
#'
#' @return An object of class `c("pcoint_westerlund", "pcoint_test")`, a
#'   list with:
#'   * `statistics`: data frame with the `statistic` (Gt, Ga, Pt, Pa), its
#'     `value`, the standardised value `z`, the asymptotic `p_value`, the
#'     bootstrap `p_boot` (`NA` without bootstrap) and the `alternative`.
#'   * `units`: one row per unit with its number of periods `T`,
#'     observations `n` and regressors `k` in the regression, selected
#'     `lags` and `leads`, bandwidth `lrwindow`, `alpha`
#'     (\eqn{\hat\alpha_i}), `se_alpha`, `t_alpha`, `alpha1`
#'     (\eqn{\hat\alpha_i(1)}) and the regression standard error `sigma`.
#'   * `coefficients`: matrix of the units' ECM coefficients (`NA` for
#'     terms outside a unit's model), named in Stata's time-series
#'     notation (`L.x` for \eqn{x_{t-1}}, `L2D.y` for \eqn{\Delta y_{t-2}},
#'     `F1D.x` for \eqn{\Delta x_{t+1}}).
#'   * `panel`: the pooled estimate `alpha`, its standard error `se`, the
#'     normaliser `t_bar` of Pa and the orders `lags` and `leads` of the
#'     panel regressions when they are common to all units.
#'   * `boot`: with a bootstrap, the number of replications `B`, the `B x 4`
#'     matrix of statistics `stats` (unless `keep_boot = FALSE`), the
#'     `critical` values, the number of resampled periods and the state of
#'     the random number generator before the draws.
#'   * `settings`: every setting used, after resolving `compat` and
#'     `control`.
#'   * `N`, `K`, `case` (1 to 3), `deterministic`, `call`, `formula`.
#'
#' @references
#' Westerlund, J. (2007). Testing for error correction in panel data.
#' *Oxford Bulletin of Economics and Statistics* 69(6), 709-748.
#' \doi{10.1111/j.1468-0084.2007.00477.x}
#'
#' Persyn, D. and Westerlund, J. (2008). Error-correction-based
#' cointegration tests for panel data. *Stata Journal* 8(2), 232-241.
#'
#' Ng, S. and Perron, P. (2005). A note on the selection of time series
#' models. *Oxford Bulletin of Economics and Statistics* 67(1), 115-134.
#'
#' @seealso [westerlund_control()], [summary.pcoint_westerlund()],
#'   [plot.pcoint_westerlund()], [hce]
#' @export
#' @examples
#' # Health care expenditure and GDP, 20 OECD countries (Westerlund 2007, Section 7)
#' pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
#'                   deterministic = "trend", lags = c(1, 3), leads = c(0, 3))
#'
#' # Stata's xtwest, as in Persyn and Westerlund (2008)
#' pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
#'                   deterministic = "trend", lags = 1, leads = 1, lrwindow = 3,
#'                   compat = "xtwest")
#'
#' \donttest{
#' # Bootstrap p-values robust to cross-sectional dependence
#' set.seed(1)
#' fit <- pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
#'                          deterministic = "trend", lags = 1, leads = 1,
#'                          bootstrap = 199, n_threads = 2)
#' fit
#' summary(fit)
#' plot(fit)
#' }
pcoint_westerlund <- function(formula, data, index = NULL,
                              deterministic = c("constant", "trend", "none"),
                              lags, leads = 0, lrwindow = NULL,
                              criterion = c("aic", "bic"), bootstrap = 0L,
                              compat = c("none", "xtwest", "xtwest_westerlund"),
                              control = westerlund_control(),
                              n_threads = getOption("panelcointr.threads", 1L)) {
  call <- match.call()
  deterministic <- match.arg(deterministic)
  criterion_given <- !missing(criterion)
  criterion <- match.arg(criterion)
  compat <- match.arg(compat)
  if (missing(lags)) {
    abort("`lags` is required: a single lag order, or `c(min, max)` to select ",
          "it by information criterion.")
  }
  lags <- order_range(lags, "lags")
  leads <- order_range(leads, "leads")
  if (!is.null(lrwindow) && (length(lrwindow) != 1L || !is_whole(lrwindow) || lrwindow < 0)) {
    abort("`lrwindow` must be NULL or a single non-negative whole number.")
  }
  if (length(bootstrap) != 1L || !is_whole(bootstrap) || bootstrap < 0) {
    abort("`bootstrap` must be a single non-negative whole number.")
  }
  if (length(n_threads) != 1L || !is_whole(n_threads) || n_threads < 1) {
    abort("`n_threads` must be a single positive whole number.")
  }
  if (!inherits(control, "westerlund_control")) {
    abort("`control` must be created by `westerlund_control()`.")
  }
  if (compat == "xtwest_westerlund" && criterion_given) {
    warn("`criterion` is ignored with `compat = \"xtwest_westerlund\"`, ",
         "which uses the information criterion of xtwest's `westerlund` option.")
  }
  s <- resolve_westerlund_settings(compat, control, deterministic, criterion)

  panel <- prepare_panel(formula, data, index)
  N <- length(panel$T)
  K <- ncol(panel$X)
  labels <- as.character(panel$ids)
  if (compat == "xtwest_westerlund") {
    if (deterministic == "none") {
      abort("`compat = \"xtwest_westerlund\"` needs a constant, as xtwest's `westerlund` ",
            "option: use `deterministic = \"constant\"` or `\"trend\"`.")
    }
    if (K > 1L) {
      abort("`compat = \"xtwest_westerlund\"` allows a single regressor, as xtwest's ",
            "`westerlund` option.")
    }
  }

  # Bandwidths and sample sizes.
  if (!is.null(lrwindow)) s$lrwindow <- as.integer(lrwindow)
  M <- if (is.null(s$lrwindow)) {
    as.integer(floor(4 * (panel$T / 100)^(2 / 9)))
  } else {
    rep(s$lrwindow, N)
  }
  s$lags <- lags
  s$leads <- leads
  dd <- match(deterministic, c("none", "constant", "trend")) - 1L
  k_max <- dd + K + lags[2L] + K * (lags[2L] + leads[2L] + 1L) + 1L
  need <- lags[2L] + leads[2L] + 1L + k_max + s$min_df
  short <- panel$T < need
  if (any(short)) {
    abort("Too few observations for ", lags[2L], " lag", if (lags[2L] != 1L) "s", " and ",
          leads[2L], " lead", if (leads[2L] != 1L) "s", ": every unit needs at least ", need,
          " periods, but ", sum(short), " unit", if (sum(short) > 1L) "s have" else " has",
          " fewer: ", enumerate(paste0(labels[short], " (T = ", panel$T[short], ")")), ".")
  }

  engine <- westerlund_engine_settings(s)
  ones <- rep(1L, N)
  res <- westerlund_rs(
    panel$y, as.vector(panel$X), K, panel$offsets, labels,
    lags[1L] * ones, lags[2L] * ones, leads[1L] * ones, leads[2L] * ones, M, engine
  )

  stat_names <- c("Gt", "Ga", "Pt", "Pa")
  stats <- stats::setNames(res$statistics, stat_names)
  moments <- westerlund_moments(deterministic, K, s$moments)
  if (is.null(moments)) {
    warn("No asymptotic moments are tabulated for ", K, " regressors (Westerlund 2007, ",
         "Table 1, covers 1 to 6): z-values and asymptotic p-values are NA. ",
         "Use `bootstrap` for inference.")
  }
  s$moment_values <- moments
  z <- westerlund_z(stats, N, moments)

  boot <- NULL
  p_boot <- rep(NA_real_, 4L)
  if (bootstrap > 0) {
    boot <- westerlund_bootstrap(
      panel, labels, res$units$lags, res$units$leads, lags, leads, M, engine, s,
      as.integer(bootstrap), as.integer(n_threads), stats
    )
    p_boot <- boot$p_value
  }

  units <- data.frame(
    id = panel$ids, T = panel$T, n = as.integer(res$units$n), k = as.integer(res$units$k),
    lags = as.integer(res$units$lags), leads = as.integer(res$units$leads), lrwindow = M,
    alpha = res$units$alpha, se_alpha = res$units$se_alpha,
    t_alpha = res$units$alpha / res$units$se_alpha, alpha1 = res$units$alpha1,
    sigma = res$units$sigma
  )
  width <- length(res$coef_names)
  coefficients <- matrix(res$coefficients, nrow = N, ncol = width, byrow = TRUE,
                         dimnames = list(labels, ecm_coef_names(dd, panel$y_name, panel$x_names,
                                                                lags[2L], leads[2L])))
  coefficients[is.nan(coefficients)] <- NA_real_
  common <- s$panel_lags == "mean" || (length(unique(units$lags)) == 1L &&
                                         length(unique(units$leads)) == 1L)

  new_pcoint_test(list(
    statistics = data.frame(
      statistic = stat_names, value = unname(stats), z = unname(z),
      p_value = unname(stats::pnorm(z)), p_boot = p_boot,
      alternative = c("alpha_i < 0 for some i", "alpha_i < 0 for some i",
                      "alpha_i = alpha < 0 for all i", "alpha_i = alpha < 0 for all i")
    ),
    units = units,
    coefficients = coefficients,
    panel = list(
      alpha = res$panel$alpha, se = res$panel$se, t_bar = res$panel$t_bar,
      lags = if (common) as.integer(res$panel$p_bar) else NA_integer_,
      leads = if (common) as.integer(res$panel$q_bar) else NA_integer_
    ),
    boot = boot,
    settings = s,
    N = N, K = K, case = dd + 1L, deterministic = deterministic,
    call = call, formula = formula, y_name = panel$y_name, x_names = panel$x_names,
    index = panel$index
  ), "pcoint_westerlund")
}

# A single order or a range c(min, max), as an integer vector c(min, max).
order_range <- function(x, name) {
  if (!length(x) %in% 1:2 || !is_whole(x) || any(x < 0)) {
    abort("`", name, "` must be a non-negative whole number or a range `c(min, max)`.")
  }
  x <- as.integer(round(x))
  if (length(x) == 1L) c(x, x) else sort(x)
}

# Names of the ECM coefficients in the layout of the Rust engine:
# deterministics, x_{t-1}, Delta y_{t-1..t-p}, Delta x_{t+q..t-p}, y_{t-1}.
ecm_coef_names <- function(dd, y, x, p, q) {
  dx <- function(l) {
    prefix <- if (l == 0L) "D." else if (l > 0L) paste0("L", l, "D.") else paste0("F", -l, "D.")
    paste0(prefix, x)
  }
  c(
    c("constant", "trend")[seq_len(dd)],
    paste0("L.", x),
    if (p > 0L) paste0("L", seq_len(p), "D.", y),
    unlist(lapply(-q:p, dx)),
    paste0("L.", y)
  )
}
