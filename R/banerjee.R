#' Banerjee and Carrion-i-Silvestre (2015) panel cointegration tests with
#' structural breaks and common factors
#'
#' Tests the null hypothesis of no cointegration in a panel whose units are
#' linked by unobserved common factors and may have structural breaks in the
#' level, the trend and the cointegrating vector (Banerjee and
#' Carrion-i-Silvestre 2015). The cross-section dependence is removed by
#' principal components, as in the PANIC approach of Bai and Ng (2004); the
#' panel statistic pools ADF tests of the defactored errors, and the common
#' factors are tested separately.
#'
#' @section Procedure:
#' For each unit, \eqn{\Delta y_{it}}{Delta y_it} is regressed on the first
#' differences of the model's deterministic terms and of the regressors
#' (Section 2.2 of the paper). The residuals of all units form an
#' \eqn{n \times N}{n x N} matrix from which \eqn{r}{r} common factors are
#' estimated by principal components; \eqn{r}{r} is chosen by the IC_p1
#' criterion of Bai and Ng (2002) with at most `max_factors` factors unless
#' it is fixed. IC_p1 is the criterion of Bai and Ng (2004), which the paper
#' cites for its "panel BIC", and it reproduces the factor choices that
#' Banerjee and Carrion-i-Silvestre (2017) report for their production
#' function data. The defactored residuals are cumulated and each unit's ADF
#' t-ratio \eqn{t_i}{t_i} is computed from a regression without
#' deterministic terms. The panel statistic is
#' \deqn{Z = \sqrt{N} (\bar t - \Theta) / \sqrt{\Psi},}{Z = sqrt(N) (tbar - Theta) / sqrt(Psi),}
#' asymptotically standard normal under the null; it is left-tailed and the
#' p-value is \eqn{\Phi(Z)}{pnorm(Z)}. \eqn{\Theta}{Theta} and \eqn{\Psi}{Psi}
#' come from the paper's Remark 1 (Models 1, 2, 4, 5) and Table II
#' (Models 3 and 6, one or two breaks on the 0.1 grid); other break
#' fractions are simulated.
#'
#' A rejection means that the idiosyncratic errors are stationary for some
#' units. The observed variables cointegrate only if, in addition, no common
#' factor is I(1): with one factor this is checked by an ADF test, with
#' several by the sequential MQ tests of Bai and Ng (2004), which estimate
#' the number of common stochastic trends `r1`. If the idiosyncratic test
#' rejects but `r1 > 0`, the variables cointegrate only together with the
#' non-stationary common factors (Section 6.3.2 of the paper).
#'
#' @section Models:
#' | `model` | Deterministic part | Cointegrating vector | Breaks |
#' |---|---|---|---|
#' | 1 | constant, level shifts | stable | none, known or estimated |
#' | 2 | constant, trend, level shifts | stable | none, known or estimated |
#' | 3 | constant, trend, level and trend shifts | stable | known, common |
#' | 4 | constant, level shifts | shifts | known or estimated |
#' | 5 | constant, trend, level shifts | shifts | known or estimated |
#' | 6 | constant, trend, level and trend shifts | shifts | known, common |
#'
#' In Models 4 to 6 the cointegrating vector changes at the same dates as
#' the deterministic terms.
#'
#' @section Breaks:
#' A break date is the **last period of the old regime**: `breaks = 1998`
#' means that the new regime starts in the period after 1998.
#' * `breaks = NULL`: no breaks (Models 1 and 2 only).
#' * A vector of time values: known breaks, common to all units.
#' * A data frame with the unit identifiers in its first column and the
#'   break dates in its second, or a list named by unit: known
#'   unit-specific breaks, possibly a different number per unit (Models 1,
#'   2, 4 and 5).
#' * `breaks = 1` or `2`: the number of breaks to estimate in each unit
#'   (Models 1, 2, 4 and 5). With `break_type = "heterogeneous"` the dates
#'   are unit-specific and estimated by the iterative least-squares
#'   procedure of Section 4.2.1, which alternates between the break dates
#'   and the common factors. With `"homogeneous"` they are common and
#'   chosen to minimise \eqn{Z}{Z} (Section 4.2.2); the infimum statistic is
#'   asymptotically standard normal, but the normal p-value is
#'   anti-conservative in finite samples, so the paper's Table III critical
#'   values are reported too. Every regime spans at least
#'   \eqn{\max(\lfloor T \cdot trim \rfloor, k + 2)}{max(floor(T trim), k + 2)}
#'   periods, `k` being the number of regime-specific slopes. In Model 2
#'   estimated break dates are consistent only for large breaks (Section
#'   6.3.2 of the paper). In Model 5 with a large shift of the cointegrating
#'   vector, the infimum statistic over-rejects severely in simulations of
#'   the paper's own design (Section 5), because the statistic at misplaced
#'   dates is shifted to the left; heterogeneous estimation or known dates
#'   are then preferable (see `inst/validation/RESULTS.md`).
#'
#' @section Endogenous regressors:
#' `dols = c(p, q)` adds the lags \eqn{1..p}{1..p} and leads
#' \eqn{1..q}{1..q} of \eqn{\Delta^2 x_t}{Delta^2 x_t}, and its current
#' value, to the first-difference regressions (dynamic OLS; working-paper
#' version, Appendix A.6), which restores the null distributions when the
#' regressors are not strictly exogenous. The orders are common to all units.
#'
#' @section Assumptions:
#' The panel must be balanced, since the principal components need a common
#' time axis. Theorem 2 of the paper needs the idiosyncratic errors to be
#' independent across units and \eqn{N/T \to 0}{N/T -> 0}; the CD statistic
#' of Pesaran (2015) on the ADF residuals is reported as a check of the
#' first condition.
#'
#' @param formula A formula `y ~ x1 + x2` with the levels of the I(1)
#'   variables. The intercept is ignored: the deterministic terms follow
#'   from `model`.
#' @param data A data frame (including tibbles, data.tables and
#'   `plm::pdata.frame`s) in long format.
#' @param index The names of the unit and time variables, as in
#'   [pcoint_westerlund()].
#' @param model The model of the paper, 1 to 6 (see Models). Required.
#' @param breaks Break specification (see Breaks).
#' @param break_type `"heterogeneous"` or `"homogeneous"` estimated breaks;
#'   only with `breaks = 1` or `2`.
#' @param factors Number of common factors: a criterion of Bai and Ng (2002),
#'   `"ic1"` (default), `"ic2"`, `"ic3"`, `"pc1"`, `"pc2"`, `"pc3"` or
#'   `"bic3"`, or a fixed number (0 for none).
#' @param max_factors Largest number of factors considered by the criterion
#'   (default 6, as in the paper, lowered to \eqn{\min(N, n) - 1}{min(N, n) - 1}
#'   for small panels, `n` being the number of first-difference
#'   observations).
#' @param lags ADF lag order: `NULL` (default) selects it between 0 and the
#'   rule of [banerjee_control()] by `lag_criterion`; a single number fixes
#'   it; `c(min, max)` selects it in that range.
#' @param dols Leads and lags of \eqn{\Delta^2 x}{Delta^2 x}: 0 (default)
#'   for none, a single number for both, or `c(lags, leads)`.
#' @param trim Trimming fraction for estimated breaks.
#' @param control Further settings from [banerjee_control()].
#' @param n_threads Number of threads for the break-date search and the
#'   simulations; the results do not depend on it.
#'
#' @return An object of class `c("pcoint_banerjee", "pcoint_test")`, a list
#'   with:
#'   * `statistics`: one row with the statistic (`Z_c`, `Z_tau`, `Z_gamma`,
#'     or `Zinf_c` / `Zinf_tau` for estimated common breaks), its `value`
#'     \eqn{\bar t}{tbar}, `z`, the normal `p_value`, `p_boot` (`NA`) and
#'     the `alternative`.
#'   * `units`: per unit, the ADF `lags`, number of observations `n`,
#'     `t_stat`, its 5% critical value `cv_5` and `reject_5`, and the break
#'     dates and fractions.
#'   * `factors`: the number of factors \code{r}, the `criterion` and its values, the
#'     `eigenvalues`, the factors `f` (first differences), `F` (cumulated)
#'     and `F_detrended`, the `loadings`; `adf` (one factor) or `mq` (the MQ
#'     sequence); and `r1_c`, `r1_f`, the estimated numbers of I(1) factors
#'     by MQ_c and MQ_f (by the ADF test when \eqn{r = 1}{r = 1}).
#'   * `breaks`: the `type`, the number `m`, the `dates`, the `positions`
#'     (1-based index of the last period of the old regime) and `lambda`;
#'     the `profile` of the search for homogeneous breaks; the iteration
#'     log (`ssr`, `converged` and the number of factors \code{r}) for heterogeneous breaks; and the
#'     Table III critical values for Z*.
#'   * `moments`: `mean`, `var` and their `source`.
#'   * `diagnostics`: the CD statistic and its p-value, the number of
#'     periods it uses, and the share of individual rejections at 5%.
#'   * `settings`, `N`, `K`, `T`, `model`, `call`, `formula`.
#'
#' @references
#' Banerjee, A. and Carrion-i-Silvestre, J. L. (2015). Cointegration in
#' panel data with structural breaks and cross-section dependence.
#' *Journal of Applied Econometrics* 30(1), 1-23. \doi{10.1002/jae.2348}
#'
#' Banerjee, A. and Carrion-i-Silvestre, J. L. (2006). Cointegration in
#' panel data with breaks and cross-section dependence. ECB Working Paper
#' 591.
#'
#' Banerjee, A. and Carrion-i-Silvestre, J. L. (2017). Testing for panel
#' cointegration using common correlated effects estimators. Birmingham
#' Business School Discussion Paper 2017-02.
#'
#' Bai, J. and Ng, S. (2002). Determining the number of factors in
#' approximate factor models. *Econometrica* 70(1), 191-221.
#'
#' Bai, J. and Ng, S. (2004). A PANIC attack on unit roots and
#' cointegration. *Econometrica* 72(4), 1127-1177.
#'
#' Ng, S. and Perron, P. (1995). Unit root tests in ARMA models with
#' data-dependent methods for the selection of the truncation lag. *Journal
#' of the American Statistical Association* 90, 268-281.
#'
#' Pesaran, M. H. (2015). Testing weak cross-sectional dependence in large
#' panels. *Econometric Reviews* 34, 1089-1117.
#'
#' @seealso [banerjee_control()], [summary.pcoint_banerjee()]
#' @export
#' @examples
#' set.seed(1)
#' N <- 15; T <- 60
#' f <- cumsum(rnorm(T))
#' d <- do.call(rbind, lapply(seq_len(N), function(i) {
#'   x <- cumsum(rnorm(T))
#'   y <- 1 + x + rnorm(1) * f + cumsum(rnorm(T))
#'   data.frame(id = i, time = 1950 + seq_len(T), y = y, x = x)
#' }))
#' pcoint_banerjee(y ~ x, data = d, model = 1)
#' pcoint_banerjee(y ~ x, data = d, model = 5, breaks = 1980)
#' \donttest{
#' pcoint_banerjee(y ~ x, data = d, model = 2, breaks = 1, break_type = "homogeneous")
#' }
pcoint_banerjee <- function(formula, data, index = NULL, model, breaks = NULL,
                            break_type = c("heterogeneous", "homogeneous"),
                            factors = "ic1", max_factors = 6L, lags = NULL, dols = 0L,
                            trim = 0.15, control = banerjee_control(),
                            n_threads = getOption("panelcointr.threads", 1L)) {
  call <- match.call()
  max_factors_given <- !missing(max_factors)
  if (missing(model)) {
    abort("`model` is required: one of the six models of Banerjee and Carrion-i-Silvestre ",
          "(2015), 1 to 6 (see ?pcoint_banerjee).")
  }
  if (length(model) != 1L || !is_whole(model) || !model %in% 1:6) {
    abort("`model` must be a whole number between 1 and 6.")
  }
  model <- as.integer(model)
  type_given <- !missing(break_type)
  break_type <- match.arg(break_type)
  if (!inherits(control, "banerjee_control")) {
    abort("`control` must be created by `banerjee_control()`.")
  }
  if (length(n_threads) != 1L || !is_whole(n_threads) || n_threads < 1) {
    abort("`n_threads` must be a single positive whole number.")
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) || trim <= 0 || trim >= 0.5) {
    abort("`trim` must be a number between 0 and 0.5.")
  }
  if (!length(dols) %in% 1:2 || !is_whole(dols) || any(dols < 0)) {
    abort("`dols` must be a non-negative whole number or `c(lags, leads)`.")
  }
  dols <- rep(as.integer(round(dols)), length.out = 2L)
  if (length(max_factors) != 1L || !is_whole(max_factors) || max_factors < 0) {
    abort("`max_factors` must be a non-negative whole number.")
  }
  max_factors <- as.integer(max_factors)
  fac <- resolve_factors(factors)

  panel <- prepare_panel(formula, data, index)
  require_balanced(panel, "pcoint_banerjee()")
  N <- length(panel$T)
  K <- ncol(panel$X)
  T <- panel$T[1L]
  labels <- as.character(panel$ids)
  times <- panel$time[seq_len(T)]
  n <- if (sum(dols) > 0L) T - 2L - sum(dols) else T - 1L
  j <- c("c", "tau", "gamma")[c(1L, 2L, 3L, 1L, 2L, 3L)[model]]

  br <- resolve_breaks(breaks, break_type, type_given, model, panel$ids, times, dols)

  # Lags, factors and the MQ settings.
  k_rule <- if (is.null(control$max_lags)) 4L * as.integer(ceiling((min(N, T) / 100)^(1 / 4))) else control$max_lags
  lag_range <- if (is.null(lags)) c(0L, k_rule) else order_range(lags, "lags")
  if (fac$type == "fixed" && fac$r >= min(N, n)) {
    abort("`factors = ", fac$r, "` must be below min(N, n) = ", min(N, n), ".")
  }
  if (!max_factors_given) max_factors <- min(max_factors, min(N, n) - 1L)
  if (max_factors >= min(N, n)) {
    abort("`max_factors` (", max_factors, ") must be below min(N, n) = ", min(N, n),
          ", where n = ", n, " is the number of rows of the first-difference regressions.")
  }
  mq_bw <- if (is.null(control$mq_bandwidth)) 4L * as.integer(ceiling((min(N, T) / 100)^(1 / 4))) else control$mq_bandwidth
  mq_pmax <- if (is.null(control$mq_max_var)) 4L * as.integer(ceiling((T / 100)^(1 / 4))) else control$mq_max_var

  # Degrees of freedom of the largest regressions.
  m_max <- max(br$m, if (is.matrix(br$pos)) ncol(br$pos) else 0L)
  ncol_fd <- (model %in% c(2L, 3L, 5L, 6L)) + K * (1L + if (sum(dols) > 0L) sum(dols) + 1L else 0L) +
    m_max * (1L + (model %in% c(3L, 6L)) + if (model >= 4L) K else 0L)
  if (n - ncol_fd < control$min_df) {
    abort("Too few periods: the first-difference regressions have ", n, " rows and ", ncol_fd,
          " columns, leaving fewer than ", control$min_df, " residual degree(s) of freedom.")
  }
  adf_df <- n - 1L - lag_range[2L] - (lag_range[2L] + 1L)
  if (adf_df < control$min_df) {
    abort("Too few periods for ADF regressions with ", lag_range[2L], " lags: reduce `lags` ",
          "or `max_lags`.")
  }

  settings <- list(
    model = model, dols_p = dols[1L], dols_q = dols[2L],
    factors = if (fac$type == "fixed") "fixed" else fac$criterion,
    r_fixed = if (fac$type == "fixed") fac$r else 0L, max_factors = max_factors,
    lag_min = lag_range[1L], lag_max = lag_range[2L], lag_criterion = control$lag_criterion,
    tsig_cut = stats::qnorm(1 - control$tsig_level / 2), lag_sample = control$lag_sample,
    pca_scale = control$pca_scale, mq_bandwidth = mq_bw, mq_max_var = mq_pmax,
    iter_max = control$iter_max, iter_tol = control$iter_tol
  )
  y <- panel$y
  X <- as.vector(panel$X)
  tlab <- format(times)
  k_reg <- if (model %in% c(4L, 5L)) K else 0L
  h <- max(floor(trim * T), k_reg + 2L)

  breaks_info <- list(type = br$type, m = br$m)
  pos <- br$pos
  if (br$type == "heterogeneous") {
    est <- banerjee_breaks_rs(y, X, K, panel$offsets, labels, tlab, br$m, h, settings,
                              as.integer(n_threads))
    pos <- matrix(est$pos, N, br$m)
    if (!est$converged) {
      warn("The estimation of the break dates did not converge in ", control$iter_max,
           " iterations; increase `iter_max` in banerjee_control().")
    }
    breaks_info$ssr <- est$ssr
    breaks_info$iterations <- length(est$ssr)
    breaks_info$converged <- est$converged
    breaks_info$r_iteration <- est$r
  } else if (br$type == "homogeneous") {
    prof <- banerjee_profile_rs(y, X, K, panel$offsets, labels, tlab, br$m, h, settings,
                                as.integer(n_threads))
    cand <- matrix(prof$pos, ncol = br$m)
    lam <- (cand + 1) / T
    mom <- banerjee_moments(j, NULL, T, control)
    z <- sqrt(N) * (prof$tbar - mom$mean) / sqrt(mom$var)
    if (all(is.na(z))) abort("The test could not be computed at any admissible break date.")
    best <- which.min(z)
    pos <- matrix(cand[best, ], N, br$m, byrow = TRUE)
    profile <- data.frame(row.names = seq_len(nrow(cand)))
    for (b in seq_len(br$m)) profile[[break_names("break_date", br$m)[b]]] <- times[cand[, b] + 1L]
    for (b in seq_len(br$m)) profile[[break_names("lambda", br$m)[b]]] <- lam[, b]
    profile$t_bar <- prof$tbar
    profile$z <- z
    profile$r <- prof$r
    breaks_info$profile <- profile
  }

  res <- banerjee_rs(y, X, K, panel$offsets, labels, tlab, as.integer(pos), ncol(pos), settings)

  # Break dates and fractions.
  m_unit <- rowSums(pos >= 0L)
  lambda <- ifelse(pos >= 0L, (pos + 1) / T, NA_real_)
  common <- ncol(pos) > 0L && all(apply(pos, 2L, function(v) length(unique(v)) == 1L))
  lambda_common <- if (common) lambda[1L, ] else NULL
  if (ncol(pos)) {
    dates <- data.frame(id = panel$ids)
    for (b in seq_len(ncol(pos))) {
      v <- times[ifelse(pos[, b] >= 0L, pos[, b] + 1L, NA_integer_)]
      dates[[break_names("break_date", ncol(pos))[b]]] <- v
    }
    breaks_info$dates <- dates
    breaks_info$positions <- ifelse(pos >= 0L, pos + 1L, NA_integer_)
    breaks_info$lambda <- lambda
    breaks_info$common <- common
  }

  # Moments, statistic and the unit critical values.
  moments <- banerjee_moments(j, lambda_common, T, control, n_threads)
  tbar <- res$tbar
  zval <- sqrt(N) * (tbar - moments$mean) / sqrt(moments$var)
  inf <- br$type == "homogeneous"
  if (inf && br$m == 1L) breaks_info$table3 <- bcs_table3_row(j, T)
  ref_unit <- null_reference("idio", j, lambda_common, control, n_threads)
  cv5 <- ref_unit$quantile(0.05)
  units <- data.frame(id = panel$ids, lags = res$lags, n = res$nobs, t_stat = res$t,
                      cv_5 = cv5, reject_5 = res$t < cv5)
  if (ncol(pos)) {
    for (b in seq_len(ncol(pos))) {
      units[[break_names("break_date", ncol(pos))[b]]] <- breaks_info$dates[[b + 1L]]
    }
    for (b in seq_len(ncol(pos))) units[[break_names("lambda", ncol(pos))[b]]] <- lambda[, b]
  }

  factors_out <- banerjee_factor_results(res, j, lambda_common, control, fac, max_factors,
                                         times, panel$ids, n_threads)
  stat_name <- paste0(if (inf) "Zinf_" else "Z_", j)

  new_pcoint_test(list(
    statistics = data.frame(statistic = stat_name, value = tbar, z = zval,
                            p_value = stats::pnorm(zval), p_boot = NA_real_,
                            alternative = "rho_i < 1 for some i"),
    units = units,
    factors = factors_out,
    breaks = breaks_info,
    moments = moments,
    diagnostics = c(cd_test(matrix(res$resid, res$n - 1L, N)), list(share_rejected = mean(units$reject_5))),
    settings = c(settings, list(lags = lag_range, dols = dols, trim = trim, h = h,
                                tsig_level = control$tsig_level, moments = control$moments,
                                sim_reps = control$sim_reps, sim_steps = control$sim_steps,
                                mq_level = control$mq_level, n_threads = n_threads)),
    N = N, K = K, T = T, n = res$n, model = model, deterministic = j,
    call = call, formula = formula, y_name = panel$y_name, x_names = panel$x_names,
    index = panel$index
  ), "pcoint_banerjee")
}

break_names <- function(prefix, m) if (m == 1L) prefix else paste0(prefix, "_", seq_len(m))

factor_criteria <- c("bic3", "ic1", "ic2", "ic3", "pc1", "pc2", "pc3")

resolve_factors <- function(factors) {
  if (is.character(factors) && length(factors) == 1L) {
    f <- tolower(factors)
    if (!f %in% factor_criteria) {
      abort("`factors` must be a number or one of ",
            enumerate(paste0("\"", factor_criteria, "\"")), ".")
    }
    return(list(type = "criterion", criterion = f))
  }
  if (length(factors) != 1L || !is_whole(factors) || factors < 0) {
    abort("`factors` must be a non-negative whole number or the name of a criterion.")
  }
  list(type = "fixed", r = as.integer(round(factors)))
}

# Positions (0-based index of the last period of the old regime) of time
# values on the time axis `times`; NA if absent.
match_times <- function(values, times) {
  if (inherits(times, "Date")) {
    return(match(as.numeric(as.Date(values)), as.numeric(times)) - 1L)
  }
  if (inherits(times, "POSIXt")) {
    return(match(as.numeric(as.POSIXct(values)), as.numeric(times)) - 1L)
  }
  if (is.numeric(times) && is.numeric(values)) return(match(values, times) - 1L)
  match(as.character(values), as.character(times)) - 1L
}

# Interpret `breaks`: the type ("none", "known", "heterogeneous",
# "homogeneous"), the number of breaks m and, for known breaks, an N x m
# matrix of 0-based positions padded with -1.
resolve_breaks <- function(breaks, break_type, type_given, model, ids, times, dols) {
  N <- length(ids)
  T <- length(times)
  none <- matrix(integer(), N, 0L)
  if (is.null(breaks)) {
    if (model >= 3L) {
      abort("Model ", model, " needs breaks: give their dates, or the number of breaks to ",
            "estimate (Models 4 and 5).")
    }
    if (type_given) abort("`break_type` applies only to estimated breaks (`breaks = 1` or `2`).")
    return(list(type = "none", m = 0L, pos = none))
  }
  estimated <- is.numeric(breaks) && !is.object(breaks) && length(breaks) == 1L &&
    is_whole(breaks) && breaks %in% 1:2
  if (estimated) {
    if (model %in% c(3L, 6L)) {
      abort("Models 3 and 6 allow known common breaks only (the paper does not cover estimated ",
            "trend breaks): give the break dates.")
    }
    return(list(type = break_type, m = as.integer(breaks), pos = none))
  }
  if (type_given) {
    abort("`break_type` applies only to estimated breaks (`breaks = 1` or `2`); ",
          "known break dates are used as given.")
  }
  if (is.numeric(breaks) && !is.object(breaks) && length(breaks) == 1L && is_whole(breaks) &&
      breaks > 2 && !breaks %in% times) {
    abort("Only 1 or 2 breaks can be estimated; `breaks = ", breaks, "` is not a time value of ",
          "the panel either.")
  }
  # Known breaks: common vector, data frame (id, date) or list by unit.
  unit_list <- NULL
  if (is.data.frame(breaks)) {
    if (ncol(breaks) < 2L) abort("A data frame of breaks needs an id and a date column.")
    unit_list <- split(breaks[[2L]], factor(as.character(breaks[[1L]]), levels = unique(as.character(breaks[[1L]]))))
  } else if (is.list(breaks)) {
    if (is.null(names(breaks)) || any(names(breaks) == "")) abort("A list of breaks must be named by unit.")
    unit_list <- breaks
  }
  pos_of <- function(v, who) {
    if (!length(v)) return(integer())
    p <- match_times(v, times)
    if (anyNA(p)) {
      abort("Break date", if (sum(is.na(p)) > 1L) "s", " ", enumerate(format(v[is.na(p)])),
            if (!is.null(who)) paste0(" of unit ", who), " not among the time values of the panel.")
    }
    p <- sort(p)
    if (anyDuplicated(p)) abort("Duplicated break dates", if (!is.null(who)) paste0(" for unit ", who), ".")
    lo <- if (sum(dols) > 0L) 1L + dols[1L] else 0L
    hi <- if (sum(dols) > 0L) T - 2L - dols[2L] else T - 2L
    bad <- p < lo | p > hi
    if (any(bad)) {
      abort("Break date", if (sum(bad) > 1L) "s", " ", enumerate(format(times[p[bad] + 1L])),
            if (!is.null(who)) paste0(" of unit ", who),
            " leave", if (sum(bad) == 1L) "s", " no period of the new regime in the sample of the ",
            "first-difference regressions; a break is the last period of the old regime.")
    }
    as.integer(p)
  }
  if (is.null(unit_list)) {
    p <- pos_of(breaks, NULL)
    return(list(type = "known", m = length(p), pos = matrix(p, N, length(p), byrow = TRUE)))
  }
  if (model %in% c(3L, 6L)) {
    abort("Models 3 and 6 need breaks common to all units: give a vector of dates.")
  }
  key <- as.character(ids)
  missing_units <- setdiff(key, names(unit_list))
  if (length(missing_units)) {
    abort("Breaks are missing for unit", if (length(missing_units) > 1L) "s", " ",
          enumerate(missing_units), " (use an empty vector for a unit without breaks).")
  }
  extra <- setdiff(names(unit_list), key)
  if (length(extra)) abort("Breaks given for unknown unit", if (length(extra) > 1L) "s", " ", enumerate(extra), ".")
  plist <- lapply(key, function(k) pos_of(unit_list[[k]], k))
  m <- max(lengths(plist))
  pos <- matrix(-1L, N, m)
  for (i in seq_len(N)) if (length(plist[[i]])) pos[i, seq_along(plist[[i]])] <- plist[[i]]
  list(type = "known", m = m, pos = pos)
}

on_grid <- function(lambda) {
  length(lambda) && all(abs(lambda * 10 - round(lambda * 10)) < 1e-7) &&
    all(round(lambda * 10) %in% 1:9)
}

# Moments of the individual statistics: Remark 1, Table II or WP Table 9,
# otherwise simulated at the exact break fractions.
banerjee_moments <- function(j, lambda, T, control, n_threads = 1L) {
  finite <- control$moments == "finite_sample"
  interp <- function(v) {
    Tc <- min(max(T, 50), 1000)
    stats::approx(1 / bcs_wp_table9$T, v, xout = 1 / Tc)$y
  }
  if (j != "gamma") {
    if (finite) {
      tab <- bcs_wp_table9[[j]]
      return(list(mean = interp(tab["mean", ]), var = interp(tab["var", ]),
                  source = "WP Table 9"))
    }
    m <- bcs_remark1[[j]]
    return(list(mean = m[["mean"]], var = m[["var"]], source = "Remark 1"))
  }
  if (is.null(lambda)) return(list(mean = NA_real_, var = NA_real_, source = "none"))
  idx <- round(lambda * 10)
  if (finite) {
    if (length(lambda) == 1L && on_grid(lambda)) {
      return(list(mean = interp(bcs_wp_table9$gamma_mean[, idx]),
                  var = interp(bcs_wp_table9$gamma_var[, idx]), source = "WP Table 9"))
    }
    warn("Finite-sample moments are tabulated only for one break on the 0.1 grid ",
         "(WP Table 9); using asymptotic moments.")
  }
  if (on_grid(lambda) && length(lambda) == 1L) {
    return(list(mean = unname(bcs_table2_one["mean", idx]), var = unname(bcs_table2_one["var", idx]),
                source = "Table II"))
  }
  if (on_grid(lambda) && length(lambda) == 2L) {
    return(list(mean = bcs_table2_two[idx[1L], idx[2L]], var = bcs_table2_two[idx[2L], idx[1L]],
                source = "Table II"))
  }
  d <- bcs_null("idio", "gamma", lambda, 1L, control, n_threads)
  list(mean = mean(d), var = stats::var(as.vector(d)), source = "simulated")
}

bcs_table3_row <- function(j, T) {
  i <- which.min(abs(bcs_table3$T - T))
  stats::setNames(bcs_table3[[j]][i, ], paste0(100 * bcs_table3$level, "%"))
}

# Seed of the simulations of null distributions: they are deterministic
# functions of their inputs and leave R's random number generator alone.
bcs_sim_seed <- 20150101
bcs_null_cache <- new.env(parent = emptyenv())

bcs_null <- function(kind, det, lambda, q, control, n_threads = 1L) {
  steps <- control$sim_steps
  reps <- control$sim_reps
  key <- paste(kind, det, paste(sprintf("%.12f", lambda), collapse = ","), q, steps, reps, sep = "|")
  hit <- bcs_null_cache[[key]]
  if (!is.null(hit)) return(hit)
  v <- matrix(banerjee_null_rs(kind, det, as.numeric(lambda), steps, reps, as.integer(q),
                               bcs_sim_seed, as.integer(n_threads)), nrow = reps)
  assign(key, v, envir = bcs_null_cache)
  v
}

# Reference distribution of an individual statistic: the quantile tables
# in sysdata for c and tau, simulated draws otherwise. Returns functions
# `quantile(p)` and `p_value(x)`.
null_reference <- function(kind, j, lambda, control, n_threads = 1L) {
  name <- paste0(kind, "_", j)
  qt <- get0("bcs_quantiles", envir = topenv(environment(null_reference)), inherits = FALSE)
  tab <- qt[[name]]
  if (j != "gamma" && !is.null(tab)) {
    probs <- qt$probs
    return(list(
      quantile = function(p) stats::approx(probs, tab, xout = p)$y,
      p_value = function(x) stats::approx(tab, probs, xout = x, rule = 2)$y,
      source = "simulated table"
    ))
  }
  if (j == "gamma" && is.null(lambda)) {
    return(list(quantile = function(p) NA_real_, p_value = function(x) NA_real_, source = "none"))
  }
  d <- as.vector(bcs_null(if (kind == "idio") "idio" else "factor_df", j,
                          if (j == "gamma") lambda else numeric(), 1L, control, n_threads))
  list(
    quantile = function(p) unname(stats::quantile(d, p, names = FALSE)),
    p_value = function(x) mean(d <= x),
    source = "simulated"
  )
}

mq_critical <- function(j, lambda, q, level, control, n_threads = 1L) {
  lev <- format(level)
  out <- rep(NA_real_, length(q))
  src <- "Bai and Ng (2004) Table I"
  if (j != "gamma") {
    ok <- q <= 6L
    out[ok] <- bai_ng_table1[[j]][q[ok], lev]
  } else if (length(lambda) == 1L && on_grid(lambda)) {
    ok <- q <= 6L
    out[ok] <- bcs_table1[q[ok], lev, round(lambda * 10)]
    src <- "BCS Table I"
  } else {
    ok <- rep(FALSE, length(q))
  }
  if (any(!ok)) {
    d <- bcs_null("mq", j, if (j == "gamma") lambda else numeric(), max(q), control, n_threads)
    for (i in which(!ok)) out[i] <- stats::quantile(d[, q[i]], level, names = FALSE)
    src <- "simulated"
  }
  list(cv = out, source = src)
}

banerjee_factor_results <- function(res, j, lambda, control, fac, max_factors, times, ids,
                                    n_threads) {
  r <- res$r
  n <- res$n
  rows <- times[res$r0 + seq_len(n)]
  crit <- stats::setNames(res$criteria, 0:max_factors)
  out <- list(r = r, criterion = if (fac$type == "fixed") "fixed" else fac$criterion,
              criteria = crit, eigenvalues = res$eigenvalues)
  if (r > 0L) {
    f <- matrix(res$f, n, r, dimnames = list(format(rows), paste0("F", seq_len(r))))
    out$f <- f
    out$F <- apply(f, 2L, cumsum)
    if (!is.matrix(out$F)) out$F <- matrix(out$F, n, r, dimnames = dimnames(f))
    out$F_detrended <- matrix(res$fd, n, r, dimnames = dimnames(f))
    out$loadings <- matrix(res$loadings, length(ids), r,
                           dimnames = list(as.character(ids), paste0("F", seq_len(r))))
  }
  level <- control$mq_level
  if (r == 0L) {
    out$r1_c <- out$r1_f <- 0L
  } else if (r == 1L) {
    ref <- null_reference("factor", j, lambda, control, n_threads)
    t <- res$factor_adf$t
    cv <- vapply(c(0.01, 0.05, 0.10), ref$quantile, numeric(1))
    out$adf <- list(statistic = t, lags = res$factor_adf$lags, n = res$factor_adf$n,
                    critical = stats::setNames(cv, c("1%", "5%", "10%")),
                    p_value = ref$p_value(t), source = ref$source)
    out$r1_c <- out$r1_f <- if (isTRUE(ref$p_value(t) < level)) 0L else 1L
  } else {
    mq <- res$mq
    crit_v <- mq_critical(j, lambda, mq$q, level, control, n_threads)
    tab <- data.frame(q = mq$q, MQ_c = mq$mq_c, MQ_f = mq$mq_f, cv = crit_v$cv,
                      var_order = mq$var_order)
    sequential <- function(stat) {
      for (i in seq_along(tab$q)) if (!isTRUE(stat[i] < tab$cv[i])) return(tab$q[i])
      0L
    }
    out$mq <- tab
    out$mq_source <- crit_v$source
    out$r1_c <- as.integer(sequential(tab$MQ_c))
    out$r1_f <- as.integer(sequential(tab$MQ_f))
  }
  out
}

# Pesaran's CD statistic on the ADF residuals (an n x N matrix padded with
# NaN), over the periods where every unit has a residual.
cd_test <- function(resid) {
  N <- ncol(resid)
  ok <- stats::complete.cases(resid)
  Tc <- sum(ok)
  if (N < 2L || Tc < 3L) return(list(cd = NA_real_, cd_p_value = NA_real_, cd_periods = Tc))
  rho <- stats::cor(resid[ok, , drop = FALSE])
  cd <- sqrt(2 * Tc / (N * (N - 1))) * sum(rho[upper.tri(rho)])
  list(cd = cd, cd_p_value = 2 * stats::pnorm(-abs(cd)), cd_periods = Tc)
}
