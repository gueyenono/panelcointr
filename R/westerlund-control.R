#' Control settings for the Westerlund (2007) tests
#'
#' Fine-grained choices for [pcoint_westerlund()]. Each argument left at
#' `NULL` takes its value from the `compat` preset of [pcoint_westerlund()]
#' (see Details); a supplied value overrides the preset. The settings
#' actually used are stored in the `settings` element of the result.
#'
#' Defaults follow Westerlund (2007). `compat = "xtwest"` reproduces the
#' Stata command `xtwest` (Persyn and Westerlund 2008) and
#' `compat = "xtwest_westerlund"` its `westerlund` option:
#'
#' | Setting | `"none"` | `"xtwest"` | `"xtwest_westerlund"` |
#' |---|---|---|---|
#' | `t_norm` | `"effective"` | `"df"` | `"effective"` |
#' | `panel_lags` | `"individual"` | `"mean"` | `"mean"` |
#' | `ic_sample` | `"common"` | `"varying"` | `"varying"` |
#' | `lrv_deterministics` | `"regression"` | `"implied"` | `"none"` |
#' | `demean_u` | `TRUE` | `FALSE` | `FALSE` |
#'
#' The compat presets also switch settings that have no control argument,
#' to match `xtwest` exactly: ties in the information criterion go to the
#' larger model; the composite error of `xtwest` also covers the last `q`
#' periods; `xtwest` averages Ga over rows (weighting units by `T_i`) and
#' gives every unit the degrees of freedom of the mean `T_i` in Pt and Pa;
#' `"xtwest_westerlund"` uses `RSS/n` in the standard errors, restricts the
#' long-run variance of \eqn{\Delta y}{Delta y} to the estimation sample, selects lags with
#' its own criterion and uses legacy asymptotic moments. The bootstrap
#' controls are not part of the presets: with any `compat`, the bootstrap
#' is the algorithm described in [pcoint_westerlund()].
#'
#' @param t_norm Sample size \eqn{\tilde T_i}{T_i} that normalises Ga (and, averaged, Pa):
#'   `"effective"`, the number of observations in the unit's regression
#'   (the paper's recommendation); `"df"`, its residual degrees of freedom
#'   (`xtwest`); or `"T"`, the number of time periods.
#' @param panel_lags Lag and lead orders of the panel statistics:
#'   `"individual"`, each unit's own orders (the paper); or `"mean"`, the
#'   floor of the average selected orders for every unit (`xtwest`).
#' @param ic_sample Sample on which the information criterion compares
#'   lag and lead orders: `"common"`, the sample of the largest model for
#'   every candidate; or `"varying"`, each candidate's own sample
#'   (`xtwest`), which makes the choice depend on the units of `y`.
#' @param lrv_deterministics Deterministic terms removed from \eqn{\Delta y}{Delta y} before
#'   its long-run variance: `"regression"`, the residuals of \eqn{\Delta y}{Delta y} on the
#'   deterministic terms (the paper); `"implied"`, demeaning in the trend
#'   case only (`xtwest`); or `"none"`.
#' @param demean_u Demean the composite error before its long-run variance
#'   when the model has a constant. This makes the trend case invariant to
#'   linear trends in the regressors; `FALSE` is the paper's literal
#'   procedure.
#' @param boot_reselect Re-select lags and leads in every bootstrap
#'   replicate (default `TRUE`); `FALSE` keeps each unit's selected orders,
#'   which is faster.
#' @param boot_burn Number of initial bootstrap periods discarded (default 50).
#' @param boot_pvalue Bootstrap p-value: `"plus_one"`,
#'   `(1 + #{S* <= S}) / (B + 1)` (default); or `"xtwest"`, `#{S* <= S} / B`.
#' @param boot_chunk Replications per call into the compiled code; R checks
#'   for interrupts and updates progress between calls.
#' @param min_boot_periods Minimum number of calendar periods in which
#'   every unit has a first-stage residual, needed for the bootstrap.
#' @param min_df Minimum residual degrees of freedom of every unit's
#'   largest regression.
#' @param keep_boot Keep the matrix of bootstrap statistics in the result.
#'
#' @return An object of class `westerlund_control`: a list of the supplied
#'   settings.
#' @seealso [pcoint_westerlund()]
#' @export
#' @examples
#' westerlund_control(t_norm = "df", boot_burn = 100)
westerlund_control <- function(t_norm = NULL, panel_lags = NULL, ic_sample = NULL,
                               lrv_deterministics = NULL, demean_u = NULL,
                               boot_reselect = NULL, boot_burn = NULL, boot_pvalue = NULL,
                               boot_chunk = 100L, min_boot_periods = 10L, min_df = 1L,
                               keep_boot = TRUE) {
  choice <- function(x, name, options) {
    if (is.null(x)) return(NULL)
    if (!is.character(x) || length(x) != 1L || !x %in% options) {
      abort("`", name, "` must be one of ", enumerate(paste0("\"", options, "\"")), ".")
    }
    x
  }
  flag <- function(x, name) {
    if (is.null(x)) return(NULL)
    if (!is.logical(x) || length(x) != 1L || is.na(x)) abort("`", name, "` must be TRUE or FALSE.")
    x
  }
  count <- function(x, name, min) {
    if (is.null(x)) return(NULL)
    if (length(x) != 1L || !is_whole(x) || x < min) {
      abort("`", name, "` must be a whole number of at least ", min, ".")
    }
    as.integer(x)
  }
  out <- list(
    t_norm = choice(t_norm, "t_norm", c("effective", "df", "T")),
    panel_lags = choice(panel_lags, "panel_lags", c("individual", "mean")),
    ic_sample = choice(ic_sample, "ic_sample", c("common", "varying")),
    lrv_deterministics = choice(lrv_deterministics, "lrv_deterministics",
                                c("regression", "implied", "none")),
    demean_u = flag(demean_u, "demean_u"),
    boot_reselect = flag(boot_reselect, "boot_reselect"),
    boot_burn = count(boot_burn, "boot_burn", 0L),
    boot_pvalue = choice(boot_pvalue, "boot_pvalue", c("plus_one", "xtwest")),
    boot_chunk = count(boot_chunk, "boot_chunk", 1L),
    min_boot_periods = count(min_boot_periods, "min_boot_periods", 2L),
    min_df = count(min_df, "min_df", 1L),
    keep_boot = flag(keep_boot, "keep_boot")
  )
  structure(out[!vapply(out, is.null, logical(1))], class = "westerlund_control")
}

#' @export
print.westerlund_control <- function(x, ...) {
  cat("<westerlund_control>\n")
  if (!length(x)) {
    cat("  (all defaults)\n")
  } else {
    for (n in names(x)) cat("  ", n, ": ", format(x[[n]]), "\n", sep = "")
  }
  invisible(x)
}

# Presets: the user-visible settings first, then the internal ones that
# only the compat modes change.
westerlund_presets <- list(
  none = list(
    t_norm = "effective", panel_lags = "individual", ic_sample = "common",
    lrv_deterministics = "regression", demean_u = TRUE,
    ic_ties = "smaller", lrv_y_sample = "all", u_sample = "estimation",
    se_alpha = "df", sigma_df = "own", ga_weights = "units", tbar = "units",
    moments = "table1", lrwindow = NULL
  ),
  xtwest = list(
    t_norm = "df", panel_lags = "mean", ic_sample = "varying",
    lrv_deterministics = "implied", demean_u = FALSE,
    ic_ties = "larger", lrv_y_sample = "all", u_sample = "extended",
    se_alpha = "df", sigma_df = "tbar", ga_weights = "rows", tbar = "rows",
    moments = "table1", lrwindow = 2L
  ),
  xtwest_westerlund = list(
    t_norm = "effective", panel_lags = "mean", ic_sample = "varying",
    lrv_deterministics = "none", demean_u = FALSE,
    ic_ties = "larger", lrv_y_sample = "estimation", u_sample = "estimation",
    se_alpha = "n", sigma_df = "tbar", ga_weights = "rows", tbar = "rows",
    moments = "legacy", lrwindow = 2L
  )
)

westerlund_boot_defaults <- list(
  boot_reselect = TRUE, boot_burn = 50L, boot_pvalue = "plus_one",
  boot_chunk = 100L, min_boot_periods = 10L, min_df = 1L, keep_boot = TRUE
)

# Resolve explicit control values, then the compat preset, then defaults.
resolve_westerlund_settings <- function(compat, control, deterministic, criterion) {
  s <- utils::modifyList(c(westerlund_presets[[compat]], westerlund_boot_defaults),
                         unclass(control))
  s$compat <- compat
  s$deterministic <- deterministic
  s$criterion <- if (compat == "xtwest_westerlund") "westerlund" else criterion
  s
}

# The settings passed to the Rust engine.
westerlund_engine_settings <- function(s) {
  s[c("deterministic", "criterion", "ic_sample", "ic_ties", "lrv_deterministics",
      "lrv_y_sample", "u_sample", "demean_u", "t_norm", "panel_lags", "se_alpha",
      "sigma_df", "ga_weights", "tbar")]
}
