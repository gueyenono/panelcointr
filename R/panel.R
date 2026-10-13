# Shared input handling for the panel cointegration tests.

abort <- function(..., class = NULL) {
  msg <- paste0(...)
  stop(structure(
    class = c(class, "panelcointr_error", "error", "condition"),
    list(message = msg, call = NULL)
  ))
}

warn <- function(...) warning(paste0(...), call. = FALSE)

# Comma-separated list of labels, truncated after `max` entries.
enumerate <- function(x, max = 10L) {
  x <- as.character(x)
  if (length(x) > max) {
    paste0(paste(x[seq_len(max)], collapse = ", "), ", ... (", length(x) - max, " more)")
  } else {
    paste(x, collapse = ", ")
  }
}

is_whole <- function(x) is.numeric(x) && all(is.finite(x)) && all(abs(x - round(x)) < 1e-8)

# Map a time variable to integer positions. Numeric time must step by one
# within each unit and is used as is; dates, date-times, factors and
# character time are mapped to their rank among the time values present.
time_positions <- function(time, time_name) {
  if (inherits(time, c("Date", "POSIXt"))) {
    v <- as.numeric(time)
    return(match(v, sort(unique(v))))
  }
  if (is.factor(time)) {
    v <- as.integer(time)
    return(match(v, sort(unique(v))))
  }
  if (is.character(time)) {
    return(match(time, sort(unique(time))))
  }
  if (is.numeric(time)) {
    if (!is_whole(time)) {
      abort("The time variable `", time_name, "` must hold whole numbers ",
            "(or be a Date, factor or character variable).")
    }
    return(as.integer(round(time)))
  }
  abort("The time variable `", time_name, "` must be numeric, a Date, a factor or character.")
}

# Validate and stack a panel for the Rust engine.
#
# Returns a list with `y` (stacked over units), `X` (matrix, `sum(T)` x K),
# `offsets` (length N + 1, 0-based), `ids` (unit identifiers in panel
# order), `T`, `start` (0-based calendar position of each unit's first
# period), `time` (original time values of the kept rows), and the names
# of the variables.
prepare_panel <- function(formula, data, index = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    abort("`formula` must be a two-sided formula such as `y ~ x1 + x2`.")
  }
  if (!is.data.frame(data)) {
    abort("`data` must be a data frame.")
  }

  # Panel index.
  id <- time <- NULL
  if (is.null(index)) {
    if (inherits(data, "pdata.frame")) {
      idx <- attr(data, "index")
      id <- idx[[1L]]
      time <- idx[[2L]]
      index <- names(idx)[1:2]
    } else {
      if (ncol(data) < 2L) {
        abort("`data` needs an id and a time column; supply `index = c(\"id\", \"time\")`.")
      }
      index <- names(data)[1:2]
    }
  } else if (!is.character(index) || length(index) != 2L || anyNA(index)) {
    abort("`index` must name the id and time columns, e.g. `index = c(\"id\", \"time\")`.")
  }
  df <- as.data.frame(data)
  if (is.null(id)) {
    missing_cols <- setdiff(index, names(df))
    if (length(missing_cols)) {
      abort("Index column", if (length(missing_cols) > 1L) "s", " not found in `data`: ",
            enumerate(missing_cols), ".")
    }
    id <- df[[index[1L]]]
    time <- df[[index[2L]]]
  }
  id <- as.vector(id)
  if (is.factor(time)) time <- droplevels(time)
  if (anyNA(id)) abort("The id variable `", index[1L], "` has missing values.")
  if (anyNA(time)) abort("The time variable `", index[2L], "` has missing values.")

  # Variables.
  others <- df[setdiff(names(df), index)]
  tt <- stats::terms(formula, data = others)
  if (attr(tt, "response") != 1L) abort("`formula` must have a response.")
  mf <- stats::model.frame(tt, data = df, na.action = stats::na.pass)
  y_name <- deparse(formula[[2L]], width.cutoff = 500L)
  x_names <- attr(tt, "term.labels")
  if (!length(x_names)) {
    abort("`formula` needs at least one regressor on the right-hand side.")
  }
  y <- mf[[1L]]
  if (!is.numeric(y) || !is.null(dim(y))) {
    abort("The response `", y_name, "` must be a numeric vector.")
  }
  X <- vapply(x_names, function(v) {
    col <- mf[[v]]
    if (is.null(col) || !is.numeric(col) || !is.null(dim(col))) {
      abort("The regressor `", v, "` must be a numeric vector ",
            "(factors, matrices and interactions are not supported).")
    }
    as.numeric(col)
  }, numeric(nrow(mf)))
  X <- matrix(X, nrow = nrow(mf), dimnames = list(NULL, x_names))
  y <- as.numeric(y)
  if (any(is.infinite(y)) || any(is.infinite(X))) {
    abort("The data contain infinite values.")
  }

  # Sort by unit and time.
  tpos <- time_positions(time, index[2L])
  ord <- order(id, tpos)
  id <- id[ord]
  tpos <- tpos[ord]
  time <- time[ord]
  y <- y[ord]
  X <- X[ord, , drop = FALSE]
  # Sorted by unit and time, so duplicated periods are adjacent.
  n_obs <- length(id)
  dup <- c(FALSE, id[-1L] == id[-n_obs] & tpos[-1L] == tpos[-n_obs])
  if (any(dup)) {
    abort("Duplicated time periods within unit", if (length(unique(id[dup])) > 1L) "s", ": ",
          enumerate(unique(id[dup])), ".")
  }

  # Trim leading and trailing incomplete rows; reject interior gaps.
  complete <- !is.na(y) & rowSums(is.na(X)) == 0L
  units <- unique(id)
  unit_of <- match(id, units)
  first <- tapply(seq_along(id), unit_of, function(r) {
    r <- r[complete[r]]
    if (length(r)) min(r) else NA_integer_
  })
  last <- tapply(seq_along(id), unit_of, function(r) {
    r <- r[complete[r]]
    if (length(r)) max(r) else NA_integer_
  })
  empty <- is.na(first)
  if (all(empty)) abort("No unit has a complete observation.")
  if (any(empty)) {
    warn("Dropping unit", if (sum(empty) > 1L) "s", " without complete observations: ",
         enumerate(units[empty]), ".")
  }
  keep <- logical(length(id))
  for (u in which(!empty)) keep[first[u]:last[u]] <- TRUE
  interior <- keep & !complete
  if (any(interior)) {
    bad <- unique(id[interior])
    abort("Missing values inside the sample of unit", if (length(bad) > 1L) "s", " ",
          enumerate(bad), "; only leading and trailing missing values can be dropped.")
  }
  id <- id[keep]
  tpos <- tpos[keep]
  time <- time[keep]
  y <- y[keep]
  X <- X[keep, , drop = FALSE]
  units <- units[!empty]
  unit_of <- match(id, units)
  step <- c(1L, diff(tpos))
  new_unit <- c(TRUE, unit_of[-1L] != unit_of[-length(unit_of)])
  gap <- !new_unit & step != 1L
  if (any(gap)) {
    bad <- unique(id[gap])
    abort("Gaps in the time variable `", index[2L], "` within unit",
          if (length(bad) > 1L) "s", " ", enumerate(bad), ".")
  }

  T_i <- tabulate(unit_of, nbins = length(units))
  first_pos <- tpos[new_unit]
  list(
    y = y,
    X = X,
    offsets = c(0L, cumsum(T_i)),
    ids = units,
    T = T_i,
    start = as.integer(first_pos - min(first_pos)),
    time = time,
    index = index,
    y_name = y_name,
    x_names = x_names
  )
}

# Stop unless every unit of a prepared panel covers the same calendar
# window; the error lists each unit's window.
require_balanced <- function(panel, test) {
  ends <- panel$start + panel$T
  if (length(unique(panel$start)) == 1L && length(unique(ends)) == 1L) return(invisible(panel))
  first <- panel$time[panel$offsets[-length(panel$offsets)] + 1L]
  last <- panel$time[panel$offsets[-1L]]
  windows <- paste0(panel$ids, " (", format(first), " to ", format(last), ")")
  abort(test, " needs a balanced panel: after dropping leading and trailing missing ",
        "values, every unit must cover the same periods. Windows: ", enumerate(windows), ".")
}
