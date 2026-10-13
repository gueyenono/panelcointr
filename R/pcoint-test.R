# The `pcoint_test` class shared by the tests of the package.
#
# Every test returns a list of class c("pcoint_<test>", "pcoint_test") with
# at least:
#   statistics  data.frame with columns statistic, value, z, p_value,
#               p_boot and alternative;
#   N, K, call, formula, settings.

new_pcoint_test <- function(x, subclass) {
  stopifnot(is.list(x), is.data.frame(x$statistics))
  structure(x, class = c(subclass, "pcoint_test"))
}

#' @export
as.data.frame.pcoint_test <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- x$statistics
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

# Table of statistics in the layout of xtwest: Value, Z-value, P-value and,
# with a bootstrap, the robust (bootstrap) P-value.
format_statistics <- function(stats, digits = 3L) {
  num <- function(v, d) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
  tab <- data.frame(
    Statistic = stats$statistic,
    Value = num(stats$value, digits),
    `Z-value` = num(stats$z, digits),
    `P-value` = num(stats$p_value, digits),
    check.names = FALSE
  )
  if (any(!is.na(stats$p_boot))) tab$`Robust P` <- num(stats$p_boot, digits)
  tab
}

#' @export
print.pcoint_test <- function(x, ...) {
  print(format_statistics(x$statistics), row.names = FALSE, right = TRUE)
  invisible(x)
}
