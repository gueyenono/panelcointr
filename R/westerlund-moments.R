# Asymptotic moments of the Westerlund (2007) statistics.
#
# Westerlund, J. (2007), Table 1, "Simulated asymptotic moments" (Table 1 of
# the working paper METEOR RM/06/056, which xtwest.ado v1.5 hard-codes
# identically). Rows: Case 1 (no deterministics), Case 2 (constant), Case 3
# (constant and trend); columns: K = 1, ..., 6 regressors.
#   Gt: mean tilde-Theta_2,           variance tilde-Sigma_22
#   Ga: mean tilde-Theta_1,           variance tilde-Sigma_11
#   Pt: mean Theta_2 Theta_1^(-1/2),  variance varphi' Sigma varphi
#   Pa: mean Theta_2 Theta_1^(-1),    variance phi' Sigma phi
westerlund_table1 <- list(
  mean = list(
    Gt = rbind(c(-0.9763, -1.3816, -1.7093, -1.9789, -2.1985, -2.4262),
               c(-1.7776, -2.0349, -2.2332, -2.4453, -2.6462, -2.8358),
               c(-2.3664, -2.5284, -2.7040, -2.8639, -3.0146, -3.1710)),
    Ga = rbind(c(-3.8022, -5.8239, -7.8108, -9.8791, -11.7239, -13.8581),
               c(-7.1423, -9.1249, -10.9667, -12.9561, -14.9752, -17.0673),
               c(-12.0116, -13.6324, -15.5262, -17.3648, -19.2533, -21.2479)),
    Pt = rbind(c(-0.5105, -0.9370, -1.3169, -1.6167, -1.8815, -2.1256),
               c(-1.4476, -1.7131, -1.9206, -2.1484, -2.3730, -2.5765),
               c(-2.1124, -2.2876, -2.4633, -2.6275, -2.7858, -2.9537)),
    Pa = rbind(c(-1.0263, -2.4988, -4.2699, -6.1141, -8.0317, -10.0074),
               c(-4.2303, -5.8650, -7.4599, -9.3057, -11.3152, -13.3180),
               c(-8.9326, -10.4874, -12.1672, -13.8889, -15.6815, -17.6515))
  ),
  var = list(
    Gt = rbind(c(1.0823, 1.0981, 1.0489, 1.0576, 1.0351, 1.0409),
               c(0.8071, 0.8481, 0.8886, 0.9119, 0.9083, 0.9236),
               c(0.6603, 0.7070, 0.7586, 0.8228, 0.8477, 0.8599)),
    Ga = rbind(c(20.6868, 29.9016, 39.0109, 50.5741, 58.9595, 69.5967),
               c(29.6336, 39.3428, 49.4880, 58.7035, 67.9499, 79.1093),
               c(46.2420, 53.7428, 64.5591, 74.7403, 84.7990, 94.0024)),
    Pt = rbind(c(1.3624, 1.7657, 1.7177, 1.6051, 1.4935, 1.4244),
               c(0.9885, 1.0663, 1.1168, 1.1735, 1.1684, 1.1589),
               c(0.7649, 0.8137, 0.8857, 0.9985, 0.9918, 0.9898)),
    Pa = rbind(c(8.3827, 24.0223, 39.8827, 53.4518, 63.2406, 76.6757),
               c(19.7090, 31.2637, 42.9975, 57.4844, 69.4374, 81.0384),
               c(37.5948, 45.6890, 57.9985, 74.1258, 81.3934, 91.2392))
  )
)

# Moments used by xtwest.ado (lines 1201-1213) under its `westerlund`
# option, for K = 1 only, presumably from Westerlund's original GAUSS code.
# They differ from Table 1.
westerlund_legacy_moments <- list(
  constant = list(mean = c(Gt = -1.793, Ga = -7.2014, Pt = -1.4746, Pa = -4.3559),
                  var = c(Gt = 0.7904, Ga = 29.3677, Pt = 1.0262, Pa = 21.0535)),
  trend = list(mean = c(Gt = -2.356, Ga = -11.8978, Pt = -2.1128, Pa = -8.9536),
               var = c(Gt = 0.6450, Ga = 44.2471, Pt = 0.7371, Pa = 35.6802))
)

# Moments for a deterministic case and K regressors, or NULL if none are
# tabulated.
westerlund_moments <- function(deterministic, K, source = c("table1", "legacy")) {
  source <- match.arg(source)
  stats <- c("Gt", "Ga", "Pt", "Pa")
  if (source == "legacy") {
    if (K != 1L || deterministic == "none") return(NULL)
    return(westerlund_legacy_moments[[deterministic]])
  }
  if (K < 1L || K > 6L) return(NULL)
  case <- match(deterministic, c("none", "constant", "trend"))
  list(
    mean = vapply(stats, function(s) westerlund_table1$mean[[s]][case, K], numeric(1)),
    var = vapply(stats, function(s) westerlund_table1$var[[s]][case, K], numeric(1))
  )
}

# Standardised statistics (Westerlund 2007, Theorem 1). Pt is not
# multiplied by sqrt(N).
westerlund_z <- function(stats, N, moments) {
  if (is.null(moments)) return(stats::setNames(rep(NA_real_, 4L), names(stats)))
  m <- moments$mean
  v <- moments$var
  rn <- sqrt(N)
  c(
    Gt = rn * (stats[["Gt"]] - m[["Gt"]]) / sqrt(v[["Gt"]]),
    Ga = rn * (stats[["Ga"]] - m[["Ga"]]) / sqrt(v[["Ga"]]),
    Pt = (stats[["Pt"]] - rn * m[["Pt"]]) / sqrt(v[["Pt"]]),
    Pa = rn * (stats[["Pa"]] - m[["Pa"]]) / sqrt(v[["Pa"]])
  )
}
