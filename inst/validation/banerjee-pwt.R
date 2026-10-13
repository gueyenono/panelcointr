# T9 of plan/banerjee.md: golden test against Banerjee and Carrion-i-Silvestre
# (2017, Birmingham DP 2017-02, Table 16, p. 45 = PDF p. 50). The table applies
# the 2015 test to a production function,
#   y_it = alpha_i + beta_1 l_it + beta_2 k_it + u_it,
# on PWT 6.3 data for 19 OECD countries over 1951-2007 (N = 19, T = 57). Model
# 1 has no breaks. The table gives Z for r = 0..6 fixed factors, with and
# without dividing by the standard deviation before the principal components
# (the "transformed" rows, our pca_scale), and the number of I(1) factors r1 by
# MQ_c and MQ_f. The text (p. 19) gives the choice of the "panel BIC" with at
# most six factors: r = 6 untransformed, and "only one integrated common
# factor" transformed, read as r = 1. The target is the package default,
# IC_p1 (plan B4); BIC3, the default until this test, is shown for reference.
#
# The authors' data were never published. data-raw/pwt-production.R rebuilds
# the panel, with a capital stock from an assumed perpetual inventory recipe.
# The DP does not give the ADF lag rule or the moments either. The script
# therefore searches a grid:
#   capital   delta in {0.05, 0.03, 0.06, 0.07}, with the initial K/Y from
#             whole-sample or first-decade averages;
#   lags      t-sig (the package default), BIC with at most 4 lags (the
#             authors' GAUSS code, plan B6), or fixed at 0..4;
#   moments   asymptotic (the default) or finite_sample (WP Table 9; the
#             GAUSS code uses the T = 100 values, plan B9).
# Targets (plan section 8): Z to the printed two decimals, r1 and the
# criterion's choice exactly. A mismatch is documented, not forced. Z at r = 0 involves no
# factors, so it checks the data, lags and moments on their own.
#
# Needs readxl, an internet connection (for the data) and a build of
# panelcointr that has pcoint_banerjee(); a release build is faster. Run from
# the package root: Rscript inst/validation/banerjee-pwt.R

library(panelcointr)
source("data-raw/pwt-production.R")

r_all <- 0:6
target <- list(
  Z = rbind(raw = c(-3.56, -3.16, -3.09, -2.84, -2.75, -2.46, -2.53),
            scaled = c(-3.56, -3.79, -3.75, -3.76, -3.88, -4.08, -3.62)),
  r1_c = rbind(raw = c(1, 1, 3, 4, 5, 6), scaled = c(1, 2, 3, 4, 5, 6)),
  r1_f = rbind(raw = c(1, 2, 3, 4, 5, 6), scaled = c(1, 2, 3, 4, 5, 6)),
  choice = c(raw = 6, scaled = 1)
)
colnames(target$Z) <- r_all
colnames(target$r1_c) <- colnames(target$r1_f) <- r_all[-1]

capital <- expand.grid(delta = c(0.05, 0.03, 0.06, 0.07), init_years = c(NA, 10))
capital$label <- sprintf("delta %.2f, K0 from %s", capital$delta,
                         ifelse(is.na(capital$init_years), "1951-2007", "1951-1960"))
lag_settings <- c(list("t-sig" = list(lags = NULL, crit = "tsig"),
                       "BIC" = list(lags = NULL, crit = "bic")),
                  setNames(lapply(0:4, function(k) list(lags = k, crit = "tsig")),
                           paste("k =", 0:4)))
moments <- c("asymptotic", "finite_sample")

# Table 16 for one panel and one setting of the lags and moments.
table16 <- function(d, lag, mom) {
  fit <- function(factors, scale) {
    ctl <- banerjee_control(lag_criterion = lag$crit, pca_scale = scale, moments = mom)
    pcoint_banerjee(y ~ l + k, data = d, index = c("country", "year"), model = 1,
                    factors = factors, max_factors = 6L, lags = lag$lags, control = ctl)
  }
  out <- list(Z = target$Z * NA, r1_c = target$r1_c * NA, r1_f = target$r1_f * NA,
              choice = target$choice * NA, bic3 = target$choice * NA)
  for (s in c("raw", "scaled")) {
    for (r in r_all) {
      if (r == 0 && s == "scaled") {
        out$Z[s, "0"] <- out$Z["raw", "0"]
        next
      }
      res <- fit(r, s == "scaled")
      out$Z[s, as.character(r)] <- res$statistics$z
      if (r > 0) {
        out$r1_c[s, as.character(r)] <- res$factors$r1_c
        out$r1_f[s, as.character(r)] <- res$factors$r1_f
      }
    }
    out$choice[s] <- fit("ic1", s == "scaled")$factors$r
    out$bic3[s] <- fit("bic3", s == "scaled")$factors$r
  }
  out
}

# The scaled Z at r = 0 repeats the unscaled one, so 13 distinct Z.
distinct <- !(row(target$Z) == 2L & col(target$Z) == 1L)
score <- function(res) {
  dz <- abs(res$Z - target$Z)[distinct]
  c(z_match = sum(abs(round(res$Z, 2) - target$Z)[distinct] < 1e-8),
    z_mad = mean(dz), z_max = max(dz),
    z0 = res$Z["raw", "0"],
    r1_match = sum(res$r1_c == target$r1_c) + sum(res$r1_f == target$r1_f),
    choice_match = sum(res$choice == target$choice))
}

grid <- expand.grid(capital = seq_len(nrow(capital)), lags = names(lag_settings),
                    moments = moments, stringsAsFactors = FALSE)
panels <- lapply(seq_len(nrow(capital)), function(i) {
  build_pwt_production(pwt63, delta = capital$delta[i],
                       init_years = if (is.na(capital$init_years[i])) NULL else capital$init_years[i])
})
t0 <- proc.time()[["elapsed"]]
results <- lapply(seq_len(nrow(grid)), function(j) {
  table16(panels[[grid$capital[j]]], lag_settings[[grid$lags[j]]], grid$moments[j])
})
message(nrow(grid), " settings in ", round(proc.time()[["elapsed"]] - t0, 1), " s")
scores <- cbind(data.frame(capital = capital$label[grid$capital], lags = grid$lags,
                           moments = grid$moments, stringsAsFactors = FALSE),
                do.call(rbind, lapply(results, score)))
ranked <- order(-scores$z_match, scores$z_mad)

show <- function(j, title) {
  res <- results[[j]]
  cat("\n", title, ": ", scores$capital[j], "; lags ", scores$lags[j], "; moments ",
      scores$moments[j], "\n", sep = "")
  for (s in c("raw", "scaled")) {
    cat(if (s == "raw") "Not transformed" else "Transformed", "\n")
    z <- function(v) sprintf("%.2f", v)
    r1 <- function(v) c("-", as.character(v))
    tab <- rbind("Z, Table 16" = z(target$Z[s, ]), "Z, ours" = z(res$Z[s, ]),
                 "r1 MQ_c, Table 16" = r1(target$r1_c[s, ]), "r1 MQ_c, ours" = r1(res$r1_c[s, ]),
                 "r1 MQ_f, Table 16" = r1(target$r1_f[s, ]), "r1 MQ_f, ours" = r1(res$r1_f[s, ]))
    colnames(tab) <- paste0("r=", r_all)
    print(tab, quote = FALSE, right = TRUE)
    cat("IC_p1 choice: ", res$choice[s], " (text: ", target$choice[s], "); BIC3: ", res$bic3[s], "\n", sep = "")
  }
  print(scores[j, -(1:3)], row.names = FALSE, digits = 3)
}

cat("Data:", length(unique(pwt_production$country)), "countries,",
    paste(range(pwt_production$year), collapse = "-"), "\n")
cat("\nZ at r = 0 (Table 16: -3.56) with the default capital series:\n")
z0 <- with(scores[grid$capital == 1L, ], tapply(z0, list(lags, moments), function(v) round(v, 2)))
print(z0[names(lag_settings), ])

cat("\nTen best settings (by Z entries matched to 2 decimals, of 13; then mean |Z diff|):\n")
print(head(scores[ranked, ], 10), row.names = FALSE, digits = 3)

show(which(grid$capital == 1L & grid$lags == "t-sig" & grid$moments == "asymptotic"),
     "Package defaults with the default capital series")
show(ranked[1], "Best setting")

best <- scores[ranked[1], ]
pass <- best$z_match == 13 && best$r1_match == 24 && best$choice_match == 2
cat("\nT9:", if (pass) "PASS" else "NOT REPRODUCED", "- best setting matches",
    best$z_match, "of 13 Z,", best$r1_match, "of 24 r1 and", best$choice_match, "of 2 IC_p1 choices\n")
