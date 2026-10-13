# Production-function panel of Banerjee and Carrion-i-Silvestre (2017,
# Birmingham DP 2017-02, Section 5.2) for golden test T9 of plan/banerjee.md
# (their Table 16): 19 OECD countries over 1951-2007, with y = log real GDP per
# capita, l = log population and k = log real capital stock per capita. The DP
# uses "the data in Banerjee, Eberhardt and Reade (2010)", which were never
# published, so the panel is rebuilt here from Penn World Table 6.3.
#
# Sample. The countries and years are the ones the DP lists. Where its sources
# disagree, the DP wins:
# - Finland is in the DP but not among the 55 countries of Banerjee, Eberhardt
#   and Reade (2010, Birmingham DP 10-33, Appendix A). Those 55 are the 60
#   countries with complete rgdpl, POP and ki over 1951-2007, less Cyprus,
#   Finland, Japan, Nicaragua and Trinidad and Tobago; nothing in PWT 6.3
#   explains the cut.
# - Japan is left out, as in the DP, although PWT 6.3 has it complete over
#   1950-2007 (the DP says it is not). Germany starts in 1970, as the DP says.
# - Greece has no rgdpl or ki in 1950, hence the 1951 start.
# - BER (2010) quote "year 2000 I$", the base year of PWT 6.2; PWT 6.3 is in
#   2005 I$. Rebasing scales each country's series by a constant, which only
#   shifts the logs and is absorbed by the unit intercepts of Model 1.
# - The file is the 18 November 2009 revision of PWT 6.3 (first released in
#   August 2009); the authors may have used an earlier one.
#
# Capital. BER (2010) say only that capital is built from ki by the perpetual
# inventory method. The recipe used here is the one Eberhardt and Teal (2011,
# CREDIT Research Paper 11/08, Appendix A-1.2) document for PWT data, after
# Klenow and Rodriguez-Clare (1997, NBER Macroeconomics Annual, p. 86):
#   Y_t = rgdpl_t POP_t,  I_t = (ki_t / 100) Y_t,
#   K_1951 = Y_1951 (I/Y) / (g + n + delta),
#   K_t = (1 - delta) K_{t-1} + I_t for t >= 1952,
# with delta = 0.05 and I/Y, g (growth of rgdpl) and n (growth of POP) the
# country's averages. The appendix does not say over which years; the default
# is the whole sample, and init_years = m uses the first m years instead (KRC
# use the first decade). These settings are assumptions, and T9 records which
# ones it used.
#
# Usage: source("data-raw/pwt-production.R") from the package root defines
# pwt63 (the raw series), build_pwt_production() and pwt_production, the panel
# with the default settings. It needs readxl and an internet connection. Per
# plan/banerjee.md (sections 7 and 11 Q7), nothing is written to data/ until T9
# passes and the terms of use are checked. The PWT 6.3 page asks only for the
# citation: Heston, A., Summers, R. and Aten, B. (2009), Penn World Table
# Version 6.3, Center for International Comparisons of Production, Income and
# Prices at the University of Pennsylvania.

iso <- c("AUS", "AUT", "BEL", "CAN", "DNK", "FIN", "FRA", "GRC", "IRL", "ITA",
         "LUX", "NLD", "NOR", "PRT", "ESP", "SWE", "CHE", "GBR", "USA")
years <- 1951:2007

read_pwt63 <- function() {
  zip <- file.path(tempdir(), "pwt63_nov182009version.zip")
  if (!file.exists(zip)) {
    utils::download.file("https://www.rug.nl/ggdc/docs/pwt63_nov182009version.zip", zip,
                         mode = "wb", quiet = TRUE)
  }
  stopifnot(unname(tools::md5sum(zip)) == "8143c0956d1f7627a542e443e0e65101")
  xls <- utils::unzip(zip, "pwt63_w_country_names.xls", exdir = tempdir())
  raw <- readxl::read_excel(xls)
  raw <- raw[raw$isocode %in% iso & raw$year %in% years, ]
  out <- data.frame(
    country = raw$country,
    year = as.integer(raw$year),
    rgdpl = as.numeric(raw$rgdpl),
    pop = as.numeric(raw$POP),
    ki = as.numeric(raw$ki)
  )
  out <- out[order(out$country, out$year), ]
  rownames(out) <- NULL
  out
}

# raw: the output of read_pwt63(). Returns one row per country and year with
# gdp (rgdpl, 2005 I$ per person), pop (thousands), capital (2005 I$ per
# person) and their logs y, l and k.
build_pwt_production <- function(raw, delta = 0.05, init_years = NULL) {
  stopifnot(delta > 0, delta < 1, is.null(init_years) || init_years >= 2)
  out <- lapply(split(raw, raw$country), function(d) {
    d <- d[order(d$year), ]
    w <- if (is.null(init_years)) seq_len(nrow(d)) else seq_len(init_years)
    gdp <- d$rgdpl * d$pop
    inv <- d$ki / 100 * gdp
    g <- mean(diff(log(d$rgdpl[w])))
    n <- mean(diff(log(d$pop[w])))
    stopifnot(g + n + delta > 0)
    cap <- numeric(nrow(d))
    cap[1] <- gdp[1] * mean(d$ki[w] / 100) / (g + n + delta)
    for (t in seq_len(nrow(d))[-1]) cap[t] <- (1 - delta) * cap[t - 1] + inv[t]
    data.frame(country = d$country, year = d$year, gdp = d$rgdpl, pop = d$pop,
               capital = cap / d$pop)
  })
  out <- do.call(rbind, out)
  out$y <- log(out$gdp)
  out$l <- log(out$pop)
  out$k <- log(out$capital)
  rownames(out) <- NULL
  out
}

pwt63 <- read_pwt63()
stopifnot(nrow(pwt63) == length(iso) * length(years), !anyNA(pwt63),
          all(table(pwt63$country) == length(years)))

pwt_production <- build_pwt_production(pwt63)
stopifnot(nrow(pwt_production) == 1083L, !anyNA(pwt_production),
          all(pwt_production$capital > 0))
