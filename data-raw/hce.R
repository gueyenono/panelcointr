# Health care expenditure and GDP for 20 OECD countries, 1970-2001: the data
# of Westerlund (2007, Section 7), distributed with Stata's xtwest (Persyn and
# Westerlund 2008, Stata Journal 8(2), st0146) as xtwestdata.dta.
#
# Run from the package root: Rscript data-raw/hce.R

url <- "http://fmwww.bc.edu/repec/bocode/x/xtwestdata.dta"
path <- tempfile(fileext = ".dta")
utils::download.file(url, path, mode = "wb", quiet = TRUE)
raw <- haven::read_dta(path)

hce <- data.frame(
  ctr = as.integer(raw$ctr),
  year = as.integer(raw$year),
  hex = as.numeric(raw$hex),
  gdp = as.numeric(raw$gdp),
  loghex = as.numeric(raw$loghex),
  loggdp = as.numeric(raw$loggdp)
)
hce <- hce[order(hce$ctr, hce$year), ]
rownames(hce) <- NULL
stopifnot(nrow(hce) == 640L, !anyNA(hce),
          all.equal(hce$loghex, log(hce$hex), tolerance = 1e-6),
          all.equal(hce$loggdp, log(hce$gdp), tolerance = 1e-6))

usethis::use_data(hce, overwrite = TRUE)
