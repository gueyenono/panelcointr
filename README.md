# panelcointr

Panel cointegration tests for R that stay valid when the units of a panel
(countries, regions, firms) are hit by common shocks or go through
structural breaks. The computations, including the bootstraps and
simulations, are written in Rust and run in parallel.

| Function | Test | Allows for |
|---|---|---|
| `pcoint_westerlund()` | Westerlund (2007) error-correction tests Gt, Ga, Pt and Pa | cross-section dependence, through a bootstrap; also reproduces the Stata command `xtwest` |
| `pcoint_banerjee()` | Banerjee and Carrion-i-Silvestre (2015) | common factors, and structural breaks that are known or estimated (one or two per unit) |

## Installation

panelcointr is not on CRAN yet. Install it from GitHub:

```r
# install.packages("remotes")
remotes::install_github("gueyenono/panelcointr")
```

Add `build_vignettes = TRUE` to install the beginner's guides as well; this
needs the knitr and rmarkdown packages.

The package compiles Rust code when it is installed, so you need the Rust
toolchain: `cargo` and `rustc` 1.71 or later.

- **Linux and macOS:** install Rust from <https://rustup.rs> (or your
  system's package manager) and check that `cargo --version` works in a new
  terminal.
- **Windows:** install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  and Rust from <https://rustup.rs>, then add the target R uses:
  `rustup target add x86_64-pc-windows-gnu`.

The Rust libraries the package uses are bundled with it, so the build does
not download anything, and the package needs no R packages beyond base R.
So far it has been tested on Linux only.

## Example

The package ships `hce`: health care expenditure and GDP for 20 OECD
countries, 1970-2001.

```r
library(panelcointr)

# Westerlund (2007). The bootstrap p-values allow for countries that are
# hit by the same shocks; report those.
set.seed(2007)
pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
                  deterministic = "trend", lags = c(1, 3), leads = c(0, 3),
                  bootstrap = 499)

# Banerjee and Carrion-i-Silvestre (2015), Model 2 (constant and trend),
# with common factors estimated from the data.
pcoint_banerjee(loghex ~ loggdp, data = hce, index = c("ctr", "year"), model = 2)
```

## Learning more

- The beginner's guides, `vignette("westerlund", package = "panelcointr")`
  and `vignette("banerjee", package = "panelcointr")`, explain what each
  test asks, how to run it, how to read every line of its output, and what
  to check before trusting the result.
- `?pcoint_westerlund` and `?pcoint_banerjee` document every argument.
- [`inst/validation/RESULTS.md`](inst/validation/RESULTS.md) records how the
  results were checked against the papers and against other software.

## References

- Banerjee, A. and Carrion-i-Silvestre, J. L. (2015). Cointegration in panel
  data with structural breaks and cross-section dependence. *Journal of
  Applied Econometrics*, 30(1), 1-23. <https://doi.org/10.1002/jae.2348>
- Persyn, D. and Westerlund, J. (2008). Error-correction-based
  cointegration tests for panel data. *Stata Journal*, 8(2), 232-241.
- Westerlund, J. (2007). Testing for error correction in panel data.
  *Oxford Bulletin of Economics and Statistics*, 69(6), 709-748.
  <https://doi.org/10.1111/j.1468-0084.2007.00477.x>

## License

MIT. The bundled Rust libraries keep their own licenses, listed in
[`LICENSE.note`](LICENSE.note).
