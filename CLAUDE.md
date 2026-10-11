# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`panelcointr` is an R package (R >= 4.6 available locally via nix) implementing panel cointegration tests robust to cross-section dependence and structural breaks. It is at the skeleton stage: no test functions are implemented yet.

## Commands

Run from the repo root:

```sh
Rscript -e 'devtools::document()'                 # regenerate NAMESPACE and man/ from roxygen
Rscript -e 'devtools::test()'                     # run all tests
Rscript -e 'devtools::test(filter = "package")'   # run one file: tests/testthat/test-package.R
Rscript -e 'devtools::check()'                    # full R CMD check (currently 0 errors/warnings/notes)
Rscript -e 'devtools::load_all()'                 # load package for interactive use
```

## Conventions

- Roxygen2 with markdown enabled; `NAMESPACE` and `man/` are generated, so edit roxygen comments in `R/` and rerun `devtools::document()` rather than editing them directly.
- testthat 3rd edition; test files are `tests/testthat/test-<name>.R`, mirroring `R/<name>.R`.
- Add dependencies with `usethis::use_package()` so DESCRIPTION stays consistent.
- Non-package directories (`literature/`, `plan/`, `CLAUDE.md`, `.claude/`) are listed in `.Rbuildignore`; add any new non-package top-level paths there too, or `R CMD check` will note them.

## Reference literature

`literature/` holds the papers that specify the methods. Treat them as the specification: test statistics, critical values and asymptotic moments should be checked against them.

- `Westerlund-PanelCointegrationTests-2008.pdf`: Westerlund (2008, *JAE*), "Panel Cointegration Tests of the Fisher Effect" (Durbin–Hausman tests with cross-section dependence).
- `guid-984886be-…-ASSET1.0.pdf`: Westerlund (2006), "Testing for error correction in panel data", METEOR RM 056 (ECM-based tests Gt, Ga, Pt, Pa).
- `BANERJEE-COINTEGRATIONPANELDATA-2015.pdf`: Banerjee & Carrion-i-Silvestre (2015, *JAE*), "Cointegration in panel data with structural breaks and cross-section dependence".
- `bbs-dp2017-02.pdf`: Banerjee & Carrion-i-Silvestre (2017), "Testing for panel cointegration using common correlated effects estimators", Birmingham Business School DP 2017-02.

Read them with `pdftotext`, e.g. `pdftotext -f 5 -l 10 literature/<file>.pdf -`. Several PDFs start with a JSTOR or repository cover page, so PDF page numbers are offset from the papers' own.

`plan/` is for planning documents.
