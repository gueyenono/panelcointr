# Plan: `pcoint_westerlund()`: Westerlund (2007) error-correction panel cointegration tests

Status: plan, not yet implemented. Written 2026-10-11.

## 0. Summary

`pcoint_westerlund()` implements the four error-correction-based panel cointegration tests of
Westerlund (2007): the group-mean tests **Gt** and **Ga**, and the panel tests **Pt** and **Pa**. It
reports asymptotic p-values from the moments in the paper's Table 1, plus optional bootstrap p-values
that are robust to cross-sectional dependence (paper §5).

The design principles are:

1. **The paper is the specification; xtwest is a validation target.** Where Stata's `xtwest`
   departs from the paper, the default follows the paper. A `compat = "xtwest"` switch reproduces
   `xtwest` to the published digits. Every divergence is listed in §3.
2. **R handles the interface; Rust handles the numerics.** R does input handling, validation, RNG
   draws, the moment tables, p-values and S3 methods. Rust does all per-unit estimation, lag/lead
   selection, long-run variances, the statistics, and the whole bootstrap loop (the hot path).
3. **Shared infrastructure for the rest of the package.** Panel preparation, the `pcoint_test`
   result class and the Rust kernels (least squares, long-run variance, panel layout, resampling) are
   written once and reused by later tests (Westerlund 2008; Banerjee & Carrion-i-Silvestre 2015, 2017).

**Already verified.** I wrote a pure-R prototype of my reading of `xtwest.ado`
(`plan/westerlund-prototype.R`). It reproduces all three published `xtwest` runs in Persyn &
Westerlund (2008, Stata Journal) on `xtwestdata`, statistics and z-values, to the printed precision:

| xtwest call (`loghex loggdp`, `constant trend`) | Published Gt / Ga / Pt / Pa | Prototype |
|---|---|---|
| `lags(1) leads(1) lrwindow(3)` | −2.681 / −10.927 / −12.035 / −10.524 | identical |
| `lags(1 2) leads(0 1) lrwindow(2)` (AIC) | −2.736 / −11.254 / −12.859 / −11.773 | identical (Ga −11.2535, see Q14) |
| `lags(1 3) leads(0 3) lrwindow(3) westerlund` (AIC) | −4.082 / −27.702 / −12.969 / −22.470 | identical; z −9.613 / −10.626 / −4.100 / −10.119 |

The average AIC-selected lag and lead lengths also match (1.15 / 0.05 and 2.8 / 1.65). So the
algorithm description below is grounded in code that has been checked against published numbers,
not only in a reading of the paper.

---

## 1. Sources

| Source | Use |
|---|---|
| Westerlund (2006), METEOR RM/06/056 (`literature/guid-984886be-…pdf`), the working-paper version of Westerlund (2007, OBES 69(6):709–748) | Specification: model (§2), test construction (§3), asymptotics and **Table 1 moments** (§4), bootstrap (§5), Monte Carlo design (§6), empirical results (Table 7) |
| `xtwest.ado` v1.5 (1 Jul 2010) and `xtwest.hlp`, SSC (`http://fmwww.bc.edu/repec/bocode/x/`) | Reference implementation: options, defaults, numerical details, quirks (§3, Appendix A) |
| Persyn & Westerlund (2008), Stata Journal 8(2):232–241 (st0146) | Published `xtwest` output used for golden tests |
| `xtwestdata.dta` (20 OECD countries × 32 years: `ctr`, `year`, `loghex`, `loggdp`, …) | Golden-test data (see open question 1) |

Table 1 of the paper and the moment matrices hard-coded in `xtwest.ado` (lines 1217–1225) agree
entry for entry. CRAN already has an R package called `Westerlund` (v0.1.4, Oct 2026), which
describes itself as a "functional approximation". We will compare against it for information only.
It is not a dependency and not an oracle.

---

## 2. The method: precise specification

### 2.1 Model and hypotheses

For unit `i = 1..N` with contiguous observations `t = 1..T_i` and `K` regressors `x_it` (all I(1)),
the conditional error-correction model (paper eq. 7, with leads as in Remark 5) is:

```
Δy_it = δ_i' d_t + α_i y_i,t-1 + λ_i' x_i,t-1 + Σ_{j=1..p_i} α_ij Δy_i,t-j
        + Σ_{j=-q_i..p_i} γ_ij' Δx_i,t-j + e_it
```

The deterministic cases are:

| Case | `deterministic =` | `d_t` | xtwest |
|---|---|---|---|
| 1 | `"none"` | (empty) | default |
| 2 | `"constant"` | `1` | `constant` |
| 3 | `"trend"` | `(1, t)'` | `constant trend` |

The null is H0: `α_i = 0` for all i (no error correction, hence no cointegration). The alternatives:

- Group-mean tests Gt and Ga: H1g, `α_i < 0` for at least one i.
- Panel tests Pt and Pa: H1p, `α_i = α < 0` for all i.

All four statistics diverge to −∞ under the alternative, so the tests are left-tailed.

### 2.2 Per-unit regression (for a given `p`, `q`)

- Estimation sample: `S_i(p,q) = { t : p+2 ≤ t ≤ T_i − q }`, with size `n_i = T_i − p − q − 1`.
- Regressor vector, ordered with **`y_{t-1}` last** (this matters; see §6.3):
  `z_t = [ d_t, x_{t-1}', Δy_{t-1..t-p}, Δx'_{t+q..t-p}, y_{t-1} ]`.
  The number of columns is `k_i = dim(d) + K + p + K(p+q+1) + 1`.
- Estimate by OLS (Householder QR). This gives `α̂_i`, the residuals `ê_it`, and
  `s_i² = RSS_i / (n_i − k_i)`.
- `SE(α̂_i) = s_i · sqrt([(Z'Z)^{-1}]_{kk})`, the conventional standard error.
- Composite error (paper §3.1):
  `û_it = Δy_it − δ̂'d_t − α̂_i y_{t-1} − λ̂_i'x_{t-1} − Σ α̂_ij Δy_{t-j}`, which equals
  `Σ γ̂_ij'Δx_{t-j} + ê_it`, for `t ∈ S_i`.

### 2.3 Long-run variances and `α̂_i(1)`

The Bartlett-kernel estimator, with bandwidth `M_i`, is:

```
LRV(v; M) = (1/n) [ Σ_t v_t² + 2 Σ_{j=1..M} (1 − j/(M+1)) Σ_t v_t v_{t-j} ]
```

The sums run over the available pairs and `n` is the number of available observations. No
demeaning is done inside the estimator; this is `xtwest`'s `lrvar …, weighted nodemean`.

- `ω̂²_yi = LRV(Δy*_i)` over all available `t = 2..T_i`. `Δy*` is `Δy` with its implied
  deterministic part removed (decision D8): unchanged in Cases 1 and 2, demeaned in Case 3.
- `ω̂²_ui = LRV(û_i)` over `S_i`.
- `α̂_i(1) = ω̂_ui / ω̂_yi`, the semiparametric estimator of the paper.

### 2.4 Statistics

With `T̃_i` the normalising sample size (D2; default `T̃_i = n_i`, the effective number of observations):

```
Gt = (1/N) Σ_i α̂_i / SE(α̂_i)
Ga = (1/N) Σ_i T̃_i α̂_i / α̂_i(1)
```

The panel statistics use the per-unit projection errors `Δỹ_it` and `ỹ_i,t-1`, which are the residuals
of `Δy_t` and `y_{t-1}` on all other regressors (paper §3.2; Frisch–Waugh–Lovell on the same regression):

```
α̂     = [ Σ_i Σ_t ỹ²_i,t-1 ]^{-1} Σ_i Σ_t (1/α̂_i(1)) ỹ_i,t-1 Δỹ_it
S_i   = σ̂_i / α̂_i(1),   σ̂_i = s_i  (regression standard error of unit i)
S_N²  = (1/N) Σ_i S_i²
SE(α̂) = S_N / sqrt( Σ_i Σ_t ỹ²_i,t-1 )
Pt    = α̂ / SE(α̂)
Pa    = T̄ α̂,   T̄ = (1/N) Σ_i T̃_i
```

### 2.5 Standardisation and asymptotic p-values

With `(μ, v)` from Table 1 for the case and K (rows Case 1..3, columns K = 1..6):

```
Z_Gt = √N (Gt − μ_Gt) / √v_Gt        μ = Θ̃₂,         v = Σ̃₂₂
Z_Ga = √N (Ga − μ_Ga) / √v_Ga        μ = Θ̃₁,         v = Σ̃₁₁
Z_Pt = (Pt − √N μ_Pt) / √v_Pt        μ = Θ₂ Θ₁^{-1/2}, v = ϕ'Σϕ   (Pt is not multiplied by √N)
Z_Pa = √N (Pa − μ_Pa) / √v_Pa        μ = Θ₂ Θ₁^{-1},   v = φ'Σφ
p    = Φ(Z)                          (left tail)
```

The moments are stored as a literal in `R/westerlund-moments.R`, with a comment citing Table 1. For
K > 6 there are no moments: z and the asymptotic p-value are `NA`, with a warning, and the bootstrap
still works. A validation script regenerates Table 1 by simulation to confirm the literal (§8, T7).

### 2.6 Lag and lead selection

`lags` and `leads` are each either a scalar (fixed) or a `c(min, max)` range. If either is a range,
each unit's `(p_i, q_i)` is chosen over the grid by an information criterion. Default: AIC,
`n log(RSS/n) + 2k`, or BIC `n log(RSS/n) + k log n`, computed on the **common sample**
`S_i(p_max, q_max)` (Ng & Perron 2005). Ties go to the smaller model. The chosen model is then
refitted on its own full sample `S_i(p_i, q_i)`.

### 2.7 Bootstrap (paper §5, extended to leads)

1. **Null-imposed first stage.** For each i, with its selected `(p_i, q_i)`, regress
   `Δy_it` on `Δy_{t-1..t-p_i}` and `Δx_{t+q_i..t-p_i}`, adding a constant in Case 3 only (paper §5).
   This gives `α̂_ij`, `γ̂_ij` and residuals `ê_it`.
2. **Joint residual vector.** On the periods `𝒯` where every unit has a residual, form
   `w_t = (ê_1t..ê_Nt, Δx_1t'..Δx_Nt')`. Centre it by its mean over `𝒯`.
3. **Resample common time indices.** Draw them with replacement (the same index for all units),
   which preserves cross-sectional dependence and cross-unit endogeneity. **The indices are drawn in
   R** (`sample.int`), so `set.seed()` gives reproducible results whatever the thread count.
4. **Generate the bootstrap panel.**
   `u*_it = e*_it + Σ_{j=−q_i..p_i} γ̂_ij'Δx*_i,t−j`, then
   `Δy*_it = Σ_j α̂_ij Δy*_i,t−j + u*_it`, started at zero, with the first `burn` values discarded.
   `y* = cumsum(Δy*)` and `x* = cumsum(Δx*)`. Unit i keeps `T_i` observations.
5. **Recompute.** Compute Gt, Ga, Pt and Pa on `(y*, x*)` with the same settings. By default the
   lags and leads are re-selected in each replicate (control `boot_reselect`).
6. **p-values.** `p*_s = (1 + #{S*_b ≤ S}) / (B + 1)`. Also report the bootstrap 1%, 5% and 10%
   critical values.

---

## 3. Where the paper and xtwest differ: design decisions

Defaults follow the paper. `compat = "xtwest"` gives `xtwest` defaults, and
`compat = "xtwest_westerlund"` gives `xtwest`'s `westerlund` option (Westerlund's original GAUSS
conventions). Individual choices can be overridden through `westerlund_control()`.

| # | Item | Paper | xtwest (default) | **Our default** | Notes |
|---|---|---|---|---|---|
| D1 | Lags/leads for the panel statistics | individual `p_i`, `q_i` (§3.2 step 1) | `floor(mean p_i)`, `floor(mean q_i)` for all units (Q1) | **individual** | `xtwest_westerlund` also floors the mean, so Westerlund's own code did; we follow the text. Control: `panel_lags = c("individual","mean")` |
| D2 | `T̃_i` in Ga | T; recommends the effective obs per unit | residual df `n_i − k_i` (Q2) | **`n_i`** (effective obs) | `xtwest_westerlund` uses `n_i` too. Control: `t_norm = c("effective","df","T")` |
| D3 | `T̄` in Pa | T; recommends the mean effective obs | `T̄`-based df (Q2) | **mean of `n_i`** | tied to D2 |
| D4 | `SE(α̂_i)` | "conventional" | df-corrected (`westerlund`: RSS/n) | **df-corrected** | |
| D5 | `σ̂_i` in `S_i` | regression s.e. of unit i | df from the mean T for all units (Q4) | **own df** | identical in balanced panels with fixed lags |
| D6 | Averaging in Ga | simple mean over i | row-weighted, i.e. weighted by `T_i` (Q3) | **simple mean** | identical in balanced panels |
| D7 | Information criterion | AIC (sample unspecified) | Stata AIC, each candidate on its own sample; ties go to the larger model (Q5) | **AIC/BIC on the common sample; ties go to the smaller model** | comparing ICs across different samples is a known flaw (Ng & Perron 2005) |
| D8 | Deterministics removed from Δy for `ω̂_y` | "residuals of Δy on `d_t`" | demean in Case 3 only (Q7) | **implied by the model: demean in Case 3 only** | Under the model Δy has a drift only in Case 3, so detrending Δy removes a term it does not have. Control: `lrv_deterministics = c("implied","regression")` |
| D9 | Sample for `ω̂_u` | `S_i` | `S_i` plus the last `q` rows (Q8) | **`S_i`** | |
| D10 | Bootstrap first-stage deterministics | constant in Case 3 only | user's `d_t` (Q10) | **paper** | |
| D11 | Bootstrap lag orders | `p_i` | re-selected on the null regression (Q10) | **`(p_i, q_i)` from the test regression** | |
| D12 | Bootstrap statistic | "constructed exactly as its sample counterpart" | re-selects lags in each replicate | **re-select** | `boot_reselect = FALSE` is faster |
| D13 | Bootstrap p-value | quantile rule | `#{S* ≤ S}/B` | **`(1 + #)/(B + 1)`** | never exactly 0 |
| D14 | Bootstrap initialisation | generate `T + n`, discard `n`, start at zero | essentially no burn-in | **`burn = 50`** (control) | |
| D15 | Bandwidth `M_i` | MC uses `4(T/100)^{2/9}` | `lrwindow`, default 2 | **`lrwindow = NULL` → `floor(4 (T_i/100)^{2/9})` per unit** | user value used as given |
| D16 | Deterministic default | n/a | none (Case 1) | **`"constant"`** | paper Remark 4: Case 1 is inconsistent against cointegration with a constant |
| D17 | Moments | Table 1 | Table 1; separate legacy moments with `westerlund` (Q13) | **Table 1** | legacy moments only under `compat = "xtwest_westerlund"` |
| D18 | Trend regressor | `t` | raw time variable (tests), `_n` (bootstrap) (Q9) | **`t = 1..T_i`** | numerically irrelevant: a trend always comes with a constant |

`compat` applies to the statistics. It cannot reproduce `xtwest`'s bootstrap draws (different RNG),
and it does not reproduce `xtwest`'s bootstrap bugs (Q10). With `compat` set, the bootstrap uses our
algorithm on the `compat` statistics.

---

## 4. User-facing API

```r
pcoint_westerlund(
  formula,                     # y ~ x1 + x2, levels of I(1) series
  data,                        # data.frame / tibble / data.table / plm::pdata.frame
  index = NULL,                # c("id", "time"); NULL → pdata.frame index, else first two columns
  deterministic = c("constant", "trend", "none"),  # "trend" means constant + trend
  lags,                        # required: scalar (fixed) or c(min, max) (IC selection)
  leads = 0,                   # scalar or c(min, max)
  lrwindow = NULL,             # Bartlett bandwidth; NULL → floor(4 (T_i/100)^(2/9)) per unit
  criterion = c("aic", "bic"),
  bootstrap = 0L,              # number of replications; 0 = asymptotic only
  compat = c("none", "xtwest", "xtwest_westerlund"),
  control = westerlund_control(),
  n_threads = getOption("panelcointr.threads", 1L)
)

westerlund_control(
  t_norm = NULL,               # "effective" | "df" | "T"
  panel_lags = NULL,           # "individual" | "mean"
  ic_sample = NULL,            # "common" | "varying"
  lrv_deterministics = NULL,   # "implied" | "regression"
  boot_reselect = NULL,        # TRUE | FALSE
  boot_burn = NULL,            # integer
  boot_pvalue = NULL,          # "plus_one" | "xtwest"
  boot_chunk = 100L,           # replications per Rust call (interrupt and progress granularity)
  min_df = 1L,                 # minimum n_i − k_i per unit
  keep_boot = TRUE             # keep the B × 4 matrix of bootstrap statistics
)
```

- `lags` has no default, like `xtwest`, which errors without `lags()`. The SJ article stresses how
  sensitive results are to this choice, so it should be a deliberate one.
- `westerlund_control()` stores only the arguments that were actually supplied. Resolution order:
  explicit control value, then the `compat` preset, then the default. The resolved settings are
  stored in the result, so a reader can always see exactly what was computed.
- Input rules:
  - Rows are sorted by (id, time).
  - Leading and trailing `NA`s in a unit are dropped. Interior `NA`s or gaps in time are an error
    that lists the offending units.
  - Numeric or integer time must step by 1 within each unit. Dates and factors are mapped to their
    rank over the union of time values.
  - Units too short for the requested `(p_max, q_max)` (`n_i − k_i < min_df`) are an error that
    lists the units and their `T_i`.
  - Non-finite values are rejected.

**Return value**, class `c("pcoint_westerlund", "pcoint_test")`:

| Field | Content |
|---|---|
| `statistics` | data.frame: `statistic` (Gt, Ga, Pt, Pa), `value`, `z`, `p_value`, `p_boot` (NA without bootstrap), `alternative` |
| `units` | data.frame per unit: `id, T, n, k, lags, leads, lrwindow, alpha, se_alpha, t_alpha, alpha1, sigma` |
| `coefficients` | per-unit ECM coefficients (NA-padded matrix) for the mean-group summary |
| `boot` | `B × 4` statistics, critical values, `B`, seed information (if bootstrapped) |
| `settings` | resolved settings including control, compat and the moments used |
| `N`, `K`, `case`, `call`, `formula` | |

**Methods:**

- `print()` gives an `xtwest`-style table (Value, Z, P, Robust P) plus the average selected lags and
  leads.
- `summary()` adds the distribution of selected lags and leads, the unit table, and the mean-group
  ECM with implied long-run coefficients `−λ̂_i/α̂_i` (`xtwest`'s `mg`).
- `plot()` draws the bootstrap distributions with the observed statistics (base graphics, no
  dependency).
- `as.data.frame()` returns the statistics table.
- `tidy()` and `glance()` are registered for `generics` through
  `@exportS3Method generics::tidy`, with `generics` in Suggests.

```r
pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
                  deterministic = "trend", lags = c(1, 3), leads = c(0, 3), lrwindow = 3)

set.seed(1)
pcoint_westerlund(loghex ~ loggdp, data = hce, index = c("ctr", "year"),
                  deterministic = "trend", lags = 1, leads = 1, lrwindow = 3,
                  bootstrap = 999, n_threads = 2)
```

---

## 5. Architecture

```
R  pcoint_westerlund()
   ├─ prepare_panel()      [shared] validate, sort, trim, stack → y, X, offsets, ids, T_i
   ├─ resolve_settings()   compat preset + control + defaults → typed list
   ├─ westerlund_rs()      ── Rust: per-unit fits, selection, LRVs, Gt/Ga/Pt/Pa, unit table
   ├─ (bootstrap > 0)
   │    westerlund_boot_prepare_rs()  ── Rust: null-imposed first stage, centred w_t, valid periods
   │    loop over chunks:
   │       idx <- matrix(sample.int(...), chunk, L)   (R RNG)
   │       westerlund_boot_rs(state, idx, n_threads)  ── Rust: chunk × 4 statistics
   │       check for interrupts / progress
   ├─ moments lookup, z, p-values, bootstrap p-values  (R)
   └─ new_pcoint_westerlund()  (S3 constructor)
```

**What goes where.**

- **Rust** does everything that runs once per unit, once per candidate model, or once per bootstrap
  replicate. That is roughly `B × N × |grid| ×` (small dense least squares), which is where R would
  spend minutes.
- **R** keeps everything that is either trivial or about the user interface, and the RNG, so
  `set.seed()` works as R users expect.

**Data passed to Rust.**

- `y`: `f64` vector of length `Σ T_i`.
- `X`: `Σ T_i × K` column-major R matrix, read in Rust as K strided slices.
- `offsets`: integer vector of length `N + 1`.
- Per-unit integer vectors: `p_min`, `p_max`, `q_min`, `q_max`, `M_i`.
- A named settings list, parsed into a Rust struct.

Rust returns a named list of plain vectors and matrices. No R objects are created inside the numeric
core.

---

## 6. Rust engine

### 6.1 Binding framework: extendr (with rextendr)

| | extendr 0.9 / rextendr 0.5 | savvy 0.11 |
|---|---|---|
| R-side tooling | `rextendr::use_extendr()`, `use_cran_defaults()`, `vendor_pkgs()`, `document()`, `use_msrv()` | `savvy-cli` (Rust-side) |
| MSRV | 1.71 | 1.81 |
| Errors | `Result<T>` → R error; panics caught | `Result` only (no panics), explicit |
| Footprint | larger (libR-sys) | smaller |

**Recommendation: extendr.** It has the R-native workflow (devtools-style, CRAN templates and
vendoring from R), the lower MSRV (CRAN asks for testing against a cargo at least two years old),
and the larger body of CRAN examples. Our FFI surface is small (vectors in, a list out), so a later
switch to savvy would be cheap if it were ever needed.

**No linear-algebra crate.** The regressions are tiny (n ≈ 20–500, k ≈ 5–40), and a
hand-written Householder QR and Cholesky (about 200 lines, tested against `lm.fit`) beat general
libraries at this size. They also keep the vendored tarball small and the MSRV low. faer 0.24, for
example, brings a large dependency tree and a recent MSRV (its repository listed ≥ 1.84 when I
checked; confirm before relying on it). Threads use `std::thread::scope`, not rayon.
**The only dependency is `extendr-api`.**

### 6.2 Crate layout (`src/rust/`)

```
Cargo.toml              # crate-type = ["staticlib"]; rust-version pinned; deps: extendr-api only
src/lib.rs              # extendr_module!; thin wrappers only (conversion + validation + Result)
src/core/mod.rs         # pure Rust, no extendr types → unit-testable with `cargo test`
src/core/panel.rs       # Panel { y, x (strided), offsets, k } + unit views
src/core/linalg.rs      # Householder QR LS (small dense, column-major); Cholesky on Gram sub-blocks; rank checks
src/core/lrv.rs         # Bartlett LRV (O(nM))
src/core/design.rs      # ECM design builder into a reusable Workspace (no allocation in hot loops)
src/core/ecm.rs         # fit(p,q), IC grid, composite error, UnitResult
src/core/westerlund.rs  # Gt/Ga/Pt/Pa aggregation; compat variants
src/core/bootstrap.rs   # null-imposed first stage, DGP, replicate loop, thread chunking
```

`[profile.release]`: `opt-level = 3`, `lto = true`, `codegen-units = 1`. Keep `panic = "unwind"`
so extendr can turn panics into R errors.

### 6.3 Algorithms

1. **One QR per unit gives both the group and the panel ingredients.** With `y_{t-1}` as the last
   column, write `Z = QR`, and let `r = R_kk` and `c = (Q'Δy)_k`. Then:

   ```
   α̂_i = c / r
   SE(α̂_i) = s_i / |r|
   Σ ỹ²_{t-1} = r²
   Σ ỹ_{t-1} Δỹ_t = r · c        (Frisch–Waugh–Lovell)
   ```

   `xtwest` runs three regressions per unit for this; we run one. Only the compat mode with
   `panel_lags = "mean"` needs a second QR, at `(p̄, q̄)`.

2. **IC grid on the common sample through Gram sub-blocks.**
   - Build the design once for `(p_max, q_max)` on `S_i(p_max, q_max)`.
   - Standardise the columns, then form `G = Z'Z` and `g = Z'Δy` once.
   - Each candidate `(p, q)` is a column subset. Its Cholesky gives
     `RSS = Δy'Δy − ‖L⁻¹g_S‖²` in `O(k³/3)`, instead of an `O(n k²)` QR on fresh data.
   - Fall back to QR if a sub-block is not positive definite.
   - Refit the winner by QR on its own sample.

   With `compat = "xtwest"` (varying samples), each candidate gets its own QR.

3. **Bartlett LRV** in `O(nM)` over the available pairs.

4. **Bootstrap.**
   - Each thread owns its workspace buffers.
   - Replicates are split into contiguous blocks, and each block writes to its own rows of the
     output. Because the indices are pre-drawn, the results are bit-identical for any `n_threads`.
   - R calls Rust in chunks (`boot_chunk`) so it can check for interrupts and show progress.
     Worker threads never touch the R API.

5. **Numerical safety.**
   - Rank deficiency is detected with a tolerance relative to `max |R_jj|`. It is reported as an
     error naming the unit and the model, never as a silent `NA`.
   - Non-finite LRVs, or `ω̂²_y ≤ 0`, are reported per unit.

### 6.4 Threads

`n_threads` defaults to `getOption("panelcointr.threads", 1L)`. Examples, tests and vignettes use at
most 2 (CRAN policy). Only the bootstrap and multi-unit fitting are parallelised.

---

## 7. Packaging Rust for R (CRAN-compliant)

- **Scaffold.** Run `rextendr::use_extendr()` and then `rextendr::use_cran_defaults()`. Review the
  generated `configure`, `configure.win`, `src/Makevars.in`, `src/Makevars.win.in`,
  `src/entrypoint.c` and `R/extendr-wrappers.R` rather than trusting them blindly.
- **CRAN requirements** (from "Using Rust in CRAN packages", plus the general CRAN Repository
  Policy rule against writing to the user's home directory):
  - `SystemRequirements: Cargo (Rust's package manager), rustc >= <MSRV>`.
  - `configure` and `configure.win` check for `cargo` and `rustc` on `PATH` **and** in
    `~/.cargo/bin`, enforce the MSRV, never install anything themselves, and tell the user how to
    get Rust if it is missing.
  - The **rustc version is printed** in the install log before compilation starts.
  - `cargo build --jobs 2` (never the default of all cores), `--offline`, with `CARGO_HOME` set
    under the build directory so nothing is written to `~/.cargo`.
  - **Vendored crates**: `rextendr::vendor_pkgs()` produces `src/rust/vendor.tar.xz` (xz
    compression is required) and `.cargo/config.toml` replaces sources with the vendor directory.
  - **Authorship and copyright of vendored crates** go in DESCRIPTION, as a
    `person("Authors of the dependency Rust crates", role = "cph", comment = "see inst/AUTHORS")`
    entry, with `inst/AUTHORS` generated by `rextendr::write_license_note()`.
  - Test against a rustc at least two years old. In practice the MSRV is ≤ 1.75 (extendr needs
    1.71), checked in CI with `rustup toolchain install <MSRV>`.
- **Installed size**: build artefacts (`target/`) are removed after linking, and `src/rust/target`
  and `src/rust/vendor` are listed in `.Rbuildignore`.
- **Windows**: the `x86_64-pc-windows-gnu` target is needed (via rustup), and the Rtools linker
  setup comes from the rextendr templates. Binary CRAN builds mean end users need no Rust.
- **Local environment**: R 4.6.0, the R packages and rustc/cargo 1.95.0 (fenix stable) all come
  from nix via the user's home-manager flake (`~/home-manager/`). Missing tools are added there and
  deployed with `hms`, never with `install.packages()` or `rustup` (see CLAUDE.md, Environment).
  Everything else the plan uses is already installed (generics, haven, plm, data.table, tibble,
  testthat, devtools). Still missing, all present in the pinned nixpkgs:

  | Needed in | Package | nixpkgs version | Role |
  |---|---|---|---|
  | P0 | `rextendr` | 0.5.0 | scaffolding, vendoring, wrapper generation |
  | P7 | `bench` | 1.1.4 | §9 benchmarks |
  | P7 | `Westerlund` | 0.1.3 (CRAN has 0.1.4) | informal cross-check (§8) |
  | P7 | rustc/cargo at the MSRV | via the flake's existing `fenix` input (`toolchainOf`) | old-rustc check |

  The MSRV toolchain must not shadow the default `cargo`/`rustc` in `home.packages`. Expose it
  under a separate name (e.g. a `cargo-msrv-check` wrapper) or keep that check in CI only;
  decide at P7.
- **CLAUDE.md** gets the Rust commands (`rextendr::document()`, `cargo test` in `src/rust`,
  re-vendoring) once phase 0 is done.

---

## 8. Validation and testing

| ID | What | How | Tolerance |
|---|---|---|---|
| T1 | **Golden tests vs published xtwest** | the three SJ runs (§0) with `compat = "xtwest"` and `"xtwest_westerlund"`; statistics, z, average lags/leads | 1e-3 abs (printed precision; Stata stores temporaries as float, Q14) |
| T2 | **Rust vs pure-R oracle** | `tests/testthat/helper-reference.R`, ported from `plan/westerlund-prototype.R` and extended to the default variant; random panels over case × K ∈ {1,2,3} × fixed/IC × balanced/unbalanced × compat | 1e-10 relative |
| T3 | Kernels | QR least squares vs `lm.fit` (including near-collinear designs); LRV vs the direct formula; Gram-Cholesky IC selection vs brute-force QR; `cargo test` for linalg and lrv | 1e-12 |
| T4 | Invariances | trend coding (`t` vs `year`); scaling of y and x (Gt, Pt unchanged); permutation of units; shifting y by a constant in Cases 2 and 3 | 1e-10 |
| T5 | Bootstrap mechanics | same `set.seed` gives identical results; `n_threads` 1 vs 2 identical; `boot_reselect` on/off; unbalanced overlap handling | exact |
| T6 | Input validation | gaps, interior NA, short units (message lists them), K > 6 (NA asymptotic p with a warning), `pdata.frame` index, factor/Date time | exact messages |
| T7 | Monte Carlo (not run in `R CMD check`; `inst/validation/`) | (a) size for the paper's Table 2 designs (Case 1, N ∈ {10,20}, T ∈ {100,200}, φ, γ); (b) bootstrap size under cross-sectional dependence (Table 5, `b ∈ {0,1}`); (c) regenerate Table 1 moments (T = 1000, 10k draws) | within Monte Carlo error |

Also, for information only: compare with CRAN `Westerlund` on `xtwestdata`.

---

## 9. Performance targets

Single thread on this machine, measured with `bench::mark` in `inst/bench/` (not run on CRAN):

| Workload | Target |
|---|---|
| SJ example, AIC `lags(1 3) leads(0 3)`, no bootstrap | < 5 ms |
| SJ example, `bootstrap(800)`, fixed lags | < 0.5 s |
| SJ example, `bootstrap(800)`, AIC re-selection over 3 × 4 grid | < 2 s |
| N = 100, T = 200, K = 2, grid 5 × 5, B = 999, re-selection | < 60 s, near-linear speed-up with threads |

The pure-R oracle is benchmarked alongside to document the speed-up.

---

## 10. Implementation phases

Each phase ends with `devtools::check()` clean.

| Phase | Deliverable | Acceptance |
|---|---|---|
| P0 Tooling | install `rextendr` via home-manager (§7); scaffold extendr; CRAN defaults; a round-trip Rust function; vendoring; offline build | check passes with Rust; `cargo test` runs; install log shows the rustc version |
| P1 Shared R infrastructure | `prepare_panel()`, formula/index parsing, validation, `pcoint_test` constructor and print skeleton | T6 |
| P2 R oracle | `helper-reference.R` (default + compat variants) | T1 passes **on the oracle** |
| P3 Rust kernels | `linalg`, `lrv`, `design`/`Workspace` | T3 |
| P4 Rust statistics | fixed lags → IC grid → compat paths; `pcoint_westerlund()` without bootstrap | T1, T2 and T4 on the Rust path |
| P5 Bootstrap | Rust first stage and DGP, R-drawn indices, chunked loop, threads, bootstrap p-values and critical values | T5 |
| P6 Interface | print/summary/plot/tidy, mean-group ECM, roxygen docs, vignette replicating the SJ application, NEWS | docs build; vignette reproduces §0 |
| P7 Validation and release prep | T7 scripts, benchmarks (§9), old-MSRV check, win-builder/rhub, CLAUDE.md update | targets met |

---

## 11. Open questions

1. **Ship `xtwestdata` in the package?** It is 20 × 32, small, and lets the vignette and golden tests
   run offline. The Stata Journal copyright statement allows copying supporting files with
   attribution to the authors and the Journal, but the numbers come from OECD Health Data 2003.
   **Recommendation:** ship it as `data/hce` with full attribution. The alternative is for tests to
   download it and skip when offline.
2. **Defaults: paper or xtwest?** **Recommendation:** paper-faithful defaults with
   `compat = "xtwest"`, so the package is defensible against the source. The other option is
   `xtwest`-identical defaults, so Stata users see the same numbers.
3. **Is CRAN the target?** I've assumed yes, which drives the vendoring, MSRV and thread work.
4. **extendr or savvy?** **Recommendation:** extendr (§6.1).

---

## 12. References

- Westerlund, J. (2007). Testing for error correction in panel data. *Oxford Bulletin of Economics and Statistics* 69(6), 709–748. doi:10.1111/j.1468-0084.2007.00477.x. Working paper: METEOR RM/06/056 (2006).
- Persyn, D. & Westerlund, J. (2008). Error-correction-based cointegration tests for panel data. *Stata Journal* 8(2), 232–241.
- Persyn, D. `xtwest` v1.5 (2010), SSC s456941. http://fmwww.bc.edu/repec/bocode/x/xtwest.ado
- Chang, Y. (2004). Bootstrap unit root tests in panels with cross-sectional dependency. *Journal of Econometrics* 120, 263–293.
- Ng, S. & Perron, P. (2005). A note on the selection of time series models. *Oxford Bulletin of Economics and Statistics* 67(1), 115–134.
- Newey, W. K. & West, K. D. (1994). Automatic lag selection in covariance matrix estimation. *Review of Economic Studies* 61, 631–653.
- CRAN, "Using Rust in CRAN packages". https://cran.r-project.org/web/packages/using_rust.html
- extendr / rextendr: https://extendr.github.io/ · savvy: https://github.com/yutannihilation/savvy

---

## Appendix A. xtwest v1.5 quirks (line numbers refer to `xtwest.ado`)

| ID | Behaviour | Lines | Handled by |
|---|---|---|---|
| Q1 | The panel statistics use `int(mean)` of the selected lags and leads for every unit | 993–998, 1027–1077 | D1 |
| Q2 | Ga and Pa are normalised by the residual df; `westerlund` uses `n_i` and rescales SE to RSS/n | 1093–1099, 1122–1127 | D2–D4 |
| Q3 | `egen GA = mean(gai)` averages over rows, so Ga is weighted by `T_i` in unbalanced panels; Gt uses `tag()` and is unweighted | 1103–1112 | D6 |
| Q4 | `σ̂_i` uses the mean-T df for all units | 1126–1131 | D5 |
| Q5 | AIC from `estat ic`, each (p,q) on its own sample; loop runs max → min with strict `<`, so ties go to the larger model | 751–785 | D7 |
| Q6 | The `westerlund` IC formula uses the tempvar names `cons` and `tren` in a scalar expression, so they evaluate to their **first-observation** values (equal to 1 only if time starts at 1) | 764 | compat only |
| Q7 | `ω̂_y`: all available Δy, demeaned only in Case 3; `westerlund`: restricted to the estimation sample, never demeaned | 788–807, 1000–1025 | D8 |
| Q8 | `û` is also computed on the last `q` rows, outside the estimation sample (it omits the Δx terms, so leads are not needed) | 879–902 | D9 |
| Q9 | Trend is the raw time variable in the tests and `_n` in the bootstrap | 711, 388 | D18 |
| Q10 | Bootstrap: the first stage uses the user's `d_t`; lags are re-selected on the null regression; the lead-coefficient loops use the **last unit's** `currlead` (wrong when leads are AIC-selected); minimal burn-in; `expandcl` + `bsample, cluster(t)` + a shared random permutation, which amounts to common time indices drawn with replacement; p-value is `#{S* ≤ S}/B` | 347–610, 1273–1316 | D10–D14, not replicated |
| Q11 | Interior gaps are an error; minimum obs = `(p+q+1) + k + 1` (df ≥ 1) | 96–107, 164 | §4 input rules |
| Q12 | At most 6 regressors (the size of the moments table) | 43 | §2.5 (NA asymptotics beyond 6) |
| Q13 | `westerlund` uses legacy K = 1 moments, different from Table 1: with trend −2.356/0.6450, −11.8978/44.2471, −2.1128/0.7371, −8.9536/35.6802; without −1.793/0.7904, −7.2014/29.3677, −1.4746/1.0262, −4.3559/21.0535 | 1201–1213 | D17 |
| Q14 | `gen` creates float (single-precision) temporaries, so the last printed digit can differ (Ga −11.2535 prints as −11.254) | throughout | T1 tolerance |
