# Plan: `pcoint_banerjee()`: Banerjee & Carrion-i-Silvestre (2015) panel cointegration tests with structural breaks and common factors

Status: written 2026-10-12. The open questions were settled the same day: all recommendations
accepted (§11). Implemented 2026-10-12: P0-P6 (the vignette uses simulated data) and P7 except T9,
which waits for the PWT data; validation results in `inst/validation/RESULTS.md`. Choices made during
the implementation are B25-B31 (§3). Builds on the infrastructure of
`pcoint_westerlund()` (`plan/westerlund.md`).

## 0. Summary

`pcoint_banerjee()` implements the panel cointegration tests of Banerjee & Carrion-i-Silvestre
(2015, *JAE*; "BCS" below). The idea follows PANIC (Bai & Ng 2004):

- Cross-section dependence comes from unobserved common factors.
- The factors are estimated by principal components from the residuals of each unit's
  cointegrating regression **in first differences**.
- The null of no cointegration is tested on the defactored idiosyncratic errors, with a
  standardised pooled ADF statistic **Z**.
- The factors are tested separately: an ADF test if there is one factor, and the Bai–Ng MQ tests if
  there are more.
- Structural breaks in the level, the trend and the cointegrating vector are allowed. Their dates
  can be known or estimated: six models, three ways of handling breaks.

A conclusion needs both parts. The observed variables cointegrate only if Z rejects (the idiosyncratic
errors are I(0)) **and** no common factor is I(1). If Z rejects but some factors are I(1), the
variables cointegrate only together with those non-stationary factors (BCS §6.3.2).

The design principles are those of the Westerlund plan:

1. **The JAE paper is the specification.** Where it is silent, the earlier ECB working-paper version
   (BCS 2006, which has the proofs appendix the JAE version moved online) and BCS (2017) fill the
   gaps. Every remaining choice is listed in §3.
2. **R handles the interface; Rust handles the numerics.** This includes the new kernels: symmetric
   eigendecomposition, ADF with lag selection, principal components, break-date search and the
   simulation of null distributions.
3. **Shared infrastructure.**
   - Reused from Westerlund: `prepare_panel()`, the `pcoint_test` class, `linalg::Qr`,
     `linalg::cholesky` and `panel::Panel`.
   - New kernels written so the later tests can reuse them: the eigensolver, PCA and factor-number
     criteria (Westerlund 2008, BCS 2017) and the ADF routine (BCS 2017).

**There is no golden reference.**

- The paper's application uses Eurostat import unit values built as in Campa & González-Mínguez
  (2006); the data are not public.
- The only port of the authors' GAUSS code, CRAN `xtbreakcoint` 1.0.4, departs from the paper in
  most steps (Appendix A).

So, unlike Westerlund, nothing is yet known to match published output to the printed digits. The
one candidate is BCS (2017) Table 16 on PWT data (T9, §11 Q7), which depends on obtaining the data.
Validation rests on three things:

- reproducing the paper's simulated moment and critical-value tables;
- an independent pure-R oracle;
- invariance tests.

**Already verified.** A simulation of my reading of the procedure
(`plan/banerjee-moments-check.R`) reproduces the paper's null moments. The procedure as simulated:
project Δe on the model's first-difference deterministic columns, cumulate, then run an ADF regression
without deterministics. Settings: 20,000 draws of length T = 1000, k = 0, Monte Carlo s.e. ≈ 0.004–0.007
for the means.

| Configuration | Paper (Thm 2, Table II) mean / var | Simulated mean / var |
|---|---|---|
| c (Models 1, 4), no break | −0.424 / 0.964 | −0.422 / 0.961 |
| c, level break at 0.5 | same (moments do not depend on breaks) | −0.436 / 0.982 |
| τ (Models 2, 5), no break | −1.535 / 0.341 | −1.534 / 0.343 |
| τ, level break at 0.5 | same | −1.536 / 0.341 |
| γ (Models 3, 6), λ = 0.5 | −1.972 / 0.339 | −1.971 / 0.333 |
| γ, λ = 0.2 | −1.809 / 0.392 | −1.803 / 0.392 |
| γ, λ = (0.3, 0.7) | −2.315 / 0.337 | −2.305 / 0.343 |

This confirms three points of §2: the first-difference projections for each model (including the
regime-specific constants `DU` in Models 3 and 6), that the ADF regression has no deterministic terms,
and that the moments do not depend on level breaks.

---

## 1. Sources

| Source | Use |
|---|---|
| BCS (2015), *JAE* 30(1):1–23, `literature/BANERJEE-COINTEGRATIONPANELDATA-2015.pdf` (PDF page = paper page + 1, JSTOR cover) | **Specification**: models (§2, eqs 1–8), procedure (§3, eqs 9–12, MQ steps), Theorem 1, Table I (MQ critical values with one trend break), panel test (§4.1, eqs 13–14, Theorem 2, moments in Remark 1, Table II), unknown breaks (§4.2, eq. 15, Theorem 3, Table III), Monte Carlo design (§5, Table IV), application (§6, Table V) |
| BCS (2006), ECB Working Paper 591, "Cointegration in panel data with breaks and cross-section dependence", `literature/Banerjee-CarrioniSilvestre-ECBWP591-2006.pdf` (PDF page = paper page + 1) | Earlier version with the appendix: **A.6, endogenous regressors (DOLS)**; **Table 9, finite-sample moments** for T = 50…1000 (the JAE values are its T = 1000 row); Table 10 (= JAE Table I). Where it differs from the JAE version, the JAE version wins: e.g. its Table 11 Z*_τ critical values differ from JAE Table III, and its unknown-trend-break factor tests were dropped |
| BCS (2017), Birmingham DP 2017-02, `literature/bbs-dp2017-02.pdf` | Names the factor criterion ("panel BIC information criterion in Bai and Ng (2002)", r_max = 6); the "transformed" option (columns divided by their standard deviation before PCA); Table 16 applies Z to a PWT 6.3 production function (golden test T9, §11 Q7). Also the specification of the later `pcoint_banerjee_cce()` |
| Methods the paper delegates to, **in `literature/`** (§11 Q3) | Bai & Ng (2002), `Bai-Ng-NumberOfFactors-2002.pdf`: factor-number criteria, pp. 201–202 (PDF pp. 11–12). Bai & Ng (2004), `Bai-Ng-PANIC-2004.pdf`: PANIC steps and **MQ tests**, pp. 1133–1134; **MQ critical values in Table I**, p. 1136 (PDF pp. 7–10). Ng & Perron (1995), `Ng-Perron-TruncationLag-1995.pdf`: t-sig lag rule (scanned, no text layer) |
| Methods the paper delegates to, **not in `literature/`** (no free copy found; add if available) | Pesaran (2015), CD test; Bai & Perron (1998, 2003), break estimation. Neither is needed for v1. The CD statistic is the standard one (Pesaran 2004, 2015), and `plm::pcdtest()` (`plm` is already in Suggests) serves as the T3 cross-check. v1 estimates at most two breaks, by exact enumeration (B14). Also cited, not needed: Stock & Watson (1993), DOLS; Perron (1989), Model C |
| CRAN `xtbreakcoint` 1.0.4 (Mar 2026, GPL-3; credits Banerjee and Carrion-i-Silvestre for "original GAUSS code"; in the pinned nixpkgs) | **Information only.** It is not faithful to the paper (Appendix A), but its source shows two things about the GAUSS code: it cumulates with `cumsumc(De)`, and it uses T = 100 finite-sample moments (for γ, exactly the WP's Table 9 values). GPL-3: compare outputs, never copy code into this MIT package |

---

## 2. The method: precise specification

### 2.1 Model, hypotheses and the six models

For unit `i = 1..N`, `t = 1..T` (balanced; §3 B21) and `K` I(1) regressors (BCS eqs 1–8):

```
y_it = D_it + x_it' δ_it + u_it,     u_it = F_t' π_i + e_it,     (1 − ρ_i L) e_it = H_i(L) ε_it
D_it = μ_i + β_i t + Σ_j θ_ij DU_ijt + Σ_j γ_ij DT_ijt
```

**Break notation.**

- `T_b` is the last period of the old regime.
- `DU_t = 1(t > T_b)`, `DT_t = (t − T_b) 1(t > T_b)`, and the impulse is `D(T_b)_t = 1(t = T_b + 1)`.
- The break fraction is `λ = T_b / T`.
- In Models 4–6 the cointegrating vector shifts at the same dates as the deterministic part
  (footnote 2; also used in §6). Separate dates (Remark 1) are not supported in v1 (B12).

The factors `F_t` (r of them) may be I(0), I(1) or a mix. `r_1` is the number of I(1) factors.

| Model | Index j | Deterministic part | Cointegrating vector | Breaks allowed |
|---|---|---|---|---|
| 1 | c | constant, level shifts | stable | none, known or estimated; heterogeneous or common |
| 2 | τ | constant, trend, level shifts | stable | same as Model 1 |
| 3 | γ | constant, trend, level and trend shifts | stable | **known and common only** |
| 4 | c | constant, level shifts | shifts | known or estimated; heterogeneous or common |
| 5 | τ | constant, trend, level shifts | shifts | same as Model 4 |
| 6 | γ | constant, trend, level and trend shifts | shifts | **known and common only** |

Models 1 and 2 without breaks are the paper's "constant term" and "time trend" cases (Remark 3).

**Hypotheses.**

- H0: `ρ_i = 1` for all i (the idiosyncratic errors have a unit root, so there is no cointegration).
- H1: `ρ_i < 1` for some i. A rejection does not mean that every unit cointegrates (§4; Pesaran 2012).
- Z is left-tailed.

### 2.2 Per-unit projection in first differences (eqs 9–11)

The rows are `t = 2..T` (`n = T − 1`). For each unit, regress `Δy_i` on `Δx_i^d` by OLS and keep
the residuals `y*_i = M_i Δy_i`. `Δx_i^d` stacks deterministic and stochastic columns:

| Model | Deterministic columns (first differences) | Stochastic columns |
|---|---|---|
| 1 | `D(T_b,j)` impulses | `Δx` |
| 2 | `1`, impulses | `Δx` |
| 3 | `1`, impulses, `DU_j` | `Δx` |
| 4 | impulses | `Δx`, `Δx·DU_j` |
| 5 | `1`, impulses | `Δx`, `Δx·DU_j` |
| 6 | `1`, impulses, `DU_j` | `Δx`, `Δx·DU_j` |

Notes:

- The impulses also absorb the one-off term `x_{T_b}'(δ_2 − δ_1)` that `Δ(DU·x)` creates at
  `T_b + 1`.
- In Models 3 and 6 the `DU_j` columns enter `M_i` non-negligibly, so the breaks must be common to
  all units (paper p. 6). Otherwise `M_i ΔF` is no longer a common factor.
- With DOLS (§2.9) the stochastic columns also include `Δ²x_{t−j}`, `j = −q..p`.

### 2.3 Common factors (eq. 11; Bai & Ng 2002, 2004)

Let `Y*` be the `n × N` matrix `[y*_1 … y*_N]`.

- **Estimation.** `f̃` is `√n` times the eigenvectors of `Y*Y*'` for the r largest eigenvalues, so
  `f̃'f̃/n = I_r`. Then `π̃ = f̃'Y*/n` and `z̃ = Y* − f̃π̃'`. Compute with the smaller of `Y*Y*'`
  and `Y*'Y*` (§6.3).
- **Number of factors.** `r̂ ∈ {0..r_max}` minimises a Bai–Ng (2002) criterion. With μ_j the
  eigenvalues, let `V(k) = (Nn)⁻¹ Σ_{j>k} μ_j` and `σ̂² = V(r_max)`.
  - **Default BIC3**: `V(k) + k σ̂² (N + n − k) ln(Nn) / (Nn)`.
  - Also offered: `IC_p1..3` and `PC_p1..3`, with `C²_NT = min(N, n)`.
  - Checked against Bai & Ng (2002), pp. 201–202 (their T is our n). There
    `V(k) = N⁻¹ Σ_i ê_i'ê_i / T`, which equals the eigenvalue sum above, and `σ̂²` "can be replaced
    by `V(kmax)`". The `IC_p` criteria use `ln V(k)` and no `σ̂²`.
- **`r_max`.** Default 6 (§5), capped at `min(N, n) − 1`. A fixed `r` (including 0) may be given.
- **Optional standardisation.** Divide each column of `Y*` by its standard deviation before the PCA
  (BCS 2017's "transformed"; B5). The defactored `z̃` are then rescaled back. The ADF t-ratios are
  scale-invariant, so only the factor estimates change.

### 2.4 Idiosyncratic unit-root tests (eq. 12, Theorem 1.1)

**Cumulation.** `ẽ_it = Σ_{s=2..t} z̃_is` for `t = 2..T`, giving T − 1 values with no zero start (B8).

**ADF regression, without deterministic terms:**

```
Δẽ_it = α_i0 ẽ_i,t−1 + Σ_{j=1..k_i} α_ij Δẽ_i,t−j + ε_it
```

`t_i` is the conventional t-ratio of `α_i0`, using `s² = RSS/(n_eff − k_i − 1)`.

**Lag selection** (B6, B7):

- **t-sig** (Ng & Perron 1995, used by BCS §5 and §6): start at `k_max` and drop the last lag while
  its |t| < 1.645 (10% two-sided), down to `k_min` (default 0).
- Default `k_max = 4⌈(min(N, T)/100)^{1/4}⌉`.
- Candidates are compared on the common sample of the largest model. The chosen model is refitted on
  its own sample.

**Null limits** (Theorem 1.1, independent of breaks and of K under strict exogeneity):

| j | Limit of `t_i` |
|---|---|
| c | `½(W(1)² − 1) / (∫W²)^{1/2}` (Dickey–Fuller without constant) |
| τ | `−½ (∫V²)^{−1/2}`, with `V(s) = W(s) − sW(1)` |
| γ | `−½ (Σ_{j=1}^{m+1} (λ_j − λ_{j−1})² ∫_0^1 V_j(b)² db)^{−1/2}`, independent bridges, one per regime |

### 2.5 Panel statistic (eqs 13–14, Theorem 2)

```
SADF = Σ_i t_i,   t̄ = SADF / N
Z_j(λ) = (N^{−1/2} SADF − Θ_j √N) / √Ψ_j = √N (t̄ − Θ_j) / √Ψ_j  ⇒ N(0, 1),   p = Φ(Z)
```

**Moments:**

- `Θ_c = −0.424`, `Ψ_c = 0.964`; `Θ_τ = −1.535`, `Ψ_τ = 0.341` (Remark 1).
- γ: Table II for one break (λ ∈ {0.1, …, 0.9}) and two breaks (all 36 pairs). The upper triangle
  holds the mean and the lower triangle the variance.
- When λ is off the grid, or m ≥ 3, the moments are simulated (B10).
- Optional finite-sample moments: WP Table 9 (B9).

Theorem 2 needs the `e_i` to be independent across i, and `N/T → 0`. Both are documented, and the
CD diagnostic (§2.10) checks the first.

### 2.6 Common-factor tests (Theorem 1.2–1.3, MQ steps of §3)

These run only when `r̂ ≥ 1`.

**Detrending.** Let `F̃_t = Σ_{s=2..t} f̃_s` and regress `F̃` on:

- `1` (c);
- `(1, t)` (τ);
- `(1, t, DT_1..DT_m)` (γ);

giving `F̃^d`. No `DU` term is used (B18): `F̃` is continuous at the breaks because the impulses were
projected out in first differences, and the limit in Theorem 1.2 is the projection on
`χ_s = (1, s, dt_j(s))`.

**r̂ = 1.** ADF regression without deterministics on `F̃^d`, lags as in §2.4. The limit is DF_μ (c),
DF_τ (τ) or a λ-dependent distribution (γ; Perron 1989, Model C, for one break).

Bai & Ng (2004, p. 1133) instead run the factor ADF with an intercept on `F̂` itself, without a
detrending step. BCS (p. 7) describe two steps: regress `F̃` on the deterministics, then run the ADF
on `F̃^d`. We follow BCS. The two versions differ only in finite samples; both have Theorem 1.2's
limit.

**r̂ > 1: MQ tests**, run sequentially from `q = r̂` downwards.

1. `β̃⊥` = eigenvectors of the q largest eigenvalues of `T⁻² Σ F̃^d F̃^d'`, and `Ỹ_t = β̃⊥' F̃^d_t`.
2. **Non-parametric `MQ_c`.**
   - Take the residuals `ξ̃` of a VAR(1) in `Ỹ` without intercept.
   - Compute `Σ̃_1 = Σ_{j=1..J} (1 − j/(J+1)) T⁻¹ Σ_t ξ̃_{t−j} ξ̃_t'`.
   - `ṽ_c` is the smallest eigenvalue of `½[Σ(Ỹ_t Ỹ'_{t−1} + Ỹ_{t−1}Ỹ_t') − T(Σ̃_1 + Σ̃_1')] (Σ Ỹ_{t−1}Ỹ'_{t−1})⁻¹`.
   - `MQ_c = T(ṽ_c − 1)`.
3. **Parametric `MQ_f`.**
   - Fit a VAR(p) in `ΔỸ`, giving `Π̃(L)`, and filter: `ỹ_t = Π̃(L)Ỹ_t`.
   - `ṽ_f` is the smallest eigenvalue of `½[Σ(ỹ_t ỹ'_{t−1} + ỹ_{t−1}ỹ_t')] (Σ ỹ_{t−1}ỹ'_{t−1})⁻¹`.
   - `MQ_f = T(ṽ_f − 1)`.
4. If `H0: r_1 = q` is rejected, set `q ← q − 1` and repeat. Otherwise `r̂_1 = q`. If q = 1 is
   rejected, `r̂_1 = 0`.

**Settings.**

- `J = 4⌈(min(N,T)/100)^{1/4}⌉` (§5, footnote 17).
- p by BIC up to `p_max = 4⌈(T/100)^{1/4}⌉` (footnote 17).
- Level 5%.
- `T` here is the number of observations of `F̃` (B20).

**Critical values:**

- Models 1, 2, 4, 5: Bai & Ng (2004) Table I (p. 1136), using `MQ^c_{c,f}` (demeaned) for c and
  `MQ^τ_{c,f}` (detrended) for τ, for m = 1..6 at 1/5/10%. The same critical values serve `MQ_c` and
  `MQ_f`, because both have the limit `ν*(m)`.
- Models 3, 6: BCS Table I for one break with λ on the 0.1 grid; simulated otherwise (B19).

The result reports `r̂_1` from `MQ_c` and `MQ_f` separately, as BCS Table V does.

### 2.7 Breaks: none, known, estimated

- **None.** Models 1 and 2 only.
- **Known.** Given as dates.
  - Units may have different dates and different numbers of breaks (Remark 1) in Models 1, 2, 4
    and 5.
  - Models 3 and 6 need one common date vector.
  - No trimming is applied, but every regime must leave the projection of full rank.
- **Estimated**, Models 1, 2, 4 and 5 only (§4.2), with `m ∈ {1, 2}` breaks per unit (B14):
  - **Heterogeneous (§4.2.1):** unit-specific dates by the iterative procedure of §2.8. Then compute Z
    as with known breaks; the moments do not depend on the dates.
  - **Homogeneous (§4.2.2, eq. 15):** compute `Z_j(λ)` with the full procedure of §§2.2–2.5 for
    every admissible common break vector λ ∈ Λ^m (B16).
    - `Z*_j = min_λ Z_j(λ)` and `λ̂ = argmin` (ties go to the earliest date).
    - By Theorem 3, `Z*_j ⇒ N(0,1)`. Table III gives finite-sample critical values (m = 1, N = 100,
      T ∈ {50, 100, 250}; B17).
    - The factor tests and the unit table are reported at `λ̂`.
  - **Trimming** (B13): `Λ = [0.15, 0.85]` (footnote 1). Every regime spans at least
    `h = max(⌊0.15 T⌋, k_reg + 2)` periods, where `k_reg` is the number of regime-specific
    coefficients.

### 2.8 Heterogeneous break estimation (§4.2.1)

Coordinate descent on the joint sum of squares
`S = Σ_i ‖Δy_i − ΔD_i(T_b,i) − ΔX_i(T_b,i)δ_i − ΔF π_i‖²`:

1. **s = 1.** For each i, choose `T_b,i ∈ Λ^m` and `δ_i` to minimise the SSR of `Δy_i` on
   `Δx_i^d(T_b,i)`; the factors sit in the error, since ΔF has mean zero (footnote 8). Keep the
   residuals `Δũ^(1)`.
2. **PCA.** Run principal components on `Δũ^(s)` with `r` factors, giving `ΔF̃^(s)` and `π̃^(s)`. Set
   `S^(s) = Σ_i ‖Δũ_i^(s) − ΔF̃^(s)π̃_i^(s)‖²`.
3. **Re-estimate breaks.** For each i, re-estimate `T_b,i` and `δ_i` on
   `Δy_i − ΔF̃^(s)π̃_i^(s)`. The residuals `Δũ^(s+1) = Δy_i − ΔD̂_i − ΔX_iδ̂_i` exclude the factor
   term.
4. **Stop or repeat.** Stop when `(S^(s−1) − S^(s)) / S^(s−1) < tol`. Otherwise repeat from step 2.

Each step minimises `S` over one block of parameters, so with `r` held fixed `S` cannot increase and
the stopping rule is well defined (B15). The default fixes `r` by the criterion at the first PCA. The
final Z re-selects `r` on `Y*` at the estimated dates, as in the known-break case.

The paper warns that in Model 2 the break dates are consistent only for large breaks (§6.3.2). This
goes in the documentation.

### 2.9 Endogenous regressors: DOLS (Assumption G(ii); WP Appendix A.6)

Under non-strict exogeneity, the limits of §2.4 hold if the first-difference regression includes leads
and lags of the regressors' second differences: `Δx_i^d = (Δx', Δ²x'_{t−j}, j = −q..p)`.

- The orders are fixed, as in Stock & Watson (1993) (the WP adds that BIC could choose them) and are
  common to all units, so that the rows of `Y*` stay aligned.
- The rows become `t = 3 + p .. T − q`.
- The default is no DOLS terms, as in the paper's exposition (B22).

### 2.10 Diagnostics

1. **Pesaran (2015) CD test** of cross-section independence, on the ADF residuals `ε̂_it` over the
   periods common to all units: `CD = √(2/(N(N−1))) Σ_{i<j} √T_ij ρ̂_ij`, two-sided N(0,1).
   Footnote 16 cites the WCD statistic of "Pesaran (2013)", the working-paper version of
   Pesaran (2015), and asks for two-sided inference. The paper reports it in Table V to justify
   Theorem 2.
2. **Share of units with `t_i` below their own 5% critical value** (Table V, "% individual
   rejections"; B24).
3. **Iteration log** (estimated breaks): S per iteration, convergence flag, `r` used.

---

## 3. Gaps, ambiguities and design decisions

Defaults follow the JAE text. Where it is silent, the WP or BCS (2017) decide if they can; otherwise
the default is my recommendation. Every choice can be changed through `banerjee_control()` or an
argument.

| # | Item | Paper | Alternatives seen | **Our default** | Notes |
|---|---|---|---|---|---|
| B1 | Which version | JAE 2015 | WP 2006 differs in places | **JAE**; WP only where JAE is silent | e.g. Z*_τ critical values: JAE Table III, not WP Table 11 |
| B2 | Model numbering | six models (§2) | xtbreakcoint: five different models (App. A) | **paper's 1–6** | `model` is required, as `lags` is in `pcoint_westerlund()`: the choice matters |
| B3 | First-difference regression | unit by unit (`M_i`) | xtbreakcoint pools δ over units | **unit by unit** | |
| B4 | Factor-number criterion | "panel BIC in Bai and Ng (2004)" (§5) | WP and BCS 2017: "panel BIC in Bai and Ng (2002)"; xtbreakcoint: IC_p2 | **BIC3 of Bai & Ng (2002), r_max = 6** | BIC3 is the only "BIC" in Bai & Ng (2002) that uses both N and T (p. 202; formula checked, §2.3). Bai & Ng (2004) define no BIC; their simulations use `IC_p1` from the 2002 paper (p. 1145). So "Bai and Ng (2004)" in the JAE is read as a citation slip. Argument `factors` |
| B5 | Scaling before PCA | none | BCS 2017 "transformed" (divide by s.d.) | **none** | control `pca_scale` |
| B6 | ADF lag rule | t-sig; `k_max = 4ceil[min[N,T]/100]^{1/4}` (§5), `12ceil[…]` in §6 | GAUSS port: BIC with k_max = 4 | **t-sig at 10%; `k_max = 4⌈(min(N,T)/100)^{1/4}⌉`** | The printed formula is ambiguous. Readings that apply `ceil` before multiplying by 4 give 4 for every min(N, T) ≤ 100. `⌈4(·)^{1/4}⌉` also gives 4 at the paper's min(N, T) = 40, but 3 below 32; a Schwert-style `⌊4(·)^{1/4}⌋` gives 3 at 40. Control `max_lags`, `lag_criterion = "tsig" / "bic" / "aic"`, `tsig_level` |
| B7 | ADF selection sample | not stated | each candidate on its own sample | **common sample of `k_max`, refit on own sample** | as Westerlund D7. Control `lag_sample` |
| B8 | Start of the cumulation | `ẽ_t = Σ_{s=2..t} z̃_s` | prepend `ẽ_1 = 0` | **T − 1 values, no zero row** | what the GAUSS code does (`cumsumc(De)` per xtbreakcoint) |
| B9 | Moments | asymptotic (T = 1000 simulation) | WP Table 9 by T. The GAUSS code (per xtbreakcoint) uses T = 100 values: γ equals WP Table 9's T = 100 column; c and τ are close to it | **asymptotic** | Control `moments = "finite_sample"`: WP Table 9, linear in 1/T, clamped to [50, 1000]; c, τ and one-break γ only |
| B10 | γ moments off the Table II grid | Table II only (λ ∈ 0.1 grid, m ≤ 2) | xtbreakcoint: nearest grid point | **Table II if every λ is on the grid (to 1e-8) and m ≤ 2; otherwise simulate at the exact λ** with the paper's design (1000 steps, 100,000 draws) and a fixed internal seed | deterministic and leaves R's RNG untouched. With T = 100 and one break, 90% of admissible dates are off the grid; the nearest bin can be off by 0.07 in Θ, which shifts Z by up to 0.7 at N = 40 |
| B11 | Break-date convention | `T_b` = last period of the old regime | strucchange uses the same; users may think "first new period" | **paper's**: `breaks = 1998` means the new regime starts in 1999 | printed explicitly ("break after 1998") |
| B12 | Separate dates for deterministic and cointegrating breaks | allowed (Remark 1), not used | — | **same dates** (footnote 2, §6) | later extension |
| B13 | Trimming | `Λ = [0.15, 0.85]` (footnote 1) | — | **`trim = 0.15`, minimum regime length `h`** (§2.7) | |
| B14 | Number of estimated breaks | any m; application uses ≤ 2; Table II covers ≤ 2 | Bai–Perron dynamic programming for larger m | **m ∈ {1, 2}, exact global minimum by enumeration** (§6.3) | larger m: known breaks only (§11 Q8) |
| B15 | r during the iteration | not stated | xtbreakcoint re-selects every iteration | **fixed after the first PCA**; final Z re-selects | keeps the SSR monotone. Control `iter_tol = 1e-6` (relative), `iter_max = 100`, warning if not converged |
| B16 | r and lags in the homogeneous search | not stated | — | **re-selected for every candidate λ** (Z_j(λ) is defined by the full procedure) | a fixed `factors = k` or `lags = k` fixes them |
| B17 | Inference for Z* | N(0,1) (Theorem 3); Table III for one break, N = 100 | — | **p-value from N(0,1), plus Table III critical values for the nearest T** | Table III's 5% values (−1.93 to −2.24) are far below −1.645, so the asymptotic p-value is anti-conservative; print says so. §11 Q4 |
| B18 | Factor detrending in Models 3, 6 | Theorem 1.2: `χ_s = (1, s, dt_j)`; Remark 5 cites Perron (1989) Model C, which has a level dummy | — | **`(1, t, DT_j)`, no DU** | follows the theorem; F̃ has no level jumps (§2.6). Critical values simulated from the theorem's functional |
| B19 | Distributions not tabulated in the paper | — | — | **simulated by the package**: DF quantiles (no constant, μ, τ), τ and γ idiosyncratic quantiles, γ factor ADF quantiles, MQ for γ off-grid or m = 2 | Published values are used where they exist and checked against our simulations (T1) |
| B20 | MQ details | as in §2.6 | WP prints the VAR in `Ỹ`, not `ΔỸ`; Bai & Ng (2004) use `ΔỸ` | **JAE and Bai–Ng**: VAR in `ΔỸ`; `T` = number of `F̃` observations; VAR(1) residuals without intercept | When r̂ = 1, ADF only (paper) |
| B21 | Unbalanced panels | not covered | — | **error**; PCA needs a common time axis | §11 Q9 |
| B22 | DOLS terms | recommended under endogeneity; exposition assumes strict exogeneity | — | **`dols = 0`** | `dols = c(lags, leads)` adds `Δ²x` terms (§2.9) |
| B23 | Residuals used for CD | `ẽ` "for different values of k" (Table V) | — | **ADF residuals of the selected models, over common periods** | |
| B24 | Individual rejections | Table V share at 5% | xtbreakcoint: −1.95 for every model | **each model's own 5% critical value** | −1.95 is right only for c |
| B25 | DOLS argument | not stated | — | **`dols = 0`: no terms; `dols = c(p, q)` adds `Δ²x_{t−j}`, `j = −q..p`, including `j = 0`** | the usual DOLS set in first differences; a scalar sets both orders |
| B26 | Names of the break columns | — | — | **`break_date` (one break), `break_date_1`, `break_date_2`** | `break` is a reserved word in R |
| B27 | PCA inside the break iteration (§2.8) | not stated | — | **unscaled, even with `pca_scale = TRUE`** | keeps `S` monotone; the final Z uses `pca_scale` |
| B28 | `mq_level` | 5% | — | **1%, 5% or 10% only** | the levels of the published tables |
| B29 | Default `r_max` in small panels | 6 | — | **`min(6, min(N, n) − 1)`; an explicit `max_factors` that is too large is an error** | with `r_max = N − 1`, BIC3 tends to pick many factors (σ̂² = V(r_max) is small); documented |
| B30 | A single `breaks` value of 1 or 2 | — | — | **number of breaks to estimate**, even if 1 or 2 is also a time value; other values are dates | give unit-specific dates as a data frame to place a break after period 1 or 2 |
| B31 | Bai–Ng Table I misprint (q = 5, detrended, 10%: −55.286) | — | our simulation −52.61 | **kept as printed** | flagged in `R/banerjee-tables.R` and RESULTS.md; affects `mq_level = 0.10` with ≥ 5 factors only |

---

## 4. User-facing API

```r
pcoint_banerjee(
  formula,                     # y ~ x1 + x2, levels of the I(1) series
  data,                        # data.frame / tibble / data.table / plm::pdata.frame (balanced)
  index = NULL,                # as pcoint_westerlund()
  model,                       # required: 1..6, the paper's models (section 2.1)
  breaks = NULL,               # NULL: no breaks (Models 1, 2)
                               # 1 or 2: number of breaks to estimate per unit
                               # vector of time values: known breaks common to all units
                               # data.frame(id, break) or named list: known unit-specific breaks
  break_type = c("heterogeneous", "homogeneous"),   # estimated breaks only
  factors = "bic3",            # criterion ("bic3", "ic1".."ic3", "pc1".."pc3") or a fixed number (0 allowed)
  max_factors = 6L,
  lags = NULL,                 # ADF lags: NULL = t-sig from the k_max rule; scalar fixed; c(min, max)
  dols = 0L,                   # leads and lags of D2 x: scalar (same) or c(lags, leads)
  trim = 0.15,
  control = banerjee_control(),
  n_threads = getOption("panelcointr.threads", 1L)
)

banerjee_control(
  lag_criterion = "tsig",      # "tsig" | "bic" | "aic"
  tsig_level = 0.10,
  lag_sample = "common",       # "common" | "varying"
  max_lags = NULL,             # NULL: 4 * ceiling((min(N, T) / 100)^(1/4))
  pca_scale = FALSE,
  moments = "asymptotic",      # "asymptotic" | "finite_sample"
  sim_reps = 100000L,          # simulated moments and critical values (B10, B19)
  sim_steps = 1000L,
  iter_max = 100L,
  iter_tol = 1e-6,
  mq_level = 0.05,
  mq_bandwidth = NULL,         # NULL: 4 * ceiling((min(N, T) / 100)^(1/4))
  mq_max_var = NULL,           # NULL: 4 * ceiling((T / 100)^(1/4))
  min_df = 1L
)
```

**Validation of the break arguments.**

- Models 4–6 need breaks.
- Models 3 and 6 need known common breaks.
- Estimated breaks need Model 1, 2, 4 or 5 and `breaks` ∈ {1, 2}.
- Known dates must be time values of the panel, and each regime must leave the projection of full
  rank.
- `break_type` with known breaks is an error, so a vector cannot be silently misread.

**Input rules** come from `prepare_panel()`, plus three checks:

1. After the leading and trailing NAs are trimmed, all units must cover the same calendar window. An
   unbalanced panel is an error that lists the units and their windows.
2. `max_factors < min(N, n)`.
3. Each unit's largest projection and ADF regression must have at least `min_df` residual degrees of
   freedom.

**Return value**, class `c("pcoint_banerjee", "pcoint_test")`:

| Field | Content |
|---|---|
| `statistics` | as for Westerlund: one row, `Z_c` / `Z_tau` / `Z_gamma`, or `Zinf_c` / `Zinf_tau` for homogeneous estimated breaks. `value` = `t̄`, `z` = Z, `p_value` = Φ(Z), `p_boot` = NA, `alternative` = "rho_i < 1 for some i" |
| `units` | per unit: `id, lags, n, t_stat, cv_5, reject_5`, break dates (one column per break) and `lambda` |
| `factors` | `r`, criterion values for `0..r_max`, eigenvalues, `f` (n × r) and `F` (cumulated), loadings. With r = 1: `adf` (statistic, lags, critical values, p-value). With r > 1: `mq` (table of q, `MQ_c`, its critical value, `MQ_f`, its critical value, VAR order), `r1_c`, `r1_f` |
| `breaks` | `type` ("none", "known", "heterogeneous", "homogeneous"), `m`, an N × m matrix of dates (time values), λ. Homogeneous: `profile` (data.frame of candidate dates, t̄, Z, r̂). Heterogeneous: `iterations`, `converged`, `ssr` path |
| `moments` | `Θ`, `Ψ` and their source ("Remark 1", "Table II", "WP Table 9", "simulated") |
| `diagnostics` | CD statistic and p-value; share of individual rejections |
| `settings`, `N`, `K`, `T`, `model`, `call`, `formula` | |

**Methods.**

- `print()`: the statistics table (`format_statistics()`), then factors, breaks, the CD line and a
  note on interpretation (§0). For Z*, the Table III critical values are shown too.
- `summary()`: adds the unit table, the criterion values, the MQ sequence and the break-date
  distribution.
- `plot()`: the Z(λ) profile (homogeneous); otherwise the estimated factors `F̃` with the break
  dates marked.
- `as.data.frame()`, `tidy()`, `glance()`: as for Westerlund.

```r
pcoint_banerjee(y ~ x1 + x2, data = d, index = c("id", "year"), model = 2)          # no breaks
pcoint_banerjee(mp ~ er + fp, data = d, model = 5, breaks = 2)                       # BCS section 6
pcoint_banerjee(y ~ x, data = d, model = 2, breaks = 1, break_type = "homogeneous")  # BCS section 5
pcoint_banerjee(y ~ x, data = d, model = 6, breaks = 1998)                           # known common break
```

---

## 5. Architecture

```
R  pcoint_banerjee()
   ├─ prepare_panel() + require_balanced()        [shared]
   ├─ resolve_breaks()      time values ↔ positions T_b; validation against the model (§4)
   ├─ resolve_settings()    banerjee_control() + defaults → typed list
   ├─ banerjee_rs()         ── Rust:
   │     none/known:   projections → PCA + criterion → cumulate → ADF (t-sig) → t̄; factor ADF / MQ
   │     heterogeneous: iteration of §2.8 → known-break procedure at the estimates
   │     homogeneous:   procedure at every candidate λ (threads) → profile → procedure at λ̂
   ├─ moments(): Remark 1 / Table II / WP Table 9 literals, or banerjee_null_rs() (fixed seed)
   ├─ critical values and p-values: literals (BCS Tables I, III; Bai-Ng Table I) and sysdata
   │  quantile tables, or banerjee_null_rs() for γ off the grid
   ├─ CD diagnostic         (R: N × N correlations of the ADF residuals)
   └─ new_pcoint_banerjee()
```

**Data passed to Rust.**

- `y`, `X`, `offsets` as for Westerlund. The panel is balanced, so `T` is common.
- Break positions as an `N × m_max` integer matrix, −1 padded.
- `model`, the lag and factor settings, `dols`, trimming, `n_threads`.

**Data returned.** A named list of plain vectors and matrices:

- per-unit `t`, `k`, `n`;
- `t̄`;
- eigenvalues, criterion values, `r̂`, `f̃`, `π̃`;
- the ADF residuals (`n × N`, NaN-padded, for CD);
- factor ADF and MQ statistics;
- break estimates, the iteration log and the homogeneous profile.

R standardises Z, because the moments may come from tables or from simulation.

---

## 6. Rust engine

### 6.1 Dependencies and packaging

**No new crates.**

- The eigensolver, random number generator and normal sampler are hand-written: about 350 lines,
  tested against R.
- The MSRV (1.71), `Cargo.lock` and the vendored tarball are unchanged, so §7 of the Westerlund plan
  needs no new work.
- Quantile tables (B19) are generated once by `data-raw/banerjee-tables.R`, which calls the package's
  own simulator, and stored in `R/sysdata.rda`.

### 6.2 Modules (`src/rust/src/core/`)

```
eigen.rs     symmetric eigendecomposition (Householder tridiagonalisation + implicit QL, as EISPACK
             tred2/tql2); generalised symmetric A v = μ B v (B SPD) through Cholesky, for MQ
factors.rs   PCA of an n × N matrix through the smaller Gram matrix; Bai-Ng criteria; sign normalisation
adf.rs       ADF regression with or without deterministics; t-sig / IC lag selection from one QR
breaks.rs    first-difference designs per model; break-date search (prefix-sum Gram); admissible sets
banerjee.rs  procedure: projections, cumulation, unit ADF, aggregation, iteration (2.8), homogeneous grid
mq.rs        factor detrending, factor ADF, MQ_c, MQ_f
rng.rs       xoshiro256++ seeded by SplitMix64; Box-Muller normals
nulldist.rs  simulation of the Theorem 1 functionals (idiosyncratic c/τ/γ, factor DF, MQ) for given λ
```

Reused: `linalg::Qr` (all regressions), `linalg::cholesky` and `forward_solve`, and `panel::Panel`.

### 6.3 Algorithms

1. **Every t-sig candidate from one QR.**
   - Order the ADF columns `[ẽ_{t−1}, Δẽ_{t−1}, …, Δẽ_{t−k_max}]` on the common sample.
   - With 0-based column indices, the last coefficient of the nested model with k lags is
     `(Q'y)_k / R_kk`, with standard error
     `s_k / |R_kk|`, where `s_k² = (y'y − Σ_{j≤k}(Q'y)_j²) / (n − k − 1)`.
   - So one QR gives the t-statistic every t-sig step needs, as Westerlund's §6.3 item 1 did for the
     panel statistics.
   - The chosen model is refitted on its own sample with `ẽ_{t−1}` last, so `t_i = c/r` from the
     last diagonal entry.
2. **Principal components.**
   - Factorise the smaller of `Y*'Y*` (N × N) and `Y*Y*'` (n × n). All eigenvalues are needed by the
     criteria; the vectors only for the top r_max.
   - Through the N × N route: `f̃ = √n Y* V_r M_r^{−1/2}`, `π̃ = V_r M_r^{1/2} / √n` and
     `z̃ = Y*(I − V_r V_r')`.
   - Tridiagonal QL costs `O(m³)` for `m = min(N, n)`, a few ms for m = 300.
   - Signs are normalised so each eigenvector's largest-|·| element is positive. The statistics do not
     depend on the sign; this keeps the reported factors deterministic.
3. **Break search with prefix-sum Gram matrices** (heterogeneous estimation).
   - An impulse at `T_b + 1` fits that row exactly, so it is equivalent to **deleting the row**.
   - Each candidate's normal equations are therefore assembled from prefix sums of `w_t w_t'` and
     `w_t Δy_t`, with `w_t = (1, Δx_t, Δ²x terms)`: whole-sample sums for common coefficients,
     segment sums for regime-specific slopes (Models 4, 5), minus the deleted rows.
   - Assembly is `O(m k₀²)` and Cholesky `O(k³)` per candidate, after an `O(n k₀²)` pass.
   - m = 2 enumerates `O(n²/2)` candidates exactly: 2,485 for T = 123 with 15% trimming.
   - Columns are standardised; the winner is refitted by QR, with a QR fallback when a Cholesky pivot
     falls below `CHOL_TOL` (as in Westerlund's IC grid).
   - In the iteration, only the `Δy`-side sums change (the factor component is subtracted), so the
     `w`-side sums are cached per unit.
4. **Homogeneous grid.**
   - Each candidate runs the full known-break procedure: QR projections (not Gram, because the
     residuals feed the statistic), PCA and ADF.
   - Candidates are independent, so they are split across threads.
   - The result is the `t̄(λ)` and `r̂(λ)` profile; R takes the minimum of Z. Θ and Ψ do not depend
     on λ for c and τ, so the argmin of Z is the argmin of t̄.
5. **Null distributions.**
   - Each statistic is a functional of discretised Brownian motions with `sim_steps` steps:
     - idiosyncratic: Riemann sums of the bridge integrals;
     - factor ADF: the DF statistic of the detrended walk;
     - MQ: the eigenvalue problem of Theorem 1.3.
   - Replicates are split into fixed-size blocks, each with its own stream from
     `SplitMix64(seed, block)`, so the results do not depend on `n_threads`.
   - Moments and critical values computed at run time (B10, B19) use a **fixed seed**: the output is a
     deterministic function of λ, and R's RNG is not touched.
6. **Numerical safety.**
   - Rank deficiency in any projection (for example a regime with too few distinct `Δx` values, or
     known breaks too close together) is an error naming the unit and the break dates.
   - Non-finite or non-positive eigenvalues are clamped to 0 with a check that their sum matches
     `‖Y*‖²` to 1e-10 relative.
   - The MQ denominators must be positive definite.

### 6.4 Threads

`n_threads` works as for Westerlund. Parallelised:

- the homogeneous candidates;
- the per-unit break search in each iteration;
- the null simulations.

Each task writes to its own output slot, and the random streams are tied to blocks, so the results are
bit-identical for any `n_threads`.

---

## 7. Packaging and environment

- **Nothing changes in the Rust packaging** (§6.1).
- **New R files:**
  - `R/banerjee.R`, `R/banerjee-control.R`, `R/banerjee-methods.R`;
  - `R/banerjee-tables.R` (literals with page citations: Remark 1, Tables I–III, WP Table 9 and
    Bai–Ng Table I);
  - `R/sysdata.rda` (simulated quantile tables) and `data-raw/banerjee-tables.R`.
- **CLAUDE.md:** add the rule that `tests/testthat/helper-banerjee-reference.R` is the oracle and new
  settings must be added to it and to the T2 comparisons.
- **home-manager:** `xtbreakcoint` (1.0.4 in the pinned nixpkgs) is needed by the T8 script only. Add
  it to `r_packages` at P7. It does not go into DESCRIPTION.
- **Possible dataset** (§11 Q7): if T9 passes and the terms allow it, `data/` gets the PWT
  production-function panel. It is built by `data-raw/pwt-production.R`, documented in
  `R/data.R`, and covered by `LazyData: true`, which is already set.

---

## 8. Validation and testing

| ID | What | How | Tolerance |
|---|---|---|---|
| T1 | **Published tables by simulation**: (a) Remark 1 moments; (b) Table II, all 9 + 36 entries; (c) Table I, MQ with one trend break; (d) WP Table 9; (e) the asymptotic DF critical values used for the factor ADF; (f) Bai & Ng (2004) Table I | `inst/validation/banerjee-tables.R` with the paper's design (1000 steps, 100,000 draws); the literals in `R/banerjee-tables.R` are checked against it | within Monte Carlo error (z-scores of the differences, as in Westerlund T7(c)) |
| T2 | **Rust vs pure-R oracle** | `helper-banerjee-reference.R`: `lm.fit` projections, `eigen()` PCA, a loop of `lm` fits for t-sig, brute-force `lm` over every break candidate, a naive iteration and grid. Configurations: models 1–6 × breaks (none, known heterogeneous, known common, estimated heterogeneous, estimated homogeneous) × K ∈ {1, 2, 3} × r fixed or selected × lags fixed or t-sig × DOLS on/off. Panels with near-ties are discarded: |t| within 1e-6 of the t-sig cut-off, criterion ties within 1e-8, SSR ties within 1e-8 relative | 1e-10 relative |
| T3 | Kernels (`cargo test` plus R tests through internal wrappers, like `ols_rs`) | eigen vs `eigen()` (random, clustered and repeated eigenvalues, n ≤ 300); generalised eigen; PCA vs `svd()` on rotation-invariant quantities; ADF vs `lm`; one-QR t-sig vs separate fits; Gram break search vs brute force; MQ vs direct R formulas; RNG against published xoshiro256++ test vectors; CD vs `plm::pcdtest(test = "cd")` on the same residual matrix (`skip_if_not_installed("plm")`) | 1e-12 (CD: 1e-10) |
| T4 | Invariances, all models and break regimes | rescaling y and each x; permuting units; adding unit constants (all models); adding unit trends (τ, γ); adding level shifts at known break dates; adding trend shifts at known common dates (3, 6); flipping eigenvector signs; `n_threads` 1 vs 2 | 1e-10 (threads: exact) |
| T5 | Determinism | repeated calls identical; simulated moments identical across runs and thread counts; `.Random.seed` unchanged by `pcoint_banerjee()` | exact |
| T6 | Input validation | unbalanced panels; dates outside the sample or inside the trimming; Models 3 and 6 with unit-specific or estimated breaks; `breaks = 3` estimated; `max_factors ≥ min(N, n)`; collinear regimes; `break_type` with known breaks | exact messages |
| T7 | Monte Carlo (`inst/validation/`, not in `R CMD check`) | (a) size of Z under H0 for every model and break regime, BCS §5 DGP with N = 40; (b) Table IV: Models 2 and 5, one estimated common break, r = 1, the grid over σ²_F, φ, ρ and T; (c) Table III critical values (N = 100); (d) MQ selection frequencies against the §5 discussion | within Monte Carlo error |
| T8 | Information only | `xtbreakcoint` on Models 1 and 2 without breaks, with r and lags fixed; differences are expected (Appendix A) and are documented, not tested | none |
| T9 | **Golden test, if the data can be obtained** (§11 Q7) | BCS (2017) Table 16: Model 1, no breaks, PWT 6.3, N = 19, T = 57; lag setting unknown (t-sig default first, then fixed k). Z for r = 0..6, with and without `pca_scale`; `r̂_1` from `MQ_c` and `MQ_f`; the BIC3 choice. In `inst/validation/`, and also in testthat if the data ship with the package | Z to the printed 2 decimals; `r̂_1` and BIC3 exact. A mismatch is investigated and documented, not forced |

---

## 9. Performance targets

Single thread, release build, `bench::mark` in `inst/bench/banerjee.R`:

| Workload | Target |
|---|---|
| Known breaks, Model 5, m = 2, N = 63, T = 123, K = 2 (the size of the paper's application) | < 10 ms |
| Heterogeneous estimation, same panel, m = 2 | < 1 s |
| Homogeneous, Model 2, m = 1, N = 40, T = 250 | < 0.5 s |
| Homogeneous, Model 5, m = 2, N = 40, T = 100 | < 3 s; near-linear speed-up with threads |
| γ moments by simulation (100,000 × 1000) | < 2 s |
| `data-raw` quantile tables (10⁶ draws per distribution) | minutes; offline only |

---

## 10. Implementation phases

Each phase ends with `devtools::check()` clean.

| Phase | Deliverable | Acceptance |
|---|---|---|
| P0 Sources | Sources are already in `literature/` (§11 Q3). Transcribe Remark 1, Tables I–III, WP Table 9 and Bai–Ng Table I into `R/banerjee-tables.R` with page citations (Appendix B) | literals double-checked against the PDFs |
| P1 Shared R | `require_balanced()`, `resolve_breaks()` (dates ↔ positions, validation), control and settings resolution | T6 (partial) |
| P2 R oracle | `helper-banerjee-reference.R` for every model and break regime | oracle reproduces §0 and T1(a)/(b) at small scale |
| P3 Rust kernels | `eigen`, `adf`, `factors`, `breaks`, `rng`, plus internal R wrappers for testing | T3 |
| P4 Rust procedure | none/known → heterogeneous iteration → homogeneous grid → factor ADF and MQ; `pcoint_banerjee()` with table moments | T2, T4, T5 |
| P5 Null distributions | `nulldist`, `banerjee_null_rs()`, `data-raw/banerjee-tables.R` → `sysdata.rda`; γ moments off the grid; unit and factor p-values; Table III lookup | T1 |
| P6 Interface | print/summary/plot/tidy/glance, roxygen, vignette (§5 DGP with a break, or the PWT data if T9 passes), NEWS, README. Obtain the BCS (2017) PWT 6.3 panel and check its terms (§11 Q7) | docs build |
| P7 Validation | T7, T8, T9 (if the data were obtained), §9 benchmarks, `inst/validation/RESULTS.md` section, CLAUDE.md update | targets met; T9 outcome recorded either way |

---

## 11. Decisions on the open questions

All recommendations were accepted on 2026-10-12. The numbering is kept, so the `§11 Qn` references
elsewhere still point here.

1. **Which paper: the 2015 JAE article.** The 2017 CCE paper becomes a separate function
   (`pcoint_banerjee_cce()`) with its own plan, after this one. It reuses the ADF and simulation
   kernels.
2. **Function name: `pcoint_banerjee()`**, after the first author, like `pcoint_westerlund()`.
3. **Sources: added to `literature/` on 2026-10-12** (§1): ECB WP 591, Bai & Ng (2002), Bai & Ng
   (2004) and Ng & Perron (1995), plus the CLAUDE.md entries with page offsets. Two checks were
   done straight away:
   - The BIC3, `IC_p` and `PC_p` formulas match §2.3.
   - The MQ steps match §2.6 and B20. Bai & Ng's Table I (p. 1136) gives the no-break MQ critical
     values, so P0 transcribes them rather than simulating.

   Pesaran (2015) and Bai & Perron (1998) have no free copy. Neither is needed for v1 (§1).
4. **Inference for Z*: an N(0,1) p-value plus the Table III critical values (B17).** The print
   method warns that the normal p-value is anti-conservative.
   - Deferred: a finite-sample p-value from simulating the whole procedure under H0 at the user's
     N and T. Each draw repeats the grid search, which is feasible in Rust. Add it if users ask.
5. **Moments: asymptotic by default** (Remark 1, Table II).
   - `moments = "finite_sample"` uses WP Table 9 (printed p. 40, PDF p. 41).
   - For γ this reproduces what the authors' GAUSS code does; for c and τ it comes close (B9).
6. **`k_max = 4⌈(min(N,T)/100)^{1/4}⌉`** (B6). It is the closest reading of the typography and gives
   k_max = 4 in the paper's own simulations. `max_lags` overrides it.
7. **Golden data: try BCS (2017) Table 16.**
   - The test: Model 1 without breaks, `y = α_i + β_1 l + β_2 k` (logs of GDP per capita,
     population and capital per capita).
   - The data: the PWT 6.3 panel of Banerjee, Eberhardt & Reade (2010), 19 OECD countries,
     1951–2007 (N = 19, T = 57).
   - The targets:
     - Z for fixed r = 0..6, with and without standardisation (`pca_scale`);
     - `r̂_1` from `MQ_c` and `MQ_f`;
     - BIC3 picking r = 6 without standardisation. With standardisation the DP (p. 19) reports
       "only one integrated common factor", read as r̂ = 1.
   - The work:
     - P6: obtain the data and check the licence. Current PWT releases are CC BY 4.0, but 6.3
       predates that, and the capital series is constructed by the authors; check both.
     - P7: reproduce the targets as T9 (§8). Ship the data as a package dataset (§7) only if T9
       passes and the terms allow it; otherwise the vignette uses simulated data.
     - The DP does not state the ADF lag settings behind Table 16. So T9 tries the default
       t-sig rule first, then fixed k = 0..4, and records which setting matches.
8. **More than two estimated breaks: not in v1** (B14). Bai–Perron dynamic programming is exact for
   Model 4 and iterative for Models 1, 2 and 5. Add it if users ask; Bai & Perron (1998) would then
   need to go into `literature/`.
9. **Unbalanced panels: not in v1** (B21). An unbalanced panel is an error. EM-based PCA (Stock &
   Watson 2002) is the usual extension.
10. **Long-run relation (eq. 20): out of scope.** It is an estimator, not a test.

---

## 12. References

- Banerjee, A. & Carrion-i-Silvestre, J. L. (2015). Cointegration in panel data with structural breaks and cross-section dependence. *Journal of Applied Econometrics* 30(1), 1–23. doi:10.1002/jae.2348
- Banerjee, A. & Carrion-i-Silvestre, J. L. (2006). Cointegration in panel data with breaks and cross-section dependence. ECB Working Paper 591.
- Banerjee, A. & Carrion-i-Silvestre, J. L. (2017). Testing for panel cointegration using common correlated effects estimators. Birmingham Business School DP 2017-02.
- Bai, J. & Ng, S. (2002). Determining the number of factors in approximate factor models. *Econometrica* 70(1), 191–221.
- Bai, J. & Ng, S. (2004). A PANIC attack on unit roots and cointegration. *Econometrica* 72(4), 1127–1177.
- Bai, J. & Perron, P. (1998). Estimating and testing linear models with multiple structural changes. *Econometrica* 66(1), 47–78.
- Ng, S. & Perron, P. (1995). Unit root tests in ARMA models with data-dependent methods for the selection of the truncation lag. *JASA* 90, 268–281.
- Pesaran, M. H. (2012). On the interpretation of panel unit root tests. *Economics Letters* 116, 545–546.
- Pesaran, M. H. (2015). Testing weak cross-sectional dependence in large panels. *Econometric Reviews* 34, 1089–1117.
- Perron, P. (1989). The great crash, the oil price shock, and the unit root hypothesis. *Econometrica* 57(6), 1361–1401.
- Stock, J. H. & Watson, M. W. (1993). A simple estimator of cointegrating vectors in higher order integrated systems. *Econometrica* 61(4), 783–820.
- Alkhalaf, M. (2026). `xtbreakcoint`: Panel Cointegration Tests with Structural Breaks, R package 1.0.4. https://CRAN.R-project.org/package=xtbreakcoint

---

## Appendix A. `xtbreakcoint` 1.0.4 compared with BCS (2015)

Read from the CRAN source (`R/factor_estimation.R`, `R/adf_tests.R`, `R/mq_test.R`, `R/utils.R`).
It is useful for one thing: its comments show what the authors' GAUSS code does in two places (X7
and B8). It cannot serve as a reference implementation.

| ID | xtbreakcoint | BCS (2015) |
|---|---|---|
| X1 | Five models: constant, trend, levelshift, trendshift, regimeshift. "regimeshift" adds a `t − T_b` column **in first differences**, a quadratic trend in levels. No breaks in the cointegrating vector | Six models (§2.1); Models 4–6 shift the cointegrating vector |
| X2 | One **pooled** first-difference regression: a common δ for all units | Unit-by-unit projection `M_i` |
| X3 | One break per unit, found by regressing the defactored differenced residuals on `[1, step]`, which is a trend break in levels, whatever the model. "regimeshift" replaces the unit dates by their rounded mean | Break dates minimise the SSR of the model's own first-difference regression (§4.2.1); common dates by the inf statistic (§4.2.2) |
| X4 | The iteration re-runs the pooled regression **without subtracting the factors**, re-selects r every round, tolerance 1e-3 | Subtract `ΔF̃π̃` before re-estimating the breaks (§4.2.1) |
| X5 | Number of factors by IC_p2 | Panel BIC (B4) |
| X6 | ADF lags by BIC, `k_max = 4` | t-sig, `k_max` rule (B6) |
| X7 | Moments from the GAUSS code at T = 100: −0.4163/0.9834 (c), −1.5377/0.3501 (τ); one-break γ = WP Table 9's T = 100 column, nearest λ bin | Asymptotic moments (Remark 1, Table II) |
| X8 | "MQ" = per-factor ADF or PP statistics compared with −2.86 | Bai–Ng `MQ_c` and `MQ_f` (§2.6) |
| X9 | Individual rejections at −1.95 for every model | Model-specific critical values (B24) |
| X10 | Missing: known breaks, the homogeneous inf statistic, Models 3 and 6 with known common breaks, DOLS, the factor ADF for r = 1, the CD diagnostic, two breaks | §§2.6–2.10 |
| X11 | Licence GPL-3 | Compare outputs only; no code is copied into this MIT package |

## Appendix B. Tables to transcribe (P0)

| Literal | Source and location | Used for |
|---|---|---|
| Θ, Ψ for c and τ | JAE Theorem 2, Remark 1 (p. 11, PDF p. 12) | Z_c, Z_τ |
| γ moments, one break (9 λ) and two breaks (36 pairs) | JAE Table II (p. 12, PDF p. 13) | Z_γ on the grid |
| MQ(q, λ) critical values, q = 1..6, λ = 0.1..0.9, 1/5/10% | JAE Table I (p. 10, PDF p. 11) | MQ in Models 3 and 6, one break |
| Z* critical values, T = 50, 100, 250, 1/2.5/5/10% | JAE Table III (p. 13, PDF p. 14) | printed with Z* |
| Finite-sample moments, T = 50…1000 | WP 591 Table 9 (p. 40, PDF p. 41) | `moments = "finite_sample"` |
| MQ critical values without breaks, `MQ^c_{c,f}` (demeaned) and `MQ^τ_{c,f}` (detrended), m = 1..6, 1/5/10% | Bai & Ng (2004) Table I (p. 1136, PDF p. 10) | MQ in Models 1, 2, 4, 5 |
| DF quantiles (no constant, μ, τ); τ and γ idiosyncratic quantiles | simulated (`data-raw`), checked against the asymptotic DF values of MacKinnon (1996) | unit and factor p-values |
