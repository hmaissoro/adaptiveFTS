# Changelog

## adaptiveFTS 0.3.0 (in development)

### Breaking changes

- Arguments were renamed for clarity. There are no deprecation shims:
  update calls that name these arguments.

  | Function | Old | New |
  |----|----|----|
  | [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md), [`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md), [`estimate_facf()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_facf.md), [`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md) | `use_same_bw` | `common_bw` |
  | [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md), [`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md), [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md), [`estimate_cov_segment_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment_risk.md), [`estimate_facf()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_facf.md), [`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md) | `center` | `center_curves` |
  | [`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md), [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md) | `optbw` | `bw` |
  | [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md) | `optbw_s`, `optbw_t` | `bw_s`, `bw_t` |
  | [`estimate_mean_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md), [`estimate_autocov_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md) | `h` | `bw` |
  | [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md), [`estimate_empirical_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_autocov.md), [`estimate_empirical_mom()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_mom.md), [`estimate_empirical_XsXt_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_XsXt_autocov.md) | `h` | `presmooth_bw` |
  | [`estimate_empirical_XsXt_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_XsXt_autocov.md) | `lag` | `autocov_lag` |
  | [`estimate_mean_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md), [`estimate_mean_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_bw_rp.md), [`estimate_autocov_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md), [`estimate_autocov_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md) | `smooth_ker` (a function) | `kernel_name` (a string) |
  | [`estimate_mean_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_bw_rp.md), [`estimate_autocov_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md) | `Kfold` | `n_folds` |
  | [`estimate_autocov_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md), [`estimate_autocov_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md) | `optbw_mean`, `dt_mean_rp` | `bw_mean`, `mean_rp` |
  | [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md), [`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md) | `id_lag` | `id_conditioning_curve` |
  | [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md), [`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md), [`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md) | `n_cv_tikhonov`, `n_subgrid_bw` | `n_cv_curves`, `bw_subgrid_size` |
  | [`simulate_far()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md), [`simulate_fma()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fma.md) | `Mdistribution`, `tdistribution`, `tdesign`, `tcommon`, `int_grid`, `burnin` | `M_distribution`, `t_distribution`, `design`, `t_common`, `n_int_grid`, `n_burnin` |
  | [`simulate_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md) | `shift_var` | `intercept_var` |
  | [`simulate_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md), [`simulate_fBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fBm.md), [`simulate_far()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md), [`simulate_fma()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fma.md) | `L` | `L2` |

  The adaptive estimators’ output column names are unchanged. The
  Rubìn-Panaretos estimators, whose bandwidth argument was renamed `h`
  -\> `bw`, rename their bandwidth output columns to match:
  [`estimate_mean_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md),
  [`estimate_mean_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_bw_rp.md)
  and
  [`estimate_autocov_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md)
  now return a `bw` column instead of `h`, and
  [`estimate_autocov_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md)
  returns `bw_mean` instead of `optbw_mean`.

- The simulators’ `L` argument is renamed `L2`, because that is what it
  always was: it multiplies the covariance, so the returned path is
  `sqrt(L) * xi` and its increments satisfy
  `E[(X(t+d) - X(t))^2] = L * d^(2 H_t)`. The package’s model writes
  that coefficient `L_t^2`, so `L` was the *squared* Hölder constant and
  the Hölder constant itself was `sqrt(L)`. The estimator output columns
  were renamed `Lt`/`Ls` -\> `Lt2`/`Ls2` in 0.2.0 for exactly this
  reason; the generators are now consistent with them, and
  `simulate_far(L2 = 4)` is the setting whose
  [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md)
  estimate of `Lt2` is 4. **Only the name changes** — the numerical
  meaning of the argument is unchanged, so existing calls keep their
  behaviour by renaming `L` to `L2`.

- The Rubìn-Panaretos estimators take the kernel by name
  (`kernel_name = "epanechnikov"`) rather than as a function object,
  matching the adaptive estimators. The kernel functions themselves
  remain exported.

- `idcol` now defaults to `"id_curve"` in
  [`estimate_sigma()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md),
  [`estimate_empirical_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_autocov.md),
  [`estimate_empirical_mom()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_mom.md)
  and
  [`estimate_empirical_XsXt_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_XsXt_autocov.md),
  as it already did elsewhere. These four previously defaulted to
  `NULL`, which made
  [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  reject a `data.frame` input.

- `get_real_data_far_kenel()` is renamed
  [`get_real_data_far_kernel()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_real_data_far_kernel.md).

- [`.Spq_fun()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-Spq_fun.md)
  and
  [`.Qpq_fun()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-Qpq_fun.md)
  are no longer exported. They are internals of
  [`estimate_autocov_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md),
  which is the entry point to use.

- The `adaptive_meta` attribute of `autocov_est` and `autocov_risk`
  objects carries `common_bw` instead of `use_same_bw`.

- The internal C++ routines were renamed to match the R argument names
  (`optbw`/`optbw_s`/`optbw_t` -\> `bw`/`bw_s`/`bw_t`, `use_same_bw` -\>
  `common_bw`, `id_lag`/`n_subgrid_bw` -\>
  `id_conditioning_curve`/`bw_subgrid_size`). These functions are not
  exported, so this affects only code that reached into the compiled
  layer directly.

- [`simulate_far()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md)
  and
  [`simulate_fma()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fma.md)
  gain an `intercept_var` argument, placed after `L`. Callers that pass
  `far_kernel`/`fma_kernel` and the arguments after it *by position*
  must be updated; named calls are unaffected. Likewise
  [`simulate_fBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fBm.md)
  gains `intercept_var` after `L`, ahead of `tied`.

- [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  now validates its result instead of passing questionable data on to
  the estimators. It fails when the observation points fall outside
  `[0, 1]` (the domain the estimators assume), when the observation
  points or the observed values are not numeric, and when any value is
  missing; it warns when a curve carries repeated observation points.
  Data that used to flow through and yield `NaN` estimates now stops at
  the formatting step.

### New features

- The local regularity step is now tunable from every adaptive
  estimator. The estimators reach it through C++, which used to
  hard-code its settings, so
  [`estimate_mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md),
  [`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
  [`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md),
  [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
  [`estimate_cov_segment_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment_risk.md),
  [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md),
  [`estimate_facf()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_facf.md),
  [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md),
  [`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md)
  and
  [`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md)
  gain the four arguments
  [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md)
  already exposed: `presmooth_bw`, `Delta`, and the two new
  cross-validation controls `presmooth_bw_grid` and `presmooth_nsubset`.
  All default to `NULL` and reproduce the previous output bit-for-bit.

  `presmooth_bw` is worth knowing about: the bandwidth chosen by the
  regularity step is reused for the empirical moment and autocovariance
  estimators that feed the risk, so setting it governs the whole
  pre-smoothing stage rather than the regularity alone.

  Note that `center_curves` still does not reach the regularity step,
  which always centres the curves; this is now stated in the
  documentation.

- [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md)
  gains `presmooth_bw_grid` and `presmooth_nsubset`, the candidate grid
  and the number of curves used by the cross-validation that selects
  `presmooth_bw`. Both are ignored when `presmooth_bw` is supplied.
  Lowering `presmooth_nsubset` is the cheapest way to speed up the
  selection on large samples.

- [`simulate_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md)’s
  `shift_var` becomes `intercept_var` and is now exposed by
  [`simulate_fBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fBm.md),
  [`simulate_far()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md)
  and
  [`simulate_fma()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fma.md).
  It is the variance of a per-curve random Gaussian intercept added to
  the innovation, expressed relative to the innovation scale: the
  intercept has variance `L * intercept_var`, so `sqrt(intercept_var)`
  is its standard deviation as a fraction of the innovation standard
  deviation at `u = 1`. Being constant in `t`, it cancels in the
  increments and leaves the local regularity (`H_t`, `L_t`) unchanged;
  it only keeps the curves from all leaving the origin at the same
  point, since `Var(xi(u)) = u^(2 H_u)` vanishes as `u -> 0`. The
  default `intercept_var = 0` reproduces the previous output
  bit-for-bit.

- [`simulate_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md)
  and
  [`simulate_fBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fBm.md)
  now warn and ignore `intercept_var` when `tied = TRUE`: a tied-down
  path carrying an intercept is neither tied down at the origin nor an
  intercept-shifted mfBm, because the tie-down turns the intercept into
  a random ramp.

### Bug fixes

- [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md)
  did not validate a user-supplied `Delta`. The range check in the C++
  core tested the still-uninitialised local variable rather than the
  supplied value, so its two range clauses could never fire. `Delta` is
  now checked in R, alongside the three other local regularity
  arguments.

- [`simulate_fBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fBm.md)
  was missing a factor of 2 in the exponent of its covariance: it used
  `u^hurst + v^hurst - |u - v|^hurst` where fractional Brownian motion
  requires `u^(2 hurst) + v^(2 hurst) - |u - v|^(2 hurst)`. The
  generated process therefore had Hurst exponent `hurst / 2`, so
  `simulate_fBm(hurst = 0.6)` returned paths of exponent 0.3 and
  disagreed with
  [`simulate_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md)
  given a Hurst function constant at 0.6. **This changes the output**:
  `hurst` now means what it says, and code calibrated against the old
  behaviour must halve its `hurst` argument to reproduce the previous
  paths.
  [`simulate_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md),
  [`simulate_far()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md)
  and
  [`simulate_fma()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fma.md)
  were never affected — they build their covariance through
  [`.covariance_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-covariance_mfBm.md),
  which always used the correct exponent — so the packaged `data_far`
  dataset is unchanged.

- The C++ layer passed unprotected `Rcpp::wrap()` temporaries into
  `Rcpp::Nullable<arma::vec>` parameters at 19 call sites.
  `Rcpp::Nullable` stores a bare `SEXP` without protecting it, so the
  garbage collector could reclaim a wrapped bandwidth vector while the
  callee was still running — the callee only converts it after
  allocating R memory of its own. The result was either a hard error
  (`Not compatible with requested type: [target=double]`, with the
  reported type varying run to run) or, when the reclaimed node was
  reused as a numeric vector of the same length, silently wrong
  bandwidths. This affected
  [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)/[`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md),
  [`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md),
  [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md)
  and
  [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md).
  Every wrapped vector is now held in a protecting `Rcpp::NumericVector`
  for the duration of the call. Results are unchanged (bit-identical to
  the committed references).

- [`estimate_autocov_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md)
  always returned a cross-validation error of zero, so the selected
  bandwidth was simply the first of the grid. The held-out mean
  estimates were read from the wrong grid object and silently resolved
  to `NULL`, which collapsed the error sum to zero for every candidate.

- [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  mis-assigned observations when the rows of a curve were not contiguous
  in the input: the curve index was rebuilt from run lengths counted by
  value but written back in row order, scattering a curve’s observation
  points across its neighbours. Curves are now identified by value, and
  the result is always sorted by `id_curve` then `tobs`.

### Dependencies

- **Requires R \>= 4.1** (was 3.5.0), for the native `|>` pipe used in
  the vignette and the `inst/` demo scripts.
- `Rdpack` is no longer an `Imports`, and `RdMacros` is dropped. The
  four references are written directly in the `\references{}` sections,
  with their DOIs. `inst/REFERENCES.bib` remains as the bibliography
  source, and a new `inst/CITATION` provides `citation("adaptiveFTS")`.
- `Suggests` goes from twelve packages to five (`ggplot2`, `knitr`,
  `rmarkdown`, `testthat`, `tikzDevice`). `crosstalk`, `DT`, `dygraphs`,
  `ggpubr`, `magrittr`, `manipulateWidget` and `plotly` are dropped: the
  examples, the `inst/` demos, the vignette and the README now use the
  package’s own
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)/[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods, plain `ggplot2`, or base graphics.

### Documentation

- The references are updated: the estimation paper is published in the
  *Journal of Time Series Analysis* (2025, <doi:10.1111/jtsa.70006>) and
  the prediction paper is a 2026 preprint.
- Nearly every example is now runnable rather than wrapped in
  `\dontrun{}`, and the estimator examples use a subset of `data_far` so
  they stay fast.
- The vignette covers the mean, the autocovariance with one and with two
  bandwidths, the functional autocorrelation and the BLUP; those
  sections were previously empty headings.
- [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  gained runnable examples, an explicit description of the three
  accepted input layouts, and a stated output contract (columns, curve
  renumbering by order of first appearance, sorting).
- The `data`, `idcol`, `tcol` and `ycol` descriptions inherited by every
  estimator are now two lines pointing at
  [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md),
  instead of a 27-line copy of its input specification repeated on 25
  help pages.

## adaptiveFTS 0.2.0

### Breaking changes

- The Hölder-constant columns are renamed from `Lt`/`Ls` to `Lt2`/`Ls2`
  in the output of
  [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
  [`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md)/[`estimate_mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md),
  [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md)/[`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md)
  and
  [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md)/[`estimate_cov_segment_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment_risk.md),
  because the estimate is the squared constant (L_t^2 / L_s^2), not L_t
  / L_s. Update any code that referred to the `Lt`/`Ls` columns.

### New features

- Design-weighted, Tikhonov-regularised adaptive functional BLUP with an
  [`lm()`](https://rdrr.io/r/stats/lm.html)/[`predict()`](https://rdrr.io/r/stats/predict.html)-style
  interface:
  - [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)
    estimates every prediction-point-independent component and caches
    the adaptive bandwidths;
    [`predict()`](https://rdrr.io/r/stats/predict.html) (method
    [`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md))
    evaluates the predictor and, for `horizon > 1`, returns every
    intermediate multi-step-ahead prediction;
    [`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md)
    is a one-call wrapper.
  - The Tikhonov parameter is selected automatically by default
    (`tikhonov = NULL`), like `optbw` for `bw_grid`:
    [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)/[`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md)
    cross-validate over `tikhonov_grid` and return the selection as
    `tikhonov_cv`. Pass a numeric `tikhonov` to skip selection.
    `select_tikhonov_parameter(method = "cv")` exposes the selection
    directly (holdout under the common design, rolling origin under the
    independent design).
  - [`summary()`](https://rdrr.io/r/base/summary.html) methods for
    `blup_fit` and `blup` objects.
  - The numerical core runs in C++
    ([`blup_fit_cpp()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit_cpp.md),
    [`blup_predict_cpp()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_predict_cpp.md)).
- Design-density estimation for the independent-design weights:
  - [`estimate_density()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_density.md)
    — leave-one-out Parzen–Rosenblatt estimator with a least-squares
    cross-validation or fixed-bandwidth path.
  - [`get_density_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_density_optimal_bw.md)
    — selects the density bandwidth on a subset of curves (median of the
    per-curve LSCV optima).
- The adaptive estimators now return classed `data.table`s (see
  [`?adaptiveFTS_est`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS_est.md))
  so they gain [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods:
  - [`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
    [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
    [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
    [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md)
    and their `*_risk()` counterparts print a compact, design-aware
    [`summary()`](https://rdrr.io/r/base/summary.html) and draw a
    diagnostic `ggplot2` plot (mean/segment curves, a regularity panel,
    a covariance surface, and risk-vs-bandwidth curves marking the
    minimiser). The returned objects remain genuine `data.table`s, so
    existing code is unaffected. `ggplot2` stays a suggested dependency.
- Descriptive statistics:
  - [`estimate_facf()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_facf.md)
    — adaptive functional autocorrelation function (FACF),
    `rho_l = ||Gamma_l|| / integral Gamma_0(t, t) dt`, built on
    [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md)
    and handling both the common and independent designs. Returns a
    classed `fts_acf` object with
    [`summary()`](https://rdrr.io/r/base/summary.html),
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
    [`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
    (an ACF-style lag plot).
- [`save_plot_tikz()`](https://hmaissoro.github.io/adaptiveFTS/reference/save_plot_tikz.md)
  — export a `ggplot` (or any printable figure) to a standalone
  TikZ/LaTeX `.tex` file (optionally compiled to PDF), so the helper no
  longer has to be copied between projects. `tikzDevice` is a new
  suggested dependency.

### Deprecations

- [`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md)
  is deprecated in favour of
  [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)/[`predict()`](https://rdrr.io/r/stats/predict.html)
  and
  [`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md).
  It still works (with a warning) this release and will be removed in a
  future version. Note the new functions compute a different quantity
  (the design-weighted adaptive BLUP), so results are not
  interchangeable.

## adaptiveFTS 0.1.1

First CRAN-targeted release.

### Features

- Local regularity parameter estimation
  ([`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md)).
- Adaptive mean function estimation
  ([`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
  [`estimate_mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md)).
- Adaptive autocovariance / covariance function estimation
  ([`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
  [`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md),
  [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md)).
- Adaptive Best Linear Unbiased Predictor
  ([`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md)).
- Nadaraya–Watson smoothing and bandwidth selection, kernels, simulation
  of functional time series (FAR/FMA, fBm/mfBm) and Rubin–Panaretos
  estimators.
- Numerical core implemented in C++ via Rcpp/RcppArmadillo.

### Reproducibility and performance

- Local regularity estimation is now fully reproducible: the
  tie-breaking jitter that previously used Armadillo’s RNG (uncontrolled
  by [`set.seed()`](https://rdrr.io/r/base/Random.html)) is replaced by
  a deterministic offset.
- Substantial speed-ups of the C++ estimators (autocovariance risk and
  estimation, mean and covariance-segment risk) with bit-identical
  results, validated against committed regression references.

### Infrastructure

- Removed the `caret`, `fastmatrix` and `parallel` dependencies in
  favour of base R and `data.table`.
- Added a `testthat` (edition 3) test suite covering every exported
  function plus numerical regression tests.
- Added continuous integration (R-CMD-check on Linux/macOS/Windows, test
  coverage, lint, and a pkgdown site).
