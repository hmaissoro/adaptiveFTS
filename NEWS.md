# adaptiveFTS 0.3.0 (in development)

## Breaking changes

* Arguments were renamed for clarity. There are no deprecation shims: update
  calls that name these arguments.

  | Function | Old | New |
  | --- | --- | --- |
  | `estimate_autocov()`, `estimate_autocov_risk()`, `estimate_facf()`, `predict_curve()` | `use_same_bw` | `common_bw` |
  | `estimate_autocov()`, `estimate_autocov_risk()`, `estimate_cov_segment()`, `estimate_cov_segment_risk()`, `estimate_facf()`, `predict_curve()` | `center` | `center_curves` |
  | `estimate_mean()`, `estimate_cov_segment()` | `optbw` | `bw` |
  | `estimate_autocov()` | `optbw_s`, `optbw_t` | `bw_s`, `bw_t` |
  | `estimate_mean_rp()`, `estimate_autocov_rp()` | `h` | `bw` |
  | `estimate_locreg()`, `estimate_empirical_autocov()`, `estimate_empirical_mom()`, `estimate_empirical_XsXt_autocov()` | `h` | `presmooth_bw` |
  | `estimate_empirical_XsXt_autocov()` | `lag` | `autocov_lag` |
  | `estimate_mean_rp()`, `estimate_mean_bw_rp()`, `estimate_autocov_rp()`, `estimate_autocov_bw_rp()` | `smooth_ker` (a function) | `kernel_name` (a string) |
  | `estimate_mean_bw_rp()`, `estimate_autocov_bw_rp()` | `Kfold` | `n_folds` |
  | `estimate_autocov_rp()`, `estimate_autocov_bw_rp()` | `optbw_mean`, `dt_mean_rp` | `bw_mean`, `mean_rp` |
  | `blup_fit()`, `blup()` | `id_lag` | `id_conditioning_curve` |
  | `blup_fit()`, `blup()`, `select_tikhonov_parameter()` | `n_cv_tikhonov`, `n_subgrid_bw` | `n_cv_curves`, `bw_subgrid_size` |
  | `simulate_far()`, `simulate_fma()` | `Mdistribution`, `tdistribution`, `tdesign`, `tcommon`, `int_grid`, `burnin` | `M_distribution`, `t_distribution`, `design`, `t_common`, `n_int_grid`, `n_burnin` |
  | `simulate_mfBm()` | `shift_var` | `intercept_var` |
  | `simulate_mfBm()`, `simulate_fBm()`, `simulate_far()`, `simulate_fma()` | `L` | `L2` |

  The adaptive estimators' output column names are unchanged. The
  Rubìn-Panaretos estimators, whose bandwidth argument was renamed `h` -> `bw`,
  rename their bandwidth output columns to match: `estimate_mean_rp()`,
  `estimate_mean_bw_rp()` and `estimate_autocov_bw_rp()` now return a `bw`
  column instead of `h`, and `estimate_autocov_rp()` returns `bw_mean` instead
  of `optbw_mean`.

* The simulators' `L` argument is renamed `L2`, because that is what it always
  was: it multiplies the covariance, so the returned path is `sqrt(L) * xi` and
  its increments satisfy `E[(X(t+d) - X(t))^2] = L * d^(2 H_t)`. The package's
  model writes that coefficient `L_t^2`, so `L` was the *squared* Hölder
  constant and the Hölder constant itself was `sqrt(L)`. The estimator output
  columns were renamed `Lt`/`Ls` -> `Lt2`/`Ls2` in 0.2.0 for exactly this
  reason; the generators are now consistent with them, and
  `simulate_far(L2 = 4)` is the setting whose `estimate_locreg()` estimate of
  `Lt2` is 4. **Only the name changes** — the numerical meaning of the argument
  is unchanged, so existing calls keep their behaviour by renaming `L` to `L2`.
* The Rubìn-Panaretos estimators take the kernel by name (`kernel_name = "epanechnikov"`)
  rather than as a function object, matching the adaptive estimators. The kernel
  functions themselves remain exported.
* `idcol` now defaults to `"id_curve"` in `estimate_sigma()`,
  `estimate_empirical_autocov()`, `estimate_empirical_mom()` and
  `estimate_empirical_XsXt_autocov()`, as it already did elsewhere. These four
  previously defaulted to `NULL`, which made `format_data()` reject a
  `data.frame` input.
* `get_real_data_far_kenel()` is renamed `get_real_data_far_kernel()`.
* `.Spq_fun()` and `.Qpq_fun()` are no longer exported. They are internals of
  `estimate_autocov_rp()`, which is the entry point to use.
* The `adaptive_meta` attribute of `autocov_est` and `autocov_risk` objects
  carries `common_bw` instead of `use_same_bw`.
* The internal C++ routines were renamed to match the R argument names
  (`optbw`/`optbw_s`/`optbw_t` -> `bw`/`bw_s`/`bw_t`, `use_same_bw` ->
  `common_bw`, `id_lag`/`n_subgrid_bw` -> `id_conditioning_curve`/`bw_subgrid_size`).
  These functions are not exported, so this affects only code that reached into
  the compiled layer directly.
* `simulate_far()` and `simulate_fma()` gain an `intercept_var` argument, placed
  after `L`. Callers that pass `far_kernel`/`fma_kernel` and the arguments after
  it *by position* must be updated; named calls are unaffected. Likewise
  `simulate_fBm()` gains `intercept_var` after `L`, ahead of `tied`.
* `format_data()` now validates its result instead of passing questionable data
  on to the estimators. It fails when the observation points fall outside
  `[0, 1]` (the domain the estimators assume), when the observation points or
  the observed values are not numeric, and when any value is missing; it warns
  when a curve carries repeated observation points. Data that used to flow
  through and yield `NaN` estimates now stops at the formatting step.

## New features

* `simulate_mfBm()`'s `shift_var` becomes `intercept_var` and is now exposed by
  `simulate_fBm()`, `simulate_far()` and `simulate_fma()`. It is the variance of a per-curve random
  Gaussian intercept added to the innovation, expressed relative to the
  innovation scale: the intercept has variance `L * intercept_var`, so
  `sqrt(intercept_var)` is its standard deviation as a fraction of the innovation
  standard deviation at `u = 1`. Being constant in `t`, it cancels in the
  increments and leaves the local regularity (`H_t`, `L_t`) unchanged; it only
  keeps the curves from all leaving the origin at the same point, since
  `Var(xi(u)) = u^(2 H_u)` vanishes as `u -> 0`. The default `intercept_var = 0`
  reproduces the previous output bit-for-bit.
* `simulate_mfBm()` and `simulate_fBm()` now warn and ignore `intercept_var` when `tied = TRUE`: a
  tied-down path carrying an intercept is neither tied down at the origin nor an
  intercept-shifted mfBm, because the tie-down turns the intercept into a random
  ramp.

## Bug fixes

* `simulate_fBm()` was missing a factor of 2 in the exponent of its covariance:
  it used `u^hurst + v^hurst - |u - v|^hurst` where fractional Brownian motion
  requires `u^(2 hurst) + v^(2 hurst) - |u - v|^(2 hurst)`. The generated process
  therefore had Hurst exponent `hurst / 2`, so `simulate_fBm(hurst = 0.6)`
  returned paths of exponent 0.3 and disagreed with `simulate_mfBm()` given a
  Hurst function constant at 0.6. **This changes the output**: `hurst` now means
  what it says, and code calibrated against the old behaviour must halve its
  `hurst` argument to reproduce the previous paths. `simulate_mfBm()`,
  `simulate_far()` and `simulate_fma()` were never affected — they build their
  covariance through `.covariance_mfBm()`, which always used the correct
  exponent — so the packaged `data_far` dataset is unchanged.
* The C++ layer passed unprotected `Rcpp::wrap()` temporaries into
  `Rcpp::Nullable<arma::vec>` parameters at 19 call sites. `Rcpp::Nullable`
  stores a bare `SEXP` without protecting it, so the garbage collector could
  reclaim a wrapped bandwidth vector while the callee was still running — the
  callee only converts it after allocating R memory of its own. The result was
  either a hard error (`Not compatible with requested type: [target=double]`,
  with the reported type varying run to run) or, when the reclaimed node was
  reused as a numeric vector of the same length, silently wrong bandwidths.
  This affected `blup_fit()`/`select_tikhonov_parameter()`, `predict_curve()`,
  `estimate_autocov()` and `estimate_cov_segment()`. Every wrapped vector is now
  held in a protecting `Rcpp::NumericVector` for the duration of the call.
  Results are unchanged (bit-identical to the committed references).
* `estimate_autocov_bw_rp()` always returned a cross-validation error of zero,
  so the selected bandwidth was simply the first of the grid. The held-out mean
  estimates were read from the wrong grid object and silently resolved to
  `NULL`, which collapsed the error sum to zero for every candidate.
* `format_data()` mis-assigned observations when the rows of a curve were not
  contiguous in the input: the curve index was rebuilt from run lengths counted
  by value but written back in row order, scattering a curve's observation
  points across its neighbours. Curves are now identified by value, and the
  result is always sorted by `id_curve` then `tobs`.

## Dependencies

* **Requires R >= 4.1** (was 3.5.0), for the native `|>` pipe used in the
  vignette and the `inst/` demo scripts.
* `Rdpack` is no longer an `Imports`, and `RdMacros` is dropped. The four
  references are written directly in the `\references{}` sections, with their
  DOIs. `inst/REFERENCES.bib` remains as the bibliography source, and a new
  `inst/CITATION` provides `citation("adaptiveFTS")`.
* `Suggests` goes from twelve packages to five (`ggplot2`, `knitr`, `rmarkdown`,
  `testthat`, `tikzDevice`). `crosstalk`, `DT`, `dygraphs`, `ggpubr`,
  `magrittr`, `manipulateWidget` and `plotly` are dropped: the examples, the
  `inst/` demos, the vignette and the README now use the package's own
  `plot()`/`autoplot()` methods, plain `ggplot2`, or base graphics.

## Documentation

* The references are updated: the estimation paper is published in the
  *Journal of Time Series Analysis* (2025, doi:10.1111/jtsa.70006) and the
  prediction paper is a 2026 preprint.
* Nearly every example is now runnable rather than wrapped in `\dontrun{}`, and
  the estimator examples use a subset of `data_far` so they stay fast.
* The vignette covers the mean, the autocovariance with one and with two
  bandwidths, the functional autocorrelation and the BLUP; those sections were
  previously empty headings.
* `format_data()` gained runnable examples, an explicit description of the three
  accepted input layouts, and a stated output contract (columns, curve
  renumbering by order of first appearance, sorting).
* The `data`, `idcol`, `tcol` and `ycol` descriptions inherited by every
  estimator are now two lines pointing at `format_data()`, instead of a
  27-line copy of its input specification repeated on 25 help pages.

# adaptiveFTS 0.2.0

## Breaking changes

* The Hölder-constant columns are renamed from `Lt`/`Ls` to `Lt2`/`Ls2` in the
  output of `estimate_locreg()`, `estimate_mean()`/`estimate_mean_risk()`,
  `estimate_autocov()`/`estimate_autocov_risk()` and
  `estimate_cov_segment()`/`estimate_cov_segment_risk()`, because the estimate is
  the squared constant (L_t^2 / L_s^2), not L_t / L_s. Update any code that
  referred to the `Lt`/`Ls` columns.

## New features

* Design-weighted, Tikhonov-regularised adaptive functional BLUP with an
  `lm()`/`predict()`-style interface:
  * `blup_fit()` estimates every prediction-point-independent component and
    caches the adaptive bandwidths; `predict()` (method `predict.blup_fit()`)
    evaluates the predictor and, for `horizon > 1`, returns every intermediate
    multi-step-ahead prediction; `blup()` is a one-call wrapper.
  * The Tikhonov parameter is selected automatically by default (`tikhonov =
    NULL`), like `optbw` for `bw_grid`: `blup_fit()`/`blup()` cross-validate over
    `tikhonov_grid` and return the selection as `tikhonov_cv`. Pass a numeric
    `tikhonov` to skip selection. `select_tikhonov_parameter(method = "cv")`
    exposes the selection directly (holdout under the common design, rolling
    origin under the independent design).
  * `summary()` methods for `blup_fit` and `blup` objects.
  * The numerical core runs in C++ (`blup_fit_cpp()`, `blup_predict_cpp()`).
* Design-density estimation for the independent-design weights:
  * `estimate_density()` — leave-one-out Parzen–Rosenblatt estimator with a
    least-squares cross-validation or fixed-bandwidth path.
  * `get_density_optimal_bw()` — selects the density bandwidth on a subset of
    curves (median of the per-curve LSCV optima).
* The adaptive estimators now return classed `data.table`s (see
  `?adaptiveFTS_est`) so they gain `summary()`, `plot()` and
  `ggplot2::autoplot()` methods:
  * `estimate_mean()`, `estimate_locreg()`, `estimate_autocov()`,
    `estimate_cov_segment()` and their `*_risk()` counterparts print a compact,
    design-aware `summary()` and draw a diagnostic `ggplot2` plot (mean/segment
    curves, a regularity panel, a covariance surface, and risk-vs-bandwidth
    curves marking the minimiser). The returned objects remain genuine
    `data.table`s, so existing code is unaffected. `ggplot2` stays a suggested
    dependency.
* Descriptive statistics:
  * `estimate_facf()` — adaptive functional autocorrelation function (FACF),
    `rho_l = ||Gamma_l|| / integral Gamma_0(t, t) dt`, built on
    `estimate_autocov()` and handling both the common and independent designs.
    Returns a classed `fts_acf` object with `summary()`, `plot()` and
    `ggplot2::autoplot()` (an ACF-style lag plot).
* `save_plot_tikz()` — export a `ggplot` (or any printable figure) to a
  standalone TikZ/LaTeX `.tex` file (optionally compiled to PDF), so the helper
  no longer has to be copied between projects. `tikzDevice` is a new suggested
  dependency.

## Deprecations

* `predict_curve()` is deprecated in favour of `blup_fit()`/`predict()` and
  `blup()`. It still works (with a warning) this release and will be removed in
  a future version. Note the new functions compute a different quantity (the
  design-weighted adaptive BLUP), so results are not interchangeable.

# adaptiveFTS 0.1.1

First CRAN-targeted release.

## Features

* Local regularity parameter estimation (`estimate_locreg()`).
* Adaptive mean function estimation (`estimate_mean()`, `estimate_mean_risk()`).
* Adaptive autocovariance / covariance function estimation
  (`estimate_autocov()`, `estimate_autocov_risk()`, `estimate_cov_segment()`).
* Adaptive Best Linear Unbiased Predictor (`predict_curve()`).
* Nadaraya–Watson smoothing and bandwidth selection, kernels, simulation of
  functional time series (FAR/FMA, fBm/mfBm) and Rubin–Panaretos estimators.
* Numerical core implemented in C++ via Rcpp/RcppArmadillo.

## Reproducibility and performance

* Local regularity estimation is now fully reproducible: the tie-breaking jitter
  that previously used Armadillo's RNG (uncontrolled by `set.seed()`) is replaced
  by a deterministic offset.
* Substantial speed-ups of the C++ estimators (autocovariance risk and
  estimation, mean and covariance-segment risk) with bit-identical results,
  validated against committed regression references.

## Infrastructure

* Removed the `caret`, `fastmatrix` and `parallel` dependencies in favour of
  base R and `data.table`.
* Added a `testthat` (edition 3) test suite covering every exported function plus
  numerical regression tests.
* Added continuous integration (R-CMD-check on Linux/macOS/Windows, test
  coverage, lint, and a pkgdown site).
