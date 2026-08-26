# Estimate the Risk of the Autocovariance Function Estimator

Estimates the risk of the adaptive lag-\\\ell\\ autocovariance function
estimator over a grid of candidate bandwidths, for \\\ell = 0, 1,
\ldots\\ (\\\ell = 0\\ being the covariance function). Minimising it
over the grid at each pair (`s`, `t`) is what
[estimate_autocov](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md)
does to select its bandwidths.

## Usage

``` r
estimate_autocov_risk(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  lag = 1,
  bw_grid = NULL,
  common_bw = FALSE,
  center_curves = TRUE,
  presmooth_bw = NULL,
  Delta = NULL,
  presmooth_bw_grid = NULL,
  presmooth_nsubset = NULL,
  kernel_name = "epanechnikov"
)
```

## Arguments

- data:

  Raw curve observations, as a `data.table` (or `data.frame`) in long
  format, or as a `list` with one element per curve. See
  [`format_data`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  for the accepted layouts and for the `id_curve` / `tobs` / `X` columns
  they are converted to.

- idcol:

  `character(1)` or `NULL`. Name of the column holding the curve index
  when `data` is a single table. Must be `NULL` when `data` is a list of
  curves.

- tcol:

  `character(1)`. Name of the column (or vector) holding the observation
  points of the curves.

- ycol:

  `character(1)`. Name of the column (or vector) holding the values
  observed at those points.

- s:

  `vector (numeric)`. First argument of the autocovariance function: the
  points `s` of the pairs (`s`, `t`). Must have the same length as `t`.

- t:

  `vector (numeric)`. Second argument of the autocovariance function:
  the points `t` of the pairs (`s`, `t`). Must have the same length as
  `s`.

- lag:

  `integer (non-negative)`. Lag \\\ell\\ of the autocovariance;
  `lag = 0` gives the covariance function.

- bw_grid:

  `vector (numeric)`. Candidate bandwidths, from which
  [estimate_autocov](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md)
  picks the risk-minimising one for each pair. Default `NULL` builds the
  grid from the data; see Details.

- common_bw:

  `logical`. If `TRUE`, a single bandwidth is selected for both
  arguments of the autocovariance; if `FALSE` (default), one bandwidth
  per argument. See Details.

- center_curves:

  `logical`. If `TRUE` (default), the curves are centred before
  smoothing. It governs the moment and autocovariance estimators only:
  the local regularity step always centres the curves.

- presmooth_bw:

  `numeric (positive vector or scalar)`. Bandwidth of the
  Nadaraya-Watson estimator used to presmooth each curve before the
  regularity is estimated. A scalar applies the same bandwidth to every
  curve; a vector must hold one bandwidth per curve, in the order the
  curves appear in `data`. Default `NULL` selects a single bandwidth by
  cross-validation over every curve, as
  [get_nw_optimal_bw](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md)
  does.

- Delta:

  `numeric (positive)`. Length of the neighbourhood around each point of
  `t` used to estimate the local regularity. Default `NULL` sets it from
  the data; see Details.

- presmooth_bw_grid:

  `vector (numeric)`. Candidate bandwidths of the cross-validation that
  selects `presmooth_bw` when the latter is `NULL`. Default `NULL` uses
  the default grid of
  [get_nw_optimal_bw](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).
  Ignored when `presmooth_bw` is supplied.

- presmooth_nsubset:

  `integer (positive)`. Number of curves used by that cross-validation.
  Default `NULL` uses `min(70, floor(N / 2))` curves, where \\N\\ is the
  number of curves. Lower it to speed up the selection on large samples.
  Ignored when `presmooth_bw` is supplied.

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

A `data.table` with one row per (pair, bandwidth) combination and
columns:

- `s`, `t`: the arguments of the autocovariance function.

- `hs`, `ht`: the candidate bandwidths for `s` and for `t`. Identical
  when `common_bw = TRUE`.

- `PNl`: the number of curves contributing to the estimate at (`s`,
  `t`), \\P\_{N,\ell}(s,t;h_s,h_t)\\.

- `locreg_bw`: the bandwidth used to estimate the local regularity.

- `Hs`, `Ls2`: the estimated local exponent \\H_s\\ and squared Hölder
  constant \\L_s^2\\ at `s`.

- `Ht`, `Lt2`: the same quantities at `t`.

- `bias_term`, `variance_term`, `dependence_term`: the three components
  of the risk.

- `autocov_risk`: the estimated risk.

## Details

Two estimators are covered. With `common_bw = TRUE` a single bandwidth
is selected for both arguments, the one-bandwidth estimator of Maissoro,
Patilea and Vimond (2025); the returned `hs` and `ht` then hold the same
value. With `common_bw = FALSE` (default) the risk is minimised over
pairs \\(h_s, h_t)\\, the two-bandwidth estimator of Maissoro, Patilea
and Vimond (2026), which adapts to the regularity of the process at `s`
and at `t` separately. The second is the more flexible but explores the
square of the grid, so it costs noticeably more.

As for the mean, the risk splits into a bias, a variance and a
dependence term, each returned separately. The local regularity
parameters are estimated internally at `s` and at `t`.

Left to `NULL`, `bw_grid` is a 20-point geometric grid running from
\\4(N\widehat\lambda)^{-0.9}\\ to \\4(N\widehat\lambda)^{-1/3}\\, where
\\N\\ is the number of curves and \\\widehat\lambda\\ the average number
of observation points per curve.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
Functional Time Series. *arXiv preprint* arXiv:2609.xxxxx.

## See also

[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_sigma()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md).

## Examples

``` r
data("data_far")

dt_autocov_risk <- estimate_autocov_risk(
  data = data_far[data_far$id_curve <= 20, ],
  idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5), t = c(1/4, 1/2), lag = 1,
  bw_grid = seq(0.04, 0.15, length.out = 5), common_bw = TRUE,
  center_curves = TRUE, kernel_name = "epanechnikov")

# The risk-minimising bandwidth for each pair.
dt_autocov_risk[, list(hs = hs[which.min(autocov_risk)]), by = c("s", "t")]
#>        s     t    hs
#>    <num> <num> <num>
#> 1:   0.2  0.25  0.04
#> 2:   0.4  0.50  0.04
```
