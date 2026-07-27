# Estimate the functional autocorrelation function (FACF)

Computes the functional autocorrelation \\\widehat\rho\_\ell\\ for lags
\\\ell = 1, \dots, \\ `lag.max`, using the adaptive lag-\\\ell\\
autocovariance estimator
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md):
\$\$\widehat\rho\_\ell = \lVert \widehat\Gamma\_{N,\ell} \rVert \left\[
\int \widehat\Gamma\_{N,0}(t, t)\\ dt \right\]^{-1},\$\$ where
\\\lVert\cdot\rVert\\ is the \\\mathbb{L}^2\\ norm on \\\[0,1\]^2\\.
This is the functional analogue of the autocorrelation function of a
scalar time series and is useful for detecting serial dependence /
non-stationarity in a functional time series.

## Usage

``` r
estimate_facf(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  lag.max = 5L,
  t = NULL,
  n_grid = 25L,
  bw_grid = NULL,
  common_bw = FALSE,
  center_curves = TRUE,
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

- lag.max:

  `integer(1)`. Largest lag. Default `5`. Clamped to `N - 1` (with a
  warning) when it reaches the number of curves.

- t:

  `numeric` or `NULL`. Evaluation grid in \\\[0,1\]\\. When `NULL`
  (default) it is chosen design-aware: the shared observation grid under
  the common design (subsampled to at most `n_grid` points), or a
  regular grid of `n_grid` points under the independent design.

- n_grid:

  `integer(1)`. Size of the default evaluation grid when `t` is `NULL`.
  Default `25`.

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
  smoothing.

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

An object of class `fts_acf` (a classed
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html);
see
[adaptiveFTS_est](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS_est.md))
with columns:

- `lag`: the lag \\\ell\\.

- `norm`: the \\\mathbb{L}^2\\ norm
  \\\lVert\widehat\Gamma\_{N,\ell}\rVert\\.

- `facf`: the functional autocorrelation \\\widehat\rho\_\ell\\.

## Details

The lag-0 variance function \\\widehat\Gamma\_{N,0}(t,t)\\ (the
denominator) and each lag-\\\ell\\ autocovariance surface (the
numerator) are estimated adaptively with
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
which handles both the common and the independent observation design.
The integrals use the trapezoidal rule on the evaluation grid `t`.
Computing a full \\\|t\| \times \|t\|\\ surface for each lag can be
costly; keep `n_grid` modest, pass a coarser `t`, or supply `bw_grid`.

## References

Horváth, L., Rice, G. and Whipple, S. (2016). Adaptive bandwidth
selection in the long run covariance estimator of functional time
series. *Computational Statistics and Data Analysis*, 100, 676–693.
[doi:10.1016/j.csda.2014.06.008](https://doi.org/10.1016/j.csda.2014.06.008)

## See also

[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
[`autoplot.fts_acf()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md),
[`summary.fts_acf()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md).

## Examples

``` r
if (FALSE) { # \dontrun{
data("data_far")
facf <- estimate_facf(data = data_far, lag.max = 5, n_grid = 20)
summary(facf)
if (requireNamespace("ggplot2", quietly = TRUE)) ggplot2::autoplot(facf)
} # }
```
