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
  use_same_bw = FALSE,
  center = TRUE,
  kernel_name = "epanechnikov"
)
```

## Arguments

- data:

  A `data.table` (or `data.frame`), a `list` of `data.table` (or
  `data.frame`), or a `list` of `list`.

  - If `data.table`: It should contain the raw curve observations in at
    least three columns.

    - `idcol` : The name of the column containing the curve index in the
      sample. Each curve index is repeated according to the number of
      observation points.

    - `tcol` : The name of the column with observation points associated
      with each curve index.

    - `ycol` : The name of the column with observed values at each
      observation point for each curve index.

  - If `list` of `data.table`: In this case, each element in the `list`
    represents the observation data of a curve in the form of a
    `data.table` or `data.frame`. Each `data.table` contains at least
    two columns.

    - `tcol` : The name of the column with observation points for the
      curve.

    - `ycol` : The name of the column with observed values for the
      curve.

  - If `list` of `list`: In this case, `data` is a list where each
    element is the observation data of a curve, given as a `list` of two
    vectors.

    - `tcol` : The vector containing observation points for the curve.

    - `ycol` : The vector containing observed values for the curve.

- idcol:

  `character`. If `data` is given as a `data.table` or `data.frame`,
  this is the name of the column that holds the curve index. Each curve
  index is repeated according to the number of observation points. If
  `data` is a `list` of `data.table` (or `data.frame`) or a `list` of
  `list`, set `idcol = NULL`.

- tcol:

  `character`. The name of the column (or vector) containing the
  observation points for the curves.

- ycol:

  `character`. The name of the column with observed values for the
  curves.

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

  `vector (numeric)`. Bandwidth grid for selecting the optimal smoothing
  parameter for each pair (`s`, `t`). Defaults to `NULL`, which
  generates an exponential grid of \\N \lambda\\.

- use_same_bw:

  `logical`. Indicates whether the same bandwidth should be used for
  both `s` and `t`. Defaults to `FALSE`.

- center:

  `logical (TRUE or FALSE)`. Default `center = TRUE` and so the curves
  are centred when the autocovariance is estimated:
  \\\mathbb{E}(X_0(s) - \mu(s))(X\_{\ell}(t) - \mu(t))\\. Otherwise, the
  two parts \\\mathbb{E}X_0(s)X\_{\ell}(t)\\ and \\\mu(s)\mu(t)\\ will
  be estimated separately. The first part with a bandwidth obtained with
  [estimate_autocov_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md)
  and the second part with a bandwidth obtained with
  [estimate_mean_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md).

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

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

Horváth L, Rice G, Whipple S (2016). “Adaptive bandwidth selection in
the long run covariance estimator of functional time series.”
*Computational Statistics & Data Analysis*, **100**, 676–693. ISSN
0167-9473.
[doi:10.1016/j.csda.2014.06.008](https://doi.org/10.1016/j.csda.2014.06.008)
.

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
