# Estimate Empirical \\X_0(s)X\_{\ell}(t)\\ Autocovariance Function for \\\ell\\ = 0, 1, ...

This function estimates the empirical \\X_0(s)X\_{\ell}(t)\\
autocovariance function for \\\ell\\ = 0, 1, ..., used in the empirical
study of the papers Maissoro, Patilea and Vimond (2025) and Maissoro,
Patilea and Vimond (2026).

## Usage

``` r
estimate_empirical_XsXt_autocov(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  cross_lag = 1,
  autocov_lag = c(0, 1, 2),
  presmooth_bw = NULL,
  center = FALSE,
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

  `vector (numeric)`. First argument in \\X_0(s)X\_{\ell}(t)\\,
  corresponding to observation points `s` in the pair (`s`, `t`). Must
  be of the same length as `t`.

- t:

  `vector (numeric)`. Second argument in \\X_0(s)X\_{\ell}(t)\\,
  corresponding to observation points `t` in the pair (`s`, `t`). Must
  be of the same length as `s`.

- cross_lag:

  `integer (positive integer)`. The lag \\\ell\\ in
  \\X_0(s)X\_{\ell}(t)\\.

- autocov_lag:

  `vector (integer)`. Lags at which the autocovariance of the scalar
  series \\n \mapsto X_n(s)X\_{n+\ell}(t)\\ is estimated, \\\ell\\ being
  `cross_lag`. If `NULL`, only \\\mathbb{E}X_0(s)X\_{\ell}(t)\\ is
  returned.

- presmooth_bw:

  `numeric (positive vector or scalar)`. Bandwidth used to presmooth
  each curve before the estimation. A scalar applies the same bandwidth
  to every curve; a vector must hold one bandwidth per curve, in the
  order the curves appear in `data`. Default `NULL` selects it by
  cross-validation, see
  [get_nw_optimal_bw](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

- center:

  `logical`. If `TRUE`, the estimated autocovariance is centered:
  \\\mathbb{E}(X_0(s) - \mu(s))(X\_{\ell}(t) - \mu(t))\\. Defaults to
  `FALSE`, providing \\\mathbb{E}X_0(s)X\_{\ell}(t)\\.

- kernel_name:

  `string`. Kernel function for estimation; defaults to "epanechnikov".
  Supported kernels are: "epanechnikov", "biweight", "triweight",
  "tricube", "triangular", and "uniform".

## Value

A `data.table` with columns:

- s : First argument in \\X_0(s)X\_{\ell}(t)\\.

- t : Second argument in \\X_0(s)X\_{\ell}(t)\\.

- cross_lag : Lag \\\ell\\ in \\X_0(s)X\_{\ell}(t)\\.

- lag : The lags at which the autocovariance of \\X_0(s)X\_{\ell}(t)\\
  is estimated; `NA` if `autocov_lag = NULL`.

- EXsXt_cross_lag : Mean of \\X_0(s)X\_{\ell}(t)\\.

- XsXt_autocov : Autocovariance estimates of \\X_0(s)X\_{\ell}(t)\\ for
  each `autocov_lag`; `NA` if `autocov_lag = NULL`.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
Functional Time Series. *arXiv preprint* arXiv:2609.xxxxx.

## See also

[`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

## Examples

``` r
# Load data
data("data_far")

# Example 1: Estimate autocovariance without centering
dt_empirical_cov <- estimate_empirical_XsXt_autocov(
  data = data_far,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  cross_lag = 1,
  autocov_lag = c(0, 1, 2),
  presmooth_bw = 0.1,
  center = FALSE,
  kernel_name = "epanechnikov"
)
dt_empirical_cov
#>        s     t cross_lag   lag EXsXt_cross_lag XsXt_autocov
#>    <num> <num>     <num> <num>           <num>        <num>
#> 1:   0.2  0.25         1     0        58539.96    1574630.1
#> 2:   0.2  0.25         1     1        58539.96    1226200.2
#> 3:   0.2  0.25         1     2        58539.96     747558.4
#> 4:   0.4  0.50         1     0        57424.45    1635518.6
#> 5:   0.4  0.50         1     1        57424.45    1238193.8
#> 6:   0.4  0.50         1     2        57424.45     701510.2
#> 7:   0.8  0.75         1     0        57456.35    1807021.6
#> 8:   0.8  0.75         1     1        57456.35    1326555.2
#> 9:   0.8  0.75         1     2        57456.35     631896.2

# Example 2: Estimate autocovariance with centering
dt_empirical_cov_centered <- estimate_empirical_XsXt_autocov(
  data = data_far,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  cross_lag = 1,
  autocov_lag = c(0, 1, 2),
  presmooth_bw = 0.1,
  center = TRUE,
  kernel_name = "epanechnikov"
)
dt_empirical_cov_centered
#>        s     t cross_lag   lag EXsXt_cross_lag XsXt_autocov
#>    <num> <num>     <num> <num>           <num>        <num>
#> 1:   0.2  0.25         1     0        5.029341    85.900437
#> 2:   0.2  0.25         1     1        5.029341    36.263913
#> 3:   0.2  0.25         1     2        5.029341     5.043176
#> 4:   0.4  0.50         1     0        4.981649    87.306598
#> 5:   0.4  0.50         1     1        4.981649    27.211200
#> 6:   0.4  0.50         1     2        4.981649    -1.967343
#> 7:   0.8  0.75         1     0        5.174762   106.003910
#> 8:   0.8  0.75         1     1        5.174762    33.311083
#> 9:   0.8  0.75         1     2        5.174762     2.787892

```
