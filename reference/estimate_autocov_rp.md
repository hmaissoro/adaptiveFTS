# Estimate lag-\\\ell\\ (\\\ell \leq 0\\) autocovariance function using Rubín and Panaretos (2020) method

Estimate lag-\\\ell\\ (\\\ell \leq 0\\) autocovariance function using
Rubín and Panaretos (2020) method

## Usage

``` r
estimate_autocov_rp(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  lag = 1,
  h,
  optbw_mean = NULL,
  dt_mean_rp = NULL,
  smooth_ker = epanechnikov
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

- s:

  `vector (numeric)`. First argument of the autocovariance function. It
  corresponds to the observation points `s` in the pair (`s`, `t`). It
  has to be of the same length as the `t`

- t:

  `vector (numeric)`. Second argument of the autocovariance function. It
  corresponds to the observation points `t` in the pair (`s`, `t`). It
  has to be of the same length as the `s`.

- lag:

  `integer (positive integer)`. Lag of the autocovariance.

- h:

  `numeric (positive scalar)`. The bandwidth of the estimator.

- optbw_mean:

  `numeric (positive scalar)`. Optimal bandwidth for the mean function
  estimator. It is `NULL` if `dt_mean_rp` is not `NULL`.

- dt_mean_rp:

  `data.table`. It contains the estimates of the mean function at each
  observation point for each curve. The name of the curve identification
  column must be `id_curve`, the observation points column `tobs` and
  the mean estimates column `muhat_RP`. Default `dt_mean_rp = NULL` and
  so it will be estimated.

- smooth_ker:

  `function`. The kernel function of the Nadaraya-Watson estimator.
  Default `smooth_ker = epanechnikov`.

## Value

A `data.table` containing the following columns.

- s : The first argument of the autocovariance function.

- t : The second argument of the autocovariance function.

- lag : The lag of the autocovariance. It corresponds to \\\ell \leq
  0\\.

- optbw_mean : The optimal bandwidth for the mean function estimator.

- h : The bandwidth used to estimate the lag-\\\ell\\, \\\ell \leq 0\\
  autocovariance function

- autocovhat_rp : The estimates of the lag-\\\ell\\ autocovariance
  function for each (`s`, `t`) using Rubìn and Panaretos (2020) method.

## References

Rubín T, Panaretos VM (2020). “Sparsely observed functional time series:
estimation and prediction.” *Electronic Journal of Statistics*,
**14**(1), 1137 – 1210.
[doi:10.1214/20-EJS1690](https://doi.org/10.1214/20-EJS1690) .
