# Bandwidth estimation using cross-validation for the Rubín and Panaretos (2020) mean function estimator.

This function estimates the optimal bandwidth for the mean function
estimator using cross-validation, as described in Rubín and Panaretos
(2020) .

## Usage

``` r
estimate_mean_bw_rp(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  Kfold = 10,
  bw_grid = seq(0.001, 0.15, len = 45),
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

- Kfold:

  `integer (positive)`. Number of fold for the cross-validation.

- bw_grid:

  `vector (numeric)`. The bandwidth grid.

- smooth_ker:

  `function`. The kernel function of the Nadaraya-Watson estimator.
  Default `smooth_ker = epanechnikov`.

## Value

A `data.table` containing the following columns.

- h : The candidate bandwidth.

- cv_error : The estimates of the Cross-Validation error for each `h`.

## References

Rubín T, Panaretos VM (2020). “Sparsely observed functional time series:
estimation and prediction.” *Electronic Journal of Statistics*,
**14**(1), 1137 – 1210.
[doi:10.1214/20-EJS1690](https://doi.org/10.1214/20-EJS1690) .

## See also

[`estimate_mean_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate a FAR A process
dt_far <- simulate_far(N = 50, lambda = 70,
                       tdesign = "random",
                       Mdistribution = rpois,
                       tdistribution = runif,
                       tcommon = NULL,
                       hurst_fun = hurst_logistic,
                       L = 4,
                       far_kernel = get_real_data_far_kenel,
                       far_mean = get_real_data_mean,
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)

# Add noise
dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]

## Estimate the bandwidth by Cross-Validation
dt_bw_mean_rp <- estimate_mean_bw_rp(
  data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  Kfold = 10, bw_grid = seq(0.001, 0.15, len = 45),
  smooth_ker = epanechnikov)

## Plot the Cross-Validation error
dygraphs::dygraph(dt_bw_mean_rp)

## Select the best bandwidth
optbw <- dt_bw_mean_rp[, h[which.min(cv_error)]]

## Estimate the mean function
dt_mean_rp <- estimate_mean_rp(
  data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), h = optbw, smooth_ker = epanechnikov)

DT::datatable(data = dt_mean_rp[, lapply(.SD, function(X) round(X, 5))])

} # }
```
