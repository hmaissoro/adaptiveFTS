# Estimate Empirical \\X_0(s)X\_{\ell}(t)\\ Autocovariance Function for \\\ell\\ = 0, 1, ...

This function estimates the empirical \\X_0(s)X\_{\ell}(t)\\
autocovariance function for \\\ell\\ = 0, 1, ..., used in the empirical
study of the papers Maissoro et al. (2024) and Maissoro et al. (2025) .

## Usage

``` r
estimate_empirical_XsXt_autocov(
  data,
  idcol = NULL,
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  cross_lag = 1,
  lag = c(0, 1, 2),
  h = NULL,
  center = FALSE,
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

- lag:

  `vector (integer)`. Lag for the autocovariance of
  \\X_0(s)X\_{\ell}(t)\\. If `lag = NULL`, only
  \\\mathbb{E}X_0(s)X\_{\ell}(t)\\ is returned.

- h:

  `numeric (positive vector or scalar)`. Smoothing bandwidth parameter.
  Defaults to `NULL`, in which case `h` is estimated via
  Cross-Validation on a subset of curves. If `h` is a scalar, all curves
  are smoothed with the same bandwidth; if a vector, it should match the
  number of curves in `data`.

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

- lag : Lags for autocovariance estimation of \\X_0(s)X\_{\ell}(t)\\;
  contains `NA` if `lag = NULL`.

- EXsXt_cross_lag : Mean of \\X_0(s)X\_{\ell}(t)\\.

- XsXt_autocov : Autocovariance estimates of \\X_0(s)X\_{\ell}(t)\\ for
  each `lag`; contains `NA` if `lag = NULL`.

## References

Maissoro H, Patilea V, Vimond M (2024). “Adaptive estimation for Weakly
Dependent Functional Times Series.” *arXiv preprint arXiv:2403.13706*.  
  
Maissoro H, Patilea V, Vimond M (2025). “Adaptive prediction for
Functional Times Series.” *arXiv preprint arXiv:2501.xxxxx*.

## See also

[`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

## Examples

``` r
if (FALSE) { # \dontrun{
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
  lag = c(0, 1, 2),
  h = 0.1,
  center = FALSE,
  kernel_name = "epanechnikov"
)
dt_empirical_cov

# Example 2: Estimate autocovariance with centering
dt_empirical_cov_centered <- estimate_empirical_XsXt_autocov(
  data = data_far,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  cross_lag = 1,
  lag = c(0, 1, 2),
  h = 0.1,
  center = TRUE,
  kernel_name = "epanechnikov"
)
dt_empirical_cov_centered
} # }
```
