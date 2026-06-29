# Estimate Empirical Autocovariance Function

This function estimates the empirical autocovariance function used in
the empirical study section of the papers Maissoro et al. (2024) and
Maissoro et al. (2025) .

## Usage

``` r
estimate_empirical_autocov(
  data,
  idcol = NULL,
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  lag = c(0, 1, 2),
  h = NULL,
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

- t:

  `vector (numeric)`. Observation points at which we want to estimate
  the empirical autocovariance function.

- lag:

  `vector (integer)`. Lag of the autocovariance.

- h:

  `numeric (positive vector or scalar)`. The smoothing bandwidth
  parameter. Default `h = NULL` and thus it will be estimated by
  Cross-Validation on a subset of curves. If `h` is a `scalar`, then all
  curves will be smoothed with the same bandwidth. Otherwise, if `h` is
  a `vector`, its length must be equal to the number of curves in `data`
  and each element of the vector must correspond to a curve given in the
  same order as in `data`.

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

## Value

A `data.table` with three columns: `t`, `lag`, and `autocov`
corresponding to the estimated autocovariance.

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

# Estimate empirical autocovariance with a specified bandwidth
dt_empirical_autocov <- estimate_empirical_autocov(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), lag = c(1, 2), h = 0.1,
  kernel_name = "epanechnikov")
dt_empirical_autocov

# Estimate empirical autocovariance with Cross-Validation bandwidth selection
dt_empirical_autocov_cv <- estimate_empirical_autocov(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), lag = c(1, 2), h = NULL,
  kernel_name = "epanechnikov")
dt_empirical_autocov_cv
} # }
```
