# Estimate Mean Function

This function estimates the mean function of an underlying process using
the adaptive estimator described in Maissoro et al. (2024) .

## Usage

``` r
estimate_mean(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  optbw = NULL,
  bw_grid = NULL,
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

  `vector (numeric)`. Observation points where the mean function of the
  underlying process is estimated.

- optbw:

  `vector (numeric)`. Optimal bandwidth parameters for mean function
  estimation at each `t`. If `optbw = NULL` (default), it will be
  estimated using the
  [estimate_mean_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md)
  function.

- bw_grid:

  `vector (numeric)`. A bandwidth grid from which the best smoothing
  parameter is selected for each `t`. Default is `NULL`, in which case
  it is defined as an exponential grid of \\N \lambda\\.

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

## Value

A `data.table` containing the following columns:

- `t` : The observation points at which the mean function is estimated.

- `optbw` : The optimal bandwidth used to estimate the mean function at
  each `t`.

- `Ht` : Local exponent estimates for each `t`, corresponding to
  \\H_t\\.

- `Lt` : Estimates of the Hölder constant for each `t`, corresponding to
  \\L_t^2\\.

- `PN` : The number of selected curves used in the estimation for each
  `t`.

- `muhat` : Estimated values of the mean function at each `t`.

## References

Maissoro H, Patilea V, Vimond M (2024). “Adaptive estimation for Weakly
Dependent Functional Times Series.” *arXiv preprint arXiv:2403.13706*.

## See also

[`estimate_mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_sigma()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md),
[`estimate_nw()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md),
[`estimate_empirical_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_autocov.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Load data
data("data_far")

# Estimate risk function for the mean
dt_mean_risk <- estimate_mean_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = NULL,
  kernel_name = "epanechnikov"
)

# Visualize mean risk at various observation points
dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
manipulateWidget::combineWidgets(
  list = list(
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.25" = `0.25`)],
      main = "t = 0.25", xlab = "h", ylab = "Risk Function"),
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.5" = `0.5`)],
      main = "t = 0.5", xlab = "h", ylab = "Risk Function"),
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.75" = `0.75`)],
      main = "t = 0.75", xlab = "h", ylab = "Risk Function")
  ),
  nrow = 3
)

# Estimate mean function with optimal bandwidths
dt_mean <- estimate_mean(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
  kernel_name = "epanechnikov"
)

# Display rounded estimates of the mean function
DT::datatable(data = dt_mean[, lapply(.SD, function(X) round(X, 3))])
} # }
```
