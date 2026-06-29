# Estimate the Risk Function of the Mean Function

This function estimates the risk function \\R\_\mu(t;h)\\ for the mean
function estimation as described in Section 4.1 of Maissoro et al.
(2024) .

## Usage

``` r
estimate_mean_risk(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
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

- bw_grid:

  `vector (numeric)`. A bandwidth grid from which the best smoothing
  parameter is selected for each `t`. Default is `NULL`, in which case
  it is defined as an exponential grid of \\N \lambda\\.

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

## Value

A `data.table` with columns:

- `t` : The observation points where the risk function is estimated.

- `h` : The candidate bandwidth values tested.

- `PN` : The number of curves used to estimate the mean at each `t`,
  corresponding to \\P_N(t;h)\\.

- `locreg_bw` : The bandwidth used to estimate the local regularity
  parameters.

- `Ht` : Estimates of the local exponent at each `t`, corresponding to
  \\H_t\\.

- `Lt` : Estimates of the Hölder constant at each `t`, corresponding to
  \\L_t^2\\.

- `bias_term` : The bias term component of the risk function.

- `variance_term` : The variance term component of the risk function.

- `dependence_term` : The dependence term component of the risk
  function.

- `mean_risk` : The estimated risk function for the mean.

## References

Maissoro H, Patilea V, Vimond M (2024). “Adaptive estimation for Weakly
Dependent Functional Times Series.” *arXiv preprint arXiv:2403.13706*.

## See also

[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_sigma()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md),
[`estimate_nw()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md),
[`estimate_empirical_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_autocov.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Load data
data("data_far")

# Estimate the risk function for mean estimation
dt_mean_risk <- estimate_mean_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = NULL,
  kernel_name = "epanechnikov"
)

# Plot the mean risk function at different points
dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
manipulateWidget::combineWidgets(
  list = list(
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.25" = `0.25`)],
      main = "t = 0.25", xlab = "h", ylab = "Risk Function"
    ),
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.5" = `0.5`)],
      main = "t = 0.5", xlab = "h", ylab = "Risk Function"
    ),
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.75" = `0.75`)],
      main = "t = 0.75", xlab = "h", ylab = "Risk Function"
    )
  ),
  nrow = 3
)
} # }
```
