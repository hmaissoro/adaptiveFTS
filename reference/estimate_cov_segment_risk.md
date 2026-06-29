# Estimate the Risk of the Covariance Segment Function

Estimates the risk \\R\_{\Gamma_0}(t; h)\\ associated with the
covariance segment line estimation proposed by Maissoro et al. (2025) .

## Usage

``` r
estimate_cov_segment_risk(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  bw_grid = NULL,
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

- t:

  A numeric vector. Observation points where the mean function of the
  underlying process is estimated.

- bw_grid:

  A numeric vector. A bandwidth grid from which the best smoothing
  parameter is selected for each `t`. Default is `NULL`, in which case
  it is defined as an exponential grid of \\N \times \lambda\\.

- center:

  Logical. If `TRUE`, centers the data before estimation. Default is
  `TRUE`.

- kernel_name:

  Character string. Specifies the kernel function for estimation;
  default is `"epanechnikov"`. Supported kernels include:
  `"epanechnikov"`, `"biweight"`, `"triweight"`, `"tricube"`,
  `"triangular"`, and `"uniform"`.

## Value

A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) with
columns:

- `t`: The points at which the risk function is estimated.

- `h`: The candidate bandwidth.

- `PN`: The number of curves used to estimate the mean at `t`, i.e.,
  \\P_N(t;h)\\.

- `locreg_bw`: The bandwidth used to estimate the local regularity
  parameters.

- `Ht`: The estimates of the local exponent \\H_t\\.

- `Lt`: The estimates of the Hölder constant \\L_t^2\\.

- `bias_term`: The bias term of the risk function.

- `variance_term`: The variance term of the risk function.

- `dependence_term`: The dependence term of the risk function.

- `cov_segment_risk`: The estimated risk of the covariance segment
  function.

## Details

The local regularity parameters are estimated within the function using
`estimate_locreg_cpp`.

The dependence term includes contributions from both a term based on
\\\mathbb{D}(t; h_t)\\ derived from fourth-moment tensors, and an
empirical autocovariance term computed using
`estimate_empirical_XsXt_autocov_cpp`.

## References

Maissoro H, Patilea V, Vimond M (2025). “Adaptive prediction for
Functional Times Series.” *arXiv preprint arXiv:2501.xxxxx*.

## See also

[estimate_mean](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[estimate_locreg](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[estimate_sigma](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md),
[estimate_nw](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md),
[estimate_empirical_autocov](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_autocov.md)

## Examples

``` r
# Example coming soon
```
