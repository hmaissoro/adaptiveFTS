# Curve prediction using the Best Linear Unbiased Predictor (BLUP).

This function predict a curve using the adaptive Best Linear Unbiased
Predictor proposed by Maissoro et al. (2025) .

## Usage

``` r
predict_curve(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = seq(0.01, 0.99, len = 99),
  id_curve_to_predict = NULL,
  bw_grid = NULL,
  use_same_bw = FALSE,
  center = TRUE,
  correct_diagonal = TRUE,
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

  A numeric vector specifying the time points at which to predict the
  curve `id_curve_to_predict`.

- id_curve_to_predict:

  An integer specifying the index of the curve to be predicted. Default
  is `NULL`, which considers the last curve in `data`.

- bw_grid:

  A numeric vector of bandwidth grid values for selecting optimal
  bandwidth parameters for (auto)covariance estimation. Default is
  `NULL`, which sets it in the function.

- use_same_bw:

  A logical value indicating whether the same bandwidth should be used
  for the arguments `s` and `t` in the (auto)covariance estimation.
  Default is `FALSE`.

- center:

  A logical value indicating whether the data should be centered before
  estimating the (auto)covariance. Default is `TRUE`.

- correct_diagonal:

  A logical value indicating whether the diagonal of the covariances
  should be corrected. Default is `TRUE`.

- kernel_name:

  A string specifying the kernel to use for estimation. Supported values
  are `"epanechnikov"`, `"biweight"`, `"triweight"`, `"tricube"`,
  `"triangular"`, and `"uniform"`. Default is `"epanechnikov"`.

## Value

A `data.table` containing the predicted curve:

- `t` : The time points at which the curve `id_curve_to_predict` is
  predicted.

- `muhat` : The estimates of the mean function.

- `prediction` : The adaptive estimates the Best Linear Unbiased
  Predictor.

## References

Maissoro H, Patilea V, Vimond M (2025). “Adaptive prediction for
Functional Times Series.” *arXiv preprint arXiv:2501.xxxxx*.

## See also

[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
[`estimate_nw()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md).
