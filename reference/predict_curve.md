# Curve prediction using the Best Linear Unbiased Predictor (BLUP).

**Deprecated.** `predict_curve()` is deprecated and will be removed in a
future release. It reconstructs a curve from a block system conditioning
on the neighbouring curve and the target curve's own partial
observations. For the design-weighted, Tikhonov-regularised adaptive
BLUP (one-step-ahead prediction of the curve following a conditioning
curve), use
[`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)
with
[`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md),
or the one-call wrapper
[`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md).
Note that these compute a different quantity, so results are not
interchangeable.

This function predict a curve using the adaptive Best Linear Unbiased
Predictor proposed by Maissoro, Patilea and Vimond (2026).

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
  common_bw = FALSE,
  center_curves = TRUE,
  correct_diagonal = TRUE,
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

- common_bw:

  A logical value indicating whether a single bandwidth is used for both
  arguments of the (auto)covariance. Default is `FALSE`.

- center_curves:

  A logical value indicating whether the curves are centred before
  smoothing. Default is `TRUE`.

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

Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
Functional Time Series. *arXiv preprint* arXiv:2609.xxxxx.

## See also

[`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md),
[`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md),
[`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md),
[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md).
