# Predict with an adaptive functional BLUP fit

Evaluates the adaptive Best Linear Unbiased Predictor of the curve
following the conditioning curve at the requested prediction points. For
`horizon > 1` (multi-step-ahead), each intermediate curve is predicted
on the target grid `t` and fed back as the new conditioning curve, so
the `(n_0 + i)`-th prediction is built from the `(n_0 + i - 1)`-th
predicted curve, and all intermediate horizons are returned.

## Usage

``` r
# S3 method for class 'blup_fit'
predict(object, t = object$Tn0, horizon = 1L, newdata = NULL, ...)
```

## Arguments

- object:

  A `blup_fit` object.

- t:

  Numeric vector of prediction points in \\\[0, 1\]\\. Default is the
  conditioning-curve design points.

- horizon:

  Integer prediction horizon (steps ahead). Default `1`. For
  `horizon > 1`, each intermediate curve is predicted on the target grid
  `t` and fed back as the new conditioning curve; the conditioning
  quantities are recomputed for that grid using the cached adaptive
  bandwidths (works for both designs).

- newdata:

  Optional numeric vector of conditioning-curve values at the fit's
  design points (`object$Tn0`), overriding `object$Yn0`. Used internally
  for the multi-step recursion; must have length `object$Mn0`.

- ...:

  Unused; for S3 compatibility.

## Value

A `data.table` with one block of rows per horizon (`horizon * length(t)`
rows in total):

- `horizon`: the prediction step, from `1` to `horizon`.

- `t`: the prediction points.

- `muhat`: the mean estimate at `t`.

- `prediction`: the adaptive BLUP at `t` for that horizon.

## Details

During multi-step prediction the estimation data is left unchanged: the
mean, the (auto)covariance operators, the adaptive bandwidths and the
noise level are estimated once and held fixed, and each predicted curve
enters only as the conditioning values of the next step. Injecting a
predicted (denoised) curve back into the estimation sample would bias
those plug-in estimates, so it is deliberately avoided.

## See also

[`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md).
