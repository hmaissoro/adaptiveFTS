# Summarise an adaptive functional time series estimator

Compact, design-aware text summaries for the objects returned by the
adaptive estimators of adaptiveFTS. Each method prints a short report
and returns its argument invisibly. See
[adaptiveFTS_est](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS_est.md)
for the class structure and
[`autoplot.mean_est()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
for the companion plots.

## Usage

``` r
# S3 method for class 'locreg_est'
summary(object, ...)

# S3 method for class 'mean_est'
summary(object, ...)

# S3 method for class 'cov_segment_est'
summary(object, ...)

# S3 method for class 'autocov_est'
summary(object, ...)

# S3 method for class 'mean_risk'
summary(object, ...)

# S3 method for class 'cov_segment_risk'
summary(object, ...)

# S3 method for class 'autocov_risk'
summary(object, ...)

# S3 method for class 'adaptiveFTS_est'
summary(object, ...)

# S3 method for class 'fts_acf'
summary(object, ...)
```

## Arguments

- object:

  An adaptive-estimator object (see
  [adaptiveFTS_est](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS_est.md)).

- ...:

  Unused; for S3 compatibility.

## Value

`object`, invisibly.

## See also

[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md),
[adaptiveFTS_est](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS_est.md).

## Examples

``` r
if (FALSE) { # \dontrun{
data("data_far")
summary(estimate_mean(data = data_far, t = seq(0.1, 0.9, length.out = 9)))
summary(estimate_autocov(data = data_far,
                         s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 1))
} # }
```
