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
data("data_far")
dt_small <- data_far[data_far$id_curve <= 20, ]
bwg <- seq(0.04, 0.15, length.out = 5)

summary(estimate_mean(data = dt_small, t = seq(0.1, 0.9, length.out = 9),
                      bw_grid = bwg))
#> Adaptive mean function estimate
#>   Evaluation points  : 9 (t in [0.1, 0.9])
#>   Training curves    : 20
#>   Kernel             : epanechnikov
#>   Optimal bandwidth  : [0.04, 0.04]
#>   Curves used (PN)   : [20, 20]
#>   muhat              : [240, 244]
summary(estimate_autocov(data = dt_small, s = c(1/5, 2/5), t = c(1/4, 1/2),
                         lag = 1, bw_grid = bwg))
#> Adaptive autocovariance estimate (lag = 1)
#>   (s, t) pairs       : 2 (0 on the diagonal s = t)
#>   Training curves    : 20
#>   Kernel             : epanechnikov (centred: TRUE)
#>   Common bw for s, t : FALSE
#>   Bandwidth (s)      : [0.04, 0.04]
#>   Bandwidth (t)      : [0.04, 0.04]
#>   Curves used (PNl)  : [19, 19]
#>   autocov            : [-0.175, 0.936]
```
