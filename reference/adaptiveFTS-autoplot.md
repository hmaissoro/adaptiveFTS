# Plot an adaptive functional time series estimator

[ggplot2::autoplot()](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods return a ggplot object for the outputs of the adaptive
estimators; the matching
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods draw
and invisibly return it. `ggplot2` is a suggested dependency and must be
installed. See
[adaptiveFTS_est](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS_est.md)
for the class structure and
[`summary.mean_est()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
for the text summaries.

## Usage

``` r
autoplot.mean_est(object, ...)

# S3 method for class 'mean_est'
plot(x, ...)

autoplot.locreg_est(object, which = c("Ht", "Lt2"), ...)

# S3 method for class 'locreg_est'
plot(x, ...)

autoplot.cov_segment_est(object, ...)

# S3 method for class 'cov_segment_est'
plot(x, ...)

autoplot.autocov_est(object, ...)

# S3 method for class 'autocov_est'
plot(x, ...)

autoplot.mean_risk(object, ...)

# S3 method for class 'mean_risk'
plot(x, ...)

autoplot.cov_segment_risk(object, ...)

# S3 method for class 'cov_segment_risk'
plot(x, ...)

autoplot.autocov_risk(object, ...)

# S3 method for class 'autocov_risk'
plot(x, ...)

autoplot.fts_acf(object, ...)

# S3 method for class 'fts_acf'
plot(x, ...)
```

## Arguments

- object, x:

  An adaptive-estimator object (see
  [adaptiveFTS_est](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS_est.md)).

- ...:

  Passed to the corresponding `autoplot` method (for `plot`) or unused.

- which:

  For `locreg_est`, the regularity parameters to display; a subset of
  `c("Ht", "Lt2")`. Default both.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object (`plot` methods return it invisibly, after drawing).

## See also

[adaptiveFTS-summary](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md),
[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md).

## Examples

``` r
if (FALSE) { # \dontrun{
data("data_far")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::autoplot(estimate_mean(data = data_far,
                                  t = seq(0.1, 0.9, length.out = 20)))
}
} # }
```
