# Draw a fractional Brownian motion sample path.

Draw a fractional Brownian motion sample path.

## Usage

``` r
simulate_fBm(t = seq(0.2, 0.8, len = 20), hurst = 0.6, L = 1, tied = TRUE)
```

## Arguments

- t:

  `vector (float)`. Grid of points between 0 and 1 where we want to
  generate the sample path.

- hurst:

  `float (positive)`. The Hurst exponent scalar value between 0 and 1.

- L:

  `float (positive)`. Hölder constant.

- tied:

  `boolean`. If `TRUE`, the sample path is tied-down.

## Value

A `data.table` containing 2 column : `t` and `mfBm`, the sample path.

## Examples

``` r

t0 <- seq(0.2, 0.8, len = 20)
dt_fBm <- simulate_fBm(t = t0, hurst = 0.6, L = 1, tied = TRUE)
plot(x = dt_fBm$t, y = dt_fBm$fBm, type = "l", col = "red")

```
