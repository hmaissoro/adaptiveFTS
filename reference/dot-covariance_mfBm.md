# Covariance matrix of the multi-fractional Brownian Motion

Covariance matrix of the multi-fractional Brownian Motion

## Usage

``` r
.covariance_mfBm(t = seq(0.2, 0.8, len = 10), hurst_fun = hurst_logistic, ...)
```

## Arguments

- t:

  `vector (float)`. Points between 0 and 1 at which to compute the
  covariance function.

- hurst_fun:

  `function`. Hurst function. It can be
  [`hurst_arctan`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_arctan.md),
  [`hurst_linear`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_linear.md),
  [`hurst_logistic`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_logistic.md).

- ...:

  Hurst function additional arguments.

## Value

a `matrix` of `t` x `t` covariance.
