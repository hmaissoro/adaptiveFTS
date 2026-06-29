# Draw a multifractional Brownian motion sample path.

This function generates a sample path of a multifractional Brownian
motion (mfBm) based on the provided Hurst function and other parameters.

## Usage

``` r
simulate_mfBm(
  t = seq(0.2, 0.8, len = 50),
  hurst_fun = hurst_logistic,
  L = 1,
  shift_var = 0,
  tied = TRUE,
  ...
)
```

## Arguments

- t:

  `vector (float)`. Grid of points between 0 and 1 where the sample path
  will be generated.

- hurst_fun:

  `function`. Hurst function. It can be
  [`hurst_arctan`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_arctan.md),
  [`hurst_linear`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_linear.md),
  [`hurst_logistic`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_logistic.md),
  or any custom Hurst function.

- L:

  `float (positive)`. Hölder constant.

- shift_var:

  `float (positive)`. The variance of the shift Gaussian random
  variable. Default is `shift_var = 1`, meaning a normal random variable
  with mean 0 and variance 1 is added.

- tied:

  `boolean`. If `TRUE`, the sample path is tied down.

- ...:

  Additional arguments for the Hurst function.

## Value

A `data.table` containing 2 columns: `t` and `mfBm`, representing the
grid points and the corresponding values of the mfBm sample path.

## Examples

``` r
t0 <- seq(0.2, 0.8, len = 20)
dt_mfBm <- simulate_mfBm(t = t0, hurst_fun = hurst_logistic, L = 1, tied = TRUE)
plot(x = dt_mfBm$t, y = dt_mfBm$mfBm, type = "l", col = "red")

```
