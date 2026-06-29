# Linear Hurst function

Linear Hurst function that can be used to generate multifractional
Brownian motion (mfBm). See the following paper
https://doi.org/10.3390/fractalfract6020074.

## Usage

``` r
hurst_linear(t = seq(0.2, 0.8, len = 10), h_left = 0.2, h_right = 0.8)
```

## Arguments

- t:

  `vector (float)`. Points between 0 and 1 at which to evaluate the
  function.

- h_left:

  `Float`. A scalar value in the interval between 0 and 1 indicating the
  minimum of the function.

- h_right:

  `Float`. A scalar value in the interval between 0 and 1 indicating the
  maximum of the function.

## Value

A `vector (float)` corresponding to the value of the function evaluated
at `t`.

## See also

[`hurst_arctan()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_arctan.md),
[`hurst_logistic()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_logistic.md).

## Examples

``` r
t0 <- seq(0.2, 0.8, len = 10)
hlinear <- hurst_linear(t = t0)
plot(x = t0, y = hlinear, type = "b", col = "red")


```
