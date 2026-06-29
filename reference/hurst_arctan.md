# Arctan Hurst function

Arctan Hurst function that can be used to generate multifractional
Brownian motion (mfBm). See the following paper
https://doi.org/10.3390/fractalfract6020074.

## Usage

``` r
hurst_arctan(t = seq(0.2, 0.8, len = 10))
```

## Arguments

- t:

  `vector (float)`. Points between 0 and 1 at which to evaluate the
  function.

## Value

A `vector (float)` corresponding to the value of the function evaluated
at `t`.

## See also

[`hurst_linear()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_linear.md),
[`hurst_logistic()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_logistic.md).

## Examples

``` r

t0 <- seq(0.2, 0.8, len = 10)
htan <- hurst_arctan(t = t0)
plot(x = t0, y = htan, type = "b", col = "red")


```
