# Logistic Hurst function

Logistic Hurst function that can be used to generate multifractional
Brownian motion (mfBm). See the following paper
https://doi.org/10.3390/fractalfract6020074.

## Usage

``` r
hurst_logistic(
  t,
  h_left = 0.2,
  h_right = 0.8,
  slope = 30,
  change_point_position = 0.5
)
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

- slope:

  `Float (positive)`. A scalar positive value corresponding to the slope
  of the logistic function.

- change_point_position:

  `Float`. A scalar value in the interval between 0 and 1 corresponding
  ti the change point position.

## Value

A `vector (float)` corresponding to the value of the function evaluated
at `t`.

## See also

[`hurst_arctan()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_arctan.md),
[`hurst_linear()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_linear.md).

## Examples

``` r
t0 <- seq(0.2, 0.8, len = 10)
hlogistic <- hurst_logistic(t = t0, h_left = 0.2,
                            h_right = 0.8, slope = 30,
                            change_point_position = 0.5)
plot(x = t0, y = hlogistic, type = "b", col = "red")


```
