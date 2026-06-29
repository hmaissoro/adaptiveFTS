# Mean function learned from the voltage curves of the electricity

For more details see the vignette:
`vignette("hybrid-simulation-setup", package = "adaptiveFTS")`

## Usage

``` r
get_real_data_mean(t = seq(0.1, 0.9, len = 10))
```

## Arguments

- t:

  `vector (numeric)`. Points at which we want to return the mean
  function. It can be a scalar.

## Value

A `data.table` containing 2 columns.

- t : The vector or scalar `t`.

- mean : The values of the mean function evaluated at `t`.

## Examples

``` r

t0 <- seq(0.1, 0.9, len = 10)
m <- get_real_data_mean(t = t0)
plot(x = t0, y = m, type = "b", col = "red",
     xlab = "t", ylab = "mean", main = "Mean function")

```
