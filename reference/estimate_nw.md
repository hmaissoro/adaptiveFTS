# Nadaraya-Watson Kernel Estimator

Estimates the Nadaraya-Watson regression function using a specified
kernel function.

## Usage

``` r
estimate_nw(y, t, tnew, h = NULL, kernel_name = "epanechnikov")
```

## Arguments

- y:

  `vector (numeric)`. A numeric vector containing the observed values of
  the dependent variable corresponding to the observation points `t`.

- t:

  `vector (numeric)`. A numeric vector containing the observed values of
  the independent variable.

- tnew:

  `vector (numeric)`. New `t` values at which to estimate the regression
  function.

- h:

  `numeric (positive)`. The bandwidth parameter. Default is `h = NULL`,
  in which case it will be chosen using cross-validation

- kernel_name:

  `string`. A string specifying the name of the kernel function to use.
  The default is "epanechnikov". Supported kernels: "epanechnikov",
  "biweight", "triweight", "tricube", "triangular", "uniform".

## Value

A `data.table` with the following columns:

- `h` : The bandwidth used to estimate the regression function.

- `tnew` : The vector of new `t` values at which the regression function
  is estimated.

- `yhat` : A vector of estimated values of the regression function at
  `tnew`.

## See also

[`estimate_nw_bw`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw_bw.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define the true regression function
m <- function(t) 4 * sin(1.5 * pi * t)

# Generate observation points
t <- sort(runif(n = 200, min = 0, max = 1))

# Generate measurement errors
e <- rnorm(n = 200, mean = 0, sd = 0.2)

# Observed values of the regression model
y <- m(t) + e

# Plot observed points and true regression function
plot(x = t, y = y, main = "Observed points and true regression function")
lines(x = t, y = m(t), col = "red")

# Estimate optimal bandwidth using cross-validation
bw_grid <- seq(1 / (2 * length(t)), length(t) ** (-1/3), length.out = 100)
hbest <- estimate_nw_bw(y = y, t = t, bw_grid = bw_grid, kernel_name = "epanechnikov")

# Estimate regression function with the Nadaraya-Watson estimator
dt_nw <- estimate_nw(y = y, t = t, tnew = seq(0.01, 0.99, length.out = 100),
                     h = hbest, kernel_name = "epanechnikov")

# Plot estimated and true regression functions
plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
     main = "Estimated and true regression function")
lines(x = dt_nw[, tnew], y = m(dt_nw[, tnew]), col = "red")
legend(x = 0.64, y = 4.1, fill = c("blue", "red"), legend = c("Estimated m", "True m"))
} # }
```
