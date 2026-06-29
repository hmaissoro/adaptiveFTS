# Nadaraya-Watson Bandwidth Selection using Cross-Validation

Selects the optimal bandwidth for Nadaraya-Watson kernel regression
using cross-validation.

## Usage

``` r
estimate_nw_bw(y, t, bw_grid = NULL, kernel_name = "epanechnikov")
```

## Arguments

- y:

  `vector (numeric)`. A numeric vector containing the observed values of
  the dependent variable.

- t:

  `vector (numeric)`. A numeric vector containing the observed values of
  the independent variable.

- bw_grid:

  `vector (numeric)`. A grid of bandwidth values to test. Default is
  `bw_grid = NULL`, in which case an exponential grid based on the
  length of `t` will be used.

- kernel_name:

  `string`. A string specifying the name of the kernel function to use.
  The default is "epanechnikov". Supported kernels: "epanechnikov",
  "biweight", "triweight", "tricube", "triangular", "uniform".

## Value

A `numeric` value representing the optimal bandwidth that minimizes the
cross-validation error.

## Details

This function computes the optimal bandwidth for the Nadaraya-Watson
kernel regression estimator by performing cross-validation over a set of
candidate bandwidths defined by `bw_grid`. The function returns the
bandwidth that minimizes the cross-validation error.

## See also

[`estimate_nw`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md)

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

# Define a grid of candidate bandwidths
bw_grid <- seq(1 / (2 * length(t)), length(t)^(-1/3), length.out = 100)

# Estimate the best bandwidth using cross-validation
hbest <- estimate_nw_bw(y = y, t = t, bw_grid = bw_grid, kernel_name = "epanechnikov")

# Estimate the regression function using the selected bandwidth
dt_nw <- estimate_nw(y = y, t = t, tnew = seq(0.01, 0.99, length.out = 100),
                     h = hbest, kernel_name = "epanechnikov")

# Plot estimated and true regression functions
plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
     main = "Estimated and true regression function")
lines(x = dt_nw[, tnew], y = m(dt_nw[, tnew]), col = "red")
legend(x = 0.64, y = 4.1, fill = c("blue", "red"), legend = c("Estimated m", "True m"))
} # }
```
