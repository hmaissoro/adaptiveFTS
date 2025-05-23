# The model
## Let
m <- function(t) 4 * sin(1.5 * pi * t)

## Observation points
t <- runif(n = 200, min = 0, max = 1)
t <- sort(t)

## Measure error
e <- rnorm(n = 200, mean = 0, sd = 0.2)

## Regression model
y <- m(t) + e

plot(x = t, y = y, main = "Observed points and true regression function")
lines(x = t, y = m(t), type = "l", col = "red")

## Estimate the best bandwidth
h_grid <- seq(1 / (2 * length(t)), length(t) ** (- 1/3), len = 100)
hbest <- estimate_nw_bw(y = y, t = t, bw_grid = h_grid,
                        kernel_name = "epanechnikov")

## Estimate the regression function
dt_nw <- estimate_nw(y = y, t = t,
                     tnew = seq(0.01, 0.99, len = 100),
                     h = hbest, kernel_name = "epanechnikov")

plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
     main = "Estimated and true regression function.")
lines(x = dt_nw[, tnew], y = m(dt_nw[, tnew]), type = "l", col = "red")
legend(x = 0.64, y = 4.1, fill = c("blue", "red"),legend = c("Estimated m", "True m"))
