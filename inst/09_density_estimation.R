## Design-density estimation demo: the leave-one-out Parzen-Rosenblatt estimator
## (estimate_density) and the subset bandwidth selection (get_density_optimal_bw)
## used for the independent-design weights of the adaptive BLUP.

library(adaptiveFTS)
library(data.table)

data("data_far")
dt <- format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")

Tn0 <- dt[id_curve == 150, sort(unique(tobs))]

## Bandwidth by least-squares cross-validation on a single curve.
fit <- estimate_density(x = Tn0, bw_grid = exp(seq(log(0.02), log(0.5), length.out = 30)))
fit$h_star

plot(fit$bw_grid, fit$cv_curve, type = "b", log = "x",
     xlab = "bandwidth", ylab = "LSCV score", main = "Design-density LSCV")
abline(v = fit$h_star, col = "red", lty = 2)

plot(Tn0, fit$estimate, type = "o", pch = 16, cex = 0.5,
     xlab = "t", ylab = "density", main = "Leave-one-out density estimate")

## Design weights.
rho <- 1 / (length(Tn0) * pmax(fit$estimate, 1e-6))
rho <- rho / sum(rho)

## A single bandwidth selected on a subset of curves, reused via the fixed-h path.
h_sub <- get_density_optimal_bw(data = dt, nsubset = 30L)
ghat  <- estimate_density(x = Tn0, h = h_sub)$estimate
