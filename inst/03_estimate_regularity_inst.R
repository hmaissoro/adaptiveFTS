# Generate a sample of FAR(1)
## Exponent H
Hfun <- function(t) {
  hurst_logistic(t = t, h_left = 0.4, h_right = 0.8, slope = 5)
}

## Hölder constant
L <- 4

dt_far <- simulate_far(N = 200L, lambda = 100L,
                       tdesign = "random",
                       Mdistribution = rpois,
                       tdistribution = runif,
                       tcommon = NULL,
                       hurst_fun = Hfun,
                       L = L,
                       far_kernel = function(s,t) 9/4 * exp( - (t + 2 * s) ** 2),
                       far_mean = function(t) 4 * sin(1.5 * pi * t),
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)

# Estimate local regularity at
t0 <- seq(0.2, 0.8, len = 8)

## If data is a data.table or a data. frame
dt_locreg <- estimate_locreg(data = dt_far,
                             idcol = "id_curve",
                             tcol = "tobs",
                             ycol = "X",
                             t = t0,
                             Delta = NULL,
                             h = NULL,
                             smooth_ker = epanechnikov,
                             weighted = FALSE)
DT::datatable(dt_locreg)

## If data is a list of data.table (or data. frame)
list_dt_far <- lapply(unique(dt_far[, id_curve]), function(idx){
  dt_far[id_curve == idx, .(tobs, X)]
})

dt_locreg_2 <- estimate_locreg(data = list_dt_far,
                               idcol = NULL,
                               tcol = "tobs",
                               ycol = "X",
                               t = t0,
                               Delta = NULL,
                               h = NULL,
                               smooth_ker = epanechnikov,
                               weighted = FALSE)
DT::datatable(dt_locreg_2)

## If data is a list of list
list_list_far <- lapply(unique(dt_far[, id_curve]), function(idx){
  list("Obs_point" = dt_far[id_curve == idx, tobs],
       "Xobs" = dt_far[id_curve == idx, X])
})

dt_locreg_3 <- estimate_locreg(data = list_list_far,
                               idcol = NULL,
                               tcol = "Obs_point",
                               ycol = "Xobs",
                               t = t0,
                               Delta = NULL,
                               h = NULL,
                               smooth_ker = epanechnikov,
                               weighted = TRUE)
DT::datatable(dt_locreg_2)

