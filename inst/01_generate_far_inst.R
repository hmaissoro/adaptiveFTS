# Generate a sample of 2 FAR
dt_far <- simulate_far(N = 2L, lambda = 70L,
                       tdesign = "common",
                       tdistribution = NULL,
                       tcommon = seq(0.2, 0.8, len = 50),
                       hurst_fun = hurst_logistic,
                       L = 4,
                       far_kernel = function(s,t) 9/4 * exp( - (t + 2 * s) ** 2),
                       far_mean = function(t) 4 * sin(1.5 * pi * t),
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)
dt_far
