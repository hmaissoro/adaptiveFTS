# Generate a FAR A process
dt_far <- simulate_far(N = 50, lambda = 70,
                       tdesign = "random",
                       tdistribution = runif,
                       tcommon = NULL,
                       hurst_fun = hurst_logistic,
                       L = 4,
                       far_kernel = get_real_data_far_kenel,
                       far_mean = get_real_data_mean,
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)

# Add noise
dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]
data <- dt_far[ttag == "trandom"]
