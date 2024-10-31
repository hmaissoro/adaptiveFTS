source("./R/02_smoothing.R")
source("./R/03_estimate_regularity.R")
Rcpp::sourceCpp("./src/03_estimate_locreg_rcpp.cpp")

# Import the data
data('data_far')

# Estimation parameters
t0 <- seq(0.1, 0.9, len = 10)

# Estimation using Rcpp
dt_locreg_cpp <- estimate_locreg_cpp(
  data = data_far,
  t = t0, Delta = NULL, h = NULL,
  kernel_name = "epanechnikov", center = TRUE)

DT::datatable(dt_locreg_cpp)
