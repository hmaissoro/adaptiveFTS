# Current function
source("./R/02_smoothing.R")
Rcpp::sourceCpp("./src/02_smoothing_rcpp.cpp")

# The data
n <- 1000
Tn <- runif(n)
Yn <- sin(2 * pi * Tn) + rnorm(n = n, mean = 0, sd = 0.25)
plot(Tn, Yn)

# Function
estimate_nw_wrap <- function(y, t, tnew, h){
  if (! is.numeric(y) | ! is.numeric(t) | ! is.numeric(tnew) |! is.numeric(h))
    stop("The arguments 'y', 't', 'tnew', 'h' must be numeric.")
  if (length(t) != length(y) & length(t) < 2)
    stop("The arguments 'y' and 't' must have a length of at least 2 and must be of the same length.")
  if (! is.null(h) & length(h) > 1)
    stop("The bandwidth 'h' must be either NULL or saclar.")

  # Estimate Nadaraya-Watson estimator using estimate_nw_cpp function
  res <- estimate_nw_cpp(y = Yn, t = Tn, tnew = tnew, h = h)

  # convert
  res <- data.table::data.table("h" = h, "tnew" = tnew, "yhat" = c(res))
  return(res)
}

# Comparaison R vs Rcpp N-W estimator
res_nw <- estimate_nw(y = Yn, t = Tn, tnew = seq(0.1, 0.9, len = 100), h = 0.05)
res_nw_cpp <- estimate_nw_cpp(y = Yn, t = Tn, tnew = seq(0.1, 0.9, len = 100), h = 0.05)
res_nw_cpp_wrap <- estimate_nw_wrap(y = Yn, t = Tn, tnew = seq(0.1, 0.9, len = 100), h = 0.05)
all.equal(target = res_nw$yhat, current = res_nw_cpp[, 1])

res_benchmark <- microbenchmark::microbenchmark(
  res_nw <- estimate_nw(y = Yn, t = Tn, tnew = seq(0.1, 0.9, len = 100), h = 0.05),
  res_nw_cpp <- estimate_nw_cpp(y = Yn, t = Tn, tnew = seq(0.1, 0.9, len = 100), h = 0.05),
  res_nw_cpp_wrap <- estimate_nw_wrap(y = Yn, t = Tn, tnew = seq(0.1, 0.9, len = 100), h = 0.05),
  times = 50
)

print(res_benchmark, unit="relative", order="median")

# Comparaison R vs Rcpp of bandwidth selection
K <- 20
b0 <- 1 / length(Tn)
bK <- length(Tn) ** (- 1 / 3)
a <- exp((log(bK) - log(b0)) / K)
bw_grid <- b0 * a ** (seq_len(K))
rm(b0, bK, a, K) ; gc()

hbest <- estimate_nw_bw(y = Yn, t = Tn, bw_grid = bw_grid)
hbest_cpp <- estimate_nw_bw_cpp(y = Yn, t = Tn, bw_grid = bw_grid)
hbest == hbest_cpp

res_bm_bw_sel <- microbenchmark::microbenchmark(
  current_nw_bw_sel = hbest <- estimate_nw_bw(y = Yn, t = Tn, bw_grid = bw_grid),
  cpp_nw_bw_sel = hbest_cpp <- estimate_nw_bw_cpp(y = Yn, t = Tn, bw_grid = bw_grid),
  times = 50
)

print(res_bm_bw_sel, unit="relative", order="median")

# Estimate best bandwidth on a subset of curves
# Import the data
data("data_far")
get_nw_optimal_bw_cpp(data_far, bw_grid = NULL, nsubset = 30, kernel_name = "epanechnikov")

