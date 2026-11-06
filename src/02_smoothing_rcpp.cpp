#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "00_internal.h"
using namespace Rcpp;
using namespace arma;

//' Biweight kernel function
 //'
 //' Calculate the biweight kernel function.
 //'
 //' @param u A numeric vector.
 //' @return A numeric vector representing the biweight kernel function.
 //'
 //' @seealso [triweight_kernel()], [tricube_kernel()], [epanechnikov_kernel()], [triangular_kernel()], and [uniform_kernel()].
 //'
 //' @examples
 //' biweight_kernel(c(-1, 0, 1))
 //'
 //'
 // [[Rcpp::export]]
 arma::vec biweight_kernel(const arma::vec u) {
   // Check if u is not empty
   if(u.size() == 0) {
     Rcpp::stop("Input vector is empty");
   }
   arma::vec result = 15.0/16.0 * arma::pow(1 - arma::square(u), 2);
   result.elem(arma::find(abs(u) > 1)).fill(0);
   return result;
 }

 //' Triweight kernel function
 //'
 //' Calculate the Triweight kernel function.
 //'
 //' @param u A numeric vector.
 //' @return A numeric vector representing the Triweight kernel function.
 //'
 //' @seealso [biweight_kernel()], [tricube_kernel()], [epanechnikov_kernel()], [triangular_kernel()], and [uniform_kernel()].
 //'
 //' @examples
 //' triweight_kernel(c(-1, 0, 1))
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 arma::vec triweight_kernel(const arma::vec u) {
   // Check if u is not empty
   if(u.size() == 0) {
     Rcpp::stop("Input vector is empty");
   }
   arma::vec result = 35.0/32.0 * arma::pow(1 - arma::square(u), 3);
   result.elem(arma::find(abs(u) > 1)).fill(0);
   return result;
 }

 //' Tricube kernel function
 //'
 //' Calculate the Tricube kernel function.
 //'
 //' @param u A numeric vector.
 //' @return A numeric vector representing the Tricube kernel function.
 //'
 //' @seealso [biweight_kernel()], [triweight_kernel()], [epanechnikov_kernel()], [triangular_kernel()], and [uniform_kernel()].
 //'
 //' @examples
 //'
 //' tricube_kernel(c(-1, 0, 1))
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 arma::vec tricube_kernel(const arma::vec u) {
   // Check if u is not empty
   if(u.size() == 0) {
     Rcpp::stop("Input vector is empty");
   }
   arma::vec result = 70.0/81.0 * arma::pow(1 - arma::pow(arma::abs(u), 3), 3);
   result.elem(arma::find(abs(u) > 1)).fill(0);
   return result;
 }

 //' Epanechnikov kernel function
 //'
 //' Calculate the Epanechnikov kernel function.
 //'
 //' @param u A numeric vector.
 //' @return A numeric vector representing the Epanechnikov kernel function.
 //'
 //' @seealso [biweight_kernel()], [triweight_kernel()], [tricube_kernel()], [epanechnikov_kernel()], and [uniform_kernel()].
 //'
 //' @examples
 //'
 //' epanechnikov_kernel(c(-1, 0, 1))
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 arma::vec epanechnikov_kernel(const arma::vec u) {
   // Check if u is not empty
   if(u.size() == 0) {
     Rcpp::stop("Input vector is empty");
   }
   arma::vec result = 0.75 * (1 - arma::square(u));
   result.elem(arma::find(abs(u) > 1)).fill(0);
   return result;
 }

 //' Triangular kernel function
 //'
 //' Calculate the Triangular kernel function.
 //'
 //' @param u A numeric vector.
 //' @return A numeric vector representing the Triangular kernel function.
 //'
 //' @seealso [biweight_kernel()], [triweight_kernel()], [tricube_kernel()], [epanechnikov_kernel()], and [uniform_kernel()].
 //'
 //' @examples
 //'
 //' triangular_kernel(c(-1, 0, 1))
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 arma::vec triangular_kernel(const arma::vec u) {
   // Check if u is not empty
   if(u.size() == 0) {
     Rcpp::stop("Input vector is empty");
   }
   arma::vec result = (1 - arma::abs(u));
   result.elem(arma::find(abs(u) > 1)).fill(0);
   return result;
 }
 //' Uniform kernel function
 //'
 //' Calculate the Uniform kernel function.
 //'
 //' @param u A numeric vector.
 //' @return A numeric vector representing the Uniform kernel function.
 //'
 //' @seealso [biweight_kernel()], [triweight_kernel()], [tricube_kernel()], [epanechnikov_kernel()], and [triangular_kernel()].
 //'
 //' @examples
 //'
 //' triangular_kernel(c(-1, 0, 1))
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 arma::vec uniform_kernel(const arma::vec u) {
   // Check if u is not empty
   if(u.size() == 0) {
     Rcpp::stop("Input vector is empty");
   }
   arma::vec result = 0.5 * arma::ones(u.size());
   result.elem(arma::find(abs(u) > 1)).fill(0);
   return result;
 }

 // Function to select kernel based on its name
 std::function<arma::vec(const arma::vec)> select_kernel(const std::string kernel_name = "epanechnikov") {
   if (kernel_name == "epanechnikov") {
     return epanechnikov_kernel;
   } else if (kernel_name == "biweight") {
     return biweight_kernel;
   } else if (kernel_name == "triweight") {
     return triweight_kernel;
   } else if (kernel_name == "tricube") {
     return tricube_kernel;
   } else if (kernel_name == "triangular") {
     return triangular_kernel;
   } else if (kernel_name == "uniform") {
     return uniform_kernel;
   } else {
     stop("Unknown kernel name");
   }
 }

 //' Nadaraya-Watson Kernel Estimation
 //'
 //' Estimates the Nadaraya-Watson kernel estimator for regression.
 //'
 //' @param y A vector of response values.
 //' @param t A vector of predictor values corresponding to the response values.
 //' @param tnew A vector of new predictor values for which to estimate the response.
 //' @param h Numeric (positive vector or scalar). The bandwidth of the Nadaraya-Watson estimator for local regularity estimation.
 //' If \code{h} is \code{NULL}, the bandwidth will be estimated by cross-validation on a subset of curves.
 //' If \code{h} is a scalar, all points will be estimated with the same bandwidth.
 //' If \code{h} is a vector, its length must be equal to the number of points in \code{t}, with each element corresponding to a point in \code{t}.
 //' @param kernel_name A string specifying the name of the kernel function to use. The default is "epanechnikov".
 //' Supported kernels: "epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform".
 //'
 //' @return A vector of estimated response values for the new predictor values in \code{tnew}.
 //'
 //' @details This function estimates the Nadaraya-Watson kernel estimator for regression,
 //' using the specified kernel function.
 //'
 //' @examples
 //' y <- c(1, 2, 3, 4, 5)
 //' t <- c(0.1, 0.2, 0.3, 0.4, 0.5)
 //' tnew <- seq(0.1, 0.5, by = 0.1)
 //' yhat <- estimate_nw_cpp(y, t, tnew, h = 0.1)
 //' yhat <- as.vector(yhat)
 //' plot(x = tnew, y = yhat)
 //'
 // [[Rcpp::export]]
 arma::vec estimate_nw_cpp(const arma::vec y, const arma::vec t,
                           const arma::vec tnew, const arma::vec h,
                           const std::string kernel_name = "epanechnikov") {
   // Check if y and t are numeric vectors and not empty
   if (!y.is_finite() || y.size() == 0 ||
       !t.is_finite() || t.size() == 0 ||
       y.size() != t.size()) {
     stop("The vectors y and t must be non-empty numeric vectors of the same length.");
   }
   int m = tnew.size();

   // Check if h is non-negative
   // Control the presmoothing bandwidth
   arma::vec h_to_use;
   int h_size = h.size();
   if (h_size != m) {
     if (h_size != 1) {
       Rcpp::stop("'h' must be a numeric vector with length equal to the number of points in 'tnew' or a scalar value.");
     } else {
       h_to_use = arma::vec(m, fill::value(h(0)));
     }
   } else {
     h_to_use = h;
   }
   // Check if kernel_name is one of the supported kernels
   std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                 "tricube", "triangular", "uniform"};
   if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
     stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
   }

   // The design is sorted once; each estimate only visits the points that can
   // lie in the kernel support (compact on [-1, 1]).
   core::Kern kern = core::parse_kernel(kernel_name);
   std::vector<double> ts;
   std::vector<arma::uword> perm;
   core::sort_design(t, ts, perm);
   arma::vec yhat(m);
   for (int i = 0; i < m; i++)
     yhat(i) = core::nw_sorted(ts, perm, y.memptr(), tnew(i), h_to_use(i), kern);

   return yhat;
 }


 /***R
 n <- 100
 Tn <- runif(n)
 Yn <- sin(2 * pi * Tn) + rnorm(n = n, mean = 0, sd = 0.25)
 estimate_nw_cpp(y = Yn, t = Tn, tnew = seq(0.1, 0.9, len = 100), h = 0.05)
 */

// Leave-one-out cross-validation of the Nadaraya-Watson bandwidth on one curve:
// for each candidate h, the residuals (y_i - yhat_i) / (1 - K(0) / S_i), with
// S_i = sum_j K((t_j - t_i) / h), non-finite ones dropped; the first minimiser
// of the mean squared metric (Inf when no residual is finite).
double nw_bw_impl(const arma::vec& y, const arma::vec& t, const arma::vec& bw_grid, core::Kern kern) {
  const int n = t.size();
  const int bw_size = bw_grid.size();
  std::vector<double> ts;
  std::vector<arma::uword> perm;
  core::sort_design(t, ts, perm);
  const double K0 = core::kval(kern, 0.0);
  arma::vec cv_error(bw_size);
  for (int i = 0; i < bw_size; i++) {
    const double h = bw_grid(i);
    double sum2 = 0;
    int nfin = 0;
    for (int j = 0; j < n; j++) {
      double S = 0;
      double yhat = core::nw_sorted(ts, perm, y.memptr(), t(j), h, kern, &S);
      double metric = (y(j) - yhat) / (1 - K0 / S);
      if (std::isfinite(metric)) { sum2 += metric * metric; ++nfin; }
    }
    cv_error(i) = (nfin == 0) ? arma::datum::inf : sum2 / nfin;
  }
  arma::uword idx_min = arma::index_min(cv_error);
  return bw_grid(idx_min);
}

 //' Estimate Bandwidth for Nadaraya-Watson Kernel Regression by Cross-Validation
  //'
  //' This function estimates the bandwidth parameter for the Nadaraya-Watson kernel regression
  //' using cross-validation.
  //'
  //' @param y A numeric vector representing the response variable.
  //' @param t A numeric vector representing the predictor variable.
  //' @param bw_grid A numeric vector containing the grid of bandwidth values to test.
  //' @param kernel_name A character string specifying the name of the kernel function to use.
  //'        Default is "epanechnikov". Supported kernels: "epanechnikov", "biweight", "triweight",
  //'        "tricube", "triangular" and "uniform".
  //'
  //' @return A numeric value representing the estimated bandwidth parameter for Nadaraya-Watson
  //'         kernel regression.
  //'
  //' @details This function estimates the bandwidth parameter for Nadaraya-Watson kernel regression
  //'          using cross-validation. It tests different bandwidth values specified by \code{bw_grid}
  //'          and selects the one that minimizes the cross-validation error.
  //'
  //' @examples
  //' \dontrun{
  //' # Generate data for the model
  //' m <- function(t) 4 * sin(1.5 * pi * t)
  //' set.seed(123)  # for reproducibility
  //' t <- runif(n = 200, min = 0, max = 1)
  //' t <- sort(t)
  //' e <- rnorm(n = 200, mean = 0, sd = 0.2)
  //' y <- m(t) + e
  //'
  //' # Plot observed points and true regression function
  //' plot(x = t, y = y, main = "Observed points and true regression function")
  //' lines(x = t, y = m(t), type = "l", col = "red")
  //'
  //' # Estimate the best bandwidth using cross-validation
  //' bw_grid <- seq(1 / (2 * length(t)), length(t) ** (- 1/3), len = 100)
  //' hbest <- estimate_nw_bw_cpp(y = y, t = t, bw_grid = bw_grid, kernel_name = "epanechnikov")
  //'
  //' # Estimate the regression function using the selected bandwidth
  //' tnew <- seq(0.01, 0.99, len = 100)
  //' dt_nw <- estimate_nw_cpp(y = y, t = t, tnew = tnew, h = hbest, kernel_name = "epanechnikov")
  //'
  //' # Plot estimated and true regression function
  //' plot(x = tnew, y = dt_nw, type = "l", col = "blue", main = "Estimated and true regression function.")
  //' lines(x = tnew, y = m(tnew), type = "l", col = "red")
  //' legend(x = 0.64, y = 4.1, fill = c("blue", "red"), legend = c("Estimated m", "True m"))
  //'
  //' }
  //'
  //' @export
  //'
  // [[Rcpp::export]]
  double estimate_nw_bw_cpp(const arma::vec y, const arma::vec t,
                            const arma::vec bw_grid,
                            const std::string kernel_name = "epanechnikov") {
    // Check if y and t are numeric vectors and not empty
    if (!y.is_finite() || y.size() == 0 ||
        !t.is_finite() || t.size() == 0 ||
        y.size() != t.size()) {
      stop("The vectors y and t must be non-empty numeric vectors of the same length.");
    }

    // Check if bw_grid is a numeric vector and not empty
    if (!bw_grid.is_finite() || bw_grid.size() == 0) {
      stop("bw_grid must be a non-empty numeric vector.");
    }

    // Check if kernel_name is one of the supported kernels
    std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                  "tricube", "triangular", "uniform"};
    if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
      stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
    }

    return nw_bw_impl(y, t, bw_grid, core::parse_kernel(kernel_name));
  }

// Presmoothing bandwidth: median of the cross-validated NW bandwidths of the
// last nsubset curves (default min(70, floor(N / 2))).
double nw_optimal_bw_impl(const core::Curves& C, const Nullable<arma::vec> bw_grid,
                          const Nullable<int> nsubset, core::Kern kern) {
  int n_curve = C.N;
  double nsubset_val = 0;
  if (nsubset.isNotNull()) {
    nsubset_val = as<int>(nsubset);
  } else {
    nsubset_val = std::min(70.0, std::floor(n_curve * 0.5));
  }
  if (nsubset.isNotNull() && (nsubset_val > n_curve)) {
    stop("If 'nsubset' is not NULL, it must be a positive integer less than or equal to the number of curves.");
  }
  arma::vec bw_grid_used;
  if (bw_grid.isNull()) {
    double lambdahat = C.lambdahat();
    double b0 = 1 / lambdahat;
    double bK = std::pow(lambdahat, -1 / 3.0);
    bw_grid_used = arma::logspace(log10(b0), log10(bK), 15);
  } else {
    bw_grid_used = as<arma::vec>(bw_grid);
  }
  if (!bw_grid_used.is_finite() || bw_grid_used.size() == 0)
    stop("bw_grid must be a non-empty numeric vector.");
  arma::vec res_opt_bw(nsubset_val);
  for (int i = 0; i < nsubset_val; ++i) {
    const int c = n_curve - 1 - i;
    if (!C.Y[c].is_finite() || C.Y[c].n_elem == 0)
      stop("The vectors y and t must be non-empty numeric vectors of the same length.");
    res_opt_bw(i) = nw_bw_impl(C.Y[c], C.T[c], bw_grid_used, kern);
  }
  return arma::median(res_opt_bw);
}

 //' Estimate Optimal Bandwidth for Nadaraya-Watson Estimator
  //'
  //' Estimates the median of the best bandwidths obtained by Nadaraya-Watson cross-validation on a subset of curves.
  //'
  //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X". Typically, the output of the function \link{format_data}.
  //' @param bw_grid A vector of candidate bandwidths. If NULL, a default grid is used.
  //' @param nsubset An integer specifying the number of curves to use for bandwidth selection. If NULL, defaults to 30 or the number of curves if less than 30.
  //' @param kernel_name A string specifying the kernel function to use. Default is "epanechnikov".
  //'
  //' @return A double value representing the estimated optimal bandwidth.
  //'
  //' @examples
  //' \dontrun{
  //' data <- data.frame(id_curve = rep(1:10, each = 10), tobs = rep(1:10, 10), X = rnorm(100))
  //' optimal_bw <- get_nw_optimal_bw_cpp(data)
  //' }
  //' @export
  //'
  // [[Rcpp::export]]
  double get_nw_optimal_bw_cpp(const Rcpp::DataFrame data,
                               const Nullable<arma::vec> bw_grid = R_NilValue,
                               const Nullable<int> nsubset = R_NilValue,
                               const std::string kernel_name = "epanechnikov") {
    // Check if kernel_name is one of the supported kernels
    std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                  "tricube", "triangular", "uniform"};
    if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
      stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
    }
    return nw_optimal_bw_impl(core::make_curves(data), bw_grid, nsubset, core::parse_kernel(kernel_name));
  }




// :::::::::::::: Design density (Parzen-Rosenblatt) ::::::::::::
// Leave-one-out estimate at the observation points,
// sum_{j != i} K((x_i - x_j) / h) / ((n - 1) h), on the sorted points (kernel
// window), returned in the order of x.
arma::vec density_loo_impl(const arma::vec& x, double h, core::Kern kern) {
  const arma::uword n = x.n_elem;
  std::vector<double> xs;
  std::vector<arma::uword> perm;
  core::sort_design(x, xs, perm);
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i) {
    arma::uword lo, hi;
    core::window(xs, x(i), h, lo, hi);
    double sum = 0;
    for (arma::uword q = lo; q < hi; ++q) if (perm[q] != i) sum += core::kval(kern, (x(i) - xs[q]) / h);
    out(i) = sum / ((n - 1) * h);
  }
  return out;
}

//' Leave-one-out design-density estimate (C++ core)
//'
//' @param x Observation times of one curve.
//' @param h Bandwidth.
//' @param kernel_name Kernel name.
//' @return The leave-one-out estimate at each point of `x`.
//' @keywords internal
// [[Rcpp::export]]
arma::vec density_loo_cpp(const arma::vec x, const double h, const std::string kernel_name) {
  return density_loo_impl(x, h, core::parse_kernel(kernel_name));
}

//' Least-squares cross-validation scores of the design density (C++ core)
//'
//' For each bandwidth, the trapezoidal integral of the squared estimate over
//' `zgrid` minus twice the mean leave-one-out estimate.
//'
//' @param x Observation times of one curve.
//' @param zgrid Integration grid.
//' @param bw_grid Candidate bandwidths.
//' @param kernel_name Kernel name.
//' @return The LSCV score per bandwidth.
//' @keywords internal
// [[Rcpp::export]]
arma::vec density_lscv_cpp(const arma::vec x, const arma::vec zgrid, const arma::vec bw_grid,
                           const std::string kernel_name) {
  core::Kern kern = core::parse_kernel(kernel_name);
  const arma::uword n = x.n_elem, nz = zgrid.n_elem;
  std::vector<double> xs;
  std::vector<arma::uword> perm;
  core::sort_design(x, xs, perm);
  arma::vec score(bw_grid.n_elem);
  arma::vec f(nz);
  for (arma::uword b = 0; b < bw_grid.n_elem; ++b) {
    const double h = bw_grid(b);
    for (arma::uword i = 0; i < nz; ++i) {
      arma::uword lo, hi;
      core::window(xs, zgrid(i), h, lo, hi);
      double sum = 0;
      for (arma::uword q = lo; q < hi; ++q) sum += core::kval(kern, (zgrid(i) - xs[q]) / h);
      f(i) = sum / (n * h);
    }
    double integral = 0;
    for (arma::uword i = 0; i + 1 < nz; ++i)
      integral += (zgrid(i + 1) - zgrid(i)) * (f(i + 1) * f(i + 1) + f(i) * f(i)) / 2;
    score(b) = integral - 2 * arma::mean(density_loo_impl(x, h, kern));
  }
  return score;
}
