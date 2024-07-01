#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
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
 //' @export
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
 //'  //' @param h Numeric (positive vector or scalar). The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
 //' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
 //' If \code{h} is a scalar, then all curves will be smoothed with the same bandwidth.
 //' Otherwise, if \code{h} is a vector, its length must be equal to the number of curves in \code{data}
 //' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
 //' @param h h Numeric (positive vector or scalar). The bandwidth parameter for the kernel estimator.
 //' If \code{h} is a scalar, then all pointd will be estimated with the same bandwidth.
 //' Otherwise, if \code{h} is a vector, its length must be equal to the number of points in \code{t}
 //' and each element of the vector must correspond to a point given in the same order as in \code{t}.
 //' @param kernel_name A string specifying the name of the kernel function to use. Default is "epanechnikov".
 //'   Supported kernels: "epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform".
 //'
 //' @return A vector of estimated response values for the new predictor values.
 //'
 //' @details This function estimates the Nadaraya-Watson kernel estimator for regression,
 //' using a specified kernel function. It returns only only vector.
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
   if (!arma::is_finite(y) || y.size() == 0 ||
       !arma::is_finite(t) || t.size() == 0 ||
       y.size() != t.size()) {
     stop("The vectors y and t must be non-empty numeric vectors of the same length.");
   }
   int m = tnew.size();
   int n = y.size();

   // Check if h is non-negative
   // Control the presmoothing bandwidth
   arma::vec h_to_use;
   if (h.size() != m) {
     if (h.size() != 1) {
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

   // Select the kernel function
   std::function<arma::vec(const arma::vec)> kernel_func = select_kernel(kernel_name);

   // Calculate weight matrix A
   arma::mat A = arma::zeros(m, n);

   for (int i = 0; i < m; i++) {
     A.row(i) = kernel_func((t - tnew(i)) / h_to_use(i)).t();

     // Normalize each row
     A.row(i) /= sum(A.row(i));
   }

   // Compute yhat
   arma::vec yhat = A * y;

   return yhat;
 }


 /***R
 n <- 100
 Tn <- runif(n)
 Yn <- sin(2 * pi * Tn) + rnorm(n = n, mean = 0, sd = 0.25)
 estimate_nw_cpp(y = Yn, t = Tn, tnew = seq(0.1, 0.9, len = 100), h = 0.05)
 */

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
    if (!arma::is_finite(y) || y.size() == 0 ||
        !arma::is_finite(t) || t.size() == 0 ||
        y.size() != t.size()) {
      stop("The vectors y and t must be non-empty numeric vectors of the same length.");
    }

    // Check if bw_grid is a numeric vector and not empty
    if (!arma::is_finite(bw_grid) || bw_grid.size() == 0) {
      stop("bw_grid must be a non-empty numeric vector.");
    }

    // Check if kernel_name is one of the supported kernels
    std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                  "tricube", "triangular", "uniform"};
    if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
      stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
    }

    // Select the kernel function
    std::function<arma::vec(const arma::vec)> kernel_func = select_kernel(kernel_name);

    // Compute the number of observations and the size of the bandwidth grid
    int n = t.size();
    int bw_size = bw_grid.size();

    // Vector to store cross-validation errors for each bandwidth value
    arma::vec cv_error(bw_size);

    // Loop over each bandwidth value in the grid
    for (int i = 0; i < bw_size; i++) {
      // Current bandwidth value
      arma::vec hi(1, fill::value(bw_grid(i)));

      // Estimate the response values using Nadaraya-Watson kernel regression
      arma::vec yhat = estimate_nw_cpp(y, t, t, hi, kernel_name);

      // Calculate the weight matrix A
      arma::mat wmat(n, n);
      for (int j = 0; j < n; j++) {
        wmat.col(j) = kernel_func((t - t(j)) / hi(0));
      }

      // Compute the residuals divided by the effective weights
      arma::vec metric = (y - yhat) / (1 - kernel_func(arma::zeros(n)) / arma::sum(wmat, 1));

      // Remove NaNs from the metric
      metric = metric.elem(find_finite(metric));

      // Calculate the mean squared error for cross-validation
      metric = arma::pow(metric, 2);
      cv_error(i) = metric.is_empty() ? arma::datum::inf : mean(metric);
    }

    // Find the bandwidth value with the minimum cross-validation error
    arma::uword idx_min = arma::index_min(cv_error);

    double hcv = bw_grid(idx_min);

    return hcv;
  }

 //' Estimate Optimal Bandwidth for Nadaraya-Watson Estimator
  //'
  //' Estimates the median of the best bandwidths obtained by Nadaraya-Watson cross-validation on a subset of curves.
  //'
  //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X". Typically, the output of the function \link{.format_data}.
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
    // Prepare the data
    arma::mat data_mat(data.nrows(), 3);
    data_mat.col(0) = as<arma::vec>(data["id_curve"]);
    data_mat.col(1) = as<arma::vec>(data["tobs"]);
    data_mat.col(2) = as<arma::vec>(data["X"]);
    arma::vec unique_id_curve = arma::unique(data_mat.col(0));
    int n_curve = unique_id_curve.n_elem;

    // Set the number of subset of curves to be used to determine the best N-W bandwidth
    double nsubset_val = 0;
    if (nsubset.isNotNull()) {
      nsubset_val = as<int>(nsubset);
    } else {
      nsubset_val = std::min(70.0, std::floor(n_curve * 0.5));
    }
    if (nsubset.isNotNull() && (nsubset_val > n_curve)) {
      stop("If 'nsubset' is not NULL, it must be a positive integer less than or equal to the number of curves.");
    }

    // Set the bandwidth grid
    arma::vec bw_grid_used;
    if (bw_grid.isNull()) {

      //// Estimate lambda
      double lambdahat = arma::mean(hist(data_mat.col(0), unique_id_curve));
      double b0 = 1 / lambdahat;
      double bK = std::pow(lambdahat, -1 / 3.0);
      bw_grid_used = arma::logspace(log10(b0), log10(bK), 15);
    } else {
      bw_grid_used = as<arma::vec>(bw_grid);
    }

    // Init to add optimal bandwidth
    arma::vec res_opt_bw(nsubset_val);

    // Estimate the best bandwidth
    for (int i = 0; i < nsubset_val; ++i) {
      arma::uvec indices_cur = arma::find(data_mat.col(0) == unique_id_curve(n_curve - 1 - i));
      double hbest = estimate_nw_bw_cpp(data_mat(indices_cur, arma::uvec({2})), data_mat(indices_cur, arma::uvec({1})), bw_grid_used, kernel_name);
      res_opt_bw(i) = hbest;
    }

    return arma::median(res_opt_bw);
  }



