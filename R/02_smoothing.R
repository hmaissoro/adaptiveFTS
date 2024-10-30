#' Nadaraya-Watson Kernel Estimator
#'
#' Estimates the Nadaraya-Watson regression function using a specified kernel function.
#'
#' @param y \code{vector (numeric)}. A numeric vector containing the observed values of the dependent variable corresponding to the observation points \code{t}.
#' @param t \code{vector (numeric)}. A numeric vector containing the observed values of the independent variable.
#' @param tnew \code{vector (numeric)}. New \code{t} values at which to estimate the regression function.
#' @param h \code{numeric (positive)}. The bandwidth parameter, where \code{h > (2 * length(x))}.
#' Default is \code{h = NULL}, in which case it will be computed automatically.
#' @param kernel_name \code{string}. A string specifying the name of the kernel function to use. The default is "epanechnikov".
#' Supported kernels: "epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform".
#'
#' @return A \code{data.table} with the following columns:
#' \itemize{
#'    \item{\code{h} :}{ The bandwidth used to estimate the regression function.}
#'    \item{\code{tnew} :}{ The vector of new \code{t} values at which the regression function is estimated.}
#'    \item{\code{yhat} :}{ A vector of estimated values of the regression function at \code{tnew}.}
#' }
#'
#' @importFrom methods is
#' @importFrom data.table data.table
#'
#' @export
#'
#' @seealso \code{\link{estimate_nw_bw}}
#'
#' @examples
#' \dontrun{
#' # Define the true regression function
#' m <- function(t) 4 * sin(1.5 * pi * t)
#'
#' # Generate observation points
#' t <- sort(runif(n = 200, min = 0, max = 1))
#'
#' # Generate measurement errors
#' e <- rnorm(n = 200, mean = 0, sd = 0.2)
#'
#' # Observed values of the regression model
#' y <- m(t) + e
#'
#' # Plot observed points and true regression function
#' plot(x = t, y = y, main = "Observed points and true regression function")
#' lines(x = t, y = m(t), col = "red")
#'
#' # Estimate optimal bandwidth using cross-validation
#' bw_grid <- seq(1 / (2 * length(t)), length(t) ** (-1/3), length.out = 100)
#' hbest <- estimate_nw_bw(y = y, t = t, bw_grid = bw_grid, kernel_name = "epanechnikov")
#'
#' # Estimate regression function with the Nadaraya-Watson estimator
#' dt_nw <- estimate_nw(y = y, t = t, tnew = seq(0.01, 0.99, length.out = 100), h = hbest, kernel_name = "epanechnikov")
#'
#' # Plot estimated and true regression functions
#' plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
#'      main = "Estimated and true regression function")
#' lines(x = dt_nw[, tnew], y = m(dt_nw[, tnew]), col = "red")
#' legend(x = 0.64, y = 4.1, fill = c("blue", "red"), legend = c("Estimated m", "True m"))
#' }
#'
estimate_nw <- function(y, t, tnew, h = NULL, kernel_name = "epanechnikov"){
  if (! is.numeric(y) | ! is.numeric(t) | ! is.numeric(tnew) |! is.numeric(h))
    stop("The arguments 'y', 't', 'tnew', 'h' must be numeric.")
  if (length(t) != length(y) & length(t) < 2)
    stop("The arguments 'y' and 't' must have a length of at least 2 and must be of the same length.")
  if (! is.null(h) & length(h) > 1)
    stop("The bandwidth 'h' must be either NULL or saclar.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  if (is.null(h)) {
    hcv <- estimate_nw_bw(y = y, t = t,
                          bw_grid = seq(2 / n, n ** (-1/3), len = 100),
                          kernel_name = kernel_name)
  }
  h <- ifelse(is.null(h), hcv, h)

  # Estimate Nadaraya-Watson estimator using estimate_nw_cpp function
  res <- estimate_nw_cpp(y = y, t = t, tnew = tnew, h = h, kernel_name = kernel_name)

  # convert to data.table
  dt_res <- data.table::data.table("h" = h, "tnew" = tnew, "yhat" = c(res))
  return(dt_res)
}

#' Nadaraya-Watson Bandwidth Selection using Cross-Validation
#'
#' Selects the optimal bandwidth for Nadaraya-Watson kernel regression using cross-validation.
#'
#' @param y \code{vector (numeric)}. A numeric vector containing the observed values of the dependent variable.
#' @param t \code{vector (numeric)}. A numeric vector containing the observed values of the independent variable.
#' @param bw_grid \code{vector (numeric)}. A grid of bandwidth values to test. Default is \code{bw_grid = NULL}, in which case an exponential grid based on the length of \code{t} will be used.
#' @param kernel_name \code{string}. A string specifying the name of the kernel function to use. The default is "epanechnikov".
#' Supported kernels: "epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform".
#'
#' @return A \code{numeric} value representing the optimal bandwidth that minimizes the cross-validation error.
#'
#' @details This function computes the optimal bandwidth for the Nadaraya-Watson kernel regression estimator by performing cross-validation over a set of candidate bandwidths defined by \code{bw_grid}. The function returns the bandwidth that minimizes the cross-validation error.
#'
#' @export
#' @seealso \code{\link{estimate_nw}}
#'
#' @examples
#' \dontrun{
#' # Define the true regression function
#' m <- function(t) 4 * sin(1.5 * pi * t)
#'
#' # Generate observation points
#' t <- sort(runif(n = 200, min = 0, max = 1))
#'
#' # Generate measurement errors
#' e <- rnorm(n = 200, mean = 0, sd = 0.2)
#'
#' # Observed values of the regression model
#' y <- m(t) + e
#'
#' # Plot observed points and true regression function
#' plot(x = t, y = y, main = "Observed points and true regression function")
#' lines(x = t, y = m(t), col = "red")
#'
#' # Define a grid of candidate bandwidths
#' bw_grid <- seq(1 / (2 * length(t)), length(t)^(-1/3), length.out = 100)
#'
#' # Estimate the best bandwidth using cross-validation
#' hbest <- estimate_nw_bw(y = y, t = t, bw_grid = bw_grid, kernel_name = "epanechnikov")
#'
#' # Estimate the regression function using the selected bandwidth
#' dt_nw <- estimate_nw(y = y, t = t, tnew = seq(0.01, 0.99, length.out = 100), h = hbest, kernel_name = "epanechnikov")
#'
#' # Plot estimated and true regression functions
#' plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
#'      main = "Estimated and true regression function")
#' lines(x = dt_nw[, tnew], y = m(dt_nw[, tnew]), col = "red")
#' legend(x = 0.64, y = 4.1, fill = c("blue", "red"), legend = c("Estimated m", "True m"))
#' }
#'
estimate_nw_bw <- function(y, t, bw_grid = NULL, kernel_name = "epanechnikov") {
  if (! is.numeric(y) | ! is.numeric(t) |! is.numeric(bw_grid))
    stop("The arguments 'y', 't', 'tnew', 'h' must be numeric.")
  if (length(t) != length(y) & length(t) < 2)
    stop("The arguments 'y' and 't' must have a length of at least 2 and must be of the same length.")
  if (! is.numeric(bw_grid) & ! is.numeric(bw_grid))
    stop("If the bandwidth grid 'bw_grid' is not NULL, so it must be a scalar or vector of numeric.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  if (is.null(bw_grid)) {
    K <- 100
    b0 <- 1 / length(t)
    bK <- length(t) ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(b0, bK, a, K) ; gc()
  }

  # Use C++ function
  hcv <- estimate_nw_bw_cpp(y = y, t = t, bw_grid = bw_grid, kernel_name = kernel_name)

  return(hcv)
}

#' Estimate Optimal Bandwidth for Nadaraya-Watson Estimator on a Subset of Curves
#'
#' This function estimates the median optimal bandwidth for the Nadaraya-Watson kernel estimator using cross-validation on a subset of curves.
#'
#' @inheritParams format_data
#' @param bw_grid \code{vector (numeric)}. A grid of candidate bandwidth values for cross-validation. Default is \code{bw_grid = NULL}, which sets an exponential grid based on the average number of observation points per curve.
#' @param nsubset \code{integer (positive)}. The number of curves to randomly and uniformly select for bandwidth optimization. Default is \code{nsubset = NULL}, in which case an optimal bandwidth is calculated for each curve.
#' @param kernel_name \code{string}. A string specifying the name of the kernel function to use, with "epanechnikov" as the default.
#' Supported kernels: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{numeric} scalar representing the estimated optimal bandwidth as the median of the best bandwidths from the subset of curves.
#'
#' @details This function performs cross-validation to determine the optimal bandwidth for each curve in a specified subset. It returns the median of these best bandwidths as the final estimate, providing a representative bandwidth that can generalize across curves.
#'
#' @export
#'
#' @importFrom methods is
#' @import data.table
#'
#' @examples
#' \dontrun{
#' # Load the dataset
#' data(data_far)
#'
#' # Estimate the optimal bandwidth on a subset of 30 curves
#' hbest <- get_nw_optimal_bw(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'                            bw_grid = NULL, nsubset = 30, kernel_name = "epanechnikov")
#' # Display the optimal bandwidth
#' hbest
#' }
#'
#'
get_nw_optimal_bw <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                              bw_grid = NULL, nsubset = NULL, kernel_name = "epanechnikov"){
  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Control parameters
  if ((!is.null(bw_grid)) & (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1)))
    stop("If'bw_grid' is not NULL, then it must be a vector of positive values between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  if (!is.null(nsubset))
    if (any(nsubset < 0) | (length(nsubset) > 1) | any(nsubset - floor(nsubset) > 0) | any(N <= nsubset))
      stop("If 'nsubset' is not NULL, then if must be a positive integer lower than the number of curves.")

  # Define the set of curve
  if (! is.null(nsubset)) {
    sample_curves <- sample(x = 1:N, size = nsubset)
  } else {
    sample_curves <- 1:N
  }

  # Define the grid
  if (is.null(bw_grid)) {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    K <- 100
    b0 <- 2 / lambdahat
    bK <- lambdahat ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))

    rm(K, b0, bK, a) ; gc()
  }

  # Estimate the bandwidths using the C++ function
  hbest <- get_nw_optimal_bw_cpp(data = data, bw_grid = bw_grid, nsubset = nsubset, kernel_name = )

  return(hbest)
}
