#' Biweight kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [triweight()], [tricube()], [epanechnikov()], [triangular()], and [uniform()].
#'
biweight <- function(u){
  ifelse(abs(u) <= 1, (15/16) * (1 - u^2)^2, 0)
}

#' Triweight kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [tricube()], [epanechnikov()], [triangular()], and [uniform()].
#'
triweight <- function(u){
  ifelse(abs(u) <= 1, (35/32) * (1 - u^2)^3, 0)
}

#' Tricube kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [triweight()], [epanechnikov()], [triangular()], and [uniform()].
#'
tricube <- function(u){
  ifelse(abs(u) <= 1, (70/81) * (1 - abs(u)^3)^3, 0)
}

#' Epanechnikov kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [triweight()], [tricube()], [triangular()], and [uniform()].
#'
epanechnikov <- function(u){
  ifelse(abs(u) <= 1, (3/4) * (1 - u ** 2), 0)
}

#' Triangular kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [triweight()], [tricube()], [epanechnikov()], and [uniform()].
#'
triangular <- function(u){
  ifelse(abs(u) <= 1, (1 - abs(u)), 0)
}

#' Uniform kernel function
#'
#' @param u \code{numeric}. Scalar or vector of numeric values at which to evaluate the function.
#'
#' @return A scalar or vector of \code{numeric}.
#' @export
#' @seealso [biweight()], [triweight()], [tricube()], [epanechnikov()], and [triangular()].
#'
uniform <- function(u){
  ifelse(abs(u) <= 1, 1/2, 0)
}

#' Nadaraya-Watson estimator
#'
#' @param y \code{vector (numeric)}. A numeric vector containing the observed values of the independent variable corresponding to the observation points \code{t}.
#' @param t \code{vector (numeric)}. A numeric vector containing the observed values of the dependent variable.
#' @param tnew \code{vector (numeric)}. New \code{t} values at which we want to estimate the regression function.
#' @param h \code{numeric (positive)}. The bandwidth parameter such that \code{h > (2 * length(x))}.
#' Default \code{h = NULL} and such it will be computed automatically.
#' @param smooth_ker \code{function}. The kernel function of the estimator.
#'
#' @return A \code{data.table} containing
#' \itemize{
#'    \item{h :}{ The bandwidth used to estimate the regression function.}
#'    \item{inKernelSupp :}{ For each \code{t_i} in \code{tnew}, it is the number of  \code{t} between \code{t_i - h } and \code{t_i + h }.}
#'    \item{tnew :}{ The vector \code{new}.}
#'    \item{yhat :}{ The regression function's vector of estimates at \code{tnew}.}
#' }
#'
#' @importFrom methods is
#' @importFrom data.table data.table
#'
#' @export
#'
#' @seealso [estimate_nw_bw()], [epanechnikov()], [biweight()], [triweight()], [tricube()], [uniform()], etc.
#'
#' @examples
#' \dontrun{
#' # The model
#' ## Let
#' m <- function(t) 4 * sin(1.5 * pi * t)
#'
#' ## Observation points
#' t <- runif(n = 200, min = 0, max = 1)
#' t <- sort(t)
#'
#' ## Measure error
#' e <- rnorm(n = 200, mean = 0, sd = 0.2)
#'
#' ## Regression model
#' y <- m(t) + e
#'
#' plot(x = t, y = y, main = "Observed points and true regression function")
#' lines(x = t, y = m(t), type = "l", col = "red")
#'
#' ## Estimate the best bandwidth
#' h_grid <- seq(1 / (2 * length(t)), length(t) ** (- 1/3), len = 100)
#' hbest <- estimate_nw_bw(y = y, t = t,
#'                         h_grid = h_grid,
#'                         smooth_ker = epanechnikov)
#'
#' ## Estimate the regression function
#' dt_nw <- estimate_nw(y = y, t = t,
#'                      tnew = seq(0.01, 0.99, len = 100),
#'                      h = hbest, smooth_ker = epanechnikov)
#'
#' plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
#'      main = "Estimated and true regression function.")
#' lines(x = dt_nw[, tnew], y = m(dt_nw[, yhat]), type = "l", col = "red")
#' legend(x = 0.64, y = 4.1, fill = c("blue", "red"),legend = c("Estimated m", "True m"))
#'
#' }
#'
estimate_nw <- function(y, t, tnew, h = NULL, smooth_ker = epanechnikov){
  if (! is.numeric(y) | ! is.numeric(t) | ! is.numeric(tnew) |! is.numeric(h))
    stop("The arguments 'y', 't', 'tnew', 'h' must be numeric.")
  if (length(t) != length(y) & length(t) < 2)
    stop("The arguments 'y' and 't' must have a length of at least 2 and must be of the same length.")
  if (! is.null(h) & length(h) > 1)
    stop("The bandwidth 'h' must be either NULL or saclar.")
  if (! methods::is(object = smooth_ker, class2 = "function"))
    stop("'smooth_ker' must be a function.")

  m <- length(tnew)
  n <- length(y)
  A <- matrix(0, nrow = m, ncol = n)
  if (is.null(h)) {
    hcv <- estimate_nw_bw(y = y, t = t,
                          h_grid = seq(2 / n, n ** (-1/3), len = 100),
                          smooth_ker = smooth_ker)
  }
  h <- ifelse(is.null(h), hcv, h)
  A <- outer(tnew, t, function(u, v) smooth_ker((u - v) / h))

  ## Get the number of points used to estimate y for each xout
  inKernelSupp <- unlist(lapply(tnew, function(tnewi, t, h){
    sum(abs(tnewi - t) <= h)
  }, t = t, h = h))

  A <- apply(A, 1, function(x) x / sum(x))
  A <- t(A)
  yhat <- A %*% y
  dt <- data.table::data.table("h" = h, "inKernelSupp" = inKernelSupp, "tnew" = tnew, "yhat" = c(yhat))
  return(dt)
}

#' Nadaraya-Watson Bandwidth Selection using cross validation.
#'
#' @param y \code{vector (numeric)}. A numeric vector containing the observed values of the independent variable corresponding to the observation points \code{t}.
#' @param t \code{vector (numeric)}. A numeric vector containing the observed values of the dependent variable.
#' @param h_grid \code{vector (numeric)}. A grid of bandwidth to test.
#' @param smooth_ker \code{function}. The kernel function of the estimator.
#'
#' @return A \code{numeric} value corresponding to the best bandwidth.
#' @export
#' @seealso [estimate_nw()]
#'
#' @examples
#' \dontrun{
#' # The model
#' ## Let
#' m <- function(t) 4 * sin(1.5 * pi * t)
#'
#' ## Observation points
#' t <- runif(n = 200, min = 0, max = 1)
#' t <- sort(t)
#'
#' ## Measure error
#' e <- rnorm(n = 200, mean = 0, sd = 0.2)
#'
#' ## Regression model
#' y <- m(t) + e
#'
#' plot(x = t, y = y, main = "Observed points and true regression function")
#' lines(x = t, y = m(t), type = "l", col = "red")
#'
#' ## Estimate the best bandwidth
#' h_grid <- seq(1 / (2 * length(t)), length(t) ** (- 1/3), len = 100)
#' hbest <- estimate_nw_bw(y = y, t = t,
#'                         h_grid = h_grid,
#'                         smooth_ker = epanechnikov)
#'
#' ## Estimate the regression function
#' dt_nw <- estimate_nw(y = y, t = t,
#'                      tnew = seq(0.01, 0.99, len = 100),
#'                      h = hbest, smooth_ker = epanechnikov)
#'
#' plot(x = dt_nw[, tnew], y = dt_nw[, yhat], type = "l", col = "blue",
#'      main = "Estimated and true regression function.")
#' lines(x = dt_nw[, tnew], y = m(dt_nw[, tnew]), type = "l", col = "red")
#' legend(x = 0.64, y = 4.1, fill = c("blue", "red"),legend = c("Estimated m", "True m"))
#'
#' }
#'
estimate_nw_bw <- function(y, t, h_grid = seq(1 / (2 * length(t)), length(t) ** (- 1/3), len = 100),
                           smooth_ker = epanechnikov) {
  if (! is.numeric(y) | ! is.numeric(t) |! is.numeric(h_grid))
    stop("The arguments 'y', 't', 'tnew', 'h' must be numeric.")
  if (length(t) != length(y) & length(t) < 2)
    stop("The arguments 'y' and 't' must have a length of at least 2 and must be of the same length.")
  if (is.null(h_grid))
    stop("The bandwidth grid 'h_grid' must be a scalar or vector of numeric.")
  if (! methods::is(object = smooth_ker, class2 = "function"))
    stop("'smooth_ker' must be a function.")

  cv_error <- sapply(h_grid, function(hi, y, t, K){
    yhat <- estimate_nw(y = y, t = t, h = hi, tnew = t, smooth_ker = K)$yhat
    wmat <- outer(X = t, Y = t, function(u, v) K((u-v) / hi))
    metric <- (y - yhat) / (1 - K(0) / rowSums(wmat))

    # If there is only one value in the kernel support, it return NaN.
    error_hi <- mean(metric[! is.nan(metric)] ** 2)
  }, y = y, t = t, K = smooth_ker)

  # If cv_error is NaN, do take it into account
  if (any(is.nan(cv_error))) {
    h_grid <- h_grid[-which(is.nan(cv_error))]
    cv_error <- cv_error[! is.nan(cv_error)]
    hcv <- h_grid[which.min(cv_error)]
  } else {
    hcv <- h_grid[which.min(cv_error)]
  }
  return(hcv)
}


#' Estimate the the standard deviation of the observation error
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the standard deviation of the error.
#'
#' @return A data.table with two columns: \code{t} and \code{sig} corresponding to the estimated standard deviation.
#' @export
#'
#' @import data.table
#'
estimate_sigma <- function(data, idcol = NULL, tcol = "tobs", ycol = "X", t = 1/2) {
  # Format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)

  dt_sig <- data.table::rbindlist(lapply(data[, unique(id_curve)], function(idx, data, t) {
    # Compute get the index the two T_{n,k} closest to time for each X_n
    ## Get the rank i_t and j_t
    data_idx <- data[id_curve == idx]
    dif_mat <- outer(X = t, Y = data_idx[, tobs], function(ti, tobsi) abs(ti - tobsi))
    rank_mat <- apply(X = dif_mat, MARGIN = 1, function(r) rank(r))

    ## Compute [Y_{n,i_t} - Y_{n,j_t}]^2 / 2 for each curve X_n and each t
    Zn <- apply(X = rank_mat, MARGIN = 2, function(c, dt){
      first <- which(c == 1)
      second <- which(c == 2)
      Y1 <- dt[first][, X]
      Y2 <- dt[second][, X]
      Zn <- 1 / 2 * (Y1 - Y2) ** 2
    }, dt = data_idx)
    dt_res <- data.table::data.table("t" = t, "Zn" = Zn)
    rm(data_idx, Zn)
    return(dt_res)
  }, data = data, t = t))
  # Estimate sigma
  dt_sig[, sig := sqrt(mean(Zn)), by = t]
  dt_sig <- unique(dt_sig[, list(t, sig)])

  return(dt_sig)
}
