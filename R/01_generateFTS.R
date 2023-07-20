#' Arctan Hurst function
#'
#' Arctan Hurst function that can be used to generate multifractional Brownian motion (mfBm).
#' See the following paper https://doi.org/10.3390/fractalfract6020074.
#'
#' @param t \code{vector (float)}. Points between 0 and 1 at which to evaluate the function.
#'
#' @return A \code{vector (float)} corresponding to the value of the function evaluated at \code{t}.
#'
#' @export
#'
#' @seealso \code{\link{hurst_linear}}, \code{\link{hurst_logistic}}.
#'
#' @examples
#'
#' t0 <- seq(0.2, 0.8, len = 10)
#' htan <- hurst_arctan(t = t0)
#' plot(x = t0, y = htan, type = "b", col = "red")
#'
#'
hurst_arctan <- function(t = seq(0.2, 0.8, len = 10)){
  hval <- atan(t) / pi + 1/2
  return(hval)
}

#' Linear Hurst function
#'
#' Linear Hurst function that can be used to generate multifractional Brownian motion (mfBm).
#' See the following paper https://doi.org/10.3390/fractalfract6020074.
#'
#' @param t \code{vector (float)}. Points between 0 and 1 at which to evaluate the function.
#' @param h_left \code{Float}. A scalar value in the interval between 0 and 1 indicating the minimum of the function.
#' @param h_right \code{Float}. A scalar value in the interval between 0 and 1 indicating the maximum of the function.
#'
#' @return A \code{vector (float)} corresponding to the value of the function evaluated at \code{t}.
#'
#' @export
#'
#' @seealso \code{\link{hurst_arctan}}, \code{\link{hurst_logistic}}.
#'
#' @examples
#' t0 <- seq(0.2, 0.8, len = 10)
#' hlinear <- hurst_linear(t = t0)
#' plot(x = t0, y = hlinear, type = "b", col = "red")
#'
#'
hurst_linear <- function(t = seq(0.2, 0.8, len = 10), h_left = 0.2, h_right = 0.8) {
  t1 <- 1
  t0 <- 0
  a <- (h_right - h_left) / (t1 - t0)
  b <- h_right - a * t1
  hval <- pmin(a * t + b, 1)
  return(hval)
}

#' Logistic Hurst function
#'
#' Logistic Hurst function that can be used to generate multifractional Brownian motion (mfBm).
#' See the following paper https://doi.org/10.3390/fractalfract6020074.
#'
#' @param t \code{vector (float)}. Points between 0 and 1 at which to evaluate the function.
#' @param h_left \code{Float}. A scalar value in the interval between 0 and 1 indicating the minimum of the function.
#' @param h_right \code{Float}. A scalar value in the interval between 0 and 1 indicating the maximum of the function.
#' @param slope \code{Float (positive)}. A scalar positive value corresponding to the slope of the logistic function.
#' @param change_point_position \code{Float}. A scalar value in the interval between 0 and 1 corresponding ti the change point position.
#'
#' @return A \code{vector (float)} corresponding to the value of the function evaluated at \code{t}.
#'
#' @export
#'
#' @seealso \code{\link{hurst_arctan}}, \code{\link{hurst_linear}}.
#'
#' @examples
#' t0 <- seq(0.2, 0.8, len = 10)
#' hlogistic <- hurst_logistic(t = t0, h_left = 0.2,
#'                             h_right = 0.8, slope = 30,
#'                             change_point_position = 0.5)
#' plot(x = t0, y = hlogistic, type = "b", col = "red")
#'
#'
hurst_logistic <- function(t, h_left = 0.2, h_right = 0.8, slope = 30,
                           change_point_position = 0.5) {
  u <- (t - change_point_position) / (1 - 0)
  hval <- (h_right - h_left) / (1 + exp(- slope * u)) + h_left
  return(hval)
}

#' Constant D(x,y) function
#'
#' See the following paper https://doi.org/10.3390/fractalfract6020074.
#'
#' @param x \code{Float (positive)}. First argument of the function.
#' @param y \code{Float (positive)}. Second argument of the function.
#'
#' @return A positive \code{Float} corresponding to the value the function evaluate at (x,y).
#' @export
#'
.constant_d <- function(x, y) {
  a <- gamma(2 * x + 1) * gamma(2 * y + 1) * sin(pi * x) * sin(pi * y)
  b <- 2 * gamma(x + y + 1) * sin(pi * (x + y) / 2)
  val <- sqrt(a) / b
  return(val)
}

#' Title
#'
#' @param t \code{vector (float)}. Points between 0 and 1 at which to compute the covariance function.
#' @param hurst_fun \code{function}. Hurst function. It can be \code{\link{hurst_arctan}}, \code{\link{hurst_linear}}, \code{\link{hurst_logistic}}.
#' @param ... Hurst function additional arguments.
#'
#' @return a \code{matrix} of \code{t} x \code{t} covariance.
#' @export
#'
.covariance_mfBm <- function(t = seq(0.2, 0.8, len = 10), hurst_fun = hurst_logistic, ...) {
  tmp <- expand.grid(u = t, v = t)
  u <- tmp$u
  v <- tmp$v
  hu <- hurst_fun(u, ...)
  hv <- hurst_fun(v, ...)
  hh <- hu + hv
  values <- .constant_d(hu, hv) *
    (u**hh + v**hh - abs(u - v) ** hh)
  mat <- matrix(values, ncol = length(t))
  return(mat)
}


#' Draw a multifractional Brownian motion sample path.
#'
#' @param t \code{vector (float)}. Grid of points between 0 and 1 where we want to generate the sample path.
#' @param hurst_fun \code{function}. Hurst function. It can be \code{\link{hurst_arctan}}, \code{\link{hurst_linear}}, \code{\link{hurst_logistic}}.
#' @param L \code{float (positive)}. Hölder constant.
#' @param tied \code{boolean}. If \code{TRUE}, the sample path is tied-down.
#' @param ... Hurst function additional arguments.
#'
#' @return A \code{data.table} containing 2 column : \code{t} and \code{mfBm}, the sample path.
#'
#' @importFrom MASS mvrnorm
#' @importFrom data.table data.table
#'
#' @export
#'
#' @examples
#'
#' t0 <- seq(0.2, 0.8, len = 20)
#' dt_mfBm <- simulate_mfBm(t = t0, hurst_fun = hurst_logistic, L = 1, tied = TRUE)
#' plot(x = dt_mfBm$t, y = dt_mfBm$mfBm, type = "l", col = "red")
#'
simulate_mfBm <- function(t = seq(0.2, 0.8, len = 50), hurst_fun = hurst_logistic, L = 1, tied = TRUE, ...) {
  t <- sort(t)
  cov_mat <- .covariance_mfBm(t = t, hurst_fun = hurst_fun, ...)
  out <- MASS::mvrnorm(1,
                       mu = rep(0, ncol(cov_mat)),
                       Sigma = L * cov_mat)
  mfBm_path <- out - tied * t * out[length(out)]
  dt <- data.table::data.table("t" = t, mfBm = mfBm_path)
  return(dt)
}


#' Draw a fractional Brownian motion sample path.
#'
#' @param t \code{vector (float)}. Grid of points between 0 and 1 where we want to generate the sample path.
#' @param hurst \code{float (positive)}. The Hurst exponent scalar value between 0 and 1.
#' @param L \code{float (positive)}. Hölder constant.
#' @param tied \code{boolean}. If \code{TRUE}, the sample path is tied-down.
#'
#' @return A \code{data.table} containing 2 column : \code{t} and \code{mfBm}, the sample path.
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom data.table data.table
#'
#' @examples
#'
#' t0 <- seq(0.2, 0.8, len = 20)
#' dt_fBm <- simulate_fBm(t = t0, hurst = 0.6, L = 1, tied = TRUE)
#' plot(x = dt_fBm$t, y = dt_fBm$fBm, type = "l", col = "red")
#'
simulate_fBm <- function(t = seq(0.2, 0.8, len = 20), hurst = 0.6, L = 1, tied = TRUE) {
  tmp <- expand.grid(u = t, v = t)
  u <- tmp$u
  v <- tmp$v
  values <- (1 / 2) *
    (u ** hurst + v ** hurst - abs(u - v) ** hurst)
  cov_mat <- matrix(values, ncol = length(t))

  out <- MASS::mvrnorm(1,
                       mu = rep(0, ncol(cov_mat)),
                       Sigma = L * cov_mat)
  fBm_path <- out - tied * t * out[length(out)]
  dt <- data.table::data.table("t" = t, "fBm" = fBm_path)
  return(dt)
}
