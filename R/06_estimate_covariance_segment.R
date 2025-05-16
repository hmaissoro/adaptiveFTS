#' Estimate the Risk of the Covariance Segment Function
#'
#' Estimates the risk \eqn{R_{\Gamma_0}(t; h)} associated with the covariance segment line estimation
#' proposed by \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param t A numeric vector. Observation points where the mean function of the underlying process is estimated.
#' @param bw_grid A numeric vector. A bandwidth grid from which the best smoothing parameter is selected for each \code{t}.
#' Default is \code{NULL}, in which case it is defined as an exponential grid of \eqn{N \times \lambda}.
#' @param center Logical. If \code{TRUE}, centers the data before estimation. Default is \code{TRUE}.
#' @param kernel_name Character string. Specifies the kernel function for estimation; default is \code{"epanechnikov"}.
#' Supported kernels include: \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"},
#' \code{"triangular"}, and \code{"uniform"}.
#'
#' @return A \link[data.table]{data.table} with columns:
#' \itemize{
#'   \item \code{t}: The points at which the risk function is estimated.
#'   \item \code{h}: The candidate bandwidth.
#'   \item \code{PN}: The number of curves used to estimate the mean at \code{t}, i.e., \eqn{P_N(t;h)}.
#'   \item \code{locreg_bw}: The bandwidth used to estimate the local regularity parameters.
#'   \item \code{Ht}: The estimates of the local exponent \eqn{H_t}.
#'   \item \code{Lt}: The estimates of the Hölder constant \eqn{L_t^2}.
#'   \item \code{bias_term}: The bias term of the risk function.
#'   \item \code{variance_term}: The variance term of the risk function.
#'   \item \code{dependence_term}: The dependence term of the risk function.
#'   \item \code{cov_segment_risk}: The estimated risk of the covariance segment function.
#' }
#'
#' @details
#' The local regularity parameters are estimated within the function using \link{estimate_locreg_cpp}.
#'
#' The dependence term includes contributions from both a term based on \eqn{\mathbb{D}(t; h_t)} derived from fourth-moment tensors,
#' and an empirical autocovariance term computed using \link{estimate_empirical_XsXt_autocov_cpp}.
#'
#' @export
#' @seealso \link{estimate_mean}, \link{estimate_locreg}, \link{estimate_sigma},
#'          \link{estimate_nw}, \link{estimate_empirical_autocov}
#'
#' @import data.table
#' @importFrom Rdpack reprompt
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Example coming soon
#'
estimate_cov_segment_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                      t = c(1/4, 1/2, 3/4),
                                      bw_grid = NULL,
                                      center = TRUE,
                                      kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (! is.null(bw_grid)) {
    if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
      stop("If 'bw_grid' is not NULL, it must be a vector of positive values between 0 and 1.")
  } else {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    K <- 20
    b0 <- 4 * (N * lambdahat) ** (- 0.9)
    bK <- 4 * (N * lambdahat) ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(K, b0, bK, a, lambdahat) ; gc()
  }

  # Estimate risk function using C++  function
  dt_risk <- estimate_cov_segment_risk_cpp(
    data = data, t = t, bw_grid = bw_grid,
    center = center, kernel_name = kernel_name)
  dt_risk <- data.table::as.data.table(dt_risk)
  data.table::setnames(x = dt_risk,
                       new = c("t", "h", "PN", "locreg_bw", "Ht", "Lt", "bias_term",
                               "variance_term", "dependence_term", "cov_segment_risk"))
  return(dt_risk)
}


#' Estimate Covariance Segment Function for Functional Data
#'
#' Estimates the covariance segment function \eqn{\Gamma_{N,0}(t,t;h_t,h_t)} for functional data
#' using the Nadaraya–Watson estimator with a specified kernel. This is part of the methodology
#' described in \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
#'
#' @inheritParams estimate_cov_segment_risk
#' @param optbw A numeric vector. Optimal bandwidth parameters for covariance segment function estimation at each \code{t}.
#' Default is \code{NULL}, in which case it will be estimated using the \link{estimate_cov_segment_risk} function.
#'
#' @return A \link[data.table]{data.table} containing the following columns:
#' \itemize{
#'   \item{\code{t} :}{ The observation points at which the covariance segment function is estimated.}
#'   \item{\code{optbw} :}{ The optimal bandwidth used to estimate covariance segment function at each \code{t}.}
#'   \item{\code{Ht} :}{ Local exponent estimates for each \code{t}, corresponding to \eqn{H_t}.}
#'   \item{\code{Lt} :}{ Estimates of the Hölder constant for each \code{t}, corresponding to \eqn{L_t^2}.}
#'   \item{\code{PN} :}{ The number of selected curves used in the estimation for each \code{t}.}
#'   \item{\code{cov_segment_hat} :}{ Uncorrected covariance segment estimate. }
#'   \item{\code{covseg_correction} :}{ Correction term based on measurement error variance. }
#'   \item{\code{cov_segment_hat_corrected} :}{ Final corrected covariance segment estimate. }
#' }
#'
#' @export
#'
#' @seealso \link{estimate_cov_segment_risk}, \link{estimate_locreg}, \link{estimate_sigma},
#'          \link{estimate_nw}, \link{estimate_empirical_autocov}
#'
#' @import data.table
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Example coming soon
#'
estimate_cov_segment <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                 t = c(1/4, 1/2, 3/4),
                                 optbw = NULL,
                                 bw_grid = NULL,
                                 center = TRUE,
                                 kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (! is.null(bw_grid)) {
    if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
      stop("If 'bw_grid' is not NULL, it must be a vector of positive values between 0 and 1.")
  } else {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    K <- 20
    b0 <- 4 * (N * lambdahat) ** (- 0.9)
    bK <- 4 * (N * lambdahat) ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(K, b0, bK, a, lambdahat) ; gc()
  }

  # Estimate covariance segment function using C++  function
  dt_res <- estimate_cov_segment_cpp(data = data, t = t, optbw = optbw, bw_grid = bw_grid, center = center, kernel_name = kernel_name)
  dt_res <- data.table::as.data.table(dt_res)
  data.table::setnames(x = dt_res, new = c("t", "optbw", "Ht", "Lt", "PN", "cov_segment_hat",
                                           "covseg_correction", "cov_segment_hat_corrected"))
  dt_res[cov_segment_hat < covseg_correction, cov_segment_hat_corrected := cov_segment_hat]
  return(dt_res)
}
