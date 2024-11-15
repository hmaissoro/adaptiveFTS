#' Local Regularity Parameters Estimation
#'
#' This function estimates the local regularity parameters \eqn{H_t} and \eqn{L_t^2} defined in Section 3 of \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points at which to estimate the local regularity parameters of the underlying process.
#' @param Delta \code{numeric (positive)}. The length of the neighborhood around each point in \code{t} for local regularity estimation.
#' Default is \code{Delta = NULL}, in which case it will be estimated from the data.
#' @param h \code{numeric (positive vector or scalar)}. Bandwidth parameter for the Nadaraya-Watson estimator used in local regularity estimation.
#' Default is \code{h = NULL}, allowing it to be determined by cross-validation on a subset of curves.
#' If \code{h} is a scalar, the same bandwidth is applied across all curves. If it is a vector, it must match the number of curves in \code{data}, with each element corresponding to a curve.
#' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov".
#' Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#' @param center \code{logical}. If \code{TRUE}, centers the curves before regularity estimation.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item{\code{t} :}{ The observation points around which local regularity is estimated.}
#'   \item{\code{locreg_bw} :}{ The presmoothing bandwidth used for estimation.}
#'   \item{\code{Delta} :}{ The neighborhood length around each \code{t} used in local regularity estimation.}
#'   \item{\code{Nused} :}{ The number of curves contributing non-degenerate estimates around each \code{t}.}
#'   \item{\code{Ht} :}{ Local exponent estimates, denoted by \eqn{H_t}.}
#'   \item{\code{Lt} :}{ Hölder constant estimates, corresponding to \eqn{L_t^2}.}
#' }
#'
#' @export
#'
#' @import data.table
#' @importFrom Rdpack reprompt
#' @importFrom methods is
#'
#' @seealso [estimate_nw()], [estimate_nw_bw()], [simulate_far()], etc.
#'
#' @references
#'  \insertAllCited{}
#'
#' @examples
#' \dontrun{
#'  # Load data
#'  data("data_far")
#'
#'  # Define observation points for local regularity estimation
#'  t0 <- seq(0.2, 0.8, length.out = 8)
#'
#'  # Estimate local regularity parameters
#'  dt_locreg <- estimate_locreg(data = data_far,
#'                               idcol = "id_curve",
#'                               tcol = "tobs",
#'                               ycol = "X",
#'                               t = t0,
#'                               Delta = NULL,
#'                               h = NULL,
#'                               kernel_name = "epanechnikov",
#'                               center = TRUE)
#'
#'  dt_locreg
#' }
#'
estimate_locreg <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                            t = 1/2, Delta = NULL, h = NULL,
                            kernel_name = "epanechnikov", center = TRUE){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! methods::is(center, "logical"))
    stop("'center' must be a TRUE or FALSE.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)

  # Estimate local regularity using C++ function
  mat_reg <- estimate_locreg_cpp(data = data, t = t, Delta = Delta, h = h,
                                    kernel_name = kernel_name, center = center)
  dt_reg <- data.table::as.data.table(mat_reg)
  data.table::setnames(x = dt_reg, new = c("t", "locreg_bw", "Delta", "Nused", "Ht", "Lt"))

  return(dt_reg)
}
