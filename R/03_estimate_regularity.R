#' Local Regularity Parameters Estimation
#'
#' @param data \code{data.table (or data.frame)} or \code{list} of \code{data.table (or data.frame)} or \code{list} of \code{list}.
#' \itemize{
#'    \item{If \code{data.table}}{
#'        It must contain the raw binding of the curve observations with at least 3 columns.
#'        \itemize{
#'          \item{\code{idcol} :}{ The name of the column that contains the index of the curve in the sample.
#'                              Each index of a curve is repeated as many times as it has observation points.}
#'          \item{\code{tcol} :}{ The name of the column that contains the observation points associated to each curve index.}
#'          \item{\code{ycol} :}{ The name of the column that contains the observed value of a curve at each point of observation and for each index of the curve.}
#'        }
#'    }
#'    \item{If \code{list} of \code{data.table}}{
#'         In this case, each element of the given \code{list} corresponds to the observation scheme of a curve, which is given as \code{data.table} or \code{data.frame}.
#'         The data.table contains at least 2 columns.
#'         \itemize{
#'          \item{\code{tcol} :}{ The name of the column that contains the observation points associated to the curve.}
#'          \item{\code{ycol} :}{ The name of the column that contains the observed value of the curve.}
#'        }
#'    }
#'    \item{If \code{list} of \code{list}}{
#'      In the latter case, the \code{data} is a list \code{list} where each element is the observation scheme of a curve given as a \code{list} of 2 vectors.
#'      \itemize{
#'          \item{\code{tcol} :}{ The name of the vector that contains the observation points associated the curve.}
#'          \item{\code{ycol} :}{ The name of the vector that contains the observed value of the curve.}
#'        }
#'    }
#' }
#' @param idcol \code{character}. If \code{data} is given as \code{data.table} or \code{data.frame},
#' it is the name of the column that contains the index of the curve in the sample.
#' Each index of a curve is repeated as many times as it has observation points.
#' Opposite, if f \code{data} is given as \code{list} of \code{data.table (or data.frame)} of \code{list} of \code{list}, \code{idcol = NULL.}
#' @param tcol \code{character}. The name of the column (or vector) that contains the observation points associated to the curves.
#' @param ycol \code{character}. The name of the column that contains the observed value of the curves.
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the local regularity parameters of the underlying process.
#' @param Delta \code{numeric (positive)}. The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.
#' Default \code{Delta = NULL} and thus it we be estimated from the data.
#' @param h \code{numeric (positive)}. The bandwidth of the Nadaraya-Watson estimator. Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#' @param weighted \code{logical}. If \code{TRUE}, a weighted versions of the local regularity estimators are estimated.
#' Default \code{weighted = FALSE}. Recall that to estimate the local regularity at $t_i$ we need to estimate each curve at $t_i - Delta/2$, $t_i$ and t_i + Delta/2$.
#' The weighting is to consider only those curves for which there is at least one observation point in the smoothing window around $t_i - Delta/2$, $t_i$ and t_i + Delta/2$.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{t :}{ The points around which the local regularity parameters are estimated.}
#'            \item{h :}{ The presmoothing bandwidth.}
#'            \item{Delta :}{ The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.}
#'            \item{H :}{ The local exponent estimates.}
#'            \item{L :}{ The Hölder constant estimates. It corresponds to $L_t^2$.}
#'            \item{Wn :}{ The number of curve used to estimate the local regularity parameters \code{weighted = TRUE}.}
#'         }
#' @export
#'
#' @import data.table
#' @importFrom methods is
#' @importFrom stats median
#'
#' @seealso [estimate_nw()], [estimate_nw_bw()], [simulate_far()], etc.
#'
#' @examples
#' \dontrun{
#'   # Generate a sample of FAR(1)
## Exponent H
#' Hfun <- function(t) {
#'   hurst_logistic(t = t, h_left = 0.4, h_right = 0.8, slope = 5)
#' }
#'
#' ## Hölder constant
#' L <- 4
#'
#' dt_far <- simulate_far(N = 200L, lambda = 100L,
#'                        tdesign = "random",
#'                        tdistribution = runif,
#'                        tcommon = NULL,
#'                        hurst_fun = Hfun,
#'                        L = L,
#'                        far_kernel = function(s,t) 9/4 * exp( - (t + 2 * s) ** 2),
#'                        far_mean = function(t) 4 * sin(1.5 * pi * t),
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#'
#' # Estimate local regularity at
#' t0 <- seq(0.2, 0.8, len = 8)
#'
#' ## If data is a data.table or a data. frame
#' dt_locreg <- estimate_locreg(data = dt_far,
#'                              idcol = "id_curve",
#'                              tcol = "tobs",
#'                              ycol = "X",
#'                              t = t0,
#'                              Delta = NULL,
#'                              h = NULL,
#'                              smooth_ker = epanechnikov,
#'                              weighted = FALSE)
#' DT::datatable(dt_locreg)
#'
#' ## If data is a list of data.table (or data. frame)
#' list_dt_far <- lapply(unique(dt_far[, id_curve]), function(idx){
#'   dt_far[id_curve == idx, .(tobs, X)]
#' })
#'
#' dt_locreg_2 <- estimate_locreg(data = list_dt_far,
#'                                idcol = NULL,
#'                                tcol = "tobs",
#'                                ycol = "X",
#'                                t = t0,
#'                                Delta = NULL,
#'                                h = NULL,
#'                                smooth_ker = epanechnikov,
#'                                weighted = FALSE)
#' DT::datatable(dt_locreg_2)
#'
#' ## If data is a list of list
#' list_list_far <- lapply(unique(dt_far[, id_curve]), function(idx){
#'   list("Obs_point" = dt_far[id_curve == idx, tobs],
#'        "Xobs" = dt_far[id_curve == idx, X])
#' })
#'
#' dt_locreg_3 <- estimate_locreg(data = list_list_far,
#'                                idcol = NULL,
#'                                tcol = "Obs_point",
#'                                ycol = "Xobs",
#'                                t = t0,
#'                                Delta = NULL,
#'                                h = NULL,
#'                                smooth_ker = epanechnikov,
#'                                weighted = TRUE)
#' DT::datatable(dt_locreg_2)
#'
#'}
#'
#'
#'
estimate_locreg <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                            t = 1/2, Delta = NULL, h = NULL,
                            smooth_ker = epanechnikov, weighted = FALSE){
  # Control see checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (! methods::is(weighted, "logical"))
    stop("'weighted' must be TRUE or FALSE")

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)

  # Control and set arguments depending on the data
  if (! is.null(Delta)) {
    if (! (methods::is(Delta, "numeric") & data.table::between(Delta, 0, 1) & length(Delta) == 1))
      stop("'Delta' must be a numeric scalar value between 0 and 1.")
  } else {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    Delta <- 2 * exp(-log(lambdahat) ** 0.72)
  }
  if (! is.null(h)) {
    if (! (methods::is(h, "numeric") & data.table::between(h, 0, 1) & length(h) == 1))
      stop("'h' must be a numeric scalar value between 0 and 1.")
  } else {
    # If h = NULL, choose the bandwidth by CV
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    N <- length(data[, unique(id_curve)])
    if (N > 50) {
      sample_curves <- sample(x = 1:N, size = 30)
    } else {
      sample_curves <- 1:N
    }
    h <- median(unlist(lapply(sample_curves, function(i, lambdahat, data, smooth_ker){
      hmax <- lambdahat ** (-1/3)
      hmin <- 1 / (2 * lambdahat)
      hgrid <- seq(hmin, hmax, len = 100)
      hbest <- estimate_nw_bw(y = data[id_curve == i, X],
                              t = data[id_curve ==i, tobs],
                              h_grid = hgrid,
                              smooth_ker = smooth_ker)

    }, lambdahat = lambdahat, data = data, smooth_ker = smooth_ker)))
  }

  # Step 1: estimate the curve at
  t1 <- t - Delta / 2
  t3 <- t + Delta / 2
  t2 <- t
  rm(t)

  dt_smooth <- rbindlist(lapply(data[, unique(id_curve)], function(i, data, h, t1, t2, t3){
    ## smooth at one at t1, t2 and t3
    m <- length(t2)
    dt <- estimate_nw(y = data[id_curve == i, X],
                      t = data[id_curve == i, tobs],
                      tnew = c(t1, t2, t3), h = h, smooth_ker = smooth_ker)
    xtilde <- dt$mhat
    inKernelSupp <- dt$inKernelSupp

    dt_out <- data.table::data.table(
      id_curve = i, t1 = t1, xt1 = xtilde[1:m], nt1_inKerSupp = inKernelSupp[1:m],
      t2 = t2, xt2 = xtilde[(m + 1):(2 * m)], nt2_inKerSupp = inKernelSupp[(m + 1):(2 * m)],
      t3 = t3, xt3 = xtilde[(2 * m + 1):(3 * m)], nt3_inKerSupp = inKernelSupp[(2 * m + 1):(3 * m)]
    )
    rm(xtilde, inKernelSupp)

    return(dt_out)
  }, h = h, data = data, t1 = t1, t2 = t2, t3 = t3))

  # Step 2 : estimate regularity parameters

  dt_reg <- rbindlist(lapply(1:length(t2), function(i, dt_smooth, t1, t2, t3, weighted) {
    ## Extract X_1(g),...,X_N(g) where g = t1, t2 or t3
    xt1 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt1]
    xt2 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt2]
    xt3 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt3]

    ## Compute Unweighed local regularity parameters
    theta_t1_t3 <- mean((xt1 - xt3) ^ 2)
    theta_t1_t2 <- mean((xt1 - xt2) ^ 2)

    H <- (log(theta_t1_t3) - log(theta_t1_t2))  / (2 * log(2))
    ## Bound H : 0.1 <= H <= 1
    # H <- min(max(H, 0.1), 1)
    L <- theta_t1_t3 / (abs(t1[i] - t3[i])**(2 * H))

    if (weighted) {
      ## Extract the number of observations T_ni used to estimate
      ## each X_1(g),...,X_N(g) where g = t1, t2 or t3
      nt1_inKerSupp <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], nt1_inKerSupp]
      nt2_inKerSupp <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], nt2_inKerSupp]
      nt3_inKerSupp <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], nt3_inKerSupp]
      w <- as.integer((nt1_inKerSupp > 1) & (nt2_inKerSupp > 1) & (nt3_inKerSupp > 1))

      theta_t1_t3w <- sum(w * (xt1 - xt3) ^ 2, na.rm = TRUE) / sum(w)
      theta_t1_t2w <- sum(w * (xt1 - xt2) ^ 2, na.rm = TRUE) / sum(w)

      Hw <- (log(theta_t1_t3w) - log(theta_t1_t2w))  / (2 * log(2))
      Lw <- theta_t1_t3w / (abs(t1[i] - t3[i]) ** (2 * Hw))

      dt_out <- data.table(t = t2[i], H = H, L = L, Hw = Hw, Lw = Lw, Wn = sum(w))

      rm(nt1_inKerSupp, nt2_inKerSupp, nt3_inKerSupp,
         theta_t1_t3, theta_t1_t3w, theta_t1_t2, theta_t1_t2w,
         H, L, Hw, Lw, xt1, xt2, xt3)
    } else {
      dt_out <- data.table(t = t2[i], H = H, L = L)
      rm(theta_t1_t3, theta_t1_t2, H, L, xt1, xt2, xt3)
    }

    return(dt_out)
  }, dt_smooth = dt_smooth, t1 = t1, t2 = t2, t3 = t3, weighted = weighted))

  dt_reg[, c("h", "Delta") := list(h, Delta)]

  if (weighted) {
    data.table::setcolorder(x = dt_reg, neworder = c("t", "h", "Delta", "H", "Hw", "Lw", "Wn"))
  } else {
    data.table::setcolorder(x = dt_reg, neworder = c("t", "h", "Delta", "H", "L"))
  }

  return(dt_reg)
}

#' Format data for local regularity estimation
#'
#' @param data \code{data.table (or data.frame)} or \code{list} of \code{data.table (or data.frame)} or \code{list} of \code{list}.
#' \itemize{
#'    \item{If \code{data.table}}{
#'        It must contain the raw binding of the curve observations with at least 3 columns.
#'        \itemize{
#'          \item{\code{idcol} :}{ The name of the column that contains the index of the curve in the sample.
#'                              Each index of a curve is repeated as many times as it has observation points.}
#'          \item{\code{tcol} :}{ The name of the column that contains the observation points associated to each curve index.}
#'          \item{\code{ycol} :}{ The name of the column that contains the observed value of a curve at each point of observation and for each index of the curve.}
#'        }
#'    }
#'    \item{If \code{list} of \code{data.table}}{
#'         In this case, each element of the given \code{list} corresponds to the observation scheme of a curve, which is given as \code{data.table} or \code{data.frame}.
#'         The data.table contains at least 2 columns.
#'         \itemize{
#'          \item{\code{tcol} :}{ The name of the column that contains the observation points associated to the curve.}
#'          \item{\code{ycol} :}{ The name of the column that contains the observed value of the curve.}
#'        }
#'    }
#'    \item{If \code{list} of \code{list}}{
#'      In the latter case, the \code{data} is a list \code{list} where each element is the observation scheme of a curve given as a \code{list} of 2 vectors.
#'      \itemize{
#'          \item{\code{tcol} :}{ The name of the vector that contains the observation points associated the curve.}
#'          \item{\code{ycol} :}{ The name of the vector that contains the observed value of the curve.}
#'        }
#'    }
#' }
#' @param idcol \code{character}. If \code{data} is given as \code{data.table} or \code{data.frame},
#' it is the name of the column that contains the index of the curve in the sample.
#' Each index of a curve is repeated as many times as it has observation points.
#' Opposite, if f \code{data} is given as \code{list} of \code{data.table (or data.frame)} of \code{list} of \code{list}, \code{idcol = NULL.}
#' @param tcol \code{character}. The name of the column (or vector) that contains the observation points associated to the curves.
#' @param ycol \code{character}. The name of the column that contains the observed value of the curves.
#'
#' @return A \code{data.table} containing 3 columns.
#'          \itemize{
#'            \item{id_curve :}{ The index of the curve.}
#'            \item{tobs :}{ The observation points associated to each curve \code{id_curve}.}
#'            \item{X :}{ The observed values of the curve associated to \code{id_curve} at \code{tobs} observation points.}
#'         }
#'
#' @import data.table
#' @importFrom methods is
#'
.format_data <- function(data, idcol = NULL, tcol = "tobs", ycol = "X"){
  # Check if data is a data.table or data.frame
  is_dt_or_df <- methods::is(data, "data.table") | methods::is(data, "data.frame")

  # Check if data is a list of data.table (or data.frame)
  is_list_dt_or_df <- methods::is(data, "list") &
    all(unlist(lapply(data, function(element){
      methods::is(element,"data.table") | methods::is(element,"data.frame")
    })))

  # Check if data is a list of list
  is_list_of_list <- methods::is(data, "list") &
    all(unlist(lapply(data, function(element){
      methods::is(element,"list")
    })))

  if (! (is_dt_or_df | is_list_dt_or_df | is_list_of_list))
    stop("'data' must of class data.table (or data.frame) or a list of data.table (or data.table) or a list of list.")

  if (is_dt_or_df) {
    if (is.null(idcol))
      stop("If the class of 'data' is data.table (or data.frame), 'idcol' need to be specifyed.")
    if (! all(c(idcol, tcol, ycol) %in% colnames(data))){
      stop("The specified column name 'idcol' or 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::as.data.table(data)
      data.table::setnames(x = data, old = c(idcol, tcol, ycol), new = c("id_curve", "tobs", "X"))
      data <- data[, list(id_curve, tobs, X)]
    }
  } else if (is_list_dt_or_df) {
    if (! is.null(idcol))
      stop("If 'data' is a list of data.table (or data.table) or a list of list, 'idcol' must be NULL.")
    check_colname <- all(unlist(lapply(data, function(element, tcol, ycol){
      all(c(tcol, ycol) %in% colnames(element))
    }, tcol = tcol, ycol = ycol)))
    if (! check_colname) {
      stop("The specified column name 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::rbindlist(lapply(1:length(data), function(i, tcol, ycol){
        data[[i]] <- data.table::as.data.table(data[[i]])
        data.table::setnames(x = data[[i]], old = c(tcol, ycol), new = c("tobs", "X"))
        dt <- data.table::data.table("id_curve" = i, data[[i]][, list(tobs, X)])
      }, tcol = tcol, ycol = ycol))
    }

  } else if (is_list_of_list) {
    if (! is.null(idcol))
      stop("If 'data' is a list of data.table (or data.table) or a list of list, 'idcol' must be NULL.")
    check_vecname <- all(unlist(lapply(data, function(element, tcol, ycol){
      all(c(tcol, ycol) %in% names(element))
    }, tcol = tcol, ycol = ycol)))
    if (! check_vecname) {
      stop("The specified vector name 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::rbindlist(lapply(1:length(data), function(i, tcol, ycol){
        data[[i]] <- data.table::as.data.table(data[[i]])
        data.table::setnames(x = data[[i]], old = c(tcol, ycol), new = c("tobs", "X"))
        dt <- data.table::data.table("id_curve" = i, data[[i]][, list(tobs, X)])
      }, tcol = tcol, ycol = ycol))
    }
  } else {
    NA
  }
  return(data)
}
