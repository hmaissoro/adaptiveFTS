#' Select the kernel function used by the density estimator
#'
#' Maps a kernel name to the corresponding (vectorised) kernel evaluated by the
#' package's compiled kernels. Only the kernels shared with the rest of the
#' package are supported.
#'
#' @param kernel_name A string giving the kernel name.
#'
#' @return A function taking a numeric vector `u` and returning the kernel
#'   values `K(u)`.
#' @keywords internal
.select_density_kernel <- function(kernel_name) {
  return(switch(
    kernel_name,
    epanechnikov = epanechnikov_kernel,
    biweight     = biweight_kernel,
    triweight    = triweight_kernel,
    tricube      = tricube_kernel,
    triangular   = triangular_kernel,
    uniform      = uniform_kernel,
    stop("Unsupported kernel name. Choose from: epanechnikov, biweight, ",
         "triweight, tricube, triangular, uniform.")
  ))
}

#' Trapezoidal integration
#'
#' Base-R trapezoidal rule, used to avoid a dependency on `pracma::trapz`.
#'
#' @param x,y Numeric vectors of equal length giving the abscissa and ordinate.
#'
#' @return The trapezoidal approximation of \eqn{\int y \, dx}.
#' @keywords internal
.trapz <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(0)
  return(sum((x[-1L] - x[-n]) * (y[-1L] + y[-n]) / 2))
}

#' Leave-one-out Parzen-Rosenblatt estimate at the observation points
#'
#' @param x Observation times of a single curve.
#' @param hval Bandwidth.
#' @param kern Kernel function.
#' @return The vector \eqn{\widehat g_{h}^{(-i)}(x_i)}, in the order of `x`.
#' @keywords internal
.density_loo_estimate <- function(x, hval, kern) {
  n <- length(x)
  return(vapply(seq_len(n), function(i) {
    sum(kern((x[i] - x[-i]) / hval)) / ((n - 1) * hval)
  }, numeric(1)))
}

#' Least-squares cross-validation score of the design-density bandwidth
#'
#' @param hval Bandwidth.
#' @param x Observation times of a single curve.
#' @param kern Kernel function.
#' @param lower,upper Bounds of the domain.
#' @return The LSCV score at `hval`.
#' @keywords internal
.density_cv_score <- function(hval, x, kern, lower, upper) {
  n <- length(x)
  zgrid <- seq(lower, upper, length.out = 500L)
  fvals <- vapply(zgrid, function(zi) sum(kern((zi - x) / hval)) / (n * hval), numeric(1))
  integral_term <- .trapz(zgrid, fvals ** 2)
  loo_term <- .density_loo_estimate(x, hval, kern)
  return(integral_term - 2 * mean(loo_term))
}

#' Leave-one-out Parzen-Rosenblatt density estimator
#'
#' Computes the leave-one-out Parzen-Rosenblatt estimator of the design density
#' evaluated at the observation points themselves, as required by the
#' independent-design weights of the adaptive BLUP.
#'
#' @details
#' For a single curve with observation times \eqn{T_{n_0,i}},
#' \eqn{1 \le i \le M_{n_0}}, the estimator at the observation points is
#' \deqn{\widehat g_{n_0,i}(T_{n_0,i}) =
#'       \frac{1}{(M_{n_0}-1)\,h} \sum_{j \ne i}
#'       K\!\left(\frac{T_{n_0,i} - T_{n_0,j}}{h}\right),}
#' feeding the independent-design weights \eqn{\varrho_{n_0,i} =
#' \{M_{n_0}\,\widehat g_{n_0,i}(T_{n_0,i})\}^{-1}}. The bandwidth \eqn{h} is
#' either supplied directly (via `h`) or selected by least-squares
#' (Rudemo-Bowman) cross-validation over `bw_grid`, which minimises an unbiased
#' estimate (up to an \eqn{h}-independent constant) of the density MISE:
#' \deqn{CV(h) = \int \widehat g_h(x)^2 dx - \frac{2}{M} \sum_i \widehat
#'       g_h^{(-i)}(T_{n_0,i}).}
#'
#' @param x Numeric vector of observation times of a *single* curve
#'   \eqn{T_{n_0,\cdot}} (not pooled across curves).
#' @param h Optional numeric scalar bandwidth. If supplied, the least-squares
#'   cross-validation is skipped and this bandwidth is used directly. This is
#'   the path used by the adaptive BLUP, where a single bandwidth selected once
#'   on a subset of curves (see [get_density_optimal_bw()]) is reused for every
#'   curve. Default is `NULL`, in which case the bandwidth is selected by LSCV.
#' @param bw_grid Numeric vector of candidate bandwidths for the LSCV search
#'   (used only when `h` is `NULL`). If `NULL` (default), a fixed log-spaced
#'   grid close to `exp(seq(log(0.01), log(0.3), length.out = 30))` is used. For
#'   Monte Carlo studies, pass a fixed grid so the cross-validation is
#'   comparable across replications.
#' @param kernel_name Kernel name: "epanechnikov" (default), "biweight",
#'   "triweight", "tricube", "triangular", or "uniform".
#' @param lower,upper Bounds of the domain \eqn{I}; default (0, 1].
#'
#' @return A list with:
#'   \itemize{
#'     \item `h_star`: the bandwidth used.
#'     \item `bw_grid`: the candidate bandwidth grid (as supplied or defaulted).
#'     \item `cv_curve`: the LSCV score per candidate bandwidth, or `NULL` when
#'       `h` is supplied.
#'     \item `kernel_name`: the kernel used.
#'     \item `estimate`: the vector \eqn{\widehat g_{n_0,i}(T_{n_0,i})}, one value
#'       per observation point, in the order of `x`.
#'   }
#'
#' @export
estimate_density <- function(x, h = NULL, bw_grid = NULL,
                             kernel_name = "epanechnikov",
                             lower = 0, upper = 1) {
  n <- length(x)
  if (n < 2L)
    stop("'x' must contain at least two observation points.")

  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )
  kern <- .select_density_kernel(kernel_name)

  if (!is.null(h)) {
    if (!(methods::is(h, "numeric") && length(h) == 1L && h > 0))
      stop("'h' must be a single positive numeric value.")
    return(list(
      h_star = h,
      bw_grid = bw_grid,
      cv_curve = NULL,
      kernel_name = kernel_name,
      estimate = .density_loo_estimate(x, h, kern)
    ))
  }

  if (is.null(bw_grid)) {
    bw_grid <- exp(seq(log(0.01), log(0.3), length.out = 30))
  }

  cv_curve <- vapply(bw_grid, function(hval) {
    tryCatch(.density_cv_score(hval, x, kern, lower, upper), error = function(e) NA_real_)
  }, numeric(1))
  if (!any(is.finite(cv_curve)))
    stop("CV score could not be computed for any h in bw_grid.")

  h_star <- bw_grid[which.min(cv_curve)]
  if (h_star %in% range(bw_grid))
    warning("h_star is at a bw_grid boundary; consider widening the grid.")

  return(list(
    h_star = h_star,
    bw_grid = bw_grid,
    cv_curve = cv_curve,
    kernel_name = kernel_name,
    estimate = .density_loo_estimate(x, h_star, kern)
  ))
}

#' Select the design-density bandwidth on a subset of curves
#'
#' Selects a single Parzen-Rosenblatt bandwidth for the design density by
#' least-squares cross-validation on a subset of curves, mirroring
#' [get_nw_optimal_bw()].
#'
#' @details
#' For each sampled curve the bandwidth minimising the LSCV score of
#' [estimate_density()] over `bw_grid` is computed, and the median of these
#' per-curve optima is returned. The adaptive BLUP calls this once and then
#' reuses the returned bandwidth for every curve (passing it as the `h` argument
#' of [estimate_density()]), so the design weights are consistent across curves
#' and the cross-validation is not repeated on every call.
#'
#' @inheritParams format_data
#' @param nsubset \code{integer (positive)}. The number of curves to randomly
#'   and uniformly select for bandwidth selection. Default is `NULL`, in which
#'   case every curve is used.
#' @param bw_grid Numeric vector of candidate bandwidths passed to
#'   [estimate_density()]. Default is `NULL` (a fixed log-spaced grid). For
#'   reproducible Monte Carlo studies, pass a fixed grid.
#' @param kernel_name Kernel name: "epanechnikov" (default), "biweight",
#'   "triweight", "tricube", "triangular", or "uniform".
#' @param lower,upper Bounds of the domain \eqn{I}; default (0, 1].
#'
#' @return A \code{numeric} scalar: the median of the per-curve LSCV-optimal
#'   bandwidths over the subset.
#'
#' @seealso [estimate_density()], [get_nw_optimal_bw()].
#'
#' @export
#'
#' @import data.table
#' @importFrom methods is
get_density_optimal_bw <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                   nsubset = NULL, bw_grid = NULL,
                                   kernel_name = "epanechnikov",
                                   lower = 0, upper = 1) {
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  ids <- data[, sort(unique(id_curve))]
  N <- length(ids)

  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  if (!is.null(nsubset))
    if (any(nsubset < 0) | (length(nsubset) > 1) | any(nsubset - floor(nsubset) > 0) | any(N <= nsubset))
      stop("If 'nsubset' is not NULL, then it must be a positive integer lower than the number of curves.")

  sample_ids <- if (!is.null(nsubset)) sample(x = ids, size = nsubset) else ids

  h_stars <- vapply(sample_ids, function(idc) {
    tvals <- data[id_curve == idc, sort(unique(tobs))]
    tryCatch(
      suppressWarnings(estimate_density(
        x = tvals, h = NULL, bw_grid = bw_grid,
        kernel_name = kernel_name, lower = lower, upper = upper)$h_star),
      error = function(e) NA_real_
    )
  }, numeric(1))

  if (!any(is.finite(h_stars)))
    stop("The density bandwidth could not be selected on any curve of the subset.")

  return(stats::median(h_stars, na.rm = TRUE))
}
