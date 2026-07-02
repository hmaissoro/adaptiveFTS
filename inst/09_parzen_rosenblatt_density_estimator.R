#' Leave-one-out Parzen-Rosenblatt density estimator with LSCV bandwidth
#'
#' Computes the leave-one-out Parzen-Rosenblatt estimator of the design
#' density evaluated at the observation points themselves, as required by the
#' independent-design weights \eqn{\varrho_{n_0,i} = \{M_{n_0}\,\widehat
#' g_{n_0,i}(T_{n_0,i})\}^{-1}} of the adaptive BLUP. For a single curve with
#' observation times \eqn{T_{n_0,i}}, \eqn{1 \le i \le M_{n_0}},
#' \deqn{\widehat g_{n_0,i}(T_{n_0,i}) =
#'       \frac{1}{(M_{n_0}-1)\,h} \sum_{j \ne i}
#'       K\!\left(\frac{T_{n_0,i} - T_{n_0,j}}{h}\right).}
#' The bandwidth \eqn{h} is selected by least-squares (Rudemo-Bowman)
#' cross-validation, which minimises an unbiased estimate (up to an
#' \eqn{h}-independent constant) of the density MISE:
#' \deqn{CV(h) = \int \widehat g_h(x)^2 dx - \frac{2}{M} \sum_i \widehat
#'       g_h^{(-i)}(T_{n_0,i}).}
#'
#' @param x Numeric vector of observation times of a *single* curve
#'   \eqn{T_{n_0,\cdot}} (not pooled across curves).
#' @param bw_grid Numeric vector of candidate bandwidths for the LSCV search.
#'   If `NULL` (default), a fixed log-spaced grid close to
#'   `exp(seq(log(0.01), log(0.3), length.out = 30))` is used. For Monte Carlo
#'   studies, pass a fixed grid so the cross-validation is comparable across
#'   replications.
#' @param kernel_name Kernel name: "epanechnikov", "gaussian", "uniform",
#'   or "biweight".
#' @param lower,upper Bounds of the domain \eqn{I}; default (0, 1].
#'
#' @return A list with `h_star` (selected bandwidth), `bw_grid`, `cv_curve`,
#'   `kernel_name`, and `estimate` (the vector \eqn{\widehat
#'   g_{n_0,i}(T_{n_0,i})}, one value per observation point, in the order of `x`).
#'
#' @export
estimate_density <- function(x, bw_grid = NULL,
                             kernel_name = "epanechnikov",
                             lower = 0, upper = 1) {
  n <- length(x)
 
  if (is.null(bw_grid)) {
    # Default to a fixed log-spaced grid so bandwidth selection remains
    # aligned with the reference Monte Carlo setup.
    bw_grid <- exp(seq(log(0.01), log(0.3), length.out = 30))
  }
 
  kern <- switch(
    kernel_name,
    epanechnikov = epanechnikov_kernel,
    gaussian     = gaussian_kernel,
    uniform      = uniform_kernel,
    biweight     = biweight_kernel,
    stop("Unknown kernel_name: ", kernel_name)
  )
 
  cv_score <- function(h) {
    # Integral term int ghat_h(x)^2 dx via trapezoidal rule
    zgrid <- seq(lower, upper, length.out = 500L)
    fvals <- vapply(zgrid, function(zi) sum(kern((zi - x) / h)) / (n * h), numeric(1))
    integral_term <- pracma::trapz(zgrid, fvals)
 
    # Leave-one-out term
    loo_term <- vapply(seq_len(n), function(i) {
      sum(kern((x[i] - x[-i]) / h)) / ((n - 1) * h)
    }, numeric(1))
 
    integral_term - 2 * mean(loo_term)
  }
 
  cv_curve <- vapply(bw_grid, function(h) {
    tryCatch(cv_score(h), error = function(e) NA_real_)
  }, numeric(1))
  if (!any(is.finite(cv_curve)))
    stop("CV score could not be computed for any h in bw_grid.")
 
  h_star <- bw_grid[which.min(cv_curve)]
  if (h_star %in% range(bw_grid))
    warning("h_star is at a bw_grid boundary; consider widening the grid.")
 
  # Leave-one-out estimate at each observation point, using h_star
  estimate <- vapply(seq_len(n), function(i) {
    sum(kern((x[i] - x[-i]) / h_star)) / ((n - 1) * h_star)
  }, numeric(1))
 
  list(
    h_star = h_star,
    bw_grid = bw_grid,
    cv_curve = cv_curve,
    kernel_name = kernel_name,
    estimate = estimate
  )
}
 

## Example / sanity check

if (FALSE) {
  data("data_far")
  Tn0 <- data_far[id_curve == 149, sort(tobs)]
  Mn0 <- length(Tn0)
  bw_grid <- exp(seq(log(0.01), log(0.3), length.out = 30))
 
  fit <- estimate_density(
    x = Tn0, bw_grid = bw_grid,
    kernel_name = "epanechnikov",
    lower = 0, upper = 1
  )
 
  fit$h_star
 
  # CV curve diagnostic
  plot(fit$bw_grid, fit$cv_curve, type = "b", log = "x",
       xlab = "bandwidth h", ylab = "CV(h)",
       main = "Leave-one-out CV score")
  abline(v = fit$h_star, col = "red", lty = 2)
 
  # ghat_{n0,i}(T_{n0,i}) vs. the true design density at the same points
  plot(Tn0, fit$estimate, pch = 16, cex = 0.5,
       xlab = "T_{n0,i}", ylab = "density",
       main = "LOO Parzen-Rosenblatt estimate at the observation points",
       type="o", col = "blue")
  legend("topright", legend = "ghat_{n0,i}",
         col = "blue", pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2))
 
  pracma::trapz(Tn0, fit$estimate)
 
  rho <- 1 / (Mn0 * fit$estimate)
  sum(rho)
}