#' Local Regularity Estimation
#'
#' @param data \code{data.table} or \code{list} of \code{data.table} or \code{list} of \code{list}
#' @param t
#' @param k
#' @param Delta
#' @param h
#' @param kernel
#'
#' @return
#' @export
#'
#' @examples
estimate_locreg <- function(data, t = 1/2, k = 1, Delta = 1/50, h = NULL, kernel = epanechnikov){
  # INPUT
  #           data : data, sample of FTS
  #           time : where we aim to estimate the local regularity
  #              k : number of point arround t1, t2 & t3 used to estimate each curve at these locations.
  #                  For example, if for a  curve i, less then k points are used to estimate t1,
  #                  this curve will not be included in the estimation of H and L2.
  #              h : the bandwidth of the N-W estimator. Default h = NULL and thus it will be estimated by CV.
  #         kernel : the kernel function used in N-W estimator
  # OUTPUT
  #  estimate of H_time and L_time at each time

  ## Step 1:  recover the function at
  t1 <- time - Delta / 2
  t3 <- time + Delta / 2
  t2 <- time

  ### calculate the mean of the number of observation points per curve
  mu_hat <- mean(unlist(lapply(data, nrow)))

  ## If h = NULL, choose the bandwidth by CV
  if (is.null(h)){
    if (length(data) > 50) {
      sample_curves <- sample(x = 1:length(data), size = 30)
    } else {
      sample_curves <- 1:length(data)
    }
    h <- median(unlist(lapply(sample_curves, function(i, mu_hat, data, kernel){

      ## Take a maximum bandwidth (the + 0.1 is arbitrary)
      h_max <- mu_hat ** (-1/3) + 0.1
      h_vec <- seq(1 / (2 * mu_hat), h_max, len = 100)
      h <- nw_cv_bw(x = data[[i]][, t], y = data[[i]][, x], h = h_vec, kernel = kernel)

    }, mu_hat = mu_hat, data = data, kernel = kernel)))
  }

  ## smooth curves
  dt_smooth <- rbindlist(lapply(1:length(data), function(i, h, data, mu_hat, t1, t2, t3){
    ## Estimate using N-W
    m <- length(t2)
    mobj <- nw(x = data[[i]][, t], y = data[[i]][, x], xout = c(t1, t2, t3),
               h = h, kernel = kernel)
    yhat <- mobj$mhat
    inKernelSupp <- mobj$inKernelSupp

    dt_out <- data.table(index = i, t1 = t1, xt1 = yhat[1:m], nt1_inKerSupp = inKernelSupp[1:m],
                         t2 = t2, xt2 = yhat[(m + 1):(2 * m)], nt2_inKerSupp = inKernelSupp[(m + 1):(2 * m)],
                         t3 = t3, xt3 = yhat[(2 * m + 1):(3*m)], nt3_inKerSupp = inKernelSupp[(2 * m + 1):(3*m)])
    return(dt_out)
  }, h = h, data = data, mu_hat = mu_hat, t1 = t1, t2 = t2, t3 = t3))

  ## Step 2 : estimate regularity parameters

  dt_reg <- rbindlist(lapply(1:length(t2), function(i, dt_smooth, t1, t2, t3) {
    ## Extract X_1(g),...,X_N(g) where g = t1, t2 or t3
    xt1 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt1]
    xt2 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt2]
    xt3 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt3]

    ## Extract the number of observations T_ni used to estimate
    ## each X_1(g),...,X_N(g) where g = t1, t2 or t3
    nt1_inKerSupp <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], nt1_inKerSupp]
    nt2_inKerSupp <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], nt2_inKerSupp]
    nt3_inKerSupp <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], nt3_inKerSupp]

    w_t1_t3 <- ifelse(test = (nt1_inKerSupp > k) & (nt3_inKerSupp > k),
                      yes = 1, no = 0)
    w_t1_t2 <- ifelse(test = (nt1_inKerSupp > k) & (nt2_inKerSupp > k),
                      yes = 1, no = 0)
    w_t2_t3 <- ifelse(test = (nt2_inKerSupp > k) & (nt3_inKerSupp > k),
                      yes = 1, no = 0)

    ## Estimate theta
    theta_t1_t3 <- sum(w_t1_t3 * (xt1 - xt3)^2, na.rm = TRUE) / sum(w_t1_t3)
    theta_t1_t2 <- sum(w_t1_t2 * (xt1 - xt2)^2, na.rm = TRUE) / sum(w_t1_t2)
    theta_t2_t3 <- sum(w_t2_t3 * (xt3 - xt2)^2, na.rm = TRUE) / sum(w_t2_t3)

    ## Estimate H
    H1 <- (log(theta_t1_t3) - log(theta_t1_t2))  / (2 * log(2))
    # H2 <- (log(theta_t1_t3) - log(theta_t2_t3))  / (2 * log(2))
    H <- H1
    # H <- (H1 + H2) / 2

    ## Bound H : 0.1 <= H <= 1
    H <- min(max(H, 0.1), 1)

    ## Estimate L2
    L2_1 <- theta_t1_t2 / (abs(t1[i] - t2[i])**(2 * H))
    L2_2 <- theta_t1_t3 / (abs(t1[i] - t3[i])**(2 * H))
    L2_3 <- theta_t2_t3 / (abs(t2[i] - t3[i])**(2 * H))
    # L2 <- (L2_1 + L2_2 + L2_3) / 3
    # L2 <- median(L2_1 + L2_2 + L2_3)
    L2 <- L2_2
    # dt_out <- data.table(t = t2[i], H1 = H1, H2 = H2, L2_H1 = L2_H1, L2_H2 = L2_H2)
    dt_out <- data.table(t = t2[i], H = H, L = L2)

    return(dt_out)
  }, dt_smooth = dt_smooth, t1 = t1, t2 = t2, t3 = t3))

  dt_reg[, c("h", "Delta") := .(h, Delta)]
  data.table::setcolorder(x = dt_reg, neworder = c("t", "h", "Delta", "H", "L"))

  return(dt_reg)
}
