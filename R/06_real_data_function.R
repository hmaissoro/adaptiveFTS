#' Mean function learned from the voltage curves of the electricity
#'
#' For more details see the vignette:
#' \code{vignette("hybrid-simulation-setup", package = "adaptiveFTS")}
#'
#' @param t \code{vector (numeric)}. Points at which we want to return the mean function.
#' It can be a scalar.
#'
#' @return A \code{data.table} containing 2 columns.
#'          \itemize{
#'            \item{t :}{ The vector or scalar \code{t}.}
#'            \item{mean :}{ The values of the mean function evaluated at \code{t}.}
#'         }
#' @export
#'
#' @importFrom data.table data.table
#'
#' @examples
#'
#' t0 <- seq(0.1, 0.9, len = 10)
#' m <- get_real_data_mean(t = t0)
#' plot(x = t0, y = m, type = "b", col = "red",
#'      xlab = "t", ylab = "mean", main = "Mean function")
#'
get_real_data_mean <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  cost_mat <- outer(X = t, Y = 1:50, function(ti, k) sqrt(2) * cos(2 * pi * k * ti))
  sint_mat <- outer(X = t, Y = 1:50, function(ti, k) sqrt(2) * sin(2 * pi * k * ti))
  eta <- cbind(1, cost_mat, sint_mat)

  # Basis coeffient
  basis_coef <- c(
    2.424666e+02, 2.045053e-01, 1.428404e-01, 3.570635e-01,
    1.090163e-01, 2.118234e-01, 1.301123e-01, 9.508016e-02,
    1.379302e-01, 4.820783e-03, -6.806055e-02, -2.388833e-02,
    -7.403338e-02, -2.313821e-02, -5.963765e-02, -2.788351e-02,
    6.017267e-02, 4.237483e-03, 2.406135e-02, 9.984997e-03,
    2.732683e-02, -3.947718e-02, -2.744963e-03, -5.624731e-03,
    -7.232138e-02, 1.228444e-02, -2.983196e-02, -1.373429e-02,
    -5.086314e-03, -5.206650e-03, 1.983353e-02, -1.671532e-02,
    1.694785e-02, 1.663588e-02, -8.924121e-03, -1.470650e-02,
    -1.568260e-02, -3.873785e-03, -1.642147e-02, 0.000000e+00,
    1.706651e-04, 3.662220e-03, 3.599005e-03, 2.163418e-02,
    2.079180e-02, -5.805163e-03, -6.198625e-03, -3.051126e-03,
    -3.050078e-02, -2.183053e-02, -2.030866e-02, 6.524725e-01,
    1.448886e+00, -3.113890e-01, -4.050726e-01, 1.478745e-01,
    -1.578363e-01, -2.346072e-01, -5.265881e-02, -5.181928e-02,
    -8.382788e-02, -5.055031e-02, 1.471149e-01, -8.402212e-03,
    -4.316013e-02, 7.528717e-02, 3.718024e-02, -4.602782e-03,
    -4.930040e-02, -7.104138e-03, -3.485272e-02, -5.034491e-02,
    2.230170e-02, 5.058664e-02, -2.996308e-02, 0.000000e+00,
    1.773518e-02, 1.664768e-03, -5.118570e-04, 2.536951e-02,
    1.103531e-02, -3.781447e-02, -9.837124e-03, 3.219296e-03,
    -1.163841e-02, -1.604513e-02, -8.183994e-03, 3.309498e-02,
    7.700235e-03, 1.578432e-02, 5.755486e-03, -3.571603e-03,
    -1.118589e-03, -1.883942e-03, -5.265843e-03, -2.892450e-02,
    1.032219e-02, 1.451413e-02, 6.348425e-04, 1.621365e-02,
    1.322576e-02
  )

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(cost_mat, sint_mat, eta, basis_coef)
  gc()

  return(muhat[, 1])
}

#' FAR kernel learned from the voltage curves of the electricity
#'
#' For more details see the vignette:
#' \code{vignette("hybrid-simulation-setup", package = "adaptiveFTS")}
#'
#' @param s \code{numeric (positive)}. A vector or scalar value(s) between 0 and 1.
#' @param t \code{numeric (positive)}. A vector or scalar value(s) between 0 and 1.
#' @param operator_norm \code{numeric (positive)}. A scalar corresponding to the norm of the integral operator associated with this kernel function.
#'
#' @return A vector (or scalar) of \code{numeric} values corresponding to the value of the kernel function evaluated at (\code{s}, \code{t}).
#' @export
#'
#' @examples
#'
#' # get the value of the kernel at (s,t) = (0.2, 0.3)
#' kerval <- get_real_data_far_kenel(s = 0.2, t = 0.3, operator_norm = 0.5)
#' kerval
#'
get_real_data_far_kenel <- function(s = 0.2, t = 0.3, operator_norm = 0.5){

  # Basis coefficient
  # For each fixed {\eta_k(s), k = 1,...,K} and {\theta_l(t), l = 1,...,L}, we have
  # c(b_{11}, b_{12}, ..., b_{1L},
  #   b_{21}, b_{22}, ..., b_{2L},
  #   ...,
  #   b_{K1}, b_{K2}, ..., b_{KL})
  basis_coef <- c(
    2.23373729709883, -4.82923610791908, 3.32335705082156, 1.79700259321478,
    -4.72746806843297, -5.73918751872499, 2.45125518767258, 1.07367258654584,
    -5.43802938214564, -3.51254535084025, -0.983346082229677, -0.325244044730888,
    0.331031200463371, -0.0540202044257287, 0, 0.050670616459408, 0.245848419309361,
    -0.0463018677073661, 0.185132606151136, -0.161785518978632, -0.114457510952933,
    0.0820977879166196, 0.187206920919163, 0.183596326343521, 0.29799917456688,
    0.0983572083790289, 0.0365619524514346, 0, -0.054901491867189, -0.0244962109053967,
    0.191707201891362, 0.128801515419686, -0.224174908366066, -0.251566711879058,
    -0.0276348680173348, -0.0715639682939576, -0.00897146114342853, 0, 0.267121966309748,
    -0.0273091594161997, 0.067130829587406, -0.00550062556625829, 0, 0.246518722930269,
    -0.0167644005889752, 0.140334666929821, 0.0286475855004669, 0, 0.28472302948554, 0,
    -0.0604572008597386, -0.0756788959584134, 0.0111654607748143, -0.118712184576549, 0,
    0.276667058786185, 0, 0.205067486705622, 0.0654171627574894, 0.012651336931906,
    -0.0896598422692541, 0.0750983773085757, -0.0390603398328494, 0.255463915244155, 0,
    -0.0538389230934492, 0.0615073725777393, -0.133119857316778, -0.131347879232277,
    0.131319341474459, 0, 0.0316171381996287, 0.1366837259548, 0.191576162910903,
    0.0821500379023849, 0.245074852934185, -0.161386348784381, 0.57160301833627,
    -0.0895419722592298, 0.353557986549909, -0.00925462916880733, -0.0698838250503909,
    -0.497687159750738, 0.142855565715012, -0.0250687004178898, 0.445786769743267,
    0.20097497729795, -0.0059528113317605, -0.572809375110091, 0.00307262837101596,
    -0.162035292129038, -0.185062678568657, 0, 0.524409493593714, -0.297541175306963,
    0.205775044697609, -0.354321377285357, -0.302443276571187, 0.213730662847565,
    0.0141379256684164, -0.0483508443593266, -0.0391148946905164, 0.0906326014808357,
    0, 0, -0.0117502269372805, -0.0517597384977126, 0.0587158994889613, 0.289741051088531,
    0, 0.0149022741092416, -0.126251329567908, -0.0459684159223364, -0.0481786322359199,
    -0.125755965631474, 0.197938756150598, 0.0270773989392422, 0.00607347821449695,
    -0.0962640169089107, 0, 0.192629550170813
  )
  # Transform to (K, L) matrix
  basis_coef_mat <- t(matrix(data = basis_coef, ncol = 11))

  ker_values <- mapply(function(s,t, coef_mat){
    # \eta(s)
    etas <- c(1, sqrt(2) * cos(2 * pi * 1:5 * s), sqrt(2) * sin(2 * pi * 1:5 * s))

    # \theta(t)
    thetat <- c(1, sqrt(2) * cos(2 * pi * 1:5 * t), sqrt(2) * sin(2 * pi * 1:5 * t))

    # Basis function
    ker_val <- matrix(etas, nrow = 1) %*% coef_mat %*% matrix(thetat, ncol = 1)
    return(c(ker_val))
  }, s = s, t = t, MoreArgs = list(coef_mat = basis_coef_mat))

  # Normalize values using operator norm
  op_norm <- 4.588783
  op_scale <- operator_norm / op_norm
  ker_values <- ker_values * op_scale

  # clean
  rm(basis_coef_mat, basis_coef)
  gc()

  return(ker_values)
}

