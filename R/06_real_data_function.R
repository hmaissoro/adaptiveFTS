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
#' dt_mu <- get_real_data_mean(t = seq(0.1, 0.9, len = 10))
#' plot(x = dt_mu$t, y = dt_mu$mean, type = "b", col = "red",
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
  dt_res <- data.table::data.table("t" = t, "mean" = muhat[, 1])

  return(dt_res)
}

#' FAR kernel learned from the voltage curves of the electricity
#'
#' For more details see the vignette:
#' \code{vignette("hybrid-simulation-setup", package = "adaptiveFTS")}
#'
#' @param s \code{vector (numeric)}. The points corresponding to the first argument of the kernel function.
#' It can be a scalar.
#' @param t \code{vector (numeric)}. The points corresponding to the second argument of the kernel function.
#' It can be a scalar.
#' @param operator_norm \code{numeric (positive)}. A scalar corresponding to the norm of the integral operator associated with this kernel function.
#' @param return_matrix \code{boolean}. If \code{TRUE}, return a \code{length(s) x length(t)} matrix. Default, \code{TRUE}.
#'
#' @return
#'    \itemize{
#'      \item{A \code{matrix}}{ of dimension \code{length(s) x length(t)} containing the kernel function elvatued at \code{s x t}.
#'                              Colnames t1, t2, ..., tlenght(t) and rownames are s1, s2, ..., slenght(t)}
#'      \item{Or a \code{data.table} containing 3 columns.}{
#'        \itemize{
#'          \item{s :}{ The vector or scalar \code{s}.}
#'          \item{t :}{ The vector or scalar \code{t}.}
#'          \item{ker_value :}{ For each pair \code{s x t}, the value of the kernel.}
#'        }
#'      }
#'         }
#' @export
#' @importFrom data.table data.table dcast
#'
#' @examples
#' # Return a matrix
#' mtx <- get_real_data_far_kenel(s = seq(0.1, 0.9, len = 5),
#'                                t = seq(0.2, 0.8, len = 5),
#'                                operator_norm = 0.5,
#'                                return_matrix = TRUE)
#' dim(mtx)
#'
#' # Return a data.table
#' dt <- get_real_data_far_kenel(s = seq(0.1, 0.9, len = 5),
#'                                t = seq(0.2, 0.8, len = 5),
#'                                operator_norm = 0.5,
#'                                return_matrix = FALSE)
#' head(dt)
#'
get_real_data_far_kenel <- function(s = seq(0.1, 0.9, len = 10),
                                    t = seq(0.2, 0.8, len = 10),
                                    operator_norm = 0.5,
                                    return_matrix = TRUE){
  # \eta(s)
  coss_mat <- outer(X = s, Y = 1:5, function(si, l) sqrt(2) * cos(2 * pi * l * si))
  sins_mat <- outer(X = s, Y = 1:5, function(si, l) sqrt(2) * sin(2 * pi * l * si))
  eta <- cbind(1, coss_mat, sins_mat)

  # \theta(t)
  cost_mat <- outer(X = t, Y = 1:5, function(ti, l) sqrt(2) * cos(2 * pi * l * ti))
  sint_mat <- outer(X = t, Y = 1:5, function(ti, l) sqrt(2) * sin(2 * pi * l * ti))
  colnames(cost_mat) <- paste0("cos", 1:5)
  colnames(sint_mat) <- paste0("sin", 1:5)
  theta <- cbind(1, cost_mat, sint_mat)

  # Basis function
  basis_fun <- kronecker(X = eta, Y = theta)

  # Basis coefficient
  basis_coef <- c(
    2.2337372971018, -0.325244044731151, 0.187206920919236, -0.251566711879163,
    -0.0167644005889917, 0.276667058786411, 0.0615073725778918, 0.571603018336332,
    -0.572809375110059, 0.0141379256684989, 0.0149022741092313, -4.82923610792512,
    0.331031200463392, 0.183596326343559, -0.027634868017173, 0.140334666930028, 0,
    -0.133119857316959, -0.0895419722599467, 0.00307262837174549, -0.0483508443594607,
    -0.126251329568179, 3.32335705082494, -0.0540202044263184, 0.297999174567218,
    -0.0715639682940563, 0.0286475855006393, 0.205067486705984, -0.131347879232633,
    0.353557986549793, -0.162035292128456, -0.039114894690579, -0.0459684159226797,
    1.7970025932255, 0, 0.0983572083794547, -0.00897146114351682, 0, 0.0654171627576935,
    0.13131934147443, -0.0092546291688908, -0.185062678568791, 0.0906326014811482,
    -0.0481786322361082, -4.72746806843479, 0.050670616460877, 0.0365619524515104, 0,
    0.284723029485189, 0.0126513369315148, 0, -0.0698838250498888, 0, 0, -0.125755965630984,
    -5.73918751880059, 0.245848419311172, 0, 0.267121966309844, 0, -0.0896598422716199,
    0.0316171382003794, -0.49768715975112, 0.52440949359224, 0, 0.197938756150832,
    2.45125518767279, -0.0463018677075481, -0.0549014918672375, -0.027309159416357,
    -0.0604572008598179, 0.0750983773087162, 0.136683725955063, 0.142855565715424,
    -0.297541175307152, -0.011750226937323, 0.027077398939407, 1.07367258651513,
    0.18513260615123, -0.0244962109056109, 0.0671308295876077, -0.0756788959583738,
    -0.0390603398330477, 0.191576162911029, -0.0250687004190137, 0.205775044698164,
    -0.0517597384979335, 0.00607347821426663, -5.43802938211584, -0.161785518980206,
    0.191707201891598, -0.00550062556645803, 0.0111654607750166, 0.255463915245325,
    0.0821500379029488, 0.445786769743474, -0.354321377284522, 0.0587158994893179,
    -0.0962640169091532, -3.51254535085072, -0.11445751095267, 0.12880151541983, 0,
    -0.118712184576829, 0, 0.245074852935059, 0.200974977296935, -0.3024432765705,
    0.289741051089142, 0, -0.983346082144421, 0.0820977879178276, -0.224174908365937,
    0.246518722931093, 0, -0.0538389230933072, -0.161386348786155, -0.00595281133102742,
    0.213730662849265, 0, 0.192629550171012
  )

  # operator norm
  op_norm <- 4.588783
  op_scale <- operator_norm / op_norm


  # FAR kernel
  far_ker <- basis_fun %*% matrix(data = basis_coef, ncol = 1)
  dt_far_ker <- data.table::as.data.table(expand.grid("s" = s, "t" = t))
  dt_far_ker <- dt_far_ker[order(t)]
  dt_far_ker[, ker_value := far_ker[, 1]]

  # Scale operator norm
  dt_far_ker[, ker_value := op_scale * ker_value]

  if (return_matrix) {
    dt_far_ker_dcast <- data.table::dcast(data = dt_far_ker, formula = s ~ t, value.var = "ker_value")
    mat_far_ker <- as.matrix(dt_far_ker_dcast[, -1])
    colnames(mat_far_ker) <- paste0("t", 1:length(s))
    rownames(mat_far_ker) <- paste0("s", 1:length(t))
    res <- mat_far_ker
  } else{
    res <- dt_far_ker
  }

  return(res)
}
