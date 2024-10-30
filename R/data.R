#' Sample for Functional Autoregressive Process of Order 1
#'
#' This dataset contains a set of 150 curves, each with an average of 90 observation points and some added noise.
#'
#' A \code{data.table} with three columns:
#'  \describe{
#'    \item{\code{id_curve} :}{ The index of each curve.}
#'    \item{\code{tobs} :}{ The observation points for each \code{id_curve}.}
#'    \item{\code{X} :}{ The observed values of the curve at each \code{tobs} point.}
#'  }
#'
"data_far"
