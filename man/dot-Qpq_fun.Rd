% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_estimate_autocovariance.R
\name{.Qpq_fun}
\alias{.Qpq_fun}
\title{\eqn{Q_{pq}^{(\ell)}}, (\eqn{\ell \leq 0}) function. See \insertCite{rubin2020;textual}{adaptiveFTS} Equation (B.7)}
\usage{
.Qpq_fun(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = 1/4,
  t = 1/2,
  lag = 1,
  p = 1,
  q = 1,
  h,
  dt_mean_rp = NULL,
  optbw_mean = NULL,
  smooth_ker = epanechnikov
)
}
\arguments{
\item{data}{A \code{data.table} (or \code{data.frame}), a \code{list} of \code{data.table} (or \code{data.frame}), or a \code{list} of \code{list}.
\itemize{
\item{If \code{data.table}:}{
It should contain the raw curve observations in at least three columns.
\itemize{
\item{\code{idcol} :}{ The name of the column containing the curve index in the sample.
Each curve index is repeated according to the number of observation points.}
\item{\code{tcol} :}{ The name of the column with observation points associated with each curve index.}
\item{\code{ycol} :}{ The name of the column with observed values at each observation point for each curve index.}
}
}
\item{If \code{list} of \code{data.table}:}{
In this case, each element in the \code{list} represents the observation data of a curve in the form of a \code{data.table} or \code{data.frame}.
Each \code{data.table} contains at least two columns.
\itemize{
\item{\code{tcol} :}{ The name of the column with observation points for the curve.}
\item{\code{ycol} :}{ The name of the column with observed values for the curve.}
}
}
\item{If \code{list} of \code{list}:}{
In this case, \code{data} is a list where each element is the observation data of a curve, given as a \code{list} of two vectors.
\itemize{
\item{\code{tcol} :}{ The vector containing observation points for the curve.}
\item{\code{ycol} :}{ The vector containing observed values for the curve.}
}
}
}}

\item{idcol}{\code{character}. If \code{data} is given as a \code{data.table} or \code{data.frame}, this is the name of the column that holds the curve index.
Each curve index is repeated according to the number of observation points. If \code{data} is a \code{list} of \code{data.table} (or \code{data.frame}) or a \code{list} of \code{list}, set \code{idcol = NULL}.}

\item{tcol}{\code{character}. The name of the column (or vector) containing the observation points for the curves.}

\item{ycol}{\code{character}. The name of the column with observed values for the curves.}

\item{s}{\code{vector (numeric)}. First argument of the autocovariance function.
It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
It has to be of the same length as the \code{t}}

\item{t}{\code{vector (numeric)}. Second argument of the autocovariance function.
It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
It has to be of the same length as the \code{s}.}

\item{lag}{\code{integer (positive integer)}. Lag of the autocovariance.}

\item{p}{\code{numeric (integer)}. It is used as exponent.}

\item{q}{\code{numeric (integer)}. It is used as exponent.}

\item{h}{\code{numeric (positive scalar)}. The bandwidth of the estimator.}

\item{dt_mean_rp}{\code{data.table}. It contains the estimates of the mean function at each observation point for each curve.
The name of the curve identification column must be \code{id_curve}, the observation points column \code{tobs} and the mean estimates column \code{muhat_RP}.
Default \code{dt_mean_rp = NULL} and so it will be estimated.}

\item{optbw_mean}{\code{numeric (positive scalar)}. Optimal bandwidth for the mean function estimator.
It is \code{NULL} if \code{dt_mean_rp} is not \code{NULL}.}

\item{smooth_ker}{\code{function}. The kernel function of the Nadaraya-Watson estimator.
Default \code{smooth_ker = epanechnikov}.}
}
\value{
A \code{numeric} scalar.
}
\description{
\eqn{Q_{pq}^{(\ell)}}, (\eqn{\ell \leq 0}) function. See \insertCite{rubin2020;textual}{adaptiveFTS} Equation (B.7)
}
\references{
\insertAllCited{}
}
