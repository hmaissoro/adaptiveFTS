% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/09_real_data_function.R
\name{get_real_data_far_kenel}
\alias{get_real_data_far_kenel}
\title{FAR kernel learned from the voltage curves of the electricity}
\usage{
get_real_data_far_kenel(s = 0.2, t = 0.3, operator_norm = 0.5)
}
\arguments{
\item{s}{\code{numeric (positive)}. A vector or scalar value(s) between 0 and 1.}

\item{t}{\code{numeric (positive)}. A vector or scalar value(s) between 0 and 1.}

\item{operator_norm}{\code{numeric (positive)}. A scalar corresponding to the norm of the integral operator associated with this kernel function.}
}
\value{
A vector (or scalar) of \code{numeric} values corresponding to the value of the kernel function evaluated at (\code{s}, \code{t}).
}
\description{
For more details see the vignette:
\code{vignette("hybrid-simulation-setup", package = "adaptiveFTS")}
}
\examples{

# get the value of the kernel at (s,t) = (0.2, 0.3)
kerval <- get_real_data_far_kenel(s = 0.2, t = 0.3, operator_norm = 0.5)
kerval

}
