% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{ppareto1}
\alias{ppareto1}
\title{Cumulative distribution function of the Pareto 1 distribution}
\usage{
ppareto1(x, mu = 1, alpha = 1/xi, xi = 1/alpha)
}
\arguments{
\item{x}{a (positive) vector}

\item{mu}{a number (the lower bound, default is 1)}

\item{alpha}{a number (the power index)}

\item{xi}{a number (the inverse of \code{alpha})}
}
\value{
the c.d.f. of the Pareto 1 distribution at points \code{x}
}
\description{
Cumulative distribution function of the Pareto 1 distribution
}
\examples{
ppareto1(2, 1, 1.5)
}
