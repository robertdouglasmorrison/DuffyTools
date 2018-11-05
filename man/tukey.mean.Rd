\name{tukey.mean}
\alias{tukey.mean}
\alias{tukey.biweight}
\alias{tukey.logmean}
\title{
Averaging via Tukey's Biweight
}
\description{
A family of robust estimators based on Tukey's biweight
}
\usage{
tukey.biweight(x, c = 5, e = 1e-04, na.rm = FALSE)
tukey.mean(x, c = 5, e = 1e-04, na.rm = FALSE)
tukey.logmean(x, c = 5, e = 1e-04, na.rm = FALSE)
}
\arguments{
  \item{x}{
numeric vector of values to find the 'average' of
}
  \item{c,e}{
arguments to the underlying biweight algorithm
}
}
\details{
Tukey's Biweight calculates a set of weight terms for the elements of
\code{x}, based on a function of the each element of \code{x} distance from
the median of \code{x}.  Those weight are then used in the calculation of the
average of \code{x}, by either mean or log-average (geometric mean).
}
\value{
For 'biweight', a numeric vector of length \code{length(x)}, of the scalar
weights for each element of \code{x}.

Otherwise, a numeric scalar that is the 'average' of \code{x}.
}
