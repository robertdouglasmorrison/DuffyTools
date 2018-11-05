\name{gaussian}
\alias{gaussian}
\alias{fit.gaussian}
\title{
Generate and/or Fit Gaussian Data
}
\description{
Generate a vector of values that represent a Gaussian distribution.
Secondly, apply a nonlinear least squares fit to gaussian data to
determine the fitted gaussian function.
}
\usage{
gaussian(x, mean = 0, sd = 1, height = NULL, offset = 0)

fit.gaussian( x, y, fit.offset = TRUE)
}
\arguments{
  \item{x}{
vector of X locations to evaluate the gaussian at.  For 'fit.gaussian',
the locations at which the measured Y values are taken from.
}
  \item{mean}{
the center locatinon X for the generated gaussian curve
}
  \item{sd}{
the standard deviation for the generated gaussian curve
}
  \item{height}{
the maximum height of the generated gaussian curve.  By default, the curve
will have height such that its integral equals 1.
}
  \item{offset}{
the linear offset for the baseline of the generated gaussian curve
}
  \item{y}{
the raw data measurements at 'x', that are to be fit by a gaussian model.
}
  \item{fit.offset}{
logical, should the model include the offset term, to fit the baseline tails
of the raw data.
}
}
\value{
For 'gaussian', a vector of length \code{ length(x)}, that gives the amplitude
of that gaussian function evaluated at \code{x}.

For \code{fit.gaussian}, a list:
  \item{coefficients }{ a vector of named coefficients (mean, sd, height, offset)
  giving the parameters of the best fit gaussian model. }
  \item{y }{ a vector of the same length as \code{y}, giving the fitted result
  of the gaussian model. }
}
\references{
  based on code from:  Earl F. Glynn,  Stowers Inst for Medical Research
}
\author{
Bob Morrison
}
\note{
Implemented via \code{ nls}
}

