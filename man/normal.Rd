\name{normal}
\alias{normal}
\alias{fit.normal}
\alias{gaussian}
\alias{fit.gaussian}
\alias{lorentzian}
\alias{fit.lorentzian}
\alias{gumbel}
\alias{fit.gumbel}
\title{
Functions for generating and fitting named distribution functions.
}
\description{
Generate or fit by nonlinear least squares a family of classic distribution functions.
}
\usage{
normal(x, mean = 0, sd = 1, height = NULL, floor = 0)
fit.normal(x, y, start.mean = 0, start.sd = 1, start.height = NULL, start.floor = 0)

gaussian(x, center = 0, width = 1, height = NULL, floor = 0)
fit.gaussian(x, y, start.center = 0, start.width = 1, start.height = NULL, start.floor = 0)

lorentzian(x, center = 0, width = 1, height = NULL, floor = 0)
fit.lorentzian(x, y, start.center = 0, start.width = 1, start.height = NULL, start.floor = 0)

gumbel(x, center = 0, width = 1, height = NULL, floor = 0)
fit.gumbel(x, y, start.center = 0, start.width = 1, start.height = NULL, start.floor = 0)
}
\arguments{
  \item{x}{
a vector of values to evaluate or fit the function at
}
  \item{y}{
a vector of observed Y values to fit the named distribution to
}
  \item{mean, center}{
the central value for the function
}
  \item{sd, width}{
the nominal measure of the width of the distribution.  Note that this is a signed value
for gumbel distributions, affecting the direction of the asymmetric tail.
}
  \item{height}{
an optional multiplier for adjusting the magnitude of the Y values returned.  By default, the height is defined 
by the underlying function.
}
  \item{floor}{
an optional floor value for the tails of the distribution.  This has the effect of applying a linear offset
to the Y values.
}
  \item{start.center, start.width, start.height, etc.}{
for the \code{fit.} functions, optional starting estimates passed to the NLS routine.  See \code{\link{nls}}.
to the Y values.
}
}
\value{
For the curve generation functions, a vector of Y values, from evaluating the function at all values in X.

For the curve fitting functions, a list containing:
  \item{y }{ a vector of best fit values of Y, at each location in X.}
  \item{mean,center }{ the best fit curve parameter of the distribution's central value.}
  \item{sd,width }{ the best fit curve parameter of the distribution's width.}
  \item{height }{ the best fit curve parameter of the distribution's height.}
  \item{floor }{ if fitted, the best fit curve parameter of the distribution's floor 
  	(linear offset of the function tails)}
}
\author{
Bob Morrison
}
\note{
While the functions are implemented differently, 'normal' an 'gaussian' are effectively the same.
}
