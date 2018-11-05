\name{logmean}
\alias{logmean}
\alias{logMeanPlusSD}
\alias{logMedianPlusSD}
\title{
Log-Average or Geometric Mean
}
\description{
Compute the log-average or geometric mean of a numeric vector
}
\usage{
logmean(x, na.rm = FALSE)

logMeanPlusSD(x, na.rm = FALSE)
logMedianPlusSD(x, na.rm = FALSE)
}
\arguments{
  \item{x}{
numeric vector of positive values
}
  \item{nSD}{
numeric scalar, the number of standard deviations to add
}
}
\details{
Implements the geometric mean, by log transforming the values in \code{x}, calculating
that mean, and then exponentiating the result.  This has better numerical stability then
then using the N-th root of the power series of \code{x}.

The 'PlusSD' forms are used to thresholds of statistical significance, by adding an 
arbitrary number of standard deviation of \code{log(x)}, to the calculation.

As log transform is undefined for zeros and negative values, some level of data cleaning
is implemented to tolerate invalid data.  Specifically, negatives are set to NA and thus
ignored.  If all values are zero, zero is returned.  If a subset of values are zero, they are 
replaced with the smallest non-zero value prior to log transform and averaging.
}
\value{
a numeric scalar
}
\seealso{
Contrast with \code{\link{mean}} and \code{\link{sqrtmean}}.
}
