\name{sqrtmean}
\alias{sqrtmean}
\title{
Generalized Mean for P=0.5
}
\description{
Compute the square of the mean of the square root of X
}
\usage{
sqrtmean(x, na.rm = FALSE)
}
\arguments{
  \item{x}{
numeric vector of positive values
}
}
\details{
Implements the generalized mean, for the case of power P=0.5.
This average will always fall between the value of the arithmatic mean \code{\link{mean}}
and the value of the geometric mean \code{\link{logmean}}.

In trying to average the gene rank positions over multiple transcripts or 
differential expression results, the various averaging tools can give quite
different results, depending on the distribution of values.  Many functions
that produce sorted lists of genes provide a \code{average.FUN} argument, so
you can control the notion of 'average'.

}
\value{
a numeric scalar
}
