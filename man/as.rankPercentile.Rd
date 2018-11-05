\name{as.rankPercentile}
\alias{as.rankPercentile}
\title{
Convert to Rank Percentiles
}
\description{
Convert a vector of numbers into rank percentiles.
}
\usage{
as.rankPercentile(x)
}
\arguments{
  \item{x}{
numeric vector, usually of gene expression abundance
}
}
\details{
Converts gene expression onto a uniform 0..100 scale, by expressing all
values as rank percentiles.  The largest value will become 100, the median
value 50, and the minimun value will be \code{100/length(x)}
}
\value{
numeric vector of the same length as \code{x}, where each value is the
rank percentile of the original.
}
\examples{
as.rankPercentile( rnorm( 10, mean=5, sd=2))
}
