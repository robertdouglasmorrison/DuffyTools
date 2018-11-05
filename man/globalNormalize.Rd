\name{globalNormalize}
\alias{globalNormalize}
\title{
Normalize Microarray Intensity
}
\description{
A generic normalization function with several methods
}
\usage{
globalNormalize(m, 
       intNormMethod = c("sum", "median", "none", "logmean", "baseline", "dkm", "rank"), 
       metricSet = NULL, target = NULL)
}
\arguments{
  \item{m}{
numeric matrix of expression data, with a row for each gene/oligo and a column for each sample
}
  \item{intNormMethod}{
the type of normalization to perform, i.e. which metric of \code{m} will be used
}
  \item{metricSet}{
the set of global metrics to base normalization on.  \code{NULL} uses \code{\link{getGlobalMetrics}}
}
  \item{target}{
the target value to normalize to. \code{NULL} uses the corresponding value from 
\code{\link{getGlobalMetrics}}
}
}
\value{
a numeric matrix of the same size as \code{m}, after normalization of each column.
}
\author{
Bob Morrison
}

\seealso{
  \code{\link{ duffyRMA}}
}
