\name{rankEquivIntensity}
\alias{rankEquivIntensity}
\title{
Quantile Normalize
}
\description{
Perform Quantile Normalization on the columns of a matrix
}
\usage{
rankEquivIntensity(m, columnSets = list(all = 1:ncol(m)))
}
\arguments{
  \item{m}{
numeric matrix of expression data, one row for each gene/oligo and one column for
each sample
}
  \item{columnSets}{
list of vectors of column numbers, for doing the normalization by subsets of columns
}
}
\value{
numeric matrix of the same size as \code{m}, with normalized columns.
}

\seealso{
  \code{\link{duffyRMA}}, for RMA normalization.
}
