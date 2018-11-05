\name{calcRP}
\alias{calcRP}
\title{
low level Rank Product function
}
\description{
Calculates the Rank Product from a set of gene rank probabilities
}
\usage{
calcRP(x)
}
\arguments{
  \item{x}{
numeric vector of gene rank probabilities
}
}
\details{
The low level calculation of the RP value, as described in the
Rank Product paper by Rainer Breitling.    Given an ordering of
N genes in a differential expression result, the probability of
a gene occuring at position j is given by (j/N).  With a set of
K replicates, Breitling defines the RP value for a gene as the
K'th root of the product of the probabilities from all K 
replicates.
}
\value{
numeric value of length 1, equivalent to the geometric mean of \code{x}
}
\references{
Breitling, R.  FEBS Letters 573 (2004)
}
\note{
Implemented via \code{log2(x)} for numerical stability, and further
normalized by the number of replicates K to allow comparisons 
between datasets.

As implemented, it is equivalent to \code{\link{logmean}}.  Contrast with
\code{\link{mean}} and \code{\link{sqrtmean}}.
}

\seealso{
  \code{\link{rankProduct}}, \code{\link{rankProductDiffExpress}}, for
  the higher level routines based on Rank Product
}
