\name{getDensityKernelMode}
\alias{getDensityKernelMode}
\title{
Mode of the Density Kernel
}
\description{
Calculate the mode of the density kernel function
}
\usage{
getDensityKernelMode(v)
}
\arguments{
  \item{v}{
numeric vector of values, typically gene or oligo expression abundance
}
}
\details{
This function is used in RMA normalization to estimate the amount
of non-specific hybridization on a microarray.
}
\value{
a numeric value on length 1, representing the most frequent value
in the distribution of \code{v}.  
}
\note{
Uses the 'epanechnikov' kernel.  See \code{\link{density}}.
}
\seealso{
  \code{\link{duffyRMA}}, \code{\link{DKMshift}}
} 
\examples{
  getDensityKernelMode( rnorm( 100, mean=10, sd=2))
}

