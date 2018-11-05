\name{DKMshift}
\alias{DKMshift}
\title{
Density Kernel Mode Linear Shift
}
\description{
Low level routine to apply the DKM background correction step of RMA
}
\usage{
DKMshift(m, target = getGlobalMetrics(obj)$DKM)
}
\arguments{
  \item{m}{
numeric matrix of expression data, with a row for each gene/oligo, and a
column for each sample
}
  \item{target}{
target value for the mode of the density kernel for each column
}
}
\details{
This function implements the background subtraction step of RMA.  It
prevents the creation of negative expression values, by a linear
compaction of values smaller than \code{target}.
}
\value{
a numeric matrix of the same size as \code{m}, with the background
correction linear shift applied to all columns 
}

\seealso{
  \code{\link{duffyRMA}}, for details about the RMA implementation.
  \code{\link{getDensityKernelMode}}, for calculating DKM.
  \code{\link{getGlobalMetrics}}, for various metrics of an expression matrix
including the mean DKM value. 
}
