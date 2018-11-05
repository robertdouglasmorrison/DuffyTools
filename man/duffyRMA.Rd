\name{duffyRMA}
\alias{duffyRMA}
\alias{duffyRMA.bgSubtract}
\alias{duffyRMA.qn}
\alias{duffyMagnitudeNormalize}
\title{
Robust Multichip Average (RMA)
}
\description{
Our local implementation of RMA normalization
}
\usage{
duffyRMA(m, targetBGlevel = NULL, magnitudeScale = NULL, columnSets = list(all = 1:ncol(m)))

duffyMagnitudeNormalize(m, globalSum = 1000 * nrow(m))
duffyRMA.bgSubtract(m, targetBGlevel = NULL)
duffyRMA.qn(m, columnSets = list(all = 1:ncol(m)))
}
\arguments{
  \item{m}{
numeric matrix of expression data, with one row per oligo/gene and one column per sample
}
  \item{targetBGlevel}{
background correction target value.  \code{NULL} means use 'DensityKernelMode'
}
  \item{magnitudeScale}{
mean gene intensity target value (aka Quackenbush).  \code{NULL} means no global scaling step
}
  \item{columnSets}{
a list of vectors of column numbers, to do RMA on subsets of columns of \code{m}
}
  \item{globalSum}{
total sum of gene intensity target value (aka Quackenbush)
}
}
\details{
This implementation of RMA gives more layers of flexibility.  Global scaling in the style
of Quackenbush (2002), with settable average gene intensity, is the first step.
Next, the adjustment for background hybridization levels can be set to an explicit value,
or calculated from the mode of the density kernel.  Lastly, the quantile normalization
step is run.  Each of these functions is callable separately for finer control.
}
\value{
a matrix of the same size as \code{m}, with normalized expression values
}
\references{
Irizarry, B. & Bolstad, B.  (2003) Nucleic Acids Research 31(4):e15
}
\author{
Bob Morrison
}

\seealso{
  \code{\link{getGlobalMetrics}}, for various metrics of an expression matrix. 
  \code{\link{DKMshift}}, for details on applying the density kernel mode 
  background correction. 
  \code{\link{rankEquivalentIntensity}}, for simple quantile normalization.
}
