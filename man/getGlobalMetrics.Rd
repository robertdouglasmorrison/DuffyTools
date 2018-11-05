\name{getGlobalMetrics}
\alias{getGlobalMetrics}
\title{
Global Metrics of Expression Data
}
\description{
Get various statistical measurement metrics about a matrix of expression data,
typically of microarray oligo/gene intensities
}
\usage{
getGlobalMetrics(m)
}
\arguments{
  \item{m}{
a numeric matrix of expression values, with a row for each gene/oligo and a column
for each sample
}
}
\details{
This function measures several metrics about the data, all of which may be
useful for the various types of data normalization common to microarrays.
}
\value{
A list of numeric scalars.  For each scalar, the mean of that metric measured
separately for each sample, is the value reported.
  \item{mean }{ the average mean of each column of \code{m}}
  \item{median }{ the average median of each column of \code{m}}
  \item{floor }{ the average floor of each column of \code{m}}
  \item{max }{ the average max of each column of \code{m}}
  \item{min }{ the average min of each column of \code{m}}
  \item{logmean }{ the average logmean (geometric mean) of each column of \code{m}. 
  See \code{\link{logmean}}}
  \item{baseline }{ the average baseline of each column of \code{m}.  This is
  measured using intergenic and negative control probes when available, and 
  by DKM if not.}
  \item{DKM }{ the average DKM (Density Kernel Mode) of each column of \code{m}.
  See \code{\link{getDensityKernelMode}}}
}

\seealso{
  \code{\link{duffyRMA}}
}
