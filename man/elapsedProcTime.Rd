\name{elapsedProcTime}
\alias{elapsedProcTime}
\title{
Provide processing timing information
}
\description{
A utility to describe elapsed processing time in a 'units centric' way.
}
\usage{
elapsedProcTime(t1, t2, N = 1, what = "Read")
}
\arguments{
  \item{t1, t2}{
the start and finish time objects, as returned by \code{\link{proc.time}}
}
  \item{N}{
integer specifying how many \code{what} have been processed in the 
time interval from \code{t1} to \code{t2}
}
  \item{what}{
character string label for the type of item processed.  Used to label the output
}
}
\value{
A list of length 4, giving the elasped time in various ways:
\item{Total.Elapsed.Time}{ Elapsed time.  See \code{\link{proc.time}} }
\item{Computation.Time }{ Total processor time.  Can be greater than elapsed time when using more than one processor.}
\item{Time.per.'what'}{ Total processor time divided by 'N'}
\item{Raw.Seconds}{ Elapsed 'Wall time'}
}
