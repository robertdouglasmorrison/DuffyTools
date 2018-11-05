\name{multicore.setup}
\alias{multicore.setup}
\alias{multicore.samplePairs}
\alias{multicore.tapply}
\alias{multicore.lapply}
\title{
Multicore Processing Wrappers
}
\description{
A set of wrapper functions for multicore processing, using the
\code{multicore} package.
}
\usage{
multicore.setup(max.cores = NULL)

multicore.samplePairs(sampleSet, symmetric.pairs = FALSE)

multicore.tapply(x, INDEX, FUN, ...)
multicore.lapply(x, FUN, ...)
}
\arguments{
  \item{max.cores}{
integer specifying how many cores to request. \code{NULL} means use all
available
}
  \item{sampleSet}{
character vector (or list of character vectors), usually SampleIDs, that
need to be turned into a list of all pairs from \code{sampleSet}, for a
subsequent multicore call to a function that expects two arguments, such
as a differential expression operation.  When \code{sampleSet} is a list,
only pairs from within each element are generated.
}
  \item{symmetric.pairs}{
logical, are symmetric pairs needed.  This is for cases where the function to be
called is sensitive to the order of its 2 arguments (i.e. 
\code{ f(x1,x2) != f(x2,x1)}.  If true, more pairs are generated.
}
  \item{x, INDEX, FUN, ...}{
the usual arguments passed on to the respective 'apply' function
}
}
\details{
These are simple wrappers around the \code{multicore} routines, to facilitate
parcelling computations that involve entire chromosomes or samples out to
separate cores.  
}
\value{
For \code{multicore.setup}, the number of cores actually available.  The function
tries to load package \code{multicore}, and if successful, sets the number of
cores to be used.  This value can be subsequently returned by \code{getOption("cores")}.

For \code{multicore.samplePairs}, a list of vectors (each of length 2), suitable for
a call to \code{multicore.lapply} for a function expecting a vector of length 2 as its
argument.

For \code{multicore.tapply} and \code{multicore.lapply}, the same type result as their
primitive functions, after being gathered from the separate cores.
}
\note{
All multicore functions turn off prescheduling (\code{ mc.preschedule=FALSE}), because
we assume the processing time of each piece is highly variable.
}

\examples{
multicore.samplePairs( list( c("Tom", "Dick", "Harry"), c("Jack","Jill")), symmetric.pairs=TRUE)
}
