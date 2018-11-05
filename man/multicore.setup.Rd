\name{multicore.setup}
\alias{multicore.setup}
\alias{multicore.currentCoreCount}
\alias{multicore.totalCoreCount}
\alias{multicore.samplePairs}
\alias{multicore.tapply}
\alias{multicore.lapply}
\alias{multicore.by}
\title{
Multicore Processing Wrappers
}
\description{
A set of wrapper functions for multicore processing, using the
\code{parallel} package.
}
\usage{
multicore.setup(max.cores = NULL)

multicore.currentCoreCount()
multicore.totalCoreCount()

multicore.samplePairs(sampleSet, symmetric.pairs = FALSE)

multicore.tapply(x, INDEX, FUN, ...)
multicore.lapply(x, FUN, ...)
multicore.by(data, INDICES, FUN, ...)
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
  \item{data, INDICES, FUN, ...}{
the usual arguments passed on to the 'by' function
}
}
\details{
These are simple wrappers around the \code{parallel} routines, to facilitate
parcelling computations that involve entire chromosomes or samples out to
separate cores.  
}
\value{
For \code{multicore.setup}, the maximum of the number of cores actually available
or the number of cores requested.  The function
tries to load package \code{parallel}, and if successful, sets the number of
cores to be used.  This value can be subsequently returned by \code{multicore.currentCoreCount}.

For \code{multicore.currentCoreCount}, integer telling the number of cores currently
requested via \code{multicore.setup}.

For \code{multicore.totalCoreCount}, integer telling the number of cores actually
present on this compute node. If \code{multicore.setup} has not yet been called, it
does so and requests all cores available.

For \code{multicore.samplePairs}, a list of vectors (each of length 2), suitable for
a call to \code{multicore.lapply} for a function expecting a vector of length 2 as its
argument.

For \code{multicore.tapply} and \code{multicore.lapply}, the same type result as their
primitive functions, after being gathered from the separate cores.
}
\note{
All multicore functions turn off prescheduling (\code{ mc.preschedule=FALSE}), because
we assume the processing time of each piece is highly variable.

All multicore functions turn off recursion (\code{ mc.allow.recursive=FALSE}), because
we assume all available cores will be used by the parent call.
}

\examples{
multicore.samplePairs( list( c("Tom", "Dick", "Harry"), c("Jack","Jill")), symmetric.pairs=TRUE)
}
