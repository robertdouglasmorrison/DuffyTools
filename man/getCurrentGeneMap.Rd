\name{getCurrentGeneMap}
\alias{getCurrentGeneMap}
\alias{getCurrentExonMap}
\alias{getCurrentSeqMap}
\alias{getCurrentRrnaMap}
\title{
Get the Current GeneMap
}
\description{
Get the annotation tables that define the current species
}
\usage{
getCurrentGeneMap()
getCurrentExonMap()
getCurrentSeqMap()
getCurrentRrnaMap()
}
\details{
See \code{\link{MapSets}} for details about the species annotation tables.
}
\value{
All these functions return one data frame, of the specified annotation
details for the current species.  
}
\note{
These tables must be treated as read-only copies of the annotation.  There are
no corresponding 'setCurrentMap' functions.  To alter any details, you must modify the
underlying MapSet.
}

\seealso{
  \code{\link{getCurrentSpecies}} {to get the SpeciesID of the current species.}
  \code{\link{setCurrentSpecies}} {to change to a different species.}
}
