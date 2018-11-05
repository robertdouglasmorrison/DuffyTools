\name{getCurrentTarget}
\alias{getAllTargets}
\alias{getCurrentTarget}
\alias{getCurrentTargetSpecies}
\alias{getCurrentTargetFilePrefix}
\title{
Get Current Target
}
\description{
Get details about the current target of organism(s)
}
\usage{
getAllTargets()
getCurrentTarget()
getCurrentTargetSpecies()
getCurrentTargetFilePrefix()
}
\details{
These functions retrieve details about the currently defined target organism(s).
See \code{\link{Targets}} for a more detailed overview.
}
\value{
For 'getAllTargets', a table of TargetIDs, SpeciesIDs, FilePrefixes, and Species Text annotation details.
For 'getCurrentTarget', a table with one row, of TargetID, SpeciesIDs, FilePrefixes, and Species Text annotation details.
For 'getCurrentTargetSpecies', a vector of this target's SpeciesIDs.
For 'getCurrentTargetFilePrefix', a vector of this target's file prefixes.
}

\seealso{
  \code{\link{exportTargets}}, for writing the targets table to a disk file.
  \code{\link{setCurrentTarget}}, to change the current target. 
}
