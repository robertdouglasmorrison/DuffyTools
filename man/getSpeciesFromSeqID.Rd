\name{getSpeciesFromSeqID}
\alias{getSpeciesFromSeqID}
\alias{getOtherSpeciesFilePrefix}
\title{
Species Identifier Converters
}
\description{
Helper functions to convert between species specific indentifiers
}
\usage{
getSpeciesFromSeqID(seqID)
getOtherSpeciesFilePrefix(speciesID)
}
\arguments{
  \item{seqID}{
character vector of SeqIDs
}
  \item{speciesID}{
character vector of SpeciesIDs
}
}
\value{
For \code{getSpeciesFromSeqID}, a character vector of SpeciesID, of the same length 
as \code{seqID}.  Invalid or not found SeqIDs will be \code{NA}.  Only species defined
by the current target are searched.

For \code{getOtherSpeciesFilePrefix}, a character vector of file prefixes.  This works
without regard for which species is current.
}
