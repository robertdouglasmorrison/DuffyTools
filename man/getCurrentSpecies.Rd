\name{getCurrentSpecies}
\alias{getAllSpecies}
\alias{getAllSpeciesFilePrefixes}
\alias{getCurrentSpecies}
\alias{getCurrentSpeciesFilePrefix}
\alias{getSpeciesText}
\title{
Get Current Species
}
\description{
Get details about the current species
}
\usage{
getAllSpecies()
getAllSpeciesFilePrefixes()
getCurrentSpecies()
getCurrentSpeciesFilePrefix()
getSpeciesText( speciesID)
}
\details{
These functions retrieve details about the currently active species.  Keep in
mind that a target can have more than one species, but only a single species
is actually current at one time.  See \code{\link{Species}} for a more detailed 
overview.  See \code{\link{MapSets}} for details on species annotations in 
general.
}
\value{
For 'getAllSpecies', a character vector of SpeciesIDs.
For 'getAllSpeciesFilePrefixes', a named character vector of file prefixes, where
the names are the corresponding SpeciesIDs.
For 'getCurrentSpecies', a single SpeciesID for the currently active species.
For 'getCurrentSpeciesFilePrefix', a single character string of this species's file prefix.
For 'getSpeciesText', a text string description of that species annotation details.
}

\seealso{
  \code{\link{setCurrentSpecies}} {to change the current species.}
  \code{\link{getAllTargets}} {to see the full list of all targets and species currently defined.}
}
