\name{setCurrentSpecies}
\alias{setCurrentSpecies}
\title{
Change the Current Species
}
\description{
Change to a different species (organism) for subsequent analysis,
thereby changing all annotation maps to this species.
}
\usage{
setCurrentSpecies(speciesID = NULL, seqID = NULL)
}
\arguments{
  \item{speciesID}{
specify the new species by SpeciesID
}
  \item{seqID}{
-or-  specify the new species by any of its SeqIDs
}
}
\details{
This is the only method to change annotations between the
various species that the Duffy Lab packages operate on.  See \code{\link{MapSets}}
for full details.
}
\value{
Returns the SpeciesID of the now-current species as a character string.  The 
main effect of the function is to change all annotation maps, so subsequent calls
to the \code{getCurrent...Map()} family of functions will return annotation
details about this now current species.
}

\seealso{
  \code{\link{setCurrentTarget}}, to change the set of organisms under study.
}
\examples{
  setCurrentSpecies( "Hs_grc")
  head( getCurrentGeneMap())
  setCurrentSpecies( "Pf3D7")
}
