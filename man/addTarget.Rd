\name{addTarget}
\alias{addTarget}
\title{
Add a Target Organism(s) to DuffyTools Environment
}
\description{
Defines a new Target as a set of one or more organisms
}
\usage{
addTarget(targetID, speciesSet = targetID, 
	  prefixSet = speciesSet, mapset.path = NULL)
}
\arguments{
  \item{targetID}{
character string, giving the name for this new target
}
  \item{speciesSet}{
character string, giving the comma separated list of speciesIDs
}
  \item{prefixSet}{
character string, giving the comma separated list of filename prefixes
}
  \item{mapset.path}{
folder name for search location for user-defined MapSets
}
}
\details{
A target defines the set of organisms to be used for subsequent analyses.  The list of
speciesID(s) and prefix(s) tell the program which MapSet(s) and organism-specific filenames
to use.
}
\value{
If successfully loaded, the name of the target as a character string
}

\seealso{
 \code{\link{MapSets}}, for an overview of organism annotation
 \code{\link{getCurrentTarget}}, \code{\link{getCurrentTargetSpecies}}, \code{\link{getCurrentTargetFilePrefix}}, for details about the current target
 \code{\link{getAllTargets}}, for the table of all currently known targets
}
\examples{
	addTarget( "HsPf", "Pf3D7,Hs_grc", "Pf,Hs")
}
