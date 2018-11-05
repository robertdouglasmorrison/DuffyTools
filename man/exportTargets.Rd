\name{exportTargets}
\alias{exportTargets}
\alias{importTargets}
\title{
Export or Import all Targets
}
\description{
Write the set of all known organism targets to a text file, or
reload the set of known organisms from a text file.
}
\usage{
exportTargets(fileout = "DuffyTools.Targets.txt")
importTargets(filein = "DuffyTools.Targets.txt")
}
\arguments{
  \item{file}{
filename for the text file of target details
}
}
\details{
This function exports the table of currently defined organism targets to
an editable text file, or does the reverse to import targets from a text file.
  See \code{\link{Targets}} for a more
detailed overview.
}
\value{
For 'export', writes one file of defined targets.
For 'import', reads the file and resets the family of known targets
}

\seealso{
  \code{\link{getCurrentTarget}}, to get the current targetID. 
  \code{\link{setCurrentTarget}}, to change the current target. 
  \code{\link{getCurrentTargetSpecies}}, to get the set of speciesIDs for this target. 
  \code{\link{getCurrentTargetFilePrefix}}, to get the set of file prefixes for this target. 
}
