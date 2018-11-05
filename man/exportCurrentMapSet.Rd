\name{exportCurrentMapSet}
\alias{exportCurrentMapSet}
\title{
Export the Current Annotation MapSet
}
\description{
Write the current organism MapSet out as a set of text files. 
}
\usage{
exportCurrentMapSet(path = ".")
}
\arguments{
  \item{path}{
the destination folder for the newly created text files 
}
}
\details{
This function exports the complete annotation for the one current
species, to a set of editable text files.  See \code{\link{MapSets}}
for a more detailed overview.
}
\value{
None.  Writes files of sequences, gene, exons, etc., to disk.
}

\seealso{
  \code{\link{importMapSet}}, for the inverse capability, to create a MapSet
  from a set of text files.
  \code{\link{getCurrentSpecies}}, to see the current species. 
  \code{\link{setCurrentSpecies}}, to change the current species (and MapSet).
}
