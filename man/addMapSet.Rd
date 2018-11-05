\name{addMapSet}
\alias{addMapSet}
\title{
Add Organism Annotation MapSet to DuffyTools environment
}
\description{
Loads an additional organism annotation to the set of knowm organisms.
}
\usage{
addMapSet(mapset)
}
\arguments{
  \item{mapset}{
a MapSet object, as created by 'importMapSet'
}
}
\details{
A MapSet is the complete description of an organism's annotation, covering
chromosomes, gene, and exons.  Adding a MapSet for an already loaded 
organism overwrites the previous annotation.

See \code{\link{MapSets}} for an overview of organism annotations.
}
\value{
If successfully loaded, the speciesID of this new organism.
}

\seealso{
  \code{\link{addTarget}}, for combining multiple organisms 
  \code{\link{exportCurrentMapSet}}, for saving an annotation to text files 
  \code{\link{importMapSet}}, for bundling annotation text files back into a MapSet object
}
