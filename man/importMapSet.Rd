\name{importMapSet}
\alias{importMapSet}
\title{
Import Annotation MapSet Files
}
\description{
Create a MapSet object from a set of text files. 
}
\usage{
importMapSet(path = ".")
}
\arguments{
  \item{path}{
the folder that holds the annotation map text files.
}
}
\details{
This function rebuilds an annotation MapSet from a set of map text files.
See \code{\link{MapSets}} for a more detailed overview.
}
\value{
A MapSet object (returned invisibly).  Also, a compressed MapSet.rda file 
gets written to \code{path}, and the MapSet gets added to the set of defined species.
}

\seealso{
  \code{\link{addMapSet}}, for adding new MapSets to the set of defined species.
  \code{\link{exportCurrentMapSet}}, for the inverse capability, to write a MapSet to disk.
}
