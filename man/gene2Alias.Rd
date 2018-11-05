\name{gene2Alias}
\alias{gene2Alias}
\title{
  Convert from a GeneID to a set of known aliases.
}
\description{
  Looks up known aliases from official GeneIDs.
}
\usage{
gene2Alias(genes, coreSubsetOnly = TRUE)
}
\arguments{
  \item{genes}{
a character vector of GeneIDs, valid for the current SpeciesID.
}
  \item{coreSubsetOnly}{
logical, should only the core subset of most common aliases be returned, or the complete
set of every known alias.
}
}
\details{
Implements the inverse of the \code{\link{alias2Gene}} function.
}
\value{
a character vector, of the same length as \code{genes}, with all found aliases for each gene
concatenated into a single comma separated string.
}
\author{
Bob Morrison
}
\seealso{
  See Also \code{\link{alias2Gene}}
}
\examples{
gene2Alias( "PF3D7_1200600")
}
