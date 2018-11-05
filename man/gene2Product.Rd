\name{gene2Product}
\alias{gene2Product}
\alias{gene2ProductAllSpecies}
\title{
Gene Product Description
}
\description{
Lookup the description of a gene product from the GeneID
}
\usage{
gene2Product(gNames)

gene2ProductAllSpecies(gNames, hints = NULL)
}
\arguments{
  \item{gNames}{
character vector of GeneIDs
}
  \item{hints}{
character vector of gene product descriptions from a secondary source
}
}
\details{
Finds and extracts the PRODUCT term from the appropriate GeneMap row, for
each GeneID.  
}
\value{
character vector of length \code{length(gNames)}, containing the gene product
descriptions.  Any GeneIDs that are not found return '""'.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
\code{gene2Product} is fastest, but only searches the current GeneMap (i.e. one
species only).  \code{gene2ProductAllSpecies} is slower, searching all MapSets 
of species in the current target.
}

\examples{
  setCurrentSpecies( "Pf3D7")
  gene2Product( c("PFA0620c", "PFL0030c"))
}
