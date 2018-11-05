\name{alias2Gene}
\alias{alias2Gene}
\alias{aliasGREP}
\title{
Convert gene aliases to standardized gene identifiers.
}
\description{
Maps (potentially) obsolete gene names and nicknames to the current standard annotation.
}
\usage{
alias2Gene(genes)
aliasGREP(genes)
}
\arguments{
  \item{genes}{
A character vector of gene names and/or nicknames.  In the case of \code{aliasGREP}, a vector
of character patterns to search for.
}
}
\details{
Uses the file "data/<speciesFilePrefix>.AliasTable.rda" from the DuffyTools installation
folder as the reference source of alias details.
Gene names that are not known or already standard nomenclature are returned unchanged.  
Standard names for Plasmodium species are from PlasmoDB.
}
\value{
A vector the same length as \code{genes}, with the most current standard gene 
name for each gene.
}
\note{
It is expected (but not enforced) that all geneID entries in the alias table
exactly match the GeneIDs in the GeneMap for that species.  The alias table used
is determined by the current species.
}
\seealso{
  \code{\link{gene2Alias}}
}
\examples{
  alias2Gene( c( "pfs25", "var2csa"))
}
