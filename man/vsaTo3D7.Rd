\name{vsaTo3D7}
\alias{vsaTo3D7}
\alias{vsaGeneProduct}
\title{
Variable Surface Antigen Genes
}
\description{
Specialized Plasmodium falciparum gene functions for variable surface antigens.
}
\usage{
vsaTo3D7(genes)

vsaGeneProduct(genes)
}
\arguments{
  \item{genes}{
character vector of VSA GeneIDs
}
}
\details{
These functions return details that are specific to the families of
variable surface antigens in the analysis of P.falciparum datasets.
}
\value{
For \code{vsaTo3D7}, the closest GeneIDs in the 3D7 annotation, base
on BLAST similarity of the protein sequences.

For \code{vsaGeneProduct}, a gene product description that encapsulates
the closest 3D7 GeneID and var gene grouping.
}
\references{
Based on the 'VarDom' database:   http://www.cbs.dtu.dk/services/VarDom/
}
\examples{
  vsaGeneProduct( c( "DD2var06", "IT4var04"))
}
