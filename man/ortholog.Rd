\name{ortholog}
\alias{ortholog}
\alias{orthologSet}
\title{
Gene Orthologs
}
\description{
Convert GeneIDs between organisms
}
\usage{
ortholog(genes, from = "PF3D7", to = "PY17X")
orthologSet(genes, from = "PF3D7", to = "PY17X")
allOrthologs(genes)
orthologTransform( x, from, to = "Pf3D7", geneColumn = "GENE_ID",
		intensityColumn = "RPKM_M", productColumn="PRODUCT",
		symmetricOnly = FALSE, sep = "\t")
}
\arguments{
  \item{genes}{
Character vector of GeneIDs, from the species specified by \code{from}
}
  \item{from}{
Name of the current species
}
  \item{to}{
Name of the alternate species, to convert the GeneIDs to
}
  \item{x}{ a filename or a data.frame
}
}
\details{
Uses a data frame object of curated GeneIDs.  Orthologs are currated from PlasmoDB,
OrthoMCL, and other sources.
}
\value{
For \code{ortholog()}, a character vector of the same length as \code{genes}, holding the one best 
GeneID in the other species for each given gene.  Genes with no ortholog, and genes not found in 
the \code{from} species, are returned as empty strings.

For \code{allOrtholog()}, a list of the same length as \code{genes}, holding a vector of
all GeneIDs in all other species, for each given gene.  Genes with no ortholog in 
any other species are returned as empty vectors.
 
For \code{orthologSet()}, a character vector of all found orthologs to any of the given genes.  Note that 
since orthologs are not in general a one-to-one matching, the returned vector may be longer or shorter than
the input vector of genes.

In general, the set of returned orthologs are not guaranteed to be in the gene annotation of the \code{to}
species.  
}
\examples{
  ortholog( "PF3D7_0102800", from="PF3D7", to="PY17X")
  orthologSet( "PF3D7_0102800", from="PF3D7", to="PY17X")
}
