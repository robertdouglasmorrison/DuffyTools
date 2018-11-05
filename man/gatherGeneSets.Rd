\name{gatherGeneSets}
\alias{gatherGeneSets}
\title{
Gather Named Sets of Genes
}
\description{
A flexible way to gather named collections of gene sets for pathway enrichment and
differential expression.
}
\usage{
gatherGeneSets(geneSets, descriptor = "GeneSets", mode = c("combined", "separate"))
}
\arguments{
  \item{geneSets}{
a character vector of object names of pre-defined gene sets, as from pathway definitions;
or a list of vectors of gene names
}
  \item{descriptor}{
if \code{geneSets} is a list, the prefix to attach to each list element.
}
  \item{mode}{
if \code{geneSets} is a character vector, controls whether all geneSets are merged into
one large combined set.
}
}
\details{
This is a low level function shared by \code{link{geneSetEnrichment}} and other pathways 
tools.  It turns named pathway definition objects into the list(s) of genes needed for
pathway analysis.

When given a character vector of pathway object names, this relies on the current species
\code{\link{getCurrentSpecies}} to load the appropriate files.
}
\value{
a LIST with objects
  \item{n }{Number of lists of geneSets returned}
  \item{descriptor }{descriptor name for each list}
  \item{geneSets }{a list of geneSets, each ready for use by a pathway enrichment algorithm}
}
\seealso{
  See Also \code{\link{geneSetEnrichment}}, \code{link{pipe.GeneSetAnalysis}}
}
