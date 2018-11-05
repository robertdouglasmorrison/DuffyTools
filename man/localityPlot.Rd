\name{localityPlot}
\alias{localityPlot}
\title{
Plot Gene Locations on Chromosomes
}
\description{
Shows the geographic layout of gene subsets on the chromosomes of a genome
}
\usage{
localityPlot(geneSet, mode = c("separate", "combined"), color = 2, label = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneSet}{
character vector of GeneIDs, or a named list of such vectors
}
  \item{mode}{
draw N separate chromosomes or one 'combined' superimposed meta-chromosome
}
  \item{color}{
colors for each set of genes in \code{geneSet}
}
  \item{label}{
optional text for plot main title
}
}
\details{
This plot shows the visual locality of any gene subsets on the genome.  Currently, it
only supports one species at a time.
}
\value{
None.  Draws a plot.
}

\seealso{
  \code{\link{setCurrentSpecies}}, to change the current species.
}
