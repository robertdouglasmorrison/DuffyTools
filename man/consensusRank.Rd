\name{consensusRank}
\alias{consensusRank}
\title{
Consensus Ranking of Genes
}
\description{
Calculate the consensus ordering for genes from a set of transcripts
or differential expression (DE) files
}
\usage{
consensusRank(  fileSet, fids = basename(fileSet), 
		colorSet = 2:(length(fileSet) + 1), 
		geneColumn = "GENE_ID", 
		valueColumn = "LOG_2_FOLD", speciesID = "Pf3D7", 
		keepNonGenes = FALSE, weights = 1, altGeneMap = NULL, 
		annotationColumns = "PRODUCT", label = "", cex = 1, useLog = FALSE)
}
\arguments{
  \item{fileSet}{
characer vector of explicit filenames of transcripts or DE results
}
  \item{fids}{
sample IDs for each file
}
  \item{colorSet}{
colors for each file
}
  \item{geneColumn}{
the column name for the GeneIDs
}
  \item{valueColumn}{
the column name for the expression or DE value to use
}
  \item{speciesID}{
speciesID for these files
}
  \item{keepNonGenes}{
logical, should all rows or just real genes be used
}
  \item{weights}{
the relative weight to give to each file
}
  \item{altGeneMap}{
alternative GeneMap to use instead of the current MapSet
}
  \item{annotationColumns}{
characer vector of extra columns to carry forward from the 
current GeneMap to the output result
}
  \item{label}{
optional text for the plot heading
}
  \item{cex}{
plot point expansion factor
}
  \item{useLog}{
logical, should the Y axis be drawn in log scale
}
}
\details{
Calculate the consensus from a family of result files, either transcriptomes or
DE fold change results.  Measures the average rank and value (expression or fold
change) for every gene, and creates a new result that represents the (weighted) 
average of the given files.
}
\value{
a data frame of the genes in common to all files, sorted by rank and value.
The first 2 elements are the genes and values, with column names \code{geneColumn} 
and \code{valueColumn}.

The third column is the average RANK, by square-root mean as in \code{\link{metaRank}}

Lastly, any other columns from \code{annotationColumns} are carried forward from 
the GeneMap.
}
\note{
Only genes in common to ALL input files are used, thus the output data frame
may have less rows then the input files.
}
