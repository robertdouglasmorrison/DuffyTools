\name{significantGenes}
\alias{significantGenes}
\title{
Select Significant Genes from a Differential Expression Results file.
}
\description{
Select the subset of genes that meet both minimum fold change and maximum P-value criteria,
from a file that wwas produced by a Differential Expression tool.
}
\usage{
significantGenes(file, geneColumn = "GENE_ID", foldColumn = "LOG2FOLD", pvalueColumn = "PVALUE", 
		keepIntergenics = FALSE, min.fold = 1, max.pvalue = 0.05, min.genes = NULL, 
		shortNames = TRUE, sep = "\t", verbose = TRUE)
}
\arguments{
  \item{file}{
character string of one full path filename.
}
  \item{geneColumn}{
character string of one column name, that contains the GeneIDs.
}
  \item{foldColumn}{
character string of one column name, that contains the fold change values.
}
  \item{pvalueColumn}{
character string of one column name, that contains the significance values.
}
  \item{keepIntergenics}{
logical, should non-genes be excluded from consideration.
}
  \item{min.fold}{
numeric value for the inclusion cutoff based on fold change.
}
  \item{max.pvalue}{
numeric value for the inclusion cutoff based on significance.
}
  \item{min.genes}{
optional numeric value, giving the minimum number of genes to return.
}
  \item{shortNames}{
%%     ~~Describe \code{shortNames} here~~
logical, should the returned GeneIDs be shortened to their gene symbol only.
}
  \item{sep}{
file field separator, passed to \code{\link{read.table}}.
}
}
\note{
This function is intended to be generic with respect to the various differential expression
tools.  But as some DE tools use slightly different column naming conventions, care should be taken 
to assure the correct column names are used for each DE tool.
}
\value{
A data frame of significant genes:
  \item{GENE_ID }{the identifiers of the significant genes}
  \item{PRODUCT }{their gene product terms}
  \item{LOG2FOLD }{their fold change values}
  \item{PVALUE }{their significance values}
}
\author{
Bob Morrison
}
