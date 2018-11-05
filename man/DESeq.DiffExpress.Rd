\name{DESeq.DiffExpress}
\alias{DESeq.DiffExpress}
\alias{EdgeR.DiffExpress}
\alias{RankProduct.DiffExpress}
\alias{SAM.DiffExpress}
\title{
Differential Expression by Group functions
}
\description{
Wrapper functions for calling various Differential Expression (DE) tools.
}
\usage{
DESeq.DiffExpress( fnames, fids, m = NULL, groupSet, targetGroup = sort(groupSet)[1], 
		geneColumn = "GENE_ID", intensityColumn = "READS_M", keepIntergenics = FALSE, 
		minimumRPKM = 1, missingGenes = "fill", extraColumn = NULL, 
		average.FUN = sqrtmean, wt.folds = 1, wt.pvalues = 1, wt.dists = 1)

DESeq.DiffExpress( fnames, fids, m = NULL, groupSet, targetGroup = sort(groupSet)[1], 
		geneColumn = "GENE_ID", intensityColumn = "READS_M", keepIntergenics = FALSE, 
		minimumRPKM = 1, missingGenes = "fill", extraColumn = NULL, 
		average.FUN = sqrtmean, wt.folds = 1, wt.pvalues = 1, wt.dists = 1)

RankProduct.DiffExpress( fnames, fids, m = NULL, groupSet, targetGroup = sort(groupSet)[1], 
		geneColumn = "GENE_ID", intensityColumn = "READS_M", keepIntergenics = FALSE, 
		offset = 0, missingGenes = "fill", extraColumn = NULL, 
		average.FUN = sqrtmean, wt.folds = 1, wt.pvalues = 1, wt.dists = 1)

SAM.DiffExpress( fnames, fids, m = NULL, groupSet, targetGroup = sort(groupSet)[1], 
		geneColumn = "GENE_ID", intensityColumn = "READS_M", keepIntergenics = FALSE, 
		minimumIntensity = NULL, missingGenes = "fill", extraColumn = NULL, 
		average.FUN = sqrtmean, wt.folds = 1, wt.pvalues = 1, wt.dists = 1,
		needJitter = FALSE)
}
\arguments{
  \item{fnames}{
a character vector of full pathnames to transcriptome files of gene expression
}
  \item{fids}{
a character vector of SampleIDs, one for each file
}
  \item{m}{
an optional matrix of gene expression, as an alternative to separate files for each sample.
In this case, the GeneIDs are taken from the rownames of \code{m} and the SampleIDs are
taken from the column names of \code{m}
}
  \item{groupSet}{
a character vector of group identifiers, one for each file, that categories that sample
}
  \item{targetGroup}{
a character string specifying the one group to be compared against all other groups
}
  \item{geneColumn}{
the name of the column in the transcriptome files that hold the GeneIDs
}
  \item{intensityColumn}{
the name of the column in the transcriptome files that hold the gene expression values
}
  \item{keepIntergenics}{
logical, should the named non-standard intergenic regions be kept in the final results
}
  \item{minimumRPKM, offset, minimumIntensity}{
linear offset added to all expression values, before taking log2 of the ratios, to prevent
divide by zero and false fold change magnitudes from very low expression values
}
  \item{missingGenes}{
how to handle missing genes when combining all files of transcriptomes.  
See \code{\link{expressionFileSetToMatrix}}
}
  \item{extraColumn}{
a character vector of additional column names to be propegated from the transcriptome files
to the resulting differential expression table
}
  \item{average.FUN}{
the math function for averaging gene expression values within each group
}
  \item{wt.folds, wt.pvalues, wt.dists}{
weighting terms for ordering the final results of differential expression.  Passed to 
function \code{\link{diffExpressRankOrder}}
}
  \item{needJitter}{
SAM can break if there is no variance in a group of expression values, as may happen in
some discrete data like proteome peptide abundance.  This logical can force the expression to
have some variance
}
}
\details{
This family of functions dispatch various differential expression tools, all of which return a
table of all genes, sorted into differential expression order such that genes most up-regulated 
in group \code{targetGroup} as compared to all other groups will be at the top, and most 
down-regulated genes at the bottom.
}
\value{
a data frame of differential results, with a row for every gene in the transcriptomes, sorted
into differential expression rank order.  Column include:
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Bob Morrison
}
\seealso{
\code{\link{metaRanks}} a low level tool for combining result files into a consensus ranking.
\code{\link{metaResults}} a wrapper function for combining results from all 5 DE functons into 
a final consensus differential expression table.
}
