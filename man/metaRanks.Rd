\name{metaRanks}
\alias{metaRanks}
\alias{metaRank.data.frames}
\alias{metaRank2html}
\title{
Meta Ranking of Transcriptome Genes
}
\description{
Combine several transcriptome or differential expression results to
build one meta result
}
\usage{
metaRanks(fnames, fids, weightset = rep(1, length(fileset)), 
		geneColumn = "GENE_ID", valueColumn="LOG_2_FOLD",
		pvalueColumn = "PVALUE", productColumn="PRODUCT", sep="\t",
		rank.average.FUN = sqrtmean, value.average.FUN = mean,
		keepIntergenics = FALSE,
		missingGenes = c("drop","fill","na"), missingValue = 0,
		naDropPercent = 0.5, nFDRsimulations=0)

metaRank.data.frames(df.list, weightset = rep(1, length(fileset)), 
		geneColumn = "GENE_ID", valueColumn="LOG_2_FOLD",
		pvalueColumn = "PVALUE", productColumn="PRODUCT", 
		rank.average.FUN = logmean, value.average.FUN = mean,
		missingGenes = c("drop","fill","na"), missingValue = 0,
		naDropPercent = 0.5)

metaRank2html( tbl, fileout = "metaRank.html", title = "", maxRows = 100,
		valueColumn = "LOG2FOLD", ...)
}
\arguments{
  \item{fnames}{
character vector of full path names to existing transcript or DE files
}
  \item{fids}{
character vector of SampleIDs for each file
}
  \item{df.list}{
list of data frame objects, that are all tables of gene abundance information.
This is just an alternative format to the file-based inputs.
}
  \item{weightset}{
numeric vector of weights for each file
}
  \item{geneColumn}{
column name that contains the GeneIDs
}
  \item{valueColumn}{
column name that contains the fold change or expression value
}
  \item{pvalueColumn}{
column name that contains the fold change or expression value
}
  \item{productColumn}{
column name that contains the fold change or expression value
}
  \item{rank.average.FUN}{
averaging function for combining the rank position of each gene
}
  \item{value.average.FUN}{
averaging function for combining the data value terms for each gene
}
  \item{keepIntergenics}{
logical, keep the explicit 'non-genes' or drop them from the result
}
  \item{missingGenes}{
mode for dealing with genes that are not found in every data set.  
\code{drop} drops the entire row/gene from the output.  
\code{fill} keeps the row/gene; and assigns the largest (worst) rank to
every missing gene.
\code{na} flags the gene's rank as \code{NA}, and then drops entire rows having
genes missing from more than \code{ naDropPercent \%} of the datasets.
}
  \item{missingValue}{
value to use for the value term for missing genes
}
  \item{naDropPercent}{
minimum percentage of missing observations, to cause a gene to be dropped
completely from the results
}
  \item{nFDRsimulations}{
number of trials of randomized trials of permuted ranks, to estimate the
false positive rate
}
  \item{tbl}{
data frame result from metaRank functions, to be rendered as an HTML file
}
  \item{fileout}{
Full path name for the newly created HTML file.
}
  \item{title}{
Character string on length 1, to be used as header and title in the HTML file
}
  \item{maxRows}{
Maximum number of rows to include in the HTML file
}
  \item{...}{
Other arguments, passed to \code{\link{table2html}}
}
}
\details{
This function combines several results files, while 'averaging' the results for
each gene from all the separate files.  All files are assumed to have the genes in
a consistent ranked order, by either expression (transcripts), fold change 
(differential expression - DE), etc.  For each gene, the average rank position and the 
average observed value are calculated, and that gene's rank position in each
input file is reported.  Final gene order is by average rank.
}
\value{
A data frame that represents the average of the input files, along with each
gene's rank position in every file. 

Also, one last column \code{FP_RATE} of the simulated false discovery rates,
which estimate how likely it was to see an average rank as high as that gene's,
by chance.

Lastly, a correlation estimate of the pairwise Spearman rank correlation is measured
and saved in the global workspace with the object name of \code{metaRankCC}.
}
\author{
Bob Morrison
}
\note{
The false discovery rate step is not currently provided in the 'data.frame' version.

The final ranking can be sensitive to the choice of rank averaging function. The default
of \code{\link{sqrtmean}} is a compromise between linear and geometric averaging.
}
