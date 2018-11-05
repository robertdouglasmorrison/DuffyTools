\name{rankProductDiffExpress}
\alias{rankProduct}
\alias{rankProductDiffExpress}
\title{
Rank Product across multiple Samples
}
\description{
Use the Rank Product method to find differentially expressed genes,
or prioritize the order a set of genes.
}
\usage{
rankProduct(rankM, nSimulations = 500)

rankProductDiffExpress(fnames, groupSet, targetGroup = groupSet[1], 
		geneColumn = "GENE_ID", intensityColumn = "INTENSITY", 
		productColumn = "PRODUCT", offset = 0, keepIntergenics = FALSE, 
		average.FUN = logmean, poolSet = rep(1, length(fnames)), 
		nSimulations = 500, missingGenes = c("drop", "fill"))
}
\arguments{
  \item{rankM}{
numeric matrix of gene ranks, with GeneIDs as the rownames and SampleIDs as the
column names
}
  \item{fnames}{
character vector of full pathnames to existing transcriptome files
}
  \item{groupSet}{
character vector of GroupIDs or conditions, to categorize the transcripts
}
  \item{targetGroup}{
the one GroupID to be the chosen subset, to compare all other groups against.
This is the group that is being tested for up-regulation.
}
  \item{geneColumn}{
column name of the column of GeneIDs
}
  \item{intensityColumn}{
column name of the column of intensity values
}
  \item{productColumn}{
column name of the column that has gene product descriptions
}
  \item{offset}{
a linear offset to add to all intensity values to prevent divide by zero
and/or extreme fold change ratios
}
  \item{keepIntergenics}{
logical, explicity keep the non-genes, or drop them from consideration
}
  \item{average.FUN}{
the averaging function for combining gene intensities within subset groups,
gene ranks, and gene RP values
}
  \item{poolSet}{
numeric vector of sample pools or tiers.  For restricting 2-sample DE tests
to samples from comparable tiers. See details.
}
  \item{nSimulations}{
number of simulations of random permutations of the data, for calculating
false positives rates.
}
  \item{missingGenes}{
method for dealing with genes that are not present in every transcript file.
Either drop entire gene rows, or fill in with minimum observed intensity.
}
}
\details{
This function implements the Rank Product algorithm of Breitling, et.al.  By 
performing all possible 2-sample DE comparisons and ranking genes by fold change,
this calculates a family of rank positions for each gene.  Turning those ranks
into probabilities of differential expression, the algorithm assigns a measure
called Rank Product (RP), as the likelihood a gene could be that high in rank across
that many DE comparisons.

For the simple case of \code{rankProduct}, given a matrix of ranks, the 
algorithm just measures RP and estimates the false positive rates.

By default, each sample is compared to all other samples that are not from
its group.  If more restrictions are warranted, the \code{poolSet} argument 
can be used to assign a pool or tier to each sample; whereby only samples from
the same pool but coming from different groups go forward into the 2-sample tests.
}
\value{
For \code{rankProduct}, a data frame with RP values, average ranks, and false positive 
rates for each gene are returned, in the same row order as the input matrix.

For \code{rankProductDiffExpress}, a data frame of consensus gene differential expression, 
sorted by RP value, with columns:
  \item{GENE_ID }{ the genes, sorted from most up-regulated to most down-regulated}
  \item{PRODUCT }{ the gene product descriptions}
  \item{LOG_2_FOLD }{ the average fold change for each gene}
  \item{RP_VALUE }{ the Rank Product value.  See \code{\link{calcRP}}}
  \item{AVG_RANK }{ the average rank position over all 2-sample DE tests, for each gene}
  \item{AVG_SET1 }{ the average gene intensity over all samples in group \code{targetGroup}}
  \item{AVG_SET2 }{ the average gene intensity over all samples in the other groups}
  \item{E_VALUE }{ the expected number of genes to have an RP value this good, by chance}
  \item{FP_RATE }{ the rate of false positive DE genes, given this RP value}
}
\references{
Rainer Breitling, et.al.  FEBS Letters 573 (2004)
}
\author{
Bob Morrison
}
\note{
Typically, this function is called once for each group, to get all possible DE comparisons
between the various groups.  While the function explicitly measures up-regulation, by reversing
the order of the rows of the result, you get the answer for down-regulation.
}
