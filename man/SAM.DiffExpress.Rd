\name{SAM.DiffExpress}
\alias{SAM.DiffExpress}
\title{
Use SAM for Differential Expression
}
\description{
Use the SAM method to find differentially expressed genes from a
set of transcriptome files.
}
\usage{
SAM.DiffExpress(fnames, fids, groupset, targetGroup = sort(groupset)[1], 
		geneColumn = "GENE_ID", intensityColumn = "INTENSITY", 
		keepIntergenics = FALSE, minimumIntensity = NULL, 
		missingGenes = "fill", 
		wt.folds = 1, wt.pvalues = 1, wt.dists = 1, ...)
}
\arguments{
  \item{fnames}{
character vector of full pathnames to existing transcriptome files
}
  \item{fids}{
character vector of SampleIDs
}
  \item{groupset}{
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
  \item{keepIntergenics}{
logical, explicity keep the non-genes, or drop them from consideration
}
  \item{minimumIntensity}{
the 'fudge factor' called 'S0' for SAM, that represents a minimum
intensity added to all denominators to prevent division by zero and 
effectively scale the distances 'D' produced by SAM
}
  \item{missingGenes}{
method for dealing with genes that are not present in every transcript file.
Either drop entire gene rows, or fill in with minimum observed intensity.
}
  \item{wt.folds, wt.pvalues, wt.dists}{
weight terms for influencing the final ordering of DE genes.  See
\code{\link{SAM.diffExpressRankOrder}}
}
}
\details{
This function implements the SAM algorithm of Tusher, et.al., as implemented
in R package 'siggenes' by H. Schwender.
}
\value{
A data frame of consensus gene differential expression, 
sorted by consensus of up-regulation, with columns:
  \item{GENE_ID }{ the genes, sorted from most up-regulated to most down-regulated}
  \item{PRODUCT }{ the gene product descriptions}
  \item{LOG_2_FOLD }{ the average fold change for each gene}
  \item{P_VALUE }{ the P-value for each gene}
  \item{Q_VALUE }{ the Q-value for each gene}
  \item{AVG_SET1 }{ the average gene intensity over all samples in group \code{targetGroup}}
  \item{AVG_SET2 }{ the average gene intensity over all samples in the other groups}
}
\references{
V. Tusher, et.al.,  PNAS 98 (2001)
}
\author{
Bob Morrison
}
\note{
Typically, this function is called once for each group, to get all possible DE comparisons
between the various groups.  While the function explicitly measures up-regulation, by reversing
the order of the rows of the result, you get the answer for down-regulation.
}

\seealso{
  \code{\link{rankProductDiffExpress}}, for the Rank Product method.
}
