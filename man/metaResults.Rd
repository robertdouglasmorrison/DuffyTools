\name{metaResults}
\alias{metaResults}
\alias{metaResultsToHTML}
\title{
MetaRank Summary of Differential Expression Results 
}
\description{
Combine DE results from several DuffyTools into a meta ranking
of most up-regulated genes.
}
\usage{
metaResults(targetGroup, results.path = "results", speciesID = getCurrentSpecies(),
	geneColumn = "GENE_ID", subfolderName = "Group", 
	tools = c("RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq"), 
	altGeneMapLabel = NULL,
	rank.average.FUN = sqrtmean, value.average.FUN = mean, 
	keepIntergenics = FALSE, 
	topFolders = NULL, other.DE.files = NULL,
	missingGenes=c("drop","fill","na"), nFDRsimulations=0,
	direction = c("UP","DOWN"))

metaResultsToHTML( m, fileout, title = "", maxRows = 100,
		linkColumnNames = NULL, ...)
}
\arguments{
  \item{targetGroup}{
character string of the group of samples being tested for up-regulation.  See details.
}
  \item{results.path}{
folder that contains all the subfolders of DE files, as made by the various
DE tools.
}
  \item{speciesID}{
SpeciesID of the one species to be considered.
}
  \item{geneColumn}{
name of the column that contains the GeneIDs.
}
  \item{subfolderName}{
name of the subfolder for this set of analyses.  This folder must contain 
all the DE results from each of the various DE tools.
}
  \item{tools}{
character vector of the DE tools that were used, to be combined into meta results.
}
  \item{rank.average.FUN}{
function to use for combining the gene rankings from all DE tools into a meta rank.
}
  \item{value.average.FUN}{
function to use for combining the gene expression values from all DE tools into an average value.
}
  \item{keepIntergenics}{
logical, use only the true genes or keep all.
}
  \item{topFolders}{
a character vector (with names) of the various folders to be combined.  It is
expected that each of these locations contains a \code{results.path} folder.
\code{NULL} means that \code{results.path} is relative to the current working
directory.
}
  \item{other.DE.files}{
an optional character vector (with names) of full pathnames to additional DE files to include
in the meta rankings.
}
  \item{missingGenes}{
how to deal with genes not in all files.  See \code{\link{metaRanks}}.
}
  \item{nFDRsimulations}{
how many random simulations to perform for estimating the false discovery rate (FDR).
}
  \item{m}{
data frame result from metaResults functions, to be rendered as an HTML file
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
  \item{linkColumnNames}{
Character vector of column names, that will be converted to hyperlinks.  See \code{\link{table2html}}
}
  \item{...}{
Other arguments, passed to \code{\link{table2html}}
}
}
\details{
This tool combines multiple differential expression results from similar or
identical datasets that were evaluated by different DE tools, to generate a
consensus or meta result.  The various DE tools, such as \code{RoundRobin}, 
\code{RankProduct}, and \code{SAM}, all use different methods and can yield
somewhat inconsistent results; and this tool combines them all into
one consensus meta result.

The \code{targetGroup} is one of the GroupIDs from the specified grouping
column of the annotation file for that analysis.  Since each DE tool uses
a common naming convention for its results, all need files are locatable by
the triplet \code{subfolderName, speciesID, targetGroup}.
}
\value{
a data frame with columns for GeneIDs, products, and:
  \item{LOG2FOLD }{ the average fold change, as the arithmatic mean over all input DE files}
  \item{AVG_PVALUE }{ the average P-value, as the geometric mean over all input DE files}
  \item{AVG_RANK }{ the average rank for each gene, from all the input files, 
  using \code{average.FUN} for calculating the average}

  Also, one column for each input DE file, showing the rank of each gene in that file,
  with column names taken from the DE tools.

  If FDR simulation was done, then 2 more columns:
  \item{E_VALUE}{ the expected number of genes having an average rank as high as the given 
  ranks.  Expressed on the range 0 to N, where N is the number of genes in the results.}
  \item{FP_RATE}{ the rate of false positives, as estimated by random permutations of
  the given ranks.  Expressed on the range 0 to 1, tells the likelihood that a gene
  could have an AVG_RANK that good purely by chance}
}
\author{
Bob Morrison
}

\seealso{
\code{\link{RankProduct}}, \code{\link{SAM}}, \code{\link{RoundRobin}}, for
the various DE tools.

\code{\link{metaResultsToHTML}}, for direct creation of the HTML form of the
the resulting data frame

\code{\link{metaRanks}}, for averaging of gene rank positions.
}
