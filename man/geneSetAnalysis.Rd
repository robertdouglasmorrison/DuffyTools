\name{geneSetDensity}
\alias{geneSetDensity}
\title{
Pathway and Gene Set Analysis using Density Curves
}
\description{
Look for differentially expressed (DE) Pathways and Gene Sets from differential
gene expression results.
}
\usage{
geneSetDensity(deList, geneSets, speciesID = getCurrentSpecies(), colorset = c(2:(length(deList) + 1)), 
	optionsFile = "Options.txt", results.path = NULL, folderName = "", 
	toolName = c("MetaResults", "RoundRobin", "RankProduct", "SAM", "DESeq", "EdgeR"), 
	descriptor = "GeneSet", minGenesPerSet = 2, trimGenesToGeneMap = TRUE,
	geneMapColumn = if (speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
	cutPvalue = 0.05, cutRankShift = NULL, cutFold = 0.1, cutFDR = 0.05,
	makePlots = TRUE, PLOT.FUN = NULL, makeGeneTables = TRUE,
	doFDR = TRUE, NgeneSets = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{deList}{
a named list of differential expression (DE) data frames (2 or more), containing the DE results 
for each condition.  The names of the list should be the names of the conditions.
}
  \item{geneSets}{
either a character vector of names of pre-defined GeneSet objects, or a list of vectors of GeneIDs.
}
  \item{speciesID}{
the species (organism) being evaluated
}
  \item{colorset}{
a vector of colors for displaying each condition
}
  \item{results.path}{
top level folder of results, defaults to the value in the Options table
}
  \item{folderName}{
the name of the subfolder of DE comparisons to use
}
  \item{toolName}{
the name of the DE tool used to generate the differential gene expression results.
}
  \item{descriptor}{
a character string for naming the final folder of created analyses.
}
  \item{minGenesPerSet}{
the minimum number of genes that a GeneSet can contain.  Very low count sets are ignored.
}
  \item{trimGenesToGeneMap}{
logical, removes extra gene terms from any gene set references that are not in the current gene map.
}
  \item{geneMapColumn}{
the name of the column in the current species' GeneMap that contains the GeneIDs found in the
GeneSets.  Some species use full notation GeneIDs and others uses simpler gene symbols.
}
  \item{cutPvalue}{
maximum P-value to include a GeneSet in the final results.  Lowering the P-value will produce less results.
}
  \item{cutRankShift}{
minimum rank shift to include a GeneSet in the final results.  Raising the ranks shift will produce less results.
}
  \item{cutFold}{
minimum fold change to include a GeneSet in the final results.  Raising the fold change will produce less results.
}
  \item{cutFDR}{
maximum FDR to include a GeneSet in the final results.  Lowering the FDR will produce less results. Only
used when \code{doFDR=TRUE}.
}
  \item{doFDR}{
logical, should a False Discovery Rate be calculated, and added to all results.  Rather slow.
}
  \item{makePlots}{
logical, should plot images for each GeneSet be created
}
  \item{PLOT.FUN}{
optional function to be called for making gene set plots.  Must take one argument that is a vector of gene IDs.
}
  \item{makeGeneTables}{
logical, should small HTML files of gene-to-GeneSet membership tables be created
}
}
\details{
This function analyzes how pre-defined sets of genes are differentially expressed between 2 or more conditions.
The methodology is to treat each gene as a point mass, located at its rank order location in each of the 
condition's DE results, and then summarize all the genes in each GeneSet as a density mass for the set.  A gene set
is called differentially expressed in a condition if its density mass is significantly shifted with respect to that
same gene set in all other conditions.
}
\value{
A folder of results, where the final name/location is a combination of 
{results.path, toolName, speciesID, folderName, descriptor}.

In that folder will be a family of HTML files for each condition, one HTML and .txt file ocvering all
GeneSets and conditions, and a subfolder of all plots and hyperlinked gene/module files.

Briefly, for each condition, there are HTML files for the most up-regulated and down-regulated GeneSets,
with P-values and density rank shifts for each significantly DE GeneSet.  Within each HTML file are links
to images of density plots for each GeneSet, and to gene membership tables to examine the details about genes
in that GeneSet.
}
\author{
Bob Morrison
}
\seealso{
\code{\link{geneSetEnrichment}} for a similar function based on enrichment instead of mass density.
\code{\link{pipe.GeneSetAnalysis}} for a high-lever wrapper function using SampleIDs.
}
