\name{pipe.GeneSetAnalysis}
\alias{pipe.GeneSetAnalysis}
\title{
Pathway and Gene Module Analysis
}
\description{
Look for differentially expressed (DE) Pathways and Gene Modules from differential
gene expression results.
}
\usage{
pipe.GeneSetAnalysis(sampleIDset, speciesID = "Pf3D7", annotationFile = "Annotation.txt", 
	optionsFile = "Options.txt", results.path = NULL, folderName = "", 
	toolName = c("MetaResults", "RoundRobin", "RankProduct", "SAM", "DESeq", "EdgeR"), 
	geneMapColumn = if (speciesID %in% c("Hs_grc", "MacMu", "Mmu_grc")) "NAME" else "GENE_ID", 
	groupColumn = "Group", colorColumn = "Color", 
	geneSets = c("GO.BiologicalProcess", "GO.MolecularFunction", "GO.CellularComponent", 
	"KEGG.Pathways", "MetabolicPathways", "GeneProduct", "PBMC.GeneModules", "Blood.GeneModules"), 
	descriptor = "GeneSets", minGenesPerSet = if (speciesID %in% c("Hs_grc", "MacMu", "Mmu_grc")) 5 else 2, 
	mode = c("combined", "separate"), cutPvalue = 0.01, verbose = T)
}
\arguments{
  \item{sampleIDset}{
character vector of 2 or more SampleIDs.  This selects which groups/conditions are to be analysed.
}
  \item{speciesID}{
the species (organism) being evaluated
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
  \item{geneMapColumn}{
the name of the column in the current species' GeneMap that contains the GeneIDs found in the
GeneSets.  Some species use full notation GeneIDs and others uses simpler gene symbols.
}
  \item{groupColumn}{
the name of one column in the Annotation table, that specifies the group/condition for each sample.
}
  \item{colorColumn}{
the name of one column in the Annotation table, that specifies the color to use for each group/condition.
}
  \item{geneSets}{
either a character vector of names of pre-defined GeneSet objects, or a list of vectors of GeneIDs.
}
  \item{descriptor}{
a character string for naming the final folder of created analyses.
}
  \item{minGenesPerSet}{
the minimum number of genes that a GeneSet can contain.  Very low count sets are ignored.
}
  \item{mode}{
how to treat the various \code{geneSets}, whether to combine them all into one folder of joint analysis,
or to evaluate each GeneSet separately in its own subfolder.
}
  \item{cutPvalue}{
maximum P-value to include a GeneSet in the final results.  Lowering the P-value will produce less results.
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
\code{\link{pipe.GeneSetEnrichment}} for a similar function based on enrichment instead of mass density.
\code{\link{geneSetAnalysis}} for the lower-lever function using lists instead of SampleIDs.
}
