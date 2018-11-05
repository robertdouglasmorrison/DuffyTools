\name{pipe.GeneSetEnrichment}
\alias{pipe.GeneSetEnrichment}
\title{
Gene Module and Pathway Enrichment
}
\description{
Submit Differential Expression (DE) Results to Pathway and Gene Module Enrichment Analysis
}
\usage{
pipe.GeneSetEnrichment(sampleIDset, speciesID = "Pf3D7", annotationFile = "Annotation.txt", 
		optionsFile = "Options.txt", results.path = NULL, folderName = "", 
		toolName = c("MetaResults", "RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq"), 
		geneColumn = if (speciesID %in% c("Hs_grc", "MacMu", "Mmu_grc")) "GENE_NAME" else "GENE_ID", 
		groupColumn = "Group", 
		geneSets = c("GO.BiologicalProcess", "GO.MolecularFunction", "GO.CellularComponent", 
		"KEGG.Pathways", "MetabolicPathways", "GeneProduct", "PBMC.GeneModules", "Blood.GeneModules"), 
		descriptor = "Enrichment", maxPvalue = 0.05, wt.enrich = 1, wt.pvalue = 2, verbose = TRUE)
}
\arguments{
  \item{sampleIDset}{
character vector of one or more SampleIDs
}
  \item{speciesID}{
the one SpeciesID of results to operate on
}
  \item{results.path}{
top level folder holding all results, by default taken from the Options Table
}
  \item{folderName}{
the named subfolder of DE results to operate on, without any species prefix
}
  \item{toolName}{
the one DE tool of results to operate on.  Note that the quartet {results.path, toolName, speciesID, folderName} all combine
to specify the exact absolute path to the DE results files that will be used as input.
}
  \item{geneColumn}{
the name of the column in the various DE results files that contains the 
GeneID terms contained in the \code{geneSets}.  Note that this can vary
from species to species
}
  \item{groupColumn}{
the name of the column in the Annotation Table that contains the names of the conditions that were used for
differential expression (DE) analysis.  Note that this determines the specific filenames of DE results that
will be used as input.
}
  \item{geneSets}{
Either a chracter vector of pre-defined gene sets, or a list of gene sets where each elements name is a 
pathway/module and each element is a vector of GeneIDs from the given SpeciesID.
}
  \item{descriptor}{
a character string of length 1, that will be the name of the created subfolder of enrichment results
}
  \item{maxPvalue}{
P-value cutoff, to only report pathways/modules that have significant enrichment
}
  \item{wt.enrich}{
the relative weight to give to the enrichment magnitudes, for row-ordering the final results.
}
  \item{wt.pvalue}{
the relative weight to give to the enrichment P-values, for row-ordering the final results.
}
}
\details{
This function takes a folder of Differential Expression results, and searches for interesting pathways and
gene modules using the probability-based notion of enrichment.  See the basics about enrichment at 
\code{\link{enrichment}}, and a lower level intermediate function \code{\link{geneSetEnrichment}}.

Whereas the enrichment primitive requires a fixed number of genes as input, this function chooses a family of 
gene counts, calls enrichment for each, and summarizes the results as per \code{\link{metaRank.data.frames}}.
}
\value{
An entire subfolder of results is created, underneath the folder that contains all the results used
as input.  For each group named by \code{groupColumn}, an HTML file of enriched (up-regulated)
pathways and/or gene modules that are significantly over-represented in the DE results for that group of
samples.
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{pipe.GeneSetAnalysis}} for a similar idea of finding DE pathways and gene modules, but based 
  on a completely different heuristic.
}
