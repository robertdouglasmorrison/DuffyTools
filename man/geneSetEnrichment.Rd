\name{geneSetEnrichment}
\alias{geneSetEnrichment}
\title{
Gene Pathway Enrichment
}
\description{
Find the pathways enriched for in a set of genes.
}
\usage{
geneSetEnrichment(genes, geneSets = defaultGeneSets(), 
		upOnly = TRUE, minEnrich = 2, maxPvalue = 0.05, wt.enrich = 1, wt.pvalue = 2,
		geneUniverse = NULL, reportGenes = FALSE, dropPathNameLinks = FALSE, verbose = TRUE)
}
\arguments{
  \item{genes}{
a vector of genes (GeneIDs for most parasite species,  Gene Symbols for most mammals),
that represent the genes of interest, perhaps from differential expression analysis.
}
  \item{geneSets}{
a named pathway definitions object as a character vector, or a list of vectors of genes.
Default is all pathway sets defined for the current species.
}
  \item{upOnly}{
logical, report only over-enriched pathways or both up and down.
}
  \item{minEnrich}{
minimum enrichment to report.  Higher value reports less pathways. A value of 1.0 means not 
enriched at all.
}
  \item{maxPvalue}{
maximum P-value to report.  Lower value reports less pathways.
}
  \item{wt.enrich}{
weight to give the enrichment component in the final sorted ordering.
}
  \item{wt.pvalue}{
weight to give the P-value component in the final sorted ordering.
}
  \item{geneUniverse}{
the universe of potentially detected genes.  Defaults to all the genes in the geneMap for
the current species.  Useful when only a subset of genes are eligble.
}
  \item{reportGenes}{
logical, report just the pathway information, or also return a second data frame of genes.
}
  \item{dropPathNameLinks}{
logical, keep any hyperlinks present in the pathway name, or remove them.
}
}
\details{
Based on the hypergeometric function.  See \code{\link{enrichment}}.  Evaluates every 
named set of genes in \code{geneSets} to see how many of the elements of \code{genes} are 
present, and if that number is greater or less than expected by chance.  Tabulated for all
and then trimmed to the subset meeting the reporting criteria.
}
\value{
A list of two data frames, one called \code{pathways} and one called \code{genes}, 
when \code{reportGenes} is true.  Else a single data frame of pathways.
The pathway object is the enriched pathway groups, sorted by both enrichment and P-value.  
The most unexpected pathway groups will be at the top.  The genes object is all the genes 
present in any enriched pathway, with all its pathways memberships, sorted by GeneID.
The pathway data frame:
  \item{PathName }{name of the pathway}
  \item{Enrichment }{the enrichment is the ratio of \code{Pct_Given/Pct_Total}}
  \item{P_value }{the likelihood of getting that many given genes in a pathway of that size}
  \item{N_Total }{total number of genes (given the geneUniverse) in that pathway}
  \item{Pct_Total }{percentage of all genes that are in that pathway}
  \item{N_Given }{number of given genes in \code{genes} that were found in that pathway}
  \item{Pct_Given }{percentage of given genes that were found in that pathway}
}
\author{
  Bob Morrison
}
\note{
  The default pathway objects are distributed as part of the DuffyTools data directory.
  Not all defined species have all types of pathway files at this time.
}

\seealso{
 \code{\link{pipe.GeneSetEnrichment}}, for a higher level wrapper for evaluating 
 enrichment of multiple differential expression results.
 \code{\link{pipe.GeneSetDensity}} and \code{\link{pipe.MetaGeneSets}}, for differential expression of pathways.
}
\examples{
  setCurrentSpecies( "Hs_grc")
  geneSetEnrichment( c("ACVR1","BAMBI","CDK8","COL1A2","ENG","FSHB","GDF9","JUNB","LTBP1","MYC","SKI","SMAD5","TGFB1","UBB"), maxP=0.0001)
}
