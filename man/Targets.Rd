\name{Targets}
\alias{Targets}
\title{
Target Organisms
}
\description{
Overview of 'DuffyTools' Target Organism(s)
}
\details{
We use the notion of 'target' to mean a small set of organisms or species
that are to be treated as a single unit for functional genomics analysis.
Our lab's focus is on host parasite interactions and most of our datasets
involve measuring gene expression of both the host and the parasite in a
single experiment, be it microarrays, proteomics, or NextGen sequencing.

A target definition consists of 3 parts:  a TargetID, a set of 1 or more
SpeciesIDs, and an equally long set of species file prefixes.  In this
way, the tools can operate on the datasets as a single entity while 
simultaneously reading, using, and creating all necessary species-specific 
information.

For the trivial case of a target being a single species, it is customary
for the TargetID to exactly match the SpeciesID.  The file prefix is
typically a shortened version of the SpeciesID, with enough uniqueness
over the set of all possible species of interest.

For a detailed overview of species annotations, see \code{\link{MapSets}}. 
for understanding the subtle difference between a target and a species,
also see \code{\link{Species}}. 
}

\seealso{
  \code{\link{getCurrentTarget}}, for getting details about the current target.
  \code{\link{getAllTargets}}, for getting details about all predefined targets.
  \code{\link{setCurrentTarget}}, to change the current target. 
  \code{\link{addTarget}}, to add a new target definition.
  \code{\link{exportTargets}}, \code{\link{importTargets}}, to managing the
  set of all known targets.
}
