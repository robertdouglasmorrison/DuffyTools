\name{Species}
\alias{Species}
\title{
Species Overview
}
\description{
Overview of 'DuffyTools' notion of Species
}
\details{
We use the notion of 'species' to mean a single organism or strain
that has it's own unique genomic annotation.  Contrast this with
\code{\link{Targets}}.

Our lab's focus is on host parasite interactions and most of our datasets
involve measuring gene expression of both the host and the parasite in a
single experiment, be it microarrays, proteomics, or NextGen sequencing.
To accomplish this, we need to work with multiple species in a seamless and
coordinated manner.

A species definition consists of a MapSet -- a collection of objects that
completely describes the genomic annotation -- and 2 unique identifiers:
a SpeciesID and a species FilePrefix.  See \code{\link{MapSets}} for full 
details.  The SpeciesID is used to tell the software which MapSet 
(annotation) to use, and the FilePrefix is used for the reading and writing
of all species-specific files.  Typically, the species file prefix is a 
shortened version of the SpeciesID, unique over the set of all possible 
species of interest.

In the trivial case of operating on datasets from just a single species, it 
is customary for the SpeciesID to exactly match the TargetID.  
}

\seealso{
  \code{\link{getCurrentSpecies}}, for getting details about the current species.
  \code{\link{setCurrentSpecies}}, to change to a different species.

  \code{\link{getCurrentGeneMap}}, to get the table of genes.
  \code{\link{getCurrentExonMap}}, to get the table of coding exons.
  \code{\link{getCurrentSeqMap}}, to get the table of sequences (chromosomes).
  \code{\link{getCurrentRrnaMap}}, to get the table of ribosomal RNA and other
special case genes.
}
