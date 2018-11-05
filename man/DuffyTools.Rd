\name{DuffyTools}
\alias{DuffyTools}
\alias{DuffyTools-package}
\title{
Duffy Lab "Tools" package Overview
}
\description{
A suite of R tools from the lab of Patrick Duffy, providing a common
infrastucture for our microarray (DuffyMA), NextGen sequencings (DuffyNGS),
and proteomics (DuffyProteome) packages.
}
\details{
This package provides a variety of utilities and data manipulation
functions that are common to the Duffy Lab's R packages for functional genomics 
spanning the domains of microarrays (DuffyMA), next gen sequencing (DuffyNGS),
and proteomics (DuffyProteome).

Primary focus is on providing a universal and flexible method of working with
genomic feature annotations, from any number of organisms, in a user controllable
format.  This allows for transcriptional studies of organisms with 
evolving annotations, such as many pathogens and parasites. Further, the 
platform is explicitly designed for mixed-organism data
collection platforms such as microarrays with both host and parasite probes, or
RNA-seq of mixed genetic material such as field samples drawn from infected
individuals.

For details about annotations, including adding new organism definitions or
customizing existing annotations, see \code{\link{MapSets}}.  For details about 
setting the target organism(s) for analysis, including adding new mixed-organism
capabilities, see \code{\link{Targets}}.  For details on switching between
organisms (species) during analysis, see \code{\link{Species}}.

A second major focus is on tools that produce and operate on transcription and 
differential expression (DE) data sets.  Utilities supporting operations such as 
clustering, plotting, and comparisons are included.  Methods for detecting
DE genes based on Round Robin, Rank Product, EdgeR, DESeq, and SAM, are supported.

To support this focus, the package tries to migrate as many user settable
parameters as feasible into configuration files that are local to each experiment.
By convention, we divide these parameters into 2 categories:  ones that can
be thought of as sample specific are found in an \code{\link{AnnotationTable}}, and
ones that are common to the entire experiment are found in an \code{\link{OptionsTable}}. 
In this manner, much of the behavior of the lab's functional genomics packages can
be customized by the end user without altering the underlying R code.

This package also includes a variety of functions for mathematic operations,
file manipulation in general, and specific tools for FASTA and FASTQ formats.
Lastly, any other tools that are generic enough to be shared by our functional
genomics packages tend to be thrown into the DuffyTools package.
}
\author{
Bob Morrison  bob.morrison@seattlebiomed.org
              robert.morrison@nih.gov
}
