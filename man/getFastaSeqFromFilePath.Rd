\name{getFastaSeqFromFilePath}
\alias{getFastaSeqFromFilePath}
\title{
Genomic DNA sequence for one SeqID
}
\description{
Get the genomic DNA nucleotide sequence for one chromosome from a file or folder
}
\usage{
getFastaSeqFromFilePath(filePath, seqID, verbose = FALSE)
}
\arguments{
  \item{filePath}{
complete path name to a FASTA file or a folder of FASTA files as a character string.
}
  \item{seqID}{
character string of one SeqID sequence identifier.
}
  \item{verbose}{
logical, show status of the search progress.
}
}
\details{
A (hopefully easy) function to find and read the DNA sequence for any single
chromosome.  In the case of a small genome contained in a single file, this
finds the named chromosome descriptor in the FASTA file, and returns that
nucleotide sequence.  In the case of very large genomes with separate files for
each chromosome, this searches the given path for the appropriately named
FASTA file.  There are hardcoded assumptions that the file names and 
sequence descriptors in the FASTA file(s) correspond exacly with the SeqID terms
in the SeqMap.
}
\value{
one chromosome of DNA nucleotides, as a Biostrings \code{\link{DNAString}},
or \code{NULL} if the file/SeqID search fails.
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{loadFasta}}, for reading arbitray FASTA files.
  \code{\link{MapSets}}, for details on synchronizing SpeciesIDs, SeqIDs, and
  FASTA descriptors.
}
