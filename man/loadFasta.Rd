\name{loadFasta}
\alias{loadFasta}
\alias{writeFasta}
\alias{writeLongFasta}
\title{
Read or write a FASTA File
}
\description{
Read/write a FASTA file from/to disk
}
\usage{
loadFasta( file = "file.fasta", verbose = TRUE, short.desc = TRUE)
writeFasta( fasta, file = NULL, line.width = 80)
writeLongFasta( desc, seq, file = NULL)
}
\arguments{
  \item{file}{
complete path name to an existing FASTA file for \code{loadFasta}, or new file name for writing
}
  \item{short.desc}{
logical, should the descriptor text be truncated at the first blank character
}
  \item{fasta}{
a fasta object, with 'desc' and 'seq' named elements
}
  \item{desc}{
character vector of descriptors
}
  \item{seq}{
character vector of sequences
}
  \item{line.width}{
integer value, the maximum line wrap length for written FASTA files
}
}
\details{
loadFasta is a wrapper around the Biostrings function 'readBStringSet'.  

The 'LongFasta' is a faster version that bypasses the line wrapping step, writing the entire sequence on a single line
}
\value{
For \code{loadFasta}, a list of length 2, with:
  \item{desc }{character vector of descriptors, from the '>' lines of the FASTA file}
  \item{seq }{character vector of sequences}
}

\seealso{
  \code{\link{as.Fasta}} for making fasta objects
  \code{\link{getFastaSeqFromFilePath}}, for getting the genomic sequence for a specific
  chrmosome.
}
