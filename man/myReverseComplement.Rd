\name{myReverseComplement}
\alias{myReverseComplement}
\alias{myReverse}
\title{
Reverse Complement of DNA nucleotides
}
\description{
Reverse Complement of DNA nucleotides
}
\usage{
myReverseComplement(dna)
myReverse(dna)
}
\arguments{
  \item{dna}{
character string of DNA bases
}
}
\details{
Reverse compliment of DNA.  Similar to the same functions in package Biostrings,
but operates on the R 'character' data type.
}
\value{
a character string, that is the reverse complement (or just the reverse) of
\code{dna}.
}
\examples{
  myReverseComplement( "ATGAAACCCGGGTTT")
}
