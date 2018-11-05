\name{AAtoCodonOptimizedDNA}
\alias{AAtoCodonOptimizedDNA}
\title{
  Convert Amino Acid sequence to Codon Optimized DNA
}
\description{
Do the reverse translation from amino acids back to DNA nucleotides,
using a reference DNA sequence as a source of codon usage frequency.
}
\usage{
AAtoCodonOptimizedDNA(aa, dnaRef)
}
\arguments{
  \item{aa}{
A character string of valid amino acids.
}
  \item{dnaRef}{
A character string of DNA, in frame, usually large, that represents the typical 
codon usage frequency.
}
}
\details{
Since the mapping from AA to DNA is not unique, the user provides a chunk of DNA reference,
that will be converted to AA and used to model the frequency distribution for randomly selecting
DNA codons from each amino acid.
}
\value{
a character string, of length \code{ nchar(aa) * 3}, of DNA that is a valid reverse translation
of the given amino acid string.
}
\author{
Bob Morrison
}
\note{
Passes the frequency distribution of codons found in \code{dnaRef} to \code{sample}, to randomly
select codons with a usage frequency similar to the given reference DNA.
}
\seealso{
See Also \code{\link{DNAtoAA}}, \code{\link{getCodonMap}}
}
