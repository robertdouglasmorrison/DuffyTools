\name{getCodonMap}
\alias{getCodonMap}
\title{
Get DNA Codon & Amino Acid Map
}
\description{
Get the conversion table of DNA codon triplets to amino acids
}
\usage{
getCodonMap()
}
\value{
A data frame with 4 columns:
  \item{DNA }{the DNA triplet codons}
  \item{RNA }{the corresponding RNA triplet codons}
  \item{AA }{the single letter amino acid code}
  \item{AA_Long }{the 3-letter amino acid code}
}

\seealso{
  \code{\link{DNAtoAA}}, \code{\link{DNAtoBestPeptide}}, \code{\link{DNAtoFrameShiftingPeptides}}, 
  \code{\link{AAtoCodonOptimizedDNA}}, 
  for converting between nucleotide strings and protein sequences.
}
