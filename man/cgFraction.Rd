\name{cgFraction}
\alias{cgFraction}
\title{ 
Percentage of C and G bases in a DNA sequence 
}
\description{
Computes the percentage of C and G bases in a DNA sequence 
}
\usage{
cgFraction(dna)
}
\arguments{
  \item{dna}{ 
a vector of character strings representing nucleic DNA sequences
}
}
\value{
  A vector of numeric values in 0 to 1, giving the percentage of C and G bases.
}
\examples{
cgFraction( "ACGTACGTAAGGCCT")
}

