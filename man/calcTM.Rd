\name{calcTM}
\alias{calcTM}
\title{
Melting Temperature (Tm) of an oligo probe
}
\description{
Calculate the Tm of a set of DNA nucleotides
}
\usage{
calcTM(probes)
}
\arguments{
  \item{probes}{
character vector of DNA nucleotide sequences
}
}
\details{
Calculates the melting temperature of the given oligos.
Any 'N's in the sequence cause a value of NA for that oligo.
}
\value{
A numeric vector of the same length as \code{probes}
}
\references{
by Joseba Bikandi,  http://biophp.org
}
\note{
Implemented as a Perl script, so requires write and execute permissions
on the current directory.  Creates a file 'MeltingTemp.pl'.
}
