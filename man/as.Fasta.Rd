\name{as.Fasta}
\alias{as.Fasta}
\title{
Convert to FASTA list format
}
\description{
Convert vectors of descriptors and sequences to a list
}
\usage{
as.Fasta(desc, seq)
}
\arguments{
  \item{desc}{
characer vector of sequence descriptors
}
  \item{seq}{
character vector of sequences (DNA or protein)
}
}
\value{
a list of length 2, with elements:
  \item{desc }{ descriptors }
  \item{seq }{ sequences }
}

\seealso{
  \code{\link{loadFasta}}, for reading FASTA files 
  \code{\link{writeFasta}}, for writing FASTA files 
  \code{\link{toTextFasta}}, for formatting FASTA sequences to a fixed width 
  \code{\link{as.Fasta.data.frame}}, for conversion to data frame format 
}
