\name{readALNfile}
\alias{readALN}
\alias{writeALN}
\alias{readNEXUS}
\title{
Read and Write Clustal .ALN and .NXS files
}
\description{
Helper functions to read and/or write the files from Clustal, for Multiple Sequence Alignments (MSA)
}
\usage{
readALN(file)
readNEXUS(file)

writeALN( aln, outfile, line.width = 60, title, capitalize = NULL, blankMissingFlanks = FALSE, 
		fastaToo = FALSE)
}
\arguments{
  \item{file, outfile}{
full pathname to existing Clustal output file.  Or a full pathname for a file to be written to.
}
  \item{aln}{
a MSA alignment object, as returned ba MAFFT.  See \code{\link{mafft}}.
}
  \item{line.width}{
integer number of characters to write per line of output.
}
  \item{capitalize}{
Optional.  Either a logical, or a character string in {"uppper", "lower"}, controlling how the 
aligned characters are case converted or not.
}
  \item{blankMissingFlanks}{
logical, should the hyphens that denote gap characters at the very beginning and ending
flanks of the MSA, which are only due to missing sequence instead of true insertions or deletions,
be represented as blanks instead. 
}
  \item{fastaToo}{
logical, should the sequences in the \code{aln} object also be written out in FASTA format.  Note
that these FASTA sequences will preserve any gap hyphens present in the MSA into the FASTA constructs.
}
}
\details{
These function facilitate manipulation of MSA results as R objects.
}
\value{
For \code{readNEXUS}, a character matrix, with a single letter for each element.  The rownames
are the names of the sequences submitted to Clustal.

For \code{readALN}, a list:
  \item{alignment}{ a character matrix with a single letter for each element.}
  \item{consensus}{ the 'consensus string' from below the MSA elements.}
  \item{pct.conserved}{ numeric vector of the amount of sequence similarity along the MSA, 
  expressing each column position as the percentage of rows that contain the most commonly
  observed character.}

For \code{writeALN}, one or two files written to disk.
}
