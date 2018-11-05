\name{readALNfile}
\alias{readALNfile}
\alias{readNEXUSfile}
\title{
Read Clustalw2 .ALN and .NXS files
}
\description{
Helper functions to read output files from Clustalw2, for multiple alignments
}
\usage{
readALNfile(file)
readNEXUSfile(file)
}
\arguments{
  \item{file}{
full pathname to existing Clustalw2 output file
}
}
\details{
Captures the multiple alignment text as a matrix of single characters.
}
\value{
a character matrix, with a single letter for each element.  The rownames
are the names of the sequences submitted to Clustal.
}
