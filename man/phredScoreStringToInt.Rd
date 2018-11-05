\name{phredScoreStringToInt}
\alias{phredScoreStringToInt}
\title{
Convert Phred Quality Score Strings to Integers
}
\description{
Turn the cryptic Phred quality scores from FASTQ files, etc., into integers
}
\usage{
phredScoreStringToInt(txt, scoreType = "Phred33")
}
\arguments{
  \item{txt}{
character vector of Phred quality score strings
}
  \item{scoreType}{
type of offset used for the ASCII coded scores.  Currently, only
'Phred33' and 'Phred66' are supported.
}
}
\details{
Most efficient for strings of equal length.
}
\value{
an integer matrix with NR rows and NC columns, where NR = \code{length(txt)}
and NC = \code{max(nchar(txt))}.  Each element is the Phred score of one
DNA base call.
}
\examples{
  phredScoreStringToInt( c( "_ABBBB&&.#$&", "IIIIIIIIIIII"))
}
