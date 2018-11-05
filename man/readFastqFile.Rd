\name{readFastqFile}
\alias{readFastqFile}
\alias{writeFastqFile}
\title{
Read and write FASTQ Files
}
\description{
Read a FASTQ file into a data frame, or write a data frame to a FASTQ file.
}
\usage{
readFastqFile(filein, maxReads = NULL, verbose = TRUE)
writeFastqFile(x, fileout, compress = FALSE, verbose = TRUE)
}
\arguments{
  \item{filein}{
full pathname to one existing, optionally compressed, FASTQ file.  Compression and .gz file extension are
automatically detected.
}
  \item{fileout}{
full pathname for newly created, optionally compressed, FASTQ file.  Compression and .gz file extension are
automatically coordinated.
}
  \item{maxReads}{
integer count of the maximum number of reads wanted, where \code{NULL} 
means read the entire file.
}
  \item{compress}{
logical, explicitly compress the created file and adjust filename accordingly.
}
}
\value{
For \code{readFastqFile}, a data frame, with 3 columns:
  \item{READ_ID }{the ReadID for each short read}
  \item{READ_SEQ }{the nucleotide sequence of each short read}
  \item{SCORE }{the cryptic ASCII Phred score string for each short read}

For \code{writeFastqFile}, the name of the created file.
}
\note{
FASTQ data can get quite large.  This routine can easily exceed available memory for
large \code{maxReads}.
}
