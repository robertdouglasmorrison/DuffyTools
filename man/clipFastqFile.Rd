\name{clipFastqFile}
\alias{clipFastq}
\title{
Clip Bases from FASTQ reads
}
\description{
Clip (trim) bases from the 5' and/or 3' ends of reads in a FASTQ file
}
\usage{
clipFastq(filein, fileout, clip5prime = 0, clip3prime = 0)
}
\arguments{
  \item{filein}{
Existing FASTQ file of reads and quality scores
}
  \item{fileout}{
New created FASTQ file with shorter reads and quality scores
}
  \item{clip5prime}{
number of bases to trim from the 5' end (high quality, front)
}
  \item{clip3prime}{
number of bases to trim from the 3' end (low quality, tail)
}
}
\details{
Overwrites the new output file if it already exists.  Either file
may be compressed if the given file suffix ends with .gz
}
\value{
No return value.  Creates a new FASTQ file with the same number of
reads.
}
