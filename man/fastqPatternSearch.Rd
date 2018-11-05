\name{fastqPatternSearch}
\alias{fastqPatternSearch}
\title{
Find DNA Patterns in FASTQ Reads
}
\description{
Find and count the frequency of DNA patterns in the reads of a FASTQ file
}
\usage{
fastqPatternSearch(filein, patterns, max.mismatch = 0, chunkSize = 4e+05)
}
\arguments{
  \item{filein}{
name of existing FASTQ file to search
}
  \item{patterns}{
character vector of DNA nucleotide patterns to search for
}
  \item{max.mismatch}{
maximum number of mismatching nucleotides allowed in a 'match'
}
  \item{chunkSize}{
buffer size, in number of reads
}
}
\details{
Uses the Biostrings package.  Useful for counting the number of reads
containing a known DNA sequence, such as adapter end primers.
}
\value{
numeric vector of length \code{length(patterns)}, giving the number of reads
that contain each pattern.
}

\seealso{
  \code{\link{matchPattern}}, in the Biostrings package.
}
