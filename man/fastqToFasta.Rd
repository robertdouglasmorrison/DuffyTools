\name{fastqToFasta}
\alias{fastqToFasta}
\title{
Convert FASTQ to FASTA
}
\description{
Convert a FASTQ file into a FASTA format file with quality score summary
}
\usage{
fastqToFasta(filein, fileout, Qscores = TRUE)
}
\arguments{
  \item{filein}{
existing file in FASTQ format
}
  \item{fileout}{
newly created file in FASTA format
}
  \item{Qscores}{
logical flag for adding summary quality score info to the ReadIDs
}
}
\value{
None.  Writes a new FASTA file to disk.
}
