\name{getFileLineCount}
\alias{getFileLineCount}
\alias{FileLineCount}
\alias{quickFileLineCountLookup}
\alias{quickFileLineCountRecord}
\title{
Manange Line and Read Counts for Large Files
}
\description{
Routines to keep track of the number of lines and/or short reads
in large files, typically NextGen sequencing data
}
\usage{
getFileLineCount(f, sampleID = "", verbose = TRUE, mode = "full", 
		 what = c("lineCount", "readCount"))

quickFileLineCountLookup(f, sampleID = "", what = c("lineCount", "readCount"))
quickFileLineCountRecord(f, sampleID = "", lineCount, readCount = lineCount)
}
\arguments{
  \item{f}{
the name of one file, that we want to find or record line/read count information about
}
  \item{sampleID}{
an optional SampleID, to specify which FileLineCount file to use
}
  \item{mode}{
character string, denotes if we need the actual exact number of lines, or just
enough to confirm that the file is not empty.  See details.
}
  \item{what}{
the type of count, either actual lines of text, or the number of reads/sequences
}
  \item{lineCount}{
the number of actual lines of text
}
  \item{readCount}{
the number of actual reads/sequences
}
}
\details{
These routine help manage the very large files typical of NextGen data.  They help
facilitate file buffering, and are helpful for generating summary statistics about
alignment results.

The 'get' and 'lookup' routines return one numeric value, giving the size of the file
in 'what' units.  When 'mode' is not 'full', it returns the count from the first buffer.
For files that are not found, or other errors, zero is returned.
The 'record' routine updates the FileLineCount file to the true counts, and is typically
called after a computationally intense step, to let other later functions know how big
the various files are without having to explicitly count records.
}
\note{
The counts are stored in files called '.FileLineCount.<SampleID>'.  If abnormal or
unexpected behaviors occur with a particular sample, try deleting its FileLineCount file;
a new one will get re-created over time.
}
