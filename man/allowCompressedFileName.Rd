\name{allowCompressedFileName}
\alias{allowCompressedFileName}
\title{
Allow .gz filename extensions
}
\description{
Find a named file, allowing for the possibility it was 'gzip' compressed
}
\usage{
allowCompressedFileName(filename)
}
\arguments{
  \item{filename}{
full pathname of one existing file, with or without .gz suffix
}
}
\details{
Confirms the existance of the named file for reading.  If not found,
searches for the file with a compression suffix of .gz or .bz2
}
\value{
the filename, with appropriate suffix attached if needed.
}
