\name{openCompressedFile}
\alias{openCompressedFile}
\title{
Open a Compressed File
}
\description{
Open a (possibly) compressed file, that was compressed 
by either 'gzip' or 'bzip2'
}
\usage{
openCompressedFile(filename, open = "r")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
character string, the full pathname to one existing file
}
  \item{open}{
mode for opening, passed on to 'gzopen' etc.
}
}
\details{
A simple wrapper function to open a possibly compressed file, without
regard for which compression tool had been used.
}
\value{
a connection to the oped file.  If an error, the error as returned by
'file', 'gzfile', etc.
}
