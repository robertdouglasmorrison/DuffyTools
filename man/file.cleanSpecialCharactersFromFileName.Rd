\name{file.cleanSpecialCharactersFromFileName}
\alias{file.cleanSpecialCharactersFromFileName}
\title{
Clean Special Character from Filenames
}
\description{
Clean non-universally acceptable special characters out of a filename,
to guarantee the file will be successfully created and readable later with
the exact same filename.
}
\usage{
file.cleanSpecialCharactersFromFileName(fnames)
}
\arguments{
  \item{fnames}{
character vector of potential filenames
}
}
\details{
This function is designed to make sure that filenames with embedded GeneIDs
always use characters allowable on any operating system.
}
\value{
character vector of (potentially modified) filenames, of the same length
as \code{fnames}, with any special characters converted to hyphens.
}
 
\examples{
  file.cleanSpecialCharactersFromFileName( "MAL8_18s:rRNA.png")
}
