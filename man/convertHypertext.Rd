\name{convertHypertext}
\alias{convertHypertext}
\title{
Convert Hypertext special characters
}
\description{
Covert Hypertext special characters back to standard ASCII
}
\usage{
convertHypertext(txt)
}
\arguments{
  \item{txt}{
character string with embedded hypertext special characters
}
}
\value{
the character string, with hypertext converted back to standard ASCII
}
\examples{
  convertHypertext( "Hello\%20World")
}
