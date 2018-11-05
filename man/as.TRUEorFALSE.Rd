\name{as.TRUEorFALSE}
\alias{as.TRUEorFALSE}
\title{
Resolve 'anything' to a logical value
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
as.TRUEorFALSE(str)
}
\arguments{
  \item{str}{
almost any value, but typically a character string of length 1
}
}
\details{
Resolve the given value to a logical, using common ideas of T/F, Y/N, 0/1, etc.
Allows automatic environment variable substituion, via \code{\link{env.sub}}
}
\value{
logical value of length 1
}
\examples{
as.TRUEorFALSE( pi)
as.TRUEorFALSE( "$HOME")
}
