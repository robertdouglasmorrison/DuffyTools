\name{env.sub}
\alias{env.sub}
\title{
Environment Variable Substitution
}
\description{
Perform environment variable substitution on a character string
}
\usage{
env.sub(x)
}
\arguments{
  \item{x}{
character string of length 1
}
}
\details{
Performs shell-like substitution of environment variables in a character
string.  Expects variable syntax in the Unix format:  ${MY_VAR}

Undefined variable names are left unexpanded.
}
\value{
a character string with all '$' escaped environment variables replaced
by their value
}
\seealso{
  \code{\link{Sys.getenv}}
}
\examples{
  file.info ( env.sub( "${HOME}/.Rprofile"))
}
