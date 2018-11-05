\name{as.percent}
\alias{as.percent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format Numbers as Percentages 
}
\description{
 Format numeric values as character strings with (optional) trailing percent sign.

}
\usage{
as.percent(x, digits = 2, big.value = 1, percentSign = TRUE, sep = " ")
}
\arguments{
  \item{x}{
vector of numeric values, typically in the range 0..'big.value'
}
  \item{digits}{
desired number of significant digits after the decimal point
}
  \item{big.value}{
largest expected value.  All values get scaled so that 'big.value'
would be returned as '100 \%' 
}
  \item{percentSign}{
logical, should the trailing percent sign be appended
}
  \item{sep}{
character string to insert between the numeric percentage and the percent sign.
}
}
\details{
  Used to represent arbitrary values as a percentage in character format, where the
  value to be represented as 100 \% is given as 'big.value'.  Implemented via 'formatC'.
}
\value{
  A character vector of the same length as 'x', scaled such that 'big.value' equals
  100, with 'digits' significant digits and a trailing ' \%' appended.
}
\examples{
  as.percent( c( 7,21, 4), big.value=30, digits=2)
}
