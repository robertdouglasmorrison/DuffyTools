\name{oneDtable}
\alias{oneDtable}
\title{
Table from a 1-D vector
}
\description{
A slightly faster Table function for 1-D Vectors
}
\usage{
oneDtable(x)
}
\arguments{
  \item{x}{
Vector of values to be tabled
}
}
\details{
A simplified, slightly faster alternative to \code{\link{table}}, for making 1-D tables from large
vectors
}
\value{
a table, with names for the unique values of \code{x}, with values being the number of times
each value of \code{x} was present in \code{x}
}
