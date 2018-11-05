\name{mergeTables}
\alias{mergeTables}
\alias{mergeIntegerTables}
\title{
Merge 1-D Tables
}
\description{
Merge two 1-D tables of counts
}
\usage{
mergeTables(t1, t2)
mergeIntegerTables(t1, t2)
}
\arguments{
  \item{t1, t2}{
named vectors of counts, as from \code{base::table}
}
}
\details{
These functions combine tabular data.  The 'Integer' form requires that the
names be interpretable as integer values.  
}
\value{
a named vector of class 'table', where the named elements are the union of
the names of the two input tables, and their values are the sum of counts.
}
