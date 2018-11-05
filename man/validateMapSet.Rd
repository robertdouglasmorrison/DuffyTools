\name{validateMapSet}
\alias{validateMapSet}
\title{
MapSet Self-Consistency Validation
}
\description{
A low level function to assure that the annotation maps within
a MapSet are self and cross consistent
}
\usage{
validateMapSet(mapset, checkOverlaps = FALSE)
}
\arguments{
  \item{mapset}{
a MapSet object, as returned by \code{\link{importMapSet}}
}
  \item{checkOverlaps}{
logical, extra test to flag gene and/or exon boundary overlaps
as potential annotation errors
}
}
\value{
Fails with a stopping error if maps are not consistent.
}
