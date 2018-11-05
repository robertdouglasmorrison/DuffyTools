\name{shortGeneName}
\alias{shortGeneName}
\title{
Shorten GeneIDs for Simplicity
}
\description{
Depending on the species, GeneIDs can get rather cumbersome to 
meet the uniqueness requirement.  This function trims GeneIDs down
a leading subset of terms for easier viewing and display.
}
\usage{
shortGeneName(genes, split = ":", keep = 2)
}
\arguments{
  \item{genes}{
character vector of GeneIDs
}
  \item{split}{
split argument for dividing the GeneID into terms.  Passed to \code{\link{strsplit}}
}
  \item{keep}{
integer number of leading terms to keep for each GeneID
}
}
\value{
a character vector of shortened GeneIDs.  These may no longer map uniquely to
exactly one gene, and are typically used for their easier appearance.
}
\seealso{
  Contrast with \code{\link{getHumanIDterms}}
}
\examples{
  shortGeneName( c("WASH7P:GI653635:01:14362","FAM138A:GI645520:01:34611"))
}
