\name{getHumanIDterms}
\alias{getHumanIDterms}
\title{
Common Human GeneID Terms
}
\description{
Gets useful human genome gene identifiers from GeneIDs
}
\usage{
getHumanIDterms(geneIDs)
}
\arguments{
  \item{geneIDs}{
character vector of GeneIDs
}
}
\details{
The rule about unique GeneIDs can make for some 'less than universal' gene
identifiers.  This function makes it easy to supplement a results data frame
with more universal gene terms.  Usually invoked as: \preformatted{
  tbl <- cbind( getHumanIDterms(tbl$GENE_ID), tbl)
}
}
\value{
A list with 2 elements:
  \item{GENE_NAME }{ the common gene symbol or nickname ortion of the GeneIDs}
  \item{ENTREZ_ID }{ the EntrezID number portion of the GeneIDs}
}
\note{
This function assumes that the common gene symbol and EntrezID are extractable
substrings of the given GeneIDs, in the pattern used for the human GeneIDs.
}
\seealso{
  Contrast with \code{\link{shortGeneName}}
}
\examples{
  getHumanIDterms( c("WASH7P:GI653635:01:14362","FAM138A:GI645520:01:34611"))
}
