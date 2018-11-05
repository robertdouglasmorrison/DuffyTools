\name{fastSeqPos2GenePos}
\alias{fastSeqPos2GenePos}
\title{
Sequence Positon to GeneID and Details
}
\description{
Map chromosomal sequences positions to their gene identifiers and positions
}
\usage{
fastSeqPos2GenePos(seqid, seqbase)
}
\arguments{
  \item{seqid}{
character vector of SeqID names
}
  \item{seqbase}{
integer vector of chromosomal sequence locations
}
}
\details{
Given a set of SeqIDs and locations, like from an alignment program's output,
map each of these to the appropriate GeneID, location within the gene, and that
gene's location within its organism's GeneMap.
}
\value{
A list with 4 elements:
  \item{GENE_ID }{ character vector of GeneIDs}
  \item{GENE_POSITION }{ integer vector of gene-centric locations (relative to
  the forward strand, regardless of the gene's actual coding strand)}
  \item{SPECIES }{ character vector of organism SpeciesIDs}
  \item{GMAP_PTR }{ integer vector of row pointers into each appropriate GeneMap}
}
\note{
This routine only searches MapSets for organisms in the current target.
}

\seealso{
  \code{\link{MapSets}}, for a detailed overview of organism annotation.
}
\examples{
  fastSeqPos2GenePos( c("Pf3D7_02", "Hs_grc_11", "Pf3D7_12"), c( 118000, 12345678, 55000))
}
