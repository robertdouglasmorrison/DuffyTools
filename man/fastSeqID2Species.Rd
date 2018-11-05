\name{fastSeqID2Species}
\alias{fastSeqID2Species}
\title{
SeqID to SpeciesID Lookup
}
\description{
Lookup the SpeciesID and location in SeqMap for a set of SeqIDs
}
\usage{
fastSeqID2Species(seqid)
}
\arguments{
  \item{seqid}{
character vector of SeqIDs
}
}
\details{
Searches all SeqMaps in all defined MapSets, to find the SpeciesID
and SeqMap location for the given set of sequence IDs.
}
\value{
A list with 2 elements:
  \item{SPECIES}{ a character vector of the SpeciesID for each SeqID}
  \item{SMAP_PTR }{ a integer vector of pointers into the appropriate sequence
  map for each SeqID}
}

\note{
An error occurs if any \code{seqid} are not found in any of the defined MapSet annotation. 
}

\seealso{
  \code{\link{MapSets}}, for a detailed overview of organism annotation.
}
\examples{
  fastSeqID2Species( c("Pf3D7_02", "Hs_grc_11", "LmjF.36", "MT_H37RV_V2"))
}
