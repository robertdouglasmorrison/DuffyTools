\name{getSamplePairIDs}
\alias{getSamplePairIDs}
\alias{originalSamplePairID}
\title{
SampleIDs for Paired End Datasets
}
\description{
Utilities to map between paired end specific files and their parent SampleID.
}
\usage{
getSamplePairIDs(sampleID, annotationFile = "Annotation.txt")

originalSamplePairID(sampleID, annotationFile = "Annotation.txt")
}
\arguments{
  \item{sampleID}{
a SampleID, that may or may not have a paired end suffix (e.g. \code{"_1", "_2"})
}
  \item{annotationFile}{
the sample annotation file containing the details about this sample
}
}
\details{
With paired end data, there are some pipeline tools that operate on each paired end
file half separately.  These need the pair suffix present.  Most other tools just want
the 'parent' SampleID with no suffix.  These utilities map back and forth between those
two formats.
}
\value{
For \code{getSamplePairIDs()}, a vector of at most length 2, with pair specific SampleIDs
for this sample.  If the dataset is not paired end data, the given SampleID is returned
as is.

For \code{originalSamplePairID}, a single SampleID with any pair suffix removed.
}
