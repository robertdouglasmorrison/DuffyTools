\name{expressionCluster}
\alias{expressionCluster}
\title{
Clustering of Gene Expression Data
}
\description{
Perform cluster analysis on gene expression data, usually transcriptomes.
}
\usage{
expressionCluster(x, useLog = TRUE, normalize = TRUE, FUN = diana)
}
\arguments{
  \item{x}{
numeric matrix of gene expression data, with a row for each gene/oligo and
a column for each sample
}
  \item{useLog}{
logical, should the data be log transformed before clustering
}
  \item{normalize}{
logical, should the data be normalize before clustering, by scaling to make all
columns have equal total expression.
}
  \item{FUN}{
clustering function to use
}
}
\details{
Uses the \code{cluster} package, to perform a hierarchical clustering of the
given samples.  
}
\value{
a cluster object of type dependent on the \code{FUN} provided, suitable for
plotting, etc.
}

\seealso{
  \code{\link{agnes}}, \code{\link{diana}}, for clustering methods.  
  \code{\link{expressionFileSetToMatrix}}, for turning sets of transcript files into one matrix.
}
