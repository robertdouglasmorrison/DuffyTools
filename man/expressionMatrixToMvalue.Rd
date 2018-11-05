\name{expressionMatrixToMvalue}
\alias{expressionMatrixToMvalue}
\title{
Turn a matrix of Gene Expression into M-values
}
\description{
Translate raw gene abundance data into log2 fold change data commonly called M-values
}
\usage{
expressionMatrixToMvalue(x, average.FUN = median, minIntensity = 0, small.offset = 1)
}
\arguments{
  \item{x}{
Matrix of gene expression data, as form \code{\link{expressionFileSetToMatrix}}
}
  \item{average.FUN}{
The averaging function to apply to each row of \code{m}, to get the denominator for 
dividing the row elements
}
  \item{minIntensity}{
The minimum allowed value for elements of \code{m}, to put a floor on low expression values
}
  \item{small.offset}{
The linear offset to add to all row elements, before division and log2 transform, to 
prevent divide-by-zero results and to minimize extreme fold changes from very low expression values
}
}
\details{
After M-value transformation, those row elements that are greater than the average expression will be positive,
while elements lower than the average expression will be negative values.  A value of 1.0 means 2-fold higher
expression and a value of -2.0 means 4-fold lower expression.
}
\value{
A matrix of M-value fold change data, of the same dimensions and names as \code{m}, suitable for heatmaps, etc.
}
\note{
This function is called by 'expressionHeatmap()'
}
